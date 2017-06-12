/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TestExamples.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "BaseLib/BidirectionalMap.h"

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/LinearEquation/DenseLinearEquation.h"

#include "GeoLib/Rectangle.h"

#include "MeshLib/Core/IMesh.h"

#include "DiscreteLib/Serial/DiscreteSystem.h"
#include "DiscreteLib/Core/IElemenetWiseLinearEquationLocalAssembler.h"

#include "NumLib/Function/TXFunction.h"

#include "FemLib/Function/FemFunction.h"
#include "FemLib/BC/DirichletBC2FEM.h"
//#include "FemLib/BC/FemDirichletBCMethod.h"
#include "FemLib/BC/NeumannBC2FEM.h"
#include "FemLib/Post/Extrapolation.h"


class GWFemTest
{
public:
    typedef FemLib::FemNodalFunctionScalar<DiscreteLib::DiscreteSystem>::type MyNodalFunctionScalar;
    typedef FemLib::FEMIntegrationPointFunctionVector<DiscreteLib::DiscreteSystem>::type MyIntegrationPointFunctionVector;
    GeoLib::Rectangle *rec;
    DiscreteLib::DiscreteSystem *dis;
    MeshLib::IMesh *msh;
    NumLib::ITXFunction *_K;
    MyNodalFunctionScalar *head;
    MyIntegrationPointFunctionVector *vel;
    std::vector<size_t> vec_bc1_nodes;
    std::vector<double> vec_bc1_vals;
    std::vector<size_t> vec_bc2_nodes;
    std::vector<double> vec_bc2_vals;
    FemLib::LagrangeFeObjectContainer* _feObjects;

    ~GWFemTest()
    {
        delete rec;
        delete _feObjects;
        delete head;
        delete vel;
        delete dis;
    }

    void define(MeshLib::IMesh *msh, NumLib::ITXFunction *K=0)
    {
        //#Define a problem
        //geometry
        rec = new GeoLib::Rectangle(GeoLib::Point(0.0, 0.0, 0.0),  GeoLib::Point(2.0, 2.0, 0.0));
        const GeoLib::Polyline &poly_left = rec->getLeft();
        const GeoLib::Polyline &poly_right = rec->getRight();
        //mesh
        this->msh = msh;
        dis = new DiscreteLib::DiscreteSystem(msh);
        _feObjects = new FemLib::LagrangeFeObjectContainer(msh);
        //discretization
        head = new MyNodalFunctionScalar();
        head->initialize(*dis, FemLib::PolynomialOrder::Linear);
        head->setFeObjectContainer(_feObjects);
        vel = new MyIntegrationPointFunctionVector();
        vel->initialize(dis);
        //bc
        NumLib::TXFunctionConstant f_bc1(.0);
        FemLib::DirichletBC2FEM bc1(*msh, poly_right, f_bc1, vec_bc1_nodes, vec_bc1_vals);
        NumLib::TXFunctionConstant f_bc2(-1e-5);
        FemLib::NeumannBC2FEM bc2(*msh, 0, *_feObjects, poly_left, f_bc2, vec_bc2_nodes, vec_bc2_vals);
        // mat
        _K = (K!=0) ? K : new NumLib::TXFunctionConstant(1.e-11);
    }

    static void calculateHead(GWFemTest &gw)
    {
        const MeshLib::IMesh *msh = gw.msh;
//        for (size_t i=0; i<gw.vec_bc1.size(); i++) gw.vec_bc1[i]->setup();
//        for (size_t i=0; i<gw.vec_bc2.size(); i++) gw.vec_bc2[i]->setup();
        // global EQS
        const size_t n_dof = msh->getNumberOfNodes();
        MathLib::EigenDenseLinearEquation eqs;
        eqs.create(n_dof);

//        MathLib::EigenDenseLinearEquation::MatrixType* globalA = eqs.getA();
        double* globalRHS = eqs.getRHS();

        //assembly
        FemLib::IFeObjectContainer* feObjects = gw.head->getFeObjectContainer();
        MathLib::LocalMatrix localK;
        std::vector<size_t> e_node_id_list;
        double gp_x[3], real_x[3];
        for (size_t i_e=0; i_e<msh->getNumberOfElements(); i_e++) {
            MeshLib::IElement *e = msh->getElement(i_e);
            FemLib::IFiniteElement *fe = feObjects->getFeObject(*e);
            const size_t &n_dof = fe->getNumberOfVariables();
            localK = MathLib::LocalMatrix::Zero(n_dof, n_dof);
            FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
            for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
                q->getSamplingPoint(j, gp_x);
                fe->computeBasisFunctions(gp_x);
                fe->getRealCoordinates(real_x);
                double k = .0;
                gw._K->eval(real_x, k);
                MathLib::LocalMatrix mat_k(1,1);
                mat_k(0,0) = k;
                fe->integrateDWxDN(j, mat_k, localK);
            }
            e->getNodeIDList(e_node_id_list);
            eqs.addAsub(e_node_id_list, localK);
        }

        //outputLinearEQS(*globalA, globalRHS);

        //apply BC
        for (size_t i=0; i<gw.vec_bc2_nodes.size(); i++) {
            globalRHS[gw.vec_bc2_nodes[i]] -= gw.vec_bc2_vals[i];
        }
        //outputLinearEQS(globalA, globalRHS);
//        FemLib::DiagonalizeMethod diag;
//        diag.apply(gw.vec_bc1_nodes, gw.vec_bc1_vals, eqs);
        eqs.setKnownX(gw.vec_bc1_nodes, gw.vec_bc1_vals);
        //for (size_t i=0; i<gw.vec_bc1.size(); i++) gw.vec_bc1[i]->apply(eqs);
        //outputLinearEQS(*globalA, globalRHS);

        //solve
        eqs.solve();

        //update head
        gw.head->setNodalValues(eqs.getX(), 0, eqs.getDimension());
    }

    static void calculateVelocity(GWFemTest &gw)
    {
        const MeshLib::IMesh *msh = gw.msh;
        FemLib::IFeObjectContainer* feObjects = gw.head->getFeObjectContainer();
        //calculate vel (vel=f(h))
        for (size_t i_e=0; i_e<msh->getNumberOfElements(); i_e++) {
            MeshLib::IElement* e = msh->getElement(i_e);
            FemLib::IFiniteElement *fe = feObjects->getFeObject(*e);
            MathLib::LocalVector local_h(e->getNumberOfNodes());
            for (size_t j=0; j<e->getNumberOfNodes(); j++)
                local_h[j] = gw.head->getValue(e->getNodeID(j));
            // for each integration points
            FemLib::IFemNumericalIntegration *integral = fe->getIntegrationMethod();
            double r[2] = {};
            const size_t n_gp = integral->getNumberOfSamplingPoints();
            gw.vel->setNumberOfIntegationPoints(i_e, n_gp);
            std::vector<double> xi(e->getNumberOfNodes());
            std::vector<double> yi(e->getNumberOfNodes());
            for (size_t i=0; i<e->getNumberOfNodes(); i++) {
                const GeoLib::Point* pt = msh->getNodeCoordinatesRef(e->getNodeID(i));
                xi[i] = (*pt)[0];
                yi[i] = (*pt)[1];
            }
            for (size_t ip=0; ip<n_gp; ip++) {
                MathLib::LocalVector q(2);
                q[0] = .0;
                q[1] = .0;
                integral->getSamplingPoint(ip, r);
                fe->computeBasisFunctions(r);
                const MathLib::LocalMatrix *dN = fe->getGradBasisFunction();
                //MathLib::Matrix<double>*N = fe->getBasisFunction();
                std::vector<double> xx(2);
                fe->getRealCoordinates(&xx[0]);

                double k;
                gw._K->eval(&xx[0], k);
                //dN->axpy(-k, &local_h[0], .0, &q[0]); //TODO  q = - K * dN * local_h;
                q = -k * (*dN) * local_h;
                gw.vel->setIntegrationPointValue(i_e, ip, q);
            }
        }
    }
};



struct DiscreteExample1
{
    std::vector<size_t> list_dirichlet_bc_id;
    std::vector<double> list_dirichlet_bc_value;
    static const size_t dim_eqs = 9;
    std::vector<double> exH;

    DiscreteExample1()
    {
        size_t int_dirichlet_bc_id[] = {2,5,8,0,3,6};
        list_dirichlet_bc_id.assign(int_dirichlet_bc_id, int_dirichlet_bc_id+6);
        list_dirichlet_bc_value.resize(6);
        fill(list_dirichlet_bc_value.begin(), list_dirichlet_bc_value.begin()+3, .0);
        fill(list_dirichlet_bc_value.begin()+3, list_dirichlet_bc_value.end(), 1.0);
        exH.resize(9);
        for (size_t i=0; i<9; i++) {
            if (i%3==0) exH[i] = 1.0;
            if (i%3==1) exH[i] = 0.5;
            if (i%3==2) exH[i] = 0.;
        }
    }

    void setLocalDirichletBC(const BaseLib::BidirectionalMap<size_t, size_t> &map_global2localNodeId, std::vector<size_t> &local_dirichlet_bc_id, std::vector<double> &local_dirichlet_bc_value)
    {
        for (size_t i=0; i<list_dirichlet_bc_id.size(); i++) {
            if (map_global2localNodeId.countInA(list_dirichlet_bc_id[i])>0) {
                size_t local_id = map_global2localNodeId.mapAtoB(list_dirichlet_bc_id[i]);
                local_dirichlet_bc_id.push_back(local_id);
                local_dirichlet_bc_value.push_back(list_dirichlet_bc_value[i]);
            }
        }
    }


    class TestElementAssembler : public DiscreteLib::IElemenetWiseLinearEquationLocalAssembler
    {
        MathLib::Matrix<double> _m;
    public:
        TestElementAssembler()
        {
            _m.resize(4,4);
            _m(0,0) = 4.0; _m(0,1) = -1.0; _m(0,2) = -2.0; _m(0,3) = -1.0; 
            _m(1,1) = 4.0; _m(1,2) = -1.0; _m(1,3) = -2.0;
            _m(2,2) = 4.0; _m(2,3) = -1.0;
            _m(3,3) = 4.0;
            for (size_t i=0; i<4; i++)
                for (size_t j=0; j<i; j++) _m(i,j) = _m(j,i);
            _m *= 1.e-11/6.0;
        }
        void assembly(MeshLib::IElement &/*e*/, MathLib::DenseLinearEquation &eqs)
        {
            (*eqs.getA()) = _m;
        }
    };
};
