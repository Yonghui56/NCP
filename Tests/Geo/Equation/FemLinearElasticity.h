/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemLinearElasticity.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/TransientAssembler/IElementWiseTransientLinearEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientResidualLocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "Tests/Geo/Material/PorousMedia.h"
#include "Tests/Geo/Material/Solid.h"

namespace Geo
{

template <class T_MATRIX>
inline void setNu_Matrix(const size_t dim, const size_t nnodes, const T_MATRIX &N, T_MATRIX &matN)
{
    matN *= .0;
    for (size_t in=0; in<nnodes; in++) {
        size_t offset = in*dim;
        for (size_t k=0; k<dim; k++)
            matN(k,offset+k) = N(in);
    }
}

template <class T_MATRIX>
inline void setB_Matrix(const size_t dim, const size_t nnodes, const T_MATRIX &dN, T_MATRIX &matB)
{
    matB *= .0;
    if (dim==2) {
        for (size_t in=0; in<nnodes; in++) {
            const size_t offset = in*dim;
            const double dshape_dx = dN(0, in);
            const double dshape_dy = dN(1, in);
            matB(0,offset) = dshape_dx;
            matB(1,offset+1) = dshape_dy;
            matB(3,offset) = dshape_dy;
            matB(3,offset+1) = dshape_dx;
        }
    } else if (dim==3) {
        for (size_t in=0; in<nnodes; in++) {
            const size_t offset = in*dim;
            const double dshape_dx = dN(0, in);
            const double dshape_dy = dN(1, in);
            const double dshape_dz = dN(2, in);
            matB(0,offset) = dshape_dx;
            matB(1,offset+1) = dshape_dy;
            matB(2,offset+2) = dshape_dz;
            matB(3,offset) = dshape_dy;
            matB(3,offset+1) = dshape_dx;
            matB(4,offset) = dshape_dz;
            matB(4,offset+2) = dshape_dx;
            matB(5,offset+1) = dshape_dz;
            matB(5,offset+2) = dshape_dy;
        }
    }
}

template <class T_MATRIX>
inline void setB_Matrix(const size_t dim, double dshape_dx, double dshape_dy, double dshape_dz, T_MATRIX &B_matrix)
{
    switch(dim)
    {
    case 2:
        // B_11, dN/dx
        (*B_matrix)(0,0) = dshape_dx;
        // B_12, 0.0
        (*B_matrix)(0,1) = 0.0;
        // B_21, 0.0
        (*B_matrix)(1,0) = 0.0;
        // B_22, dN/dy
        (*B_matrix)(1,1) = dshape_dy;
        // B_31, 0.0
        (*B_matrix)(2,0) = 0.0;
        // B_32, 0.0
        (*B_matrix)(2,1) = 0.0;
        // B_41, dN/dy
        (*B_matrix)(3,0) = dshape_dy;
        // B_42, dN/dx
        (*B_matrix)(3,1) = dshape_dx;

        break;
    case 3:
        // B_11, dN/dx
        (*B_matrix)(0,0) = dshape_dx;
        // B_22, dN/dy
        (*B_matrix)(1,1) = dshape_dy;
        // B_33, dN/dz
        (*B_matrix)(2,2) = dshape_dz;
        //
        // B_41, dN/dy
        (*B_matrix)(3,0) = dshape_dy;
        // B_42, dN/dx
        (*B_matrix)(3,1) = dshape_dx;
        //
        // B_51, dN/dz
        (*B_matrix)(4,0) = dshape_dz;
        // B_53, dN/dx
        (*B_matrix)(4,2) = dshape_dx;
        //
        // B_62, dN/dz
        (*B_matrix)(5,1) = dshape_dz;
        // B_63, dN/dy
        (*B_matrix)(5,2) = dshape_dy;

        break;
    }
}

template <class T_MATRIX>
inline void setB_Matrix4axisymmetry(const size_t dim, double radius, double dshape_dx, double dshape_dy, T_MATRIX &B_matrix)
{
    // B_11, dN/dx
    (*B_matrix)(0,0) = dshape_dx;
    // B_12, 0.0
    (*B_matrix)(0,1) = 0.0;
    // B_21, N/r
    (*B_matrix)(1,0) = dshape_dx / radius;
    // B_22, 0.0
    (*B_matrix)(1,1) = 0.0;
    // B_31, 0.0
    (*B_matrix)(2,0) = 0.0;
    // B_32, dN/dz
    (*B_matrix)(2,1) = dshape_dy;
    // B_41, dN/dy
    (*B_matrix)(3,0) = dshape_dy;
    // B_42, dN/dx
    (*B_matrix)(3,1) = dshape_dx;
}

#if 1
class FemLinearElasticLinearLocalAssembler: public NumLib::IElementWiseTransientLinearEQSLocalAssembler
{
private:
    PorousMedia* _pm;
    FemLib::LagrangeFeObjectContainer* _feObjects;
public:
    typedef MathLib::LocalVector LocalVectorType;
    typedef MathLib::LocalMatrix LocalMatrixType;

    FemLinearElasticLinearLocalAssembler(FemLib::LagrangeFeObjectContainer &feObjects, PorousMedia &pm)
    : _pm(&pm), _feObjects(&feObjects)
    {
    };

    virtual ~FemLinearElasticLinearLocalAssembler() {};

    virtual void assembly(const NumLib::TimeStep &/*time*/,  const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &/*localDofManager*/, const LocalVectorType &/*local_u_n1*/, const LocalVectorType &/*local_u_n*/, MathLib::LocalEquation &eqs)
    {
        FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);

        const size_t dim = e.getDimension();
        const size_t n_strain_components = (dim==2 ? 4 : 6);
        const size_t nnodes = e.getNumberOfNodes();
        Solid *solidphase = _pm->solidphase;

        // set D
        LocalMatrixType matD = LocalMatrixType::Zero(n_strain_components, n_strain_components);
        //matD *= .0;
        double nv = solidphase->poisson_ratio;
        double E = solidphase->Youngs_modulus;
        double Lambda, G, K;
        calculateLameConstant(nv, E, Lambda, G, K);
        setElasticConsitutiveTensor(dim, Lambda, G, matD);

        // body force
        LocalVectorType body_force = LocalVectorType::Zero(dim);
        bool hasGravity = false;
        if (hasGravity)
            body_force[dim-1] = solidphase->density * 9.81;

        //
        LocalMatrixType &localK = *eqs.getA();
        LocalVectorType &localRHS = *eqs.getRHSAsVec();
        LocalMatrixType matB = LocalMatrixType::Zero(n_strain_components, nnodes*dim);
        LocalMatrixType matN = LocalMatrixType::Zero(dim, nnodes*dim);
        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
            q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);
            double fac = fe->getDetJ() * q->getWeight(j);


            // set N,B
            LocalMatrixType &N = *fe->getBasisFunction();
            LocalMatrixType &dN = *fe->getGradBasisFunction();
            setNu_Matrix(dim, nnodes, N, matN);
            setB_Matrix(dim, nnodes, dN, matB);

            // K += B^T * D * B
            localK.noalias() += fac * matB.transpose() * matD * matB;

            // RHS
            localRHS.noalias() += fac * matN.transpose() * body_force;
        }
    }
};

class FemLinearElasticResidualLocalAssembler : public NumLib::IElementWiseTransientResidualLocalAssembler
{
private:
    PorousMedia* _pm;
    FemLib::LagrangeFeObjectContainer* _feObjects;
public:
    FemLinearElasticResidualLocalAssembler(FemLib::LagrangeFeObjectContainer &feObjects, PorousMedia &pm)
    : _pm(&pm), _feObjects(&feObjects)
    {
    };

    virtual ~FemLinearElasticResidualLocalAssembler() {};

    /// assemble a local residual for the given element
    /// @param time            time step
    /// @param e            element
    /// @param local_u_n1    guess of current time step value
    /// @param local_u_n    previous time step value
    /// @param local_r        local residual
    virtual void assembly(const NumLib::TimeStep &/*time*/,  const MeshLib::IElement &/*e*/, const DiscreteLib::DofEquationIdTable &/*localDofManager*/, const MathLib::LocalVector &/*local_u_n1*/, const MathLib::LocalVector &/*local_u_n*/, MathLib::LocalVector &/*local_r*/)
    {

    }
};

class FemLinearElasticJacobianLocalAssembler: public NumLib::IElementWiseTransientJacobianLocalAssembler
{
private:
    PorousMedia* _pm;
    FemLib::LagrangeFeObjectContainer* _feObjects;
public:
    FemLinearElasticJacobianLocalAssembler(FemLib::LagrangeFeObjectContainer &feObjects, PorousMedia &pm)
    : _pm(&pm), _feObjects(&feObjects)
    {
    };

    void assembly(const NumLib::TimeStep &/*time*/, const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &/*localDofManager*/, const MathLib::LocalVector &/*u1*/, const MathLib::LocalVector &/*u0*/,  MathLib::LocalMatrix &/*localJ*/)
    {
        FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);

        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
            q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);

//            double k;
//            _pm->hydraulic_conductivity->eval(real_x, k);
//            LocalMatrix matK(1,1);
//
//            //fe->integrateWxN(_pm->storage, localM);
//            fe->integrateDWxDN(j, k, localJ);
        }
    }
};

#endif

} //end
