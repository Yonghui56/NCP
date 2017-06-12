/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FdmGw1D.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "NumLib/Function/DiscreteDataConvergenceCheck.h"
#include "Tests/Geo/Material/PorousMedia.h"
#include "IStencilWiseTransientLinearEQSLocalAssembler.h"
#include "FdmIVBVProblem.h"
#include "SingleStepFDM.h"
#include "FdmFunction.h"

class FdmGw1DLocalAssembler: public FdmLib::IStencilWiseTransientLinearEQSLocalAssembler
{
private:
    Geo::PorousMedia* _pm;
    double _h;
    double _h2;
public:
    FdmGw1DLocalAssembler(double h, Geo::PorousMedia &pm)
    : _pm(&pm), _h(h), _h2(h*h)
    {
    };

    virtual void assembly(const NumLib::TimeStep &/*time*/,  FdmLib::IStencil &s, const LocalVectorType &/*local_u_n1*/, const LocalVectorType &/*local_u_n*/, LocalEquationType &eqs)
    {
        //const double dt = time.getTimeStepSize();
        const size_t center_point_id = s.getCentralNodeID();
        //double storage = .0;
        //_pm->storage->eval(0, storage);
        double k = .0;
        _pm->hydraulic_conductivity->eval(0, k);

        if (center_point_id==0) { //Neuman BC
            eqs.addA(0, 0, k/_h);
            const std::vector<size_t> &neighbor_points = s.getSurroundingNodes();
            for(size_t i=0; i<neighbor_points.size(); i++)
            {
                eqs.addA(0, i+1, -k/_h);
            }
//        } else if (center_point_id == 4) {
        } else if (center_point_id == 20) { // Dirichlet BC
            eqs.addA(0, 0, 1.0);
        } else {
            eqs.addA(0, 0, 2.0*k/_h2);
            const std::vector<size_t> &neighbor_points = s.getSurroundingNodes();

            for(size_t i=0; i<neighbor_points.size(); i++)
            {
              eqs.addA(0, i+1, -k/_h2);
            }
        }

//        std::cout << center_point_id << ":"; eqs.printout();

    }
};

typedef FdmLib::FdmIVBVProblem
        <
        FdmGw1DLocalAssembler
        > GWFdmProblem;



template <
    class T_LINEAR_SOLVER
    >
class FunctionHead : public NumLib::AbstractTransientMonolithicSystem
{
    enum Out { Head=0 };
public:
    typedef FdmLib::SingleStepFDM
            <
                GWFdmProblem,
                T_LINEAR_SOLVER
            > SolutionForHead;

    FunctionHead()
    : _problem(NULL), _solHead(NULL), _rec(NULL)
    {
        AbstractTransientMonolithicSystem::resizeOutputParameter(1);
    };

    void define(DiscreteLib::DiscreteSystem* dis, GWFdmProblem* problem, BaseLib::Options &option)
    {
        //solution algorithm
        _solHead = new SolutionForHead(dis, problem);
        //_solHead->getTimeODEAssembler()->setTheta(1.0);
        typename SolutionForHead::LinearSolverType* linear_solver = _solHead->getLinearEquationSolver();
        linear_solver->setOption(option);
        this->setOutput(Head, problem->getIC(0));
    }

    int solveTimeStep(const NumLib::TimeStep &time)
    {
        _solHead->solveTimeStep(time);
        setOutput(Head, _solHead->getCurrentSolution(0));
        return 0;
    }

    double suggestNext(const NumLib::TimeStep &time_current) { return _solHead->suggestNext(time_current); }

    bool isAwake(const NumLib::TimeStep &time) { return _solHead->isAwake(time);  }

    void accept(const NumLib::TimeStep &time)
    {
        _solHead->accept(time);

        //std::cout << "Head=" << std::endl;
        //_solHead->getCurrentSolution(0)->printout();
    };

    NumLib::DiscreteDataConvergenceCheck _checker;
    virtual NumLib::IConvergenceCheck* getConvergenceChecker()
    {
        return &_checker;
    }

private:
    GWFdmProblem* _problem;
    SolutionForHead* _solHead;
    GeoLib::Rectangle *_rec;

    DISALLOW_COPY_AND_ASSIGN(FunctionHead);
};

class FunctionFdmVelocity
    : public NumLib::AbstractTransientMonolithicSystem
{
    enum In { Head=0 };
    enum Out { Velocity=0 };
public:

    FunctionFdmVelocity()
    : _dis(NULL), _vel(NULL), _K(NULL)
    {
        AbstractTransientMonolithicSystem::resizeInputParameter(1);
        AbstractTransientMonolithicSystem::resizeOutputParameter(1);
    };

    void define(DiscreteLib::DiscreteSystem &dis, Geo::PorousMedia &pm)
    {
        _dis = &dis;
        _K = pm.hydraulic_conductivity;
        _vel = new FdmLib::FdmCellVectorFunction(dis.getMesh()->getNumberOfElements());
        //this->setOutput(Velocity, _vel);
    }

    int solveTimeStep(const NumLib::TimeStep &/*time*/)
    {
        const MeshLib::IMesh *msh = _dis->getMesh();
        FdmLib::FdmFunctionScalar *head = (FdmLib::FdmFunctionScalar*)getInput(Head);

        //calculate vel (vel=f(h))
        for (size_t i_e=0; i_e<msh->getNumberOfElements(); i_e++) {
            MeshLib::IElement* e = msh->getElement(i_e);
            std::vector<double> local_h(e->getNumberOfNodes());
            for (size_t j=0; j<e->getNumberOfNodes(); j++)
                local_h[j] = head->getValue(e->getNodeID(j));
            // for each integration points
            std::vector<double> xi(e->getNumberOfNodes());
            std::vector<double> yi(e->getNumberOfNodes());
            for (size_t i=0; i<e->getNumberOfNodes(); i++) {
                const GeoLib::Point* pt = msh->getNodeCoordinatesRef(e->getNodeID(i));
                xi[i] = (*pt)[0];
                yi[i] = (*pt)[1];
            }
            MathLib::LocalVector q(2);
            q[0] = .0;
            q[1] = .0;
            double k;
            _K->eval(&xi[1], k);
            q[0] = -k* (local_h[1] - local_h[0]) / (xi[1]-xi[0]);
            _vel->setValue(i_e, q);
        }

        setOutput(Velocity, _vel);
        return 0;
    }

    double suggestNext(const NumLib::TimeStep &/*time_current*/)
    {
        return .0;
    }

    bool isAwake(const NumLib::TimeStep &/*time*/)
    {
        return true;
    }

    void accept(const NumLib::TimeStep &/*time*/)
    {
        //std::cout << "Velocity=" << std::endl;
        //_vel->printout();
    };

    NumLib::DiscreteDataConvergenceCheck _checker;
    virtual NumLib::IConvergenceCheck* getConvergenceChecker()
    {
        return &_checker;
    }

private:
    DiscreteLib::DiscreteSystem* _dis;
    FdmLib::FdmCellVectorFunction* _vel;
    NumLib::ITXFunction *_K;

    DISALLOW_COPY_AND_ASSIGN(FunctionFdmVelocity);
};
