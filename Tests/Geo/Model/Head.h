/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Head.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MathLib/LinAlg/VectorNorms.h"
#include "GeoLib/Rectangle.h"
#include "DiscreteLib/Serial/DiscreteSystem.h"
#include "FemLib/Function/FemFunction.h"
#include "NumLib/Function/Function.h"
#include "NumLib/Function/DiscreteDataConvergenceCheck.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerResidualLocalAssembler.h"
#include "SolutionLib/Fem/FemIVBVProblem.h"
#include "SolutionLib/Fem/SingleStepFEM.h"

#include "Tests/Geo/Equation/FemGroundwaterFlow.h"
#include "Tests/Geo/Material/PorousMedia.h"

using namespace GeoLib;
using namespace MathLib;
using namespace FemLib;
using namespace NumLib;
using namespace MeshLib;
using namespace SolutionLib;
using namespace DiscreteLib;

namespace Geo
{




template <
    class T_LINEAR_SOLVER
    >
class FunctionHead : public AbstractTransientMonolithicSystem
{
    enum Out { Head=0 };
public:

    typedef TemplateFemEquation<
            DiscreteSystem,
            T_LINEAR_SOLVER,
            Geo::GroundwaterFlowTimeODELocalAssembler<ElementWiseTimeEulerEQSLocalAssembler>,
            Geo::GroundwaterFlowTimeODELocalAssembler<ElementWiseTimeEulerResidualLocalAssembler>,
            Geo::GroundwaterFlowJacobianLocalAssembler
            >
            GWFemEquation;

    typedef FemIVBVProblem< DiscreteSystem, GWFemEquation > GWFemProblem;

    typedef SingleStepFEM
            <
                GWFemProblem,
                T_LINEAR_SOLVER
            > SolutionForHead;

    FunctionHead() 
    : _problem(NULL), _solHead(NULL), _rec(NULL), _feObjects(NULL)
    {
        AbstractTransientMonolithicSystem::resizeOutputParameter(1);
    };

    void define(DiscreteSystem* dis, GWFemProblem* problem, BaseLib::Options &option)
    {
        //MeshLib::IMesh *msh = dis.getMesh();
        //size_t nnodes = msh->getNumberOfNodes();
        //solution algorithm
        _solHead = new SolutionForHead(dis, problem);
        //_solHead->getTimeODEAssembler()->setTheta(1.0);
        typename SolutionForHead::LinearSolverType* linear_solver = _solHead->getLinearEquationSolver();
        linear_solver->setOption(option);
        _solHead->getNonlinearSolver()->setOption(option);
        this->setOutput(Head, _solHead->getCurrentSolution(0));
    }

    int solveTimeStep(const TimeStep &time)
    {
        _solHead->solveTimeStep(time);
        setOutput(Head, _solHead->getCurrentSolution(0));
        return 0;
    }

    double suggestNext(const TimeStep &time_current) { return _solHead->suggestNext(time_current); }

    bool isAwake(const TimeStep &time) { return _solHead->isAwake(time);  }

    void accept(const TimeStep &time) { _solHead->accept(time); };



    //std::cout << "Head=" << std::endl;
    //_solHead->getCurrentSolution(0)->printout();
    NumLib::DiscreteDataConvergenceCheck _checker;
    virtual NumLib::IConvergenceCheck* getConvergenceChecker()
    {
        return &_checker;
    }

private:
    GWFemProblem* _problem;
    SolutionForHead* _solHead;
    GeoLib::Rectangle *_rec;
    //FemNodalFunctionScalar *_head;
    LagrangeFeObjectContainer* _feObjects;

    DISALLOW_COPY_AND_ASSIGN(FunctionHead);
};


} //end


