/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Concentration.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "NumLib/Function/Function.h"
#include "NumLib/Function/DiscreteDataConvergenceCheck.h"
#include "MathLib/LinAlg/VectorNorms.h"
#include "GeoLib/Rectangle.h"
#include "DiscreteLib/Serial/DiscreteSystem.h"
#include "FemLib/Function/FemFunction.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerResidualLocalAssembler.h"
#include "SolutionLib/Fem/FemIVBVProblem.h"
#include "SolutionLib/Fem/SingleStepFEM.h"

#include "Tests/Geo/Equation/FemMassTransport.h"
#include "Tests/Geo/Material/PorousMedia.h"
#include "Tests/Geo/Material/Compound.h"

namespace Geo
{



template <
    class T_LINEAR_SOLVER
    >
class FunctionConcentration : public NumLib::AbstractTransientMonolithicSystem
{
    enum In { Velocity=0 };
    enum Out { Concentration = 0 };

public:
    typedef TemplateFemEquation<
            DiscreteLib::DiscreteSystem,
            T_LINEAR_SOLVER,
            Geo::MassTransportTimeODELocalAssembler<NumLib::ElementWiseTimeEulerEQSLocalAssembler>,
            Geo::MassTransportTimeODELocalAssembler<NumLib::ElementWiseTimeEulerResidualLocalAssembler>,
            Geo::MassTransportJacobianLocalAssembler
            >
            MassFemEquation;

    typedef FemIVBVProblem< DiscreteLib::DiscreteSystem, MassFemEquation > MassFemProblem;
    typedef SolutionLib::SingleStepFEM
            <
                MassFemProblem,
                T_LINEAR_SOLVER
            > SolutionForConc;

    FunctionConcentration() 
    {
        AbstractTransientMonolithicSystem::resizeInputParameter(1);
        AbstractTransientMonolithicSystem::resizeOutputParameter(1);
    };

    void define(DiscreteLib::DiscreteSystem* dis, MassFemProblem* problem, BaseLib::Options &option)
    {
        //solution algorithm
        _solConc = new SolutionForConc(dis, problem);
        //_solConc->getTimeODEAssembler()->setTheta(1.0);
        typename SolutionForConc::LinearSolverType* linear_solver = _solConc->getLinearEquationSolver();
        linear_solver->setOption(option);
        _solConc->getNonlinearSolver()->setOption(option);
        this->setOutput(Concentration, _solConc->getCurrentSolution(0));
    }

    int solveTimeStep(const NumLib::TimeStep &time)
    {
        //input
        const NumLib::ITXFunction *vel = this->getInput<NumLib::ITXFunction>(Velocity);
        _solConc->getProblem()->getEquation()->getLinearAssembler()->initialize(vel);
        _solConc->getProblem()->getEquation()->getResidualAssembler()->initialize(vel);
        _solConc->getProblem()->getEquation()->getJacobianAssembler()->initialize(vel);

        // solve
        _solConc->solveTimeStep(time);

        // output
        setOutput(Concentration, _solConc->getCurrentSolution(0));

        return 0;
    }

    double suggestNext(const NumLib::TimeStep &time_current) { return _solConc->suggestNext(time_current); }

    bool isAwake(const NumLib::TimeStep &time) { return _solConc->isAwake(time);  }

    void accept(const NumLib::TimeStep &time)
    {
        _solConc->accept(time);

        //std::cout << "Concentration=" << std::endl;
        //_solConc->getCurrentSolution(0)->printout();
    };

    NumLib::DiscreteDataConvergenceCheck _checker;
    virtual NumLib::IConvergenceCheck* getConvergenceChecker()
    {
        return &_checker;
    }

private:
    DISALLOW_COPY_AND_ASSIGN(FunctionConcentration);

    SolutionForConc* _solConc;
    //FemNodalFunctionScalar *_conc;
    FemLib::LagrangeFeObjectContainer* _feObjects;
};

} //end
