/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Concentration.h
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerResidualLocalAssembler.h"
#include "NumLib/Function/DiscreteDataConvergenceCheck.h"
#include "SolutionLib/Fem/FemIVBVProblem.h"
#include "SolutionLib/Fem/TemplateFemEquation.h"
#include "SolutionLib/Fem/SingleStepFEM.h"
#include "ProcessLib/AbstractTransientProcess.h"

#include "MaterialLib/Compound.h"
#include "RichardsFlowTimeODELocalAssember.h"
#include "RichardsFlowJacobianLocalAssembler.h"



/**
 *
 */
template <class T_DISCRETE_SYSTEM, class T_LINEAR_SOLVER>
class FunctionRichards
: public ProcessLib::AbstractTransientProcess
{
public:
    enum In {};
    enum Out { Pressure = 0 };

    typedef T_DISCRETE_SYSTEM MyDiscreteSystem;
    typedef T_LINEAR_SOLVER MyLinearSolver;
    // local assembler // Todo this assembler must be changed later, after we setup some for Richards flow
    typedef RichardsFlowTimeODELocalAssembler<NumLib::ElementWiseTimeEulerEQSLocalAssembler> MyLinearAssemblerType;
    typedef RichardsFlowTimeODELocalAssembler<NumLib::ElementWiseTimeEulerResidualLocalAssembler> MyResidualAssemblerType;
    typedef RichardsFlowJacobianLocalAssembler MyJacobianAssemblerType;
    // Equation definition
    typedef SolutionLib::TemplateFemEquation<
            MyDiscreteSystem,
            MyLinearSolver,
            MyLinearAssemblerType,
            MyResidualAssemblerType,
            MyJacobianAssemblerType
            >
            MyEquationType;
    // IVBV problem definition
    typedef SolutionLib::FemIVBVProblem<
            MyDiscreteSystem,
            MyEquationType
            > MyProblemType;
    // Solution algorithm definition
    typedef SolutionLib::SingleStepFEM
            <
            MyProblemType,
            MyLinearSolver
            > MySolutionType;
	//Constructor
    FunctionRichards() 
        : AbstractTransientProcess("RICHARDS_FLOW", 0, 1),
          _problem(0), _solution(0), _feObjects(0)
    {
        // set default parameter name
        //ProcessLib::AbstractTransientProcess::setInputParameterName(Velocity, "Velocity");
        ProcessLib::AbstractTransientProcess::setOutputParameterName(Pressure, "Pressure");
    };

	//DESTRUCTOR
    virtual ~FunctionRichards()
    {
        BaseLib::releaseObject(_problem, _solution, _feObjects);
    };

    /// initialize this process
    virtual bool initialize(const BaseLib::Options &option);

    /// finalize but nothing to do here
    virtual void finalize() {};

    virtual NumLib::IConvergenceCheck* getConvergenceChecker()
    {
        return &_checker;
    }

protected:
    virtual void initializeTimeStep(const NumLib::TimeStep &time);

    virtual void updateOutputParameter(const NumLib::TimeStep &time);

    virtual MySolutionType* getSolution() {return _solution;};

    virtual void output(const NumLib::TimeStep &time);

private:
    DISALLOW_COPY_AND_ASSIGN(FunctionRichards);

private:
    MyProblemType* _problem;
    MySolutionType* _solution;
    FemLib::LagrangeFeObjectContainer* _feObjects;
    //MaterialLib::Compound* _compound;
    NumLib::DiscreteDataConvergenceCheck _checker;
};

#include "FunctionRichards.hpp"

