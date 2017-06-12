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
#include "MassTransportTimeODELocalAssembler.h"
#include "MassTransportJacobianLocalAssembler.h"

/**
 *
 */
template <class T_DISCRETE_SYSTEM, class T_LINEAR_SOLVER>
class FunctionConcentration
: public ProcessLib::AbstractTransientProcess
{
public:
    enum In { Velocity=0 };
    enum Out { Concentration = 0 };

    typedef T_DISCRETE_SYSTEM MyDiscreteSystem;
    typedef T_LINEAR_SOLVER MyLinearSolver;
    // local assembler
    typedef MassTransportTimeODELocalAssembler<NumLib::ElementWiseTimeEulerEQSLocalAssembler> MyLinearAssemblerType;
    typedef MassTransportTimeODELocalAssembler<NumLib::ElementWiseTimeEulerResidualLocalAssembler> MyResidualAssemblerType;
    typedef MassTransportJacobianLocalAssembler MyJacobianAssemblerType;
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

    FunctionConcentration() 
        : AbstractTransientProcess("MASS_TRANSPORT", 1, 1),
          _problem(0), _solution(0), _feObjects(0), _compound(0)
    {
        // set default parameter name
        ProcessLib::AbstractTransientProcess::setInputParameterName(Velocity, "Velocity");
        ProcessLib::AbstractTransientProcess::setOutputParameterName(Concentration, "Concentration");
    };

    virtual ~FunctionConcentration()
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
    DISALLOW_COPY_AND_ASSIGN(FunctionConcentration);

private:
    MyProblemType* _problem;
    MySolutionType* _solution;
    FemLib::LagrangeFeObjectContainer* _feObjects;
    MaterialLib::Compound* _compound;
    NumLib::DiscreteDataConvergenceCheck _checker;
};

#include "Concentration.hpp"

