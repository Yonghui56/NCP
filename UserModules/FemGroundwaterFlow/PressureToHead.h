/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FunctionPressureToHead.h
 *
 * Created on 2012-08-23 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "FemLib/Function/FemFunction.h"
#include "NumLib/Function/Function.h"
#include "NumLib/Function/DiscreteDataConvergenceCheck.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "ProcessLib/AbstractTimeIndependentProcess.h"

template <class T_DISCRETE_SYSTEM>
class FunctionPressureToHead
    : public ProcessLib::AbstractTimeIndependentProcess
{
public:
    enum In { Pressure=0 };
    enum Out { Head=0 };

    typedef T_DISCRETE_SYSTEM MyDiscreteSystem;
    typedef typename FemLib::FemNodalFunctionScalar<MyDiscreteSystem>::type MyNodalFunctionScalar;

    FunctionPressureToHead()
    : ProcessLib::AbstractTimeIndependentProcess("PRESSURE_TO_HEAD", 1, 1), _dis(NULL), _h(NULL)
    {
        // set default parameter name
        ProcessLib::AbstractTimeIndependentProcess::setInputParameterName(Pressure, "Pressure");
        ProcessLib::AbstractTimeIndependentProcess::setOutputParameterName(Head, "Head");
    };

    virtual ~FunctionPressureToHead()
    {
        BaseLib::releaseObject(_h);
    };


    virtual bool initialize(const BaseLib::Options &op);

    virtual void finalize() {};

    int solveTimeStep(const NumLib::TimeStep &/*time*/);

    virtual void finalizeTimeStep(const NumLib::TimeStep &/*time*/);

    ///
    virtual NumLib::IConvergenceCheck* getConvergenceChecker() { return &_checker; };

private:
    DiscreteLib::IDiscreteSystem* _dis;
    MyNodalFunctionScalar* _h;
    NumLib::DiscreteDataConvergenceCheck _checker;

    DISALLOW_COPY_AND_ASSIGN(FunctionPressureToHead);
};

#include "PressureToHead.hpp"



