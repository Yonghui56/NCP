/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementVelocity.h
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "FemLib/Function/FemFunction.h"
#include "NumLib/Function/Function.h"
#include "NumLib/Function/DiscreteDataConvergenceCheck.h"
#include "NumLib/Function/TXWrapped3DVectorFunction.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "ProcessLib/AbstractTimeIndependentProcess.h"

template <class T_DISCRETE_SYSTEM>
class FunctionHeadToElementVelocity
    : public ProcessLib::Process
{
public:
    enum In { Head=0 };
    enum Out { Velocity=0 };

    typedef T_DISCRETE_SYSTEM MyDiscreteSystem;
    typedef typename FemLib::FemNodalFunctionScalar<MyDiscreteSystem>::type MyNodalFunctionScalar;
    typedef typename FemLib::FEMIntegrationPointFunctionVector<MyDiscreteSystem>::type MyIntegrationPointFunctionVector;
    typedef typename NumLib::TXWrapped3DVectorFunction<MyIntegrationPointFunctionVector> My3DIntegrationPointFunctionVector;

    FunctionHeadToElementVelocity() 
    : ProcessLib::Process("HEAD_TO_ELEMENT_VELOCITY", 1, 1),
      _dis(NULL), _vel(NULL), _feObjects(NULL), _vel_3d(NULL), _tim(NULL)
    {
        // set default parameter name
        ProcessLib::Process::setInputParameterName(Head, "Head");
        ProcessLib::Process::setOutputParameterName(Velocity, "Velocity");
    };

    virtual ~FunctionHeadToElementVelocity()
    {
        BaseLib::releaseObject(_feObjects, _vel, _vel_3d);
    };


    virtual bool initialize(const BaseLib::Options &op);

    virtual void finalize() {};

    double suggestNext(const NumLib::TimeStep &time_current) { return _tim->getNext(time_current.getTime()); }

    bool isAwake(const NumLib::TimeStep &time) { return time.getTime()==suggestNext(time.getPreviousTime());  }

    int solveTimeStep(const NumLib::TimeStep &/*time*/);

    virtual bool accept(const NumLib::TimeStep &time)
    {
        return _tim->accept(time.getTime());
    }

    virtual void finalizeTimeStep(const NumLib::TimeStep &/*time*/);

    ///
    virtual NumLib::IConvergenceCheck* getConvergenceChecker() { return &_checker; };


private:
    MyDiscreteSystem* _dis;
    MyIntegrationPointFunctionVector* _vel;
    FemLib::LagrangeFeObjectContainer* _feObjects;
    NumLib::DiscreteDataConvergenceCheck _checker;
    My3DIntegrationPointFunctionVector* _vel_3d;
    NumLib::ITimeStepFunction* _tim;

    DISALLOW_COPY_AND_ASSIGN(FunctionHeadToElementVelocity);
};

#include "HeadToElementVelocity.hpp"



