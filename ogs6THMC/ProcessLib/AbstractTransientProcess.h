/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file AbstractTransientProcess.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "logog.hpp"

#include "BaseLib/Options.h"
#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"
#include "SolutionLib/Core/AbstractTimeSteppingAlgorithm.h"
#include "Process.h"

namespace ProcessLib
{

/**
 * \brief Implementation of Process (ITransientSystem) class for monolithic system
 *
 */
class AbstractTransientProcess : public Process
{
public:
    ///
    AbstractTransientProcess(const std::string &pcs_name, size_t n_in_parameters, size_t n_out_parameters)
    : Process(pcs_name, n_in_parameters, n_out_parameters)
    {
    }

    ///
    virtual ~AbstractTransientProcess() {};

    ///
    virtual int solveTimeStep(const NumLib::TimeStep &time)
    {
        INFO("Solving %s...", getProcessName().c_str());
        initializeTimeStep(time);
        getSolution()->solveTimeStep(time);
        postSolutionAlgorithm(time);
        updateOutputParameter(time);
        return 0;
    }

    /// 
    virtual double suggestNext(const NumLib::TimeStep &time_current) 
    {
        return getSolution()->suggestNext(time_current); 
    }

    ///
    virtual bool isAwake(const NumLib::TimeStep &time) 
    { 
        return getSolution()->isAwake(time);  
    }

    virtual bool accept(const NumLib::TimeStep &time)
    {
        return getSolution()->accept(time);
    }

    ///
    virtual void finalizeTimeStep(const NumLib::TimeStep &time)
    {
        postTimeStep(time);
        output(time);
        getSolution()->finalizeTimeStep(time);
    };

protected:
    ///
    virtual SolutionLib::AbstractTimeSteppingAlgorithm* getSolution() = 0;
    ///
    virtual void initializeTimeStep(const NumLib::TimeStep &/*time*/) {};
    ///
    virtual void postSolutionAlgorithm(const NumLib::TimeStep &/*time*/) {};
    ///
    virtual void updateOutputParameter(const NumLib::TimeStep &time) = 0;
    ///
    virtual void postTimeStep(const NumLib::TimeStep &/*time*/) {};
    ///
    virtual void output(const NumLib::TimeStep &time) = 0;
};

} //end ProcessLib

