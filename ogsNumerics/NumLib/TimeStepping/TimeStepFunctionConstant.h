/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TimeStepFunctionConstant.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "TimeUnit.h"
#include "AbstractTimeStepFunction.h"

namespace NumLib
{


class TimeStepFunctionConstant : public AbstractTimeStepFunction
{
private:
    double _dt;
public:
    TimeStepFunctionConstant(double t0, double tn, double dt, TimeUnit::type t_unit=TimeUnit::Second) : AbstractTimeStepFunction(t0, tn, t_unit)
    {
        _dt = TimeUnit::convertToSec(t_unit, dt);
    };
    ITimeStepFunction* clone()
    {
        TimeStepFunctionConstant *obj = new TimeStepFunctionConstant(getBeginning(), getEnd(), _dt);
        return obj;
    }

protected:
    double suggestNext(double) 
    {
        double t = getPrevious() + _dt;
        if (t > getEnd())
            return getEnd();
        else
            return t;
    };
};


}
