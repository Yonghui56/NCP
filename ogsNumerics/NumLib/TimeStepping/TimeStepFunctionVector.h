/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TimeStepFunctionVector.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "TimeUnit.h"
#include "AbstractTimeStepFunction.h"

namespace NumLib
{


class TimeStepFunctionVector : public AbstractTimeStepFunction
{
public:
    TimeStepFunctionVector(double t0, double tn, const std::vector<double> &time_step_vector, TimeUnit::type t_unit=TimeUnit::Second)
    : AbstractTimeStepFunction(t0, tn, t_unit)
    {
        _dt_vector.resize(time_step_vector.size());
        double t_sum = 0;
        for (size_t i=0; i<time_step_vector.size(); i++) {
            _dt_vector[i] = TimeUnit::convertToSec(t_unit, time_step_vector[i]);
            t_sum += _dt_vector[i];
        }

        if (t_sum < AbstractTimeStepFunction::getEnd())
            AbstractTimeStepFunction::setEnd(t_sum);
    };

    TimeStepFunctionVector* clone()
    {
        TimeStepFunctionVector *obj = new TimeStepFunctionVector(getBeginning(), getEnd(), _dt_vector);
        return obj;
    }

protected:
    double suggestNext(double)
    {
        double t = getPrevious() + _dt_vector[AbstractTimeStepFunction::getStep()];
        if (t > getEnd())
            return getEnd();
        else
            return t;
    };

private:
    std::vector<double> _dt_vector;
};


}
