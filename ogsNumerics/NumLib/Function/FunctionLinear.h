/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FunctionLinear.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "IFunction.h"
#include "MathLib/Interpolation/LinearInterpolation.h"

namespace NumLib
{

class FunctionLinear1D : public TemplateFunction<double, double>
{
public:
    FunctionLinear1D(MathLib::LinearInterpolation* linear) {
        _linear = linear;
    };
    virtual void eval(const double& x, double &v)
    {
        v = _linear->getValue(x);
    };
    virtual void eval_slope(const double &x, double &slope)
    {
        slope = _linear->getSlope(x); 
    };
    virtual TemplateFunction<double,double>* clone() const
    {
        FunctionLinear1D* obj = new FunctionLinear1D(_linear);
        return obj;
    }
private:
    MathLib::LinearInterpolation *_linear;
};

} //end
