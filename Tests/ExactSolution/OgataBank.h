/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file OgataBank.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MathLib/MathTools.h"

inline double analyticalOgataBank(double x, double t, double v, double alpha)
{
    double u = MathLib::erfc((x-v*t)/std::sqrt(4*alpha*t))+std::exp(v*x/alpha)*MathLib::erfc((x+v*t)/std::sqrt(4*alpha*t));
    u *= 0.5;
    return u;
}

