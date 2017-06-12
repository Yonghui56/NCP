/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Convergence.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#pragma once

#include <cmath>
#include <algorithm>

#include "logog.hpp"

#include "MathLib/LinAlg/VectorNorms.h"


namespace MathLib
{

/**
 * \brief Check norm1 of residual
 */
class NRErrorNorm1Residual
{
public:
    template<class T_D0>
    inline double error(const T_D0* r, const T_D0*, const T_D0*) const
    {
        if (r==0) return .0;
        return norm1(*r, r->size());
    }
};

template<>
inline double NRErrorNorm1Residual::error(const double* r, const double*, const double*) const
{
    return fabs(*r);
}

}
