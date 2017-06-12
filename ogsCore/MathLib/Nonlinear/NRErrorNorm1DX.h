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

#include "MathLib/LinAlg/VectorNorms.h"


namespace MathLib
{

/**
 * \brief Check norm1 of dx
 */
class NRErrorNorm1DX
{
public:
    template<class T_D0>
    inline double error(const T_D0*, const T_D0* dx, const T_D0*) const
    {
        if (dx==0) return .0;
        return norm1(*dx, dx->size());
    }
};

template<>
inline double NRErrorNorm1DX::error(const double*, const double* dx, const double*) const
{
    return fabs(*dx);
}

}
