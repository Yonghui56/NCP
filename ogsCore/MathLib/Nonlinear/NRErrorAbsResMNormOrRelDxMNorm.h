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
 * \brief Check norm-max of residual and dx. Relative error is used for dx.
 */
class NRErrorAbsResMNormOrRelDxMNorm
{
public:
    NRErrorAbsResMNormOrRelDxMNorm() : _itr_count(0) {};

    template<class T_D0>
    inline double error(const T_D0* r, const T_D0* dx, const T_D0* x)
    {
        double abs_mnorm_r = norm_max(*r, r->size());
        double abs_mnorm_x = norm_max(*x, x->size());
        double abs_mnrom_dx = norm_max(*dx, dx->size());
        return calc(abs_mnorm_r, abs_mnorm_x, abs_mnrom_dx);
    }

private:
    inline double calc(double abs_mnorm_r, double abs_mnorm_x, double abs_mnrom_dx)
    {
        double rel_mnorm_dx = .0;
        if (abs_mnorm_x!=.0) rel_mnorm_dx = abs_mnrom_dx/abs_mnorm_x;

        INFO("-> %d: ||r||_inf=%1.3e, ||dx||_inf=%1.3e, ||x||_inf=%1.3e, ||dx||/||x||=%1.3e", _itr_count, abs_mnorm_r, abs_mnrom_dx, abs_mnorm_x, rel_mnorm_dx);

        _itr_count++;

        if (abs_mnorm_x == .0) {
            return abs_mnorm_r;
        } else {
            // return std::max(abs_mnorm_r, rel_mnorm_dx);
            // HS temperary testing...
            return abs_mnorm_r;
        }
    }

private:
    size_t _itr_count;
};

template<>
inline double NRErrorAbsResMNormOrRelDxMNorm::error(const double* r, const double* dx, const double* x)
{
    double abs_mnorm_r = fabs(*r);
    double abs_mnorm_x = fabs(*x);
    double abs_mnrom_dx = fabs(*dx);
    return calc(abs_mnorm_r, abs_mnorm_x, abs_mnrom_dx);
}

}
