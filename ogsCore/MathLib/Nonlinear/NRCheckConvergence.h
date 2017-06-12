/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NRCheckConvergence.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#pragma once

#include <cmath>
#include <algorithm>

namespace MathLib
{

/**
 * \brief Convergence check for Newton iterations
 *
 * \tparam T_D0     Value type
 * \tparam T_ERROR  Error type
 */
template<class T_D0, class T_ERROR>
class NRCheckConvergence
{
public:
    /// Default setting
    NRCheckConvergence() : _tolerance(1.e-6), _error(.0) {};

    /// 
    explicit NRCheckConvergence(double err) : _tolerance(err), _error(.0) {};

    /**
     * return if the solution convers or not
     *
     * \param r         residual
     * \param dx        solution increment
     * \param x_new     new solution
     */
    bool check(const T_D0* r, const T_D0* dx, const T_D0* x_new)
    {
        _error = calc.error(r, dx, x_new);
        return (fabs(_error) < _tolerance);
    }

    /// return calculated error
    double getError() const {return _error;};

    /// return current tolerance
    double getTolerance() const {return _tolerance;};

private:
    double _tolerance;
    double _error;
    T_ERROR calc;
};


}
