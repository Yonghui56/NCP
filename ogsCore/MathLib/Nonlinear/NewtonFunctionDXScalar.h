/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NewtonFunctionDXScalar.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#pragma once

namespace MathLib
{

/**
 * \brief Scalar function evaluating solution increment for Newton method 
 *
 * \tparam F_JACOBIAN   A function evaluating jacobian value
 */
template<class F_JACOBIAN>
class NewtonFunctionDXScalar
{
public:
    explicit NewtonFunctionDXScalar(F_JACOBIAN &f) : _f_j(&f) {};

    void eval(const double &x, const double &r, double &dx)
    {
        double j;
        _f_j->eval(x, j);
        if (j==0) return; //TODO error
        dx = -r / j;
    }

private:
    F_JACOBIAN* _f_j;
};


} //end

