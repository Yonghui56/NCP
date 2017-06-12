/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Picard.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MathLib/Nonlinear/Picard.h"
#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "INonlinearSolver.h"

namespace NumLib
{

/**
 * \brief Picard
 */
template <class T_DIS_SYS, class F_LINEAR>
class Picard : public INonlinearSolver
{
public:
    typedef T_DIS_SYS MyDiscreteSystem;

    Picard(MyDiscreteSystem* dis_sys, F_LINEAR* linear_f) : _linear_f(linear_f), _dis_sys(dis_sys)
    {
        _x_old = _dx = 0;
    };

    virtual ~Picard()
    {
        _dis_sys->deleteVector(_x_old);
        _dis_sys->deleteVector(_dx);
    };

    virtual void solve(const VectorType &x_0, VectorType &x_new)
    {
        if (_x_old==0) {
            _x_old = _dis_sys->template createVector<double>(x_0.size());
            _dx = _dis_sys->template createVector<double>(x_0.size());
        }

        MathLib::PicardMethod picard;
        MathLib::NRCheckConvergence<VectorType, MathLib::NRErrorNorm1DX> check(getOption().error_tolerance);
        picard.solve(*_linear_f, x_0, x_new, *_x_old, *_dx, getOption().max_iteration, &check);
    }

    virtual void recordLog(BaseLib::Options& /*opt*/) {};

private:
    F_LINEAR* _linear_f;
    MyDiscreteSystem* _dis_sys;
    VectorType* _x_old;
    VectorType* _dx;
};

} //end

