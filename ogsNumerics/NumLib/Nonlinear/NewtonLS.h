/**
* Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.com)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.com/LICENSE.txt
*
*
* \file NewtonLS.h
*
* Created on 2015-03-31 by Haibing Shao
*/

#pragma once

#include "BaseLib/CodingTools.h"
#include "MathLib/Nonlinear/NewtonLineSearch.h"
#include "MathLib/Nonlinear/NRIterationStepInitializerDummy.h"
#include "DiscreteLib/Core/IDiscreteSystem.h"

#include "INonlinearSolver.h"

namespace NumLib
{

    /**
    * \brief Newton Line Search solver for discrete systems
    */
    template <
        class T_DIS_SYS,
        class F_R,
        class F_DX,
        class T_STEP_INIT = MathLib::NRIterationStepInitializerDummy
    >
    class NewtonLineSearch : public INonlinearSolver
    {
    public:
        typedef T_DIS_SYS MyDiscreteSystem;

        NewtonLineSearch(MyDiscreteSystem* dis_sys, F_R* f_r, F_DX* f_dx)
            : _dis_sys(dis_sys), _f_r(f_r), _f_dx(f_dx), _nl_step_init(NULL)
        {
            _r = _dx = 0;
            _x_tmp = 0;
        };

        NewtonLineSearch(MyDiscreteSystem* dis_sys, F_R* f_r, F_DX* f_dx, T_STEP_INIT* nl_step_init)
            : _dis_sys(dis_sys), _f_r(f_r), _f_dx(f_dx), _nl_step_init(nl_step_init)
        {
            _r = _dx = 0;
            _x_tmp = 0; 
        };

        virtual ~NewtonLineSearch()
        {
            _dis_sys->deleteVector(_r);
            _dis_sys->deleteVector(_dx);
            _dis_sys->deleteVector(_x_tmp);
            BaseLib::releaseObject(_nl_step_init);
        };

        virtual void solve(const VectorType &x_0, VectorType &x_new)
        {
            if (_r == 0) {
                _r = _dis_sys->template createVector<double>(x_0.size());
                _dx = _dis_sys->template createVector<double>(x_0.size());
                _x_tmp = _dis_sys->template createVector<double>(x_0.size());
            }
            MathLib::NRCheckConvergence<VectorType, MathLib::NRErrorAbsResMNormOrRelDxMNorm> check(_option.error_tolerance);
            _nr.solve(*_f_r, *_f_dx, x_0, *_x_tmp, x_new, *_r, *_dx, _option.max_iteration, &check, _nl_step_init);
        }

        virtual void recordLog(BaseLib::Options& opt)
        {
            BaseLib::OptionGroup* opNR = opt.getSubGroup("NewtonLineSearch");
            if (opNR == nullptr) {
                opNR = opt.addSubGroup("NewtonLineSearch");
            }
            opNR->setOptionAsNum<size_t>("iterations", _nr.getNIterations());
            opNR->setOptionAsNum<double>("error", _nr.getError());
        }

    private:
        DISALLOW_COPY_AND_ASSIGN(NewtonLineSearch);

    private:
        MathLib::NewtonLineSearchMethod _nr;
        MyDiscreteSystem* _dis_sys;
        F_R* _f_r;
        F_DX* _f_dx;
        VectorType* _r;
        VectorType* _dx;
        VectorType* _x_tmp; 
        T_STEP_INIT* _nl_step_init;
    };



} //end