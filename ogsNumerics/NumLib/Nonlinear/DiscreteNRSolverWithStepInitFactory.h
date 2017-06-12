/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DiscreteNRSolverWithStepInitFactory.h
 *
 * Created on 2012-10-04 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/Options.h"
#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "DiscreteLib/Core/IDiscreteVector.h"
#include "NonlinearSolver.h"
#include "NonlinearSolverOption.h"

namespace NumLib
{

/**
 * \brief Factory class generating a discrete nonlinear solver
 */
template <class T_NR_STEP_INIT>
class DiscreteNRSolverWithStepInitFactory
{
public:
    DiscreteNRSolverWithStepInitFactory(T_NR_STEP_INIT *step_init) 
        : _step_init(step_init)
    {
    }

    template <class T_DIS_SYS, class F_LINEAR, class F_R, class F_DX>
    INonlinearSolver* create(const NonlinerSolverOption &nl_option, T_DIS_SYS* dis_sys,F_LINEAR* f_l, F_R* f_r, F_DX* f_dx)
    {
        INonlinearSolver* solver = NULL;
        switch (nl_option.solver_type)
        {
        case NonlinerSolverOption::LINEAR:
            solver = new Linear<F_LINEAR>(f_l);
            break;
        case NonlinerSolverOption::PICARD:
            solver = new Picard<T_DIS_SYS, F_LINEAR>(dis_sys, f_l);
            break;
        case NonlinerSolverOption::NEWTON:
            solver = new NewtonRaphson<T_DIS_SYS, F_R, F_DX, T_NR_STEP_INIT>(dis_sys, f_r, f_dx, _step_init);
            break;
		case NonlinerSolverOption::NEWTON_LINESEARCH:
			solver = new NewtonLineSearch<T_DIS_SYS, F_R, F_DX, T_NR_STEP_INIT>(dis_sys, f_r, f_dx, _step_init);
			break;
        default:
            break;
        }
        solver->setOption(nl_option);
        return solver;
    }

private:
    T_NR_STEP_INIT* _step_init;
};

} //
