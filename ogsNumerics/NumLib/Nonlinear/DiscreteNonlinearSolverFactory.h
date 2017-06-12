/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DiscreteNonlinearSolverFactory.h
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
class DiscreteNonlinearSolverFactory
{
public:
    template <class T_DIS_SYS, class F_LINEAR, class F_R, class F_DX>
    INonlinearSolver* create(const NonlinerSolverOption &nl_option, T_DIS_SYS* dis_sys,F_LINEAR* f_l, F_R* f_r, F_DX* f_dx)
    {
        INonlinearSolver* solver = 0;
        switch (nl_option.solver_type)
        {
        case NonlinerSolverOption::LINEAR:
            solver = new Linear<F_LINEAR>(f_l);
            break;
        case NonlinerSolverOption::PICARD:
            solver = new Picard<T_DIS_SYS, F_LINEAR>(dis_sys, f_l);
            break;
        case NonlinerSolverOption::NEWTON:
            solver = new NewtonRaphson<T_DIS_SYS, F_R, F_DX>(dis_sys, f_r, f_dx);
            break;
		case NonlinerSolverOption::NEWTON_LINESEARCH:
			solver = new NewtonLineSearch<T_DIS_SYS, F_R, F_DX>(dis_sys, f_r, f_dx);
			break;
        default:
            break;
        }
        solver->setOption(nl_option);
        return solver;
    }
};

} //
