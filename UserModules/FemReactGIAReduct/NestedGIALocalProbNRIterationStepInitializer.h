/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NestedLocalProbNRIterationStepInitializer.h
 *
 * This class is used in Newton-Raphson iterations to include the 
 * calcualation of local ODE problems
 *
 * after the file DiscreteNonlinearSolverFactory.h
 * 
 * Created on 2012-10-05 by Haibing Shao
 */

#ifndef NESTED_GIA_Local_Prob_NR_ITERATION_STEP_INITIALIZER_H
#define NESTED_GIA_Local_Prob_NR_ITERATION_STEP_INITIALIZER_H


#include "NumLib/TransientAssembler/ElementWiseTimeEulerResidualLocalAssembler.h"

template <class T_USER_DATA>
class NestedGIALocalProbNRIterationStepInitializer
{
public:
    NestedGIALocalProbNRIterationStepInitializer(T_USER_DATA & user_data)
      : _user_data(user_data) {};

	template<class T_X, class F_RESIDUALS, class F_DX>
    void pre_process(const T_X &/*dx*/, const T_X &/*x_new*/, F_RESIDUALS & f_residuals, F_DX &/*f_dx*/)
    {
    };

    template<class T_X, class F_RESIDUALS, class F_DX>
    void post_process(const T_X &/*dx*/, const T_X & x_new, F_RESIDUALS &/*f_residuals*/, F_DX &/*f_dx*/)
    {
        // very important is, now we have x_new, 
        // the data in the FunctionConcentrations needs to be updated, 
        // so that the in the next Newton iteration, concentrations are calculated
        // based on the new x_new values
    	_user_data.update_xi_global_cur_nodal_values( x_new );

    	// call the local problem
    	_user_data.calc_nodal_local_problem(_user_data.getTimeStep()->getTimeStepSize(), 1.0E-8, 1.0E-20, 150);
    };

private:
    /**
      * pointer to user data
      */ 
    T_USER_DATA & _user_data;
};

#endif  // end of ifndef
