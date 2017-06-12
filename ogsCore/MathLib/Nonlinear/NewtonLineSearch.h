/**
* Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.com)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.com/LICENSE.txt
*
*
* \file NewtonRaphson.h
*
* Created on 2015-03-31 by Haibing Shao
*/

#pragma once

#include "logog.hpp"

#include <iostream>
#include <fstream>
#include "MathLib/LinAlg/LinearEquation/DenseLinearEquation.h"
#include "NRCheckConvergence.h"
#include "NRErrorAbsResMNormOrRelDxMNorm.h"
#include "NRErrorAbsResMNorm.h"
#include "NewtonFunctionDXVector.h"
#include "NewtonFunctionDXScalar.h"
#include "NRIterationStepInitializerDummy.h"

namespace MathLib
{

    /**
    * \brief Newton Line Search method
    */
    class NewtonLineSearchMethod
    {
    public:
        NewtonLineSearchMethod() : _log_itr_count(0), _log_error(.0) {};

        /// general solver
        /// @tparam F_RESIDUALS
        /// @tparam F_DX
        /// @tparam T_VALUE
        /// @tparam T_CONVERGENCE
        template<class F_RESIDUALS, class F_DX, class T_VALUE, class T_CONVERGENCE, class T_PREPOST>
        int solve(F_RESIDUALS &f_residuals, F_DX &f_dx, const T_VALUE &x0, T_VALUE &x_pre, T_VALUE &x_new, T_VALUE &r, T_VALUE &dx, size_t max_itr_count = 100, T_CONVERGENCE* convergence = NULL, T_PREPOST* pre_post = NULL)
        {
            std::size_t j(0); 
            double d_norm(9999999.9), d1_norm(9999999.9);
            const double N_L = 1000 / 0.018;
            const double H = 7.65e-6;
            const double damping_factor(1.0); 
            std::size_t n_nodes = x0.size() / 3;
            T_CONVERGENCE _default_convergence;
            if (convergence == NULL) convergence = &_default_convergence;
            r = .0;
            x_new = x0;
            x_pre = x0;
            INFO("Newton Line Search iteration started!");
            // HS: verifying my theory..............................
            if (pre_post)
                pre_post->pre_process(dx, x_pre, f_residuals, f_dx);
            // .....................................................

            bool converged = false;
            size_t itr_cnt = 0;
            f_residuals.eval(x_pre, r);
            converged = convergence->check(&r, &dx, &x_pre);
            d_norm = convergence->getError();
            if (!converged) {
                for (itr_cnt = 0; itr_cnt<max_itr_count; itr_cnt++) {
                    // preprocessing
                    if (pre_post)
                        pre_post->pre_process(dx, x_pre, f_residuals, f_dx);
                    // Jacobian
                    f_dx.eval(x_pre, r, dx);
                    // x increment
                    x_new = x_pre; 
					dx *= damping_factor;
                    x_new -= dx;

                    std::size_t m_flag = 1;
                    ///modified newton regarding phase appearance
                    ///detect the phase change zone
                    ///update new x_pre
                    /// based on Rachford-Rich
                    ///after the first Newton Iteration

                    for (std::size_t n = 0; n < n_nodes; n++) {
                        double l_pre = (1 - x_pre[3 * n + 1]) / (1 - x_pre[3 * n] * H / N_L);
                        double l_new = (1 - x_new[3 * n + 1]) / (1 - x_new[3 * n] * H / N_L);
                        double const N_G = x_pre[3 * n] / 8.314 / 295.15;
                        double test = (1 - l_new)*N_L / (N_G*l_new - N_L*l_new + N_L);
                        if (abs(x_pre[3 * n + 2])<1e-9 && l_new > 0 && l_new<=1 && l_pre>1) {
                            INFO("This point is going to switch to two-phase zone");
                            INFO("the phase molar fraction is %f", l_new);
                            x_pre[3 * n + 2] = (1 - l_new)*N_L / (N_G*l_new - N_L*l_new + N_L);
                            m_flag = 0;
                        }
                    }

                    if (m_flag == 0) {
                        if (pre_post)
                            pre_post->post_process(dx, x_pre, f_residuals, f_dx);
                        INFO("Newton-Raphson Line Search modification");
                        f_residuals.eval(x_pre, r);
                        itr_cnt -= 1;
                        continue;
                    }
                    // post processing
                    if (pre_post)
                        pre_post->post_process(dx, x_new, f_residuals, f_dx);
                    // update residual
                    f_residuals.eval(x_new, r);

                    // the line search operations
                    j = 0; 
                    while (j < 20)
                    {
						
                        // calculate the norm of r -> d1_norm
                        converged = convergence->check(&r, &dx, &x_new);
                        d1_norm = convergence->getError();

                        if (d1_norm < d_norm)
                        {
                            break; 
                        }
                        else
                        {
							INFO("Global Line Search begins!");
							x_new = x_pre; 
                            dx *= 0.5;
                            x_new  -= dx;//+= dx;//
                        }

                        // post processing
                        if (pre_post)
                            pre_post->post_process(dx, x_new, f_residuals, f_dx);
                        f_residuals.eval(x_new, r); 
                        j++;
                    }

                    d_norm = d1_norm; 
                    x_pre = x_new; 

                    // printout(itr_cnt, x_new, r, dx);
                   if (converged)
                        break;
                }
            }

            // store log
            this->_log_itr_count = itr_cnt;
            this->_log_error = convergence->getError();

            INFO("------------------------------------------------------------------");
            INFO("*** Newton Line Search nonlinear solver computation result");
            if (max_itr_count == 1) {
                INFO("status    : iteration not required");
            }
            else {
                INFO("status    : %s", (converged ? "converged" : "***ERROR - DID NOT CONVERGE!"));
            }
            INFO("iteration : %d/%d", itr_cnt, max_itr_count);
            INFO("residuals : %1.3e (tolerance=%1.3e)", convergence->getError(), convergence->getTolerance());
            INFO("------------------------------------------------------------------");

            if (converged) return 0;
            std::cout << "->*** Warning: the iterations didn't converge." << std::endl;

            return -1; //not converged
        }

        /// solve scalar problems
        template<class F_RESIDUALS, class F_JACOBIAN>
        int solve(F_RESIDUALS &f_residuals, F_JACOBIAN &f_jac, const double &x0, double &x_new, double error = 1e-6, size_t max_itr_count = 100)
        {
            NewtonFunctionDXScalar<F_JACOBIAN> f_dx(f_jac);
            double r, dx;
            typedef NRCheckConvergence<double, NRErrorAbsResMNormOrRelDxMNorm> MyConvergence;
            MyConvergence check(error);
            return solve<   F_RESIDUALS,
                NewtonFunctionDXScalar<F_JACOBIAN>,
                double,
                MyConvergence,
                NRIterationStepInitializerDummy
            >(f_residuals, f_dx, x0, x_new, r, dx, max_itr_count, &check);
        }

        /// solve vector problems using a direct linear solver
        template<class F_RESIDUALS, class F_JACOBIAN, class T_V, class T_CONVERGENCE>
        int solve(F_RESIDUALS &f_residuals, F_JACOBIAN &f_jac, const T_V &x0, T_V &x_new, T_CONVERGENCE* check_error = 0, size_t max_itr_count = 100)
        {
            const size_t n = x0.size();
            T_V r(n), dx(n);
            MathLib::DenseLinearEquation dense;
            dense.create(n);
            NewtonFunctionDXVector<F_JACOBIAN, MathLib::DenseLinearEquation> f_dx(f_jac, dense);
            return solve<   F_RESIDUALS,
                NewtonFunctionDXVector<F_JACOBIAN, MathLib::DenseLinearEquation>,
                T_V,
                T_CONVERGENCE,
                NRIterationStepInitializerDummy
            >(f_residuals, f_dx, x0, x_new, r, dx, max_itr_count, check_error);
        }

        /// solve vector problems using a direct linear solver
        template<class F_RESIDUALS, class F_JACOBIAN, class T_V>
        int solve(F_RESIDUALS &f_residuals, F_JACOBIAN &f_jac, const T_V &x0, T_V &x_new, double error = 1e-6, size_t max_itr_count = 100)
        {
            NRCheckConvergence<T_V, NRErrorAbsResMNormOrRelDxMNorm> check(error);
            return solve(f_residuals, f_jac, x0, x_new, &check, max_itr_count);
        }

        /// get the number of iterations computed
        size_t getNIterations() const { return _log_itr_count; };

        /// get the final error
        double getError() const { return _log_error; };

    private:
        template<class T_VALUE>
        inline void printout(size_t i, T_VALUE &x_new, T_VALUE &r, T_VALUE &dx)
        {
            std::cout << "-> " << i << ": ";
            std::cout << "r=(";
            for (size_t i = 0; i<dx.size(); i++) std::cout << r[i] << " ";
            std::cout << std::endl;
#if 1
            std::cout << "-> " << i << ": x=(";
            for (size_t i = 0; i<x_new.size(); i++) std::cout << x_new[i] << " ";
            std::cout << "), r=(";
            for (size_t i = 0; i<dx.size(); i++) std::cout << r[i] << " ";
            std::cout << "), dx=(";
            for (size_t i = 0; i<dx.size(); i++) std::cout << dx[i] << " ";
            std::cout << ")" << std::endl;
#endif
        }

    private:
        size_t _log_itr_count;
        double _log_error;
    };

#if 1
    // template specialization
    template<>
    inline void NewtonLineSearchMethod::printout(size_t i, double &x_new, double &r, double &dx)
    {
        std::cout << "-> " << i << ": x=" << x_new << ", r=" << r << ", dx=" << dx << std::endl;
    }
#endif

} //end