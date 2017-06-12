/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file BisectionMethod.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#pragma once

#include "NRCheckConvergence.h"
#include "NRErrorNorm1Residual.h"

namespace MathLib
{

/**
 * \brief Bisection nonlinear solver
 */
class BisectionMethod
{
public:
    /**
     * solve the given nonlinear problem
     *
     * \tparam F_PROBLEM        Function class evaluating residuals
     * \tparam T_VALUE          Value type
     * \tparam T_CONVERGENCE    Convergence check class
     *
     * \param f_residuals
     * \param x_a
     * \param x_c
     * \param x_b
     * \param r_a
     * \param r_b
     * \param r_c
     * \param max_itr_count
     * \param convergence
     */
    template<class F_PROBLEM, class T_VALUE, class T_CONVERGENCE>
    int solve(F_PROBLEM &f_residuals,  T_VALUE &x_a, T_VALUE &x_c, T_VALUE &x_b, T_VALUE &r_a, T_VALUE &r_b, T_VALUE &r_c, size_t max_itr_count=100, T_CONVERGENCE* convergence=NULL)
    {
        T_CONVERGENCE _default_convergence;
        if (convergence==0) convergence = &_default_convergence;

        f_residuals.eval(x_a, r_a);
        f_residuals.eval(x_c, r_c);

        bool converged = false;
        std::cout << "Bisection iteration started!" << std::endl;
        for (size_t i=0; i<max_itr_count; i++) {
            x_b = x_a+x_c;
            x_b *= 0.5;
            f_residuals.eval(x_b, r_b);
            printout(i, x_a, x_b, x_c, r_a, r_b, r_c);
            //TODO check singularity. abs(r_c-r_a) should be smaller than the previous one.
            if (convergence->check(&r_b, 0, &x_b)) {
                converged = true;
                break;
            }
            if (r_a*r_b<.0) {
                x_c = x_b;
                r_c = r_b;
            } else if (r_b*r_c<.0) {
                x_a = x_b;
                r_a = r_b;
            } else {
                std::cout << "->*** Warning: neither f(a)*f(b) or f(b)*f(c) is negative." << std::endl;
                break;
            }
        }

        if (converged) return 0;
        std::cout << "->*** Warning: the iterations didn't converge." << std::endl;

        return -1; //not converged
    }

    /// solve scalar problems
    template<class F_PROBLEM>
    int solve(F_PROBLEM &fun, double &x_a, double &x_c, double &x_new, double error=1e-6, size_t max_itr_count=100)
    {
        double u[3];
        NRCheckConvergence<double,NRErrorNorm1Residual> check(error);
        return solve(fun, x_a, x_c, x_new, u[0], u[1], u[2], max_itr_count, &check);
    }

private:
    template<class T_VALUE>
    inline void printout(size_t i,  T_VALUE &x_a, T_VALUE &x_b, T_VALUE &x_c, T_VALUE &u_a, T_VALUE &u_b, T_VALUE &u_c);
//    {
//        std::cout << "-> " << i <<": x=(";
//        for (size_t i=0; i<x_new.size(); i++) std::cout << x_new[i] << " ";
//        std::cout << "), dx=(";
//        for (size_t i=0; i<dx.size(); i++) std::cout << dx[i] << " ";
//        std::cout << ")" << std::endl;
//    }
};

template<>
inline void BisectionMethod::printout(size_t i,  double &x_a, double &x_b, double &x_c, double &u_a, double &u_b, double &u_c)
{
    std::cout << "-> " << i <<": a=" << x_a << ", b=" << x_b << ", c=" << x_c << ", f(a)=" << u_a << ", f(b)=" << u_b << ", f(c)=" << u_c << std::endl;
}

} //end
