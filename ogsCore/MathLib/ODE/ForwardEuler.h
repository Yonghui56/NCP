/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ForwardEuler.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#pragma once

namespace MathLib
{

/**
 * \brief Forward Euler ODE solver
 */
class ForwardEuler
{
public:
    /**
     * solve dX/dt = F(t,X)
     *
     * \param f     Pointer to a function F(t,X)
     * \param t     Pointer to time step size array
     * \param n     the number of time steps
     * \param y     Pointer to time step value array
     */
    void solve(double (*f)(double, double), double *t, int n, double *y) const {
        for (int i=0; i<n-1; i++) {
            double dt = t[i+1] - t[i];
            y[i+1] = solve(f, y[i], t[i], dt);
        }
    };

    /**
     * solve dX/dt = F(t,X)
     *
     * \param f     Pointer to a function F(t,X)
     * \param y0    Previous time step value
     * \param t0    Previous time step
     * \param dt    Time step size
     */
    double solve(double (*f)(double, double), double y0, double t0, double dt) const {
        double y = y0 + dt*f(t0, y0);
        return y;
    };
};

}
