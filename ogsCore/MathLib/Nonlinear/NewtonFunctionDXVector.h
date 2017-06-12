/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NewtonRaphson.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#pragma once

namespace MathLib
{

/**
 * \brief Vector function evaluating solution increment for Newton method 
 *
 * \tparam F_JACOBIAN   A function evaluating jacobian value
 */
template<class F_JACOBIAN, class T_LINEAR_SOLVER>
class NewtonFunctionDXVector
{
public:
    typedef typename T_LINEAR_SOLVER::MatrixType JacobianType;

    NewtonFunctionDXVector(F_JACOBIAN &f, T_LINEAR_SOLVER &linear_solver)
        : _f_j(&f), _linear_solver(&linear_solver) {};

    template<class T_VALUE>
    void eval(const T_VALUE &x, const T_VALUE &r, T_VALUE &dx)
    {
        const size_t n_rows = _linear_solver->getDimension();
        JacobianType* jac = _linear_solver->getA();
        _f_j->eval(x, *jac);
        //if (!check_jac(jac)) return false;

        // rhs = -r
        for (size_t i=0; i<n_rows; ++i)
            _linear_solver->setRHS(i, -1.*r[i]);

        // solve J dx = -r
        _linear_solver->solve();

        // get dx
        double *u = _linear_solver->getX();
        for (size_t i=0; i<n_rows; ++i)
            dx[i] = u[i];
    }

private:
    inline bool check_jac(JacobianType &jac)
    {
        const size_t n = _linear_solver->getDimension();
        for (size_t i=0; i<n; i++)
            for (size_t j=0; j<n; j++)
                if (jac(i,j)!=.0) return true;
        return false;
    }

private:
    F_JACOBIAN* _f_j;
    T_LINEAR_SOLVER* _linear_solver;
};

} //end

