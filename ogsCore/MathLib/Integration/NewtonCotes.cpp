/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NewtonCotes.cpp
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#include "NewtonCotes.h"

namespace MathLib
{

double NewtonCotes::middpoint(double (*fun)(double), double a, double b, std::size_t n) {
    double h = (b-a) / (double)n;
    double val = .0;
    double x = a + 0.5*h;
    for (std::size_t i=0; i<n; i++) {
        val += fun(x);
        x += h;
    }
    val *= h;
    return val;
};

double NewtonCotes::trapezoid(double (*fun)(double), double a, double b, std::size_t n) {
    double h = (b-a) / (double)n;
    double val = 0.5*(fun(a) + fun(b));
    double x = a;
    for (std::size_t i=0; i<n-1; i++) {
        x += h;
        val += fun(x);
    }
    val *= h;
    return val;
};

double NewtonCotes::simpson(double (*fun)(double), double a, double b, std::size_t n) {
    double h = (b-a) / (double)n;
    double val = fun(a) + 4*fun(b-0.5*h) + fun(b);
    double x = a;
    for (std::size_t i=0; i<n-1; i++) {
        val += 4*fun(x+0.5*h) + 2*fun(x+h);
        x += h;
    }
    val = val * h / 6.0;
    return val;
};

}
