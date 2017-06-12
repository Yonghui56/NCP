/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NewtonCotes.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#pragma once

#include <cstddef>

namespace MathLib
{

/**
 * \brief Newton Cotes quadrature method
 *
 */
class NewtonCotes
{
public:
    /// midpoint integration
    double middpoint(double (*fun)(double), double a, double b, std::size_t n);

    /// trapezoid integration
    double trapezoid(double (*fun)(double), double a, double b, std::size_t n);

    /// simpson integration
    double simpson(double (*fun)(double), double a, double b, std::size_t n);
};

}

