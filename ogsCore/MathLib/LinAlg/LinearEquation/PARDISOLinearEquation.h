/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file PARDISOLinearEquation.h
 *
 * Original work is by Chan-Hee Park
 * Moved the work here by Norihiro Watanabe on 2012-08-23
 */

#ifdef USE_PARDISO

#pragma once

#include "AbstractCRSLinearEquation.h"

namespace MathLib
{

/**
 * \brief Linear equation class using PARDISO solver
 *
 */
class PARDISOLinearEquation
: public AbstractCRSLinearEquation<signed>
{
public:
    virtual ~PARDISOLinearEquation() {};

    virtual void setOption(const BaseLib::Options &/*option*/) {};

protected:
    virtual void solveEqs(CRSMatrix<double, signed> *A, double *rhs, double *x);

};

}

#endif

