/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DenseLinearEquations.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#pragma once

#include <valarray>
#include <vector>

#include "AbstractDenseLinearEquation.h"


namespace MathLib
{

/**
 * \brief Linear equation based on dense matrix
 */
class DenseLinearEquation : public AbstractDenseLinearEquation
{
public:
    void initialize() {};
    void finalize() {};

    void setOption(const BaseLib::Options &option);

    void solve();
};


}
