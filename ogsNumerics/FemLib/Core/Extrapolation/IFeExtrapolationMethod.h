/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IFeExtrapolationMethod.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "MathLib/DataType.h"

namespace FemLib
{

class IFiniteElement;

/**
 * \brief Extrapolation of integration point values to nodal values
 */
class IFeExtrapolationMethod
{
public:
    virtual ~IFeExtrapolationMethod() {};
    virtual void extrapolate(IFiniteElement &fe, const std::vector<MathLib::LocalVector> &gp_values, std::vector<MathLib::LocalVector> &nodal_values) = 0;
};

} // end
