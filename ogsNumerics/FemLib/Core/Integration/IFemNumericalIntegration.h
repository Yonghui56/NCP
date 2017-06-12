/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IFemNumericalIntegration.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MeshLib/Core/IElement.h"

namespace FemLib
{

/**
 * \brief Interface for integration methods used in FEM
 */
class IFemNumericalIntegration
{
public:
    virtual ~IFemNumericalIntegration() {};
    virtual void initialize(MeshLib::IElement &e, size_t n_sampl_level) = 0;
    virtual size_t getNumberOfSamplingPoints() const = 0;
    virtual void getSamplingPoint(size_t igp, double *) const = 0;
    virtual double getWeight(size_t igp) const = 0;
};

}
