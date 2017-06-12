/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemIntegrationAnalytical.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <algorithm>
#include "MeshLib/Core/IElement.h"
#include "IFemNumericalIntegration.h"

namespace FemLib
{

/**
 * \brief Dummy class for analytical integration methods which are not using sampling.
 */
class FemIntegrationAnalytical : public IFemNumericalIntegration
{
public:
    virtual ~FemIntegrationAnalytical() {};

    virtual void initialize(MeshLib::IElement &/*e*/, size_t)
    {
    };

    virtual size_t getNumberOfSamplingPoints() const {return 1;};

    virtual void getSamplingPoint(size_t, double*) const 
    {
    };
    
    virtual double getWeight(size_t) const {return 1.0;};
};

}
