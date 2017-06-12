/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemIntegrationGaussLine.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MathLib/Integration/GaussLegendre.h"
#include "MeshLib/Core/IElement.h"

#include "AbstractFemIntegrationGaussBase.h"

namespace FemLib
{

class FemIntegrationGaussLine : public AbstractFemIntegrationGaussBase
{
public:
    void getSamplingPoint(size_t igp, double* x) const
    {
        x[0] = MathLib::GaussLegendre::getPoint(getSamplingLevel(), igp);
    }

    double getWeight(size_t igp) const
    {
        return MathLib::GaussLegendre::getWeight(getSamplingLevel(), igp);
    }

private:
    size_t getTotalNumberOfSamplingPoints(MeshLib::IElement&, size_t n_sampl_level) const
    {
        return n_sampl_level;
    }
};

}
