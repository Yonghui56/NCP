/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemIntegrationGaussQuad.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MathLib/Integration/GaussLegendre.h"
#include "MeshLib/Core/IElement.h"

#include "AbstractFemIntegrationGaussBase.h"

namespace FemLib
{

class FemIntegrationGaussQuad : public AbstractFemIntegrationGaussBase
{
public:
    void getSamplingPoint(size_t igp, double* x) const
    {
        size_t nGauss = getSamplingLevel();
        size_t gp_r = igp/nGauss;
        size_t gp_s = igp%nGauss;
        x[0] = MathLib::GaussLegendre::getPoint(nGauss, gp_r);
        x[1] = MathLib::GaussLegendre::getPoint(nGauss, gp_s);
    }

    double getWeight(size_t igp) const
    {
        size_t nGauss = getSamplingLevel();
        size_t gp_r = igp/nGauss;
        size_t gp_s = igp%nGauss;
        return  MathLib::GaussLegendre::getWeight(nGauss, gp_r)*MathLib::GaussLegendre::getWeight(nGauss, gp_s);
    }

private:
    size_t getTotalNumberOfSamplingPoints(MeshLib::IElement&, size_t n_sampl_level) const
    {
        return n_sampl_level*n_sampl_level;
    }
};


}
