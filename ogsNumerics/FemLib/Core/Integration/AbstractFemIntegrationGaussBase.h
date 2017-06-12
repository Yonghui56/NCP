/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file AbstractFemIntegrationGaussBase.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MeshLib/Core/IElement.h"

#include "IFemNumericalIntegration.h"

namespace FemLib
{

/**
 * \brief Gauss–Legendre integration class
 * 
 * 
 */
class AbstractFemIntegrationGaussBase : public IFemNumericalIntegration
{
public:
    AbstractFemIntegrationGaussBase()
    {
        _n_sampl_level = 0;
        _n_sampl_pt = 0;
    }

    void initialize(MeshLib::IElement &e, size_t n_sampl_level)
    {
        _n_sampl_level = n_sampl_level;
        _n_sampl_pt = getTotalNumberOfSamplingPoints(e, n_sampl_level);
    };

    size_t getSamplingLevel() const {return _n_sampl_level;};
    size_t getNumberOfSamplingPoints() const {return _n_sampl_pt;};

private:
    size_t _n_sampl_level;
    size_t _n_sampl_pt;

    virtual size_t getTotalNumberOfSamplingPoints(MeshLib::IElement &e, size_t n_sampl_level) const = 0;

};

}
