/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file PorousMedia.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"

#include "Solid.h"

namespace NumLib
{
class ITXFunction;
}

namespace Geo
{

struct PorousMedia
{
    NumLib::ITXFunction* hydraulic_conductivity;
    NumLib::ITXFunction* porosity;
    NumLib::ITXFunction* storage;
    Solid* solidphase;

    PorousMedia()
    {
        BaseLib::zeroObject(hydraulic_conductivity, porosity, storage, solidphase);
    }
    ~PorousMedia()
    {
        BaseLib::releaseObject(hydraulic_conductivity, porosity, storage);
    }
};

} //end
 
