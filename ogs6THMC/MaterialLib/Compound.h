/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Compound.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "NumLib/Function/TXFunction.h"

namespace MaterialLib
{

struct Compound
{
    /**
      * name of the component
      */
    std::string name;

    /**
      * whether this component is moving
      */
    bool is_mobile;

    /**
      * charge of the component, in the unit of equivalent
      */
    double charge; 

    std::string comp_type; 

    NumLib::ITXFunction* molecular_diffusion;

    Compound()
    : is_mobile(true), molecular_diffusion(nullptr), 
      charge(0.0)
    {
    }
    ~Compound()
    {
        BaseLib::releaseObject(molecular_diffusion);
    }

};

} //end

