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

#ifndef COMPOUND_H_
#define COMPOUND_H_

#include "BaseLib/CodingTools.h"

namespace NumLib
{
class ITXFunction;
}

namespace Geo
{

struct Compound
{
    NumLib::ITXFunction* molecular_diffusion;

    Compound()
    {
        BaseLib::zeroObject(molecular_diffusion);
    }
    ~Compound()
    {
        BaseLib::releaseObject(molecular_diffusion);
    }

};

} //end

#endif /* COMPOUND_H_ */
