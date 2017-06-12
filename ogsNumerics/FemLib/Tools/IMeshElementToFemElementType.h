/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IMeshElementToFemElementType.h
 *
 * Created on 2012-11-27 by Norihiro Watanabe
 */

#pragma once

#include "MeshLib/Core/IElement.h"

namespace FemLib
{

class IMeshElementToFemElementType
{
public:
    virtual ~IMeshElementToFemElementType() {};

    virtual int getFeType(const MeshLib::IElement& e, const int order) const = 0;
};


}
