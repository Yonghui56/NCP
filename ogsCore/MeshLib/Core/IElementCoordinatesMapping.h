/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementCoordinatesMapping.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "GeoLib/Point.h"

namespace MeshLib
{

/**
 * \brief Interface of methods mapping coordinates of elements.
 *
 *
 */
class IElementCoordinatesMapping
{
public:
    ///
    virtual ~IElementCoordinatesMapping() {};
    
    /// get mapped coordinates of the node
    virtual GeoLib::Point* getNodePoint(size_t node_id) = 0;
};

}
