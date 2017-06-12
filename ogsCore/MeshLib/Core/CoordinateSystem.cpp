/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file CoordinateSystem.cpp
 *
 * Created on 2012-09-19 by Norihiro Watanabe
 */

#include "CoordinateSystem.h"

#include "BaseLib/CodingTools.h"

namespace MeshLib
{

/// has z dimension
bool CoordinateSystem::hasX() const {
    switch (_type) {
    case CoordinateSystemType::X:
    case CoordinateSystemType::XY:
    case CoordinateSystemType::XZ:
    case CoordinateSystemType::XYZ:
        return true;
    default:
        return false;
    }
}

size_t CoordinateSystem::getIndexOfX() const {
    if (hasX())
        return 0;
    else
        return BaseLib::index_npos;
}

/// has z dimension
bool CoordinateSystem::hasY() const {
    switch (_type) {
    case CoordinateSystemType::Y:
    case CoordinateSystemType::XY:
    case CoordinateSystemType::YZ:
    case CoordinateSystemType::XYZ:
        return true;
    default:
        return false;
    }
}

size_t CoordinateSystem::getIndexOfY() const {
    if (hasY()) {
        if (hasX())
            return 1;
        else
            return 0;
    } else {
        return BaseLib::index_npos;
    }
}

/// has z dimension
bool CoordinateSystem::hasZ() const {
    switch (_type) {
    case CoordinateSystemType::Z:
    case CoordinateSystemType::XZ:
    case CoordinateSystemType::YZ:
    case CoordinateSystemType::XYZ:
        return true;
    default:
        return false;
    }
}

size_t CoordinateSystem::getIndexOfZ() const {
    if (hasZ()) {
        return getDimension() - 1;
    } else {
        return BaseLib::index_npos;
    }
}

} // end
