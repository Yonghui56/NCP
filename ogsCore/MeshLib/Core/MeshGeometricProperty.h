/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file MeshGeometricProperty
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <limits>

#include "GeoLib/AxisAlignedBoundingBox.h"
#include "CoordinateSystem.h"

namespace MeshLib
{

/**
 * \brief Mesh geometric property
 *
 * This class contains the following geometric information of a particular mesh
 * - coordinate system
 * - boundign box
 * - minimum edge length
 */
class MeshGeometricProperty
{
public:
    ///
    MeshGeometricProperty()
    {
        _min_edge_len = std::numeric_limits<double>::epsilon();
    }
    
    ///
    double getMinEdgeLength() const
    {
        return _min_edge_len;
    }

    ///
    void setMinEdgeLength(double l)
    {
        _min_edge_len = l;
    }

    ///
    const CoordinateSystem& getCoordinateSystem() const
    {
        return _coord;
    }

    ///
    void setCoordinateSystem(const CoordinateSystemType::type &coord)
    {
        _coord.setType(coord);
    }

    ///
    const GeoLib::AxisAlignedBoundingBox& getBoundingBox() const
    {
        return _aabox;
    }

    ///
    GeoLib::AxisAlignedBoundingBox& getBoundingBox()
    {
        return _aabox;
    }

private:
    double _min_edge_len;
    CoordinateSystem _coord;
    GeoLib::AxisAlignedBoundingBox _aabox;
};

}// end namespace

