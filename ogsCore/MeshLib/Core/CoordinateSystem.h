/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file CoordinateSystem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <cassert>
#include <cstddef>

namespace MeshLib
{

/**
 * \brief Coordinate system type
 */
struct CoordinateSystemType
{
    enum type {
        X = 10,
        Y = 11,
        Z = 12,
        XY = 21,
        XZ = 22,
        YZ = 23,
        XYZ = 32,
        INVALID = -1
    };
};

/**
 * \brief Coordinate systems
 *
 *
 */
class CoordinateSystem
{
public:

    ///
    CoordinateSystem() {
        _type = CoordinateSystemType::INVALID;
    }

    ///
    explicit CoordinateSystem(CoordinateSystemType::type coord) {
        _type = coord;
    }

    ///
    void setType(CoordinateSystemType::type coord) {
        _type = coord;
    }

    /// get this coordinate type
    CoordinateSystemType::type getType() const {
        return _type;
    }

    /// get dimension size
    size_t getDimension() const {
        return _type / 10;
    }

    /// has z dimension
    bool hasX() const;

    /// get index of x components
    size_t getIndexOfX() const;

    /// has z dimension
    bool hasY() const;

    /// get index of y components
    size_t getIndexOfY() const;

    /// has z dimension
    bool hasZ() const;

    /// get index of z components
    size_t getIndexOfZ() const;

private:
    CoordinateSystemType::type _type; 
};

}
