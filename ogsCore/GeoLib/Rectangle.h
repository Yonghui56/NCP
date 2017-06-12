/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Rectangle.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "BaseLib/CodingTools.h"

#include "Point.h"
#include "Polyline.h"

namespace GeoLib
{

/**
 * \brief Rectangle shape
 *
 */
class Rectangle
{
public:
    /// Create this shape by specifying the lower left and lower right points
    Rectangle(const Point &pt_lower_left, const Point &pt_upper_right);

    ///
    virtual ~Rectangle();

    ///
    const Polyline& getDomainPolyline() const;

    ///
    const Polyline& getLeft() const;

    ///
    const Polyline& getRight() const;

    ///
    const Polyline& getTop() const;

    ///
    const Polyline& getBottom() const;

private:
    std::vector<Point*> _pnt_vec;
    Polyline* _poly_domain;
    std::vector<Polyline*> _poly_vec;
};

}
