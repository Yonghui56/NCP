/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Line.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "Point.h"
#include "Polyline.h"

namespace GeoLib
{

/**
 * \brief Line shape composed of two points and one polyline
 *
 */
class Line
{
public:
    /**
     * Line shape can be created by two points
     */
    Line(const Point &pt1, const Point &pt2);

    ///
    virtual ~Line();

    /// return polyline
    const Polyline& getPolyline() const;

    /// return point 1
    const Point& getPoint1() const;

    /// return point 2
    const Point& getPoint2() const;

private:
    std::vector<Point*> _pnt_vec;
    Polyline* _poly;
};

}
