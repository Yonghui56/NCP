/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Line.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "Line.h"

namespace GeoLib
{

Line::Line(const Point &pt1, const Point &pt2)
{
    _pnt_vec.push_back(new Point(pt1));
    _pnt_vec.push_back(new Point(pt2));
    _poly = new Polyline(_pnt_vec);
    _poly->addPoint(0);
    _poly->addPoint(1);
}

Line::~Line() {
    BaseLib::releaseObjectsInStdVector(_pnt_vec);
    BaseLib::releaseObject(_poly);
}

const Polyline& Line::getPolyline() const
{
    return *_poly;
}

const Point& Line::getPoint1() const
{
    return *_pnt_vec[0];
}

const Point& Line::getPoint2() const
{
    return *_pnt_vec[1];
}


}
