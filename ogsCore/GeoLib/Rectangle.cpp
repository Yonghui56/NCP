/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Rectangle.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "Rectangle.h"

namespace GeoLib
{

Rectangle::Rectangle(const Point &pt_lower_left, const Point &pt_upper_right) {
    const double z = pt_lower_left[2];
    _pnt_vec.push_back(new Point(pt_lower_left));
    _pnt_vec.push_back(new Point(pt_upper_right[0], pt_lower_left[1], z));
    _pnt_vec.push_back(new Point(pt_upper_right));
    _pnt_vec.push_back(new Point(pt_lower_left[0], pt_upper_right[1], z));
    _poly_domain = new Polyline(_pnt_vec);
    _poly_domain->addPoint(0);
    _poly_domain->addPoint(1);
    _poly_domain->addPoint(2);
    _poly_domain->addPoint(3);
    Polyline *poly_left = new Polyline(_pnt_vec);
    poly_left->addPoint(0);
    poly_left->addPoint(3);
    Polyline *poly_right = new Polyline(_pnt_vec);
    poly_right->addPoint(1);
    poly_right->addPoint(2);
    Polyline *poly_top = new Polyline(_pnt_vec);
    poly_top->addPoint(2);
    poly_top->addPoint(3);
    Polyline *poly_bottom = new Polyline(_pnt_vec);
    poly_bottom->addPoint(0);
    poly_bottom->addPoint(1);
    _poly_vec.push_back(poly_left);
    _poly_vec.push_back(poly_right);
    _poly_vec.push_back(poly_top);
    _poly_vec.push_back(poly_bottom);
};

Rectangle::~Rectangle() {
    BaseLib::releaseObjectsInStdVector(_pnt_vec);
    BaseLib::releaseObjectsInStdVector(_poly_vec);
    BaseLib::releaseObject(_poly_domain);
}

const Polyline& Rectangle::getDomainPolyline() const 
{
    return *_poly_domain;
};

const Polyline& Rectangle::getLeft() const
{
    return *_poly_vec[0];
};

const Polyline& Rectangle::getRight() const
{
    return *_poly_vec[1];
};

const Polyline& Rectangle::getTop() const
{
    return *_poly_vec[2];
};

const Polyline& Rectangle::getBottom() const
{
    return *_poly_vec[3];
};


}
