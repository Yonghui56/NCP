/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Node.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "GeoLib/Point.h"

namespace MeshLib
{

/**
 * \brief Mesh node class
 */
class Node
{
public:
    ///
    Node() : _node_id(0) {};

    ///
    Node (size_t id, const GeoLib::Point &x) : _node_id(id), _x(x) {};

    ///
    size_t getNodeID() const { return _node_id;};
    ///
    void setNodeID(const size_t &id) { _node_id = id;};

    ///
    const GeoLib::Point* getX() const  { return &_x; };
    ///
    void setX(const GeoLib::Point &x) { _x = x; };

private:
    size_t _node_id;
    GeoLib::Point _x;
};

} // end namespace

