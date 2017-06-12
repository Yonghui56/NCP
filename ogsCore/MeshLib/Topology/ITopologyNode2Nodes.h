/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ITopologyNode2Nodes.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <cstddef>
#include <set>

namespace MeshLib
{
/**
 * \brief Interface to node-node connectivity
 */ 
class ITopologyNode2Nodes
{
public:
    virtual ~ITopologyNode2Nodes() {};

    /// get a list of connected nodes
    virtual const std::set<size_t>& getConnectedNodes(size_t node_id) const = 0;

    /// get a total number of nodes
    virtual size_t getNumberOfNodes() const = 0;
};


}
