/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TopologyNode2NodesConnectedByEdges.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <cstddef>
#include <vector>
#include <set>
#include "ITopologyNode2Nodes.h"

namespace MeshLib
{
class IMesh;


/**
 * \brief Mesh topological data: node <-> connected nodes vie edges
 */
class TopologyNode2NodesConnectedByEdges : public ITopologyNode2Nodes
{
public:
    TopologyNode2NodesConnectedByEdges(const IMesh& msh);

    virtual ~TopologyNode2NodesConnectedByEdges() {};

    const std::set<size_t>& getConnectedNodes(size_t node_id) const 
    {
        return _node2conn_nodes[node_id];
    }

    size_t getNumberOfNodes() const 
    {
        return _node2conn_nodes.size();
    }

private:
    std::vector<std::set<size_t> > _node2conn_nodes;
};

}
