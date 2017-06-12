/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TopologySequentialNodes2Elements.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <cstddef>
#include <vector>
#include <map>
#include <set>

namespace MeshLib
{
class IMesh;

/**
 * \brief Mesh topological data: node <-> connected elements
 */    
class TopologySequentialNodes2Elements
{
public:
    TopologySequentialNodes2Elements(const IMesh &msh);

    const std::vector<size_t>& getConnectedElements(size_t node_id) const {
        return _node2conn_eles[node_id];
    }

private:
    std::vector<std::vector<size_t> > _node2conn_eles;
};

}
