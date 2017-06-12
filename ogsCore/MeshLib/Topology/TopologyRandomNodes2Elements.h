/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TopologyRandomNodes2Elements.h
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

class TopologyRandomNodes2Elements
{
public:
    TopologyRandomNodes2Elements(const IMesh &msh, const std::vector<size_t> &nodes);

    const std::vector<size_t>* getConnectedElements(size_t node_id) const {
        std::map<size_t, std::vector<size_t> >::const_iterator itr = _map_node2conn_eles.find(node_id);
        if (itr!=_map_node2conn_eles.end())
            return &itr->second;
        return 0;
    }

private:
    std::map<size_t, std::vector<size_t> > _map_node2conn_eles;
};

}
