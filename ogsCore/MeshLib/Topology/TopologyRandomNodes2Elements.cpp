/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TopologyRandomNodes2Elements.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "TopologyRandomNodes2Elements.h"

#include <algorithm>

#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Core/IElement.h"

namespace MeshLib
{

TopologyRandomNodes2Elements::TopologyRandomNodes2Elements(const IMesh &msh, const std::vector<size_t> &nodes) 
{
    for (size_t i=0; i<msh.getNumberOfElements(); i++) {
        IElement* e = msh.getElement(i);
        for (size_t j=0; j<e->getNumberOfNodes(); j++) {
            if (std::find(nodes.begin(), nodes.end(), e->getNodeID(j))!=nodes.end()) {
                _map_node2conn_eles[e->getID()].push_back(i);
                break;
            }
        }
    }
}

} // end namespace

