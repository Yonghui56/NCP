/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TopologyNode2NodesConnectedByElements.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "TopologyNode2NodesConnectedByElements.h"

#include <algorithm>

#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Core/IElement.h"

namespace MeshLib
{

TopologyNode2NodesConnectedByElements::TopologyNode2NodesConnectedByElements(const IMesh &msh) 
{
    _node2conn_nodes.resize(msh.getNumberOfNodes());
    const size_t nr_ele = msh.getNumberOfElements();
    for (size_t i=0; i<nr_ele; i++) {
        const IElement* e = msh.getElement(i);
        const size_t e_nnodes = e->getNumberOfNodes();
        for (size_t j=0; j<e_nnodes; j++) {
            std::set<size_t> &set_conn_nodes = _node2conn_nodes[e->getNodeID(j)];
            for (size_t k=0; k<e_nnodes; k++) {
                if (j==k) continue;
                set_conn_nodes.insert(e->getNodeID(k));
            }
        }
    }
}

} // end namespace

