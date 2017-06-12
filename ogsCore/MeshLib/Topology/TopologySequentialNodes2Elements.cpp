/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TopologySequentialNodes2Elements.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "TopologySequentialNodes2Elements.h"

#include <algorithm>

#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Core/IElement.h"

namespace MeshLib
{

TopologySequentialNodes2Elements::TopologySequentialNodes2Elements(const IMesh &msh) 
{
    _node2conn_eles.resize(msh.getNumberOfNodes());
    const size_t nr_ele = msh.getNumberOfElements();
    for (size_t i=0; i<nr_ele; i++) {
        const IElement* e = msh.getElement(i);
        for (size_t j=0; j<e->getNumberOfNodes(); j++) {
            _node2conn_eles[e->getNodeID(j)].push_back(i);
        }
    }
}

} // end namespace

