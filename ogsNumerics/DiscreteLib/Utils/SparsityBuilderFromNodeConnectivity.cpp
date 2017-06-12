/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SparsityBuilderFromNodeConnectivity.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "SparsityBuilderFromNodeConnectivity.h"

#include "MeshLib/Topology/TopologyNode2NodesConnectedByElements.h"

#include "DiscreteLib/DDC/SubDomain.h"
#include "SparsityTool.h"

namespace DiscreteLib
{

SparsityBuilderFromNodeConnectivity::SparsityBuilderFromNodeConnectivity(const MeshLib::IMesh &msh, DofEquationIdTable &dofManager, MathLib::RowMajorSparsity &sparse)
{
    msh.setCurrentOrder(msh.getMaxiumOrder());
    MeshLib::TopologyNode2NodesConnectedByElements topo_node2nodes(msh);
    SparsityTool::createRowMajorSparsityFromNodeConnectivity(topo_node2nodes, msh.getID(), dofManager, sparse);
}

SparsityBuilderFromNodeConnectivity::SparsityBuilderFromNodeConnectivity(SubDomain &ddc_dom, const MeshLib::IMesh &msh, DofEquationIdTable &dofManager, MathLib::RowMajorSparsity &sparse)
{
    MeshLib::TopologyNode2NodesConnectedByElements topo_node2nodes(msh);
    SparsityTool::createRowMajorSparsityFromNodeConnectivity(topo_node2nodes, *ddc_dom.getGlobalLocalIdMap(), msh.getID(), dofManager, sparse);
}


}
