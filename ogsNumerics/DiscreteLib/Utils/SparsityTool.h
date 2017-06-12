/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SparsityBuilder.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MathLib/LinAlg/Sparse/Sparsity.h"

#include "MeshLib/Topology/ITopologyNode2Nodes.h"

#include "DiscreteLib/Utils/DofEquationIdTable.h"
#include "DiscreteLib/DDC/IGlobaLocalMappingTable.h"

namespace DiscreteLib
{

namespace SparsityTool
{
/**
 * create row-major sparsity based on node connectivities
 */
void createRowMajorSparsityFromNodeConnectivity(const MeshLib::ITopologyNode2Nodes &topo_node2nodes, MathLib::RowMajorSparsity &row_major_entries);

void createRowMajorSparsityFromNodeConnectivity(const MeshLib::ITopologyNode2Nodes &topo_node2nodes, size_t mesh_id, DofEquationIdTable &local_dofTable, MathLib::RowMajorSparsity &row_major_entries);

void createRowMajorSparsityFromNodeConnectivity(const MeshLib::ITopologyNode2Nodes &local_topo_node2nodes, IGlobaLocalMappingTable &pt_mapping, size_t mesh_id, DofEquationIdTable &global_dofTable, MathLib::RowMajorSparsity &row_major_entries);

/**
 * create row-major sparsity for multiple DOFs with a single mesh
 */
void createRowMajorSparsityForMultipleDOFs(const MeshLib::ITopologyNode2Nodes &topo_node2nodes, size_t n_dof, MathLib::RowMajorSparsity &row_major_entries);

}; // SparsityTool
}; //DiscreteLib
