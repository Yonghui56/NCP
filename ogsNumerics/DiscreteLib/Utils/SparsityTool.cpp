/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SparsityTool.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "SparsityTool.h"

#include <set>


namespace DiscreteLib
{

namespace SparsityTool
{
/**
 * create row-major sparsity based on node connectivities
 */
void createRowMajorSparsityFromNodeConnectivity(const MeshLib::ITopologyNode2Nodes &topo_node2nodes, MathLib::RowMajorSparsity &row_major_entries)
{
    const size_t n_nodes = topo_node2nodes.getNumberOfNodes();
    row_major_entries.resize(n_nodes);

    for (size_t i=0; i<n_nodes; i++) {
        // search connected nodes
        std::set<size_t> &setConnection = row_major_entries[i];
        setConnection.insert(i);
        const std::set<size_t> connected_nodes = topo_node2nodes.getConnectedNodes(i);
        for (std::set<size_t>::const_iterator it=connected_nodes.begin(); it!=connected_nodes.end(); ++it) {
            setConnection.insert(*it);
        }
    }
}

#if 0
void createRowMajorSparsityFromNodeConnectivity(const MeshLib::ITopologyNode2Nodes &topo_node2nodes, size_t mesh_id, DofEquationIdTable &local_dofTable, MathLib::RowMajorSparsity &row_major_entries)
{
    row_major_entries.resize(local_dofTable.getTotalNumberOfActiveDoFs());

    const size_t n_nodes = topo_node2nodes.getNumberOfNodes();
    const size_t n_var = local_dofTable.getNumberOfVariables();
    for (size_t i=0; i<n_nodes; i++) {
        const std::set<size_t> connected_nodes = topo_node2nodes.getConnectedNodes(i);
        for (size_t j=0; j<n_var; j++) {
            // for each DoF(var,pt)
            std::set<size_t> &setConnection = row_major_entries[i*n_var+j];
            // add all DoFs defined in this point
            for (size_t l=0; l<n_var; l++) {
                size_t col_id = local_dofTable.mapEqsID(l, mesh_id, i);
                setConnection.insert(col_id);
            }
            // add DoFs defined in connected nodes
            for (std::set<size_t>::const_iterator it=connected_nodes.begin(); it!=connected_nodes.end(); it++) {
                for (size_t l=0; l<n_var; l++) {
                    size_t col_id = local_dofTable.mapEqsID(l, mesh_id, *it);
                    setConnection.insert(col_id);
                }
            }
        }
    }
}
#else
void createRowMajorSparsityFromNodeConnectivity(const MeshLib::ITopologyNode2Nodes &topo_node2nodes_max_order, size_t mesh_id, DofEquationIdTable &local_dofTable, MathLib::RowMajorSparsity &row_major_entries)
{
    const size_t n_rows = local_dofTable.getTotalNumberOfActiveDoFs();
    row_major_entries.resize(n_rows);

    size_t n_nodes_max = topo_node2nodes_max_order.getNumberOfNodes();

    const size_t n_var = local_dofTable.getNumberOfVariables();
    for (size_t pt_id=0; pt_id<n_nodes_max; pt_id++) {
        for (size_t var_id=0; var_id<n_var; var_id++) {
//            if (!local_dofTable.getPointEquationIdTable(var_id, mesh_id)->hasKey(pt_id))
//                continue;
//            size_t row_id = pt_id*n_var + var_id;
            size_t row_id = local_dofTable.mapEqsID(var_id, mesh_id, pt_id);
            if (row_id == BaseLib::index_npos) continue;
            // for each DoF(var,pt)
            std::set<size_t> &setConnection = row_major_entries[row_id];
            // add all DoFs defined in this point
            for (size_t l=0; l<n_var; l++) {
                size_t col_id = local_dofTable.mapEqsID(l, mesh_id, pt_id);
                if (col_id != BaseLib::index_npos)
                    setConnection.insert(col_id);
            }
            // add DoFs defined in connected nodes
            const std::set<size_t> connected_nodes = topo_node2nodes_max_order.getConnectedNodes(pt_id);
            for (std::set<size_t>::const_iterator it=connected_nodes.begin(); it!=connected_nodes.end(); ++it) {
                for (size_t l=0; l<n_var; l++) {
                    size_t col_id = local_dofTable.mapEqsID(l, mesh_id, *it);
                    if (col_id != BaseLib::index_npos)
                        setConnection.insert(col_id);
                }
            }
        }
    }
}
#endif

void createRowMajorSparsityFromNodeConnectivity(const MeshLib::ITopologyNode2Nodes &local_topo_node2nodes, IGlobaLocalMappingTable &pt_mapping, size_t mesh_id, DofEquationIdTable &global_dofTable, MathLib::RowMajorSparsity &row_major_entries)
{
    row_major_entries.resize(global_dofTable.getTotalNumberOfActiveDoFs());

    const size_t n_nodes = local_topo_node2nodes.getNumberOfNodes();
    for (size_t i=0; i<n_nodes; i++) {
        const std::set<size_t> connected_nodes = local_topo_node2nodes.getConnectedNodes(i);
        const size_t global_i = pt_mapping.local2global(i);
        for (size_t j=0; j<global_dofTable.getNumberOfVariables(); j++) {
            // for each DoF(var,pt)
            //long row_id = global_dofTable.mapEqsID(j, mesh_id, global_i);
            std::set<size_t> &setConnection = row_major_entries[i];
            // add all DoFs defined in this point
            for (size_t l=0; l<global_dofTable.getNumberOfVariables(); l++) {
                size_t col_id = global_dofTable.mapEqsID(l, mesh_id, global_i);
                setConnection.insert(col_id);
            }
            // add DoFs defined in connected nodes
            for (std::set<size_t>::const_iterator it=connected_nodes.begin(); it!=connected_nodes.end(); ++it) {
                const size_t global_pt_id = pt_mapping.local2global(*it);
                for (size_t l=0; l<global_dofTable.getNumberOfVariables(); l++) {
                    size_t col_id = global_dofTable.mapEqsID(l, mesh_id, global_pt_id);
                    setConnection.insert(col_id);
                }
            }
        }
    }
}

/**
 * create row-major sparsity for multiple DOFs with a single mesh
 */
void createRowMajorSparsityForMultipleDOFs(const MeshLib::ITopologyNode2Nodes &topo_node2nodes, size_t n_dof, MathLib::RowMajorSparsity &row_major_entries)
{
    const size_t n_nodes = topo_node2nodes.getNumberOfNodes();
    row_major_entries.resize(n_nodes*n_dof);

    MathLib::RowMajorSparsity sparsity_msh;
    createRowMajorSparsityFromNodeConnectivity(topo_node2nodes, sparsity_msh);

    for (size_t i_dof=0; i_dof<n_dof; i_dof++) {
        for (size_t i=0; i<n_nodes; i++) {
            std::set<size_t> &msh_row_entries = sparsity_msh[i];
            std::set<size_t> &row_entries = row_major_entries[i];
            for (std::set<size_t>::iterator it=msh_row_entries.begin(); it!=msh_row_entries.end(); ++it) {
                size_t nod_id = *it;
                for (size_t j=0; j<n_dof; j++) {
                    size_t eqs_id = 0;
                    if (0) {
                        eqs_id = nod_id*n_dof + j;
                    } else {
                        eqs_id = i + j*n_nodes;
                    }
                    row_entries.insert(eqs_id);
                }
            }
        }
    }
}
};

}
