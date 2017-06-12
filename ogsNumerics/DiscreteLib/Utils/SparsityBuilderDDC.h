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

#include <set>

#include "MathLib/LinAlg/Sparse/SparseTableCRS.h"

#include "MeshLib/Topology/TopologyNode2NodesConnectedByElements.h"

#include "DiscreteLib/Utils/DofEquationIdTable.h"
#include "DiscreteLib/DDC/DecomposedDomain.h"
#include "DiscreteLib/DDC/SubDomain.h"
#include "DiscreteLib/DDC/IGlobaLocalMappingTable.h"
#include "SparsityTool.h"

namespace DiscreteLib
{
//class SparsityBuilderFromNodeConnectivityDDC
//{
//public:
//    SparsityBuilderFromNodeConnectivityDDC(SubDomain &ddc_dom, MeshLib::IMesh &msh, DofEquationIdTable &dofManager, MathLib::RowMajorSparsity &sparse)
//    {
//        MeshLib::TopologyNode2NodesConnectedByElements topo_node2nodes(&msh);
//        SparsityTool::createRowMajorSparsityFromNodeConnectivity(topo_node2nodes, *ddc_dom.getGlobalLocalIdMap(), msh.getID(), dofManager, sparse);
//    }
//};

class SparsityBuilderFromNodeConnectivityWithGhostDoFs
{
public:
    SparsityBuilderFromNodeConnectivityWithGhostDoFs(MeshLib::IMesh &msh, DofEquationIdTable &dofManager, MathLib::RowMajorSparsity &sparse)
    {
        MeshLib::TopologyNode2NodesConnectedByElements topo_node2nodes(msh);
        if (dofManager.getNumberOfVariables()==1) {
            SparsityTool::createRowMajorSparsityFromNodeConnectivity(topo_node2nodes, sparse);
        } else {
            SparsityTool::createRowMajorSparsityForMultipleDOFs(topo_node2nodes, dofManager.getNumberOfVariables(), sparse);
        }
        size_t n_rows = dofManager.getTotalNumberOfActiveDoFsWithoutGhost();
        sparse.erase(sparse.begin()+n_rows, sparse.end());
    }
};

//class SparsityBuilderFromNodeConnectivityWithInactiveDoFs
//{
//public:
//    SparsityBuilderFromNodeConnectivityWithInactiveDoFs(MeshLib::IMesh &msh, DofEquationIdTable &dofManager, MathLib::RowMajorSparsity &sparse)
//    {
//        MeshLib::TopologyNode2NodesConnectedByElements topo_node2nodes(&msh);
//        if (dofManager.getNumberOfVariables()==1) {
//            SparsityBuilder::createRowMajorSparsityFromNodeConnectivity(topo_node2nodes, *dofManager.getVariableDoF(0), sparse);
//        } else {
//            SparsityBuilder::createRowMajorSparsityForMultipleDOFs(topo_node2nodes, dofManager.getNumberOfVariables(), sparse);
//        }
//    }
//};

//class SparsityBuilderFromLocalSparsity
//{
//public:
//    SparsityBuilderFromLocalSparsity(std::vector<MathLib::RowMajorSparsity*> &list_local_sparse, DofEquationIdTable &dofManager, MathLib::RowMajorSparsity &sparse)
//    {
//
//    }
//};

template <class T_SPARSITY_BUILDER>
class SparsityBuilderFromDDC
{
public:
    SparsityBuilderFromDDC(DecomposedDomain &ddc_global, DofEquationIdTable &dofManager, MathLib::RowMajorSparsity &sparse)
    {
        const size_t n_dom = ddc_global.getNumberOfSubDomains();
        const size_t n_dofs = dofManager.getTotalNumberOfActiveDoFs();
        sparse.resize(n_dofs);

        std::vector<MathLib::RowMajorSparsity> local_sparse(ddc_global.getNumberOfSubDomains());
        for (size_t i=0; i<n_dom; i++) {
            SubDomain* dom = ddc_global.getSubDomain(i);
            T_SPARSITY_BUILDER builder(*dom, *dom->getLoalMesh(), dofManager, local_sparse[i]);
        }

        for (size_t i=0; i<n_dom; i++) {
            SubDomain* dom = ddc_global.getSubDomain(i);
            IGlobaLocalMappingTable* mapping = dom->getGlobalLocalIdMap();
            MathLib::RowMajorSparsity &local_sp = local_sparse[i];
            for (size_t j=0; j<local_sp.size(); j++) {
                const std::set<size_t> &local_row = local_sp[j];
                if (local_row.size()==0) continue;
                const size_t i_gloal = mapping->local2global(j);
                std::set<size_t> &global_row = sparse[i_gloal];
                for (std::set<size_t>::const_iterator itr=local_row.begin(); itr!=local_row.end(); ++itr) {
                    global_row.insert(*itr);
                    //global_row.insert(mapping->local2global(*itr));
                }
            }
        }

        //size_t i_row = 0;
        //for (size_t i=0; i<local_sparse.size(); i++) {
        //    size_t offset = 0;
        //    MathLib::RowMajorSparsity &local_sp = local_sparse[i];
        //    for (size_t j=0; j<local_sp.size(); j++) {
        //        std::set<size_t> &setConnection = sparse[i_row++];
        //        std::set<size_t> &local_conn = local_sp[j];
        //        for (std::set<size_t>::iterator itr = local_conn.begin(); itr!=local_conn.end(); ++itr) {
        //            setConnection.insert(offset + *itr);
        //        }
        //    }
        //    offset += local_sp.size();
        //}
    }
};

}
