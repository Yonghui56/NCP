/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementWiseVectorUpdater.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include "MeshLib/Core/IElement.h"
#include "DiscreteLib/Utils/DofEquationIdTable.h"

namespace DiscreteLib
{

/**
 * \brief Element-wise vector updater
 *
 * \tparam T_VALUE  vector data type
 * \tparam T_LOCAL  local assembler
 */
template <class T_VALUE, class T_LOCAL>
class ElementWiseVectorUpdater
{
public:
    typedef T_LOCAL LocalAssemblerType;
    typedef IDiscreteVectorAssembler<T_VALUE>::VectorType GlobalVectorType;

    /**
     *
     * @param msh
     * @param dofManager
     * @param a
     */
    ElementWiseVectorUpdater(MeshLib::IMesh* msh, LocalAssemblerType* a)
    : _msh(msh), _e_assembler(a)
    {

    }

    /**
     *
     * @param e
     * @param globalVec
     */
    void update(const MeshLib::IElement &e, const DofEquationIdTable &dofManager, GlobalVectorType &globalVec)
    {
        LocalVector localVec;
        std::vector<size_t> ele_node_ids, ele_node_size_order;
        std::vector<long> local_dofmap_row;//, local_dofmap_column;

        std::vector<T_VALUE> local_u_n;
        // get dof map
        e.getNodeIDList(e.getMaximumOrder(), ele_node_ids);
        e.getListOfNumberOfNodesForAllOrders(ele_node_size_order);
        dofManager.mapEqsID(_msh->getID(), ele_node_ids, local_dofmap_row); //TODO order
        // local assembly
        localVec.resize(local_dofmap_row.size(), .0);
        _e_assembler->assembly(*e, localVec);
        // update global
        globalVec.addSubvector(local_dofmap_row, &localVec[0]);
    }

private:
    MeshLib::IMesh* _msh;
    LocalAssemblerType* _e_assembler;
};

} //end
