/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TransientElementWiseVectorUpdater.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "MeshLib/Core/IElement.h"
#include "DiscreteLib/Utils/DofEquationIdTable.h"
#include "DiscreteLib/Utils/Tools.h"

namespace NumLib
{

template <class T_LOCAL>
class TransientElementWiseVectorUpdater
{
public:
    typedef T_LOCAL LocalAssemblerType;
    typedef DiscreteLib::IDiscreteVector<double> GlobalVector;

    TransientElementWiseVectorUpdater(const TimeStep* time, MeshLib::IMesh* msh, DiscreteLib::DofEquationIdTable* /*dofManager*/, const GlobalVector* u0, const GlobalVector* u1, LocalAssemblerType* a)
    : _msh(msh), _transient_e_assembler(a), _timestep(time), _vec_u0(u0), _vec_u1(u1)
    {

    }

    TransientElementWiseVectorUpdater(const TransientElementWiseVectorUpdater<LocalAssemblerType> &obj)
    {
        _msh = obj._msh;
        _transient_e_assembler = obj._transient_e_assembler;
        _timestep = obj._timestep;
        _vec_u0 = obj._vec_u0;
        _vec_u1 = obj._vec_u1;
    }

    TransientElementWiseVectorUpdater<LocalAssemblerType> &operator=(const TransientElementWiseVectorUpdater<LocalAssemblerType> &obj)
    {
        return * new TransientElementWiseVectorUpdater<LocalAssemblerType>(obj);
    }

    void update(const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &dofManager, GlobalVector &globalVec)
    {
        std::vector<size_t> ele_node_ids, ele_node_size_order;
        std::vector<size_t> local_dofmap_row, local_dofmap_column;
        MathLib::LocalVector localVec;
        MathLib::LocalVector local_u_n1;
        MathLib::LocalVector local_u_n;

        e.getNodeIDList(e.getMaximumOrder(), ele_node_ids);
        e.getListOfNumberOfNodesForAllOrders(ele_node_size_order);
        dofManager.mapEqsIDreduced(_msh->getID(), ele_node_ids, local_dofmap_row, local_dofmap_column);

        // previous and current results
        DiscreteLib::getLocalVector(local_dofmap_row, *_vec_u1, local_u_n1);
        DiscreteLib::getLocalVector(local_dofmap_row, *_vec_u0, local_u_n);
        
        DiscreteLib::DofEquationIdTable localDofMap;
        dofManager.createLocalMappingTable(_msh->getID(), ele_node_ids, localDofMap);
        
        // local assembly
        localVec = MathLib::LocalVector::Zero(local_dofmap_row.size());
        _transient_e_assembler->assembly(*_timestep, e, localDofMap, local_u_n1, local_u_n, localVec);
        // update global

		globalVec.addSubvector(local_dofmap_row, &localVec[0]);// the address of localVec[0]
		
    }

private:
    MeshLib::IMesh* _msh;
    LocalAssemblerType* _transient_e_assembler;
    const TimeStep* _timestep;
    const GlobalVector* _vec_u0;
    const GlobalVector* _vec_u1;
};

} //end
