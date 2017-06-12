/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TransientElementWiseMatrixUpdater.h
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
class TransientElementWiseMatrixUpdater
{
public:
    typedef T_LOCAL LocalAssemblerType;
    typedef DiscreteLib::IDiscreteVector<double> GlobalVector;

    TransientElementWiseMatrixUpdater(const TimeStep* time, MeshLib::IMesh* msh, const GlobalVector* u0, const GlobalVector* u1, LocalAssemblerType* a)
    : _msh(msh), _transient_e_assembler(a), _timestep(time), _vec_u0(u0), _vec_u1(u1)
    {

    }

    TransientElementWiseMatrixUpdater(const TransientElementWiseMatrixUpdater<LocalAssemblerType> &obj)
    {
        _msh = obj._msh;
        _transient_e_assembler = obj._transient_e_assembler;
        _timestep = obj._timestep;
        _vec_u0 = obj._vec_u0;
        _vec_u1 = obj._vec_u1;
    }

    TransientElementWiseMatrixUpdater<LocalAssemblerType> &operator=(const TransientElementWiseMatrixUpdater<LocalAssemblerType> &obj)
    {
        return * new TransientElementWiseMatrixUpdater<LocalAssemblerType>(obj);
    }

    template <class T_LINEAR_EQS>
    void update(const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &dofManager, T_LINEAR_EQS &eqs)
    {
        std::vector<size_t> ele_node_ids, ele_node_size_order;
        std::vector<size_t> local_dofmap_row, local_dofmap_column;
        MathLib::LocalEquation localEQS;
        MathLib::LocalVector local_u_n1;
        MathLib::LocalVector local_u_n;

        // get dof map
        e.getNodeIDList(e.getMaximumOrder(), ele_node_ids);
        dofManager.mapEqsIDreduced(_msh->getID(), ele_node_ids, local_dofmap_row, local_dofmap_column);
        // previous and current results
        DiscreteLib::getLocalVector(local_dofmap_row, *_vec_u1, local_u_n1);
        DiscreteLib::getLocalVector(local_dofmap_row, *_vec_u0, local_u_n);
        
        DiscreteLib::DofEquationIdTable localDofMap;
        dofManager.createLocalMappingTable(_msh->getID(), ele_node_ids, localDofMap);

        // local assembly
        localEQS.create(local_dofmap_row.size());
        _transient_e_assembler->assembly(*_timestep, e, localDofMap, local_u_n1, local_u_n, *localEQS.getA());
		//_transient_e_assembler->assembly_exact(*_timestep, e, localDofMap, local_u_n1, local_u_n, *localEQS.getA());
        // update global
        eqs.addAsub(local_dofmap_row, local_dofmap_column, *localEQS.getA());
    }

private:
    MeshLib::IMesh* _msh;
    LocalAssemblerType* _transient_e_assembler;
    const TimeStep* _timestep;
    const GlobalVector* _vec_u0;
    const GlobalVector* _vec_u1;
};

} //end
