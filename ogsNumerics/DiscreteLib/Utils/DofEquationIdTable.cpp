/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DofEquationIdTable.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "DofEquationIdTable.h"
#include <limits>

namespace DiscreteLib
{

void DofEquationIdTable::construct(long offset)
{
    if (_numbering_type==DofNumberingType::BY_VARIABLE) {
        //order by dof
        long eqs_id = offset;
        for (size_t i=0; i<_map_var2dof.size(); i++) {
            std::map<size_t, IEquationIdStorage*> &obj = _map_var2dof[i];
            for (std::map<size_t, IEquationIdStorage*>::iterator itr = obj.begin(); itr!=obj.end(); ++itr) {
                IEquationIdStorage* pt2eq = itr->second;
                eqs_id = pt2eq->setAll(eqs_id);
            }
        }
        _total_dofs_without_ghosts = eqs_id; //TODO
        _total_dofs = eqs_id;
    } else {
        //order by discrete points
        long eqs_id = offset;
        if (_is_seq_address) {
            assert(_ghost_pt.size()==0);
            for (std::map<size_t, std::vector<size_t> >::iterator itr=_map_msh2var.begin(); itr!=_map_msh2var.end(); ++itr) {
                size_t mesh_id = itr->first;
                std::vector<size_t> &list_var = itr->second;
                for (size_t i=0; i<list_var.size(); i++) {
                    IEquationIdStorage* pt2eq = getPointEquationIdTable(list_var[i], mesh_id);
                    eqs_id = offset + i;
                    eqs_id = pt2eq->setAll(eqs_id, list_var.size());
                }
            }
            _total_dofs_without_ghosts = eqs_id; //TODO
            _total_dofs = eqs_id;
        } else { //random address
            for (std::map<size_t, std::vector<size_t> >::iterator itr=_map_msh2var.begin(); itr!=_map_msh2var.end(); ++itr) {
                size_t mesh_id = itr->first;
                std::vector<size_t> &list_var = itr->second;
                // get range of point id
                size_t i_min = std::numeric_limits<size_t>::max();
                size_t i_max = 0;
                for (size_t i=0; i<list_var.size(); i++) {
                    IEquationIdStorage* pt2eq = getPointEquationIdTable(list_var[i], mesh_id);
                    size_t tmp_min, tmp_max;
                    pt2eq->key_range(tmp_min, tmp_max);
                    i_min = std::min(i_min, tmp_min);
                    i_max = std::max(i_max, tmp_max);
                }
                // first, real nodes
                for (size_t i=i_min; i<i_max+1; i++) {
                    if (isGhostPoint(mesh_id, i)) continue; //skip ghost
                    for (size_t j=0; j<list_var.size(); j++) {
                        IEquationIdStorage* pt2eq = getPointEquationIdTable(list_var[j], mesh_id);
                        if (pt2eq->hasKey(i) && pt2eq->isActive(i)) {
                            pt2eq->set(i, eqs_id++);
                        }
                    }
                }
                _total_dofs_without_ghosts = eqs_id;
                // next ghost nodes
                for (size_t i=i_min; i<i_max+1; i++) {
                    if (!isGhostPoint(mesh_id, i)) continue; //skip real
                    for (size_t j=0; j<list_var.size(); j++) {
                        IEquationIdStorage* pt2eq = getPointEquationIdTable(list_var[j], mesh_id);
                        if (pt2eq->hasKey(i) && pt2eq->isActive(i)) {
                            pt2eq->set(i, eqs_id++);
                        }
                    }
                }
            }
            _total_dofs = eqs_id;
        }
    }
}

//void DofEquationIdTable::mapEqsID(const std::vector<size_t> &ele_node_ids, const std::vector<size_t> &ele_node_size_order, std::vector<long> &local_dofmap) const
//{
//    const size_t n_dof = getNumberOfVariables();
//    if (n_dof<2) {
//        const DofMap *dofMap = getVariableDoF(0);
//        dofMap->getListOfEqsID(ele_node_ids, ele_node_size_order[dofMap->getOrder()-1], local_dofmap);
//    } else {
//        for (size_t j=0; j<n_dof; j++) {
//            std::vector<long> tmp_map;
//            const DofMap *dofMap = getVariableDoF(j);
//            dofMap->getListOfEqsID(ele_node_ids, ele_node_size_order[dofMap->getOrder()-1], tmp_map);
//            local_dofmap.insert(local_dofmap.end(), tmp_map.begin(), tmp_map.end());
//        }
//    }
//}

//size_t DofEquationIdTable::getTotalNumberOfDOFs(const std::vector<size_t> &list_vec_entry_id, const std::vector<size_t> &list_vec_size_for_order) const
//{
//    size_t s = 0;
//    //const size_t n_dof = getNumberOfVariables();
//    //for (size_t i=0; i<n_dof; i++) {
//    //    const DofMap *dofMap = getVariableDoF(i);
//    //    s += list_vec_size_for_order[dofMap->getOrder()-1];
//    //}
//    return s;
//}


} //end
