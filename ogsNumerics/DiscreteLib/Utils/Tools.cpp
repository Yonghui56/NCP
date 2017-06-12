/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Tools.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "Tools.h"

namespace DiscreteLib
{


/// create a subset of vector u corresponding to the given vector index
void getLocalVector2(const DofEquationIdTable &dofManager, const std::vector<size_t> &list_vec_entry_id, const std::vector<size_t> &list_vec_size_for_order, const std::vector<IDiscreteVector<double>*> &list_multiple_u, std::vector<double> &local_u)
{
    local_u.clear();
    const size_t n_var = dofManager.getNumberOfVariables();
    for (size_t i=0; i<n_var; i++) {
        size_t var_order = 1; //TODO
        const size_t n_dof = list_vec_size_for_order[var_order-1];
        const DiscreteLib::IDiscreteVector<double> &var_u = *list_multiple_u[i];
        for (size_t j=0; j<n_dof; j++) {
            local_u.push_back(var_u[list_vec_entry_id[j]]);
        }
    }
    //const size_t n_dof = getNumberOfVariables();
    //for (size_t i=0; i<n_dof; i++) {
    //    const DofMap *dofMap = getVariableDoF(i);
    //    const size_t vec_size = list_vec_size_for_order[dofMap->getOrder()-1];
    //    const DiscreteVector<double> &dof_u = *list_multiple_u[i];
    //    for (size_t j=0; j<vec_size; j++) {
    //        local_u.push_back(dof_u[list_vec_entry_id[j]]);
    //    }
    //}
}

/// create a subset of vector u corresponding to the given vector index
void getLocalVector(const DofEquationIdTable &dofManager, const std::vector<size_t> &list_vec_entry_id, const std::vector<size_t> &list_vec_size_for_order, const std::vector<IDiscreteVector<double>*> &list_multiple_u, MathLib::LocalVector &local_u)
{
    std::vector<double> temp_v;
    getLocalVector2(dofManager, list_vec_entry_id, list_vec_size_for_order, list_multiple_u, temp_v);
    local_u.resize(temp_v.size());
    for (size_t i=0; i<temp_v.size(); i++)
        local_u[i] = temp_v[i];
}

/**
 * extract entries and create subset
 *
 * @param list_vec_entry_id     a list of entry id
 * @param global_u              global vector
 * @param local_u               subset
 */
void getLocalVector(const std::vector<size_t> &list_vec_entry_id, const IDiscreteVector<double> &global_u, MathLib::LocalVector &local_u)
{
    size_t valid_entry_cnt = 0;
    for (size_t i=0; i<list_vec_entry_id.size(); i++) {
        if (list_vec_entry_id[i] != BaseLib::index_npos)
            valid_entry_cnt++;
    }
    local_u.resize(valid_entry_cnt);
    valid_entry_cnt = 0;
    for (size_t i=0; i<list_vec_entry_id.size(); i++) {
        if (list_vec_entry_id[i] != BaseLib::index_npos)
            local_u[valid_entry_cnt++] = global_u[list_vec_entry_id[i]];
    }
}

void getLocalVector(const std::vector<size_t> &list_vec_entry_id, const MathLib::LocalVector &global_u, MathLib::LocalVector &local_u)
{
    size_t valid_entry_cnt = 0;
    for (size_t i=0; i<list_vec_entry_id.size(); i++) {
        if (list_vec_entry_id[i] != BaseLib::index_npos)
            valid_entry_cnt++;
    }
    local_u.resize(valid_entry_cnt);
    valid_entry_cnt = 0;
    for (size_t i=0; i<list_vec_entry_id.size(); i++) {
        if (list_vec_entry_id[i] != BaseLib::index_npos)
            local_u[valid_entry_cnt++] = global_u[list_vec_entry_id[i]];
    }
}

void setGlobalVector(const DofEquationIdTable &dofManager, size_t var_id, size_t mesh_id, const IDiscreteVector<double> &u, IDiscreteVector<double> &global_vec)
{
    for (size_t i=u.getRangeBegin(); i<u.getRangeEnd(); i++) {
        size_t eqs_id = dofManager.mapEqsID(var_id, mesh_id, i);
        global_vec[eqs_id] = u[i];
    }
}

void setLocalVector(const DofEquationIdTable &dofManager, size_t var_id, size_t mesh_id, const IDiscreteVector<double> &global_vec, IDiscreteVector<double> &u)
{
    for (size_t i=u.getRangeBegin(); i<u.getRangeEnd(); i++) {
        size_t eqs_id = dofManager.mapEqsID(var_id, mesh_id, i);
        u[i] = global_vec[eqs_id];
    }
}

void convertToEqsValues(const DiscreteLib::DofEquationIdTable &eqs_map, size_t var_id, size_t msh_id, const std::vector<size_t> &list_node_id, const std::vector<double> &list_node_values, std::vector<size_t> &list_eqs_id, std::vector<double> &list_eqs_val)
{
    for (size_t j=0; j<list_node_id.size(); j++) {
        size_t pt_id = list_node_id[j];
        if (eqs_map.isActiveDoF(var_id, msh_id, pt_id)) {
            size_t eqs_id = eqs_map.mapEqsID(var_id, msh_id, pt_id);
            double bc_val = list_node_values[j];
            list_eqs_id.push_back(eqs_id);
            list_eqs_val.push_back(bc_val);
        }
    }
}

} //end

