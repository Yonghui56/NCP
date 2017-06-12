/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Tools.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include "DiscreteLib/Serial/DiscreteVector.h"
#include "MathLib/DataType.h"
#include "DiscreteLib/Utils/DofEquationIdTable.h"

namespace DiscreteLib
{

/// create a subset of vector u corresponding to the given vector index
void getLocalVector2(const DofEquationIdTable &dofManager, const std::vector<size_t> &list_vec_entry_id, const std::vector<size_t> &list_vec_size_for_order, const std::vector<IDiscreteVector<double>*> &list_multiple_u, std::vector<double> &local_u);

/// create a subset of vector u corresponding to the given vector index
void getLocalVector(const DofEquationIdTable &dofManager, const std::vector<size_t> &list_vec_entry_id, const std::vector<size_t> &list_vec_size_for_order, const std::vector<IDiscreteVector<double>*> &list_multiple_u, MathLib::LocalVector &local_u);

void getLocalVector(const std::vector<size_t> &list_vec_entry_id, const IDiscreteVector<double> &global_u, MathLib::LocalVector &local_u);

void getLocalVector(const std::vector<size_t> &list_vec_entry_id, const MathLib::LocalVector &global_u, MathLib::LocalVector &local_u);

/// set
void setGlobalVector(const DofEquationIdTable &dofManager, size_t var_id, size_t mesh_id, const IDiscreteVector<double> &u, IDiscreteVector<double> &global_vec);

void setLocalVector(const DofEquationIdTable &dofManager, size_t var_id, size_t mesh_id, const IDiscreteVector<double> &global_vec, IDiscreteVector<double> &u);

void convertToEqsValues(const DiscreteLib::DofEquationIdTable &eqs_map, size_t var_id, size_t msh_id, const std::vector<size_t> &list_node_id, const std::vector<double> &list_node_values, std::vector<size_t> &list_eqs_id, std::vector<double> &list_eqs_val);
} //end

