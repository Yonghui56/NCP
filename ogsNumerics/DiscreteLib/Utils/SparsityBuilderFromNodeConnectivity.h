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
#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/Utils/DofEquationIdTable.h"

namespace DiscreteLib
{
class SubDomain;

class SparsityBuilderFromNodeConnectivity
{
public:
    SparsityBuilderFromNodeConnectivity(const MeshLib::IMesh &msh, DofEquationIdTable &dofManager, MathLib::RowMajorSparsity &sparse);
    SparsityBuilderFromNodeConnectivity(SubDomain &ddc_dom, const MeshLib::IMesh &msh, DofEquationIdTable &dofManager, MathLib::RowMajorSparsity &sparse);
};

}
