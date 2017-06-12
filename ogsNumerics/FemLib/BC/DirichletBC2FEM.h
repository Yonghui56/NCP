/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DirichletBC2FEM.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */


#pragma once

#include <vector>

#include "GeoLib/GeoObject.h"
#include "MeshLib/Core/IMesh.h"
#include "NumLib/Function/TXFunction.h"

namespace FemLib
{

/**
 * \brief A class constructing DirichletBC data for FEM
 */
class DirichletBC2FEM
{
public:
    ///
    DirichletBC2FEM(const MeshLib::IMesh &msh, const GeoLib::GeoObject &geo, const NumLib::ITXFunction &bc_func, std::vector<size_t> &vec_nodes, std::vector<double> &vec_values);
};


}
