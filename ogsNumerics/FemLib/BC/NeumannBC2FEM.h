/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NeumannBC2FEM.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "GeoLib/GeoObject.h"
#include "MeshLib/Core/IMesh.h"
#include "NumLib/Function/TXFunction.h"
#include "FemLib/Tools/IFeObjectContainer.h"

namespace FemLib
{

/**
 * Neumann BC
 */
class NeumannBC2FEM
{
public:
    /// 
    NeumannBC2FEM(const MeshLib::IMesh &msh, const double &current_time, IFeObjectContainer &feObjects, const GeoLib::GeoObject &_geo, const NumLib::ITXFunction &_bc_func, std::vector<size_t> &_vec_nodes, std::vector<double> &_vec_values);
};



}
