/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file MeshGenerator.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <memory>
#include "BaseLib/BidirectionalMap.h"

#include "MeshLib/Core/Element.h"
#include "MeshLib/Core/UnstructuredMesh.h"
#include "MeshLib/Core/StructuredMesh.h"

namespace MeshLib
{

namespace MeshGenerator
{

/**
 * generate a line mesh
 *
 * \param length
 * \param subdivision
 * \param origin_x
 * \param origin_y
 * \param origin_z
 */
MeshLib::UnstructuredMesh* generateLineMesh(const double length, const size_t subdivision, const double origin_x, const double origin_y, const double origin_z);

/**
 * generate a reguar quad mesh
 *
 * \param length
 * \param subdivision
 * \param origin_x
 * \param origin_y
 * \param origin_z
 */
MeshLib::UnstructuredMesh* generateRegularQuadMesh(const double length, const size_t subdivision, const double origin_x, const double origin_y, const double origin_z);

/**
 * generate a reguar quad mesh but structured
 *
 * \param length
 * \param subdivision
 * \param origin_x
 * \param origin_y
 * \param origin_z
 */
StructuredMesh<ElementShape::QUAD>* generateStructuredRegularQuadMesh(const double length, const size_t subdivision, const double origin_x, const double origin_y, const double origin_z);

/**
 * generate a sub mesh from the original mesh
 *
 * \param src               Original mesh
 * \param list_e            List of elements making a new mesh
 * \param dest              New mesh
 * \param map_global2local  Mapping of element id between two meshes
 */
void generateSubMesh(const MeshLib::IMesh &src, const std::vector<size_t> &list_e, MeshLib::IMesh* &dest, BaseLib::BidirectionalMap<size_t, size_t> &map_global2local);

/// generate higher order mesh
void generateHigherOrderUnstrucuredMesh(UnstructuredMesh &msh, size_t order);


};

}
