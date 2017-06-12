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

#include "GeoLib/GeoObject.h"
#include "GeoLib/Polyline.h"
#include "GeoLib/Surface.h"

#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Core/UnstructuredMesh.h"
#include "MeshLib/Core/StructuredMesh.h"

namespace MeshLib
{

/// seeks all nodes located on a given object and returns a list of the found node pointers.
void findNodesOnGeometry(IMesh const* msh, GeoLib::GeoObject const* obj, std::vector<size_t> *vec_nodes);

/// find nodes on a given point
void findNodesOnPoint(IMesh const* msh, GeoLib::Point const* point, std::vector<size_t> *vec_nodes);

/// seeks all nodes located on a given polyline and returns a list of the found node pointers.
void findNodesOnPolyline(IMesh const* msh, GeoLib::Polyline const* poly, std::vector<size_t> *vec_nodes);

/// seeks all nodes located on a given surface and returns a list of the found node pointers.
void findNodesOnSurface(IMesh const* msh, GeoLib::Surface const* surf, std::vector<size_t> *vec_nodes);

/// find boundary elements located on given geometries
void findBoundaryElementsOnGeometry(IMesh * msh, GeoLib::GeoObject const* obj, std::vector<IElement*> *vec_eles);

/// find edge elements on a given polyline
void findEdgeElementsOnPolyline(IMesh * msh, GeoLib::Polyline const* poly, std::vector<IElement*> *vec_edges_on_poly);



/// get a list of elements connected to one of give nodes
void findConnectedElements(IMesh const* msh, const std::vector<size_t> &nodes, std::vector<size_t> &connected_elements);

IElement* createEdgeElement(IMesh& msh, IElement &e, size_t edge_id);

IElement* findEdgeElement(const std::vector<IElement*> &edges, std::vector<size_t> &vec_edge_nodes);

void findEdgeElements(IMesh& msh, IElement &e, std::vector<IElement*> &edges);

/// get a list of edge elements of given elements
void createEdgeElements(IMesh * msh, const std::vector<size_t> &selected_ele, std::vector<IElement*> &edges);

/// get a list of edge elements of given elements
void createEdgeElements(IMesh * msh);

/// calculate geometric properties of an unstructured mesh
void calculateMeshGeometricProperties(UnstructuredMesh &msh);

/// set element coordinates mapping
void setMeshElementCoordinatesMapping(IMesh &msh);

/// calculate minimum edge length in a given mesh
double calculateMeshMinimumEdgeLength(UnstructuredMesh &msh);


}// end namespace
