/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file MeshNodesAlongSurface.cpp
 *
 * Created on 2012-11-14 by Norihiro Watanabe
 */

#include "MeshNodesAlongSurface.h"

// c++
#include <algorithm>

// Base
#include "BaseLib/quicksort.h"
// Math
#include "MathLib/Vector.h"
#include "MathLib/MathTools.h"
// Mesh
#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Core/Node.h"
#include "MeshLib/Core/MeshGeometricProperty.h"


namespace MeshLib
{

MeshNodesAlongSurface::MeshNodesAlongSurface(GeoLib::Surface const* const sfc, IMesh const* mesh) :
    _sfc(sfc), _mesh(mesh), _linear_nodes (0)
{
    //std::vector<Node*> const& mesh_nodes (mesh->getNodeVector());
    const double min_edge_length (mesh->getGeometricProperty()->getMinEdgeLength() * 0.1);

    //size_t n_linear_order_nodes (mesh->getNumberOfNodes ());
    const size_t n_nodes (mesh->getNumberOfNodes());

    const_cast<GeoLib::Surface*>(sfc)->initSurfaceGrid();

    for (size_t j(0); j < n_nodes; j++) {
        const GeoLib::Point* nod_x = mesh->getNodeCoordinatesRef(j);
        if (sfc->isPntInBV(nod_x->getData(), min_edge_length)) {
            if (sfc->isPntInSfc(nod_x->getData(), min_edge_length)) {
                _msh_node_ids.push_back(j);
            }
        }
    }

}

const std::vector<size_t>& MeshNodesAlongSurface::getNodeIDs () const
{
    return _msh_node_ids;
}

const GeoLib::Surface* MeshNodesAlongSurface::getSurface () const
{
    return _sfc;
}

size_t MeshNodesAlongSurface::getNumberOfLinearNodes () const
{
    return _linear_nodes;
}

//const std::vector<double>& MeshNodesAlongSurface::getDistOfProjNodeFromPlyStart() const
//{
//    return _dist_of_proj_node_from_ply_start;
//}

} // end namespace
