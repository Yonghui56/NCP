/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file MeshNodesAlongPolyline.h
 *
 * Created on 2010-08-09 by Thomas Fischer
 */

#ifndef MESHNODESALONGPOLYLINE_H_
#define MESHNODESALONGPOLYLINE_H_

// GEOLIB
#include "GeoLib/Polyline.h"

// Mesh
#include "MeshLib/Core/IMesh.h"


namespace MeshLib
{
/**
 * This class computes the ids of the mesh nodes along a polyline.
 *
 * The mesh nodes are sorted as follow:
 * [ ... ids of sorted linear nodes ... | ... ids of unsorted higher order nodes ]
 */
class MeshNodesAlongPolyline
{
public:
    MeshNodesAlongPolyline(GeoLib::Polyline const* const ply, IMesh const* mesh);
    const std::vector<size_t>& getNodeIDs () const;
    const GeoLib::Polyline* getPolyline () const;
    size_t getNumberOfLinearNodes () const;
    std::vector<double> const & getDistOfProjNodeFromPlyStart() const;

private:
    const GeoLib::Polyline* _ply;
    const IMesh* _mesh;
    size_t _linear_nodes;
    std::vector<size_t> _msh_node_ids;
    std::vector<double> _dist_of_proj_node_from_ply_start;
};
}

#endif /* MESHNODESALONGPOLYLINE_H_ */
