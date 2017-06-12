/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file MeshNodesAlongSurface.h
 *
 * Created on 2012-11-14 by Norihiro Watanabe
 */

#ifndef MESHNODESALONGSURFACE_H_
#define MESHNODESALONGSURFACE_H_

// GEOLIB
#include "GeoLib/Surface.h"

// Mesh
#include "MeshLib/Core/IMesh.h"


namespace MeshLib
{
/**
 * This class computes the ids of the mesh nodes along a surface.
 *
 * The mesh nodes are sorted as follow:
 * [ ... ids of sorted linear nodes ... | ... ids of unsorted higher order nodes ]
 */
class MeshNodesAlongSurface
{
public:
    MeshNodesAlongSurface(GeoLib::Surface const* const sfc, IMesh const* mesh);
    const std::vector<size_t>& getNodeIDs () const;
    const GeoLib::Surface* getSurface () const;
    size_t getNumberOfLinearNodes () const;
    std::vector<double> const & getDistOfProjNodeFromPlyStart() const;

private:
    const GeoLib::Surface* _sfc;
    const IMesh* _mesh;
    size_t _linear_nodes;
    std::vector<size_t> _msh_node_ids;
    //std::vector<double> _dist_of_proj_node_from_ply_start;
};
}

#endif /* MESHNODESALONGSURFACE_H_ */
