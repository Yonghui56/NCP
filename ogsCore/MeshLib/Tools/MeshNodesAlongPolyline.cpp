/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file MeshNodesAlongPolyline.cpp
 *
 * Created on 2010-08-09 by Thomas Fischer
 */

#include "MeshNodesAlongPolyline.h"

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

MeshNodesAlongPolyline::MeshNodesAlongPolyline(GeoLib::Polyline const* const ply, IMesh const* mesh) :
    _ply(ply), _mesh(mesh), _linear_nodes (0)
{
    //std::vector<Node*> const& mesh_nodes (mesh->getNodeVector());
    double min_edge_length (mesh->getGeometricProperty()->getMinEdgeLength() * 0.1);

    //size_t n_linear_order_nodes (mesh->getNumberOfNodes ());
    size_t n_nodes (mesh->getNumberOfNodes());

    //std::vector<size_t> msh_node_higher_order_ids;
    //std::vector<double> dist_of_proj_higher_order_node_from_ply_start;

    // repeat until at least one relevant node was found
    while (_msh_node_ids.empty())
    {
        // loop over all line segments of the polyline
        for (size_t k = 0; k < ply->getNumberOfPoints() - 1; k++)
        {
            double act_length_of_ply(ply->getLength(k));
            // loop over all nodes
            for (size_t j = 0; j < n_nodes; j++)
            {
                double dist, lambda;

                // is the orthogonal projection of the j-th node to the
                // line g(lambda) = _ply->getPoint(k) + lambda * (_ply->getPoint(k+1) - _ply->getPoint(k))
                // at the k-th line segment of the polyline, i.e. 0 <= lambda <= 1?
                const GeoLib::Point* nod_x = mesh->getNodeCoordinatesRef(j);
                if (MathLib::calcProjPntToLineAndDists(nod_x->getData(), (_ply->getPoint(k))->getData(), (_ply->getPoint(k + 1))->getData(), lambda, dist) 
                        < min_edge_length)
                    if (0 <= lambda && lambda <= 1)
                    {
                        //if (mesh_nodes[j]->getNodeID() < n_linear_order_nodes)
                        //{
                            // check if node id is already in the vector
                            if (std::find (_msh_node_ids.begin(), _msh_node_ids.end(), j)
                                == _msh_node_ids.end())
                            {
                                _msh_node_ids.push_back(j);
                                _dist_of_proj_node_from_ply_start.push_back (act_length_of_ply + dist);
                                _linear_nodes++;
                            }
                        //}
                        //else
                        //// check if node id is already in the vector
                        //if (std::find (msh_node_higher_order_ids.begin(), msh_node_higher_order_ids.end(), mesh_nodes[j]->getNodeID()) ==  msh_node_higher_order_ids.end())
                        //{
                        //    msh_node_higher_order_ids.push_back (mesh_nodes[j]->getNodeID());
                        //    dist_of_proj_higher_order_node_from_ply_start.push_back (act_length_of_ply + dist);
                        //}
                    } // end if lambda
            } // end node loop
        } // end line segment loop

        if (_msh_node_ids.empty())
            min_edge_length *= 2.0;
    }

    // sort the (linear) nodes along the polyline according to their distances
    BaseLib::Quicksort<double> (_dist_of_proj_node_from_ply_start, 0, _dist_of_proj_node_from_ply_start.size(), _msh_node_ids);

//#ifndef NDEBUG
//    std::cout << "*****" << std::endl;
//    std::cout << "distances of linear nodes along polyline " << *ply <<
//    "(min_edge_length = " << min_edge_length << "): " << std::endl;
//    for (size_t k(0); k < _dist_of_proj_node_from_ply_start.size(); k++)
//        std::cout << "\t" << _msh_node_ids[k] << " " <<
//        _dist_of_proj_node_from_ply_start[k] << std::endl;
//    //std::cout << "number of linear nodes along polyline " << ply << ": " <<
//    //_dist_of_proj_node_from_ply_start.size()
//    //          << ", number of higher order nodes: " << msh_node_higher_order_ids.size() <<
//    //std::endl;
//#endif
    //// assign/append higher order nodes at the end of vector _msh_node_ids
    //for (size_t k(0); k < msh_node_higher_order_ids.size(); k++)
    //    _msh_node_ids.push_back (msh_node_higher_order_ids[k]);
    //// append distances for higher order nodes at the end of vector _dist_of_proj_node_from_ply_start
    //for (size_t k(0); k < dist_of_proj_higher_order_node_from_ply_start.size(); k++)
    //    _dist_of_proj_node_from_ply_start.push_back (dist_of_proj_higher_order_node_from_ply_start[k]);
}

const std::vector<size_t>& MeshNodesAlongPolyline::getNodeIDs () const
{
    return _msh_node_ids;
}

const GeoLib::Polyline* MeshNodesAlongPolyline::getPolyline () const
{
    return _ply;
}

size_t MeshNodesAlongPolyline::getNumberOfLinearNodes () const
{
    return _linear_nodes;
}

const std::vector<double>& MeshNodesAlongPolyline::getDistOfProjNodeFromPlyStart() const
{
    return _dist_of_proj_node_from_ply_start;
}

} // end namespace
