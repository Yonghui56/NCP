/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Tools.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "Tools.h"

#include <algorithm>
#include <exception>

#include "logog.hpp"

#include "GeoLib/GeoType.h"

#include "MeshLib/Core/ElementFactory.h"
#include "MeshLib/Tools/MeshNodesAlongPolyline.h"
#include "MeshLib/Tools/MeshNodesAlongSurface.h"
#include "MeshLib/Core/ElementCoordinatesInvariant.h"
#include "MeshLib/Core/ElementCoordinatesMappingLocal.h"
//#include "MeshLib/Topology/Topology.h"

namespace MeshLib
{

/// seeks all nodes located on a given polyline and returns a list of the found node pointers.
void findNodesOnPolyline(IMesh const* msh, GeoLib::Polyline const* poly, std::vector<size_t> *vec_nodes)
{
    MeshNodesAlongPolyline obj(poly, msh);
    std::vector<size_t> vec_node_id = obj.getNodeIDs();
    vec_nodes->assign(vec_node_id.begin(), vec_node_id.end());
};

void findNodesOnPoint(IMesh const* msh, GeoLib::Point const* point, std::vector<size_t> *vec_nodes)
{
    for (size_t i=0; i<msh->getNumberOfNodes(); ++i) {
        if (*msh->getNodeCoordinatesRef(i) == *point) {
            vec_nodes->push_back(i);
            break;
        }
    }
};

void findNodesOnSurface(IMesh const* msh, GeoLib::Surface const* sfc, std::vector<size_t> *vec_nodes)
{
    MeshNodesAlongSurface obj(sfc, msh);
    std::vector<size_t> vec_node_id = obj.getNodeIDs();
    vec_nodes->assign(vec_node_id.begin(), vec_node_id.end());
}

///
void findNodesOnGeometry(IMesh const* msh, GeoLib::GeoObject const* obj, std::vector<size_t> *vec_nodes)
{
    switch (obj->getGeoType()) {
        case GeoLib::POINT:
            findNodesOnPoint(msh, static_cast<GeoLib::Point const*>(obj), vec_nodes);
            break;
        case GeoLib::POLYLINE:
            findNodesOnPolyline(msh, static_cast<GeoLib::Polyline const*>(obj), vec_nodes);
            break;
        case GeoLib::SURFACE:
            findNodesOnSurface(msh, static_cast<GeoLib::Surface const*>(obj), vec_nodes);
            break;
        case GeoLib::GEODOMAIN:
            vec_nodes->resize(msh->getNumberOfNodes());
            for (size_t i=0; i<vec_nodes->size(); i++)
                (*vec_nodes)[i] = i;
            break;
        default:
            throw "This geo type is not supported in MeshLib::findNodesOnGeometry()";
            break;
    }
};

/// get a list of elements connected to one of give nodes
void findConnectedElements(IMesh const* msh, const std::vector<size_t> &nodes, std::vector<size_t> &connected_elements)
{
    for (size_t i=0; i<msh->getNumberOfElements(); i++) {
        IElement* e = msh->getElement(i);
        for (size_t j=0; j<e->getNumberOfNodes(); j++) {
            if (std::find(nodes.begin(), nodes.end(), e->getNodeID(j))!=nodes.end()) {
                connected_elements.push_back(e->getID());
                break;
            }
        }
    }
};

/// find edge elements on a given polyline
void findEdgeElementsOnPolyline(IMesh * msh, GeoLib::Polyline const* poly, std::vector<IElement*> *vec_edges_on_poly)
{
    // get a list of nodes on the polyline
    std::vector<size_t> nodes_on_poly;
    findNodesOnPolyline(msh, poly, &nodes_on_poly);
    // get a list of elements having the nodes
    std::vector<size_t> elements_near_poly;
    findConnectedElements(msh, nodes_on_poly, elements_near_poly);
    // get a list of edges made of the nodes
    std::vector<IElement*> selected_edges;
    createEdgeElements(msh, elements_near_poly, selected_edges);
    const size_t nr_edges = selected_edges.size();
    for (size_t i=0; i<nr_edges; i++) {
        IElement *edge_e = selected_edges[i];
        // check
        size_t cnt_match = 0;
        for (size_t j=0; j<edge_e->getNumberOfNodes(); j++) {
            //const INode* edge_nod = msh->getNode(edge_e->getNodeID(j));
            if (std::find(nodes_on_poly.begin(), nodes_on_poly.end(), edge_e->getNodeID(j)) != nodes_on_poly.end())
                cnt_match++;
            else
                break;
        }
        // update the list
        if (cnt_match==edge_e->getNumberOfNodes())
            vec_edges_on_poly->push_back(edge_e);
    }
};

/// find boundary elements located on given geometries
void findBoundaryElementsOnGeometry(IMesh * msh, GeoLib::GeoObject const* obj, std::vector<IElement*> *vec_eles)
{
    switch (obj->getGeoType()) {
        case GeoLib::POLYLINE:
            findEdgeElementsOnPolyline(msh, static_cast<GeoLib::Polyline const*>(obj), vec_eles);
            break;
        default:
            throw "This geo type is not supported in MeshLib::findNodesOnGeometry()";
            break;
    }
};

IElement* createEdgeElement(IMesh& msh, IElement &e, size_t j)
{
    std::vector<size_t> vec_edge_nodes;
    IElement* edge = ElemenetFactory::createNewElement(e.getEdgeType(j)); 
    e.getNodeIDsOfEdges(j, vec_edge_nodes);
    for (size_t k=0; k<vec_edge_nodes.size(); k++) {
        edge->setNodeID(k, vec_edge_nodes[k]);
    }
    msh.addEdgeElement(edge);
    return edge;
}

IElement* findEdgeElement(const std::vector<IElement*> &edges, std::vector<size_t> &vec_edge_nodes)
{
    std::sort(vec_edge_nodes.begin(), vec_edge_nodes.end());

    IElement *e_edge = 0;
    for (std::vector<IElement*>::const_iterator itr=edges.begin(); itr!=edges.end(); ++itr) {
        if ((*itr)->hasNodeIds(vec_edge_nodes)) {
            e_edge = *itr;
            break;
        }
    }
    return e_edge;
}

void findEdgeElements(IMesh& msh, IElement &e, std::vector<IElement*> &edges)
{
    std::vector<size_t> vec_edge_nodes;
    const size_t n_edges = e.getNumberOfEdges();
    for (size_t j=0; j<n_edges; j++) {
        //check if already exists
        e.getNodeIDsOfEdges(j, vec_edge_nodes);
        IElement* edge = findEdgeElement(edges, vec_edge_nodes);
        if (edge!=0) {
            if (e.getEdge(j)==0)
                e.setEdge(j, edge);
            continue;
        }

        edge = e.getEdge(j);
        if (edge==0) {
            //if new, create 
            if (edge==0) {
                edge = createEdgeElement(msh, e, j);
            }
            e.setEdge(j, edge);
        }
        edges.push_back(edge);

    }
}

/// get a list of edge elements of given elements
void createEdgeElements(IMesh * msh, const std::vector<size_t> &selected_ele, std::vector<IElement*> &edges)
{
    for (size_t i=0; i<selected_ele.size(); i++) {
        IElement *e = msh->getElement(selected_ele[i]);
        findEdgeElements(*msh, *e, edges);
    }
};

void createEdgeElements(IMesh * msh)
{
    std::vector<IElement*> edges;
    for (size_t i=0; i<msh->getNumberOfElements(); i++) {
        IElement *e = msh->getElement(i);
        findEdgeElements(*msh, *e, edges);
    }
};



MeshLib::CoordinateSystemType::type getCoordinateSystemFromBoundingBox(const GeoLib::AxisAlignedBoundingBox &bbox)
{
    GeoLib::Point pt_diff = bbox.getMaxPoint() - bbox.getMinPoint();
    MeshLib::CoordinateSystemType::type coords;
    bool hasX = fabs(pt_diff[0]) > .0;
    bool hasY = fabs(pt_diff[1]) > .0;
    bool hasZ = fabs(pt_diff[2]) > .0;

    if (hasX) {
        if (hasY) {
            if (hasZ) {
                coords = CoordinateSystemType::XYZ;
            } else {
                coords = CoordinateSystemType::XY;
            }
        } else if (hasZ) {
            coords = CoordinateSystemType::XZ;
        } else {
            coords = CoordinateSystemType::X;
        }
    } else if (hasY) {
        if (hasZ) {
            coords = CoordinateSystemType::YZ;
        } else {
            coords = CoordinateSystemType::Y;
        }
    } else if (hasZ) {
        coords = CoordinateSystemType::Z;
    } else {
        coords = CoordinateSystemType::INVALID;
    }

    return coords;
}

double calculateMeshMinimumEdgeLength(UnstructuredMesh &msh)
{
    std::vector<size_t> vec_edge_nodes;
    double min_edge_len = std::numeric_limits<double>::max();
    const size_t n_ele = msh.getNumberOfElements();
    for (size_t i=0; i<n_ele; i++) {
        MeshLib::IElement* e = msh.getElement(i);
        for (size_t j=0; j<e->getNumberOfEdges(); j++) {
            e->getNodeIDsOfEdges(j, vec_edge_nodes);
            assert (vec_edge_nodes.size() == 2);

            double edge_len = std::sqrt(GeoLib::sqrDist(msh.getNodeCoordinatesRef(vec_edge_nodes[0]), msh.getNodeCoordinatesRef(vec_edge_nodes[1])));
            min_edge_len = std::min(min_edge_len, edge_len);
        }
    }

    return min_edge_len;
}

void calculateMeshGeometricProperties(UnstructuredMesh &msh)
{
    MeshGeometricProperty* geo_prop = msh.getGeometricProperty();
    //double tol = std::numeric_limits<double>::epsilon();

    // coordinate systems
    geo_prop->setCoordinateSystem(getCoordinateSystemFromBoundingBox(geo_prop->getBoundingBox()));

    //
//    GeoLib::Point pt_diff = geo_prop->getBoundingBox().getMaxPoint() - geo_prop->getBoundingBox().getMinPoint();
//    double max_len = std::max(pt_diff[0], pt_diff[1]);
//    max_len = std::max(max_len, pt_diff[2]);
//    double min_edge_len = max_len / msh.getNumberOfNodes();
    double min_edge_len = calculateMeshMinimumEdgeLength(msh);
    geo_prop->setMinEdgeLength(min_edge_len);

    INFO("-> calculate mesh geometric properties");
    INFO("* min. edge length = %f", min_edge_len);

}

void setMeshElementCoordinatesMapping(IMesh &msh)
{
    const size_t msh_dim = msh.getDimension();
    for (size_t i=0; i<msh.getNumberOfElements(); i++) {
        IElement* e = msh.getElement(i);
        if (e->getMappedCoordinates()==NULL) {
            MeshLib::IElementCoordinatesMapping* ele_map;
            size_t ele_dim = e->getDimension();
            assert(msh_dim >= ele_dim);
            if (msh_dim == ele_dim) {
                ele_map = new ElementCoordinatesInvariant(&msh, e);
            } else {
                ele_map = new ElementCoordinatesMappingLocal(&msh, *e, msh.getGeometricProperty()->getCoordinateSystem());
            }
            e->setMappedCoordinates(ele_map);
        }

    }
};

} // end namespace
