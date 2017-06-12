/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file MeshGenerator.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "MeshGenerator.h"

#include <memory>
#include <set>
#include <cstddef>
#include "GeoLib/Point.h"
#include "MeshLib/Topology/TopologySequentialNodes2Elements.h"
#include "MeshLib/Tools/Tools.h"

namespace MeshLib
{

MeshLib::UnstructuredMesh* MeshGenerator::generateLineMesh(const double length, const size_t subdivision, const double origin_x, const double origin_y, const double origin_z) 
{
    MeshLib::UnstructuredMesh* msh = new MeshLib::UnstructuredMesh(MeshLib::CoordinateSystemType::X);

    const double unit_length = length / subdivision;
    const size_t n_nodes_per_axis = subdivision+1;
    const size_t n_eles = subdivision;
    //const size_t n_nodes = subdivision+1;

    //nodes
    size_t node_id(0);
    for (size_t i_z=0; i_z<n_nodes_per_axis; i_z++) {
        const double x = unit_length*i_z;
        GeoLib::Point p(x+origin_x, origin_y, origin_z);
        msh->setNodeCoordinates(node_id++, p);
        //msh->setNodeCoordinates(node_id++, GeoLib::Point(x+origin_x, origin_y, origin_z));
    }

    //elements
    //size_t ele_id(0);
    for (size_t i_z=0; i_z<n_eles; i_z++) {
        Line *e = new Line();
        e->setNodeID(0, i_z);
        e->setNodeID(1, i_z+1);
        msh->addElement(e);
    }

    return msh;
}

MeshLib::UnstructuredMesh* MeshGenerator::generateRegularQuadMesh(const double length, const size_t subdivision, const double origin_x, const double origin_y, const double origin_z) 
{

    MeshLib::UnstructuredMesh* msh = new MeshLib::UnstructuredMesh(MeshLib::CoordinateSystemType::XY);

    //size_t n_eles = static_cast<size_t>(pow(static_cast<double>(subdivision), 2));
    //size_t n_nodes = static_cast<size_t>(pow(static_cast<double>(subdivision+1), 2));
    const double unit_length = length / subdivision;
    const size_t n_nodes_per_axis = subdivision+1;

    //nodes
    size_t node_id(0);
    const double z = origin_z;
    for (size_t j_y=0; j_y<n_nodes_per_axis; j_y++) {
        const double y = unit_length*j_y + origin_y;
        for (size_t k_x=0; k_x<n_nodes_per_axis; k_x++) {
            const double x = unit_length*k_x + origin_x;
            GeoLib::Point p(x, y, z);
            msh->setNodeCoordinates(node_id++, p);
            //msh->setNodeCoordinates(node_id++, GeoLib::Point(x, y, z));
        }
    }

    //elements
    //size_t ele_id(0);
    for (size_t j=0; j<subdivision; j++) {
        const size_t offset_y1 = j*n_nodes_per_axis;
        const size_t offset_y2 = (j+1)*n_nodes_per_axis;
        for (size_t k=0; k<subdivision; k++) {
            Quadrirateral *e = new Quadrirateral();
            e->setNodeID(0, offset_y1+k);
            e->setNodeID(1, offset_y1+k+1);
            e->setNodeID(2, offset_y2+k+1);
            e->setNodeID(3, offset_y2+k);
            //for (size_t l=0; l<4; l++)
            //    e->setNode(l, msh->getNode(e->getNodeID(l)));
            msh->addElement(e);
        }
    }

    return msh;

};

StructuredMesh<ElementShape::QUAD>* MeshGenerator::generateStructuredRegularQuadMesh(const double length, const size_t subdivision, const double origin_x, const double origin_y, const double origin_z) 
{
    GeoLib::Point org(origin_x, origin_y, origin_z);
    GeoLib::Point pt_length(length, length, .0);
    size_t unit_div[3] = {subdivision, subdivision, 1};

    return new StructuredMesh<ElementShape::QUAD>(MeshLib::CoordinateSystemType::XY, org, pt_length.getData(), unit_div);

};

void MeshGenerator::generateSubMesh(const MeshLib::IMesh &src, const std::vector<size_t> &list_e, MeshLib::IMesh* &dest, BaseLib::BidirectionalMap<size_t, size_t> &map_global2local)
{
    MeshLib::UnstructuredMesh* new_msh = new MeshLib::UnstructuredMesh(src.getGeometricProperty()->getCoordinateSystem().getType());
    std::set<size_t> list_nodes_subset;
    for (size_t i=0; i<list_e.size(); i++) {
        const MeshLib::IElement *e = src.getElement(list_e[i]);
        for (size_t j=0; j<e->getNumberOfNodes(); j++) {
            list_nodes_subset.insert(e->getNodeID(j));
        }
    }

    for (std::set<size_t>::iterator itr=list_nodes_subset.begin(); itr!=list_nodes_subset.end(); ++itr) {
        GeoLib::Point p(src.getNodeCoordinates(*itr));
        size_t new_node_id = new_msh->addNode(p);
        //size_t new_node_id = new_msh->addNode(src.getNodeCoordinates(*itr));
        map_global2local.insert(*itr, new_node_id);
    }

    for (size_t i=0; i<list_e.size(); i++) {
        const MeshLib::IElement *e = src.getElement(list_e[i]);
        MeshLib::IElement *new_e = e->clone();
        for (size_t j=0; j<e->getNumberOfNodes(); j++) {
            new_e->setNodeID(j, map_global2local.mapAtoB(e->getNodeID(j)));
        }
        new_msh->addElement(new_e);
    }
    dest = new_msh;
}

/**
 * generate higher order mesh
 */
void MeshGenerator::generateHigherOrderUnstrucuredMesh(UnstructuredMesh &msh, size_t order)
{
    assert(order<3);

    TopologySequentialNodes2Elements nod2ele(msh);

    if (msh.getNumberOfEdges()==0) {
        createEdgeElements(&msh);
    }

    // make all edges higher order
    for (size_t i=0; i<msh.getNumberOfEdges(); i++) {
        IElement* e = msh.getEdgeElement(i);
        e->setMaximumOrder(order);
        if (order==2) {
            double const* const pnt0(msh.getNodeCoordinatesRef(e->getNodeID(0))->getData());
            double const* const pnt1(msh.getNodeCoordinatesRef(e->getNodeID(1))->getData());
            // add midpoint
            GeoLib::Point p(.5 * (pnt0[0] + pnt1[0]), 0.5 * (pnt0[1] + pnt1[1]), 0.5 * (pnt0[2] + pnt1[2]));
            size_t nod_id = msh.addNode(p, 2);
            //size_t nod_id = msh.addNode(GeoLib::Point(.5 * (pnt0[0] + pnt1[0]), 0.5 * (pnt0[1] + pnt1[1]), 0.5 * (pnt0[2] + pnt1[2])), 2);
            e->setNodeID(2, nod_id);
        } else {
            //
            std::cout << "***Error: order (" << order << ") is not supported in generateHigherOrderUnstrucuredMesh()" << std::endl;
        }
    }

    // set the new node ids to all elements
    for (size_t i=0; i<msh.getNumberOfElements(); i++) {
        IElement* e = msh.getElement(i);
        size_t e_nnodes1 = e->getNumberOfNodes(1);
        e->setMaximumOrder(order);
        // for each edge
        for (size_t j=0; j<e->getNumberOfEdges(); j++) {
            IElement *edge = e->getEdge(j);
            size_t new_nod = edge->getNodeID(2);
            size_t local_id = e_nnodes1 + j;
            e->setNodeID(local_id, new_nod);
        }

        // Quad 9
        if (e->getShapeType() == ElementShape::QUAD || e->getNumberOfNodes(2)==9)
        {
            double x0, y0, z0;
            x0 = y0 = z0 = 0.0;
            std::vector<size_t> list_e_nodes;
            e->getNodeIDList(1, list_e_nodes);
            for (size_t i = 0; i < list_e_nodes.size(); i++) // Nodes
            {
                const GeoLib::Point* pt = msh.getNodeCoordinatesRef(list_e_nodes[i]);
                double const* const pnt_i(pt->getData());
                x0 += pnt_i[0];
                y0 += pnt_i[1];
                z0 += pnt_i[2];
            }
            x0 /= (double) e_nnodes1;
            y0 /= (double) e_nnodes1;
            z0 /= (double) e_nnodes1;
            GeoLib::Point p(x0,y0,z0);
            size_t nodid = msh.addNode(p,2);
            //size_t nodid = msh.addNode(GeoLib::Point(x0,y0,z0),2);
            e->setNodeID(8, nodid);
        }
    }
}

}
