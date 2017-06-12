/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IElement.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma  once

#include <vector>

#include "GeoLib/Point.h"

#include "MeshLib/Core/Node.h"

namespace MeshLib
{

class IElementCoordinatesMapping;

/// List of element shape types
struct ElementShape
{
    enum type {
        LINE = 1,
        QUAD = 2,
        HEXAHEDRON = 3,
        TRIANGLE = 4,
        TETRAHEDRON = 5,
        PRISM = 6,
        PYRAMID = 7,
        INVALID = -1
    };
};

/**
 * \brief Interface of element classes
 *
 * Any elements have the following attributes
 * - element id
 * - group id
 * - shape type
 * - dimension size
 * - node index
 * - faces
 * - edges
 * - coordinates mapping
 * - order
 */
class IElement
{
public:
    ///
    virtual ~IElement() {};
    
    //--- init/finalize ---
    /// reset element data 
    virtual void reset() = 0;
    /// return a clone of this object
    virtual IElement* clone() const = 0;

    //--- element topology ---
    /// return the shape type
    virtual ElementShape::type getShapeType() const = 0;
    /// return intrinsic dimensions of this element
    virtual size_t getDimension() const = 0;
    /// return the number of nodes under current order
    virtual size_t getNumberOfNodes() const = 0;
    /// get the number of element faces
    virtual size_t getNumberOfFaces() const = 0;
    /// get the number of element edges
    virtual size_t getNumberOfEdges() const = 0;

    //--- element attributes ---
    /// return this element id
    virtual size_t getID() const = 0;
    /// set this element id
    virtual void setID(size_t id) = 0;

    /// return the group id of this element
    virtual size_t getGroupID() const = 0;
    /// set the group if of this element
    virtual void setGroupID(size_t id) = 0;

    /// get mapped coordinates
    virtual IElementCoordinatesMapping* getMappedCoordinates() const = 0;
    /// set mapped coordinates
    virtual void setMappedCoordinates(IElementCoordinatesMapping* mapping) = 0;

    /// return node id
    virtual size_t getNodeID(size_t local_node_id) const = 0;
    /// set node id
    virtual void setNodeID(const size_t &local_node_id, const size_t &node_id) = 0;
    /// get a list of node ids
    virtual void getNodeIDList( std::vector<size_t> &e_node_id_list ) const = 0;
    /// return if this element has the given list of nodes 
    virtual bool hasNodeIds(const std::vector<size_t> &node_ids) const = 0;

    /// get an edge element
    virtual IElement* getEdge(size_t edge_id) const = 0;
    /// set an edge element
    virtual void setEdge(size_t edge_id, IElement* e) = 0;
    /// get node ids of the edge element
    virtual void getNodeIDsOfEdges(size_t edge_id, std::vector<size_t> &vec_node_ids) const = 0;
    /// get edge element type
    virtual ElementShape::type getEdgeType(size_t edge_id) const = 0;

    //--- element order ---
    /// set max. order
    virtual void setMaximumOrder(size_t) = 0;
    /// get max. order
    virtual size_t getMaximumOrder() const = 0;
    /// set current order
    virtual void setCurrentOrder(size_t) const = 0;
    /// get current order
    virtual size_t getCurrentOrder() const = 0;
    /// return the number of nodes under the given order
    virtual size_t getNumberOfNodes(size_t order) const = 0;
    /// get a list of node id
    virtual void getNodeIDList( size_t order, std::vector<size_t> &e_node_id_list ) const = 0;
    /// get a list of node id of the specified edge
    virtual void getNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_node_ids) const = 0;
    /// get a list of node size with orders
    virtual void getListOfNumberOfNodesForAllOrders(std::vector<size_t> &vec) const = 0;
};

}
