/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IMesh.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "GeoLib/Point.h"
#include "IElement.h"

namespace MeshLib
{
class CoordinateSystem;
class Node;
class MeshGeometricProperty;

/**
 * \brief Interface to mesh classes
 *
 */
class IMeshBase
{
public:
    ///
    virtual ~IMeshBase() {};

    //--- mesh properties ---
    /// set mesh id
    virtual void setID(size_t id) = 0;
    /// return mesh id
    virtual size_t getID() const = 0;
    /// return mesh dimension
    virtual size_t getDimension() const = 0;
    /// get mesh geometric property
    virtual const MeshGeometricProperty* getGeometricProperty() const = 0;
    /// return if this mesh is axisymmetric or not
    virtual bool isAxisymmetric() const = 0;
    /// set axisymmetric
    virtual void setAxisymmetric(bool flag) = 0;

    //--- nodes ---
    /// get the number of nodes
    virtual size_t getNumberOfNodes() const = 0;
    /// get a node object
    virtual Node* getNode( size_t id ) const = 0;
    /// get a point object
    virtual const GeoLib::Point* getNodeCoordinatesRef(size_t id) const = 0;
    /// get a point object
    virtual GeoLib::Point getNodeCoordinates(size_t id) const = 0;
    /// get a list of points in the given element
    virtual void getListOfNodeCoordinates(const std::vector<size_t> &vec_node_id, std::vector<GeoLib::Point> &vec_pt) const = 0;

    //--- elements ---
    /// get the number of elements
    virtual size_t getNumberOfElements() const = 0;
    /// get the number of elements of the given type
    virtual size_t getNumberOfElements(ElementShape::type ele_type) const = 0;
    /// get an element
    virtual IElement* getElement( size_t element_id ) const = 0;

    //--- edges ---
    /// get the number of edges
    virtual size_t getNumberOfEdges() const = 0;
    /// add a new element
    virtual void addEdgeElement(IElement*) = 0;
    /// get a edge element
    virtual IElement* getEdgeElement(size_t edge_id) = 0;

    //--- faces ---
    
    //--- setup ---
    /// construct geometric data of this mesh
    virtual void constructGeometricProperty() = 0;
};

/**
 * \brief Interface for mixed order mesh
 */
class IMesh : public IMeshBase
{
public:
    ///
    virtual ~IMesh() {};
    
    /// get the maximum possible order
    virtual size_t getMaxiumOrder() const = 0;
    
    /// set current order
    virtual void setCurrentOrder(size_t order) const = 0;
    
    /// get current order
    virtual size_t getCurrentOrder() const = 0;
    
    /// get the number of nodes with the current order
    virtual size_t getNumberOfNodes() const = 0; // i don't understand why this redefinition is needed

    /// get the number of nodes with the given order
    virtual size_t getNumberOfNodes(size_t order) const = 0;

    /// get the number of all possible nodes 
    size_t getNumberOfTotalNodes() const
    {
        return getNumberOfNodes(getMaxiumOrder());
    };
};

}
