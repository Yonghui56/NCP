/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DecomposedDomain.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MeshLib/Core/IMesh.h"
#include "DecomposedDomain.h"
#include "SubDomain.h"

namespace DiscreteLib
{

/**
 * \brief 
 * 
 */
class DecomposedMesh : public MeshLib::IMesh
{
public:
    /**
     * 
     * @param id
     * @param ddc_dom
     */
    DecomposedMesh(size_t id, DecomposedDomain* ddc_dom)
    : _id(id), _ddc_dom(ddc_dom), _msh0(NULL)
    {
        if (ddc_dom->getNumberOfSubDomains()>0) {
            _msh0 = ddc_dom->getSubDomain(0)->getLoalMesh();
        }
    }

    /**
     * 
     */
    virtual ~DecomposedMesh()
    {
    }

    //--- mesh properties ---
    /// set mesh id
    virtual void setID(size_t id) {_id = id;};
    /// return mesh id
    virtual size_t getID() const { return _id;};
    /// return mesh dimension
    virtual size_t getDimension() const { return _msh0->getDimension();};
    /// get mesh geometric property
    virtual const MeshLib::MeshGeometricProperty* getGeometricProperty() const
    {
        return _msh0->getGeometricProperty();
    };
    /// return if this mesh is axisymmetric or not
    virtual bool isAxisymmetric() const { return _msh0->isAxisymmetric();};
    /// set axisymmetric
    virtual void setAxisymmetric(bool flag) { _msh0->setAxisymmetric(flag);};

    //--- nodes ---
    /// get the number of nodes
    virtual size_t getNumberOfNodes() const {return _ddc_dom->getNumberOfNodes();};
    /// get a node object
    virtual MeshLib::Node* getNode( size_t id ) const  {return _msh0->getNode(id);};
    /// get a point object
    virtual const GeoLib::Point* getNodeCoordinatesRef(size_t id) const {return _msh0->getNodeCoordinatesRef(id);};
    /// get a point object
    virtual GeoLib::Point getNodeCoordinates(size_t id) const {return _msh0->getNodeCoordinates(id);};
    /// get a list of points in the given element
    virtual void getListOfNodeCoordinates(const std::vector<size_t> &vec_node_id, std::vector<GeoLib::Point> &vec_pt) const
    {
        _msh0->getListOfNodeCoordinates(vec_node_id, vec_pt);
    }

    //--- elements ---
    /// get the number of elements
    virtual size_t getNumberOfElements() const {return _ddc_dom->getNumberOfElements();};
    /// get the number of elements of the given type
    virtual size_t getNumberOfElements(MeshLib::ElementShape::type ele_type) const {return _msh0->getNumberOfElements(ele_type);};
    /// get an element
    virtual MeshLib::IElement* getElement( size_t element_id ) const {return _msh0->getElement(element_id);};

    //--- edges ---
    /// get the number of edges
    virtual size_t getNumberOfEdges() const {return _msh0->getNumberOfNodes();};
    /// add a new element
    virtual void addEdgeElement(MeshLib::IElement*) {};
    /// get a edge element
    virtual MeshLib::IElement* getEdgeElement(size_t edge_id) {return _msh0->getEdgeElement(edge_id);};

    //--- faces ---

    //--- setup ---
    /// construct geometric data of this mesh
    virtual void constructGeometricProperty()
    {
        _msh0->constructGeometricProperty();
    };
    
    ///
    DecomposedDomain* getDecomposedDomain() const {return _ddc_dom;};

    /// get the maximum possible order
    virtual size_t getMaxiumOrder() const {return _msh0->getMaxiumOrder();};

    /// set current order
    virtual void setCurrentOrder(size_t order) const {return _msh0->setCurrentOrder(order);};

    /// get current order
    virtual size_t getCurrentOrder() const {return _msh0->getCurrentOrder();};

    /// get the number of nodes with the given order
    virtual size_t getNumberOfNodes(size_t order) const {return _msh0->getNumberOfNodes(order);};

private:
    size_t _id;
    DecomposedDomain* _ddc_dom;
    MeshLib::IMesh* _msh0;
};


};
