/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file UnstructuredMesh.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <map>

#include "IMesh.h"
#include "IElement.h"
#include "Node.h"
#include "MeshGeometricProperty.h"


namespace MeshLib
{

/**
 * \brief Unstructured mesh class
 *
 *
 */
class UnstructuredMesh : public IMesh
{
public:
    UnstructuredMesh() : _msh_id(0), _isAxisymmetric(false), _order(1)
    {
    };

    explicit UnstructuredMesh(const CoordinateSystemType::type coord)
    : _msh_id(0), _isAxisymmetric(false), _order(1)
    {
        _geo_prop.setCoordinateSystem(coord);
    };

    virtual ~UnstructuredMesh();

    //--- attributes ---
    void setID(size_t id) {_msh_id = id;};
    virtual size_t getID() const {return _msh_id;};

    size_t getDimension() const
    {
        return _geo_prop.getCoordinateSystem().getDimension();
    };

    /// get mesh geometric property
    const MeshGeometricProperty* getGeometricProperty() const { return &_geo_prop; }
    MeshGeometricProperty* getGeometricProperty() { return &_geo_prop; }

    ///
    bool isAxisymmetric() const { return _isAxisymmetric; }
    void setAxisymmetric(bool flag) { _isAxisymmetric = flag; }

    void setCurrentOrder(size_t order) const { _order = order; };
    size_t getCurrentOrder() const { return _order; };

    size_t getMaxiumOrder() const { return _map_order_nnodes.size(); }


    //--- nodes ---
    size_t getNumberOfNodes(size_t order) const
    {
        if (order-1 < _map_order_nnodes.size()) {
            size_t sum = 0;
            for (size_t i=0; i<order; ++i) {
                std::map<size_t, size_t>::const_iterator itr = _map_order_nnodes.find(i+1);
                sum += itr->second;
            }
            return sum;
        }
        return 0;
    }
    
    size_t getNumberOfNodes() const { return getNumberOfNodes(_order); }

    Node* getNode( size_t id ) const { return _list_nodes[id]; };

    const std::vector<Node*>& getNodeVector() const { return _list_nodes; };

    size_t addNode( GeoLib::Point &x, size_t order=1 ); 

    size_t setNodeCoordinates( size_t node_id, GeoLib::Point &x );

    const GeoLib::Point*  getNodeCoordinatesRef(size_t id) const
    {
        return getNode(id)->getX();
    }

    GeoLib::Point getNodeCoordinates(size_t id) const
    {
        return *getNode(id)->getX();
    }

    void getListOfNodeCoordinates(const std::vector<size_t> &vec_node_id, std::vector<GeoLib::Point> &vec_pt) const;


    //--- elements ---
    ///
    size_t getNumberOfElements() const { return _list_elements.size(); };

    size_t getNumberOfElements(ElementShape::type ele_type) const;

    IElement* getElement( size_t element_id ) const 
    {
        assert(element_id<_list_elements.size());
        IElement* e = _list_elements[element_id];
        e->setCurrentOrder(_order);
        return e;
    };

    size_t addElement( IElement *e);

    //--- edges ---
    void addEdgeElement(IElement *e);
    size_t getNumberOfEdges() const;
    IElement* getEdgeElement(size_t edge_id);

    //--- setup ---
    /// construct geometric data of this mesh
    virtual void constructGeometricProperty();

private:
    size_t _msh_id;
    std::vector<Node*> _list_nodes;
    std::vector<IElement*> _list_elements;
    std::map<ElementShape::type, size_t> _n_eles_type;
    std::vector<IElement*> _list_edge_elements;
    MeshGeometricProperty _geo_prop;
    bool _isAxisymmetric;
    mutable size_t _order;
    std::map<size_t, size_t> _map_order_nnodes;
};

}
