/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file StructuredMesh.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <map>
#include <algorithm>

#include "BaseLib/CodingTools.h"
#include "GeoLib/Point.h"

#include "IMesh.h"
#include "IElement.h"
#include "ElementFactory.h"
#include "MeshGeometricProperty.h"

namespace MeshLib
{

/**
 * \brief a structured mesh class
 * 
 * \tparam T_ELEMENT_TYPE element type
 */
template <ElementShape::type T_ELEMENT_TYPE>
class StructuredMesh : public IMesh
{
public:

    ///
    StructuredMesh(CoordinateSystemType::type coord, GeoLib::Point &org_pt, const double* len, const size_t* division)
    : _msh_id(0), _origin(org_pt), _e(0), _nod(0)
    {
        _order = 1;
        _isAxisymmetric = false;
        _geo_prop.setCoordinateSystem(coord);
        const size_t dim = getDimension();
        for (size_t i=0; i<dim; i++)
            _length[i] = len[i];
        for (size_t i=0; i<dim; i++)
            _number_of_elements_per_dimension[i] = division[i];
        for (size_t i=0; i<dim; i++)
            _spacing[i] = _length[i] / division[i];

        double d = _spacing[0];
        for (size_t i=1; i<getDimension(); i++)
            d = std::min(d, _spacing[i]);
        _geo_prop.setMinEdgeLength(d);

        construct();
    }

    ///
    virtual ~StructuredMesh()
    {
        BaseLib::releaseObject(_e, _nod);
        BaseLib::releaseObjectsInStdVector(_list_edge_elements);
    }

    //--- attributes ---
    /// 
    void setID(size_t id) {_msh_id = id;};
    virtual size_t getID() const {return _msh_id;};

    /// 
    size_t getDimension() const 
    {
        return _geo_prop.getCoordinateSystem().getDimension();
    }

    /// 
    const MeshGeometricProperty* getGeometricProperty() const { return &_geo_prop; }

    ///
    void setAxisymmetric(bool flag) { _isAxisymmetric = flag; }
    bool isAxisymmetric() const { return _isAxisymmetric; }

    /// 
    void setMaxiumOrder(size_t order)
    {
        assert(order<3);

        size_t n_edges = _number_of_elements_per_dimension[0]*(_number_of_elements_per_dimension[1]+1) + _number_of_elements_per_dimension[1]*(_number_of_elements_per_dimension[0]+1);

        if (order==2) {
            _map_order_nnodes[2] = _map_order_nnodes[1] + n_edges + _n_ele;
        }
    }
    size_t getMaxiumOrder() const { return _map_order_nnodes.size(); }

    ///
    void setCurrentOrder(size_t order) const { _order = order; }
    size_t getCurrentOrder() const { return _order; }
    
    
    //--- nodes ---
    /// get the number of nodes
    size_t getNumberOfNodes() const { return getNumberOfNodes(_order); }
    size_t getNumberOfNodes(size_t order) const
    {
        if (order-1 < _map_order_nnodes.size()) {
            std::map<size_t, size_t>::const_iterator itr = _map_order_nnodes.find(order);
            return itr->second;
        }
        return 0;
    }
    /// get a node object
    Node* getNode( size_t id ) const;
    /// get a point object
    const GeoLib::Point* getNodeCoordinatesRef(size_t id) const 
    {
        getNodeCoordinates(id, const_cast<GeoLib::Point*>(_nod->getX()));
        return _nod->getX();
    }
    /// get a point object
    GeoLib::Point getNodeCoordinates(size_t id) const 
    {
        return *getNodeCoordinatesRef(id);
    }
    ///
    void getListOfNodeCoordinates(const std::vector<size_t> &vec_node_id, std::vector<GeoLib::Point> &vec_pt) const 
    {
        for (size_t i=0; i<vec_node_id.size(); i++)
            vec_pt.push_back(*this->getNodeCoordinatesRef(vec_node_id[i]));
    }
    
   
    //--- elements ---
    /// get the number of elements
    size_t getNumberOfElements() const { return _n_ele; }
    /// get the number of elements of the given type
    size_t getNumberOfElements(ElementShape::type ele_type) const {
        if (ele_type==T_ELEMENT_TYPE) return _n_ele;
        else return 0;
    }
    /// get an element
    IElement* getElement( size_t element_id ) const;

    //--- edges ---
    /// add a new element
    void addEdgeElement(IElement* e) { _list_edge_elements.push_back(e); }

    size_t getNumberOfEdges() const { return _list_edge_elements.size(); }

    IElement* getEdgeElement(size_t edge_id) {
        if (edge_id<_list_edge_elements.size())
            return _list_edge_elements[edge_id];
        return 0;
    }

    //--- faces ---

    // --- setup ---
    virtual void constructGeometricProperty() {}
    

private:
    void construct();
    void getNodeCoordinates(size_t id, GeoLib::Point* p) const;
    
private:
    size_t _msh_id;
    GeoLib::Point _origin;
    double _length[3];
    double  _spacing[3];
    size_t  _number_of_elements_per_dimension[3];
    size_t  _number_of_nodes_per_dimension[3];
    std::map<size_t, size_t> _map_order_nnodes;
    IElement* _e;
    Node* _nod;
    size_t _n_ele;
    std::vector<IElement*> _list_edge_elements;
    MeshGeometricProperty _geo_prop;
    bool _isAxisymmetric;
    mutable size_t _order;
};

template<> void StructuredMesh<ElementShape::QUAD>::construct();
template<> IElement* StructuredMesh<ElementShape::QUAD>::getElement( size_t element_id ) const;
template<> Node* StructuredMesh<ElementShape::QUAD>::getNode( size_t id ) const; 
template<> void StructuredMesh<ElementShape::QUAD>::getNodeCoordinates(size_t id,  GeoLib::Point* p) const;

}
