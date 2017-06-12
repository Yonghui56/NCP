/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file StructuredMesh.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "StructuredMesh.h"
#include "Element.h"

namespace MeshLib
{

template<> void StructuredMesh<ElementShape::QUAD>::construct()
{
    const size_t dim = this->getDimension();

    _n_ele = 1;
    for (size_t i=0; i<dim; i++)
        _n_ele *= _number_of_elements_per_dimension[i];

    for (size_t i=0; i<dim; i++)
        _number_of_nodes_per_dimension[i] = _number_of_elements_per_dimension[i] + 1;

    size_t _nod_size = 1;
    for (size_t i=0; i<dim; i++)
        _nod_size *= _number_of_nodes_per_dimension[i];
    _map_order_nnodes[1] = _nod_size;

    _e = new Quadrirateral();
    _nod = new Node();
}


template<> IElement* StructuredMesh<ElementShape::QUAD>::getElement( size_t element_id ) const
{
    //set e
    _e->reset();
    _e->setID(element_id);
    _e->setMaximumOrder(getMaxiumOrder());
    _e->setCurrentOrder(_order);
    const size_t x_j = element_id / _number_of_elements_per_dimension[0];
    const size_t offset_y1 = x_j*_number_of_nodes_per_dimension[0];
    const size_t offset_y2 = (x_j+1)*_number_of_nodes_per_dimension[0];
    const size_t k = element_id % _number_of_elements_per_dimension[0];
    for (size_t i=0; i<4; i++) {
        _e->setNodeID(0, offset_y1+k);
        _e->setNodeID(1, offset_y1+k+1);
        _e->setNodeID(2, offset_y2+k+1);
        _e->setNodeID(3, offset_y2+k);
    }
    if (getMaxiumOrder()==2) {
        const size_t nnodes1 = getNumberOfNodes(1);
        const size_t nnodes2 = getNumberOfNodes(2);
        const size_t offset2_y1 = nnodes1 + x_j*(_number_of_elements_per_dimension[0] + _number_of_nodes_per_dimension[0]);
        const size_t offset2_y2 = offset2_y1 + _number_of_elements_per_dimension[0];
        const size_t offset2_y3 = offset2_y2 + _number_of_nodes_per_dimension[0];
        _e->setNodeID(4, offset2_y1+k);
        _e->setNodeID(5, offset2_y2+k+1);
        _e->setNodeID(6, offset2_y3+k);
        _e->setNodeID(7, offset2_y2+k);
        _e->setNodeID(8, nnodes2-_n_ele+element_id);
    }
    return _e;
};

template<> Node* StructuredMesh<ElementShape::QUAD>::getNode( size_t id ) const 
{
    _nod->setNodeID(id);
    getNodeCoordinates(id, const_cast<GeoLib::Point*>(_nod->getX()));

    return _nod;
};

template<> void StructuredMesh<ElementShape::QUAD>::getNodeCoordinates(size_t id,  GeoLib::Point* p) const
{
    double *pt = const_cast<double*>(p->getData());
    const size_t nnodes1 = getNumberOfNodes(1);
    if (id < nnodes1) {
        size_t k_x = id % _number_of_nodes_per_dimension[0]; 
        size_t j_y = id / _number_of_nodes_per_dimension[1]; 
        pt[0] = _spacing[0]*k_x + _origin[0];
        pt[1] = _spacing[1]*j_y + _origin[1];
    } else {
        const size_t nnodes2 = getNumberOfNodes(2);
        if (id < nnodes2 - _n_ele) {
            const size_t id1 = id - nnodes1;
            const size_t n_row_block = _number_of_nodes_per_dimension[0] + _number_of_elements_per_dimension[0];
            size_t block_id = id1 / n_row_block;
            size_t local_id = id1 % n_row_block;
            bool isRow1 = (0==local_id / _number_of_elements_per_dimension[0]);
            if (isRow1) {
                pt[0] = _spacing[0]*local_id + 0.5*_spacing[0] + _origin[0];
                pt[1] = _spacing[1]*block_id + _origin[1];
            } else {
                pt[0] = _spacing[0]*(local_id-_number_of_elements_per_dimension[0]) + _origin[0];
                pt[1] = _spacing[1]*block_id + 0.5*_spacing[1] + _origin[1];
            }
        } else {
            const size_t e_id = id - (nnodes2 - _n_ele);
            const size_t j_y = e_id / _number_of_elements_per_dimension[0];
            const size_t k_x = e_id % _number_of_elements_per_dimension[0];
            pt[0] = _spacing[0]*k_x + 0.5*_spacing[0] + _origin[0];
            pt[1] = _spacing[1]*j_y + 0.5*_spacing[1] + _origin[1];
        }

    }
    pt[2] = _origin[2];
}

}
