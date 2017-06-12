/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TemplateElementTopology.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */


#pragma once

#include "IElement.h"

namespace MeshLib
{

/**
 * \brief Template class for element topological rule
 * 
 * \tparam TYPE             element shape type
 * \tparam DIMENSION        element dimension
 * \tparam NUMBER_OF_FACES  the number of faces
 * \tparam NUMER_OF_EDGES   the number of edges
 * \tparam NUMBER_OF_NODES1 the number of nodes of the 1st order element
 * \tparam NUMBER_OF_NODES2 the number of nodes of the 2nd order element
 */
template <
    ElementShape::type TYPE, 
    size_t DIMENSION, 
    size_t NUMBER_OF_FACES, 
    size_t NUMER_OF_EDGES, 
    size_t NUMBER_OF_NODES1, 
    size_t NUMBER_OF_NODES2
    >
class TemplateElementTopology
{
public:
    /// return shape type
    static ElementShape::type getShapeType() {return TYPE;};
    
    /// return element dimension
    static size_t getDimension() {return DIMENSION;};

    /// return the number of element faces
    static size_t getNumberOfFaces() {return NUMBER_OF_FACES;};

    /// return the number of element edges
    static size_t getNumberOfEdges() {return NUMER_OF_EDGES;};
    
    /// return the number of element nodes
    static size_t getNumberOfNodes(size_t order)
    {
        switch (order) {
        case 1: return NUMBER_OF_NODES1;
        case 2: return NUMBER_OF_NODES2;
        }
        return 0;
    }
};

} //end namespace

