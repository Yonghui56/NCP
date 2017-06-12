/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementTopology.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */


#pragma once

#include "IElement.h"
#include "TemplateElementTopology.h"

namespace MeshLib
{

class LineTopology : public TemplateElementTopology<ElementShape::LINE, 1, 0, 1, 2, 3>
{
public:
    static ElementShape::type getEdgeElementType(size_t) { return ElementShape::LINE; }
    static void getLocalNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_local_node_ids);
};

class TriangleTopology : public TemplateElementTopology<ElementShape::TRIANGLE, 2, 0, 3, 3, 6>
{
public:
    static ElementShape::type getEdgeElementType(size_t) { return ElementShape::LINE; }
    static void getLocalNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_local_node_ids);
};

class Quad8Topology : public TemplateElementTopology<ElementShape::QUAD, 2, 0, 4, 4, 8>
{
public:
    static ElementShape::type getEdgeElementType(size_t) { return ElementShape::LINE; }
    static void getLocalNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_local_node_ids);
};

class Quad9Topology : public TemplateElementTopology<ElementShape::QUAD, 2, 0, 4, 4, 9>
{
public:
    static ElementShape::type getEdgeElementType(size_t) { return ElementShape::LINE; }
    static void getLocalNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_local_node_ids);
};

class TetraTopology : public TemplateElementTopology<ElementShape::TETRAHEDRON, 3, 4, 6, 4, 10>
{
public:
    static ElementShape::type getEdgeElementType(size_t) { return ElementShape::LINE; }
    static void getLocalNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_local_node_ids);
    static ElementShape::type getFaceElementType(size_t) { return ElementShape::TRIANGLE; }
    static void getLocalNodeIDsOfFaceElement(size_t order, size_t face_id, std::vector<size_t> &vec_local_node_ids);
};

class HexTopology : public TemplateElementTopology<ElementShape::HEXAHEDRON, 3, 6, 12, 8, 20>
{
public:
    static ElementShape::type getEdgeElementType(size_t) { return ElementShape::LINE; }
    static void getLocalNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_local_node_ids);
    static ElementShape::type getFaceElementType(size_t) { return ElementShape::QUAD; }
    static void getLocalNodeIDsOfFaceElement(size_t order, size_t face_id, std::vector<size_t> &vec_local_node_ids);
};

class PrismTopology : public TemplateElementTopology<ElementShape::PRISM, 3, 5, 9, 6, 15>
{
public:
    static ElementShape::type getEdgeElementType(size_t) { return ElementShape::LINE; }
    static void getLocalNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_local_node_ids);
    static ElementShape::type getFaceElementType(size_t face_id)
    { 
        if (face_id<2) return ElementShape::TRIANGLE;
        else return ElementShape::QUAD;
    }
    static void getLocalNodeIDsOfFaceElement(size_t order, size_t face_id, std::vector<size_t> &vec_local_node_ids);
};

class PyramidTopology : public TemplateElementTopology<ElementShape::PYRAMID, 3, 5, 8, 5, 13>
{
public:
    static ElementShape::type getEdgeElementType(size_t) { return ElementShape::LINE; }
    static void getLocalNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_local_node_ids);
    static ElementShape::type getFaceElementType(size_t face_id)
    { 
        if (face_id==0) return ElementShape::QUAD;
        else return ElementShape::TRIANGLE;
    }
    static void getLocalNodeIDsOfFaceElement(size_t order, size_t face_id, std::vector<size_t> &vec_local_node_ids);
};

} //end namespace

