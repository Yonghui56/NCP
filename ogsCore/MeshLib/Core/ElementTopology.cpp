/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementTopology.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "ElementTopology.h"

#include <vector>

namespace MeshLib
{

void LineTopology::getLocalNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_local_node_ids)
{
    vec_local_node_ids.resize(2); 
    vec_local_node_ids[0] = edge_id;
    vec_local_node_ids[1] = edge_id+1;
    if (order==2)
        vec_local_node_ids[2] = edge_id+2;
}

void TriangleTopology::getLocalNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_local_node_ids)
{
    size_t n_nodes = LineTopology::getNumberOfNodes(order);
    vec_local_node_ids.resize(n_nodes); 
    vec_local_node_ids[0] = edge_id;
    vec_local_node_ids[1] = (edge_id+1)%3;
    if (order==2)
        vec_local_node_ids[2] = edge_id+3;
}

void Quad8Topology::getLocalNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_local_node_ids)
{
    size_t n_nodes = LineTopology::getNumberOfNodes(order);
    vec_local_node_ids.resize(n_nodes); 
    vec_local_node_ids[0] = edge_id;
    vec_local_node_ids[1] = (edge_id+1)%4;
    if (order==2)
        vec_local_node_ids[2] = edge_id+4;
}

void Quad9Topology::getLocalNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_local_node_ids)
{
    size_t n_nodes = LineTopology::getNumberOfNodes(order);
    vec_local_node_ids.resize(n_nodes); 
    vec_local_node_ids[0] = edge_id;
    vec_local_node_ids[1] = (edge_id+1)%4;
    if (order==2)
        vec_local_node_ids[2] = edge_id+4;
}

void TetraTopology::getLocalNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_local_node_ids)
{
    size_t n_nodes = LineTopology::getNumberOfNodes(order);
    vec_local_node_ids.resize(n_nodes); 
    if (edge_id<3) {
        vec_local_node_ids[0] = edge_id;
        vec_local_node_ids[1] = (edge_id+1)%3;
    } else {
        vec_local_node_ids[0] = 3;
        vec_local_node_ids[1] = (edge_id+1)%3;
    }
    if (order==2)
        vec_local_node_ids[2] = edge_id+4;
}

void TetraTopology::getLocalNodeIDsOfFaceElement(size_t order, size_t face_id, std::vector<size_t> &vec_local_node_ids)
{
    size_t n_nodes = TriangleTopology::getNumberOfNodes(order);
    vec_local_node_ids.resize(n_nodes); 
    switch(face_id)
    {
    case 0:
        vec_local_node_ids[0] = 1;
        vec_local_node_ids[1] = 2;
        vec_local_node_ids[2] = 3;
        if(order==2)
        {
            vec_local_node_ids[3] = 5;
            vec_local_node_ids[4] = 8;
            vec_local_node_ids[5] = 7;
        }
        break;
    case 1:
        vec_local_node_ids[0] = 3;
        vec_local_node_ids[1] = 2;
        vec_local_node_ids[2] = 0;
        if(order==2)
        {
            vec_local_node_ids[3] = 8;
            vec_local_node_ids[4] = 6;
            vec_local_node_ids[5] = 9;
        }
        break;
    case 2:
        vec_local_node_ids[0] = 1;
        vec_local_node_ids[1] = 3;
        vec_local_node_ids[2] = 0;
        if(order==2)
        {
            vec_local_node_ids[3] = 7;
            vec_local_node_ids[4] = 9;
            vec_local_node_ids[5] = 4;
        }
        break;
    case 3:
        vec_local_node_ids[0] = 0;
        vec_local_node_ids[1] = 2;
        vec_local_node_ids[2] = 1;
        if(order==2)
        {
            vec_local_node_ids[3] = 6;
            vec_local_node_ids[4] = 5;
            vec_local_node_ids[5] = 4;
        }
        break;
    }
}

void HexTopology::getLocalNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_local_node_ids)
{
    size_t n_nodes = LineTopology::getNumberOfNodes(order);
    vec_local_node_ids.resize(n_nodes); 
    if (edge_id<8) {
        vec_local_node_ids[0] = edge_id;
        vec_local_node_ids[1] = (edge_id + 1) % 4 + 4 * static_cast<int>(edge_id / 4);
    } else {
        vec_local_node_ids[0] = edge_id % 4;
        vec_local_node_ids[1] = edge_id % 4 + 4;
    }
    if (order==2)
        vec_local_node_ids[2] = edge_id+8;
}

void HexTopology::getLocalNodeIDsOfFaceElement(size_t order, size_t face_id, std::vector<size_t> &vec_local_node_ids)
{
    size_t n_nodes = Quad8Topology::getNumberOfNodes(order);
    vec_local_node_ids.resize(n_nodes); 
    switch(face_id)
    {
    case 0:
        for(size_t k = 0; k < 4; k++)
            vec_local_node_ids[k] = k;
        if(order==2)
            for(size_t k = 0; k < 4; k++)
                vec_local_node_ids[k + 4] = k + 8;
        break;
    case 1:
        for(size_t k = 0; k < 4; k++)
            vec_local_node_ids[k] = k + 4;
        if(order==2)
            for(size_t k = 0; k < 4; k++)
                vec_local_node_ids[k + 4] = k + 12;
        break;
    case 2:
        vec_local_node_ids[0] = 0;
        vec_local_node_ids[1] = 4;
        vec_local_node_ids[2] = 5;
        vec_local_node_ids[3] = 1;
        if(order==2)
        {
            vec_local_node_ids[4] = 16;
            vec_local_node_ids[5] = 12;
            vec_local_node_ids[6] = 17;
            vec_local_node_ids[7] = 8;
        }
        break;
    case 3:
        vec_local_node_ids[0] = 1;
        vec_local_node_ids[1] = 5;
        vec_local_node_ids[2] = 6;
        vec_local_node_ids[3] = 2;
        if(order==2)
        {
            vec_local_node_ids[4] = 17;
            vec_local_node_ids[5] = 13;
            vec_local_node_ids[6] = 18;
            vec_local_node_ids[7] = 9;
        }

        break;
    case 4:
        vec_local_node_ids[0] = 2;
        vec_local_node_ids[1] = 6;
        vec_local_node_ids[2] = 7;
        vec_local_node_ids[3] = 3;
        if(order==2)
        {
            vec_local_node_ids[4] = 18;
            vec_local_node_ids[5] = 14;
            vec_local_node_ids[6] = 19;
            vec_local_node_ids[7] = 10;
        }
        break;
    case 5:
        vec_local_node_ids[0] = 0;
        vec_local_node_ids[1] = 3;
        vec_local_node_ids[2] = 7;
        vec_local_node_ids[3] = 4;
        if(order==2)
        {
            vec_local_node_ids[4] = 11;
            vec_local_node_ids[5] = 19;
            vec_local_node_ids[6] = 15;
            vec_local_node_ids[7] = 16;
        }
        break;
    }
}

void PrismTopology::getLocalNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_local_node_ids)
{
    size_t n_nodes = LineTopology::getNumberOfNodes(order);
    vec_local_node_ids.resize(n_nodes); 
    if (edge_id<6) {
        vec_local_node_ids[0] = edge_id;
        vec_local_node_ids[1] = (edge_id + 1) % 3 + 3 * (int)(edge_id / 3);
    } else {
        vec_local_node_ids[0] = edge_id%3;
        vec_local_node_ids[1] = edge_id % 3 + 3;
    }
    if (order==2)
        vec_local_node_ids[2] = edge_id+6;
}

void PrismTopology::getLocalNodeIDsOfFaceElement(size_t order, size_t face_id, std::vector<size_t> &vec_local_node_ids)
{
    size_t n_nodes = 0;
    if (face_id<2) {
        n_nodes = TriangleTopology::getNumberOfNodes(order);
    } else {
        n_nodes = Quad8Topology::getNumberOfNodes(order);
    }
    vec_local_node_ids.resize(n_nodes); 
    switch(face_id)
    {
    case 0:
        for(size_t k = 0; k < 3; k++)
            vec_local_node_ids[k] = k;
        if(order==2)
        {
            for(size_t k = 0; k < 3; k++)
                vec_local_node_ids[k + 3] = k + 6;
        }
        break;
    case 1:
        for(size_t k = 0; k < 3; k++)
            vec_local_node_ids[k] = k + 3;
        if(order==2)
        {
            for(size_t k = 0; k < 3; k++)
                vec_local_node_ids[k + 3] = k + 9;
        }
        break;
    case 2:
        vec_local_node_ids[0] = 1;
        vec_local_node_ids[1] = 2;
        vec_local_node_ids[2] = 5;
        vec_local_node_ids[3] = 4;
        if(order==2)
        {
            vec_local_node_ids[4] = 7;
            vec_local_node_ids[5] = 14;
            vec_local_node_ids[6] = 10;
            vec_local_node_ids[7] = 13;
        }
        break;
    case 3:
        vec_local_node_ids[0] = 5;
        vec_local_node_ids[1] = 2;
        vec_local_node_ids[2] = 0;
        vec_local_node_ids[3] = 3;
        if(order==2)
        {
            vec_local_node_ids[4] = 14;
            vec_local_node_ids[5] =  8;
            vec_local_node_ids[6] = 12;
            vec_local_node_ids[7] = 10;
        }
        break;
    case 4:
        vec_local_node_ids[0] = 0;
        vec_local_node_ids[1] = 1;
        vec_local_node_ids[2] = 4;
        vec_local_node_ids[3] = 3;
        if(order==2)
        {
            vec_local_node_ids[4] = 6;
            vec_local_node_ids[5] = 13;
            vec_local_node_ids[6] = 9;
            vec_local_node_ids[7] = 12;
        }
        break;
    }
}

void PyramidTopology::getLocalNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_local_node_ids)
{
    size_t n_nodes = LineTopology::getNumberOfNodes(order);
    vec_local_node_ids.resize(n_nodes); 
    if (edge_id < 4) {
        vec_local_node_ids[0] = edge_id;
        vec_local_node_ids[1] = (edge_id+1)%4;
    } else {
        vec_local_node_ids[0] = edge_id % 4;
        vec_local_node_ids[1] = 4;
    }
    if (order==2)
        vec_local_node_ids[2] = edge_id+4;
}

void PyramidTopology::getLocalNodeIDsOfFaceElement(size_t order, size_t face_id, std::vector<size_t> &vec_local_node_ids)
{
    size_t n_nodes = 0;
    if (face_id==0) {
        n_nodes = Quad8Topology::getNumberOfNodes(order);
    } else {
        n_nodes = TriangleTopology::getNumberOfNodes(order);
    }
    vec_local_node_ids.resize(n_nodes); 
    switch(face_id)
    {
    case 0:
        for(size_t k = 0; k < 4; k++)
            vec_local_node_ids[k] = k;
        if(order==2)
        {
            for(size_t k = 0; k < 4; k++)
                vec_local_node_ids[k + 4] = k + 5;
        }
        break;
    case 1:
        vec_local_node_ids[0] = 0;
        vec_local_node_ids[1] = 1;
        vec_local_node_ids[2] = 4;
        if(order==2)
        {
            vec_local_node_ids[3] = 5;
            vec_local_node_ids[4] = 10;
            vec_local_node_ids[5] = 9;
        }
        break;
    case 2:
        vec_local_node_ids[0] = 1;
        vec_local_node_ids[1] = 2;
        vec_local_node_ids[2] = 4;
        if(order==2)
        {
            vec_local_node_ids[3] = 6;
            vec_local_node_ids[4] = 11;
            vec_local_node_ids[5] = 10;
        }
        break;
    case 3:
        vec_local_node_ids[0] = 2;
        vec_local_node_ids[1] = 3;
        vec_local_node_ids[2] = 4;
        if(order==2)
        {
            vec_local_node_ids[3] = 7;
            vec_local_node_ids[4] = 12;
            vec_local_node_ids[5] = 11;
        }
        break;
    case 4:
        vec_local_node_ids[0] = 3;
        vec_local_node_ids[1] = 0;
        vec_local_node_ids[2] = 4;
        if(order==2)
        {
            vec_local_node_ids[3] = 8;
            vec_local_node_ids[4] = 9;
            vec_local_node_ids[5] = 12;
        }
        break;
    }
}


} // end namespace
