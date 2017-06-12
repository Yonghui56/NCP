/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file MeshIoOgs5.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "MeshIoOgs5.h"

#include "MeshLib/Core/UnstructuredMesh.h"
#include "MeshLib/Core/ElementFactory.h"
#include "MeshLib/Tools/Tools.h"

namespace MeshIoOgs5
{

const std::string convertElementType2String(const MeshLib::ElementShape::type t)
{
    if (t == MeshLib::ElementShape::LINE)            return "line";
    if (t == MeshLib::ElementShape::QUAD)            return "quad";
    if (t == MeshLib::ElementShape::HEXAHEDRON)    return "hex";
    if (t == MeshLib::ElementShape::TRIANGLE)        return "tri";
    if (t == MeshLib::ElementShape::TETRAHEDRON)    return "tet";
    if (t == MeshLib::ElementShape::PRISM)        return "pris";
    if (t == MeshLib::ElementShape::PYRAMID)        return "pyra";
    return "none";
};

MeshLib::ElementShape::type convertString2ElementType(const std::string &s)
{
    if (s.compare("line") == 0) return MeshLib::ElementShape::LINE;
    if (s.compare("quad") == 0) return MeshLib::ElementShape::QUAD;
    if (s.compare("hex")  == 0) return MeshLib::ElementShape::HEXAHEDRON;
    if (s.compare("tri")  == 0) return MeshLib::ElementShape::TRIANGLE;
    if (s.compare("tet")  == 0) return MeshLib::ElementShape::TETRAHEDRON;
    if (s.compare("pris") == 0) return MeshLib::ElementShape::PRISM;
    if (s.compare("pyra") == 0) return MeshLib::ElementShape::PYRAMID;
    return MeshLib::ElementShape::INVALID;
};

void readMesh(std::string const& fileName, std::vector<MeshLib::IMesh*> &vec_mesh) 
{

    std::ifstream if_file (fileName.data(),std::ios::in);
    if (!if_file.is_open())
    {
      std::cout << "CFEMesh::FEMRead() - Could not open file...\n";
      return;
    }
    
    std::cout << "MSHRead:  ASCII file" << std::endl;
    std::string line_string ("");
    
    MeshLib::UnstructuredMesh *msh = NULL;
    
    while (!if_file.eof())
    {
      getline(if_file, line_string);
    
      if (line_string.find("#STOP") != std::string::npos)
        break;
    
      if (line_string.find("#FEM_MSH")!=std::string::npos) {
        msh = new MeshLib::UnstructuredMesh();
        vec_mesh.push_back(msh);
        continue;
      } 
      
      if (!msh) continue;
    
      if (line_string.find("$NODES") != std::string::npos)
      {
        double x, y, z;
        size_t no_nodes, idx;
        if_file >> no_nodes >> std::ws;
        std::string s;
        std::ios::pos_type position = if_file.tellg();
        for (size_t i = 0; i < no_nodes; i++)
        {
          if_file >> idx >> x >> y >> z;
          GeoLib::Point p(x, y, z);
          msh->setNodeCoordinates(idx, p);
          position = if_file.tellg();
          if_file >> s;
          if_file.seekg(position, std::ios::beg);
          if_file >> std::ws;
        }
      }
      else  if (line_string.find("$ELEMENTS") != std::string::npos)
      {
        size_t no_elements, idx, group_id, idummy;
        if_file >> no_elements >> std::ws;
        for (size_t i = 0; i < no_elements; i++)
        {
          std::string buffer("");
          std::string name("");
          if_file >> idx >> group_id;
          if_file >> buffer;
    
          if (buffer.find("-1") != std::string::npos) {
            //int grid_adaptation = strtol(buffer.data(), NULL, 0);
            if_file >> name;
          } else {
            name = buffer;
          }
    
          MeshLib::ElementShape::type ele_type = convertString2ElementType(name);
          MeshLib::IElement *newElem = MeshLib::ElemenetFactory::createNewElement(ele_type);
          newElem->setID(idx);
          newElem->setGroupID(group_id);
          for (size_t j=0; j < newElem->getNumberOfNodes(); j++) {
            if_file >> idummy;
            newElem->setNodeID(j, idummy);
          }
          msh->addElement(newElem);
        }
      }

    }
    
    if_file.close();
};

void write(MeshLib::IMesh const * /*mesh*/, std::string const& /*fileName*/)
{
    std::cout << "Not implemented yet." << std::endl;
};
}
