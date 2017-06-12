/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file MeshIoLegacyVtk.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "MeshIoLegacyVtk.h"

#include "MeshLib/Core/IMesh.h"

using namespace MeshLib;

void MeshIOLegacyVtk::WriteAsciiFile(const std::string &str_file_name, const MeshLib::IMesh &mesh)
{
    std::fstream vtk_file (str_file_name.c_str(), std::ios::out);
    vtk_file.setf(std::ios::scientific,std::ios::floatfield);
    vtk_file.precision(12);
    if (!vtk_file.good())
    {
        std::cout << "Could not open file for writing: " << str_file_name << std::endl;
        return;
    }
    vtk_file.seekg(0L,std::ios::beg);
    WriteVTKHeader(vtk_file);
    WriteVTKPoints(vtk_file, mesh);
    WriteVTKCells(vtk_file, mesh);
    WriteVTKCellDataHeader(vtk_file, mesh);
    WriteVTKMaterialID(vtk_file, mesh);

    vtk_file.close();
};

void MeshIOLegacyVtk::WriteVTKHeader(std::fstream &vtk_file)
{
    vtk_file << "# vtk DataFile Version 3.0" << std::endl;
    vtk_file << "Unstructured Grid from OpenGeoSys" << std::endl;
    vtk_file << "ASCII"  << std::endl;
    vtk_file << "DATASET UNSTRUCTURED_GRID"  << std::endl;
}


void MeshIOLegacyVtk::WriteVTKPoints(std::fstream &vtk_file, const MeshLib::IMesh &mesh)
{
    const size_t n_nodes(mesh.getNumberOfNodes());
    vtk_file << "POINTS "<< n_nodes << " double" << std::endl;

    for(size_t i = 0; i < n_nodes ; i++)
    {
        const GeoLib::Point *pt = mesh.getNodeCoordinatesRef(i);
        vtk_file << (*pt)[0] << " " << (*pt)[1] << " " << (*pt)[2] << std::endl;
    }
}

void MeshIOLegacyVtk::WriteVTKCells(std::fstream &vtk_file, const MeshLib::IMesh &mesh)
{
    size_t numCells = mesh.getNumberOfElements();

    // count overall length of element vector
    long numAllPoints =0;
    for(size_t i=0; i < numCells; i++)
    {
        IElement* ele = mesh.getElement(i);
        numAllPoints = numAllPoints + (ele->getNumberOfNodes()) + 1;
    }

    // write elements
    vtk_file << "CELLS " << numCells << " " << numAllPoints << std::endl;
    for(size_t i=0; i < numCells; i++)
    {
        IElement* ele = mesh.getElement(i);

        // Write number of points per cell
        switch(ele->getShapeType())
        {
        case ElementShape::LINE:
            vtk_file << "2"; break;
        case ElementShape::QUAD:
            vtk_file << "4"; break;
        case ElementShape::HEXAHEDRON:
            vtk_file << "8"; break;
        case ElementShape::TRIANGLE:
            vtk_file << "3"; break;
        case ElementShape::TETRAHEDRON:
            vtk_file << "4"; break;
        case ElementShape::PRISM:
            vtk_file << "6"; break;
        case ElementShape::PYRAMID:
            vtk_file << "5"; break;
        default:
            std::cerr << "COutput::WriteVTKElementData MshElemType not handled" << std::endl;
            break;
        }

        for(size_t j = 0; j < ele->getNumberOfNodes(); j++)
            vtk_file << " " << ele->getNodeID(j);

        vtk_file << std::endl;
    }
    vtk_file << std::endl;

    // write cell types
    vtk_file << "CELL_TYPES " << numCells << std::endl;

    for(size_t i=0; i < numCells; i++)
    {
        IElement* ele = mesh.getElement(i);

        // Write vtk cell type number (see vtkCellType.h)
        switch(ele->getShapeType())
        {
        case ElementShape::LINE:
            vtk_file << "3" << std::endl; break;
        case ElementShape::QUAD:
            vtk_file << "9" << std::endl; break;
        case ElementShape::HEXAHEDRON:
            vtk_file << "12" << std::endl; break;
        case ElementShape::TRIANGLE:
            vtk_file << "5" << std::endl; break;
        case ElementShape::TETRAHEDRON:
            vtk_file << "10" << std::endl; break;
        case ElementShape::PRISM: // VTK_WEDGE
            vtk_file << "13" << std::endl; break;
        case ElementShape::PYRAMID:
            vtk_file << "14" << std::endl; break;
        default:
            std::cerr << "COutput::WriteVTKElementData MshElemType not handled" << std::endl;
            break;
        }
    }
    vtk_file << std::endl;
}

void MeshIOLegacyVtk::WriteVTKCellDataHeader(std::fstream &vtk_file, const IMesh &mesh)
{
    vtk_file << "CELL_DATA " << mesh.getNumberOfElements() << std::endl;
}

void MeshIOLegacyVtk::WriteVTKMaterialID(std::fstream &vtk_file, const IMesh &mesh)
{
    vtk_file << "SCALARS " << "GroupID" << " int 1" << std::endl;
    vtk_file << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < mesh.getNumberOfElements(); i++) {
        IElement* e = mesh.getElement(i);
        vtk_file << e->getGroupID() << std::endl;
    }
}





