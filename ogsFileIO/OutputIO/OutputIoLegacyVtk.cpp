/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file OutputIoLegacyVtk.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "OutputIoLegacyVtk.h"

#include "MeshIO/MeshIoLegacyVtk.h"

using namespace MeshLib;

void OutputIOLegacyVtk::WriteAsciiFile(std::string str_file_name, const IMesh &mesh, int time_step_number, double simulation_time, std::vector<NodalScalarValue> &nodalScalarValues) {
    std::fstream vtk_file (str_file_name.c_str(), std::ios::out);
    vtk_file.setf(std::ios::scientific,std::ios::floatfield);
    vtk_file.precision(12);
    if (!vtk_file.good())
    {
        std::cout << "Could not open file for writing: " << str_file_name << std::endl;
        return;
    }
    vtk_file.seekg(0L,std::ios::beg);
    WriteVTKHeader(vtk_file, time_step_number, simulation_time, nodalScalarValues);
    MeshIOLegacyVtk::WriteVTKPoints(vtk_file, mesh);
    MeshIOLegacyVtk::WriteVTKCells(vtk_file, mesh);
    WriteVTKDataArrays(vtk_file, mesh, nodalScalarValues);

    vtk_file.close();
};

void OutputIOLegacyVtk::WriteVTKHeader(std::fstream &vtk_file, int time_step_number, double simulation_time, std::vector<NodalScalarValue> &/*nodalScalarValues*/)
{
    MeshIOLegacyVtk::WriteVTKHeader(vtk_file);
    // time information
    // see http://www.visitusers.org/index.php?title=Time_and_Cycle_in_VTK_files
    vtk_file << "FIELD TimesAndCycles 2" << std::endl;
    vtk_file << "TIME 1 1 double" << std::endl;
    vtk_file << simulation_time << std::endl;
    vtk_file << "CYLCE 1 1 long" << std::endl;
    vtk_file << time_step_number << std::endl;
}

void OutputIOLegacyVtk::WriteVTKDataArrays(std::fstream &vtk_file, const IMesh &mesh, std::vector<NodalScalarValue> &nodalScalarValues)
{
    const size_t n_nodes = mesh.getNumberOfNodes();
    vtk_file << "POINT_DATA " << n_nodes << std::endl;
    for (size_t i=0; i<nodalScalarValues.size(); i++) {
        vtk_file << "SCALARS " << nodalScalarValues.at(i).first << " double 1" << std::endl;
        vtk_file << "LOOKUP_TABLE default" << std::endl;
        for (size_t j = 0; j<n_nodes; j++)
        {
            vtk_file << nodalScalarValues.at(i).second[j] << std::endl;
        }
    }
}
