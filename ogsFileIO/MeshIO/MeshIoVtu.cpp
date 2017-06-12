/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file MeshIoVtu.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "MeshIoVtu.h"

#include <fstream>

#include "BaseLib/SystemTools.h"
#include "BaseLib/FileTools.h"

using namespace std;
using namespace MeshLib;

MeshIoVtu::MeshIoVtu(bool binary_mode)
: VtuWriter(binary_mode)
{
}

bool MeshIoVtu::WriteXMLUnstructuredGrid(const std::string &vtkfile,
                                    MeshLib::IMesh &msh)
{
    //-------------------------------------------------------------------------
    //# Setup file stream
    //-------------------------------------------------------------------------
    std::fstream fin;
    if (_useBinary)
        fin.open(vtkfile.data(), std::ios::out | std::ios::binary);
    else
        fin.open(vtkfile.data(), std::ios::out);

    if (!fin.good())
    {
        std::cout << "***Warning: Cannot open the output file, " << vtkfile << std::endl;
        return false;
    }

    if (!_useBinary)
    {
        fin.setf(std::ios::scientific,std::ios::floatfield);
        fin.precision(12);
    }

    //-------------------------------------------------------------------------
    //# Output
    //-------------------------------------------------------------------------
    long offset = 0;

    string str_format;
    if (!_useBinary)
        str_format = "ascii";
    else
        str_format = "appended";
    bool data_out = !_useBinary;

    //# Header
    fin << "<?xml version=\"1.0\"?>" << std::endl;
    fin << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"";
    if (!_useBinary || _isLittleEndian)
        fin << " byte_order=\"LittleEndian\"";
    else
        fin << " byte_order=\"BigEndian\"";
    fin << ">" << std::endl;
    //  fin << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">"  << endl;

    //# Unstructured Grid information
    fin << "  <UnstructuredGrid>" << endl;
    fin << "    <Piece NumberOfPoints=\"" << msh.getNumberOfNodes() <<
    "\" NumberOfCells=\"" << msh.getNumberOfElements() << "\">" << std::endl;
    //....................................................................
    // Nodes
    fin << "      <Points>" << std::endl;
    writeDataArrayHeader(fin, type_Double, "", 3, str_format, offset);
    writePoints(fin, data_out, msh, offset);
    writeDataArrayFooter(fin);
    fin << "      </Points>" << endl;
    //....................................................................
    // Elements
    fin << "      <Cells>" << endl;
    //connectivity
    writeDataArrayHeader(fin, type_Long, "connectivity", 0, str_format, offset);
    long sum_ele_components = 0;
    writeCellConnectivity(fin, data_out, msh, offset, sum_ele_components);
    writeDataArrayFooter(fin);
    //offset
    writeDataArrayHeader(fin, type_Long, "offsets", 0, str_format, offset);
    writeCellOffset(fin, data_out, msh, offset);
    writeDataArrayFooter(fin);
    //type
    writeDataArrayHeader(fin, type_UChar, "types", 0, str_format, offset);
    writeCellType(fin, data_out, msh, offset);
    writeDataArrayFooter(fin);
    fin << "      </Cells>" << endl;
    //....................................................................
    // Element values
    fin << "      <CellData>" << endl;
    writeElementGroupID(fin, data_out, msh, offset);
    fin << "      </CellData>" << endl;
    fin << "    </Piece>" << endl;
    fin << "  </UnstructuredGrid>" << endl;

    //======================================================================
    // Raw data (for binary mode)
    if (_useBinary)
    {
        fin << "  <AppendedData encoding=\"raw\">" << endl;
        fin << "    _";

        //Node
        this->writePoints(fin, true, msh, offset);
        //Element
        //conncectivity
        this->writeCellConnectivity(fin, true, msh, offset, sum_ele_components);
        //offset
        this->writeCellOffset(fin, true, msh, offset);
        //type
        this->writeCellType(fin, true, msh, offset);
        // Nodal values
        //this->WriteNodalValue(fin, true, out, msh, offset);
        // Elemental values
        this->writeElementGroupID(fin, true,  msh, offset);

        fin << endl;
        fin << "  </AppendedData>" << endl;
    }

    fin << "</VTKFile>" << endl;
    fin.close();

    return true;
}


