/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file MeshIoVtu.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include <vector>

#include "CommonIO/VtuWriter.h"

class MeshIoVtu : public VtuWriter
{
public:
    explicit MeshIoVtu(bool binary_mode);
    virtual ~MeshIoVtu(void){}

public:
    bool WriteXMLUnstructuredGrid(const std::string &vtkfile,
                                  MeshLib::IMesh &msh);
};
