/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file MeshIoOgs5.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once
#include <iostream>
#include <fstream>
#include <string>

#include "MeshLib/Core/IMesh.h"

namespace MeshIoOgs5
{
void readMesh(std::string const& fileName, std::vector<MeshLib::IMesh*> &vec_mesh);
void writeMesh(MeshLib::IMesh const * mesh, std::string const& fileName);
const std::string convertElementType2String(const MeshLib::ElementShape::type t);
MeshLib::ElementShape::type convertString2ElementType(const std::string &s);
}

