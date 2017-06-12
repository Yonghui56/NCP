/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Ogs4GeoIO.h
 *
 * Created on 2010-01-14 by Thomas Fischer
 */


#ifndef OGSIOVER4_H_
#define OGSIOVER4_H_

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// forward declaration
namespace GeoLib
{
class GEOObjects;
}

namespace Ogs4GeoIO
{
/** I/O - routines for the OGS-4 gli file format */

/** method reads geometric objects from file in gli format */
bool readGLIFileV4 (const std::string& fname, GeoLib::GEOObjects* geo, std::string& unique_name, std::vector<std::string>& errors);

void writeGLIFileV4 (const std::string& fname,
                     const std::string& proj_name,
                     const GeoLib::GEOObjects& geo);

void writeAllDataToGLIFileV4 (const std::string& fname, const GeoLib::GEOObjects& geo);
} // end namespace

#endif /* OGSIOVER4_H_ */
