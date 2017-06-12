/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file GeoType.h
 *
 * Created on 2010-06-17 by Thomas Fischer
 */

#ifndef GEOTYPE_H_
#define GEOTYPE_H_

#include <string>

namespace GeoLib 
{

/**
 * \ingroup GeoLib
 */

enum GEOTYPE 
{
    INVALID = 0,
    POINT,     //!< POINT
    POLYLINE,  //!< POLYLINE
    POLYGON,   //!< POLYGON  
    SURFACE,   //!< SURFACE
    VOLUME,    //!< VOLUME
    GEODOMAIN  //!< DOMAIN (ENTIRE SPACE)
};

GEOTYPE convertGeoType (const std::string& geo_type_str);

std::string convertGeoTypeToString (GEOTYPE geo_type);

} // end namespace GEOLIB

#endif /* GEOTYPE_H_ */
