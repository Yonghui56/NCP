/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file GeoDomain.h
 *
 * Created on 2012-09-23 by Norihiro Watanabe
 */

#include "GeoObject.h"

namespace GeoLib
{

/**
 * \brief Domain class representing entire space
 */
class GeoDomain : public GeoObject
{
public:
    ///
    virtual ~GeoDomain() {}

    /// return the geometric type
    virtual GEOTYPE getGeoType() const {return GeoLib::GEODOMAIN;};
};
    
} //end
