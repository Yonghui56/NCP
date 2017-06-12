/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IClonable.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

namespace NumLib
{

class IClonable
{
public:
    virtual ~IClonable() {};
    virtual IClonable* clone() const = 0;
};

}
