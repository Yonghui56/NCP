/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IDiscreteObject.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <cstddef>

namespace DiscreteLib
{

/**
 * \brief Object in discrete systems
 */
class IDiscreteObject
{
public:
    ///
    IDiscreteObject() : _obj_id(0) {};
    ///
    virtual ~IDiscreteObject() {};
    ///
    size_t getObjectID() const {return _obj_id;};
    ///
    void setObjectID(size_t i) {_obj_id = i;};
    
private:
    size_t _obj_id;
};

} // end
