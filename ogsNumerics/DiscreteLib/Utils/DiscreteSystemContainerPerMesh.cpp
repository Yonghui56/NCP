/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DiscreteSystemContainerPerMesh.cpp
 *
 * Created on 2012-08-17 by Norihiro Watanabe
 */

#include "DiscreteSystemContainerPerMesh.h"

namespace DiscreteLib
{

DiscreteSystemContainerPerMesh* DiscreteSystemContainerPerMesh::_obj = 0;

DiscreteSystemContainerPerMesh* DiscreteSystemContainerPerMesh::getInstance()
{
    if (_obj==0) _obj = new DiscreteSystemContainerPerMesh();
    return _obj;
}


}
