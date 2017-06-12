/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemElementalFunction.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <algorithm>

#include "MeshLib/Core/IMesh.h"

#include "NumLib/Function/Function.h"

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Tools/FemElementObjectContainer.h"

namespace FemLib
{



/**
 * \brief Template class for FEM element-based functions
 */
template<typename Tvalue>
class TemplateFEMElementalFunction : public NumLib::TemplateFunction<GeoLib::Point,Tvalue>
{
public:
    void eval(const GeoLib::Point &pt, Tvalue &v) {
        v = _ele_values[0];
    };
    Tvalue& getValue(size_t id) const {
        return _ele_values[id];
    };
private:
    Tvalue* _ele_values;
};


}



