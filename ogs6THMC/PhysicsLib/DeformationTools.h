/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DeformationTools.h
 *
 * Created on 2012-11-30 by Norihiro Watanabe
 */


#pragma once

#include <string>
#include "MathLib/DataType.h"

inline size_t getNumberOfStrainComponents(size_t dim) {return (dim==2 ? 4 : 6);};

inline MathLib::LocalMatrix get_m(const size_t dim)
{
    MathLib::LocalMatrix m = MathLib::LocalMatrix::Zero(getNumberOfStrainComponents(dim),1);
    for (size_t i=0; i<dim; i++)
        m(i,0) = 1.0;
    return m;
}

inline const std::string getStressStrainComponentPostfix(size_t i)
{
    switch (i) {
    case 0: return "_XX";
    case 1: return "_YY";
    case 2: return "_ZZ";
    case 3: return "_XY";
    case 4: return "_XZ";
    case 5: return "_YZ";
    }
    return "";
}

inline const std::string getDisplacementComponentPostfix(size_t i)
{
    switch (i) {
    case 0: return "_X";
    case 1: return "_Y";
    case 2: return "_Z";
    }
    return "";
}

template<class MyVariable>
inline MyVariable* getDisplacementComponentVariable(MyVariable *u_x, MyVariable* u_y,
        MyVariable* u_z, const std::string &var_name)
{
    if (var_name.find("_X") != std::string::npos) {
        return u_x;
    } else if (var_name.find("_Y") != std::string::npos) {
        return u_y;
    } else {
        return u_z;
    }
}

