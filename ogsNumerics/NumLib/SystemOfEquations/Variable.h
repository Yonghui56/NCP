/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Variable.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <string>

namespace NumLib
{

struct Variable
{
    size_t id;
    size_t n_dof;
    std::string name;

    Variable(size_t i, size_t n, std::string str) : id(i), n_dof(n), name(str) {};
};


} //end

