/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SparseTableCRS.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <set>
#include <cstddef>

namespace MathLib 
{

/**
 * Row-major sparsity class
 */
typedef std::vector<std::set<size_t> > RowMajorSparsity;

} // end namespace MathLib

