/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DataType.h
 *
 * Created on 2012-08-31 by Norihiro Watanabe
 */

#pragma once

#include <Eigen>
#include "MathLib/LinAlg/LinearEquation/EigenDenseLinearEquation.h"

namespace MathLib
{
/// Local linear equation type
typedef MathLib::EigenDenseLinearEquation LocalEquation;
/// Local dense matrix type (row-majored)
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> LocalMatrix;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SparseMatrix;
/// Local dense vector type
typedef Eigen::VectorXd LocalVector;

}
