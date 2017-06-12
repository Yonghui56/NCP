/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LinearSolverList.h
 *
 * Created on 2012-08-18 by Norihiro Watanabe
 */


#pragma once

#include "MathLib/LinAlg/LinearEquation/DenseLinearEquation.h"
#include "MathLib/LinAlg/LinearEquation/SparseLinearEquation.h"
#include "MathLib/LinAlg/LinearEquation/EigenDenseLinearEquation.h"
#ifdef USE_LIS
#include "MathLib/LinAlg/LinearEquation/LisLinearEquation.h"
#endif
#ifdef USE_PARDISO
#include "MathLib/LinAlg/LinearEquation/PARDISOLinearEquation.h"
#endif
#ifdef USE_PETSC
#include "MathLib/LinAlg/LinearEquation/PETScLinearEquation.h"
#endif

