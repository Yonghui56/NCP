/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DiscreteSystemList.h
 *
 * Created on 2012-08-18 by Norihiro Watanabe
 */


#pragma once

#include "DiscreteLib/Serial/DiscreteSystem.h"
#include "DiscreteLib/OpenMP/OMPDiscreteSystem.h"
#include "DiscreteLib/SerialNodeDdc/SerialNodeDdcSharedDiscreteSystem.h"
#include "DiscreteLib/SerialNodeDdc/SerialNodeDdcDistributedDiscreteSystem.h"
#if defined(USE_LIS) && defined(USE_MPI)
#include "DiscreteLib/lis/LisDiscreteSystem.h"
#include "DiscreteLib/lis/LisMPILinearSystem.h"
#endif

