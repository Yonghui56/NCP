/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file GeoProcessBuilder.cpp
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#include "GeoProcessBuilder.h"

#include "ProcessLib/ProcessFactoryImpl.h"
#include "UserModules/ProcessList.h"
#include "DiscreteSystemList.h"
#include "LinearSolverList.h"

ProcessLib::AbstractProcessBuilder* GeoProcessBuilder::_obj = 0;

ProcessLib::AbstractProcessBuilder* GeoProcessBuilder::getInstance()
{
    if (_obj==0) _obj = new GeoProcessBuilder();
    return _obj;
}


GeoProcessBuilder::GeoProcessBuilder()
{
#include "UserModules/ProcessReg.h"
}
