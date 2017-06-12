/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SimulationInfo.cpp
 *
 * Created on 2012-07-05 by Norihiro Watanabe
 */

#include "SimulationInfo.h"

#include <cstdio>
#include <iostream>

#include "logog.hpp"

#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"

#include "Configure.h"
#include "BuildInfo.h"

namespace ogs6
{

void SimulationInfo::outputHeader ( void )
{
//    char buf[128];
//
//    const int len = 47;

    INFO("");
    INFO("          ###################################################");
    INFO("          ##                                               ##");
    INFO("          ##              OpenGeoSys-Project 6             ##");
#ifdef USE_LIS
    INFO("          ## %s ##", BaseLib::bothPadding("powered by LIS",45).c_str());
#endif
    INFO("          ##                                               ##");
    INFO("          ##   Contributors                                ##");
    INFO("          ##   * Helmholtz Centre for Environmental        ##");
    INFO("          ##     Research - UFZ                            ##");
    INFO("          ##   * TU Dresden                                ##");
    INFO("          ##   * University of Kiel                        ##");
    INFO("          ##   * University of Edinburgh                   ##");
    INFO("          ##   * University of Tuebingen (ZAG)             ##");
    INFO("          ##   * Federal Institute for Geosciences         ##");
    INFO("          ##     and Natural Resources (BGR)               ##");
    INFO("          ##   * Helmholtz Centre Potsdam GFZ              ##");
    INFO("          ##     German Research Centre for Geosciences    ##");
    INFO("          ##                                               ##");
    INFO("          ##   Program version                             ##");
    INFO("          ##   * Version: %s ##", BaseLib::rightPadding(OGS_VERSION, 32).c_str());
    INFO("          ##   * Date   : %s ##", BaseLib::rightPadding(OGS_DATE, 32).c_str());
#ifdef GIT_COMMIT_INFO
    INFO("          ##   * Rev.   :                                  ##");
    INFO("          ##     %s ##", BaseLib::rightPadding(GIT_COMMIT_INFO, 41).c_str());
#endif
    INFO("          ##                                               ##");
    INFO("          ###################################################");
    INFO("");
    INFO("");
}

SimulationInfo::SimulationInfo(const std::string &project_path, const std::string &output_dir_path)
{
    this->setProjectPath(project_path);
    if (output_dir_path.length()>0) {
        _output_dir = output_dir_path;
    } else {
        _output_dir = _project_dir;
    }
}


void SimulationInfo::setProjectPath(const std::string& path)
{
    _project_path = path;
    _project_dir = BaseLib::getFileDirecotryPath(path);
    _project_name = BaseLib::getFileBaseName(path);
};

} //end
