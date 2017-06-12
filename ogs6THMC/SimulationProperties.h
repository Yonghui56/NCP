/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SimulationProperties.h
 *
 * Created on 2012-07-05 by Norihiro Watanabe
 */

#pragma once

#include <string>

#include "logog.hpp"
#include "tinyxml2.h"

#include "BaseLib/Options.h"
#include "Ogs5ToOgs6.h"


/*

xml
    ProcessData
        Process, type=gw, name=gw1, time_group_id=0
            file
        Process, type=element_velocity, name=ev
            file
    TimeData
        TimeGroup, id=0, start=0, end=30, unit=DAY, control=STEPS
    coupling
        P
            M, pcs=gw1
                out, name=p
            M, pcs=ev
                in, name=p
                out, name=v

 */

inline bool readSimulationProperties(const std::string &xml_file, BaseLib::Options &properties)
{
    tinyxml2::XMLDocument doc;
    int ret = doc.LoadFile(xml_file.c_str());
    if (ret != tinyxml2::XML_NO_ERROR) {
        LOGOG_CERR << "Error in reading a XML file " << xml_file << " with error code " << ret << std::endl;
        return false;
    }
    
    // quick check
    if (doc.FirstChildElement("ProcessData")==0) {
        LOGOG_CERR << "Error in reading a XML file: ProcessData tag is not found." << std::endl;
        return false;
    }
    if (doc.FirstChildElement("TimeData")==0) {
        LOGOG_CERR << "Error in reading a XML file: TimeData tag is not found." << std::endl;
        return false;
    }
    
    // read data
    tinyxml2::XMLElement* e;
    BaseLib::Options* optProcessData = properties.addSubGroup("ProcessData");
    for (e = doc.FirstChildElement("ProcessData")->FirstChildElement("Process"); e != 0; e->NextSibling()) {
        const char* pcs_type = e->Attribute("type");
        if (pcs_type==0) {
            ERR("Error in reading a XML file: attribute \"type\" in Process tag must be specified.");
            return false;
        }
        const char* pcs_name = e->Attribute("name");
        if (pcs_name == 0)
            pcs_name = pcs_type;
        int tim_id = -1;
        e->QueryIntAttribute("time_group_id", &tim_id);
        
        
        BaseLib::Options* optProcess = optProcessData->addSubGroup(pcs_name);
        optProcess->addOption("type", pcs_type);
        optProcess->addOptionAsNum("time_group_id", tim_id);
    }

    BaseLib::Options* optTimeData = properties.addSubGroup("TimeData");
    for (e = doc.FirstChildElement("TimeData")->FirstChildElement("TimeGroup"); e != 0; e->NextSibling()) {
        int tim_id = e->IntAttribute("id");
        double tim_start = e->DoubleAttribute("start");
        double tim_end = e->DoubleAttribute("end");
        const char* tim_unit = e->Attribute("unit");
        const char* tim_control = e->Attribute("control");
        
        
        BaseLib::Options* optTimeGroup = optTimeData->addSubGroup("TimeData");
        optTimeGroup->addOptionAsNum("time_group_id", tim_id);
    }

    
    return true;
}



#if 0
//inline void readOgs5(const std::string &input_path)
inline int readOgs5(SimulationInfo &sim_info)
{
    // fem
    ogs5::Ogs5FemData ogs5femdata;
    ogs5femdata.read(sim_info.getProjectPath());

    // gli
    std::string geo_unique_name;
    std::vector<std::string> geo_errors;
    GeoLib::GEOObjects geo;
    Ogs4GeoIO::readGLIFileV4(sim_info.getProjectPath()+".gli", &geo, geo_unique_name, geo_errors);

    // mesh
    std::vector<MeshLib::IMesh*> msh_vector;
    Ogs5MeshIO::readMesh(sim_info.getProjectPath()+".msh", msh_vector);
    if (msh_vector.size()==0) {
        LOGOG_CERR << " Error: Cannot find mesh " << sim_info.getProjectPath()+".msh" << std::endl;
        return -1;
    }

    // ddc

    //-------------------------------------------------------------------------
    // Setup simulation
    //-------------------------------------------------------------------------
    MeshLib::IMesh* msh = msh_vector[0];
    DiscreteLib::DiscreteSystem dis(*msh);
    // - create pcs
    Ogs5ToOgs6::convert(ogs5femdata, geo, geo_unique_name, dis);
    // - ddc
}
#endif

