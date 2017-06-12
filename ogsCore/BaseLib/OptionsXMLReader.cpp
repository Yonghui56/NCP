/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file OptionsXMLReader.cpp
 *
 * Created on 2012-07-17 by Norihiro Watanabe
 */

#include "OptionsXMLReader.h"

#include <map>
#include <vector>

#include "logog.hpp"
#include "tinyxml2.h"

#include "StringTools.h"

namespace BaseLib
{

void addXMLtoOptions(tinyxml2::XMLElement* e_root, BaseLib::OptionGroup &properties)
{
    BaseLib::OptionGroup* optRoot = properties.addSubGroup(e_root->Name());

    // attributes
    for (const tinyxml2::XMLAttribute* att = e_root->FirstAttribute(); att != 0; att = att->Next())
    {
        optRoot->addOption(att->Name(), att->Value());
    }

    // element
    for (tinyxml2::XMLElement* e=e_root->FirstChildElement(); e!=0; e=e->NextSiblingElement())
    {
        if (e->FirstChildElement()!=0 || e->FirstAttribute()!=0) {
            addXMLtoOptions(e, *optRoot);
        } else if (e->GetText() != 0){
            optRoot->addOption(e->Name(), e->GetText());
        }
    }
}

bool addXMLtoOptions(const std::string &xml_file, BaseLib::Options &properties)
{
    tinyxml2::XMLDocument doc;
    int ret = doc.LoadFile(xml_file.c_str());
    if (ret != tinyxml2::XML_NO_ERROR) {
        LOGOG_CERR << "Error in reading a XML file " << xml_file << " with error code " << ret << std::endl;
        return false;
    }

    // read data
    addXMLtoOptions(doc.RootElement(), properties);


    return true;
}

} //end
