/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemVariableBuilder.cpp
 *
 * Created on 2012-10-19 by Norihiro Watanabe
 */

#include "FemVariableBuilder.h"

#include <vector>
#include <cassert>
#include "logog.hpp"
#include "NumLib/Function/TXFunctionBuilder.h"
#include "SolutionLib/Fem/IFemNeumannBC.h"
#include "SolutionLib/Fem/FemSourceTerm.h"

void FemVariableBuilder::doit(const std::string &given_var_name, const BaseLib::Options &option, const MeshLib::IMesh* msh, const GeoLib::GEOObjects *geo, const std::string &geo_unique_name, FemLib::IFeObjectContainer* _feObjects, SolutionLib::FemVariable* var) const
{
    // IC
    SolutionLib::FemIC* var_ic = new SolutionLib::FemIC(msh);
    const BaseLib::Options* opICList = option.getSubGroup("ICList");
    if (opICList!=nullptr) {
        std::vector<const BaseLib::Options*> vec_opIC = opICList->getSubGroupList("IC");
        for (size_t i=0; i<vec_opIC.size(); i++)
        {
            const BaseLib::Options* opIC = vec_opIC[i];
            std::string var_name = opIC->getOption("Variable");
            if (var_name.compare(given_var_name)!=0) continue;
            std::string geo_type = opIC->getOption("GeometryType");
            std::string geo_name = opIC->getOption("GeometryName");
            const GeoLib::GeoObject* geo_obj = geo->searchGeoByName(geo_unique_name, geo_type, geo_name);
            assert(opIC->hasOption("DistributionType"));
            NumLib::ITXFunction* f_ic = NumLib::TXFunctionBuilder::create(*opIC);
            var_ic->addDistribution(geo_obj, f_ic);
        }
    }
    var->setIC(var_ic);
    // BC
    const BaseLib::Options* opBCList = option.getSubGroup("BCList");
    if (opBCList!=nullptr) {
        std::vector<const BaseLib::Options*> vec_opBC = opBCList->getSubGroupList("BC");
        for (size_t i=0; i<vec_opBC.size(); i++)
        {
            const BaseLib::Options* opBC = vec_opBC[i];
            std::string var_name = opBC->getOption("Variable");
            if (var_name.compare(given_var_name)!=0) continue;
            std::string geo_type = opBC->getOption("GeometryType");
            std::string geo_name = opBC->getOption("GeometryName");
            const GeoLib::GeoObject* geo_obj = geo->searchGeoByName(geo_unique_name, geo_type, geo_name);
            assert(opBC->hasOption("DistributionType"));
            NumLib::ITXFunction* f_bc = NumLib::TXFunctionBuilder::create(*opBC);
            var->addDirichletBC(new SolutionLib::FemDirichletBC(msh, geo_obj, f_bc));
        }
    }
    // ST
    const BaseLib::Options* opSTList = option.getSubGroup("STList");
    if (opSTList != nullptr) {
        std::vector<const BaseLib::Options*> vec_opST = opSTList->getSubGroupList("ST");
        for (size_t i=0; i<vec_opST.size(); i++)
        {
            const BaseLib::Options* opST = vec_opST[i];
            std::string var_name = opST->getOption("Variable");
            if (var_name.compare(given_var_name)!=0) continue;
            std::string geo_type = opST->getOption("GeometryType");
            std::string geo_name = opST->getOption("GeometryName");
            const GeoLib::GeoObject* geo_obj = geo->searchGeoByName(geo_unique_name, geo_type, geo_name);
            std::string st_type = opST->getOption("STType");
            assert(opST->hasOption("DistributionType"));
            NumLib::ITXFunction* f_st = NumLib::TXFunctionBuilder::create(*opST);
            if (st_type.compare("NEUMANN")==0) {
                // user set inflow as positive sign but internally negative
                NumLib::ITXFunction* f_st_old = f_st;
                NumLib::TXFunctionConstant constantFunction(-1.);
                f_st = new NumLib::TXCompositFunction
                <
                    NumLib::ITXFunction, NumLib::TXFunctionConstant,
                    NumLib::Multiplication
                >(f_st, &constantFunction);
                delete f_st_old;
            }
            SolutionLib::IFemNeumannBC *femSt = 0;
            if (st_type.compare("NEUMANN")==0) {
                femSt = new SolutionLib::FemNeumannBC(msh, _feObjects, geo_obj, f_st);
            } else if (st_type.compare("SOURCESINK")==0) {
                femSt = new SolutionLib::FemSourceTerm(msh, geo_obj, f_st);
            }
            delete f_st;
            var->addNeumannBC(femSt);
        }
    }
}
