/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file CouplingStrucutreBuilder.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <string>

#include "BaseLib/Options.h"

#include "ICoupledProblem.h"
#include "MonolithicProblem.h"
#include "PartitionedProblem.h"
#include "Algorithm/PartitionedAlgorithmFactory.h"


namespace NumLib
{

template <class T_I, class T_M, class T_P, class T_ALGORITHM>
class TemplateCouplingStrucutreBuilder
{
    std::vector<T_M*> _vec_m;
    std::vector<std::string> _vec_m_name;
public:
    TemplateCouplingStrucutreBuilder() {};

    template <class T_EQS_FACTORY>
    T_I* build(const BaseLib::Options *option, T_EQS_FACTORY &eqs_fac)
    {
        const BaseLib::Options* op_cpl = option->getSubGroup("coupling");
        if (op_cpl==0) {
            INFO("tag<coupling> not found.");
            return NULL;
        }
        if (op_cpl->hasSubGroup("M")) {
            const BaseLib::Options* op_sub = op_cpl->getSubGroup("M");
            T_M *sys = buildMonolithicSystem(op_sub, eqs_fac);
            return sys;
        } else if (op_cpl->hasSubGroup("P")) {
            const BaseLib::Options* op_sub = op_cpl->getSubGroup("P");
            T_P *sys = buildPartitionedSystem(op_sub, eqs_fac);
            return sys;
        }
        return NULL;
    }

    std::vector<T_M*>& getListOfMonolithicSystem()
    {
        return _vec_m;
    }

    std::vector<std::string>& getListOfMonolithicSystemName()
    {
        return _vec_m_name;
    }

private:
    template <class T_EQS_FACTORY>
    T_M* buildMonolithicSystem(const BaseLib::Options *option, T_EQS_FACTORY &eqs_fac)
    {
        std::string eqs_type = option->getOption("type");
        if (eqs_type.empty()) {
            ERR("***Error in TemplateCouplingStrucutreBuilder: attribute \"type\" is not found.");
            return NULL;
        }
        T_M* eqs = eqs_fac.create(eqs_type);
        if (eqs==NULL) {
            ERR("***Error in TemplateCouplingStrucutreBuilder: Monolithic system type <%s> is not found.", eqs_type.c_str());
            return NULL;
        }
        if (option->hasOption("in")) {
            std::vector<std::string> in_names = option->getOptionList<std::string>("in");
            for (size_t i=0; i<in_names.size(); i++) {
                eqs->setInputParameterName(i, in_names[i]);
            }
        }
        if (option->hasOption("out")) {
            std::vector<std::string> out_names = option->getOptionList<std::string>("out");
            for (size_t i=0; i<out_names.size(); i++) {
                eqs->setOutputParameterName(i, out_names[i]);
            }
        }
        _vec_m.push_back(eqs);
        std::string eqs_name = option->getOption("name");
        _vec_m_name.push_back(eqs_name);
        return eqs;
    }

    template <class T_EQS_FACTORY>
    T_P* buildPartitionedSystem(const BaseLib::Options *option, T_EQS_FACTORY &eqs_fac)
    {
        T_P* part = new T_P();
        //para
        if (option->hasOption("in")) {
            std::vector<std::string> in_names = option->getOptionList<std::string>("in");
            part->resizeInputParameter(in_names.size());
            for (size_t i=0; i<in_names.size(); i++) {
                part->setInputParameterName(i, in_names[i]);
            }
        }
        if (option->hasOption("out")) {
            std::vector<std::string> out_names = option->getOptionList<std::string>("out");
            part->resizeOutputParameter(out_names.size());
            for (size_t i=0; i<out_names.size(); i++) {
                part->setOutputParameterName(i, out_names[i]);
            }
        }
        //alg
        size_t max_itr = option->getOptionAsNum<size_t>("max_itr");
        double epsilon = option->getOptionAsNum<double>("epsilon");
        part->setAlgorithm(*T_ALGORITHM::create(option->getOption("algorithm"), max_itr, epsilon));
//        IPartitionedAlgorithm* alg = T_ALGORITHM::create(option->getOption("algorithm"), checker, max_itr, epsilon);
//        if (alg!=0) {
//        }
        //problems
        const BaseLib::Options* op_problems = option->getSubGroup("problems");
        for (BaseLib::Options::const_iterator itr=op_problems->begin(); itr!=op_problems->end(); ++itr) {
            std::string str = itr->first;
            BaseLib::Options* op_sub = static_cast<BaseLib::Options*>(itr->second);
            T_I* sys = 0;
            if (str.compare("M")==0) {
                sys = buildMonolithicSystem(op_sub, eqs_fac);
                if (sys!=0) part->addProblem(*sys);
            } else if (str.compare("P")==0) {
                sys = buildPartitionedSystem(op_sub, eqs_fac);
                if (sys!=0) part->addProblem(*sys, true);
            }
        }
        part->connectParameters();
        return part;
    }
};

typedef class TemplateCouplingStrucutreBuilder<ICoupledSystem,TemplateSteadyMonolithicSystem,PartitionedProblem,PartitionedAlgorithmFactory> CouplingStrucutreBuilder;


} //end
