/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file CouplingStrucutreBuilder4SysEqs.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <string>

#include "BaseLib/Options.h"

#include "NumLib/Coupling/ICoupledProblem.h"
#include "NumLib/Coupling/MonolithicProblem.h"
#include "NumLib/Coupling/PartitionedProblem.h"
#include "NumLib/Coupling/Algorithm/PartitionedAlgorithmFactory.h"
#include "MyCouplingEQS.h"

namespace NumLib
{

template <class T_I, class T_M, class T_P, class T_ALGORITHM>
class TemplateCouplingStrucutreBuilder4SysEqs
{
public:
    template <class T_EQS_FACTORY>
    T_I* build(const BaseLib::Options *option, T_EQS_FACTORY eqs_fac)
    {
        const BaseLib::Options* op_cpl = option->getSubGroup("coupling");
        if (op_cpl->begin()==op_cpl->end()) return 0;

        std::string str = op_cpl->begin()->first;
        BaseLib::Options* op_sub = static_cast<BaseLib::Options*>(op_cpl->begin()->second);
        if (str.find("M")==0) {
            //const Base::Options* op_sub = op_cpl->getSubGroup("M");
            T_M *sys = buildMonolithicSystem(op_sub, eqs_fac);
            return sys;
        } else if (str.find("P")==0) {
            //const Base::Options* op_sub = op_cpl->getSubGroup("P");
            T_P *sys = buildPartitionedSystem(op_sub, eqs_fac);
            return sys;
        }

        //if (op_cpl->hasSubGroup("M")) {
        //    const Base::Options* op_sub = op_cpl->getSubGroup("M");
        //    T_M *sys = buildMonolithicSystem(op_sub, eqs_fac);
        //    return sys;
        //} else if (op_cpl->hasSubGroup("P")) {
        //    const Base::Options* op_sub = op_cpl->getSubGroup("P");
        //    T_P *sys = buildPartitionedSystem(op_sub, eqs_fac, check_fack);
        //    return sys;
        //}
        return 0;
    }

private:
    template <class T_EQS_FACTORY>
    T_M* buildMonolithicSystem(const BaseLib::Options *option, T_EQS_FACTORY &eqs_fac)
    {
        std::vector<std::string> list_var_name = option->getOptionList<std::string>("variable");
        T_M* eqs = eqs_fac.create(list_var_name);
        return eqs;
    }

    template <class T_EQS_FACTORY>
    T_P* buildPartitionedSystem(const BaseLib::Options *option, T_EQS_FACTORY &eqs_fac)
    {
        T_P* part = new T_P();
        //alg
        size_t max_itr = option->getOptionAsNum<size_t>("max_itr");
        double epsilon = option->getOptionAsNum<double>("epsilon");
        //IConvergenceCheck* checker = check_fac.create(option->getOption("convergence"));
        part->setAlgorithm(*T_ALGORITHM::create(option->getOption("algorithm"), max_itr, epsilon));
        //problems
        const BaseLib::Options* op_problems = option->getSubGroup("problems");
        std::vector<T_I*> list_subproblems;
        for (BaseLib::Options::const_iterator itr=op_problems->begin(); itr!=op_problems->end(); ++itr) {
            std::string str = itr->first;
            BaseLib::Options* op_sub = static_cast<BaseLib::Options*>(itr->second);
            T_I* sys = 0;
            if (str.find("M")==0) {
                sys = buildMonolithicSystem(op_sub, eqs_fac);
            } else if (str.compare("P")==0) {
                sys = buildPartitionedSystem(op_sub, eqs_fac);
            }
            if (sys!=0) {
                part->addProblem(*sys);
            }
        }
        eqs_fac.setPartitionedProblem(part);
        return part;
    }
};

template <class T_CONVERGENCE>
struct CouplingStrucutreBuilder4SysEqs
{
    typedef TemplateCouplingStrucutreBuilder4SysEqs<ICoupledSystem,MyCouplingEQS<T_CONVERGENCE>,PartitionedProblem,PartitionedAlgorithmFactory> type;

};

//typedef class TemplateCouplingStrucutreBuilder4SysEqs<ICoupledSystem,MyCouplingEQS,PartitionedProblem,PartitionedAlgorithmFactory> CouplingStrucutreBuilder4SysEqs;


} //end
