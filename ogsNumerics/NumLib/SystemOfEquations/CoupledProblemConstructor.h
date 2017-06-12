/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file CoupledProblemConstructor.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "NumLib/Coupling/PartitionedProblem.h"
#include "SystemOfEquations.h"
#include "MyCouplingEQS.h"

namespace NumLib
{

template <class T_CONVERGENCE>
class CoupledProblemConstructor
{
public:
    ICoupledSystem* build(SystemOfEquations &sys, std::vector<std::vector<Variable*> > &list_active_vars, std::vector<IPartitionedAlgorithm* > &list_part_alg, std::vector<MyCouplingEQS<T_CONVERGENCE>*> &list_sub_problem)
    {
        for (size_t i=0; i<list_active_vars.size(); i++) {
            MyCouplingEQS<T_CONVERGENCE>* eqs1 = createPartitionedProblem(sys, list_active_vars[i]);
            list_sub_problem.push_back(eqs1);
        }

        if (list_sub_problem.size()==1) {
            return list_sub_problem[0];
        } else {
            PartitionedProblem *part1 = new PartitionedProblem();
            part1->resizeOutputParameter(sys.getNumberOfVariables());
            for (size_t i=0; i<sys.getNumberOfVariables(); i++)
                part1->setOutputParameterName(i, sys.getVariable(i)->name);
            for (size_t i=0; i<list_sub_problem.size(); i++)
                part1->addProblem(*list_sub_problem[i]);
            part1->connectParameters();
            part1->setAlgorithm(*list_part_alg[0]);

            return part1;
        }
    }

    MyCouplingEQS<T_CONVERGENCE>* createPartitionedProblem(SystemOfEquations &sys, const std::vector<Variable*> &active_vars)
    {
        std::vector<LinearEquation*> active_eqs;
        for (size_t i=0; i<active_vars.size(); i++) {
            for (size_t j=0; j<sys.getNumberOfEquations(); j++) {
                if (active_vars[i]->id == sys.getEquation(j)->primary_variable->id)
                    active_eqs.push_back(sys.getEquation(j));
            }
        }
        MyCouplingEQS<T_CONVERGENCE>* eqs = new MyCouplingEQS<T_CONVERGENCE>(active_vars, active_eqs);
        return eqs;
    }

    MyCouplingEQS<T_CONVERGENCE>* createMonolithicProblem(SystemOfEquations &sys)
    {
        MyCouplingEQS<T_CONVERGENCE>* eqs = new MyCouplingEQS<T_CONVERGENCE>(*sys.getListOfVariables(), *sys.getListOfEquations());
        return eqs;
    }
};

template <class T_CONVERGENCE>
class CoupledProblemFactory
{
    SystemOfEquations* _sys;
    std::vector<MyCouplingEQS<T_CONVERGENCE>*>* _list_sub_problem;
    std::vector<typename MyCouplingEQS<T_CONVERGENCE>::ArrayType*>* _ini_para;
public:
    CoupledProblemFactory(SystemOfEquations &sys, std::vector<typename MyCouplingEQS<T_CONVERGENCE>::ArrayType*> &ini_para, std::vector<MyCouplingEQS<T_CONVERGENCE>*> &list_sub_problem)
    {
        _sys = &sys;
        _list_sub_problem = &list_sub_problem;
        _ini_para = &ini_para;
    }

    // this function should be called after sub problems are added
    void setPartitionedProblem(PartitionedProblem* part)
    {
        // find a list of input/output parameter names
        std::set<std::string> set_out_names;
        std::set<std::string> set_in_names;
        std::vector<std::string> out_names;
        for (size_t i=0; i<part->getNumberOfSubProblems(); i++) {
            ICoupledSystem* prob = part->getProblem(i);
            for (size_t j=0; j<prob->getNumberOfOutputParameterNames(); j++) {
                const std::string &name = prob->getOutputParameterName(j);
                if (set_out_names.count(name)==0) {
                    set_out_names.insert(name);
                    out_names.push_back(name);
                }
            }
            for (size_t j=0; j<prob->getNumberOfInputParameterNames(); j++) {
                const std::string &name = prob->getInputParameterName(j);
                set_in_names.insert(name);
            }
        }

        // select input parameters which are not output of any sub problems
        std::vector<std::string> in_names;
        for (std::set<std::string>::iterator itr=set_in_names.begin(); itr!=set_in_names.end(); ++itr)
        {
            if (set_out_names.count(*itr)==0)
                in_names.push_back(*itr);
        }

        // configure
        part->resizeInputParameter(in_names.size());
        for (size_t i=0; i<in_names.size(); i++)
            part->setInputParameterName(i, in_names[i]);
        part->resizeOutputParameter(out_names.size());
        for (size_t i=0; i<out_names.size(); i++)
            part->setOutputParameterName(i, out_names[i]);

        part->connectParameters();
    }

    MyCouplingEQS<T_CONVERGENCE>* create(const std::vector<std::string> &list_var_name)
    {
        std::vector<Variable*> list_active_vars;
        for (size_t i=0; i<list_var_name.size(); i++) {
            for (size_t j=0; j<_sys->getNumberOfVariables(); j++) {
                Variable* var = _sys->getVariable(j);
                if (var->name.compare(list_var_name[i])==0) {
                    list_active_vars.push_back(var);
                }
            }
        }

        std::vector<LinearEquation*> active_eqs;
        for (size_t i=0; i<list_active_vars.size(); i++) {
            for (size_t j=0; j<_sys->getNumberOfEquations(); j++) {
                if (list_active_vars[i]->id == _sys->getEquation(j)->primary_variable->id)
                    active_eqs.push_back(_sys->getEquation(j));
            }
        }
        MyCouplingEQS<T_CONVERGENCE>* eqs = new MyCouplingEQS<T_CONVERGENCE>(list_active_vars, active_eqs);
        eqs->setInitial(*_ini_para);
        _list_sub_problem->push_back(eqs);

        return eqs;
    }
};
} //end
