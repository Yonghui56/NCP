/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TemplateDiscreteNonlinearSolver.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <string>

#include "BaseLib/Options.h"
#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "DiscreteLib/Core/IDiscreteVector.h"
#include "NonlinearSolver.h"
#include "NonlinearSolverOption.h"
#include "DiscreteNonlinearSolverFactory.h"

namespace NumLib
{

template <
    class T_DIS_SYS, 
    class F_LINEAR, 
    class F_R, 
    class F_DX, 
    class T_NL_FACTORY=DiscreteNonlinearSolverFactory
>
class TemplateDiscreteNonlinearSolver
{
public:
    typedef T_DIS_SYS MyDiscreteSystem;
    typedef DiscreteLib::IDiscreteVector<double> VectorType;

    TemplateDiscreteNonlinearSolver(MyDiscreteSystem* dis_sys, F_LINEAR* f_l, F_R* f_r, F_DX* f_dx)
    : _dis_sys(dis_sys), _f_l(f_l), _f_r(f_r), _f_dx(f_dx), _solver(NULL), _nl_factory(new T_NL_FACTORY), _log(nullptr)
    {
    }

    TemplateDiscreteNonlinearSolver(MyDiscreteSystem* dis_sys, F_LINEAR* f_l, F_R* f_r, F_DX* f_dx, T_NL_FACTORY* nl_factory)
    : _dis_sys(dis_sys), _f_l(f_l), _f_r(f_r), _f_dx(f_dx), _solver(NULL), _nl_factory(nl_factory), _log(nullptr)
    {
    }


    virtual ~TemplateDiscreteNonlinearSolver()
    {
        BaseLib::releaseObject(_solver, _nl_factory);
    }

    void setOption(const BaseLib::Options &option)
    {
        const BaseLib::Options *op = option.getSubGroup("NonlinearSolver");
        if (op==0) return;

        if (op->hasOption("solver_type"))
            _option.solver_type = _option.getSolverType(op->getOption("solver_type"));
        if (op->hasOption("error_tolerance"))
            _option.error_tolerance = op->getOptionAsNum<double>("error_tolerance");
        if (op->hasOption("max_iteration_step"))
            _option.max_iteration = op->getOptionAsNum<int>("max_iteration_step");
    }

    void setOption(const NonlinerSolverOption &option) { _option = option; }

    NonlinerSolverOption &getOption() const { return _option; }

    void setLog(BaseLib::Options &log) { _log = &log; }

    const BaseLib::Options* getLog() const {return _log;};

    void solve(const VectorType &x_0, VectorType &x_new)
    {
        if (_solver==0)
            _solver = _nl_factory->create(_option, _dis_sys, _f_l, _f_r, _f_dx);
        _solver->solve(x_0, x_new);

        if (_log!=nullptr)
            _solver->recordLog(*_log);
    }

private:
    DISALLOW_COPY_AND_ASSIGN(TemplateDiscreteNonlinearSolver);

private:
    NonlinerSolverOption _option;
    MyDiscreteSystem* _dis_sys;
    F_LINEAR* _f_l;
    F_R* _f_r;
    F_DX* _f_dx;
    INonlinearSolver* _solver;
    T_NL_FACTORY* _nl_factory;
    BaseLib::Options* _log;
};

} //end

