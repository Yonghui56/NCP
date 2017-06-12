/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemLinearEQS.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "NumLib/Function/IFunction.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "NumLib/TransientAssembler/TransientElementWiseLinearEquationUpdater.h"
#include "FemLib/Function/FemFunction.h"
#include "FemVariable.h"
#include "FemDirichletBC.h"
#include "IFemNeumannBC.h"

#include "SolutionLib/DataType.h"

namespace SolutionLib
{

/**
 * \brief Template class for transient linear FEM functions
 *
 * - FEM variables
 * - Local assembler
 * - Linear EQS
 * - Current time step
 * - Previous time step value
 *
 * \tparam T_ASSEMBLER
 */
template <
    class T_DIS_SYS,
    class T_LINEAR_SOLVER,
    class T_LOCAL_ASSEMBLER
    >
class TemplateTransientLinearFEMFunction
    : public NumLib::TemplateFunction<SolutionVector, SolutionVector>
{
public:
    typedef T_DIS_SYS MyDiscreteSystem;
    typedef T_LINEAR_SOLVER LinearSolverType;
    typedef FemVariable MyVariable;
    typedef T_LOCAL_ASSEMBLER UserLocalAssembler;
    typedef typename NumLib::TransientElementWiseLinearEquationUpdater<UserLocalAssembler> MyUpdater;
    typedef typename T_DIS_SYS::template MyLinearEquationAssembler<MyUpdater,LinearSolverType>::type MyGlobalAssembler;

    /// constructor
    /// \param list_var
    /// \param assembler
    /// \param linear_eqs    Linear equation
    TemplateTransientLinearFEMFunction(MyDiscreteSystem* sys, const std::vector<MyVariable*> &list_var, UserLocalAssembler* assembler, DiscreteLib::IDiscreteLinearEquation* linear_eqs)
        : _sys(sys), _local_assembler(assembler),  _linear_eqs(linear_eqs),
          _t_n1(0), _u_n0(0), _list_var(list_var)
    {
    };

    ///
    virtual ~TemplateTransientLinearFEMFunction() {};

    ///
    NumLib::TemplateFunction<SolutionVector,SolutionVector>* clone() const
    {
        return new TemplateTransientLinearFEMFunction<
                MyDiscreteSystem,
                LinearSolverType,
                UserLocalAssembler
                    >(_sys, _list_var, _local_assembler, _linear_eqs);
    }

    /// reset property
    void reset(const NumLib::TimeStep* t, const SolutionVector* u_n0)
    {
        this->_t_n1 = t;
        this->_u_n0 = u_n0;
    };

    /// solve linear equations discretized with FEM
    /// @param u_k0    initial guess
    /// @param u_k1    new results
    void eval(const SolutionVector &u_k0, SolutionVector &u_k1);

private:
    MyDiscreteSystem* _sys;
    UserLocalAssembler* _local_assembler;
    DiscreteLib::IDiscreteLinearEquation* _linear_eqs;
    const NumLib::TimeStep* _t_n1;
    const SolutionVector* _u_n0;
    std::vector<MyVariable*> _list_var;
};


template <class T1, class T2, class T3>
void TemplateTransientLinearFEMFunction<T1, T2, T3>::eval(const SolutionVector &u_k0, SolutionVector &u_k1)
{
    // input, output
    const NumLib::TimeStep &t_n1 = *this->_t_n1;
    const SolutionVector* u_n = this->_u_n0;

    _linear_eqs->initialize();

    // setup BC
    for (size_t i=0; i<_list_var.size(); i++) {
        MyVariable* var = _list_var[i];
        std::vector<size_t> var_bc_id;
        std::vector<double> var_bc_val;
        for (size_t j=0; j<var->getNumberOfDirichletBC(); j++) {
            FemDirichletBC* bc1 = var->getDirichletBC(j);
            bc1->setup(var->getCurrentOrder());
            var_bc_id.insert(var_bc_id.end(), bc1->getListOfBCNodes().begin(), bc1->getListOfBCNodes().end());
            var_bc_val.insert(var_bc_val.end(), bc1->getListOfBCValues().begin(), bc1->getListOfBCValues().end());
        }
        _linear_eqs->setPrescribedDoF(i, var_bc_id, var_bc_val);
        for (size_t j=0; j<var->getNumberOfNeumannBC(); j++) {
            IFemNeumannBC* st = var->getNeumannBC(j);
            st->initCurrentTime(t_n1.getTime());
            st->setup(var->getCurrentOrder());
        }
    }

    // assembly
    MeshLib::IMesh* msh = _sys->getMesh();
    MyUpdater updater(&t_n1, msh, u_n, &u_k0, _local_assembler);
    MyGlobalAssembler assembler(&updater);
    _linear_eqs->construct(assembler);

    //apply BC1,2
    for (size_t i=0; i<_list_var.size(); i++) {
        MyVariable* var = _list_var[i];
        for (size_t j=0; j<var->getNumberOfNeumannBC(); j++) {
            IFemNeumannBC* bc2 = var->getNeumannBC(j);
            _linear_eqs->addRHS(i, bc2->getListOfBCNodes(), bc2->getListOfBCValues(), -1.0);
        }
    }

    // solve
    _linear_eqs->setX(u_k0);
    _linear_eqs->solve();
    _linear_eqs->getX(u_k1);
}

} //end
