/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemDxEQS.h
 *
 * Created on 2013-08-16 by Yonghui Huang & Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <fstream>


#include "MathLib/DataType.h"
#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "NumLib/Function/IFunction.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "NumLib/TransientAssembler/TransientElementWiseMatrixUpdater.h"
#include "FemLib/Function/FemFunction.h"
#include "SolutionLib/DataType.h"
#include "SolutionLib/Fem/FemDirichletBC.h"
#include "SolutionLib/Fem/FemNeumannBC.h"
#include "SolutionLib/Fem/FemVariable.h"
#include "IterativeLinearSolvers"
#include "time.h"
/**
 * \brief Template class for transient linear FEM functions
 *
 * \tparam T_SPACE_ASSEMBLER
 */
template <class T_DIS_SYS,
          class T_LINEAR_SOLVER,
          class T_FUNCTION_DATA,
		  class T_LOCAL_JACOBIAN_ASSEMBLER>
class TemplateTransientDxFEMFunction_Global_Complementary
{
public:
    typedef T_DIS_SYS MyDiscreteSystem;
    typedef typename T_DIS_SYS::template MyLinearEquation<T_LINEAR_SOLVER,DiscreteLib::SparsityBuilderFromNodeConnectivity>::type MyLinaerEQS;
    typedef SolutionLib::FemVariable MyVariable;
    typedef T_LINEAR_SOLVER LinearSolverType;
	typedef T_LOCAL_JACOBIAN_ASSEMBLER UserLocalJacobianAssembler;
    typedef typename FemLib::FemNodalFunctionScalar<MyDiscreteSystem>::type MyNodalFunctionScalar;
	typedef typename NumLib::TransientElementWiseMatrixUpdater<UserLocalJacobianAssembler> MyUpdater;
	typedef typename T_DIS_SYS::template MyLinearEquationAssembler<MyUpdater, LinearSolverType>::type MyGlobalAssembler;
    /// constructor
    /// @param problem        Fem problem
    /// @param linear_eqs    Discrete linear equation
	TemplateTransientDxFEMFunction_Global_Complementary(MeshLib::IMesh* msh,
    		                                  std::vector<MyVariable*> &list_var,
    		                                  MyLinaerEQS* linear_eqs,
    		                                  DiscreteLib::DofEquationIdTable* dofManager,
											  T_FUNCTION_DATA* func_data,
											  UserLocalJacobianAssembler* asssembler)
        : _msh(msh),  _linear_eqs(linear_eqs),
		_t_n1(0), _u_n0(0), _list_var(list_var), _dofManager(dofManager), _function_data(func_data),
		_local_assembler(asssembler)
    {
    };

    ///
	virtual ~TemplateTransientDxFEMFunction_Global_Complementary()
    {
		_function_data = NULL;
    };

    ///
    NumLib::TemplateFunction<SolutionLib::SolutionVector,SolutionLib::SolutionVector>* clone() const
    {
        return new TemplateTransientDxFEMFunction_Global_Complementary<
                    T_DIS_SYS,
                    T_LINEAR_SOLVER,
                    T_FUNCTION_DATA,
					T_LOCAL_JACOBIAN_ASSEMBLER
		>(_list_var, _linear_eqs, _dofManager, _function_data, _local_assembler);
    }

    /// reset property
    void reset(const NumLib::TimeStep* t, const SolutionLib::SolutionVector* u_n0)
    {
        this->_t_n1 = t;
        this->_u_n0 = u_n0;
    };

    /// solve Newton-Raphson equations
    /// @param u_n1    current solution
    /// @param r       residual
    /// @param du      increment
    void eval(const SolutionLib::SolutionVector &u_n1, const SolutionLib::SolutionVector &r, SolutionLib::SolutionVector &du);

private:
	void GlobalJacobianAssembler(T_FUNCTION_DATA* _function_data,
    							 LinearSolverType & eqsJacobian_global);
	

private:
    MeshLib::IMesh* _msh;
    DiscreteLib::DofEquationIdTable* _dofManager;
    MyLinaerEQS* _linear_eqs;
    const NumLib::TimeStep* _t_n1;
    const SolutionLib::SolutionVector* _u_n0;
    std::vector<MyVariable*> _list_var;
	T_FUNCTION_DATA* _function_data;
    //MyDiscreteSystem* _dis;
    FemLib::IFiniteElement* _fe;
	UserLocalJacobianAssembler *_local_assembler;

	
};


template <class T1, class T2, class T3, class T4>
void TemplateTransientDxFEMFunction_Global_Complementary<T1, T2, T3, T4>::eval(const SolutionLib::SolutionVector &u_n1, const SolutionLib::SolutionVector &r, SolutionLib::SolutionVector &du)
{
	
    // input, output
    const NumLib::TimeStep &t_n1 = *this->_t_n1;
	const SolutionLib::SolutionVector* u_n = this->_u_n0;
	size_t num_node = _msh->getNumberOfNodes();
	size_t num_var = _list_var.size();
	// const SolutionLib::SolutionVector* u_n = this->_u_n0;

    _linear_eqs->initialize();

    // setup BC1 for solution increment
    for (size_t i=0; i<_list_var.size(); i++) {
        MyVariable* var = _list_var[i];
        std::vector<size_t> var_bc_id;
        std::vector<double> var_bc_val;
        for (size_t j=0; j<var->getNumberOfDirichletBC(); j++) {
            SolutionLib::FemDirichletBC* bc1 = var->getDirichletBC(j);
            bc1->setup(var->getCurrentOrder());
            std::vector<double> bc_value_for_dx(bc1->getListOfBCNodes().size(), .0);
            for (size_t k=0; k<bc_value_for_dx.size(); k++) {
                size_t id = bc1->getListOfBCNodes()[k];
                double v =  bc1->getListOfBCValues()[k];
                size_t eqs_id = _linear_eqs->getDofMapManger()->mapEqsID(i, _msh->getID(), id);
                bc_value_for_dx[k] = v - u_n1[eqs_id]; // dx = (bc value) - (current value)
            }
            var_bc_id.insert(var_bc_id.end(), bc1->getListOfBCNodes().begin(), bc1->getListOfBCNodes().end());
            var_bc_val.insert(var_bc_val.end(), bc_value_for_dx.begin(), bc_value_for_dx.end());

        }
        _linear_eqs->setPrescribedDoF(i, var_bc_id, var_bc_val);
    }
	
    LinearSolverType* eqsJacobian_global = _linear_eqs->getLinearSolver();
    // Assembling global jacobian
	MyUpdater updater(&t_n1, _msh, u_n, &u_n1, _local_assembler);
	//_local_assembler->get_function_data()->get_vec_Res();
	MyGlobalAssembler assembler(&updater);
	_linear_eqs->construct(assembler);
	//mat_Jacobian = _linear_eqs.getA();
    GlobalJacobianAssembler(_function_data, *eqsJacobian_global);
    // set residual
    _linear_eqs->addRHS(r, -1.0);

    // solve
	_linear_eqs->solve();
    _linear_eqs->getX(du);

}

template <class T_DIS_SYS, class T_LINEAR_SOLVER, class T_FUNCTION_DATA, class T_LOCAL_JACOBIAN_ASSEMBLER>
void TemplateTransientDxFEMFunction_Global_Complementary<T_DIS_SYS, T_LINEAR_SOLVER, T_FUNCTION_DATA, T_LOCAL_JACOBIAN_ASSEMBLER>::GlobalJacobianAssembler(T_FUNCTION_DATA* _function_data, LinearSolverType & eqsJacobian_global)
{
	const size_t nnodes = _msh->getNumberOfNodes();
	auto& mat_Jacobian = _function_data->get_mat_Jacob();
	for (size_t node_idx_row = 0; node_idx_row < nnodes; node_idx_row++){
		for (size_t node_idx_col = 0; node_idx_col < 3*nnodes; node_idx_col++){
			double const tmp =
			    mat_Jacobian.coeff(2 * nnodes + node_idx_row, node_idx_col);
			eqsJacobian_global.addA(3*node_idx_row+2, node_idx_col, -tmp);
		}
	}
}


