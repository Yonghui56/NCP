/**
* Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.com/LICENSE.txt
*
*
* \file SingleStepCMP_PressureForm.h
*
* Created on 2015-03-12 by Yonghui HUANG
*/

#ifndef SINGLE_STEP_CMP_PRESSURE_FORM_H 
#define SINGLE_STEP_CMP_PRESSURE_FORM_H 

#include <vector>
#include <map>
#include "logog.hpp"

#include "BaseLib/CodingTools.h"
#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "DiscreteLib/Core/IDiscreteLinearEquation.h"
#include "DiscreteLib/Utils/SparsityBuilderFromNodeConnectivity.h"
#include "FemLib/Function/FemFunction.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/Nonlinear/TemplateDiscreteNonlinearSolver.h"

#include "SolutionLib/DataType.h"
#include "SolutionLib/Core/AbstractTimeSteppingAlgorithm.h"
#include "SolutionLib/Fem/TemplateFemEquation.h"
#include "SolutionLib/Fem/FemDirichletBC.h"
#include "SolutionLib/Fem/FemNeumannBC.h"


namespace SolutionLib
{

/**
* \brief Solution algorithm for kinetic reaction reduction problems using FEM with single time stepping method
*
* - previous and current time step values
* - ST values
* - Linear equation
* - Residual equation
* - Dx equation
* - Linear EQS
* - Linear solver
*
* \tparam T_USER_FEM_PROBLEM      the FEM problem class
* \tparam T_LINEAR_SOLVER         Linear equation solver class
*/
template <
	class T_USER_FUNCTION_DATA,
	class T_USER_FEM_PROBLEM,
	class T_USER_NON_LINEAR_PROBLEM,
	class T_USER_NON_LINEAR_SOLUTION
>
class SingleStepCMP_PressureForm
	: public AbstractTimeSteppingAlgorithm
{
public:

	typedef T_USER_FUNCTION_DATA       UserFunctionData;
	typedef T_USER_FEM_PROBLEM         UserFemProblem;
	typedef T_USER_NON_LINEAR_PROBLEM  UserNonLinearProblem;
	typedef T_USER_NON_LINEAR_SOLUTION UserNonLinearSolution;

	typedef typename UserFemProblem::MyDiscreteSystem MyDiscreteSystem;
	typedef typename UserFemProblem::MyVariable MyVariable;
	typedef typename FemLib::FemNodalFunctionScalar<MyDiscreteSystem>::type MyNodalFunctionScalar;

	/**
	* constructor, including the following tasks
	* - initialize solution vectors
	* - set up DoFs and equation index
	* - prepare linear and nonlinear equations
	* - prepare linear and nonlinear solutions
	*/
	SingleStepCMP_PressureForm(MyDiscreteSystem*      dis,
                             UserFemProblem*        problem,
							 UserNonLinearSolution* nlin_solution, 
                             UserFunctionData*      function_data);

	/**
	* destructor, reclaiming the memory
	*/
	virtual ~SingleStepCMP_PressureForm()
	{
		/*
		BaseLib::releaseObject(_problem);

		BaseLib::releaseObjectsInStdVector(_vec_u_n1);
		_discrete_system->deleteVector(_x_n0);
		_discrete_system->deleteVector(_x_n1);
		_discrete_system->deleteVector(_x_n1_0);
		_discrete_system->deleteVector(_x_st);
		BaseLib::releaseObject(_discrete_system);

		BaseLib::releaseObject(_feObjects);
		BaseLib::releaseObjectsInStdVector(_linear_problem);
		BaseLib::releaseObjectsInStdVector(_lin_solutions);
		BaseLib::releaseObject(_function_data);

		BaseLib::releaseObjectsInStdMap(_bc_info);
		*/
	}

	/**
	* solve the time step
	*/
	int solveTimeStep(const NumLib::TimeStep &t_n1);

	/**
	* called when this solution is accepted
	*/
	virtual bool accept(const NumLib::TimeStep &t)
	{
		// this solution itself
		if (!AbstractTimeSteppingAlgorithm::accept(t)) return false;

		return true;
	}

	/**
	* called when this solution is accepted
	*/
	virtual void finalizeTimeStep(const NumLib::TimeStep &t)
	{
		// this solution itself
		AbstractTimeSteppingAlgorithm::finalizeTimeStep(t);
		*_x_n0 = *_x_n1; //copy current value to previous value
		_nlin_solution->finalizeTimeStep(t);
	}

	/**
	* get the corresponding solution pointer
	*/
	MyNodalFunctionScalar* getCurrentSolution(size_t var_id) { return _vec_u_n1[var_id]; }

	/**
	* get the pointer to the problem class
	*/
	UserFemProblem* getProblem() { return _problem; };

	/**
	* get the domain of freedom equation ID table
	*/
	DiscreteLib::DofEquationIdTable* getDofEquationIdTable() { return &_dofManager; };

private:
	DISALLOW_COPY_AND_ASSIGN(SingleStepCMP_PressureForm);

private:
	/**
	* DOF table
	*/
	DiscreteLib::DofEquationIdTable _dofManager;

	/**
	* pointer to the problem
	*/
	UserFemProblem* _problem;

	/**
	  * pointer to non-linear solution
	  */
	UserNonLinearSolution*           _nlin_solution;

	/**
	  * pointer to discretization
	  */
	MyDiscreteSystem *_discrete_system;


	std::vector<MyNodalFunctionScalar*> _vec_u_n1;
	SolutionVector *_x_n0;
	SolutionVector *_x_n1_0;
	SolutionVector *_x_n1;
	SolutionVector *_x_st;

	/**
	* pointer to FEM object
	*/
	FemLib::LagrangeFeObjectContainer* _feObjects;

	/**
	* pointer to the function concentrations
	*/
	UserFunctionData*                _function_data;

};

template <
	class T_USER_FUNCTION_DATA,
	class T_USER_FEM_PROBLEM,
	class T_USER_NON_LINEAR_PROBLEM,
	class T_USER_NON_LINEAR_SOLUTION
>
SingleStepCMP_PressureForm<T_USER_FUNCTION_DATA, T_USER_FEM_PROBLEM, T_USER_NON_LINEAR_PROBLEM, T_USER_NON_LINEAR_SOLUTION>
::SingleStepCMP_PressureForm(MyDiscreteSystem*                dis,
                           UserFemProblem*                  problem,
						   UserNonLinearSolution*           nlin_solution,
						   UserFunctionData*                function_data)
: AbstractTimeSteppingAlgorithm(*problem->getTimeSteppingFunction()),
_problem(problem), _nlin_solution(nlin_solution), _discrete_system(dis), _function_data(function_data)
{

	INFO("->Setting up a solution algorithm for CMP_PressureForm problem.");

	size_t i, j, i_var;
	const size_t n_var = problem->getNumberOfVariables();
	MeshLib::IMesh* msh = dis->getMesh();

	// create dof map
	for (size_t i = 0; i<n_var; i++) {
		MyVariable* var = problem->getVariable(i);
		size_t n_dof_per_var = msh->getNumberOfNodes(var->getCurrentOrder());
		_dofManager.addVariableDoFs(msh->getID(), 0, n_dof_per_var);
		INFO("* Variable %d: name=%s, order=%d, n_dof=%d", i, var->getName().c_str(), var->getCurrentOrder(), n_dof_per_var);
	}
	_dofManager.construct();
	const size_t n_total_dofs = _dofManager.getTotalNumberOfActiveDoFs();
	INFO("* Total number of DoFs = %d", n_total_dofs);

	_feObjects = new FemLib::LagrangeFeObjectContainer(msh);

	// set up initial condition
	_vec_u_n1.resize(n_var, 0);
	for (size_t i = 0; i<n_var; i++) {
		MyVariable* femVar = problem->getVariable(i);
		FemIC* femIC = femVar->getIC();
		MyNodalFunctionScalar* u0 = new MyNodalFunctionScalar();
		u0->initialize(*dis, femVar->getCurrentOrder(), 0);
		//u0->setFeObjectContainer(_feObjects);
		femIC->setup(*u0);
		u0->setFeObjectContainer(_feObjects);
		_vec_u_n1[i] = u0; //u0->clone()
	}

	// initialize vectors for solution and ST
	_x_n0 = dis->template createVector<double>(n_total_dofs);
	_x_n1 = dis->template createVector<double>(n_total_dofs);
	_x_n1_0 = dis->template createVector<double>(n_total_dofs);
	_x_st = dis->template createVector<double>(n_total_dofs);

	// copy values of each variable to one solution vector
	for (i = 0; i<n_var; i++) {
		SolutionVector* vec_var = _vec_u_n1[i]->getDiscreteData();
		DiscreteLib::setGlobalVector(_dofManager, i, msh->getID(), *vec_var, *_x_n0);
	}

	// setup functions
	std::vector<MyVariable*> list_var(n_var);
	for (i = 0; i<n_var; i++)
		list_var[i] = problem->getVariable(i);

	// getting the boundary conditions of concentrations for all components, 
	//const size_t msh_id = _discrete_system->getMesh()->getID();
	std::vector<size_t> list_bc1_eqs_id;
	std::vector<double> list_bc1_val;
	for (i_var = 0; i_var < n_var; i_var++) {
		MyVariable* var = _problem->getVariable(i_var);
		for (i = 0; i < var->getNumberOfDirichletBC(); i++) {
			SolutionLib::FemDirichletBC *bc1 = var->getDirichletBC(i);
			bc1->setup(var->getCurrentOrder());
			std::vector<size_t> &list_bc_nodes = bc1->getListOfBCNodes();
			std::vector<double> &list_bc_values = bc1->getListOfBCValues();

			// now loop over this vector
			for (j = 0; j < list_bc_nodes.size(); j++)
			{
				size_t node_id = list_bc_nodes[j];
				double node_value = list_bc_values[j];
				
				// set BC values
				if (i_var == 0)
					_function_data->set_P_node_values(node_id, node_value);
				else if (i_var == 1)
					_function_data->set_X_node_values(node_id, node_value);

			}  // end of for j
		}  // end of for i
	}  // end of for i_var

};

template <class T_USER_FUNCTION_DATA,
          class T_USER_FEM_PROBLEM,
		  class T_USER_NON_LINEAR_PROBLEM,
		  class T_USER_NON_LINEAR_SOLUTION>
int SingleStepCMP_PressureForm<T_USER_FUNCTION_DATA, T_USER_FEM_PROBLEM, T_USER_NON_LINEAR_PROBLEM, T_USER_NON_LINEAR_SOLUTION>
::solveTimeStep(const NumLib::TimeStep &t_n1)
{
	// print solving information
	INFO("--Solving nonlinear equation for the CMP_PressureForm process. \n");
	
	_nlin_solution->solveTimeStep(t_n1);

    // update the P and X vector
	_function_data->set_P(_nlin_solution->getCurrentSolution(0));
    _function_data->set_X(_nlin_solution->getCurrentSolution(1));
	
    // calcuate the equilibrium reaction system on each node
	_function_data->calc_nodal_eos_sys(t_n1.getTimeStepSize());



	return 0;
}


}

#endif  // end of ifndef
