/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SingleStepReactOPS.h
 *
 * Created on 2013-03-26 by Haibing Shao
 */

#ifndef SINGLE_STEP_REACT_OPS_H
#define SINGLE_STEP_REACT_OPS_H

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

#include "ConcNodeInfo.h"

// class ReductionKinNodeInfo; 

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
    class T_USER_LINEAR_PROBLEM,
    class T_USER_LINEAR_SOLUTION,
    class T_USER_NON_LINEAR_PROBLEM, 
	class T_USER_NON_LINEAR_SOLUTION 
    >
class SingleStepReactOPS
    : public AbstractTimeSteppingAlgorithm
{
public:

	typedef T_USER_FUNCTION_DATA       UserFunctionData; 
    typedef T_USER_FEM_PROBLEM         UserFemProblem;
    typedef T_USER_LINEAR_PROBLEM      UserLinearProblem; 
	typedef T_USER_LINEAR_SOLUTION     UserLinearSolution;
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
    SingleStepReactOPS(MyDiscreteSystem*                    dis, 
		                   UserFemProblem*                  problem, 
						   UserFunctionData*                function_data, 
                           std::vector<UserLinearProblem*>  linear_problem,
						   std::vector<UserLinearSolution*> linear_solutions);

    /**
      * destructor, reclaiming the memory
      */ 
    virtual ~SingleStepReactOPS()
    {
        BaseLib::releaseObject(_problem); 
        
        BaseLib::releaseObjectsInStdVector( _vec_u_n1 ); 
        _discrete_system->deleteVector(_x_n0);
        _discrete_system->deleteVector(_x_n1);
        _discrete_system->deleteVector(_x_n1_0);
        _discrete_system->deleteVector(_x_st);
        BaseLib::releaseObject(_discrete_system);

        BaseLib::releaseObject(_feObjects);
        BaseLib::releaseObjectsInStdVector( _linear_problem );
        BaseLib::releaseObjectsInStdVector( _lin_solutions );
        BaseLib::releaseObject(_function_data);

        BaseLib::releaseObjectsInStdMap( _bc_info ); 
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
        // call all linear solutions
        for (size_t i=0; i < _lin_solutions.size(); i++ )
            if (!_lin_solutions[i]->accept(t)) return false;
        return true;
    }

    /**
      * called when this solution is accepted
      */
    virtual void finalizeTimeStep(const NumLib::TimeStep &t)
    {
        // this solution itself
        AbstractTimeSteppingAlgorithm::finalizeTimeStep(t);
        // call all linear solutions
        for (size_t i=0; i < _lin_solutions.size(); i++ )
            _lin_solutions[i]->finalizeTimeStep(t);
    }

    /** 
      * get the corresponding solution pointer
      */ 
    MyNodalFunctionScalar* getCurrentSolution(size_t var_id) { return _vec_u_n1[var_id]; }

    /**
      * get the pointer to the problem class
      */ 
    UserFemProblem* getProblem() {return _problem;};

    /** 
      * get the domain of freedom equation ID table
      */ 
    DiscreteLib::DofEquationIdTable* getDofEquationIdTable() {return &_dofManager;};

    /**
      * tell whether a particular node is a boundary node
      */ 
	bool isBCNode(size_t node_idx); 

private:
    DISALLOW_COPY_AND_ASSIGN(SingleStepReactOPS);

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
      * pointer to linear problems
      */ 
    std::vector<UserLinearProblem*>  _linear_problem;

    /**
      * pointer to linear solutions
      */ 
	std::vector<UserLinearSolution*> _lin_solutions; 

    /**
      * pointer to the function concentrations
      */ 
	UserFunctionData*                _function_data; 

    /**
      * boundary condition values
      */ 
	std::map<size_t, ConcNodeInfo*>  _bc_info; 
};

template <
	class T_USER_FUNCTION_DATA, 
    class T_USER_FEM_PROBLEM,
    class T_USER_LINEAR_PROBLEM,
    class T_USER_LINEAR_SOLUTION,
    class T_USER_NON_LINEAR_PROBLEM, 
	class T_USER_NON_LINEAR_SOLUTION 
    >
SingleStepReactOPS<T_USER_FUNCTION_DATA, T_USER_FEM_PROBLEM, T_USER_LINEAR_PROBLEM, T_USER_LINEAR_SOLUTION, T_USER_NON_LINEAR_PROBLEM, T_USER_NON_LINEAR_SOLUTION>
    ::SingleStepReactOPS(MyDiscreteSystem*                dis, 
                         UserFemProblem*                  problem, 
						 UserFunctionData*                function_data, 
                         std::vector<UserLinearProblem*>  linear_problem,
						 std::vector<UserLinearSolution*> linear_solutions)
    : AbstractTimeSteppingAlgorithm(*problem->getTimeSteppingFunction()),
      _problem(problem), _discrete_system(dis),
      _linear_problem(linear_problem), _lin_solutions(linear_solutions), _function_data(function_data)
{
    INFO("->Setting up a solution algorithm for Reactive Transport OPS problem.");

    size_t i, j, i_var; 
    const size_t n_var = problem->getNumberOfVariables();
    MeshLib::IMesh* msh = dis->getMesh();

    // create dof map
    for (size_t i=0; i<n_var; i++) {
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
    for (size_t i=0; i<n_var; i++) {
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
    for ( i=0; i<n_var; i++) {
        SolutionVector* vec_var = _vec_u_n1[i]->getDiscreteData();
        DiscreteLib::setGlobalVector(_dofManager, i, msh->getID(), *vec_var, *_x_n0);
    }

    // setup functions
    std::vector<MyVariable*> list_var(n_var);
    for ( i=0; i<n_var; i++) 
		list_var[i] = problem->getVariable(i);

	// getting the boundary conditions of concentrations for all components, 
    //const size_t msh_id = _discrete_system->getMesh()->getID();
    std::vector<size_t> list_bc1_eqs_id;
    std::vector<double> list_bc1_val;
    for ( i_var=0; i_var < n_var; i_var++ ) {
        MyVariable* var = _problem->getVariable(i_var);
        for ( i=0; i < var->getNumberOfDirichletBC(); i++ ) {
            SolutionLib::FemDirichletBC *bc1 = var->getDirichletBC(i);
            bc1->setup(var->getCurrentOrder());
            std::vector<size_t> &list_bc_nodes = bc1->getListOfBCNodes();
            std::vector<double> &list_bc_values = bc1->getListOfBCValues();
            
	        // now loop over this vector
		    for ( j=0; j < list_bc_nodes.size(); j++ )
		    {
				size_t node_id = list_bc_nodes[j]; 
				double node_value = list_bc_values[j]; 
                std::map<size_t, ConcNodeInfo*>::iterator my_bc; 
				my_bc = _bc_info.find( node_id ); 
                // now check whether this node has already been in the BC info
                if ( my_bc != _bc_info.end() )
                {   // this node is already there
					my_bc->second->set_comp_conc(i_var, node_value); 	             
                }
                else  // create a new structure and fill it in 
				{
					my_bc = _bc_info.begin(); 
					ConcNodeInfo* bc_node = new ConcNodeInfo( node_id, 
						                                      this->_function_data->get_n_Comp(),
                                                              this->_problem->getEqReactSys() ); 
					bc_node->set_comp_conc( i_var, node_value ); 
					_bc_info.insert( my_bc, std::pair<size_t, ConcNodeInfo*>(node_id, bc_node) ); 
				}  // end of if else

				// set BC concentration values
				_function_data->set_BC_conc_node_values(node_id, i_var, node_value);
	        }  // end of for j
        }  // end of for i
    }  // end of for i_var


	// loop over all the boundary nodes, and
	// transform these concentrations to eta and xi values
	std::map<size_t, ConcNodeInfo*>::iterator bc_node_it;
	std::vector<size_t> vec_bc_node_idx;
	std::vector<std::vector<double>> vec_node_conc_values;

	for (i = 0; i < _linear_problem.size(); i++)
	{
		std::vector<double> vec_conc;
		vec_node_conc_values.push_back(vec_conc);
	}

	for (bc_node_it = _bc_info.begin(); bc_node_it != _bc_info.end(); bc_node_it++)
	{
		size_t node_idx = bc_node_it->second->get_node_id();
		vec_bc_node_idx.push_back(node_idx);

		for (i = 0; i < _linear_problem.size(); i++)
		{
			double conc_value = bc_node_it->second->get_comp_conc(i);
			vec_node_conc_values[i].push_back(conc_value);
		}
	}

	// imposing BC for eta
	for (i = 0; i < _linear_problem.size(); i++)
	{
		_linear_problem[i]->getVariable(0)->addDirichletBC(new SolutionLib::FemDirichletBC(vec_bc_node_idx, vec_node_conc_values[i]));
	}

};

template <
	class T_USER_FUNCTION_DATA, 
    class T_USER_FEM_PROBLEM,
    class T_USER_LINEAR_PROBLEM,
    class T_USER_LINEAR_SOLUTION,
    class T_USER_NON_LINEAR_PROBLEM, 
	class T_USER_NON_LINEAR_SOLUTION 
    >
int SingleStepReactOPS<T_USER_FUNCTION_DATA, T_USER_FEM_PROBLEM, T_USER_LINEAR_PROBLEM, T_USER_LINEAR_SOLUTION, T_USER_NON_LINEAR_PROBLEM, T_USER_NON_LINEAR_SOLUTION>
    ::solveTimeStep(const NumLib::TimeStep &t_n1)
{
	size_t i; 

	// solving linear problems one after the other
	for ( i=0; i < _lin_solutions.size(); i++)
	{
        // read in the concentration value as initial condition
        _lin_solutions[i]->overwriteLastTimeStepValue( _function_data->get_concentrations(i) ); 

		// print solving information
        INFO("--Solving linear equation for the transport of Comp #%s. \n",  this->_linear_problem[i]->getVariable(0)->getName().c_str() ); 
		_lin_solutions[i]->solveTimeStep( t_n1 );
        
        _function_data->set_concentrations( i, _lin_solutions[i]->getCurrentSolution(0) ); 
	}
	// calcuate the equilibrium reaction system on each node
	_function_data->calc_nodal_react_sys(t_n1.getTimeStepSize());

    return 0;
}


template <
	class T_USER_FUNCTION_DATA, 
    class T_USER_FEM_PROBLEM,
    class T_USER_LINEAR_PROBLEM,
    class T_USER_LINEAR_SOLUTION,
    class T_USER_NON_LINEAR_PROBLEM, 
	class T_USER_NON_LINEAR_SOLUTION 
    >
bool SingleStepReactOPS<T_USER_FUNCTION_DATA, T_USER_FEM_PROBLEM, T_USER_LINEAR_PROBLEM, T_USER_LINEAR_SOLUTION, T_USER_NON_LINEAR_PROBLEM, T_USER_NON_LINEAR_SOLUTION>
    :: isBCNode(size_t node_idx)
{    
	if ( _bc_info.find( node_idx ) != _bc_info.end() )
		return true;  // is BC node
	else
		return false;  // not BC node
}

}

#endif  // end of ifndef
