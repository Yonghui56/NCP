/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Concentration.h
 *
 * Created on 2012-09-06 by Haibing Shao
 */

#include "logog.hpp"

#include "MathLib/DataType.h"
#include "DiscreteLib/Utils/DiscreteSystemContainerPerMesh.h"
#include "FemLib/Function/FemNodalFunction.h"
#include "NumLib/Function/TXFunctionBuilder.h"
#include "NumLib/Function/TXFunctionDirect.h"
#include "OutputIO/OutputBuilder.h"
#include "OutputIO/OutputTimingBuilder.h"
#include "SolutionLib/Fem/FemSourceTerm.h"
#include "Ogs6FemData.h"
#include "MathLib/ODE/RungeKutta4.h"

template <class T1, class T2>
bool FunctionOPSConc<T1,T2>::initialize(const BaseLib::Options &option)
{
    size_t i;  // index
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    _msh_id = option.getOptionAsNum<size_t>("MeshID");
    size_t time_id = option.getOptionAsNum<size_t>("TimeGroupID");
    NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];
        
    //mesh and FE objects
    MeshLib::IMesh* msh = femData->list_mesh[_msh_id];
    MyDiscreteSystem* dis = 0;
    dis = DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(msh);
    _feObjects = new FemLib::LagrangeFeObjectContainer(msh);

	// initialize the equilibrium reaction system
	_local_eq_react_sys = femData->m_EqReactSys;
	_local_kin_react_sys = femData->m_KinReactSys;

    // make sure the reduction scheme is already initialized. 
    if (!(this->_local_eq_react_sys->IsInitialized()) && !(this->_local_kin_react_sys->IsInitialized()))
	{
		// error msg
	    ERR("While initialize the OPS Reactive Transport Process, the chemEqReactSysActivity class has not been correctly initialized! ");
		// then stop the program
		exit(1);
	}

    // getting the number of components
	if (this->_local_eq_react_sys->get_n_Eq_React() > 0 && this->_local_eq_react_sys->get_n_Eq_React() > 0)
	{
		_n_Comp = std::max(this->_local_eq_react_sys->get_n_Comp(), this->_local_kin_react_sys->get_n_Comp());
		_n_Comp_mob = std::max(this->_local_eq_react_sys->get_n_Comp_mob(), this->_local_kin_react_sys->get_n_Comp_mob());
	}
	else if (this->_local_eq_react_sys->get_n_Eq_React() > 0)
	{
		_n_Comp = this->_local_eq_react_sys->get_n_Comp();
		_n_Comp_mob = this->_local_eq_react_sys->get_n_Comp_mob();
	}
	else if (this->_local_kin_react_sys->get_n_Kin_React() > 0)
	{
		_n_Comp = this->_local_kin_react_sys->get_n_Comp();
		_n_Comp_mob = this->_local_kin_react_sys->get_n_Comp_mob(); 
	}

	// set concentrations of all components as output
	for ( i=0; i < _n_Comp; i++ )
		this->setOutputParameterName( i, femData->map_ChemComp[i]->second->get_name() ); 

	// reactive transport problem with OPS
	_problem = new MyReactOPSProblemType( dis, _local_eq_react_sys ); 
    _problem->setTimeSteppingFunction(*tim);  // applying the same time stepping function for all linear problems

    // creating concentrations vector
	MyNodalFunctionScalar* tmp_conc;
	// first basis component
    for ( i=0; i < _n_Comp; i++ )
	{
		if (femData->map_ChemComp[i]->second->getCompType() == ogsChem::BASIS_COMP)
		{
			MyVariableConc* comp = _problem->addVariable(femData->map_ChemComp[i]->second->get_name());
			FemVariableBuilder var_builder;
			var_builder.doit(femData->map_ChemComp[i]->second->get_name(), option, msh, femData->geo, femData->geo_unique_name, _feObjects, comp);
			SolutionLib::FemIC* femIC = _problem->getVariable(i)->getIC();
			tmp_conc = new MyNodalFunctionScalar();
			if (femIC)
			{
				// FemIC vector is not empty
				tmp_conc->initialize(*dis, _problem->getVariable(i)->getCurrentOrder(), 1.0e-30);
				femIC->setup(*tmp_conc);
			}
			else
			{
				// FemIC vector is empty
				// initialize the vector with zeros
				tmp_conc->initialize(*dis, _problem->getVariable(i)->getCurrentOrder(), 1.0e-30);
			}
			_concentrations.push_back(tmp_conc);
		}
    }
	// then secondary mobile
	for (i = 0; i < _n_Comp; i++)
	{
		if (femData->map_ChemComp[i]->second->getCompType() == ogsChem::AQ_PHASE_COMP)
		{
			MyVariableConc* comp = _problem->addVariable(femData->map_ChemComp[i]->second->get_name());
			FemVariableBuilder var_builder;
			var_builder.doit(femData->map_ChemComp[i]->second->get_name(), option, msh, femData->geo, femData->geo_unique_name, _feObjects, comp);
			SolutionLib::FemIC* femIC = _problem->getVariable(i)->getIC();
			tmp_conc = new MyNodalFunctionScalar();
			if (femIC)
			{
				// FemIC vector is not empty
				tmp_conc->initialize(*dis, _problem->getVariable(i)->getCurrentOrder(), 1.0e-30);
				femIC->setup(*tmp_conc);
			}
			else
			{
				// FemIC vector is empty
				// initialize the vector with zeros
				tmp_conc->initialize(*dis, _problem->getVariable(i)->getCurrentOrder(), 1.0e-30);
			}
			_concentrations.push_back(tmp_conc);
		}
	}
	// then secondary sorption
	for (i = 0; i < _n_Comp; i++)
	{
		if (femData->map_ChemComp[i]->second->getCompType() == ogsChem::SORPTION_COMP)
		{
			MyVariableConc* comp = _problem->addVariable(femData->map_ChemComp[i]->second->get_name());
			FemVariableBuilder var_builder;
			var_builder.doit(femData->map_ChemComp[i]->second->get_name(), option, msh, femData->geo, femData->geo_unique_name, _feObjects, comp);
			SolutionLib::FemIC* femIC = _problem->getVariable(i)->getIC();
			tmp_conc = new MyNodalFunctionScalar();
			if (femIC)
			{
				// FemIC vector is not empty
				tmp_conc->initialize(*dis, _problem->getVariable(i)->getCurrentOrder(), 1.0e-30);
				femIC->setup(*tmp_conc);
			}
			else
			{
				// FemIC vector is empty
				// initialize the vector with zeros
				tmp_conc->initialize(*dis, _problem->getVariable(i)->getCurrentOrder(), 1.0e-30);
			}
			_concentrations.push_back(tmp_conc);
		}
	}
	// then secondary mineral
	for (i = 0; i < _n_Comp; i++)
	{
		if (femData->map_ChemComp[i]->second->getCompType() == ogsChem::MIN_PHASE_COMP)
		{
			MyVariableConc* comp = _problem->addVariable(femData->map_ChemComp[i]->second->get_name());
			FemVariableBuilder var_builder;
			var_builder.doit(femData->map_ChemComp[i]->second->get_name(), option, msh, femData->geo, femData->geo_unique_name, _feObjects, comp);
			SolutionLib::FemIC* femIC = _problem->getVariable(i)->getIC();
			tmp_conc = new MyNodalFunctionScalar();
			if (femIC)
			{
				// FemIC vector is not empty
				tmp_conc->initialize(*dis, _problem->getVariable(i)->getCurrentOrder(), 1.0e-30);
				femIC->setup(*tmp_conc);
			}
			else
			{
				// FemIC vector is empty
				// initialize the vector with zeros
				tmp_conc->initialize(*dis, _problem->getVariable(i)->getCurrentOrder(), 1.0e-30);
			}
			_concentrations.push_back(tmp_conc);
		}
	}
	// then kinetic immob
	for (i = 0; i < _n_Comp; i++)
	{
		if (femData->map_ChemComp[i]->second->getCompType() == ogsChem::KIN_COMP)
		{
			MyVariableConc* comp = _problem->addVariable(femData->map_ChemComp[i]->second->get_name());
			FemVariableBuilder var_builder;
			var_builder.doit(femData->map_ChemComp[i]->second->get_name(), option, msh, femData->geo, femData->geo_unique_name, _feObjects, comp);
			SolutionLib::FemIC* femIC = _problem->getVariable(i)->getIC();
			tmp_conc = new MyNodalFunctionScalar();
			if (femIC)
			{
				// FemIC vector is not empty
				tmp_conc->initialize(*dis, _problem->getVariable(i)->getCurrentOrder(), 1.0e-30);
				femIC->setup(*tmp_conc);
			}
			else
			{
				// FemIC vector is empty
				// initialize the vector with zeros
				tmp_conc->initialize(*dis, _problem->getVariable(i)->getCurrentOrder(), 1.0e-30);
			}
			_concentrations.push_back(tmp_conc);
		}
	}

	// adding variables into the linear problem
	// for the linear transport problem, variables are c_mobha
	for (i = 0; i < _n_Comp_mob; i++)
	{
		// set up problem
		MyLinearTransportProblemType*  linear_problem = new MyLinearTransportProblemType(dis);

		// linear assemblers
		MyLinearAssemblerType* linear_assembler = new MyLinearAssemblerType(_feObjects);
		MyLinearResidualAssemblerType* linear_r_assembler = new MyLinearResidualAssemblerType(_feObjects);
		MyLinearJacobianAssemblerType* linear_j_eqs = new MyLinearJacobianAssemblerType(_feObjects);

		MyLinearEquationType* linear_eqs = linear_problem->createEquation();
		linear_eqs->initialize(linear_assembler, linear_r_assembler, linear_j_eqs);
		linear_problem->setTimeSteppingFunction(*tim);

		// set up variables
		// in this case, the variables are concentrations of all mobile chemical components 
		MyVariableConc* comp_conc = linear_problem->addVariable(femData->map_ChemComp[i]->second->get_name());
		_linear_problems.push_back(linear_problem);
	}

	// set IC for linear transport components
	for (i = 0; i < _n_Comp_mob; i++)
	{
		SolutionLib::FemIC* linear_femIC = new SolutionLib::FemIC(msh);
		linear_femIC->addDistribution(femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>(_concentrations[i]->getDiscreteData()));
		_linear_problems[i]->getVariable(0)->setIC(linear_femIC);  
	}

	// set up linear solutions
	for (i = 0; i < _n_Comp_mob; i++)
	{
		MyLinearSolutionType* linear_solution = new MyLinearSolutionType(dis, this->_linear_problems[i]);
		MyLinearSolver* linear_solver = linear_solution->getLinearEquationSolver();
		const BaseLib::Options* optNum = option.getSubGroup("Numerics");
		linear_solver->setOption(*optNum);
		this->_linear_solutions.push_back(linear_solution);
	}

	// set up solution
	_solution = new MyReactOPSSolution(dis, _problem, this, _linear_problems, _linear_solutions);

	// run equilibrium reactions on each node that is not on boundary
	this->calc_nodal_react_sys(0.0);

	// set initial output parameter
	for (i = 0; i<_concentrations.size(); i++) {
		OutputVariableInfo var1(this->getOutputParameterName(i), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _concentrations[i]);
		femData->outController.setOutput(var1.name, var1);
	}

    return true;
}

template <class T1, class T2>
void FunctionOPSConc<T1, T2>::initializeTimeStep(const NumLib::TimeStep &/*time*/)
{
	size_t i; 
    const NumLib::ITXFunction *vel = this->getInput<NumLib::ITXFunction>(Velocity);

	// set velocity for linear problem
	for ( i=0; i < _linear_problems.size(); i++ ) {
		_linear_problems[i]->getEquation()->getLinearAssembler()->setVelocity(vel);
		_linear_problems[i]->getEquation()->getResidualAssembler()->setVelocity(vel);
		_linear_problems[i]->getEquation()->getJacobianAssembler()->setVelocity(vel);
	}
}

template <class T1, class T2>
void FunctionOPSConc<T1, T2>::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{ 

}

template <class T1, class T2>
void FunctionOPSConc<T1, T2>::output(const NumLib::TimeStep &/*time*/)
{
    size_t i;

    // update data for output
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    this->_msh_id = this->_problem->getDiscreteSystem()->getMesh()->getID(); 
	// set the new output
	for (i=0; i<_concentrations.size(); i++) {
		OutputVariableInfo var1(this->getOutputParameterName(i), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _concentrations[i]);
        femData->outController.setOutput(var1.name, var1);
    }

}

template <class T1, class T2>
void FunctionOPSConc<T1, T2>::set_BC_conc_node_values(std::size_t node_idx, std::size_t i_var, double node_value)
{
	_concentrations[i_var]->setValue(node_idx, node_value);
}

template <class T1, class T2>
void FunctionOPSConc<T1, T2>::calc_nodal_react_sys(double dt)
{
    size_t i, node_idx, err_node_count; 
    size_t result = 1; 
	
	// initialize the local vector
	LocalVector loc_conc; 
	loc_conc = LocalVector::Zero( _n_Comp ); 

	// signal, solving local equilibrium reaction system
	INFO("--Solving local chemical reaction system...");

    // clear the counter
    err_node_count = 0; 
    // loop over all the nodes
    for (node_idx = _concentrations[0]->getDiscreteData()->getRangeBegin();
	     node_idx < _concentrations[0]->getDiscreteData()->getRangeEnd(); 
		 node_idx++ )
	{
		// skip the boundary nodes
		if ( ! this->_solution->isBCNode(node_idx) )
		{
			// on each node, get the right start value
			// get the concentration vector
			for (i = 0; i < _n_Comp; i++)
			{
				double tmp_conc; 
				tmp_conc = this->_concentrations[i]->getValue(node_idx);
				if (tmp_conc < 0.0)
					tmp_conc = 1.0e-30; 
				loc_conc[i] = tmp_conc; 
			}

			// solve the local equilibrium system
			if (this->_local_eq_react_sys->get_n_Eq_React() > 0 )
				this->_local_eq_react_sys->solve_EqSys_Newton( loc_conc, 
                                                               result, 
                                                               node_idx, 
		                                                       1.0e-8, 
                                                               1.0e-12,
                                                               50 ); 

			// solving kinetic reactions. 
			if (this->_local_kin_react_sys->get_n_Kin_React() > 0)
				this->_local_kin_react_sys->solve_KinSys(loc_conc,
				                                         result,
				                                         node_idx,
				                                         1.0e-12,
				                                         1.0e-12,
				                                         dt,
				                                         50); 
        	// if the iteration converged. 
			if (result == 0)
			{
				// collect the solved concentrations
				for (i = 0; i < _n_Comp; i++)
					this->_concentrations[i]->setValue(node_idx, loc_conc[i]);
			}
            else 
                err_node_count++; 
		} // end of if
        else 
        {
            // if it is boundary nodes
            // then do nothing, using the originally set concentrations
        }

	}  // end of for node_idx

    if ( err_node_count == 0 )
    {
        INFO("--Solution of chemical system all converged...");
    }
    else
    {
        INFO("--Solution of chemical system failed on %i nodes...", err_node_count);
    }
}

