/**
* Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.com/LICENSE.txt
*
*
* \file FunctionCMP_2P2C.hpp
*
* Created on 2014-05-12 by Yonghui HUANG
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
#include "EOS_GlobalComplementaryForm.h"

template <class T1, class T2>
bool FunctionCMP_GlobalComplementaryForm<T1, T2>::initialize(const BaseLib::Options &option)
{

	Ogs6FemData* femData = Ogs6FemData::getInstance();
	_msh_id = option.getOptionAsNum<size_t>("MeshID");
	size_t time_id = option.getOptionAsNum<int>("TimeGroupID");
	NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];

	//mesh and FE objects
	MeshLib::IMesh* msh = femData->list_mesh[_msh_id];
	// get the number of elements 
	size_t n_ele = msh->getNumberOfElements();
	size_t n_nodes = msh->getNumberOfNodes();
	MyDiscreteSystem* dis = 0;
	dis = DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(msh);
	_feObjects = new FemLib::LagrangeFeObjectContainer(msh);

	// set names of the output parameters
	this->setOutputParameterName(0, "LIQUID_PRESSURE");
	this->setOutputParameterName(1, "MOLAR_FRACTION");
	this->setOutputParameterName(2, "SATURATION");
	this->setOutputParameterName(3, "Gas_Pressure");
	this->setOutputParameterName(4, "Capillary_Pressure");
	this->setOutputParameterName(5, "Mass_Density_L_H");
	this->setOutputParameterName(6, "Mass_Density_G_H");
	// also secondary variables
	// TODO: set seconary variable names also as output parameters

	// create the MyCMPPressureForm problem
	_problem = new MyCMPGlobalComplementaryFormProblemType(dis);
	_problem->setTimeSteppingFunction(*tim);  // applying the time stepping function
	
	// creating mean pressure vector
	MyVariableCMPGlobalComplementaryForm* Liquid_Pressure = _problem->addVariable("LIQUID_PRESSURE");
	FemVariableBuilder var_builder;
	var_builder.doit("LIQUID_PRESSURE", option, msh, femData->geo, femData->geo_unique_name, _feObjects, Liquid_Pressure);
	SolutionLib::FemIC* femIC = _problem->getVariable(0)->getIC();
	_P = new MyNodalFunctionScalar();
	if (femIC)
	{
		// FemIC vector is not empty
		_P->initialize(*dis, _problem->getVariable(0)->getCurrentOrder(), 0.0);
		femIC->setup(*_P);
	}
	else
	{
		// FemIC vector is empty
		// initialize the vector with zeros
		_P->initialize(*dis, _problem->getVariable(0)->getCurrentOrder(), 0.0);
	}

	// creating molar fraction

	MyVariableCMPGlobalComplementaryForm* total_mass_density = _problem->addVariable("MOLAR_FRACTION");
	var_builder.doit("MOLAR_FRACTION", option, msh, femData->geo, femData->geo_unique_name, _feObjects, total_mass_density);
	femIC = _problem->getVariable(1)->getIC();
	_X = new MyNodalFunctionScalar();
	if (femIC)
	{
		// FemIC vector is not empty
		_X->initialize(*dis, _problem->getVariable(1)->getCurrentOrder(), 0.0);
		femIC->setup(*_X);
	}
	else
	{
		// FemIC vector is empty
		// initialize the vector with zeros
		_X->initialize(*dis, _problem->getVariable(1)->getCurrentOrder(), 0.0);
	}

	MyVariableCMPGlobalComplementaryForm* saturation = _problem->addVariable("SATURATION");
	var_builder.doit("SATURATION", option, msh, femData->geo, femData->geo_unique_name, _feObjects, saturation);
	femIC = _problem->getVariable(2)->getIC();
	_S = new MyNodalFunctionScalar();
	if (femIC)
	{
		// FemIC vector is not empty
		_S->initialize(*dis, _problem->getVariable(2)->getCurrentOrder(), 0.0);
		femIC->setup(*_S);
	}
	else
	{
		// FemIC vector is empty
		// initialize the vector with zeros
		_S->initialize(*dis, _problem->getVariable(2)->getCurrentOrder(), 0.0);
	}

	// initialize the local EOS problem
	
	//ogsChem::LocalVector _output1;
	//_output1 = ogsChem::LocalVector::Zero(_LP_EOS->N);

	_PC = new MyNodalFunctionScalar();	
	_PC->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
	//_PC->setValue(0, 0.0); //PC_ini
	_PL = new MyNodalFunctionScalar();
	_PL->initialize(*dis, FemLib::PolynomialOrder::Linear, 1e+6);//PL_ini
	_PG = new MyNodalFunctionScalar();
	_PG->initialize(*dis, FemLib::PolynomialOrder::Linear, 1e+6);
	_rho_L_h = new MyNodalFunctionScalar();//MASS density of hydrogen in liquid phase
	_rho_L_h->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
	_rho_G_h = new MyNodalFunctionScalar();//mass density of hydrogen in gas phase
	_rho_G_h->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
	
	_dPcdSg = new MyNodalFunctionScalar();
	_dPcdSg->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
	

	/**
	* initialize the vector of tempVar
	*/
	_mat_Jacob.resize(
	    3 * n_nodes,
	    3 * n_nodes);  // Resizing also initializes the matrix with zero

	_vec_Res = LocalVector::Zero(3 * n_nodes);

	//initialize the local M matrix and K matrix and partially Jacobian matrix
	for (size_t index = 0; index < n_ele; index++)
	{
		MathLib::LocalMatrix test1 = MathLib::LocalMatrix::Zero(9, 9);
		_elem_M_matrix.push_back(test1);
		MathLib::LocalMatrix test2 = MathLib::LocalMatrix::Zero(9, 9);
		_elem_K_matrix.push_back(test2);
		MathLib::LocalMatrix test3 = MathLib::LocalMatrix::Zero(9, 9);
		_elem_J_matrix.push_back(test3);
	}
	

	// linear assemblers
	MyNonLinearAssemblerType* nonlin_assembler = new MyNonLinearAssemblerType(_feObjects, this, msh->getGeometricProperty()->getCoordinateSystem());
	MyNonLinearResidualAssemblerType* non_linear_r_assembler = new MyNonLinearResidualAssemblerType(_feObjects, this, msh->getGeometricProperty()->getCoordinateSystem());
	MyNonLinearJacobianAssemblerType* non_linear_j_assembler = new MyNonLinearJacobianAssemblerType(_feObjects, this, msh->getGeometricProperty()->getCoordinateSystem());

	// define nonlinear problem
	_non_linear_problem = new MyNonLinearCMPGlobalComplementaryFORMProblemType(dis);
	_non_linear_eqs = _non_linear_problem->createEquation();
	_non_linear_eqs->initialize(nonlin_assembler, non_linear_r_assembler, non_linear_j_assembler);
	//_non_linear_eqs->initialize(nonlin_assembler,non_linear_r_assembler, non_linear_j_assembler);
	_non_linear_problem->setTimeSteppingFunction(*tim);
	_non_linear_problem->addVariable("LIQUID_PRESSURE");
	_non_linear_problem->addVariable("MOLAR_FRACTION");
	_non_linear_problem->addVariable("SATURATION");

	SolutionLib::FemIC* P_ic = new SolutionLib::FemIC(msh);
	P_ic->addDistribution(femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>(_P->getDiscreteData()));
	var_builder.doit("LIQUID_PRESSURE", option, msh, femData->geo, femData->geo_unique_name, _feObjects, _non_linear_problem->getVariable(0));
	_non_linear_problem->getVariable(0)->setIC(P_ic);

	SolutionLib::FemIC* X_ic = new SolutionLib::FemIC(msh);
	X_ic->addDistribution(femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>(_X->getDiscreteData()));
	var_builder.doit("MOLAR_FRACTION", option, msh, femData->geo, femData->geo_unique_name, _feObjects, _non_linear_problem->getVariable(1));
	_non_linear_problem->getVariable(1)->setIC(X_ic);

	SolutionLib::FemIC* S_ic = new SolutionLib::FemIC(msh);
	S_ic->addDistribution(femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>(_S->getDiscreteData()));
	var_builder.doit("SATURATION", option, msh, femData->geo, femData->geo_unique_name, _feObjects, _non_linear_problem->getVariable(2));
	_non_linear_problem->getVariable(2)->setIC(S_ic);

	// set up non-linear solution
	myNRIterator = new MyNRIterationStepInitializer(non_linear_r_assembler, non_linear_j_assembler);
	myNSolverFactory = new MyDiscreteNonlinearSolverFactory(myNRIterator);
	this->_non_linear_solution = new MyNonLinearSolutionType(dis, this->_non_linear_problem, this,myNSolverFactory);
	this->_non_linear_solution->getDofEquationIdTable()->setNumberingType(DiscreteLib::DofNumberingType::BY_POINT);  // global order
	this->_non_linear_solution->getDofEquationIdTable()->setLocalNumberingType(DiscreteLib::DofNumberingType::BY_VARIABLE);  // local order
	const BaseLib::Options* optNum = option.getSubGroup("Numerics");

	// linear solver
	MyLinearSolver* linear_solver = this->_non_linear_solution->getLinearEquationSolver();
	linear_solver->setOption(*optNum);
	// set nonlinear solver options
	this->_non_linear_solution->getNonlinearSolver()->setOption(*optNum);
    // set theta----------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (!optNum->hasOption("TimeTheta"))
    {
        ERR("Time theta setting not found!!!");
        exit(1);
    }
    else
    {
        double tmp_theta(optNum->getOptionAsNum<double>("TimeTheta"));
        non_linear_j_assembler->setTheta(tmp_theta);
        non_linear_r_assembler->setTheta(tmp_theta);
        _non_linear_eqs->getLinearAssembler()->setTheta(tmp_theta);
    }
	


	// get the nonlinear solution dof manager
	this->_nl_sol_dofManager = this->_non_linear_solution->getDofEquationIdTable();
	// set up solution
	_solution = new MyCMPGlobalComplementaryFormSolution(dis, _problem, _non_linear_solution, this);
	/**
	*Calculate the secondary variables on each node
	*/
	this->calc_nodal_eos_sys(0.0);
	

	// set initial output parameter
	OutputVariableInfo var1(this->getOutputParameterName(0), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _P);
	femData->outController.setOutput(var1.name, var1);
	OutputVariableInfo var2(this->getOutputParameterName(1), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _X);
	femData->outController.setOutput(var2.name, var2);
	OutputVariableInfo var3(this->getOutputParameterName(2), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _S);
	femData->outController.setOutput(var3.name, var3);

	OutputVariableInfo var_sec_1(this->getOutputParameterName(3), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _PG);
	femData->outController.setOutput(var_sec_1.name, var_sec_1);

	OutputVariableInfo var_sec_2(this->getOutputParameterName(4), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _PC);
	femData->outController.setOutput(var_sec_2.name, var_sec_2);

	OutputVariableInfo var_sec_3(this->getOutputParameterName(5), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _rho_L_h);
	femData->outController.setOutput(var_sec_3.name, var_sec_3);

	OutputVariableInfo var_sec_4(this->getOutputParameterName(6), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _rho_G_h);
	femData->outController.setOutput(var_sec_4.name, var_sec_4);

	return true;
}

template <class T1, class T2>
void FunctionCMP_GlobalComplementaryForm<T1, T2>::initializeTimeStep(const NumLib::TimeStep &/*time*/)
{
	/*
	size_t i;
	//const NumLib::ITXFunction *vel = this->getInput<NumLib::ITXFunction>(Velocity);

	// set velocity for linear problem
	for (i = 0; i < _linear_problems.size(); i++) {
		_linear_problems[i]->getEquation()->getLinearAssembler()->setVelocity(vel);
		_linear_problems[i]->getEquation()->getResidualAssembler()->setVelocity(vel);
		_linear_problems[i]->getEquation()->getJacobianAssembler()->setVelocity(vel);
	}
	*/
}

template <class T1, class T2>
void FunctionCMP_GlobalComplementaryForm<T1, T2>::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{

}

template <class T1, class T2>
void FunctionCMP_GlobalComplementaryForm<T1, T2>::output(const NumLib::TimeStep &/*time*/)
{
	
	//update data for output
	Ogs6FemData* femData = Ogs6FemData::getInstance();
	this->_msh_id = this->_problem->getDiscreteSystem()->getMesh()->getID();
	// set the new output
	// we have 2 primary variables
	OutputVariableInfo var1("LIQUID_PRESSURE", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _P);
	femData->outController.setOutput(var1.name, var1);
	OutputVariableInfo var2("MOLAR_FRACTION", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _X);
	femData->outController.setOutput(var2.name, var2);
	OutputVariableInfo var3("SATURATION", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _S);
	femData->outController.setOutput(var3.name, var3);
	// and 8 secondary variables
	// add all seconary variables as well. 
	OutputVariableInfo var_Sec_1("Gas_Pressure", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _PG);
	femData->outController.setOutput(var_Sec_1.name, var_Sec_1);

	OutputVariableInfo var_Sec_2("Capillary_Pressure", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _PC);
	femData->outController.setOutput(var_Sec_2.name, var_Sec_2);



}


template <class T1, class T2>
void FunctionCMP_GlobalComplementaryForm<T1, T2>::calc_nodal_eos_sys(double dt)
{
	// solve the EOS system on each node
	std::size_t node_id(0);
	std::size_t num_nodes(0);
	num_nodes = _P->getDiscreteData()->getRangeEnd();
	double EPS = 1E-7;
	MathLib::LocalMatrix matSecDer = MathLib::LocalMatrix::Zero(3, 2);
	MathLib::LocalVector vecTempVar = MathLib::LocalVector::Zero(2);
	MathLib::LocalMatrix test = MathLib::LocalMatrix::Zero(1, 1);

	ogsChem::LocalVector output = ogsChem::LocalVector::Zero(4);
	ogsChem::LocalVector output_ini = output;
	ogsChem::LocalVector output_dP = ogsChem::LocalVector::Zero(3);
	ogsChem::LocalVector output_dX = ogsChem::LocalVector::Zero(3);
	ogsChem::LocalVector input = ogsChem::LocalVector::Zero(3);
	
	INFO("Calculating EOS on each node...");

	// loop for each node
	for (node_id = _P->getDiscreteData()->getRangeBegin();
		node_id < _P->getDiscreteData()->getRangeEnd();
		node_id++)
	{
		input(0) = _P->getValue(node_id);
		input(1) = _X->getValue(node_id);
		input(2) = _S->getValue(node_id);

		
        double const x_L_h = input(1);
        double const pg = input(0);
        double const sg = input(2);
        double const N_L = 1000 / 0.018;
        double const H = 7.65e-6;
        double const var = pg*H / N_L - x_L_h;
        double const var_hold = pg*H / N_L;
        _vec_Res(2 * num_nodes + node_id) = -std::min(sg, pg*H/N_L- x_L_h);
 		if (input(2) <= (pg*H / N_L - x_L_h)){//then rho_L^h=C_h*PG
			//Calc each entry of the mass matrix
			_mat_Jacob.coeffRef(2 * num_nodes + node_id, 3 * node_id) = 0.0;
			_mat_Jacob.coeffRef(2 * num_nodes + node_id, 3 * node_id + 1) = 0.0;
			_mat_Jacob.coeffRef(2 * num_nodes + node_id, 3 * node_id+2) = 1.0;
		}
		else if (input(2) > (pg*H / N_L - x_L_h)){
			_mat_Jacob.coeffRef(node_id + 2 * num_nodes, 3 * node_id) = H / N_L;
			_mat_Jacob.coeffRef(node_id + 2 * num_nodes, 3 * node_id + 1) = -1;
			_mat_Jacob.coeffRef(node_id + 2 * num_nodes, 3 * node_id + 2) = 0.0;
		}
		// end
	}

		
	

}
