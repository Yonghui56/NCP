/**
* Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.com/LICENSE.txt
*
*
* \file FunctionCMP_CapPresFormular.hpp
*
* Created on 2015-03-12 by Yonghui HUANG
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
#include "EOS_CapPresFormular.h"

template <class T1, class T2>
bool FunctionCMP_CapPresFormular<T1, T2>::initialize(const BaseLib::Options &option)
{

	Ogs6FemData* femData = Ogs6FemData::getInstance();
	_msh_id = option.getOptionAsNum<size_t>("MeshID");
	size_t time_id = option.getOptionAsNum<size_t>("TimeGroupID");
	NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];

	//mesh and FE objects
	MeshLib::IMesh* msh = femData->list_mesh[_msh_id];
	// get the number of elements 
	size_t n_ele = msh->getNumberOfElements();
	MyDiscreteSystem* dis = 0;
	dis = DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(msh);
	_feObjects = new FemLib::LagrangeFeObjectContainer(msh);

	// set names of the output parameters
	this->setOutputParameterName(0, "NONWET_PRESSURE");
	this->setOutputParameterName(1, "CAPILLARY_PRESSURE");
	this->setOutputParameterName(2, "Saturation");
	this->setOutputParameterName(3, "Liquid_Pressure");
	// also secondary variables
	// TODO: set seconary variable names also as output parameters

	// create the MyCMPCapPresFormular problem
	_problem = new MyCMPCapPresFormularProblemType(dis);
	_problem->setTimeSteppingFunction(*tim);  // applying the time stepping function
	
	// creating mean pressure vector
	MyVariableCMPCapPresFormular* mean_pressure = _problem->addVariable("NONWET_PRESSURE");
	FemVariableBuilder var_builder;
	var_builder.doit("NONWET_PRESSURE", option, msh, femData->geo, femData->geo_unique_name, _feObjects, mean_pressure);
	SolutionLib::FemIC* femIC = _problem->getVariable(0)->getIC();
	_PN = new MyNodalFunctionScalar();
	if (femIC)
	{
		// FemIC vector is not empty
		_PN->initialize(*dis, _problem->getVariable(0)->getCurrentOrder(), 0.0);
		femIC->setup(*_PN);
	}
	else
	{
		// FemIC vector is empty
		// initialize the vector with zeros
		_PN->initialize(*dis, _problem->getVariable(0)->getCurrentOrder(), 0.0);
	}

	// creating molar fraction
	MyVariableCMPCapPresFormular* molar_fraction = _problem->addVariable("CAPILLARY_PRESSURE");
	var_builder.doit("CAPILLARY_PRESSURE", option, msh, femData->geo, femData->geo_unique_name, _feObjects, molar_fraction);
	femIC = _problem->getVariable(1)->getIC();
	_PC = new MyNodalFunctionScalar();
	if (femIC)
	{
		// FemIC vector is not empty
		_PC->initialize(*dis, _problem->getVariable(1)->getCurrentOrder(), 0.0);
		femIC->setup(*_PC);
	}
	else
	{
		// FemIC vector is empty
		// initialize the vector with zeros
		_PC->initialize(*dis, _problem->getVariable(1)->getCurrentOrder(), 0.0);
	}

	// initialize the local EOS problem
	//_LP_EOS = new LocalProblem_EOS();
	//ogsChem::LocalVector _output1;
	//_output1 = ogsChem::LocalVector::Zero(_LP_EOS->N);

	/*_PC = new MyNodalFunctionScalar();	
	_PC->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
	//_PC->setValue(0, 0.0); //PC_ini
	_PL = new MyNodalFunctionScalar();
	_PL->initialize(*dis, FemLib::PolynomialOrder::Linear, 100.0);//PL_ini
	_PG = new MyNodalFunctionScalar();
	_PG->initialize(*dis, FemLib::PolynomialOrder::Linear, 100.0);
	_X_m = new MyNodalFunctionScalar();
	_X_m->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0381488377805587);
	_X_M = new MyNodalFunctionScalar();
	_X_M->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.190206132157118);
	_NG = new MyNodalFunctionScalar();
	_NG->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
	_NL = new MyNodalFunctionScalar();
	_NL->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
	_S = new MyNodalFunctionScalar();
	_S->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.00749928306032062);*/
	_S = new MyIntegrationPointFunctionVector();
	_S->initialize(dis);
	_PL = new MyIntegrationPointFunctionVector();
	_PL->initialize(dis);
	//_S->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0); 
	
	_dPc = new MyIntegrationPointFunctionVector();
	/*for (size_t idx = 0; idx < n_ele; idx++)
	{
		_S->setValue(idx, 0.0);
		_dPc->setValue(idx, 0.0);
	}
	*/
	_dPc->initialize(dis);
	//_dPc->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);

	_mat_secDer = new MyNodalFunctionMatrix();// define a matrix to store all the derivatives of the secondary variables 
	MathLib::LocalMatrix tmp_mat = MathLib::LocalMatrix::Zero(8, 2); //the size I will modify later
	_mat_secDer->initialize(*dis, FemLib::PolynomialOrder::Linear, tmp_mat);
	
	//initialize the local M matrix and K matrix
	for (size_t index = 0; index < n_ele; index++)
	{
		MathLib::LocalMatrix test1 = MathLib::LocalMatrix::Zero(6, 6);
		_elem_M_matrix.push_back(test1);
		MathLib::LocalMatrix test2 = MathLib::LocalMatrix::Zero(6, 6);
		_elem_K_matrix.push_back(test2);
	}
	
	// linear assemblers
	MyNonLinearAssemblerType* nonlin_assembler = new MyNonLinearAssemblerType(_feObjects, this, msh->getGeometricProperty()->getCoordinateSystem());
	MyNonLinearResidualAssemblerType* non_linear_r_assembler = new MyNonLinearResidualAssemblerType(_feObjects, this, msh->getGeometricProperty()->getCoordinateSystem());
	MyNonLinearJacobianAssemblerType* non_linear_j_assembler = new MyNonLinearJacobianAssemblerType(_feObjects, this,msh->getGeometricProperty()->getCoordinateSystem());

	// define nonlinear problem
	_non_linear_problem = new MyNonLinearCMPCapPresFormularProblemType(dis);
	_non_linear_eqs = _non_linear_problem->createEquation();
	_non_linear_eqs->initialize(nonlin_assembler, non_linear_r_assembler, non_linear_j_assembler);
	_non_linear_problem->setTimeSteppingFunction(*tim);
	_non_linear_problem->addVariable("NONWET_PRESSURE");
	_non_linear_problem->addVariable("CAPILLARY_PRESSURE");

	SolutionLib::FemIC* PN_ic = new SolutionLib::FemIC(msh);
	PN_ic->addDistribution(femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>(_PN->getDiscreteData()));
	var_builder.doit("NONWET_PRESSURE", option, msh, femData->geo, femData->geo_unique_name, _feObjects, _non_linear_problem->getVariable(0));
	_non_linear_problem->getVariable(0)->setIC(PN_ic);

	SolutionLib::FemIC* PC_ic = new SolutionLib::FemIC(msh);
	PC_ic->addDistribution(femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>(_PC->getDiscreteData()));
	var_builder.doit("CAPILLARY_PRESSURE", option, msh, femData->geo, femData->geo_unique_name, _feObjects, _non_linear_problem->getVariable(1));
	_non_linear_problem->getVariable(1)->setIC(PC_ic);

	// set up non-linear solution
	myNRIterator = new MyNRIterationStepInitializer(non_linear_r_assembler, non_linear_j_assembler);
	myNSolverFactory = new MyDiscreteNonlinearSolverFactory(myNRIterator);
	this->_non_linear_solution = new MyNonLinearSolutionType(dis, this->_non_linear_problem, myNSolverFactory);
	this->_non_linear_solution->getDofEquationIdTable()->setNumberingType(DiscreteLib::DofNumberingType::BY_POINT);  // global order
	this->_non_linear_solution->getDofEquationIdTable()->setLocalNumberingType(DiscreteLib::DofNumberingType::BY_VARIABLE);  // local order
	const BaseLib::Options* optNum = option.getSubGroup("Numerics");

	// linear solver
	MyLinearSolver* linear_solver = this->_non_linear_solution->getLinearEquationSolver();
	linear_solver->setOption(*optNum);
	// set nonlinear solver options
	this->_non_linear_solution->getNonlinearSolver()->setOption(*optNum);
    // set theta
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
	_solution = new MyCMPCapPresFormularSolution(dis, _problem, _non_linear_solution, this);
	/**
	*Calculate the secondary variables on each node
	*/
	this->calc_nodal_eos_sys(0.0);
	

	// set initial output parameter
	OutputVariableInfo var1(this->getOutputParameterName(0), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _PN);
	femData->outController.setOutput(var1.name, var1);
	OutputVariableInfo var2(this->getOutputParameterName(1), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _PC);
	femData->outController.setOutput(var2.name, var2);
	OutputVariableInfo var_Sec_1(this->getOutputParameterName(2), _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _S);
	femData->outController.setOutput(var_Sec_1.name, var_Sec_1);
	
	return true;
}

template <class T1, class T2>
void FunctionCMP_CapPresFormular<T1, T2>::initializeTimeStep(const NumLib::TimeStep &/*time*/)
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
void FunctionCMP_CapPresFormular<T1, T2>::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{

}

template <class T1, class T2>
void FunctionCMP_CapPresFormular<T1, T2>::output(const NumLib::TimeStep &/*time*/)
{
	
	//update data for output
	Ogs6FemData* femData = Ogs6FemData::getInstance();
	this->_msh_id = this->_problem->getDiscreteSystem()->getMesh()->getID();
	// set the new output
	// we have 2 primary variables
	OutputVariableInfo var1("NONWET_PRESSURE", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _PN);
	femData->outController.setOutput(var1.name, var1);
	OutputVariableInfo var2("CAPILLARY_PRESSURE", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _PC);
	femData->outController.setOutput(var2.name, var2);
	// and 8 secondary variables
	// add all seconary variables as well. 
	//OutputVariableInfo var_Sec_1("Saturation", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _S);
	//femData->outController.setOutput(var_Sec_1.name, var_Sec_1);
	OutputVariableInfo var_Sec_1("Saturation", _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _S);
	femData->outController.setOutput(var_Sec_1.name, var_Sec_1);
	//OutputVariableInfo var_Sec_2("Derived Capillary Pressure", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _dPc);
	OutputVariableInfo var_Sec_3("Liquid_Pressure", _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _PL);
	femData->outController.setOutput(var_Sec_3.name, var_Sec_3);
	
	//femData->outController.setOutput(var_Sec_2.name, var_Sec_2);
	OutputVariableInfo var_Sec_2("Derived Capillary Pressure", _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _dPc);
	femData->outController.setOutput(var_Sec_2.name, var_Sec_2);
	
}


template <class T1, class T2>
void FunctionCMP_CapPresFormular<T1, T2>::calc_nodal_eos_sys(double dt = 0.0)
{
	/*std::size_t node_id(0);
	//MathLib::LocalMatrix matSecDer = MathLib::LocalMatrix::Zero(8, 2);
	MathLib::LocalVector vecSecDer = MathLib::LocalVector::Zero(1);

	ogsChem::LocalVector output = ogsChem::LocalVector::Zero(2);    
	
	ogsChem::LocalVector input = ogsChem::LocalVector::Zero(2);

	INFO("Calculating EOS on each node...");

	for (node_id = _P->getDiscreteData()->getRangeBegin();
		node_id < _P->getDiscreteData()->getRangeEnd();
		node_id++)
	{
		input(0) = _P->getValue(node_id);
		input(1) = _PC->getValue(node_id);

		//output(0) = _S->getValue(node_id);

		_EOS->set_env_condition(input);//apply the primary variables as the input variables for Eos
		output(0) = _EOS->getSat_byPC();
		output(1) = _EOS->get_diriv_PC();
		_S->setValue(node_id, output(0));
		_dPc->setValue(node_id, output(1));

		//_LP_EOS->solve(input, output);
	}
	*/
	

}


