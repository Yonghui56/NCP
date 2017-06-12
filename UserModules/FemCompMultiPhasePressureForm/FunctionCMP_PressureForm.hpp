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
#include "EOS_PressureForm.h"

template <class T1, class T2>
bool FunctionCMP_PressureForm<T1, T2>::initialize(const BaseLib::Options &option)
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
	this->setOutputParameterName(0, "LIQUID_PRESSURE");
	this->setOutputParameterName(1, "MOLAR_FRACTION");
	this->setOutputParameterName(2, "Saturation");
	this->setOutputParameterName(3, "Gas_Pressure");
	this->setOutputParameterName(4, "Capillary_Pressure");
	// also secondary variables
	// TODO: set seconary variable names also as output parameters

	// create the MyCMPPressureForm problem
	_problem = new MyCMPPressureFormProblemType(dis);
	_problem->setTimeSteppingFunction(*tim);  // applying the time stepping function
	
	// creating mean pressure vector
	MyVariableCMPPressureForm* mean_pressure = _problem->addVariable("LIQUID_PRESSURE");
	FemVariableBuilder var_builder;
	var_builder.doit("LIQUID_PRESSURE", option, msh, femData->geo, femData->geo_unique_name, _feObjects, mean_pressure);
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
	MyVariableCMPPressureForm* molar_fraction = _problem->addVariable("MOLAR_FRACTION");
	var_builder.doit("MOLAR_FRACTION", option, msh, femData->geo, femData->geo_unique_name, _feObjects, molar_fraction);
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

	// initialize the local EOS problem
	//_LP_EOS = new LocalProblem_EOS();
	//ogsChem::LocalVector _output1;
	//_output1 = ogsChem::LocalVector::Zero(_LP_EOS->N);
	/*
	* Initialize of the secondary variables on each element
	*/
	_S = new MyIntegrationPointFunctionVector();
	_S->initialize(dis);
	_PG = new MyIntegrationPointFunctionVector();
	_PG->initialize(dis);
	_PC = new MyIntegrationPointFunctionVector();
	_PC->initialize(dis);
	_dPc = new MyIntegrationPointFunctionVector();
	_dPc->initialize(dis);

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
	MyNonLinearJacobianAssemblerType* non_linear_j_assembler = new MyNonLinearJacobianAssemblerType(_feObjects, this, msh->getGeometricProperty()->getCoordinateSystem());

	// define nonlinear problem
	_non_linear_problem = new MyNonLinearCMPPRESSUREFORMProblemType(dis);
	_non_linear_eqs = _non_linear_problem->createEquation();
	_non_linear_eqs->initialize(nonlin_assembler, non_linear_r_assembler, non_linear_j_assembler);
	_non_linear_problem->setTimeSteppingFunction(*tim);
	_non_linear_problem->addVariable("LIQUID_PRESSURE");
	_non_linear_problem->addVariable("MOLAR_FRACTION");

	SolutionLib::FemIC* P_ic = new SolutionLib::FemIC(msh);
	P_ic->addDistribution(femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>(_P->getDiscreteData()));
	var_builder.doit("LIQUID_PRESSURE", option, msh, femData->geo, femData->geo_unique_name, _feObjects, _non_linear_problem->getVariable(0));
	_non_linear_problem->getVariable(0)->setIC(P_ic);

	SolutionLib::FemIC* X_ic = new SolutionLib::FemIC(msh);
	X_ic->addDistribution(femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>(_X->getDiscreteData()));
	var_builder.doit("MOLAR_FRACTION", option, msh, femData->geo, femData->geo_unique_name, _feObjects, _non_linear_problem->getVariable(1));
	_non_linear_problem->getVariable(1)->setIC(X_ic);

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
	_solution = new MyCMPPressureFormSolution(dis, _problem, _non_linear_solution, this);
	/**
	*Calculate the secondary variables on each node
	*/
	this->calc_nodal_eos_sys(0.0);
	

	// set initial output parameter
	OutputVariableInfo var1(this->getOutputParameterName(0), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _P);
	femData->outController.setOutput(var1.name, var1);
	OutputVariableInfo var2(this->getOutputParameterName(1), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _X);
	femData->outController.setOutput(var2.name, var2);
	// set initial output parameter for secondary variable 
	// these are on each element
	OutputVariableInfo var_Sec_1(this->getOutputParameterName(2), _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _S);
	femData->outController.setOutput(var_Sec_1.name, var_Sec_1);

	OutputVariableInfo var_Sec_2(this->getOutputParameterName(3), _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _PG);
	femData->outController.setOutput(var_Sec_2.name, var_Sec_2);

	OutputVariableInfo var_Sec_3(this->getOutputParameterName(4), _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _PC);
	femData->outController.setOutput(var_Sec_3.name, var_Sec_3);

	return true;
}

template <class T1, class T2>
void FunctionCMP_PressureForm<T1, T2>::initializeTimeStep(const NumLib::TimeStep &/*time*/)
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
void FunctionCMP_PressureForm<T1, T2>::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{

}

template <class T1, class T2>
void FunctionCMP_PressureForm<T1, T2>::output(const NumLib::TimeStep &/*time*/)
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
	
	//we have several secondary variables the values of which are defined on each elements
	OutputVariableInfo var_Sec_1("Saturation", _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _S);
	femData->outController.setOutput(var_Sec_1.name, var_Sec_1);

	OutputVariableInfo var_Sec_2("Gas_Pressure", _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _PG);
	femData->outController.setOutput(var_Sec_2.name, var_Sec_2);

	OutputVariableInfo var_Sec_3("Capillary_Pressure", _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _PC);
	femData->outController.setOutput(var_Sec_3.name, var_Sec_3);

	
}


template <class T1, class T2>
void FunctionCMP_PressureForm<T1, T2>::calc_nodal_eos_sys(double dt = 0.0)
{

}


