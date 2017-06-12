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
#include "EOS_TotalDensityForm.h"

template <class T1, class T2>
bool FunctionCMP_TotalDensityForm<T1, T2>::initialize(const BaseLib::Options &option)
{

	Ogs6FemData* femData = Ogs6FemData::getInstance();
	_msh_id = option.getOptionAsNum<size_t>("MeshID");
	size_t time_id = option.getOptionAsNum<size_t>("TimeGroupID");
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
	this->setOutputParameterName(0, "MEAN_PRESSURE");
	this->setOutputParameterName(1, "TOTAL_MASS_DENSITY");

	this->setOutputParameterName(2, "Saturation");
	this->setOutputParameterName(3, "Liquid_Pressure");
	this->setOutputParameterName(4, "Gas_Pressure");
	this->setOutputParameterName(5, "Capillary_Pressure");
	this->setOutputParameterName(6, "Mass_Density_L_H");
	this->setOutputParameterName(7, "Mass_Density_G_H");
	this->setOutputParameterName(8, "Mass_Density_G_W");
	// also secondary variables
	// TODO: set seconary variable names also as output parameters

	// create the MyCMPPressureForm problem
	_problem = new MyCMPTotalDensityFormProblemType(dis);
	_problem->setTimeSteppingFunction(*tim);  // applying the time stepping function
	
	// creating mean pressure vector
	MyVariableCMPTotalDensityForm* mean_pressure = _problem->addVariable("MEAN_PRESSURE");
	FemVariableBuilder var_builder;
	var_builder.doit("MEAN_PRESSURE", option, msh, femData->geo, femData->geo_unique_name, _feObjects, mean_pressure);
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

	MyVariableCMPTotalDensityForm* total_mass_density = _problem->addVariable("TOTAL_MASS_DENSITY");
	var_builder.doit("TOTAL_MASS_DENSITY", option, msh, femData->geo, femData->geo_unique_name, _feObjects, total_mass_density);
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
	_LP_EOS = new LocalProblem_EOS_TotalDensity();
	//ogsChem::LocalVector _output1;
	//_output1 = ogsChem::LocalVector::Zero(_LP_EOS->N);
	_S = new MyIntegrationPointFunctionVector();
	_S->initialize(dis);
	_PL = new MyIntegrationPointFunctionVector();
	_PL->initialize(dis);
	_PG = new  MyIntegrationPointFunctionVector();
	_PG->initialize(dis);
	_PC = new MyIntegrationPointFunctionVector();
	_PC->initialize(dis);
	_rho_L_h = new MyIntegrationPointFunctionVector();
	_rho_L_h->initialize(dis);
	_rho_G_h = new MyIntegrationPointFunctionVector();
	_rho_G_h->initialize(dis);
	_rho_G_w = new  MyIntegrationPointFunctionVector();
	_rho_G_w->initialize(dis);
	MathLib::LocalVector tmp = MathLib::LocalVector::Zero(1);
	for (size_t ele_id = 0; ele_id < n_ele; ele_id++)
	{
		_S->setNumberOfIntegationPoints(ele_id, 4);	
		_PL->setNumberOfIntegationPoints(ele_id, 4);
		_PG->setNumberOfIntegationPoints(ele_id, 4);
		_PC->setNumberOfIntegationPoints(ele_id, 4);
		_rho_G_h->setNumberOfIntegationPoints(ele_id, 4);
		_rho_G_w->setNumberOfIntegationPoints(ele_id, 4);
		_rho_L_h->setNumberOfIntegationPoints(ele_id, 4);
		for (size_t jj = 0; jj < 4; jj++)
		{
			_S->setIntegrationPointValue(ele_id, jj, tmp);
			_PL->setIntegrationPointValue(ele_id, jj, tmp);
			_PG -> setIntegrationPointValue(ele_id, jj, tmp);
			_PC->setIntegrationPointValue(ele_id, jj, tmp);
			_rho_G_h->setIntegrationPointValue(ele_id, jj, tmp);
			_rho_L_h->setIntegrationPointValue(ele_id, jj, tmp);
			_rho_G_w->setIntegrationPointValue(ele_id, jj, tmp);
		}
	}
	
	_mat_secDer = new MyNodalFunctionMatrix();// define a matrix to store all the derivatives of the secondary variables 
	MathLib::LocalMatrix tmp_mat = MathLib::LocalMatrix::Zero(3, 2); //the size I will modify later
	_mat_secDer->initialize(*dis, FemLib::PolynomialOrder::Linear, tmp_mat);
	/**
	* initialize the vector of tempVar
	*/
	
	_vec_tempVar = new MyNodalFunctionVector();
	MathLib::LocalVector tmp_vec = MathLib::LocalVector::Zero(2);
	_vec_tempVar->initialize(*dis, FemLib::PolynomialOrder::Linear, tmp_vec);


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
	_non_linear_problem = new MyNonLinearCMPTOTALDENSITYFORMProblemType(dis);
	_non_linear_eqs = _non_linear_problem->createEquation();
	_non_linear_eqs->initialize(nonlin_assembler, non_linear_r_assembler, non_linear_j_assembler);
	_non_linear_problem->setTimeSteppingFunction(*tim);
	_non_linear_problem->addVariable("MEAN_PRESSURE");
	_non_linear_problem->addVariable("TOTAL_MASS_DENSITY");

	SolutionLib::FemIC* P_ic = new SolutionLib::FemIC(msh);
	P_ic->addDistribution(femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>(_P->getDiscreteData()));
	var_builder.doit("MEAN_PRESSURE", option, msh, femData->geo, femData->geo_unique_name, _feObjects, _non_linear_problem->getVariable(0));
	_non_linear_problem->getVariable(0)->setIC(P_ic);

	SolutionLib::FemIC* X_ic = new SolutionLib::FemIC(msh);
	X_ic->addDistribution(femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>(_X->getDiscreteData()));
	var_builder.doit("TOTAL_MASS_DENSITY", option, msh, femData->geo, femData->geo_unique_name, _feObjects, _non_linear_problem->getVariable(1));
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
	_solution = new MyCMPTotalDensityFormSolution(dis, _problem, _non_linear_solution, this);
	/**
	*Calculate the secondary variables on each node
	*/
	this->calc_nodal_eos_sys(0.0);
	

	// set initial output parameter
	OutputVariableInfo var1(this->getOutputParameterName(0), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _P);
	femData->outController.setOutput(var1.name, var1);
	OutputVariableInfo var2(this->getOutputParameterName(1), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _X);
	femData->outController.setOutput(var2.name, var2);
	// 
	OutputVariableInfo var_Sec_1(this->getOutputParameterName(2), _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _S);
	femData->outController.setOutput(var_Sec_1.name, var_Sec_1);

	OutputVariableInfo var_Sec_2(this->getOutputParameterName(3), _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _PL);
	femData->outController.setOutput(var_Sec_2.name, var_Sec_2);

	OutputVariableInfo var_Sec_3(this->getOutputParameterName(4), _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _PG);
	femData->outController.setOutput(var_Sec_3.name, var_Sec_3);

	OutputVariableInfo var_Sec_4(this->getOutputParameterName(5), _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _PC);
	femData->outController.setOutput(var_Sec_4.name, var_Sec_4);
	
	OutputVariableInfo var_Sec_5(this->getOutputParameterName(6), _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _rho_L_h);
	femData->outController.setOutput(var_Sec_5.name, var_Sec_5);

	OutputVariableInfo var_Sec_6(this->getOutputParameterName(7), _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _rho_G_h);
	femData->outController.setOutput(var_Sec_6.name, var_Sec_6);

	OutputVariableInfo var_Sec_7(this->getOutputParameterName(8), _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _rho_G_w);
	femData->outController.setOutput(var_Sec_7.name, var_Sec_7);
	return true;
}

template <class T1, class T2>
void FunctionCMP_TotalDensityForm<T1, T2>::initializeTimeStep(const NumLib::TimeStep &/*time*/)
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
void FunctionCMP_TotalDensityForm<T1, T2>::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{

}

template <class T1, class T2>
void FunctionCMP_TotalDensityForm<T1, T2>::output(const NumLib::TimeStep &/*time*/)
{
	
	//update data for output
	Ogs6FemData* femData = Ogs6FemData::getInstance();
	this->_msh_id = this->_problem->getDiscreteSystem()->getMesh()->getID();
	// set the new output
	// we have 2 primary variables
	OutputVariableInfo var1("MEAN_PRESSURE", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _P);
	femData->outController.setOutput(var1.name, var1);
	OutputVariableInfo var2("TOTAL_MASS_DENSITY", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _X);
	femData->outController.setOutput(var2.name, var2);
	// and 8 secondary variables
	// add all seconary variables as well. 
	//we have several secondary variables the values of which are defined on each elements
	
	OutputVariableInfo var_Sec_1("Saturation", _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _S);
	femData->outController.setOutput(var_Sec_1.name, var_Sec_1);

	OutputVariableInfo var_Sec_2("Liquid_Pressure", _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _PL);
	femData->outController.setOutput(var_Sec_2.name, var_Sec_2);

	OutputVariableInfo var_Sec_3("Gas_Pressure", _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _PG);
	femData->outController.setOutput(var_Sec_3.name, var_Sec_3);

	OutputVariableInfo var_Sec_4("Capillary_Pressure", _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _PC);
	femData->outController.setOutput(var_Sec_4.name, var_Sec_4);

	OutputVariableInfo var_Sec_5("Mass_Density_L_H", _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _rho_L_h);
	femData->outController.setOutput(var_Sec_5.name, var_Sec_5);

	OutputVariableInfo var_Sec_6("Mass_Density_G_H", _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _rho_G_h);
	femData->outController.setOutput(var_Sec_6.name, var_Sec_6);

	OutputVariableInfo var_Sec_7("Mass_Density_G_W", _msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _rho_G_w);
	femData->outController.setOutput(var_Sec_7.name, var_Sec_7);

}


template <class T1, class T2>
void FunctionCMP_TotalDensityForm<T1, T2>::calc_nodal_eos_sys(double dt = 0.0)
{
	

	
	/*
	// solve the EOS system on each node
	std::size_t node_id(0);
	double EPS = 1E-7;
	MathLib::LocalMatrix matSecDer = MathLib::LocalMatrix::Zero(3, 2);
	MathLib::LocalVector vecTempVar = MathLib::LocalVector::Zero(2);


	ogsChem::LocalVector output = ogsChem::LocalVector::Zero(10);
	ogsChem::LocalVector output_ini = output;
	ogsChem::LocalVector output_dP = ogsChem::LocalVector::Zero(3);
	ogsChem::LocalVector output_dX = ogsChem::LocalVector::Zero(3);
	ogsChem::LocalVector input = ogsChem::LocalVector::Zero(2);
	
	INFO("Calculating EOS on each node...");
	
	// loop for each node
	for (node_id = _P->getDiscreteData()->getRangeBegin();
		node_id < _P->getDiscreteData()->getRangeEnd();
		node_id++)
	{
		input(0) = _P->getValue(node_id);
		input(1) = _X->getValue(node_id);
		//This part is for standard newton iteration of local problem
		//solve EOS 
		//This part is for complementary condition
		output(0) = _S->getValue(node_id);
		output(1) = _rho_L_h->getValue(node_id);
		output(2) = _rho_G_h->getValue(node_id);
		output(3) = _PC->getValue(node_id);
		output(4) = _PL->getValue(node_id);
		output(5) = _PG->getValue(node_id);
		output(6) = _rho_G_w->getValue(node_id);
		output(7) = _dPGh_dPg->getValue(node_id);
		output(8) = _dPcdSg->getValue(node_id);
		output(9) = _Func_C->getValue(node_id);
		_EOS->set_env_condition(input);
		_LP_EOS->solve(input, output);
		
		//store the output value
		output_ini = output;

		_S->setValue(node_id, output(0));

		_rho_L_h->setValue(node_id, output(1));
		_rho_G_h->setValue(node_id, output(2));
		_PC->setValue(node_id, output(3));
		_PL->setValue(node_id, output(4));
		_PG->setValue(node_id, output(5));
		_rho_G_w->setValue(node_id, output(6));
		_dPGh_dPg->setValue(node_id, output(7));
		_dPcdSg->setValue(node_id, output(8));
		_Func_C->setValue(node_id, output(9));
		

		vecTempVar(0) = _EOS->Func_Omega(output(0));
		vecTempVar(1) = _EOS->Func_InterM(output(0));

		this->_vec_tempVar->setValue(node_id, vecTempVar);
		//set the values of derivatives
		if (this->m_EOS_EVL_METHOD == EOS_EVL_FIN_DIFF)
		{
			//CALCULATE THE DERIVATIVE OF THE SMALL VALUE ON P
			input(0) = input(0) + EPS;
			// directly use the output value as the initial guess of the next calculation
			_LP_EOS->solve(input, output);
			//store the output value
			output_dP = output;

			////CALCULATE THE DERIVATIVE OF THE SMALL VALUE ON X
			input(0) = _P->getValue(node_id);
			input(1) = input(1) + EPS;//add a small value on X
			_LP_EOS->solve(input, output);
			//store the output value
			output_dX = output;


			// here is the derivative operations. 
			matSecDer.col(0) = (output_dP - output_ini) / EPS;
			matSecDer.col(1) = (output_dX - output_ini) / EPS;
			//std::cout << matSecDer << std::endl;
		}
		else if (this->m_EOS_EVL_METHOD == EOS_EVL_ANALYTICAL)
		{
			//Compute the derivatives using the analytical way
			output_dP(0) = _EOS->Deriv_dSgdP(output(0));//dSgdP
			output_dX(0) = _EOS->Deriv_dSgdX(output(0));//dSgdX

			output_dP(1) = _EOS->Deriv_drhoLh_dP(output(0), vecTempVar(0));//
			output_dX(1) = _EOS->Deriv_drhoLh_dX(output(0), vecTempVar(1));

			output_dP(2) = _EOS->Deriv_drhoGh_dP(output(0), vecTempVar(0));
			output_dX(2) = _EOS->Deriv_drhoGh_dX(output(0), vecTempVar(1));

			matSecDer.col(0) = output_dP;
			matSecDer.col(1) = output_dX;
			//std::cout << matSecDer << std::endl;

		}

		this->_mat_secDer->setValue(node_id, matSecDer);

		// end
	}
	*/	
	

}


