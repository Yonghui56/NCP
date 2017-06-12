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
#include "LocalProblem_EOS.h"

template <class T1, class T2>
bool FunctionCMP_2P2C<T1, T2>::initialize(const BaseLib::Options &option)
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
	this->setOutputParameterName(0, "MEAN_PRESSURE");
	this->setOutputParameterName(1, "MOLAR_FRACTION");
	// also secondary variables
	// TODO: set seconary variable names also as output parameters

	// create the MyCompMultiPhase problem
	_problem = new MyCompMultiPhaseProblemType(dis);
	_problem->setTimeSteppingFunction(*tim);  // applying the time stepping function
	
	// creating mean pressure vector
	MyVariableCompMultiPhase* mean_pressure = _problem->addVariable("MEAN_PRESSURE");
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
	MyVariableCompMultiPhase* molar_fraction = _problem->addVariable("MOLAR_FRACTION");
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
	_LP_EOS = new LocalProblem_EOS();
	//ogsChem::LocalVector _output1;
	//_output1 = ogsChem::LocalVector::Zero(_LP_EOS->N);

	_PC = new MyNodalFunctionScalar();	
	_PC->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
	//_PC->setValue(0, 0.0); //PC_ini
	_PL = new MyNodalFunctionScalar();
	_PL->initialize(*dis, FemLib::PolynomialOrder::Linear, 100.0);//PL_ini
	_PG = new MyNodalFunctionScalar();
	_PG->initialize(*dis, FemLib::PolynomialOrder::Linear, 100.0);
	_X_m = new MyNodalFunctionScalar();
	_X_m->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
	_X_M = new MyNodalFunctionScalar();
	_X_M->initialize(*dis, FemLib::PolynomialOrder::Linear,0.0);
	_NG = new MyNodalFunctionScalar();
	_NG->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
	_NL = new MyNodalFunctionScalar();
	_NL->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
	_S = new MyNodalFunctionScalar();
	_S->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
	
	_mat_secDer = new MyNodalFunctionMatrix();
	MathLib::LocalMatrix tmp_mat = MathLib::LocalMatrix::Zero(8, 2); 
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
	MyNonLinearAssemblerType* nonlin_assembler = new MyNonLinearAssemblerType(_feObjects, this);
	MyNonLinearResidualAssemblerType* non_linear_r_assembler = new MyNonLinearResidualAssemblerType(_feObjects, this);
	MyNonLinearJacobianAssemblerType* non_linear_j_assembler = new MyNonLinearJacobianAssemblerType(_feObjects, this);

	// define nonlinear problem
	_non_linear_problem = new MyNonLinearCompMultiPhaseProblemType(dis);
	_non_linear_eqs = _non_linear_problem->createEquation();
	_non_linear_eqs->initialize(nonlin_assembler, non_linear_r_assembler, non_linear_j_assembler);
	_non_linear_problem->setTimeSteppingFunction(*tim);
	_non_linear_problem->addVariable("MEAN_PRESSURE");
	_non_linear_problem->addVariable("MOLAR_FRACTION");

	SolutionLib::FemIC* P_ic = new SolutionLib::FemIC(msh);
	P_ic->addDistribution(femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>(_P->getDiscreteData()));
	var_builder.doit("MEAN_PRESSURE", option, msh, femData->geo, femData->geo_unique_name, _feObjects, _non_linear_problem->getVariable(0));
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
	_solution = new MyCompMultiPhaseSolution(dis, _problem, _non_linear_solution, this);
	/**
	*Calculate the secondary variables on each node
	*/
	this->calc_nodal_eos_sys(0.0);
	

	// set initial output parameter
	OutputVariableInfo var1(this->getOutputParameterName(0), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _P);
	femData->outController.setOutput(var1.name, var1);
	OutputVariableInfo var2(this->getOutputParameterName(1), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _X);
	femData->outController.setOutput(var2.name, var2);

	return true;
}

template <class T1, class T2>
void FunctionCMP_2P2C<T1, T2>::initializeTimeStep(const NumLib::TimeStep &/*time*/)
{
	/*size_t i;
	//const NumLib::ITXFunction *vel = this->getInput<NumLib::ITXFunction>(Velocity);

	// set velocity for linear problem
	for (i = 0; i < _linear_problems.size(); i++) {
		_linear_problems[i]->getEquation()->getLinearAssembler()->setVelocity(vel);
		_linear_problems[i]->getEquation()->getResidualAssembler()->setVelocity(vel);
		_linear_problems[i]->getEquation()->getJacobianAssembler()->setVelocity(vel);
	}*/
}

template <class T1, class T2>
void FunctionCMP_2P2C<T1, T2>::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{

}

template <class T1, class T2>
void FunctionCMP_2P2C<T1, T2>::output(const NumLib::TimeStep &/*time*/)
{
	// update data for output
	Ogs6FemData* femData = Ogs6FemData::getInstance();
	this->_msh_id = this->_problem->getDiscreteSystem()->getMesh()->getID();
	// set the new output
	// we have 2 primary variables
	OutputVariableInfo var1("MEAN_PRESSURE", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _P);
	femData->outController.setOutput(var1.name, var1);
	OutputVariableInfo var2("MOLAR_FRACTION", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _X);
	femData->outController.setOutput(var2.name, var2);
	// and 8 secondary variables
	// add all seconary variables as well. 
	OutputVariableInfo var_Sec_1("Capillary Pressure", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _PC);
	femData->outController.setOutput(var_Sec_1.name, var_Sec_1);
	OutputVariableInfo var_Sec_2("Liquid_Phase_Pressure", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _PL);
	femData->outController.setOutput(var_Sec_2.name, var_Sec_2);
	OutputVariableInfo var_Sec_3("Gas_Phase_Pressure", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _PG);
	femData->outController.setOutput(var_Sec_3.name, var_Sec_3);
	OutputVariableInfo var_Sec_4("Molar_Fraction_Light_component_Liquid_Phase", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _X_m);
	femData->outController.setOutput(var_Sec_4.name, var_Sec_4);
	OutputVariableInfo var_Sec_5("Molar_Fraction_Light_component_Gas_Phase", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _X_M);
	femData->outController.setOutput(var_Sec_5.name, var_Sec_5);
	OutputVariableInfo var_Sec_6("Molar_Density_Gas_Phase", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _NG);
	femData->outController.setOutput(var_Sec_6.name, var_Sec_6);
	OutputVariableInfo var_Sec_7("Molar_Density_Liquid_Phase", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _NL);
	femData->outController.setOutput(var_Sec_7.name, var_Sec_7);
	OutputVariableInfo var_Sec_8("Saturation", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _S);
	femData->outController.setOutput(var_Sec_8.name, var_Sec_8);
}


template <class T1, class T2>
void FunctionCMP_2P2C<T1, T2>::calc_nodal_eos_sys(double dt = 0.0)
{
	// solve the EOS system on each node
	std::size_t node_id(0); 
	double EPS = 1E-7;
	MathLib::LocalMatrix matSecDer = MathLib::LocalMatrix::Zero(8, 2); 

	ogsChem::LocalVector output = ogsChem::LocalVector::Zero(8);
	ogsChem::LocalVector output_ini = output;
	ogsChem::LocalVector output_dP = output;
	ogsChem::LocalVector output_dX = output;
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
		/*output(0) = _PC->getValue(node_id);
		output(1) = _PL->getValue(node_id);
		output(2) = _PG->getValue(node_id);
		output(3) = _X_m->getValue(node_id);
		output(4) = _X_M->getValue(node_id);
		output(5) = _NG->getValue(node_id);
		output(6) = _NL->getValue(node_id);
		output(7) = _S->getValue(node_id);*/
		//This part is for complementary condition
		output(0) = _S->getValue(node_id);
		output(1) = _X_m->getValue(node_id);
		output(2) = _X_M->getValue(node_id);
		output(3) = _PC->getValue(node_id);
		output(4) = _PL->getValue(node_id);
		output(5) = _PG->getValue(node_id);
		output(6) = _NL->getValue(node_id);
		output(7) = _NG->getValue(node_id);
		_LP_EOS->solve(input, output);
		//store the output value
		output_ini = output;
		//pass the value
		/*_PC->setValue(node_id, output(0));
		_PL->setValue(node_id, output(1));
		_PG->setValue(node_id, output(2));
		_X_m->setValue(node_id, output(3));
		_X_M->setValue(node_id, output(4));
		_NG->setValue(node_id, output(5));
		_NL->setValue(node_id, output(6));
		_S->setValue(node_id, output(7));*/
		_S->setValue(node_id, output(0));
		_X_m->setValue(node_id, output(1));
		_X_M->setValue(node_id, output(2));
		_PC->setValue(node_id, output(3));
		_PL->setValue(node_id, output(4));
		_PG->setValue(node_id, output(5));
		_NL->setValue(node_id, output(6));
		_NG->setValue(node_id, output(7));
		
		std::size_t flag_phase = this->_EOS->Characteristic_Function_xi(output(0));

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


			// TODO: here is the derivative operations. 
			matSecDer.col(0) = (output_dP - output_ini) / EPS;
			matSecDer.col(1) = (output_dX - output_ini) / EPS;
			//std::cout << matSecDer << std::endl;
		}
		else if (this->m_EOS_EVL_METHOD == EOS_EVL_ANALYTICAL)
		{
			//Compute the derivatives using the analytical way
			output_dP(0) = this->_EOS->Calculate_dSgdP(input(0), input(1), output(0), output(4), output(5));//dsgdp
			output_dP(3) = this->_EOS->Calculate_dPCdP(output_dP(0));//dpcdp
			output_dP(4) = this->_EOS->Calculate_dPLdP(output(0), output(3), output_dP(0), output_dP(3));//dpldp
			output_dP(5) = this->_EOS->Calculate_dPGdP(output(0), output(3), output_dP(0), output_dP(3));//dpgdp
			
			output_dP(1) = this->_EOS->Calculate_dX_mdP(output(5), output_dP(5));//dXmdp
			output_dP(2) = this->_EOS->Calculate_dXMdP(output(5), output_dP(5));//dXMdp
			output_dP(6) = this->_EOS->Calculate_dN_LdP(output(1), output(4), output_dP(1), output_dP(4));//dNLdp
			output_dP(7) = this->_EOS->Calculate_dN_GdP(output_dP(5));//dNGdp


			output_dX(0) = this->_EOS->Calculate_dSgdX(input(0), input(1), output(0), output(4), output(5));//dsgdx
			output_dX(3) = this->_EOS->Calculate_dPCdX(output_dX(0));//dpcdx
			output_dX(4) = this->_EOS->Calculate_dPLdX(output(0), output(3), output_dX(0), output_dX(3));//dpldx
			output_dX(5) = this->_EOS->Calculate_dPGdX(output(0), output(3), output_dX(0), output_dX(3));//dpgdx

			output_dX(1) = this->_EOS->Calculate_dX_mdX(output(5), output_dX(5));//dXmdx
			output_dX(2) = this->_EOS->Calculate_dXMdX(output(5), output_dX(5));//dXMdx
			output_dX(6) = this->_EOS->Calculate_dN_LdX(output(1), output(4), output_dX(1), output_dX(4));//dNLdx
			output_dX(7) = this->_EOS->Calculate_dN_GdX(output_dX(5));//dNGdx


			matSecDer.col(0) = output_dP;
			matSecDer.col(1) = output_dX;
			//std::cout << matSecDer << std::endl;
			
		}
		
		this->_mat_secDer->setValue(node_id, matSecDer); 

		// end
	}

}


