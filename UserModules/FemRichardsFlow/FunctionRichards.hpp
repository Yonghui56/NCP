/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Concentration.cpp
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */


//#include "Concentration.h"

#include "logog.hpp"

#include "DiscreteLib/Utils/DiscreteSystemContainerPerMesh.h"
#include "FemLib/Function/FemNodalFunction.h"
#include "NumLib/Function/TXFunctionBuilder.h"
#include "OutputIO/OutputBuilder.h"
#include "OutputIO/OutputTimingBuilder.h"
#include "SolutionLib/Fem/FemSourceTerm.h"
#include "Ogs6FemData.h"
#include "FemVariableBuilder.h"

template <class T1, class T2>
bool FunctionRichards<T1,T2>::initialize(const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    size_t msh_id = option.getOptionAsNum<size_t>("MeshID");
    size_t time_id = option.getOptionAsNum<size_t>("TimeGroupID");
    NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];

    //mesh and FE objects
    MeshLib::IMesh* msh = femData->list_mesh[msh_id];
    MyDiscreteSystem* dis = 0;
    dis = DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(msh);
    _feObjects = new FemLib::LagrangeFeObjectContainer(msh);

    ////Compound
    //size_t compound_id = option.getOptionAsNum<size_t>("CompoundID");
    //if (femData->list_compound.size() < compound_id) {
    //    ERR("***Error: compound data not found (compound id=%d)", compound_id);
    //    return false;
    //}
    //_compound = femData->list_compound[compound_id];

    // local assemblers
    MyLinearAssemblerType* linear_assembler = new MyLinearAssemblerType(_feObjects, msh->getGeometricProperty()->getCoordinateSystem());
    MyResidualAssemblerType* r_assembler = new MyResidualAssemblerType(_feObjects, msh->getGeometricProperty()->getCoordinateSystem());
    MyJacobianAssemblerType* j_eqs = new MyJacobianAssemblerType(_feObjects);
	
    // set up problem
    _problem = new MyProblemType(dis);
    MyEquationType* eqs = _problem->createEquation();
    eqs->initialize(linear_assembler, r_assembler, j_eqs);
    _problem->setTimeSteppingFunction(*tim);
    // set up variable
    typename MyProblemType::MyVariable* pressure_w = _problem->addVariable("PRESSURE1");
    FemVariableBuilder varBuilder;
    varBuilder.doit("PRESSURE1", option, msh, femData->geo, femData->geo_unique_name, _feObjects, pressure_w);

    // set up solution
    _solution = new MySolutionType(dis, _problem);
    MyLinearSolver* linear_solver = _solution->getLinearEquationSolver();
    const BaseLib::Options* optNum = option.getSubGroup("Numerics");
    linear_solver->setOption(*optNum);
    _solution->getNonlinearSolver()->setOption(*optNum);

    // set initial output
    OutputVariableInfo var(this->getOutputParameterName(Pressure), msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));
    //OutputVariableInfo var(this->getOutputParameterName(Pressure), _problem->getDiscreteSystem()->getMesh()->getID(), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));

	femData->outController.setOutput(var.name, var); 

    // initial output parameter
    this->setOutput(Pressure, _solution->getCurrentSolution(pressure_w->getID()));

    return true;
}

template <class T1, class T2>
void FunctionRichards<T1, T2>::initializeTimeStep(const NumLib::TimeStep &/*time*/)
{
  
}

template <class T1, class T2>
void FunctionRichards<T1, T2>::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{
    setOutput(Pressure, _solution->getCurrentSolution(0));
}

template <class T1, class T2>
void FunctionRichards<T1, T2>::output(const NumLib::TimeStep &/*time*/)
{
    //update data for output
    const size_t msh_id = _problem->getDiscreteSystem()->getMesh()->getID();
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    OutputVariableInfo var(this->getOutputParameterName(Pressure), msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));
    femData->outController.setOutput(var.name, var); 


};
