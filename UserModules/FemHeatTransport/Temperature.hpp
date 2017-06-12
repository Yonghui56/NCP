/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Temperature.cpp
 *
 * Created on 2012-10-23 by Norihiro Watanabe
 */


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
bool FunctionTemperature<T1,T2>::initialize(const BaseLib::Options &option)
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

    // equations
    MyLinearAssemblerType* linear_assembler = new MyLinearAssemblerType(_feObjects);
    MyResidualAssemblerType* r_assembler = new MyResidualAssemblerType(_feObjects);
    MyJacobianAssemblerType* j_eqs = new MyJacobianAssemblerType(_feObjects);

    // set up problem
    _problem = new MyProblemType(dis);
    MyEquationType* eqs = _problem->createEquation();
    eqs->initialize(linear_assembler, r_assembler, j_eqs);
    _problem->setTimeSteppingFunction(*tim);
    // set up variable
    typename MyProblemType::MyVariable* concentration = _problem->addVariable("concentration");
    FemVariableBuilder varBuilder;
    varBuilder.doit(this->getOutputParameterName(Temperature), option, msh, femData->geo, femData->geo_unique_name, _feObjects, concentration);

    // set up solution
    _solution = new MySolutionType(dis, _problem);
    MyLinearSolver* linear_solver = _solution->getLinearEquationSolver();
    const BaseLib::Options* optNum = option.getSubGroup("Numerics");
    linear_solver->setOption(*optNum);
    _solution->getNonlinearSolver()->setOption(*optNum);

    // set initial output
    OutputVariableInfo var(this->getOutputParameterName(Temperature), msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));
    femData->outController.setOutput(var.name, var); 

    // initial output parameter
    this->setOutput(Temperature, _solution->getCurrentSolution(concentration->getID()));

    return true;
}

template <class T1, class T2>
void FunctionTemperature<T1, T2>::initializeTimeStep(const NumLib::TimeStep &/*time*/)
{
    const NumLib::ITXFunction *vel = this->getInput<NumLib::ITXFunction>(Velocity);
    this->_problem->getEquation()->getLinearAssembler()->setVelocity(vel);
    this->_problem->getEquation()->getResidualAssembler()->setVelocity(vel);
    this->_problem->getEquation()->getJacobianAssembler()->setVelocity(vel);
}

template <class T1, class T2>
void FunctionTemperature<T1, T2>::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{
    setOutput(Temperature, _solution->getCurrentSolution(0));
}

template <class T1, class T2>
void FunctionTemperature<T1, T2>::output(const NumLib::TimeStep &/*time*/)
{
    //update data for output
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    OutputVariableInfo var(this->getOutputParameterName(Temperature), _problem->getDiscreteSystem()->getMesh()->getID(), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));
    femData->outController.setOutput(var.name, var); 
};
