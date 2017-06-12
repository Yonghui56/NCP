/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Head.hpp
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

//#include "Head.h"
#include "logog.hpp"
#include "DiscreteLib/Utils/DiscreteSystemContainerPerMesh.h"
#include "FemLib/Function/FemNodalFunction.h"
#include "NumLib/Function/TXFunctionBuilder.h"
#include "OutputIO/OutputBuilder.h"
#include "OutputIO/OutputTimingBuilder.h"
#include "SolutionLib/Fem/FemSourceTerm.h"
#include "SolutionLib/Fem/FemIC.h"
#include "Ogs6FemData.h"
#include "FemVariableBuilder.h"

template <class T1, class T2>
bool FunctionHead<T1,T2>::initialize(const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    size_t msh_id = option.getOptionAsNum<size_t>("MeshID");
    size_t time_id = option.getOptionAsNum<size_t>("TimeGroupID");

    //--------------------------------------------------------------------------
    // set up mesh and FE objects
    //--------------------------------------------------------------------------
    MeshLib::IMesh* msh = femData->list_mesh[msh_id];
    MyDiscreteSystem* dis = 0;
    dis = DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(msh);
    _feObjects = new FemLib::LagrangeFeObjectContainer(msh);

    //--------------------------------------------------------------------------
    // set up problem
    //--------------------------------------------------------------------------
    _problem = new MyProblemType(dis);
    NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];
    _problem->setTimeSteppingFunction(*tim);

    //--------------------------------------------------------------------------
    // set up variables
    //--------------------------------------------------------------------------
    MyVariable* head = _problem->addVariable("head"); //internal name
    FemVariableBuilder var_builder;
    var_builder.doit(this->getOutputParameterName(Head), option, msh, femData->geo, femData->geo_unique_name, _feObjects, head);

    //--------------------------------------------------------------------------
    // set up equations
    //--------------------------------------------------------------------------
    // local assemblers
    MyLinearAssemblerType* linear_assembler = new MyLinearAssemblerType(*_feObjects);
    MyResidualAssemblerType* r_assembler = new MyResidualAssemblerType(*_feObjects);
    MyJacobianAssemblerType* j_eqs = new MyJacobianAssemblerType(*_feObjects);
    // equations
    MyEquationType* eqs = _problem->createEquation();
    eqs->initialize(linear_assembler, r_assembler, j_eqs);

    //--------------------------------------------------------------------------
    // set up this solution algorithm
    //--------------------------------------------------------------------------
    _solution = new MySolutionType(dis, _problem);
    typename MySolutionType::LinearSolverType* linear_solver = _solution->getLinearEquationSolver();
    const BaseLib::Options* optNum = option.getSubGroup("Numerics");
    linear_solver->setOption(*optNum);
    _solution->getNonlinearSolver()->setOption(*optNum);

    //--------------------------------------------------------------------------
    // set up data for output
    //--------------------------------------------------------------------------
    OutputVariableInfo var(this->getOutputParameterName(Head), msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));
    femData->outController.setOutput(var.name, var);

    //--------------------------------------------------------------------------
    // set initial output
    //--------------------------------------------------------------------------
    this->setOutput(Head, _solution->getCurrentSolution(0));

    return true;
}

template <class T1, class T2>
void FunctionHead<T1,T2>::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{
    setOutput(Head, _solution->getCurrentSolution(0));
}

template <class T1, class T2>
void FunctionHead<T1,T2>::output(const NumLib::TimeStep &/*time*/)
{
    //update data for output
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    OutputVariableInfo var(this->getOutputParameterName(Head), _problem->getDiscreteSystem()->getMesh()->getID(), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));
    femData->outController.setOutput(var.name, var);
};
