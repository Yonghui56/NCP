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
#include "Ogs6FemData.h"
#include "FemVariableBuilder.h"

template <class T1, class T2>
bool FunctionLiquidPressure<T1,T2>::initialize(const BaseLib::Options &option)
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

    // local assemblers
    MyLinearAssemblerType* linear_assembler = new MyLinearAssemblerType(*_feObjects, msh->getGeometricProperty()->getCoordinateSystem());
    MyResidualAssemblerType* r_assembler = new MyResidualAssemblerType(*_feObjects, msh->getGeometricProperty()->getCoordinateSystem());
    MyJacobianAssemblerType* j_eqs = new MyJacobianAssemblerType(*_feObjects);

    // set up problem
    _problem = new MyProblemType(dis);
    MyEquationType* eqs = _problem->createEquation();
    eqs->initialize(linear_assembler, r_assembler, j_eqs);
    _problem->setTimeSteppingFunction(*tim);
    // set up variable
    typename MyProblemType::MyVariable* head = _problem->addVariable("pressure"); //internal name
    FemVariableBuilder varBuilder;
    varBuilder.doit(this->getOutputParameterName(Pressure), option, msh, femData->geo, femData->geo_unique_name, _feObjects, head);

    // set up solution
    _solution = new MySolutionType(dis, _problem);
    typename MySolutionType::LinearSolverType* linear_solver = _solution->getLinearEquationSolver();
    const BaseLib::Options* optNum = option.getSubGroup("Numerics");
    linear_solver->setOption(*optNum);
    _solution->getNonlinearSolver()->setOption(*optNum);

    // set initial output
    OutputVariableInfo var(this->getOutputParameterName(Pressure), _problem->getDiscreteSystem()->getMesh()->getID(), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));
    femData->outController.setOutput(var.name, var);

    // initial output parameter
    this->setOutput(Pressure, _solution->getCurrentSolution(head->getID()));

    return true;
}

template <class T1, class T2>
void FunctionLiquidPressure<T1,T2>::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{
    setOutput(Pressure, _solution->getCurrentSolution(0));
}

template <class T1, class T2>
void FunctionLiquidPressure<T1,T2>::output(const NumLib::TimeStep &/*time*/)
{
    //update data for output
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    OutputVariableInfo var(this->getOutputParameterName(Pressure), _problem->getDiscreteSystem()->getMesh()->getID(), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));
    femData->outController.setOutput(var.name, var);
};
