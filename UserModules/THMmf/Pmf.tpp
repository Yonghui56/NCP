/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Pmf.tpp
 *
 * Created on 2012-10-24 by Norihiro Watanabe
 */

#include "logog.hpp"
#include "MeshLib/Tools/Tools.h"
#include "DiscreteLib/Utils/DiscreteSystemContainerPerMesh.h"
#include "FemLib/Function/FemNodalFunction.h"
#include "FemLib/Tools/FeObjectContainerPerFeType.h"
#include "FemLib/Tools/MeshElementToFemElementType.h"
#include "NumLib/Function/TXFunctionBuilder.h"
#include "OutputIO/OutputBuilder.h"
#include "OutputIO/OutputTimingBuilder.h"
#include "SolutionLib/Fem/FemSourceTerm.h"
#include "Ogs6FemData.h"
#include "FemVariableBuilder.h"
#include "THMmfFemElementCatalogBuilder.h"
#include "THMmfMeshElement2FeTypeBuilder.h"

namespace THMmf
{

template <class T1, class T2>
bool Pmf<T1,T2>::initialize(const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    size_t msh_id = option.getOptionAsNum<size_t>("MeshID");
    size_t time_id = option.getOptionAsNum<size_t>("TimeGroupID");
    NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];

    //mesh and FE objects
    MeshLib::IMesh* msh = femData->list_mesh[msh_id];
    MeshLib::setMeshElementCoordinatesMapping(*msh);
    MyDiscreteSystem* dis = 0;
    dis = DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(msh);
    FemLib::FemElementCatalog* feCatalog = THMmfFemElementCatalogBuilder::construct();
    FemLib::MeshElementToFemElementType* eleTofe = THMmfMeshElement2FeTypeBuilder::construct(*msh, FemLib::PolynomialOrder::Quadratic, femData->list_medium);
    _feObjects = new FemLib::FeObjectContainerPerFeType(feCatalog, eleTofe, msh);

    // local assemblers
    MyLinearAssemblerType* linear_assembler = new MyLinearAssemblerType();
    MyLinearAssemblerTypeForPorousMedia* linear_assembler_pm = new MyLinearAssemblerTypeForPorousMedia(_feObjects, msh->getGeometricProperty()->getCoordinateSystem());
    MyLinearAssemblerTypeForFracture* linear_assembler_frac = new MyLinearAssemblerTypeForFracture(_feObjects, msh->getGeometricProperty()->getCoordinateSystem());
    const MeshLib::CoordinateSystem coord = msh->getGeometricProperty()->getCoordinateSystem();
    linear_assembler->addLocalAssembler(coord.getDimension(), linear_assembler_pm);
    linear_assembler->addLocalAssembler(coord.getDimension()-1, linear_assembler_frac);

    // set up problem
    _problem = new MyProblemType(dis);
    MyEquationType* eqs = _problem->createEquation();
    eqs->initialize(linear_assembler);
    _problem->setTimeSteppingFunction(*tim);
    // set up variable
    typename MyProblemType::MyVariable* head = _problem->addVariable("pressure"); //internal name
    head->setFeObjectContainer(_feObjects);
    FemVariableBuilder varBuilder;
    varBuilder.doit(this->getOutputParameterName(Pressure), option, msh, femData->geo, femData->geo_unique_name, _feObjects, head);

    // set up solution
    _solution = new MySolutionType(dis, _problem);
    typename MySolutionType::LinearSolverType* linear_solver = _solution->getLinearEquationSolver();
    const BaseLib::Options* optNum = option.getSubGroup("Numerics");
    linear_solver->setOption(*optNum);
    _solution->getNonlinearSolver()->setOption(*optNum);

    // set initial output
    OutputVariableInfo var(this->getOutputParameterName(Pressure), msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));
    femData->outController.setOutput(var.name, var);

    // initial output parameter
    this->setOutput(Pressure, _solution->getCurrentSolution(head->getID()));

    return true;
}

template <class T1, class T2>
void Pmf<T1,T2>::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{
    setOutput(Pressure, _solution->getCurrentSolution(0));
}

template <class T1, class T2>
void Pmf<T1,T2>::output(const NumLib::TimeStep &/*time*/)
{
    //update data for output
    const size_t msh_id = _problem->getDiscreteSystem()->getMesh()->getID();
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    OutputVariableInfo var(this->getOutputParameterName(Pressure), msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));
    femData->outController.setOutput(var.name, var);
};

}
