/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Tmf.tpp
 *
 * Created on 2012-10-23 by Norihiro Watanabe
 */


#include "logog.hpp"

#include "MeshLib/Tools/Tools.h"
#include "DiscreteLib/Utils/DiscreteSystemContainerPerMesh.h"
#include "FemLib/Function/FemNodalFunction.h"
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
bool Tmf<T1,T2>::initialize(const BaseLib::Options &option)
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

    // equations
    MyLinearAssemblerType* linear_assembler = new MyLinearAssemblerType();
    MyLinearAssemblerTypeForPorousMedia* linear_assembler_pm = new MyLinearAssemblerTypeForPorousMedia(_feObjects);
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
    typename MyProblemType::MyVariable* var_T = _problem->addVariable("T");
    var_T->setFeObjectContainer(_feObjects);
    FemVariableBuilder varBuilder;
    varBuilder.doit(this->getOutputParameterName(Temperature), option, msh, femData->geo, femData->geo_unique_name, _feObjects, var_T);

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
    this->setOutput(Temperature, _solution->getCurrentSolution(var_T->getID()));

    return true;
}

template <class T1, class T2>
void Tmf<T1, T2>::initializeTimeStep(const NumLib::TimeStep &/*time*/)
{
    const NumLib::ITXFunction *vel = this->getInput<NumLib::ITXFunction>(Velocity);
    MyLinearAssemblerType* dim_assembers = this->_problem->getEquation()->getLinearAssembler();
    const MeshLib::CoordinateSystem coord = this->_problem->getDiscreteSystem()->getMesh()->getGeometricProperty()->getCoordinateSystem();
    MyLinearAssemblerTypeForPorousMedia* ap = static_cast<MyLinearAssemblerTypeForPorousMedia*>(dim_assembers->getLocalAssembler(coord.getDimension()));
    MyLinearAssemblerTypeForFracture* af = static_cast<MyLinearAssemblerTypeForFracture*>(dim_assembers->getLocalAssembler(coord.getDimension()-1));
    ap->setVelocity(vel);
    af->setVelocity(vel);
}

template <class T1, class T2>
void Tmf<T1, T2>::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{
    setOutput(Temperature, _solution->getCurrentSolution(0));
}

template <class T1, class T2>
void Tmf<T1, T2>::output(const NumLib::TimeStep &/*time*/)
{
    //update data for output
    const size_t msh_id = _problem->getDiscreteSystem()->getMesh()->getID();
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    OutputVariableInfo var(this->getOutputParameterName(Temperature), msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));
    femData->outController.setOutput(var.name, var); 
};

} // end THMmf

