/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Umf.tpp
 *
 * Created on 2012-11-15 by Norihiro Watanabe
 */


#include "logog.hpp"

#include "MeshLib/Tools/MeshGenerator.h"
#include "FemLib/Function/FemNodalFunction.h"
#include "NumLib/Function/TXFunctionBuilder.h"
#include "OutputIO/OutputBuilder.h"
#include "OutputIO/OutputTimingBuilder.h"
#include "SolutionLib/Fem/FemSourceTerm.h"
#include "SolutionLib/Fem/FemVariable.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Solid.h"
#include "PhysicsLib/FemLinearElasticTools.h"
#include "Ogs6FemData.h"
#include "FemVariableBuilder.h"
#include "THMmfFemElementCatalogBuilder.h"
#include "THMmfMeshElement2FeTypeBuilder.h"

namespace THMmf
{

template <class T1, class T2>
typename Umf<T1,T2>::MyVariable* Umf<T1,T2>::getDisplacementComponent(MyVariable *u_x, MyVariable* u_y, MyVariable* u_z, const std::string &var_name)
{
    if (var_name.find("_X")!=std::string::npos) {
        return u_x;
    } else if (var_name.find("_Y")!=std::string::npos) {
        return u_y;
    } else {
        return u_z;
    }
}

template <class T1, class T2>
size_t Umf<T1,T2>::getDisplacementComponentIndex(const std::string &var_name) const
{
    if (var_name.find("_X")!=std::string::npos) {
        return 0;
    } else if (var_name.find("_Y")!=std::string::npos) {
        return 1;
    } else {
        return 2;
    }
}

template <class T1, class T2>
bool Umf<T1,T2>::initialize(const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    size_t msh_id = option.getOptionAsNum<size_t>("MeshID");
    size_t time_id = option.getOptionAsNum<size_t>("TimeGroupID");

    //--------------------------------------------------------------------------
    // set up mesh and FE objects
    //--------------------------------------------------------------------------
    MeshLib::IMesh* msh = femData->list_mesh[msh_id];
    const FemLib::PolynomialOrder::type msh_order = FemLib::PolynomialOrder::Linear;
    if (msh_order != FemLib::PolynomialOrder::Linear) {
        INFO("->generating higher order mesh...");
        MeshLib::MeshGenerator::generateHigherOrderUnstrucuredMesh(*(MeshLib::UnstructuredMesh*)msh, 2);
        INFO("* mesh id %d: order=%d, nodes=%d, elements=%d", msh_id, msh->getMaxiumOrder(), msh->getNumberOfNodes(msh->getMaxiumOrder()), msh->getNumberOfElements());
    }
    MyDiscreteSystem* dis = 0;
    dis = DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(msh);
    FemLib::FemElementCatalog* feCatalog = THMmfFemElementCatalogBuilder::construct();
    FemLib::MeshElementToFemElementType* eleTofe = THMmfMeshElement2FeTypeBuilder::construct(*msh, FemLib::PolynomialOrder::Quadratic, femData->list_medium);
    _feObjects = new FemLib::FeObjectContainerPerFeType(feCatalog, eleTofe, msh);
    const size_t dim = msh->getDimension();

    //--------------------------------------------------------------------------
    // set up problem
    //--------------------------------------------------------------------------
    _problem = new MyProblemType(dis);
    NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];
    _problem->setTimeSteppingFunction(*tim);

    //--------------------------------------------------------------------------
    // set up variables
    //--------------------------------------------------------------------------
    // definitions
    MyVariable* u_x = _problem->addVariable("u_x", msh_order);
    u_x->setFeObjectContainer(_feObjects);
    MyVariable* u_y = _problem->addVariable("u_y", msh_order);
    u_y->setFeObjectContainer(_feObjects);
    MyVariable* u_z = NULL;
    if (dim==3) {
        u_z = _problem->addVariable("u_z", msh_order);
        u_z->setFeObjectContainer(_feObjects);
    }
    FemVariableBuilder var_builder;
    var_builder.doit(this->getOutputParameterName(Displacement)+"_X1", option, msh, femData->geo, femData->geo_unique_name, _feObjects, u_x);
    var_builder.doit(this->getOutputParameterName(Displacement)+"_Y1", option, msh, femData->geo, femData->geo_unique_name, _feObjects, u_y);
    if (dim==3) {
        var_builder.doit(this->getOutputParameterName(Displacement)+"_Z1", option, msh, femData->geo, femData->geo_unique_name, _feObjects, u_z);
    }

    //--------------------------------------------------------------------------
    // set up equations
    //--------------------------------------------------------------------------
    MyEquationType* eqs = _problem->createEquation();
    std::vector<size_t> vec_orders;
    for (size_t i=0; i<_problem->getNumberOfVariables(); i++)
        vec_orders.push_back(_problem->getVariable(i)->getCurrentOrder());

    MyLinearAssemblerType* linear_assembler = new MyLinearAssemblerType();
    MyLinearAssemblerTypeForPorousMedia* linear_assembler_pm = new MyLinearAssemblerTypeForPorousMedia(_feObjects, _problem->getNumberOfVariables(), vec_orders, msh->getGeometricProperty()->getCoordinateSystem());
    //MyLinearAssemblerTypeForFracture* linear_assembler_frac = new MyLinearAssemblerTypeForFracture(*_feObjects, msh->getGeometricProperty()->getCoordinateSystem());
    const MeshLib::CoordinateSystem coord = msh->getGeometricProperty()->getCoordinateSystem();
    linear_assembler->addLocalAssembler(coord.getDimension(), linear_assembler_pm);
    //linear_assembler->addLocalAssembler(coord.getDimension()-1, linear_assembler_frac);
    eqs->initialize(linear_assembler);


    //--------------------------------------------------------------------------
    // set up this solution algorithm
    //--------------------------------------------------------------------------
    _solution = new MySolutionType(dis, _problem);
    typename MySolutionType::LinearSolverType* linear_solver = _solution->getLinearEquationSolver();
    const BaseLib::Options* optNum = option.getSubGroup("Numerics");
    linear_solver->setOption(*optNum);
    _solution->getNonlinearSolver()->setOption(*optNum);
    _solution->getDofEquationIdTable()->setNumberingType(DiscreteLib::DofNumberingType::BY_POINT);

    //--------------------------------------------------------------------------
    // set up data for output
    //--------------------------------------------------------------------------
    // create u variable which is vector
    MathLib::LocalVector tmp_u0(3);
    tmp_u0 *= .0;
    _displacement = new MyNodalFunctionVector();
    _displacement->initialize(*dis, u_x->getCurrentOrder(), tmp_u0);
    for (size_t i=0; i<_displacement->getNumberOfNodes(); i++) {
        for (size_t j=0; j<dim; j++) {
            _displacement->getValue(i)(j) = _solution->getCurrentSolution(j)->getValue(i);
        }
    }
    for (size_t i=0; i<dim; i++) {
        _vec_u_components.push_back(new NodalPointScalarWrapper(_displacement, i));
    }

    // set initial output
    OutputVariableInfo var(this->getOutputParameterName(Displacement), msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 2, _displacement);
    femData->outController.setOutput(var.name, var);
    for (size_t i=0; i<_vec_u_components.size(); i++) {
        OutputVariableInfo var1(this->getOutputParameterName(Displacement) + getDisplacementComponentPostfix(i), msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _vec_u_components[i]);
        femData->outController.setOutput(var1.name, var1);
    }

    //--------------------------------------------------------------------------
    // set initial output
    //--------------------------------------------------------------------------
    // initial output parameter
    this->setOutput(Displacement, _displacement);

    std::cout << "At the end of Umf::init(): x_max = " << femData->geo->getSurfaceVec(femData->geo_unique_name)->at(2)->getAABB().getMaxPoint()[0] << std::endl;

    return true;
}

template <class T1, class T2>
void Umf<T1,T2>::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{
    const size_t dim = _problem->getDiscreteSystem()->getMesh()->getDimension();
    for (size_t i=0; i<_displacement->getNumberOfNodes(); i++) {
        for (size_t j=0; j<dim; j++) {
            _displacement->getValue(i)(j) = _solution->getCurrentSolution(j)->getValue(i);
        }
    }
    setOutput(Displacement, _displacement);
}

template <class T1, class T2>
void Umf<T1,T2>::output(const NumLib::TimeStep &/*time*/)
{
    //update data for output
    const size_t msh_id = _problem->getDiscreteSystem()->getMesh()->getID();
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    OutputVariableInfo var(this->getOutputParameterName(Displacement), msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 2, _displacement);
    femData->outController.setOutput(var.name, var);
    for (size_t i=0; i<_vec_u_components.size(); i++) {
        OutputVariableInfo var1(this->getOutputParameterName(Displacement) + getDisplacementComponentPostfix(i), msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _vec_u_components[i]);
        femData->outController.setOutput(var1.name, var1);
    }
};

}
