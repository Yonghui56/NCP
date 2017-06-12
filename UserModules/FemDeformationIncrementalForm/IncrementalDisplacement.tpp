/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IncrementalDisplacement.tpp
 *
 * Created on 2012-11-29 by Norihiro Watanabe
 */

#include "logog.hpp"

#include "DiscreteLib/Core/IDiscreteVector.h"
#include "FemLib/Function/FemNodalFunction.h"
#include "NumLib/Function/TXFunctionBuilder.h"
#include "OutputIO/OutputBuilder.h"
#include "OutputIO/OutputTimingBuilder.h"
#include "SolutionLib/Fem/FemSourceTerm.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Solid.h"
#include "PhysicsLib/FemLinearElasticTools.h"
#include "PhysicsLib/SmallDeformationMedia.h"
#include "Ogs6FemData.h"
#include "FemVariableBuilder.h"
#include "IncrementalDeformationTools.h"


template<class T1, class T2>
bool FunctionIncrementalDisplacement<T1, T2>::initialize(
        const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    const size_t msh_id = option.getOptionAsNum<size_t>("MeshID");
    const size_t time_id = option.getOptionAsNum<size_t>("TimeGroupID");
    MeshLib::IMesh* msh = femData->list_mesh[msh_id];
    const size_t dim = msh->getDimension();
    const MeshLib::CoordinateSystem msh_coord =
            msh->getGeometricProperty()->getCoordinateSystem();
    NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];

    //--------------------------------------------------------------------------
    // set up higher-order mesh and FE type catalog
    //--------------------------------------------------------------------------
    const FemLib::PolynomialOrder::type msh_order = FemLib::PolynomialOrder::Linear; //TODO user input
    if (msh_order != FemLib::PolynomialOrder::Linear) {
        INFO("->generating higher order mesh...");
        MeshLib::MeshGenerator::generateHigherOrderUnstrucuredMesh(*(MeshLib::UnstructuredMesh*)msh, 2);
        INFO("* mesh id %d: order=%d, nodes=%d, elements=%d", msh_id, msh->getMaxiumOrder(), msh->getNumberOfNodes(msh->getMaxiumOrder()), msh->getNumberOfElements());
    }
    _feObjects = new FemLib::LagrangeFeObjectContainer(msh);

    //--------------------------------------------------------------------------
    // create a discrete system
    //--------------------------------------------------------------------------
    MyDiscreteSystem* dis =
            DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(msh);

    //--------------------------------------------------------------------------
    // set up FEM IVBV problem
    //--------------------------------------------------------------------------
    _problem = new MyProblemType(dis);
    _problem->setTimeSteppingFunction(*tim);

    //--------------------------------------------------------------------------
    // set up problem variables
    //--------------------------------------------------------------------------
    MyVariable* u_x = _problem->addVariable("u_x", msh_order);
    MyVariable* u_y = _problem->addVariable("u_y", msh_order);
    MyVariable* u_z = NULL;
    if (dim==3) {
        u_z = _problem->addVariable("u_z", msh_order);
    }
    u_x->setFeObjectContainer(_feObjects);
    u_y->setFeObjectContainer(_feObjects);
    if (dim==3) {
        u_z->setFeObjectContainer(_feObjects);
    }
    FemVariableBuilder var_builder;
    var_builder.doit(this->getOutputParameterName(Displacement) + "_X1", option,
            msh, femData->geo, femData->geo_unique_name, _feObjects, u_x);
    var_builder.doit(this->getOutputParameterName(Displacement) + "_Y1", option,
            msh, femData->geo, femData->geo_unique_name, _feObjects, u_y);
    if (dim==3) {
        var_builder.doit(this->getOutputParameterName(Displacement)+"_Z1", option,
                msh, femData->geo, femData->geo_unique_name, _feObjects, u_z);
    }

    //--------------------------------------------------------------------------
    // set up local assemblers and equations
    //--------------------------------------------------------------------------
    MyLinearAssemblerType* linear_assembler = new MyLinearAssemblerType(
            _feObjects);
    MyResidualAssemblerType* r_assembler = new MyResidualAssemblerType(
            _feObjects, msh_coord);
    MyJacobianAssemblerType* j_eqs = new MyJacobianAssemblerType(_feObjects);
    MyEquationType* eqs = _problem->createEquation();
    eqs->initialize(linear_assembler, r_assembler, j_eqs);

    //--------------------------------------------------------------------------
    // set up this solution algorithm
    //--------------------------------------------------------------------------
    _solution = new MySolutionType(dis, _problem);
    typename MySolutionType::LinearSolverType* linear_solver =
            _solution->getLinearEquationSolver();
    const BaseLib::Options* optNum = option.getSubGroup("Numerics");
    linear_solver->setOption(*optNum);
    _solution->getNonlinearSolver()->setOption(*optNum);
    DiscreteLib::DofEquationIdTable* dofMapping =
            _solution->getDofEquationIdTable();
    dofMapping->setNumberingType(DiscreteLib::DofNumberingType::BY_POINT);
    dofMapping->setLocalNumberingType(DiscreteLib::DofNumberingType::BY_POINT);

    //--------------------------------------------------------------------------
    // set up secondary variables: displacement, stress, strain
    //--------------------------------------------------------------------------
    MathLib::LocalVector tmp_u0(3);
    tmp_u0 *= .0;
    _previous_displacement = new MyNodalFunctionVector();
    _previous_displacement->initialize(*dis, u_x->getCurrentOrder(), tmp_u0);
    _current_displacement = new MyNodalFunctionVector(*_previous_displacement);
    _delta_displacement = new MyNodalFunctionVector(*_previous_displacement);
    // create strain, stress vectors
    const size_t n_strain_components = getNumberOfStrainComponents(dim);
    MathLib::LocalVector v0 = MathLib::LocalVector::Zero(n_strain_components);
    _previous_strain = new MyIntegrationPointFunctionVector();
    _previous_strain->initialize(dis, v0);
    _current_strain = new MyIntegrationPointFunctionVector(*_previous_strain);
    _previous_stress = new MyIntegrationPointFunctionVector(*_previous_strain);
    _current_stress = new MyIntegrationPointFunctionVector(*_previous_stress);


    //--------------------------------------------------------------------------
    // set up data for output
    //--------------------------------------------------------------------------
    for (size_t i=0; i<dim; i++) {
        _vec_u_components.push_back(new NodalPointScalarWrapper(_current_displacement, i));
    }
    for (size_t i=0; i<n_strain_components; i++) {
        _vec_strain_components.push_back(new IntegrationPointScalarWrapper(_current_strain, i));
        _vec_stress_components.push_back(new IntegrationPointScalarWrapper(_current_stress, i));
    }

    // set initial output
    OutputVariableInfo var(this->getOutputParameterName(Displacement), msh_id,
            OutputVariableInfo::Node, OutputVariableInfo::Real, 2,
            _current_displacement);
    femData->outController.setOutput(var.name, var);
    for (size_t i = 0; i < _vec_u_components.size(); i++) {
        OutputVariableInfo var1(
                this->getOutputParameterName(Displacement)
                        + getDisplacementComponentPostfix(i), msh_id,
                OutputVariableInfo::Node, OutputVariableInfo::Real, 1,
                _vec_u_components[i]);
        femData->outController.setOutput(var1.name, var1);
    }
    for (size_t i=0; i<n_strain_components; i++) {
        OutputVariableInfo var1(this->getOutputParameterName(Strain) + getStressStrainComponentPostfix(i), msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _vec_strain_components[i]);
        femData->outController.setOutput(var1.name, var1);
        OutputVariableInfo var2(this->getOutputParameterName(Stress) + getStressStrainComponentPostfix(i), msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _vec_stress_components[i]);
        femData->outController.setOutput(var2.name, var2);
    }

    //--------------------------------------------------------------------------
    // set initial system output
    //--------------------------------------------------------------------------
    this->setOutput(Displacement, _current_displacement);
    this->setOutput(Strain, _current_strain);
    this->setOutput(Stress, _current_stress);

    return true;
}

template<class T1, class T2>
void FunctionIncrementalDisplacement<T1, T2>::initializeTimeStep(
        const NumLib::TimeStep &/*time*/)
{
    INFO("->pre-processing before calling a solution algorithm...");
    this->_problem->getEquation()->getResidualAssembler()->setPreviousStress(_previous_stress);
}

template<class T1, class T2>
void FunctionIncrementalDisplacement<T1, T2>::postSolutionAlgorithm(const NumLib::TimeStep &/*time*/)
{
    INFO("->post-processing after calling a solution algorithm...");
    MeshLib::IMesh* msh = _problem->getDiscreteSystem()->getMesh();
    const size_t dim = msh->getDimension();
    const size_t n_u_nodes = _delta_displacement->getNumberOfNodes();
    //--------------------------------------------------------------------------
    // compute current total displacement
    //--------------------------------------------------------------------------
    for (size_t i=0; i<n_u_nodes; i++) {
        for (size_t j=0; j<dim; j++) {
            _delta_displacement->getValue(i)(j) = _solution->getCurrentSolution(j)->getValue(i);
        }
    }
    *_current_displacement->getDiscreteData() = *_previous_displacement->getDiscreteData();
    *_current_displacement->getDiscreteData() +=  *_delta_displacement->getDiscreteData();
    //--------------------------------------------------------------------------
    // compute strain, stress
    //--------------------------------------------------------------------------
    calculateStrainStress(Ogs6FemData::getInstance(), msh, _feObjects,_delta_displacement, _previous_strain, _previous_stress, _current_strain, _current_stress);
}

template<class T1, class T2>
void FunctionIncrementalDisplacement<T1, T2>::postTimeStep(const NumLib::TimeStep &/*time*/)
{
    INFO("->post-processing after solving time step...");

    // copy current value to previous value
    *_previous_displacement = *_current_displacement;
    *_previous_strain = *_current_strain;
    *_previous_stress = *_current_stress;
}

template<class T1, class T2>
void FunctionIncrementalDisplacement<T1, T2>::updateOutputParameter(
        const NumLib::TimeStep &/*time*/)
{
    //--------------------------------------------------------------------------
    // set system output
    //--------------------------------------------------------------------------
    setOutput(Displacement, _current_displacement);
    setOutput(Strain, _current_strain);
    setOutput(Stress, _current_stress);
}


template<class T1, class T2>
void FunctionIncrementalDisplacement<T1, T2>::output(
        const NumLib::TimeStep &/*time*/)
{
    const size_t msh_id = _problem->getDiscreteSystem()->getMesh()->getID();
    //update data for output
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    OutputVariableInfo var(this->getOutputParameterName(Displacement), msh_id,
            OutputVariableInfo::Node, OutputVariableInfo::Real, 2,
            _current_displacement);
    femData->outController.setOutput(var.name, var);
    for (size_t i = 0; i < _vec_u_components.size(); i++) {
        OutputVariableInfo var1(
                this->getOutputParameterName(Displacement)
                        + getDisplacementComponentPostfix(i), msh_id,
                OutputVariableInfo::Node, OutputVariableInfo::Real, 1,
                _vec_u_components[i]);
        femData->outController.setOutput(var1.name, var1);
    }
}



