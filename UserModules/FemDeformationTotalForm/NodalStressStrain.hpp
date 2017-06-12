/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NodalStressStrain.cpp
 *
 * Created on 2012-08-14 by Norihiro Watanabe
 */

//#include "NodalStressStrain.h"

#include "logog.hpp"

#include "DiscreteLib/Utils/DiscreteSystemContainerPerMesh.h"
#include "FemLib/Post/Extrapolation.h"

#include "Ogs6FemData.h"
#include "PhysicsLib/FemLinearElasticTools.h"

template <class T>
size_t FunctionNodalStressStrain<T>::getNumberOfStrainComponents() const
{
    MeshLib::IMesh* msh = _dis->getMesh();
    const size_t dim = msh->getDimension();
    return (dim==2 ? 4 : 6);
}

template <class T>
bool FunctionNodalStressStrain<T>::initialize(const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();

    size_t msh_id = option.getOptionAsNum<size_t>("MeshID");
    MeshLib::IMesh* msh = femData->list_mesh[msh_id];
    _dis = DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(msh);
    const size_t n_strain_components = getNumberOfStrainComponents();
    _feObjects = new FemLib::LagrangeFeObjectContainer(msh);

    // create strain, stress vectors
    MathLib::LocalVector v0(n_strain_components);
    v0 *= .0;

    // set initial output
    _nodal_strain = new MyNodalFunctionVector();
    _nodal_strain->initialize(*_dis, FemLib::PolynomialOrder::Linear, v0);
    _nodal_stress = new MyNodalFunctionVector();
    _nodal_stress->initialize(*_dis, FemLib::PolynomialOrder::Linear, v0);
    for (size_t i=0; i<n_strain_components; i++) {
        _vec_nodal_strain_components.push_back(new NodalPointScalarWrapper(_nodal_strain, i));
        _vec_nodal_stress_components.push_back(new NodalPointScalarWrapper(_nodal_stress, i));
    }
    for (size_t i=0; i<n_strain_components; i++) {
        OutputVariableInfo var1(this->getOutputParameterName(NodStrain) + getStressStrainComponentPostfix(i), msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _vec_nodal_strain_components[i]);
        femData->outController.setOutput(var1.name, var1);
        OutputVariableInfo var2(this->getOutputParameterName(NodStress) + getStressStrainComponentPostfix(i), msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _vec_nodal_stress_components[i]);
        femData->outController.setOutput(var2.name, var2);
    }

    // initial output parameter
    setOutput(NodStrain, _nodal_strain);
    setOutput(NodStress, _nodal_stress);

    return true;
}

template <class T>
void FunctionNodalStressStrain<T>::finalizeTimeStep(const NumLib::TimeStep &/*time*/)
{
    //update data for output
    const size_t n_strain_components = getNumberOfStrainComponents();
    const size_t msh_id = _dis->getMesh()->getID();
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    for (size_t i=0; i<n_strain_components; i++) {
        OutputVariableInfo var1(this->getOutputParameterName(NodStrain) + getStressStrainComponentPostfix(i), msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _vec_nodal_strain_components[i]);
        femData->outController.setOutput(var1.name, var1);
        OutputVariableInfo var2(this->getOutputParameterName(NodStress) + getStressStrainComponentPostfix(i), msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _vec_nodal_stress_components[i]);
        femData->outController.setOutput(var2.name, var2);
    }
};

template <class T>
int FunctionNodalStressStrain<T>::solveTimeStep(const NumLib::TimeStep &/*time*/)
{
    INFO("Solving NODAL_STRESS_STRAIN...");

    MyIntegrationPointFunctionVector* strain = (MyIntegrationPointFunctionVector*)getInput(GpStrain);
    MyIntegrationPointFunctionVector* stress = (MyIntegrationPointFunctionVector*)getInput(GpStress);

    FemLib::FemExtrapolationAverage<MyDiscreteSystem, MathLib::LocalVector> extrapo;
    _nodal_strain->setFeObjectContainer(_feObjects);
    _nodal_stress->setFeObjectContainer(_feObjects);
    extrapo.extrapolate(*strain, *_nodal_strain);
    extrapo.extrapolate(*stress, *_nodal_stress);

    setOutput(NodStrain, _nodal_strain);
    setOutput(NodStress, _nodal_stress);

    return 0;
}

