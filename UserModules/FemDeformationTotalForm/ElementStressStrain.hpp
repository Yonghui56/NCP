/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file StressStrain.cpp
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

//#include "ElementStressStrain.h"

#include "logog.hpp"

#include "DiscreteLib/Utils/DiscreteSystemContainerPerMesh.h"
#include "Ogs6FemData.h"
#include "PhysicsLib/FemLinearElasticTools.h"

template <class T>
size_t FunctionElementStressStrain<T>::getNumberOfStrainComponents() const
{
    MeshLib::IMesh* msh = _dis->getMesh();
    const size_t dim = msh->getDimension();
    return (dim==2 ? 4 : 6);
}

template <class T>
bool FunctionElementStressStrain<T>::initialize(const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();

    size_t msh_id = option.getOptionAsNum<size_t>("MeshID");
    MeshLib::IMesh* msh = femData->list_mesh[msh_id];
    _dis = DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(msh);
    const size_t n_strain_components = getNumberOfStrainComponents();

    // create strain, stress vectors
    MathLib::LocalVector v0(n_strain_components);
    v0 *= .0;
    _strain = new MyIntegrationPointFunctionVector();
    _strain->initialize(_dis, v0);
    _stress = new MyIntegrationPointFunctionVector();
    _stress->initialize(_dis, v0);
    _feObjects = new FemLib::LagrangeFeObjectContainer(msh);

    // set initial output
    for (size_t i=0; i<n_strain_components; i++) {
        _vec_strain_components.push_back(new IntegrationPointScalarWrapper(_strain, i));
        _vec_stress_components.push_back(new IntegrationPointScalarWrapper(_stress, i));
    }
    for (size_t i=0; i<n_strain_components; i++) {
        OutputVariableInfo var1(this->getOutputParameterName(Strain) + getStressStrainComponentPostfix(i), msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _vec_strain_components[i]);
        femData->outController.setOutput(var1.name, var1);
        OutputVariableInfo var2(this->getOutputParameterName(Stress) + getStressStrainComponentPostfix(i), msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _vec_stress_components[i]);
        femData->outController.setOutput(var2.name, var2);
    }

    // initial output parameter
    setOutput(Strain, _strain);
    setOutput(Stress, _stress);

    return true;
}

template <class T>
void FunctionElementStressStrain<T>::finalizeTimeStep(const NumLib::TimeStep &/*time*/)
{
    //update data for output
    const size_t n_strain_components = getNumberOfStrainComponents();
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    const size_t msh_id = _dis->getMesh()->getID();
    for (size_t i=0; i<n_strain_components; i++) {
        OutputVariableInfo var1(this->getOutputParameterName(Strain) + getStressStrainComponentPostfix(i), msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _vec_strain_components[i]);
        femData->outController.setOutput(var1.name, var1);
        OutputVariableInfo var2(this->getOutputParameterName(Stress) + getStressStrainComponentPostfix(i), msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _vec_stress_components[i]);
        femData->outController.setOutput(var2.name, var2);
    }
};

template <class T>
int FunctionElementStressStrain<T>::solveTimeStep(const NumLib::TimeStep &/*time*/)
{
    INFO("Solving ELEMENT_STRESS_STRAIN...");

    const MeshLib::IMesh *msh = _dis->getMesh();
    MyNodalFunctionVector* u = (MyNodalFunctionVector*)getInput(Displacement);
    MyIntegrationPointFunctionVector* strain = _strain;
    MyIntegrationPointFunctionVector* stress = _stress;
    FemLib::LagrangeFeObjectContainer* feObjects = _feObjects;

    const size_t dim = msh->getDimension();
    const size_t n_strain_components = getNumberOfStrainComponents();

    Ogs6FemData* femData = Ogs6FemData::getInstance();

    //calculate strain, stress
    for (size_t i_e=0; i_e<msh->getNumberOfElements(); i_e++)
    {
        // element setup
        MeshLib::IElement* e = msh->getElement(i_e);
        const size_t nnodes = e->getNumberOfNodes();
        FemLib::IFiniteElement *fe = feObjects->getFeObject(*e);
        FemLib::IFemNumericalIntegration *integral = fe->getIntegrationMethod();
        const size_t n_gp = integral->getNumberOfSamplingPoints();
        strain->setNumberOfIntegationPoints(i_e, n_gp);
        stress->setNumberOfIntegationPoints(i_e, n_gp);
        const NumLib::TXPosition e_pos(NumLib::TXPosition::Element, e->getID());

        MaterialLib::Solid *solidphase = femData->list_solid[e->getGroupID()];
        // set D
        MathLib::LocalMatrix matD = MathLib::LocalMatrix::Zero(n_strain_components, n_strain_components);
        //matD *= .0;
        MathLib::LocalMatrix nv(1,1);
        MathLib::LocalMatrix E(1,1);
        solidphase->poisson_ratio->eval(e_pos, nv);
        solidphase->Youngs_modulus->eval(e_pos, E);
        double Lambda, G, K;
        MaterialLib::calculateLameConstant(nv(0,0), E(0,0), Lambda, G, K);
        MaterialLib::setElasticConsitutiveTensor(dim, Lambda, G, matD);

        // local u
        MathLib::LocalVector local_u(dim*nnodes);
        for (size_t j=0; j<nnodes; j++) {
            const size_t node_id = e->getNodeID(j);
            for (size_t k=0; k<dim; k++)
                local_u[j*dim+k] = u->getValue(node_id)(k);
        }

        // for each integration points
        MathLib::LocalMatrix matB = MathLib::LocalMatrix::Zero(n_strain_components, nnodes*dim);
        MathLib::LocalMatrix matN = MathLib::LocalMatrix::Zero(dim, nnodes*dim);
        double r[3] = {};
        double x[3] = {};
        for (size_t ip=0; ip<n_gp; ip++) {
            integral->getSamplingPoint(ip, r);
            fe->computeBasisFunctions(r);
            MathLib::LocalMatrix &N = *fe->getBasisFunction();
            const MathLib::LocalMatrix &dN = *fe->getGradBasisFunction();
            fe->getRealCoordinates(x);

            // set N,B
            setNu_Matrix_byPoint(dim, nnodes, N, matN);
            setB_Matrix_byPoint(dim, nnodes, dN, matB);

            // strain
            MathLib::LocalVector gp_strain(n_strain_components);
            gp_strain.noalias() = matB * local_u;
            strain->setIntegrationPointValue(i_e, ip, gp_strain);

            // stress
            MathLib::LocalVector gp_stress(n_strain_components);
            gp_stress = matD * gp_strain;
            stress->setIntegrationPointValue(i_e, ip, gp_stress);

//                std::cout << "strain=\n" << gp_strain << std::endl;
//                std::cout << "D=\n" << matD << std::endl;
//                std::cout << "stress=\n" << gp_stress << std::endl;
        }
    }

    setOutput(Strain, strain);
    setOutput(Stress, stress);

    return 0;
}

