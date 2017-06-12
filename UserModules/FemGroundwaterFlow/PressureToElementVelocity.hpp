/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementVelocity.hpp
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

//#include "ElementVelocity.h"

#include "logog.hpp"
#include "MathLib/Vector.h"
#include "DiscreteLib/Utils/DiscreteSystemContainerPerMesh.h"
#include "MaterialLib/PorousMedia.h"

#include "Ogs6FemData.h"

template <class T>
bool FunctionPressureToElementVelocity<T>::initialize(const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();

    size_t msh_id = option.getOptionAsNum<size_t>("MeshID");
    MeshLib::IMesh* msh = femData->list_mesh[msh_id];
    MyDiscreteSystem* dis = 0;
    dis = DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(msh);

    _feObjects = new FemLib::LagrangeFeObjectContainer(msh);
    _dis = dis;
    _vel = new MyIntegrationPointFunctionVector();
    _vel->initialize(dis);

    // set initial output
    OutputVariableInfo var(this->getOutputParameterName(Velocity), msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 3, _vel);
    femData->outController.setOutput(var.name, var);

    // initial output parameter
    _vel_3d = new My3DIntegrationPointFunctionVector(_vel, msh->getGeometricProperty()->getCoordinateSystem());
    this->setOutput(Velocity, _vel_3d);

    return true;
}

template <class T>
void FunctionPressureToElementVelocity<T>::finalizeTimeStep(const NumLib::TimeStep &/*time*/)
{
    //std::cout << "Velocity=" << std::endl;
    //_vel->printout();
    //update data for output
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    OutputVariableInfo var(this->getOutputParameterName(Velocity), _dis->getMesh()->getID(), OutputVariableInfo::Element, OutputVariableInfo::Real, 3, _vel_3d);
    femData->outController.setOutput(var.name, var);
};

template <class T>
int FunctionPressureToElementVelocity<T>::solveTimeStep(const NumLib::TimeStep &/*time*/)
{
    INFO("Calculating Darcy velocity within elements from fluid pressure...");

    const MeshLib::IMesh *msh = _dis->getMesh();
    const MeshLib::CoordinateSystem coord = msh->getGeometricProperty()->getCoordinateSystem();
    MyNodalFunctionScalar *f_p = (MyNodalFunctionScalar*)getInput(Pressure);
    MyIntegrationPointFunctionVector *vel = _vel;;

    MaterialLib::Fluid* fluidphase = Ogs6FemData::getInstance()->list_fluid[0];

    FemLib::LagrangeFeObjectContainer* feObjects = _feObjects;
    //calculate vel (vel=f(h))
    for (size_t i_e=0; i_e<msh->getNumberOfElements(); i_e++) {
        MeshLib::IElement* e = msh->getElement(i_e);
        const NumLib::TXPosition e_pos(NumLib::TXPosition::Element, e->getID());
        MeshLib::ElementCoordinatesMappingLocal* ele_local_coord;
        ele_local_coord = (MeshLib::ElementCoordinatesMappingLocal*)e->getMappedCoordinates();
        const MathLib::LocalMatrix &matR = ele_local_coord->getRotationMatrixToOriginal();
        size_t mat_id = e->getGroupID();
        MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];
        // fluid
        double mu = .0;
        fluidphase->dynamic_viscosity->eval(e_pos, mu);
        double rho_f = .0;
        MathLib::LocalVector vec_g;
        const bool hasGravityEffect = coord.hasZ();
        if (hasGravityEffect) {
            fluidphase->density->eval(e_pos, rho_f);
            vec_g = MathLib::LocalVector::Zero(coord.getDimension());
            vec_g[coord.getIndexOfZ()] = -9.81;
        }
        // media
        double k;
        pm->permeability->eval(e_pos, k);
        double k_mu;
        k_mu = k / mu;
        MathLib::LocalMatrix local_k_mu = MathLib::LocalMatrix::Identity(e->getDimension(), e->getDimension());
        local_k_mu *= k_mu;
        MathLib::LocalMatrix global_k_mu;
        if (e->getDimension() < coord.getDimension()) {
            MathLib::LocalMatrix local2 = MathLib::LocalMatrix::Zero(coord.getDimension(), coord.getDimension());
//            local2.topLeftCorner(local_k_mu.rows(), local_k_mu.cols()) = local_k_mu;
            local2.block(0, 0, local_k_mu.rows(), local_k_mu.cols()) = local_k_mu.block(0, 0, local_k_mu.rows(), local_k_mu.cols());
            global_k_mu = matR * local2 * matR.transpose();
        } else {
            global_k_mu = local_k_mu;
        }

        FemLib::IFiniteElement *fe = feObjects->getFeObject(*e);
        MathLib::LocalVector local_p(e->getNumberOfNodes());
        for (size_t j=0; j<e->getNumberOfNodes(); j++)
            local_p[j] = f_p->getValue(e->getNodeID(j));
        // for each integration points
        FemLib::IFemNumericalIntegration *integral = fe->getIntegrationMethod();
        double r[3] = {};
        const size_t n_gp = integral->getNumberOfSamplingPoints();
        vel->setNumberOfIntegationPoints(i_e, n_gp);
        MathLib::LocalVector local_q;
        for (size_t ip=0; ip<n_gp; ip++) {
            integral->getSamplingPoint(ip, r);
            fe->computeBasisFunctions(r);
            const MathLib::LocalMatrix* dN = fe->getGradBasisFunction();
            //MathLib::LocalMatrix* N = fe->getBasisFunction();

            //TODO casting cause zero q
            //static_cast<MathLib::LocalVector>(q.head(msh->getDimension())) = (*dN) * local_p * (-1.0) * k_mu;
            local_q = - global_k_mu * (*dN) * local_p;
            if (hasGravityEffect) {
                // F += dNp^T * K * rho * gz
                local_q.noalias() += global_k_mu * rho_f * vec_g;
            }

            vel->setIntegrationPointValue(i_e, ip, local_q);
			
        }
    }
    setOutput(Velocity, vel);
    return 0;
}

