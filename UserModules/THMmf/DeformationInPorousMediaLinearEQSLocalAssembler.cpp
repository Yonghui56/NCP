/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DeformationInPorousMediaLinearEQSLocalAssembler.cpp
 *
 * Created on 2012-11-15 by Norihiro Watanabe
 */

#include "DeformationInPorousMediaLinearEQSLocalAssembler.h"

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/TransientAssembler/IElementWiseTransientLinearEQSLocalAssembler.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Solid.h"
#include "MaterialLib/Fluid.h"
#include "PhysicsLib/DeformationTools.h"
#include "PhysicsLib/FemLinearElasticTools.h"
#include "Ogs6FemData.h"

namespace THMmf
{

void DeformationInPorousMediaLinearEQSLocalAssembler::assembleComponents(
        const NumLib::TimeStep &/*timestep*/,
        const MeshLib::IElement &e, 
        const std::vector<size_t> &vec_order, 
        const std::vector<LocalVectorType> &vec_x0, 
        const std::vector<LocalVectorType> &/*vec_u1*/, 
        std::vector<std::vector<LocalMatrixType> > &vec_K,
        std::vector<LocalVectorType> &vec_F
        )
{
    const size_t dim = e.getDimension();

    assert(vec_order.size()==dim);

    const size_t u_order = vec_order[0];
    const size_t u_nnodes = vec_x0[0].rows();
    assert(u_order==vec_order[1]);
    // combine ui
    LocalVectorType u0(u_nnodes*dim);
    for (size_t i=0; i<u_nnodes; i++) {
        for (size_t j=0; j<dim; j++) {
            u0(i+j*u_nnodes) = vec_x0[j](i);
        }
    }
    // ------------------------------------------------------------------------
    // Element
    // ------------------------------------------------------------------------
    const size_t n_strain_components = getNumberOfStrainComponents(dim);
    const size_t nnodes_u = e.getNumberOfNodes(u_order);
    const NumLib::TXPosition e_pos(NumLib::TXPosition::Element, e.getID());

    // ------------------------------------------------------------------------
    // Transient
    // ------------------------------------------------------------------------
    //const double dt = timestep.getTimeStepSize();
    const double theta = 1.0;

    // ------------------------------------------------------------------------
    // Material (assuming element constant)
    // ------------------------------------------------------------------------
    size_t mat_id = e.getGroupID();
    size_t fluid_id = 0; //TODO
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    MaterialLib::PorousMedia* pm = femData->list_pm[mat_id];
    MaterialLib::Solid *solidphase = femData->list_solid[mat_id];
    MaterialLib::Fluid *fluidphase = femData->list_fluid[fluid_id];

    // solid
    double rho_s = .0;
    if (solidphase->density!=NULL)
        solidphase->density->eval(e_pos, rho_s);
    LocalMatrixType De = LocalMatrixType::Zero(n_strain_components, n_strain_components);
    MathLib::LocalMatrix nv(1,1);
    MathLib::LocalMatrix E(1,1);
    solidphase->poisson_ratio->eval(e_pos, nv);
    solidphase->Youngs_modulus->eval(e_pos, E);
    double Lambda, G, K;
    MaterialLib::calculateLameConstant(nv(0,0), E(0,0), Lambda, G, K);
    MaterialLib::setElasticConsitutiveTensor(dim, Lambda, G, De);

    // fluid
    double rho_f = .0;
    fluidphase->density->eval(e_pos, rho_f);

    // media
    double n = .0;
    pm->porosity->eval(e_pos, n);

    double pm_rho = (1.-n) * rho_s + n * rho_f;


    // ------------------------------------------------------------------------
    // Body force
    // ------------------------------------------------------------------------
    LocalVectorType body_force = LocalVectorType::Zero(dim);
    const bool hasGravity = _problem_coordinates.hasZ();
    if (hasGravity) {
        body_force[_problem_coordinates.getIndexOfZ()] = pm_rho * 9.81;
    }

    // ------------------------------------------------------------------------
    // Local component assembly
    // ------------------------------------------------------------------------
    LocalMatrixType Kuu = LocalMatrixType::Zero(nnodes_u*dim, nnodes_u*dim);
    LocalVectorType Fu = LocalVectorType::Zero(nnodes_u*dim);

    // temp matrix
    LocalMatrixType B = LocalMatrixType::Zero(n_strain_components, nnodes_u*dim);
    LocalMatrixType Nuvw = LocalMatrixType::Zero(dim, nnodes_u*dim);
    const LocalMatrixType m = get_m(dim);

    //
    e.setCurrentOrder(u_order);
    FemLib::IFiniteElement* fe_u = _feObjects->getFeObject(e);
    FemLib::IFemNumericalIntegration *q_u = fe_u->getIntegrationMethod();
    double gp_x[3], real_x[3];
    for (size_t j=0; j<q_u->getNumberOfSamplingPoints(); j++) {
        q_u->getSamplingPoint(j, gp_x);
        fe_u->computeBasisFunctions(gp_x);
        fe_u->getRealCoordinates(real_x);
        double fac_u = fe_u->getDetJ() * q_u->getWeight(j);

        //--- local component ----
        // set N,B
        LocalMatrixType &Nu = *fe_u->getBasisFunction();
        LocalMatrixType &dNu = *fe_u->getGradBasisFunction();
        setNu_Matrix_byComponent(dim, nnodes_u, Nu, Nuvw);
        setB_Matrix_byComponent(dim, nnodes_u, dNu, B);

        // K_uu += B^T * D * B
        Kuu.noalias() += fac_u * B.transpose() * De * B;

        // Fu += N^T * b
        if (hasGravity)
            Fu.noalias() += fac_u * Nuvw.transpose() * body_force;
    }
    Fu.noalias() += (theta - 1) * Kuu * u0;

    //
    for (size_t i=0; i<dim; i++) {
        for (size_t j=0; j<dim; j++) {
            vec_K[i][j] = Kuu.block(i*nnodes_u, j*nnodes_u, nnodes_u, nnodes_u);
        }
    }

    for (size_t i=0; i<dim; i++) {
        vec_F[i] = Fu.segment(i*nnodes_u, nnodes_u);
    }
}

}

