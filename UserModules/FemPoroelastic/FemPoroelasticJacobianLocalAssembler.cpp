/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemPoroelasticJacobianLocalAssembler.cpp
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#include "FemPoroelasticJacobianLocalAssembler.h"

#include "FemLib/Core/Element/IFemElement.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Solid.h"
#include "PhysicsLib/DeformationTools.h"
#include "PhysicsLib/FemLinearElasticTools.h"
#include "Ogs6FemData.h"

void FemPoroelasticJacobianLocalAssembler::assembleComponents
    (  const NumLib::TimeStep &timestep,
        const MeshLib::IElement &e,
        const std::vector<size_t> &vec_order,
        const std::vector<LocalVectorType> &/*vec_x0*/,
        const std::vector<LocalVectorType> &/*vec_x1*/,
        std::vector<std::vector<LocalMatrixType> > &vec_K
        )
{
    assert(vec_order.size()==3);

    //TODO how do you know 1st is ux and 3rd is p? who decides it?
    const size_t u_order = vec_order[0];
    assert(u_order==vec_order[1]);
    const size_t p_order = vec_order[2];
    // ------------------------------------------------------------------------
    // Element
    // ------------------------------------------------------------------------
    const size_t dim = e.getDimension();
    const size_t n_strain_components = getNumberOfStrainComponents(dim);
    const size_t nnodes_u = e.getNumberOfNodes(u_order);
    const size_t nnodes_p = e.getNumberOfNodes(p_order);
    const NumLib::TXPosition e_pos(NumLib::TXPosition::Element, e.getID());

    // ------------------------------------------------------------------------
    // Transient
    // ------------------------------------------------------------------------
    const double dt = timestep.getTimeStepSize();
    const double theta = 1.0; //TODO

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
    double mu = .0;
    fluidphase->dynamic_viscosity->eval(e_pos, mu);

    // media
    double k;
    pm->permeability->eval(e_pos, k);
//    double n = .0;
//    pm->porosity->eval(e_pos, n);
    double s = .0;
    pm->storage->eval(e_pos, s);
    double k_mu;
    k_mu = k / mu;


    // ------------------------------------------------------------------------
    // Local component assembly
    // ------------------------------------------------------------------------
    LocalMatrixType Kuu = LocalMatrixType::Zero(nnodes_u*dim, nnodes_u*dim);
    LocalMatrixType Cup = LocalMatrixType::Zero(nnodes_u*dim, nnodes_p);
    LocalMatrixType Kpp = LocalMatrixType::Zero(nnodes_p, nnodes_p);
    LocalMatrixType Mpp = LocalMatrixType::Zero(nnodes_p, nnodes_p);
    LocalMatrixType Cpu = LocalMatrixType::Zero(nnodes_p, nnodes_u*dim);

    // temp matrix
    LocalMatrixType B = LocalMatrixType::Zero(n_strain_components, nnodes_u*dim);
    //LocalMatrixType Nuvw = LocalMatrixType::Zero(dim, nnodes_u*dim);
    const LocalMatrixType m = get_m(dim);

    //
    FemLib::IFiniteElement* fe_u = _feObjects.getFeObject(e, u_order);
    FemLib::IFiniteElement* fe_p = _feObjects.getFeObject(e, p_order);
    FemLib::IFemNumericalIntegration *q_u = fe_u->getIntegrationMethod();
    double gp_x[3], real_x[3];
    for (size_t j=0; j<q_u->getNumberOfSamplingPoints(); j++) {
        q_u->getSamplingPoint(j, gp_x);
        fe_u->computeBasisFunctions(gp_x);
        fe_p->computeBasisFunctions(gp_x);
        fe_u->getRealCoordinates(real_x);
        double fac_u = fe_u->getDetJ() * q_u->getWeight(j);

        //--- local component ----
        // set N,B
        //LocalMatrixType &Nu = *fe_u->getBasisFunction();
        LocalMatrixType &dNu = *fe_u->getGradBasisFunction();
        //setNu_Matrix_byComponent(dim, nnodes_u, Nu, Nuvw);
        setB_Matrix_byComponent(dim, nnodes_u, dNu, B);
        LocalMatrixType &Np = *fe_p->getBasisFunction();

        // K_uu += B^T * D * B
        Kuu.noalias() += fac_u * B.transpose() * De * B;

        // C_up += B^T * m * Np
        Cup.noalias() += fac_u * B.transpose() * m * Np;
    }

    FemLib::IFemNumericalIntegration *q_p = fe_p->getIntegrationMethod();
    for (size_t j=0; j<q_p->getNumberOfSamplingPoints(); j++) {
        q_p->getSamplingPoint(j, gp_x);
        fe_u->computeBasisFunctions(gp_x);
        fe_p->computeBasisFunctions(gp_x);
        fe_p->getRealCoordinates(real_x);
        double fac = fe_p->getDetJ() * q_p->getWeight(j);

        //--- local component ----
        // set N,B
        LocalMatrixType &dNu = *fe_u->getGradBasisFunction();
        setB_Matrix_byComponent(dim, nnodes_u, dNu, B);
        LocalMatrixType &Np = *fe_p->getBasisFunction();
        LocalMatrixType &dNp = *fe_p->getGradBasisFunction();

        // M_pp += Np^T * S * Np
        Mpp.noalias() += fac * Np.transpose() * s * Np;

        // K_pp += dNp^T * K * dNp
        Kpp.noalias() += fac * dNp.transpose() * k_mu * dNp;

        // C_pu += Np^T * m^T * B
        Cpu.noalias() += fac * Np.transpose() * m.transpose() * B;
    }

    //
    for (size_t i=0; i<dim; i++) {
        for (size_t j=0; j<dim; j++) {
            vec_K[i][j] = Kuu.block(i*nnodes_u, j*nnodes_u, nnodes_u, nnodes_u);
        }
        vec_K[i][dim] = -1. * Cup.block(i*nnodes_u, 0, nnodes_u, nnodes_p);
        vec_K[dim][i] = 1.0/dt * Cpu.block(0, i*nnodes_u, nnodes_p, nnodes_u);
    }
    vec_K[dim][dim] = 1.0/dt * Mpp + theta * Kpp ;

}
