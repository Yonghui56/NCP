/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemPoroelasticResidualLocalAssembler.cpp
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */


#include "FemPoroelasticResidualLocalAssembler.h"

#include "FemLib/Core/Element/IFemElement.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Solid.h"
#include "PhysicsLib/DeformationTools.h"
#include "PhysicsLib/FemLinearElasticTools.h"
#include "Ogs6FemData.h"

void FemPoroelasticResidualLocalAssembler::assembleComponents
    (   const NumLib::TimeStep &timestep,
        const MeshLib::IElement &e, 
        const std::vector<size_t> &vec_order, 
        const std::vector<LocalVectorType> &vec_x0,
        const std::vector<LocalVectorType> &vec_x1,
        std::vector<LocalVectorType> &vec_r
        )
{
    assert(vec_order.size()==3);

    //TODO how do you know 1st is ux and 3rd is p? who decides it?
    const size_t id_ux = 0;
    const size_t id_uy = 1;
    const size_t id_p = 2;

    const size_t u_order = vec_order[id_ux];
    assert(u_order==vec_order[id_uy]);
    const size_t p_order = vec_order[id_p];
    const LocalVectorType &ux0 = vec_x0[id_ux];
    const LocalVectorType &uy0 = vec_x0[id_uy];
    const LocalVectorType &ux1 = vec_x1[id_ux];
    const LocalVectorType &uy1 = vec_x1[id_uy];
    // combine ux and uy
    LocalVectorType u0(ux0.rows()*2);
    LocalVectorType u1(ux0.rows()*2);
    for (int i=0; i<ux0.rows(); i++) {
        u0(i) = ux0(i);
        u0(i+ux0.rows()) = uy0(i);
        u1(i) = ux1(i);
        u1(i+ux0.rows()) = uy1(i);
    }
    const LocalVectorType &p0 = vec_x0[id_p];
    const LocalVectorType &p1 = vec_x1[id_p];
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
    double mu = .0;
    fluidphase->dynamic_viscosity->eval(e_pos, mu);
    double rho_f = .0;
    fluidphase->density->eval(e_pos, rho_f);

    // media
    double k;
    pm->permeability->eval(e_pos, k);
    double n = .0;
    pm->porosity->eval(e_pos, n);
    double s = .0;
    pm->storage->eval(e_pos, s);
    double k_mu;
    k_mu = k / mu;


    // ------------------------------------------------------------------------
    // Body force
    // ------------------------------------------------------------------------
    LocalVectorType body_force = LocalVectorType::Zero(dim);
    bool hasGravity = false;
    if (hasGravity) {
        body_force[dim-1] = rho_s * 9.81;
    }

    // ------------------------------------------------------------------------
    // Local component assembly
    // ------------------------------------------------------------------------
    LocalMatrixType Kuu = LocalMatrixType::Zero(nnodes_u*dim, nnodes_u*dim);
    LocalMatrixType Cup = LocalMatrixType::Zero(nnodes_u*dim, nnodes_p);
    LocalMatrixType Kpp = LocalMatrixType::Zero(nnodes_p, nnodes_p);
    LocalMatrixType Mpp = LocalMatrixType::Zero(nnodes_p, nnodes_p);
    LocalMatrixType Cpu = LocalMatrixType::Zero(nnodes_p, nnodes_u*dim);
    LocalVectorType Fu = LocalVectorType::Zero(nnodes_u*dim);
    LocalVectorType Fp = LocalVectorType::Zero(nnodes_p);

    // temp matrix
    LocalMatrixType B = LocalMatrixType::Zero(n_strain_components, nnodes_u*dim);
    LocalMatrixType Nuvw = LocalMatrixType::Zero(dim, nnodes_u*dim);
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
        LocalMatrixType &Nu = *fe_u->getBasisFunction();
        LocalMatrixType &dNu = *fe_u->getGradBasisFunction();
        setNu_Matrix_byComponent(dim, nnodes_u, Nu, Nuvw);
        setB_Matrix_byComponent(dim, nnodes_u, dNu, B);
        LocalMatrixType &Np = *fe_p->getBasisFunction();

        // K_uu += B^T * D * B
        Kuu.noalias() += fac_u * B.transpose() * De * B;

        // C_up += B^T * m * Np
        Cup.noalias() += fac_u * B.transpose() * m * Np;

        // Fu += N^T * b
        if (hasGravity) {
            Fu.noalias() += fac_u * Nuvw.transpose() * body_force;
        }
    }
    Fu.noalias() += (theta - 1) * Kuu * u0 + (1-theta)* Cup * p0;

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

    // Backward euler
    Fp = (1.0/dt * Mpp - (1-theta)*Kpp)* p0 + 1.0/dt * Cpu * u0;

    // r = K*u - RHS
    LocalVectorType r_u = Kuu * u1 - Cup * p1 - Fu;
    LocalVectorType r_p = 1.0/dt * Cpu * u1 + (1.0/dt * Mpp + theta * Kpp) * p1 - Fp;
//    if (e.getID()==0) {
//        std::cout << "u1=" << std::endl << u1 << std::endl;
//        std::cout << "p1=" << std::endl << p1 << std::endl;
//        std::cout << "Fp=" << std::endl << Fp << std::endl;
//        std::cout << "r_p=" << std::endl << r_p << std::endl;
//    }

    //
    for (size_t i=0; i<dim; i++) {
        vec_r[i] = r_u.segment(i*nnodes_u, nnodes_u);
    }
    vec_r[id_p] = r_p;
}
