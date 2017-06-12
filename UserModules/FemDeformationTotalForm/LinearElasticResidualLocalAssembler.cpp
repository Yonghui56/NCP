/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LinearElasticResidualLocalAssembler.cpp
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */


#include "LinearElasticResidualLocalAssembler.h"

#include "FemLib/Core/Element/IFemElement.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Solid.h"
#include "PhysicsLib/FemLinearElasticTools.h"
#include "Ogs6FemData.h"

/// assemble a local residual for the given element
void FemLinearElasticResidualLocalAssembler::assembly
    (   const NumLib::TimeStep &/*time*/,
        const MeshLib::IElement &e,
        const DiscreteLib::DofEquationIdTable &/*localDofManager*/,
        const MathLib::LocalVector &local_u_n1,
        const MathLib::LocalVector &/*local_u_n*/,
        MathLib::LocalVector &local_r    )
{

    FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
    size_t mat_id = e.getGroupID();
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    //MaterialLib::PorousMedia* pm = femData->list_pm[mat_id];
    MaterialLib::Solid *solidphase = femData->list_solid[mat_id];

    const size_t dim = e.getDimension();
    const size_t n_strain_components = (dim==2 ? 4 : 6);
    const size_t nnodes = e.getNumberOfNodes();
    const NumLib::TXPosition e_pos(NumLib::TXPosition::Element, e.getID());

    // set D
    MathLib::LocalMatrix matD = MathLib::LocalMatrix::Zero(n_strain_components, n_strain_components);
    MathLib::LocalMatrix nv(1,1);
    MathLib::LocalMatrix E(1,1);
    solidphase->poisson_ratio->eval(e_pos, nv);
    solidphase->Youngs_modulus->eval(e_pos, E);
    double Lambda, G, K;
    MaterialLib::calculateLameConstant(nv(0,0), E(0,0), Lambda, G, K);
    MaterialLib::setElasticConsitutiveTensor(dim, Lambda, G, matD);

    // body force
    MathLib::LocalVector body_force = MathLib::LocalVector::Zero(dim);
    bool hasGravity = false;
    if (hasGravity) {
        double rho_s = .0;
        solidphase->density->eval(e_pos, rho_s);
        body_force[dim-1] = rho_s * 9.81;
    }

    //
    MathLib::LocalMatrix localK = MathLib::LocalMatrix::Zero(local_r.rows(), local_r.rows());
    MathLib::LocalVector localRHS = MathLib::LocalVector::Zero(local_r.rows());
    MathLib::LocalMatrix matB = MathLib::LocalVector::Zero(n_strain_components, nnodes*dim);
    MathLib::LocalMatrix matN = MathLib::LocalVector::Zero(dim, nnodes*dim);
    FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
    double gp_x[3], real_x[3];
    for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
        q->getSamplingPoint(j, gp_x);
        fe->computeBasisFunctions(gp_x);
        fe->getRealCoordinates(real_x);
        double fac = fe->getDetJ() * q->getWeight(j);

        // set N,B
        MathLib::LocalMatrix &N = *fe->getBasisFunction();
        MathLib::LocalMatrix &dN = *fe->getGradBasisFunction();
        setNu_Matrix_byPoint(dim, nnodes, N, matN);
        setB_Matrix_byPoint(dim, nnodes, dN, matB);

        // K += B^T * D * B
        localK.noalias() += fac * matB.transpose() * matD * matB;

        // RHS
        localRHS.noalias() += fac * matN.transpose() * body_force;

    }

    // r = K*u - RHS
    local_r.noalias() = localK * local_u_n1 - localRHS;
}
