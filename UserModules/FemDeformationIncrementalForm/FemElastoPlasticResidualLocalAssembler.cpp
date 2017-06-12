/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemLinearElasticResidualLocalAssembler.cpp
 *
 * Created on 2012-11-29 by Norihiro Watanabe
 */

#include "FemElastoPlasticResidualLocalAssembler.h"

#include "FemLib/Core/Element/IFemElement.h"
#include "MaterialLib/Solid.h"
#include "PhysicsLib/FemLinearElasticTools.h"
#include "PhysicsLib/SmallDeformationMedia.h"
#include "Ogs6FemData.h"


void FemElastoPlasticResidualLocalAssembler::assembly(
        const NumLib::TimeStep &/*time*/, const MeshLib::IElement &e,
        const DiscreteLib::DofEquationIdTable &/*localDofManager*/,
        const MathLib::LocalVector &local_du_n1,
        const MathLib::LocalVector &/*local_du_n*/,
        MathLib::LocalVector &local_r)
{
    //---------------------------------------------------------------
    // configure element
    //---------------------------------------------------------------
    FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);
    const size_t mat_id = e.getGroupID();
    const size_t dim = e.getDimension();
    const size_t nnodes = e.getNumberOfNodes();
    const NumLib::TXPosition e_pos(NumLib::TXPosition::Element, e.getID());

    //---------------------------------------------------------------
    // configure material
    //---------------------------------------------------------------
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    MaterialLib::Solid *solidphase = femData->list_solid[mat_id];
    const size_t n_strain_components = getNumberOfStrainComponents(dim);
    const bool hasGravity = _problem_coordinates.hasZ();

    //---------------------------------------------------------------
    // compute elastic matrix (assuming it is invariant with gauss points)
    //---------------------------------------------------------------
    double nv, E;
    solidphase->poisson_ratio->eval(e_pos, nv);
    solidphase->Youngs_modulus->eval(e_pos, E);

    //---------------------------------------------------------------
    // compute at each integration point
    //---------------------------------------------------------------
    MathLib::LocalVector localInternalForce = MathLib::LocalVector::Zero(local_r.rows());
    MathLib::LocalMatrix matB = MathLib::LocalMatrix::Zero(n_strain_components,
            nnodes * dim);
    MathLib::LocalMatrix matN = MathLib::LocalMatrix::Zero(dim, nnodes * dim);
    FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
    double gp_x[3], real_x[3];
    for (size_t j = 0; j < q->getNumberOfSamplingPoints(); j++) {
        //---------------------------------------------------------------
        // get integration point in reference coordinates
        //---------------------------------------------------------------
        q->getSamplingPoint(j, gp_x);

        //---------------------------------------------------------------
        // compute shape functions
        //---------------------------------------------------------------
        fe->computeBasisFunctions(gp_x);
        MathLib::LocalMatrix &N = *fe->getBasisFunction();
        MathLib::LocalMatrix &dN = *fe->getGradBasisFunction();
        const double fac = fe->getDetJ() * q->getWeight(j);

        //---------------------------------------------------------------
        // compute physical coordinates of this integration point
        //---------------------------------------------------------------
        fe->getRealCoordinates(real_x);
        NumLib::TXPosition gp_pos(NumLib::TXPosition::IntegrationPoint, e.getID(), j, real_x);

        //---------------------------------------------------------------
        // compute Nu, B
        //---------------------------------------------------------------
        setNu_Matrix_byPoint(dim, nnodes, N, matN);
        setB_Matrix_byPoint(dim, nnodes, dN, matB);

        //---------------------------------------------------------------
        // get total stress
        //---------------------------------------------------------------
        // previous stress: S_n
        MathLib::LocalMatrix stress_n = MathLib::LocalMatrix::Zero(n_strain_components, 1);
        _previous_stress->eval(gp_pos, stress_n);
        // incremental stress: dS_n+1 = D B du_n+1
        PhysicsLib::SmallDeformationMedia sdMedia(dim, nv, E);
        sdMedia.setInitialStress(stress_n);
        const MathLib::LocalVector dStrain = matB * local_du_n1;
        sdMedia.incrementStrain(dStrain);
        // total stress: S_n+1 = dS + S_n
        const MathLib::LocalVector stress_n1 = sdMedia.getTotalStress();

        //---------------------------------------------------------------
        // get body force
        //---------------------------------------------------------------
        MathLib::LocalVector body_force = MathLib::LocalVector::Zero(dim);
        if (hasGravity) {
            double rho_s = .0;
            solidphase->density->eval(e_pos, rho_s);
            body_force[dim - 1] = - rho_s * 9.81;
        }

        //---------------------------------------------------------------
        // compute residual: r = int {B^T S_{n+1} + N^T b} - f_ex
        // * the external force term (f_ex) is set by source term
        // * the internal force terms (f_in) is computed here
        //---------------------------------------------------------------
        // f_in += B^T * S +N^T b
        localInternalForce.noalias() += fac * matB.transpose() * stress_n1;
        if (hasGravity) {
            localInternalForce.noalias() += fac * matN.transpose() * body_force;
        }
    }

    local_r.noalias() = localInternalForce;
}


