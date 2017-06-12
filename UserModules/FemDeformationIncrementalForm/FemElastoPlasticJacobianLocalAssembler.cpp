/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemElastoPlasticJacobianLocalAssembler.cpp
 *
 * Created on 2012-11-29 by Norihiro Watanabe
 */

#include "FemElastoPlasticJacobianLocalAssembler.h"

#include "FemLib/Core/Element/IFemElement.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Solid.h"
#include "PhysicsLib/FemLinearElasticTools.h"
#include "Ogs6FemData.h"

void FemElastoPlasticJacobianLocalAssembler::assembly(
        const NumLib::TimeStep &/*time*/, const MeshLib::IElement &e,
        const DiscreteLib::DofEquationIdTable &/*localDofManager*/,
        const MathLib::LocalVector &/*u1*/, const MathLib::LocalVector &/*u0*/,
        MathLib::LocalMatrix &localJ)
{
    FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);
    size_t mat_id = e.getGroupID();
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    //MaterialLib::PorousMedia* pm = femData->list_pm[mat_id];
    MaterialLib::Solid *solidphase = femData->list_solid[mat_id];

    const size_t dim = e.getDimension();
    const size_t n_strain_components = (dim == 2 ? 4 : 6);
    const size_t nnodes = e.getNumberOfNodes();
    const NumLib::TXPosition e_pos(NumLib::TXPosition::Element, e.getID());

    // set tangential D
    Eigen::MatrixXd matDT = Eigen::MatrixXd::Zero(n_strain_components,
            n_strain_components);
    MathLib::LocalMatrix nv(1, 1);
    MathLib::LocalMatrix E(1, 1);
    solidphase->poisson_ratio->eval(e_pos, nv);
    solidphase->Youngs_modulus->eval(e_pos, E);
    double Lambda, G, K;
    MaterialLib::calculateLameConstant(nv(0, 0), E(0, 0), Lambda, G, K);
    MaterialLib::setElasticConsitutiveTensor(dim, Lambda, G, matDT);

    // J = [B^T DT B]
    MathLib::LocalMatrix matB(n_strain_components, nnodes * dim);
    MathLib::LocalMatrix matN(dim, nnodes * dim);
    FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
    double gp_x[3], real_x[3];
    for (size_t j = 0; j < q->getNumberOfSamplingPoints(); j++) {
        q->getSamplingPoint(j, gp_x);
        fe->computeBasisFunctions(gp_x);
        fe->getRealCoordinates(real_x);
        double fac = fe->getDetJ() * q->getWeight(j);

        // set N,B
        MathLib::LocalMatrix &N = *fe->getBasisFunction();
        MathLib::LocalMatrix &dN = *fe->getGradBasisFunction();
        setNu_Matrix_byPoint(dim, nnodes, N, matN);
        setB_Matrix_byPoint(dim, nnodes, dN, matB);

        // J += B^T * D_T * B
        localJ.noalias() += fac * matB.transpose() * matDT * matB;
    }
}
