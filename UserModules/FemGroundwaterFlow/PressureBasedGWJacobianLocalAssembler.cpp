/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file GWJacobianLocalAssembler.cpp
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#include "PressureBasedGWJacobianLocalAssembler.h"

#include "NumLib/Function/TXFunction.h"
#include "Ogs6FemData.h"

void PressureBasedGWJacobianLocalAssembler::assembly(const NumLib::TimeStep &ts, const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &/*localDofManager*/, const MathLib::LocalVector &/*u1*/, const MathLib::LocalVector &/*u0*/,  MathLib::LocalMatrix &localJ)
{
    FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
    const NumLib::TXPosition e_pos(NumLib::TXPosition::Element, e.getID());
    size_t mat_id = e.getGroupID();
    MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];
    MaterialLib::Fluid* fluid = Ogs6FemData::getInstance()->list_fluid[0];
    double mu = .0;
    fluid->dynamic_viscosity->eval(e_pos, mu);
    double rho_f = .0;
    fluid->density->eval(e_pos, rho_f);

    MathLib::LocalMatrix localM = MathLib::LocalMatrix::Zero(localJ.rows(), localJ.cols());
    MathLib::LocalMatrix localK = MathLib::LocalMatrix::Zero(localJ.rows(), localJ.cols());

    FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
    double gp_x[3], real_x[3];
    for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
        q->getSamplingPoint(j, gp_x);
        fe->computeBasisFunctions(gp_x);
        fe->getRealCoordinates(real_x);
        MathLib::LocalMatrix &Np = *fe->getBasisFunction();
        MathLib::LocalMatrix &dNp = *fe->getGradBasisFunction();
        fe->getRealCoordinates(real_x);
        double fac = fe->getDetJ() * q->getWeight(j);

        double k;
        pm->hydraulic_conductivity->eval(real_x, k);
        double s;
        pm->storage->eval(real_x, s);
        double k_mu;
        k_mu = k / mu;

        // M_pp += Np^T * S * Np
        localM.noalias() += fac * Np.transpose() * s * Np;

        // K_pp += dNp^T * K * dNp
        localK.noalias() += fac * dNp.transpose() * k_mu * dNp;
    }

    double euler_theta = 1.0; //TODO where to get this parameter?
    localJ.noalias() = 1./ts.getTimeStepSize() * localM + euler_theta *localK;
}

