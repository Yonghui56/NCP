/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FlowInPorousMediaTimeODELocalAssembler.h
 *
 * Created on 2012-10-24 by Norihiro Watanabe
 */

#pragma once

#include "FemLib/Core/Element/IFemElement.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Fluid.h"
#include "Ogs6FemData.h"

/**
 * \brief Local assembly of time-ODE for liquid flow in porous media
 *
 * \tparam T_TIME_DIS   Time discretization scheme
 */
template <class T_TIME_DIS>
class FlowInPorousMediaTimeODELocalAssembler: public T_TIME_DIS
{
public:
    /**
     *
     * @param feObjects
     * @param problem_coordinates
     */
    FlowInPorousMediaTimeODELocalAssembler (
            const FemLib::IFeObjectContainer* feObjects,
            const MeshLib::CoordinateSystem &problem_coordinates
            )
    : _feObjects(feObjects->clone()), _problem_coordinates(problem_coordinates)
    {
    };

    /**
     *
     */
    virtual ~FlowInPorousMediaTimeODELocalAssembler()
    {
        BaseLib::releaseObject(_feObjects);
    };

protected:
    /**
     * assemble components of local time ODE
     *
     * @param time
     * @param e
     * @param u1
     * @param u0
     * @param localM
     * @param localK
     * @param localF
     */
    virtual void assembleODE(const NumLib::TimeStep &/*time*/, const MeshLib::IElement &e, const MathLib::LocalVector &/*u1*/, const MathLib::LocalVector &/*u0*/, MathLib::LocalMatrix &localM, MathLib::LocalMatrix &localK, MathLib::LocalVector &localF);

private:
    FemLib::IFeObjectContainer* _feObjects;
    const MeshLib::CoordinateSystem _problem_coordinates;
};


template <class T>
void FlowInPorousMediaTimeODELocalAssembler<T>::assembleODE(const NumLib::TimeStep &/*time*/, const MeshLib::IElement &e, const MathLib::LocalVector &/*u1*/, const MathLib::LocalVector &/*u0*/, MathLib::LocalMatrix &localM, MathLib::LocalMatrix &localK, MathLib::LocalVector &localF)
{
    // element information
    const size_t mat_id = e.getGroupID();

    // material data
    MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];
    MaterialLib::Fluid* fluid = Ogs6FemData::getInstance()->list_fluid[0];

    // gravity effect
    const bool hasGravityEffect = _problem_coordinates.hasZ();
    MathLib::LocalVector vec_g;
    if (hasGravityEffect) {
        vec_g = MathLib::LocalVector::Zero(_problem_coordinates.getDimension());
        vec_g[_problem_coordinates.getIndexOfZ()] = -9.81;
    }

    // numerical integration
    FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);
    FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
    double gp_x[3], real_x[3];
    for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
        // compute shape functions and parameters for integration
        q->getSamplingPoint(j, gp_x);
        fe->computeBasisFunctions(gp_x);
        MathLib::LocalMatrix &Np = *fe->getBasisFunction();
        MathLib::LocalMatrix &dNp = *fe->getGradBasisFunction();
        const double fac = fe->getDetJ() * q->getWeight(j);

        // evaluate material properties
        fe->getRealCoordinates(real_x);
        NumLib::TXPosition gp_pos(NumLib::TXPosition::IntegrationPoint, e.getID(), j, real_x);
        double mu = .0;
        double rho_f = .0;
        fluid->dynamic_viscosity->eval(gp_pos, mu);
        if (hasGravityEffect) {
            fluid->density->eval(gp_pos, rho_f);
        }
        double k = .0;
        double s = .0;
        pm->permeability->eval(gp_pos, k);
        pm->storage->eval(gp_pos, s);
        MathLib::LocalMatrix mat_k_mu = MathLib::LocalMatrix::Identity(e.getDimension(), e.getDimension());
        mat_k_mu *= k / mu;

        // M_pp += Np^T * S * Np
        localM.noalias() += fac * Np.transpose() * s * Np;

        // K_pp += dNp^T * K * dNp
        localK.noalias() += fac * dNp.transpose() * mat_k_mu * dNp;

        if (hasGravityEffect) {
            // F += dNp^T * K * rho * gz
            localF.noalias() += fac * dNp.transpose() * mat_k_mu * rho_f * vec_g;
        }
    }
}

