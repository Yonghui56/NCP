/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FlowInFractureTimeODELocalAssembler.h
 *
 * Created on 2012-10-24 by Norihiro Watanabe
 */

#pragma once

#include "FemLib/Core/Element/IFemElement.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Fluid.h"
#include "Ogs6FemData.h"

/**
 * \brief Local assembly class for time-ODE GW equation in fractures
 */
template <class T>
class FlowInFractureTimeODELocalAssembler: public T
{
public:
    /**
     *
     * @param feObjects
     * @param problem_coordinates
     */
    FlowInFractureTimeODELocalAssembler(const FemLib::IFeObjectContainer* feObjects, const MeshLib::CoordinateSystem &problem_coordinates)
    : _feObjects(feObjects->clone()), _problem_coordinates(problem_coordinates)
    {
    };

    /**
     *
     */
    virtual ~FlowInFractureTimeODELocalAssembler()
    {
        BaseLib::releaseObject(_feObjects);
    };

protected:
    /**
     *
     * @param
     * @param e
     * @param
     * @param
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
void FlowInFractureTimeODELocalAssembler<T>::assembleODE(const NumLib::TimeStep &/*time*/, const MeshLib::IElement &e, const MathLib::LocalVector &/*u1*/, const MathLib::LocalVector &/*u0*/, MathLib::LocalMatrix &localM, MathLib::LocalMatrix &localK, MathLib::LocalVector &localF)
{
    // element information
    const size_t mat_id = e.getGroupID();
    MeshLib::ElementCoordinatesMappingLocal* ele_local_coord;
    ele_local_coord = (MeshLib::ElementCoordinatesMappingLocal*)e.getMappedCoordinates();
    assert(ele_local_coord!=NULL);
    const MathLib::LocalMatrix* matR = &ele_local_coord->getRotationMatrixToOriginal();

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

    // numerical integrations
    FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);
    FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
    double gp_x[3], real_x[3];
    for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
        // compute shape functions and parameters for integration
        q->getSamplingPoint(j, gp_x);
        fe->computeBasisFunctions(gp_x);
        MathLib::LocalMatrix &Np = *fe->getBasisFunction();
        MathLib::LocalMatrix &dNp = *fe->getGradBasisFunction();
        fe->getRealCoordinates(real_x);
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
        double b = 1.0;
        double k = .0;
        double s = .0;
        pm->geo_area->eval(gp_pos, b);
        pm->permeability->eval(gp_pos, k);
        pm->storage->eval(gp_pos, s);
        MathLib::LocalMatrix local_k_mu = MathLib::LocalMatrix::Identity(e.getDimension(), e.getDimension());
        local_k_mu *= k / mu;
        MathLib::LocalMatrix local2 = MathLib::LocalMatrix::Zero(_problem_coordinates.getDimension(), _problem_coordinates.getDimension());
        local2.block(0, 0, local_k_mu.rows(), local_k_mu.cols()) = local_k_mu.block(0, 0, local_k_mu.rows(), local_k_mu.cols());
        //local2.topLeftCorner(local_k_mu.rows(), local_k_mu.cols()) = local_k_mu;
        MathLib::LocalMatrix global_k_mu = (*matR) * local2 * matR->transpose();

        // M_pp += Np^T * b * S * Np
        localM.noalias() += fac * Np.transpose() * b * s * Np;

        // K_pp += dNp^T * b * K * dNp
        localK.noalias() += fac * dNp.transpose() * b * global_k_mu * dNp;

        if (hasGravityEffect) {
            // F += dNp^T * b * K * rho * gz
            localF.noalias() += fac * dNp.transpose() * b * global_k_mu * rho_f * vec_g;
        }
    }
}
