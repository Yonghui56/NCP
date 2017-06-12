/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file PressureBasedGWTimeODELocalAssembler.h
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#pragma once

#include "FemLib/Core/Element/IFemElement.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Fluid.h"
#include "Ogs6FemData.h"

/**
 * \brief Local assembly class for time-ODE GW equation
 */
template <class T>
class PressureBasedGWTimeODELocalAssembler: public T
{
public:
    typedef MathLib::LocalVector LocalVector;
    typedef MathLib::LocalMatrix LocalMatrix;

    explicit PressureBasedGWTimeODELocalAssembler(FemLib::LagrangeFeObjectContainer &feObjects, const MeshLib::CoordinateSystem &problem_coordinates)
    : _feObjects(feObjects), _problem_coordinates(problem_coordinates)
    {
    };

    virtual ~PressureBasedGWTimeODELocalAssembler() {};

protected:
    virtual void assembleODE(const NumLib::TimeStep &/*time*/, const MeshLib::IElement &e, const LocalVector &/*u1*/, const LocalVector &/*u0*/, LocalMatrix &localM, LocalMatrix &localK, LocalVector &localF)
    {
        FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
        const NumLib::TXPosition e_pos(NumLib::TXPosition::Element, e.getID());
        const MathLib::LocalMatrix* matR = NULL;
        if (e.getDimension() < _problem_coordinates.getDimension()) {
            MeshLib::ElementCoordinatesMappingLocal* ele_local_coord;
            ele_local_coord = (MeshLib::ElementCoordinatesMappingLocal*)e.getMappedCoordinates();
            matR = &ele_local_coord->getRotationMatrixToOriginal();
        }

        size_t mat_id = e.getGroupID();
        MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];
        MaterialLib::Fluid* fluid = Ogs6FemData::getInstance()->list_fluid[0];
        double geo_area = 1.0;
        pm->geo_area->eval(e_pos, geo_area);
        double mu = .0;
        fluid->dynamic_viscosity->eval(e_pos, mu);
        double rho_f = .0;
        LocalVector vec_g;
        const bool hasGravityEffect = _problem_coordinates.hasZ();
        if (hasGravityEffect) {
            fluid->density->eval(e_pos, rho_f);
            vec_g = LocalVector::Zero(_problem_coordinates.getDimension());
            vec_g[_problem_coordinates.getIndexOfZ()] = -9.81;
        }

        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
            q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            LocalMatrix &Np = *fe->getBasisFunction();
            LocalMatrix &dNp = *fe->getGradBasisFunction();
            fe->getRealCoordinates(real_x);
            double fac = geo_area * fe->getDetJ() * q->getWeight(j);

            double k;
            pm->permeability->eval(real_x, k);
            double s;
            pm->storage->eval(real_x, s);
            double k_mu;
            k_mu = k / mu;
            MathLib::LocalMatrix local_k_mu = MathLib::LocalMatrix::Identity(e.getDimension(), e.getDimension());
            local_k_mu *= k_mu;
            MathLib::LocalMatrix global_k_mu;
            if (e.getDimension() < _problem_coordinates.getDimension()) {
                MathLib::LocalMatrix local2 = MathLib::LocalMatrix::Zero(_problem_coordinates.getDimension(), _problem_coordinates.getDimension());
                local2.block(0, 0, local_k_mu.rows(), local_k_mu.cols()) = local_k_mu.block(0, 0, local_k_mu.rows(), local_k_mu.cols());
                //local2.topLeftCorner(local_k_mu.rows(), local_k_mu.cols()) = local_k_mu;
                global_k_mu = (*matR) * local2 * matR->transpose();
            } else {
                global_k_mu = local_k_mu;
            }

            // M_pp += Np^T * S * Np
            localM.noalias() += fac * Np.transpose() * s * Np;

            // K_pp += dNp^T * K * dNp
            localK.noalias() += fac * dNp.transpose() * global_k_mu * dNp;

            if (hasGravityEffect) {
                // F += dNp^T * K * rho * gz
                localF.noalias() += fac * dNp.transpose() * global_k_mu * rho_f * vec_g;
            }
        }
    }

private:
    FemLib::LagrangeFeObjectContainer _feObjects;
    MeshLib::CoordinateSystem _problem_coordinates;
};

