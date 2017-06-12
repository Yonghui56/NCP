/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file HeatTransportJacobianLocalAssembler.h
 *
 * Created on 2012-10-23 by Norihiro Watanabe
 */

#pragma once

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Core/Integration/Integration.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "MaterialLib/PorousMedia.h"

#include "Ogs6FemData.h"

class HeatTransportJacobianLocalAssembler: public NumLib::IElementWiseTransientJacobianLocalAssembler
{
public:
    HeatTransportJacobianLocalAssembler(FemLib::LagrangeFeObjectContainer* feObjects)
        : _feObjects(*feObjects), _vel(NULL)
    {
    };

    virtual ~HeatTransportJacobianLocalAssembler() {};

    void setVelocity(const NumLib::ITXFunction *vel)
    {
        _vel = const_cast<NumLib::ITXFunction*>(vel);
    }

    void assembly(const NumLib::TimeStep &time, const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &/*localDofManager*/, const MathLib::LocalVector &/*u1*/, const MathLib::LocalVector &/*u0*/, MathLib::LocalMatrix &localJ)
    {
        FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
        const NumLib::TXPosition e_pos(NumLib::TXPosition::Element, e.getID());
        const size_t n_dim = e.getDimension();
        size_t mat_id = e.getGroupID();
        MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];
        MaterialLib::Solid* solid = Ogs6FemData::getInstance()->list_solid[mat_id];
        MaterialLib::Fluid* fluid = Ogs6FemData::getInstance()->list_fluid[0];

        MathLib::LocalMatrix localM = MathLib::LocalMatrix::Zero(localJ.rows(), localJ.cols());
        MathLib::LocalMatrix localDispersion = MathLib::LocalMatrix::Zero(localJ.rows(), localJ.cols());
        MathLib::LocalMatrix localAdvection = MathLib::LocalMatrix::Zero(localJ.rows(), localJ.cols());

        double n = .0;
        pm->porosity->eval(e_pos, n);

        double rho_f = .0;
        double cp_f = .0;
        double lambda_f = .0;
        fluid->density->eval(e_pos, rho_f);
        fluid->specific_heat->eval(e_pos, cp_f);
        fluid->thermal_conductivity->eval(e_pos, lambda_f);
        double rho_s = .0;
        double cp_s = .0;
        double lambda_s = .0;
        solid->density->eval(e_pos, rho_s);
        solid->specific_heat->eval(e_pos, cp_s);
        solid->thermal_conductivity->eval(e_pos, lambda_s);
        double rho_cp_eff = n * rho_f * cp_f + (1.-n) * rho_s * cp_s;
        double lambda_eff = n * lambda_f + (1.-n) * lambda_s;

        MathLib::LocalMatrix rho_cp(1,1);
        rho_cp(0,0) = rho_cp_eff;
        MathLib::LocalMatrix lambda(1,1);
        lambda(0,0) = lambda_eff;

        NumLib::ITXFunction::DataType v;
        _vel->eval(e_pos, v);
        NumLib::ITXFunction::DataType v2 = v.topRows(n_dim).transpose();
        v2 *= rho_cp_eff;


        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
            q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);

            fe->integrateWxN(j, rho_cp, localM);
            fe->integrateDWxDN(j, lambda, localDispersion);
            fe->integrateWxDN(j, v2, localAdvection);
        }

        double dt = time.getTimeStepSize();
        double theta = 1.0;
        localJ = (1.0/dt * localM + theta *(localDispersion + localAdvection));

        //std::cout << "M="; localM.write(std::cout); std::cout << std::endl;
        //std::cout << "L="; matDiff.write(std::cout); std::cout << std::endl;
        //std::cout << "A="; matAdv.write(std::cout); std::cout << std::endl;
    }

private:
    FemLib::LagrangeFeObjectContainer _feObjects;
    NumLib::ITXFunction* _vel;
};
