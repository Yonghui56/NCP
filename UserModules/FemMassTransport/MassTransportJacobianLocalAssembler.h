/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file MassTransportJacobianLocalAssembler.h
 *
 * Created on 2012-07-13 by Norihiro Watanabe
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
#include "MaterialLib/Compound.h"

#include "Ogs6FemData.h"

class MassTransportJacobianLocalAssembler: public NumLib::IElementWiseTransientJacobianLocalAssembler
{
public:
    MassTransportJacobianLocalAssembler(MaterialLib::Compound* cmp, FemLib::LagrangeFeObjectContainer* feObjects)
        : _cmp(cmp), _feObjects(*feObjects), _vel(NULL)
    {
    };

    virtual ~MassTransportJacobianLocalAssembler() {};

    void setVelocity(const NumLib::ITXFunction *vel)
    {
        _vel = const_cast<NumLib::ITXFunction*>(vel);
    }

    void assembly(const NumLib::TimeStep &time, const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &/*localDofManager*/, const MathLib::LocalVector &/*u1*/, const MathLib::LocalVector &/*u0*/, MathLib::LocalMatrix &localJ)
    {
        FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
        const size_t n_dim = e.getDimension();
        size_t mat_id = e.getGroupID();
        MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];
        double cmp_mol_diffusion = .0;
        _cmp->molecular_diffusion->eval(0, cmp_mol_diffusion);

        MathLib::LocalMatrix matM = MathLib::LocalMatrix::Zero(localJ.rows(), localJ.cols());
        MathLib::LocalMatrix matDiff = MathLib::LocalMatrix::Zero(localJ.rows(), localJ.cols());
        MathLib::LocalMatrix matAdv = MathLib::LocalMatrix::Zero(localJ.rows(), localJ.cols());

        MathLib::LocalMatrix poro(1,1);
        MathLib::LocalMatrix d_poro(1,1);
        NumLib::ITXFunction::DataType v;

        double gp_x[3], real_x[3];
        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
            q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);
            NumLib::TXPosition gp_pos(NumLib::TXPosition::IntegrationPoint, e.getID(), j, real_x);

            pm->porosity->eval(gp_pos, poro);
            d_poro(0,0) = cmp_mol_diffusion * poro(0,0);
            _vel->eval(gp_pos, v);
            NumLib::ITXFunction::DataType v2 = v.topRows(n_dim).transpose();

            fe->integrateWxN(j, poro, matM);
            fe->integrateDWxDN(j, d_poro, matDiff);
            fe->integrateWxDN(j, v2, matAdv);
        }

        double dt = time.getTimeStepSize();
        double theta = 1.0;
        localJ = (1.0/dt * matM + theta *(matDiff + matAdv));

        //std::cout << "M="; localM.write(std::cout); std::cout << std::endl;
        //std::cout << "L="; matDiff.write(std::cout); std::cout << std::endl;
        //std::cout << "A="; matAdv.write(std::cout); std::cout << std::endl;
    }

private:
    MaterialLib::Compound* _cmp;
    FemLib::LagrangeFeObjectContainer _feObjects;
    NumLib::ITXFunction* _vel;
};
