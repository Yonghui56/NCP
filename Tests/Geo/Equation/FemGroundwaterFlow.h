/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemGroundwaterFlow.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "FemLib/Core/Element/IFemElement.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "Tests/Geo/Material/PorousMedia.h"

namespace Geo
{

template <class T>
class GroundwaterFlowTimeODELocalAssembler: public T
{
private:
    PorousMedia* _pm;
    FemLib::LagrangeFeObjectContainer* _feObjects;
public:
    typedef MathLib::LocalVector LocalVector;
    typedef MathLib::LocalMatrix LocalMatrix;

    GroundwaterFlowTimeODELocalAssembler(FemLib::LagrangeFeObjectContainer &feObjects, PorousMedia &pm)
    : _pm(&pm), _feObjects(&feObjects)
    {
    };

    virtual ~GroundwaterFlowTimeODELocalAssembler() {};

protected:
    virtual void assembleODE(const NumLib::TimeStep &/*time*/, const MeshLib::IElement &e, const LocalVector &/*u1*/, const LocalVector &/*u0*/, LocalMatrix &/*localM*/, LocalMatrix &localK, LocalVector &/*localF*/)
    {
        FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);

        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
            q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);

            LocalMatrix k;
            _pm->hydraulic_conductivity->eval(real_x, k);

            //fe->integrateWxN(_pm->storage, localM);
            fe->integrateDWxDN(j, k, localK);
        }
    }
};

class GroundwaterFlowJacobianLocalAssembler: public NumLib::IElementWiseTransientJacobianLocalAssembler
{
private:
    PorousMedia* _pm;
    FemLib::LagrangeFeObjectContainer* _feObjects;
public:
    GroundwaterFlowJacobianLocalAssembler(FemLib::LagrangeFeObjectContainer &feObjects, PorousMedia &pm)
    : _pm(&pm), _feObjects(&feObjects)
    {
    };

    void assembly(const NumLib::TimeStep &/*time*/, const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &/*localDofManager*/, const MathLib::LocalVector &/*u1*/, const MathLib::LocalVector &/*u0*/,  MathLib::LocalMatrix &localJ)
    {
        FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);

        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
            q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);

            MathLib::LocalMatrix k;
            _pm->hydraulic_conductivity->eval(real_x, k);

            //fe->integrateWxN(_pm->storage, localM);
            fe->integrateDWxDN(j, k, localJ);
        }
    }
};


} //end
