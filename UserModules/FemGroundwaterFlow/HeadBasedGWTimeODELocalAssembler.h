/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file GWTimeODELocalAssembler.h
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#pragma once

#include "FemLib/Core/Element/IFemElement.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "MaterialLib/PorousMedia.h"
#include "Ogs6FemData.h"

/**
 * \brief Local assembly class for time-ODE GW equation
 */
template <class T>
class HeadBasedGWTimeODELocalAssembler: public T
{
public:
    typedef MathLib::LocalVector LocalVector;
    typedef MathLib::LocalMatrix LocalMatrix;

    explicit HeadBasedGWTimeODELocalAssembler(
                FemLib::LagrangeFeObjectContainer &feObjects)
    : _feObjects(feObjects)
    {};

    virtual ~HeadBasedGWTimeODELocalAssembler() {};

protected:
    virtual void assembleODE(   
                const NumLib::TimeStep &/*time*/, const MeshLib::IElement &e, 
                const LocalVector &/*u1*/, const LocalVector &/*u0*/, 
                LocalMatrix &localM, LocalMatrix &localK, LocalVector &/*localF*/
                )
    {
        FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);

        size_t mat_id = e.getGroupID();
        MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];

        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
            q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);

            LocalMatrix k;
            pm->hydraulic_conductivity->eval(real_x, k);
            LocalMatrix s;
            pm->storage->eval(real_x, s);

            fe->integrateWxN(j, s, localM);
            fe->integrateDWxDN(j, k, localK);
        }
    }

private:
    FemLib::LagrangeFeObjectContainer _feObjects;
};

