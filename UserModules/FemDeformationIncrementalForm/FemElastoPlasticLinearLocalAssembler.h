/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemElastoPlasticLinearLocalAssembler.h
 *
 * Created on 2012-11-29 by Norihiro Watanabe
 */

#pragma once

#include "MeshLib/Core/IElement.h"
#include "FemLib/Tools/IFeObjectContainer.h"
#include "NumLib/TransientAssembler/IElementWiseTransientLinearEQSLocalAssembler.h"
#include "NumLib/TimeStepping/TimeStep.h"

/**
 *
 */
class FemElastoPlasticLinearLocalAssembler:
        public NumLib::IElementWiseTransientLinearEQSLocalAssembler
{
public:

    /**
     *
     * @param feObjects
     */
    explicit FemElastoPlasticLinearLocalAssembler(
            FemLib::IFeObjectContainer* feObjects)
            : _feObjects(feObjects->clone())
    {
    }

    /**
     * Copy constructor
     * @param src
     */
    FemElastoPlasticLinearLocalAssembler(const FemElastoPlasticLinearLocalAssembler &src)
    : _feObjects(src._feObjects->clone())
    {
    }

    /**
     *
     */
    virtual ~FemElastoPlasticLinearLocalAssembler()
    {
        BaseLib::releaseObject(_feObjects);
    }

    /**
     *
     * @param time
     * @param e
     * @param localDofManager
     * @param local_u_n1
     * @param local_u_n
     * @param eqs
     */
    virtual void assembly(const NumLib::TimeStep &time,
            const MeshLib::IElement &e,
            const DiscreteLib::DofEquationIdTable &localDofManager,
            const MathLib::LocalVector &local_u_n1,
            const MathLib::LocalVector &local_u_n, MathLib::LocalEquation &eqs) OGS_DECL_OVERRIDE;

private:
    FemLib::IFeObjectContainer* _feObjects;
};

