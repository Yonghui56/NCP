/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemElastoPlasticJacobianLocalAssembler.h
 *
 * Created on 2012-11-29 by Norihiro Watanabe
 */

#pragma once

#include "FemLib/Tools/IFeObjectContainer.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"

/**
 * \brief Local Jacobian matrix assembler
 *
 * The residual \f$ r$ \f$ given as
 * \f[
 *  r = B^T d\sigma_{n+1} + B^T \sigma_n + b - f
 *    = (B^T D_T B) du_{n+1} + B^T S_n + b - f
 * \f]
 * where B is strain divergence matrix, \f$\sigma\f$ is stress, \f$b\f$ is body force,
 * and \f$f\f$ is external force. \f$n+1, n\f$ mean current time step and previous time step.
 * Its Jacobian matrix is expressed
 * \f[
 *  J = B^T D_T B
 * \f]
 */
class FemElastoPlasticJacobianLocalAssembler:
        public NumLib::IElementWiseTransientJacobianLocalAssembler
{
public:
    /**
     *
     * @param feObjects
     */
    explicit FemElastoPlasticJacobianLocalAssembler(
            const FemLib::IFeObjectContainer* feObjects)
            : _feObjects(feObjects->clone())
    {
    }

    /**
     *
     */
    virtual ~FemElastoPlasticJacobianLocalAssembler()
    {
        BaseLib::releaseObject(_feObjects);
    }

    /**
     *
     * @param time
     * @param e
     * @param localDofManager
     * @param u1
     * @param u0
     * @param localJ
     */
    void assembly(const NumLib::TimeStep &time, const MeshLib::IElement &e,
            const DiscreteLib::DofEquationIdTable &localDofManager,
            const MathLib::LocalVector &u1,
            const MathLib::LocalVector &u0, MathLib::LocalMatrix &localJ);

private:
    FemLib::IFeObjectContainer* _feObjects;
};
