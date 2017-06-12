/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemElastoPlasticResidualLocalAssembler.h
 *
 * Created on 2012-11-29 by Norihiro Watanabe
 */

#pragma once

#include "MeshLib/Core/CoordinateSystem.h"
#include "FemLib/Tools/IFeObjectContainer.h"
#include "NumLib/Function/ITXFunction.h"
#include "NumLib/TransientAssembler/IElementWiseTransientResidualLocalAssembler.h"

/**
 * \brief Local residual vector assembler
 *
 * The residual \f$ r$ \f$ given as
 * \f[
 *  r = B^T d\sigma_{n+1} + B^T \sigma_n + b - f
 *    = (B^T D_T B) du_{n+1} + B^T S_n + b - f
 * \f]
 * where B is strain divergence matrix, \f$\sigma\f$ is stress, \f$b\f$ is body force,
 * and \f$f\f$ is external force. \f$n+1, n\f$ mean current time step and previous time step.
 */
class FemElastoPlasticResidualLocalAssembler:
        public NumLib::IElementWiseTransientResidualLocalAssembler
{
public:
    /**
     *
     * @param feObjects
     * @param problem_coordinates
     */
    FemElastoPlasticResidualLocalAssembler(
            const FemLib::IFeObjectContainer* feObjects,
            const MeshLib::CoordinateSystem &problem_coordinates)
            : _feObjects(feObjects->clone()), _previous_stress(nullptr), _problem_coordinates(
                    problem_coordinates)
    {
    }

    /**
     *
     */
    virtual ~FemElastoPlasticResidualLocalAssembler()
    {
        BaseLib::releaseObject(_feObjects);
    }

    /**
     *
     * @param s
     */
    void setPreviousStress(const NumLib::ITXFunction* s)
    {
        _previous_stress = s;
    }

    /**
     * assemble a local residual for the given element
     * @param time
     * @param e
     * @param localDofManager
     * @param local_u_n1
     * @param local_u_n
     * @param local_r
     */
    virtual void assembly(const NumLib::TimeStep &time,
            const MeshLib::IElement &e,
            const DiscreteLib::DofEquationIdTable &localDofManager,
            const MathLib::LocalVector &local_u_n1,
            const MathLib::LocalVector &local_u_n,
            MathLib::LocalVector &local_r);

private:
    FemLib::IFeObjectContainer* _feObjects;
    const NumLib::ITXFunction* _previous_stress;
    const MeshLib::CoordinateSystem _problem_coordinates;
};
