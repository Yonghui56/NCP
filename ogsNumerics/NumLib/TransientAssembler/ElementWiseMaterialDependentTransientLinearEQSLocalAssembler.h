/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementWiseMaterialDependentTransientLinearEQSLocalAssembler.h
 *
 * Created on 2012-10-24 by Norihiro Watanabe
 */

#pragma once

#include <cassert>

#include "MeshLib/Core/IElement.h"
#include "IElementWiseTransientLinearEQSLocalAssembler.h"


namespace NumLib
{

/**
 * \brief Euler scheme element assembler for time ODE formulations
 *
 * @tparam  T_USER_ASSEMBLY     User-given assembler
 */
class ElementWiseMaterialDependentTransientLinearEQSLocalAssembler : public IElementWiseTransientLinearEQSLocalAssembler
{
public:
    ElementWiseMaterialDependentTransientLinearEQSLocalAssembler() {};

    virtual ~ElementWiseMaterialDependentTransientLinearEQSLocalAssembler()
    {
        BaseLib::releaseObjectsInStdVector(_vec_mat_assembler);
    };

    /**
     * add a local assembler
     * @param mat_id
     * @param assembler
     */
    void addLocalAssembler(IElementWiseTransientLinearEQSLocalAssembler* assembler)
    {
        _vec_mat_assembler.push_back(assembler);
    }

    /**
     * return a local assembler
     * @param mat_id
     * @return
     */
    IElementWiseTransientLinearEQSLocalAssembler* getLocalAssembler(size_t id)
    {
        return _vec_mat_assembler[id];
    }

    /// assemble a local linear equation for the given element
    /// @param time            time step
    /// @param e            element
    /// @param local_u_n1    guess of current time step value
    /// @param local_u_n    previous time step value
    /// @param eqs            local algebraic equation
    virtual void assembly(const TimeStep &time, const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &localDofManager, const MathLib::LocalVector &local_u_n1, const MathLib::LocalVector &local_u_n, MathLib::LocalEquation &eqs)
    {
        IElementWiseTransientLinearEQSLocalAssembler* assembler = getLocalAssembler(e.getGroupID());
        assembler->assembly(time, e, localDofManager, local_u_n1, local_u_n, eqs);
    }


private:
    std::vector<IElementWiseTransientLinearEQSLocalAssembler*> _vec_mat_assembler;
};



} //end
