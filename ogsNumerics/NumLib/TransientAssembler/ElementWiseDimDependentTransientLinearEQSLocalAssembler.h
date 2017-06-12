/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementWiseDimDependentTransientLinearEQSLocalAssembler.h
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
 * \brief Element dimension dependent linear EQS local assembler
 *
 * @tparam  T_USER_ASSEMBLY     User-given assembler
 */
class ElementWiseDimDependentTransientLinearEQSLocalAssembler
: public IElementWiseTransientLinearEQSLocalAssembler
{
public:
    ElementWiseDimDependentTransientLinearEQSLocalAssembler()
    : _vec_dim_assembler(3, nullptr) {};

    virtual ~ElementWiseDimDependentTransientLinearEQSLocalAssembler()
    {
//        BaseLib::releaseObjectsInStdVector(_vec_dim_assembler);
    };

    /**
     * add a local assembler
     * @param mat_id
     * @param assembler
     */
    void addLocalAssembler(size_t i_dim, IElementWiseTransientLinearEQSLocalAssembler* assembler)
    {
        assert(0 < i_dim && i_dim < 4);
        _vec_dim_assembler[i_dim-1] = assembler;
    }

    /**
     * return a local assembler
     * @param mat_id
     * @return
     */
    IElementWiseTransientLinearEQSLocalAssembler* getLocalAssembler(size_t i_dim)
    {
        assert(0 < i_dim && i_dim < 4);
        return _vec_dim_assembler[i_dim-1];
    }

    /// assemble a local linear equation for the given element
    /// @param time            time step
    /// @param e            element
    /// @param local_u_n1    guess of current time step value
    /// @param local_u_n    previous time step value
    /// @param eqs            local algebraic equation
    virtual void assembly(const TimeStep &time, const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &localDofManager, const MathLib::LocalVector &local_u_n1, const MathLib::LocalVector &local_u_n, MathLib::LocalEquation &eqs)
    {
        IElementWiseTransientLinearEQSLocalAssembler* assembler = getLocalAssembler(e.getDimension());
        assembler->assembly(time, e, localDofManager, local_u_n1, local_u_n, eqs);
    }


private:
    std::vector<IElementWiseTransientLinearEQSLocalAssembler*> _vec_dim_assembler;
};



} //end
