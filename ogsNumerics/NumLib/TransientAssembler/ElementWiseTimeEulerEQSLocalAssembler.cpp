/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementWiseTimeEulerEQSLocalAssembler.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "ElementWiseTimeEulerEQSLocalAssembler.h"

namespace NumLib
{

/// assemble a local linear equation for the given element
/// @param time            time step
/// @param e            element
/// @param local_u_n1    guess of current time step value
/// @param local_u_n    previous time step value
/// @param eqs            local algebraic equation
void ElementWiseTimeEulerEQSLocalAssembler::assembly(const TimeStep &time, const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &/*localDofManager*/, const MathLib::LocalVector &local_u_n1, const MathLib::LocalVector &local_u_n, MathLib::LocalEquation &eqs)
{
    const double delta_t = time.getTimeStepSize();
    const size_t n_dof = eqs.getDimension();

    MathLib::LocalMatrix M = MathLib::LocalMatrix::Zero(n_dof, n_dof);
    MathLib::LocalMatrix K = MathLib::LocalMatrix::Zero(n_dof, n_dof);
    MathLib::LocalVector F = MathLib::LocalVector::Zero(n_dof);

    assembleODE(time, e, local_u_n1, local_u_n, M, K, F);

    //std::cout << "M="; M.write(std::cout); std::cout << std::endl;
    //std::cout << "K="; K.write(std::cout); std::cout << std::endl;

    MathLib::LocalMatrix *localA = eqs.getA();
    MathLib::LocalVector *localRHS = eqs.getRHSAsVec();

    // A = 1/dt M + theta K
    (*localA) = 1.0/delta_t * M + _theta * K;
    // RHS = (1/dt M - (1-theta) K) u0 + F
    (*localRHS) = (1.0/delta_t * M - (1.-_theta) * K) * local_u_n;
    (*localRHS) += F;
}


} //end
