/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementWiseTimeEulerResidualLocalAssembler.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "ElementWiseTimeEulerResidualLocalAssembler.h"


namespace NumLib
{

void ElementWiseTimeEulerResidualLocalAssembler::assembly(const TimeStep &time, const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &/*localDofManager*/, const MathLib::LocalVector  &local_u_n1, const MathLib::LocalVector  &local_u_n, MathLib::LocalVector  &local_r)
{
    const double delta_t = time.getTimeStepSize();
    const size_t n_dof = local_r.size();
	size_t elem_id(0);
	elem_id = e.getGroupID();
    MathLib::LocalMatrix M = MathLib::LocalMatrix::Zero(n_dof, n_dof);
    MathLib::LocalMatrix K = MathLib::LocalMatrix::Zero(n_dof, n_dof);
    MathLib::LocalVector  F = MathLib::LocalVector ::Zero(n_dof);

    // get M,K,F in M du/dt + K = F
    assembleODE(time, e, local_u_n1, local_u_n, M, K, F);

    MathLib::LocalMatrix TMP_M(n_dof, n_dof);
    //MathLib::LocalVector  TMP_V(n_dof);

    // evaluate r: r = (1/dt M + theta K) * u1 -(1/dt M - (1-theta) K) * u0 - F
    // r = (1/dt M + theta K) * u1
    TMP_M = 1.0/delta_t * M + _theta * K;
	//std::cout << TMP_M << std::endl;
    local_r = TMP_M * local_u_n1;
    // r -= (1/dt M - (1-theta) K) u0
    TMP_M = 1.0/delta_t * M - (1.-_theta) * K;
	//std::cout << TMP_M << std::endl;
    //local_r.noalias() -= TMP_M * local_u_n;
	local_r.noalias() -= TMP_M * local_u_n;
	// r -= F
	local_r.noalias() -= F;
	/*if (local_u_n1.size()<9)
	{
		
	}*/
	

    
	// debugging
	//std::cout << "r=" << local_r << std::endl;
	//
#if 0
    if (e.getID()<1) {
        std::cout << "Element: " << e.getID() << std::endl;
        std::cout << "M=" << M << std::endl;
        std::cout << "K=" << K << std::endl;
        std::cout << "u1=" << local_u_n1 << std::endl;
        std::cout << "u0=" << local_u_n << std::endl;
        TMP_M = 1.0/delta_t * M + _theta * K;
        std::cout << "A=" << TMP_M << std::endl;
        TMP_M = 1.0/delta_t * M - (1.-_theta) * K;
        std::cout << "B=" << TMP_M << std::endl;
        std::cout << "r=" << local_r << std::endl;
    }

#endif
}


} //end
