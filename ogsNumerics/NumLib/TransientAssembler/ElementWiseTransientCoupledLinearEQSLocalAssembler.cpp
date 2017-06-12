/* Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementWiseTransientCoupledLinearEQSLocalAssembler.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "ElementWiseTransientCoupledLinearEQSLocalAssembler.h"

#include "DiscreteLib/Utils/Tools.h"

namespace NumLib
{

void ElementWiseTransientCoupledLinearEQSLocalAssembler::assembly(  const NumLib::TimeStep &timestep,
                        const MeshLib::IElement &e,
                        const DiscreteLib::DofEquationIdTable &localDofManager,
                        const LocalVectorType &local_u_n1,
                        const LocalVectorType &local_u_n,
                        MathLib::LocalEquation &eqs
                        )
{
    // parameters need to be passed
    const size_t n_var = _n_var;
    std::vector<size_t> &vec_order = _vec_order;
    size_t mesh_id = 0;

    //
    std::vector<std::vector<size_t> > vec_local_pos(n_var);
    for (size_t i=0; i<n_var; i++) {
        std::vector<size_t> list_nodeid;
        e.getNodeIDList(vec_order[i], list_nodeid);
        localDofManager.mapEqsID(i, mesh_id, list_nodeid, vec_local_pos[i]);
    }

    std::vector<LocalVectorType> vec_u0(n_var), vec_u1(n_var);
    for (size_t i=0; i<n_var; i++) {
        DiscreteLib::getLocalVector(vec_local_pos[i], local_u_n, vec_u0[i]);
        DiscreteLib::getLocalVector(vec_local_pos[i], local_u_n1, vec_u1[i]);
    }

    std::vector<std::vector<LocalMatrixType> > vec_K(n_var);
    std::vector<LocalVectorType> vec_F(n_var);
    for (size_t i=0; i<n_var; i++) {
        for (size_t j=0; j<n_var; j++)
            vec_K[i].push_back(LocalMatrixType::Zero(vec_local_pos[i].size(), vec_local_pos[j].size()));
        //vec_K[i].resize(n_var);
    }


    this->assembleComponents(timestep, e, vec_order, vec_u0, vec_u1, vec_K, vec_F);

    for (size_t i=0; i<n_var; i++) {
        for (size_t j=0; j<n_var; j++) {
            eqs.addAsub(vec_local_pos[i], vec_local_pos[j], vec_K[i][j]);
        }
        eqs.addRHSsub(vec_local_pos[i], vec_F[i]);
    }

}


} //end
