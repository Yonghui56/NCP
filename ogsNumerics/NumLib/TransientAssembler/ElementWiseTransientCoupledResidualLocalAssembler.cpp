/* Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementWiseTransientCoupledResidualLocalAssembler.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "ElementWiseTransientCoupledResidualLocalAssembler.h"

#include "DiscreteLib/Utils/Tools.h"

namespace NumLib
{

void ElementWiseTransientCoupledResidualLocalAssembler::assembly(  const NumLib::TimeStep &timestep,
                        const MeshLib::IElement &e,
                        const DiscreteLib::DofEquationIdTable &localDofManager,
                        const LocalVectorType &local_u_n1,
                        const LocalVectorType &local_u_n,
                        LocalVectorType &local_r
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

    std::vector<LocalVectorType> vec_r(n_var);

    this->assembleComponents(timestep, e, vec_order, vec_u0, vec_u1, vec_r);

    for (size_t i=0; i<n_var; i++) {
        for (size_t j=0; j<vec_local_pos[i].size(); j++) {
            local_r[vec_local_pos[i][j]] += vec_r[i][j];
        }
    }

//    if (e.getID()==0) {
//        std::cout << "local_u_n0=" << local_u_n << std::endl;
//        std::cout << "local_u_n1=" << local_u_n1 << std::endl;
//        std::cout << "local_r=" << local_r << std::endl;
//    }
}


} //end
