/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file StencilWiseTransientLinearEQSAssembler.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "MathLib/LinAlg/LinearEquation/ILinearEquation.h"
#include "MeshLib/Topology/TopologyNode2NodesConnectedByEdges.h"

#include "DiscreteLib/Utils/DofEquationIdTable.h"
#include "DiscreteLib/Serial/DiscreteVector.h"
#include "DiscreteLib/Core/IDiscreteLinearEquationAssembler.h"
#include "DiscreteLib/Utils/Tools.h"
#include "IStencilWiseTransientLinearEQSLocalAssembler.h"
#include "IStencil.h"

namespace MeshLib
{
class IMesh;
}

namespace NumLib
{

class TimeStep;
}

namespace FdmLib
{

/**
 * \brief Element-based discrete system assembler classes
 */
class StencilWiseTransientLinearEQSAssembler : public DiscreteLib::IDiscreteLinearEquationAssembler
{
public:
    typedef IStencilWiseTransientLinearEQSLocalAssembler::LocalEquationType LocalEquationType;
    typedef IStencilWiseTransientLinearEQSLocalAssembler::LocalMatrixType LocalMatrixType;
    typedef IStencilWiseTransientLinearEQSLocalAssembler::LocalVectorType LocalVectorType;

    /// @param time
    /// @param u0
    /// @param u1
    /// @param a
    StencilWiseTransientLinearEQSAssembler(const NumLib::TimeStep* time, const std::vector<DiscreteLib::IDiscreteVector<double>*>* u0, const std::vector<DiscreteLib::IDiscreteVector<double>*>* u1, IStencilWiseTransientLinearEQSLocalAssembler* a)
        : _transient_e_assembler(a), _timestep(time), _u0(u0), _u1(u1)
    { };


    /// Conduct the element by element assembly procedure
    ///
    /// @param msh                 Mesh
    /// @param dofManager         Dof map manager
    /// @param list_dofId         List of Dof IDs used in this problem
    /// @param eqs                 Linear equation solver
    void assembly(const MeshLib::IMesh &msh, const DiscreteLib::DofEquationIdTable &dofManager, MathLib::ILinearEquation &eqs)
    {
        const NumLib::TimeStep &time = *_timestep;
        LocalEquationType localEQS;
        std::vector<size_t> ele_node_ids, ele_node_size_order;
        std::vector<size_t> local_dofmap;
        const size_t n_nod = msh.getNumberOfNodes();
        MeshLib::TopologyNode2NodesConnectedByEdges topo(msh);

        LocalVectorType local_u_n1;
        LocalVectorType local_u_n;
        for (size_t i=0; i<n_nod; i++) {
            MeshLib::Node* nod = msh.getNode(i);
            ele_node_ids.clear();
            ele_node_size_order.clear();
            local_dofmap.clear();
            Stencil5 stencil;
            // set stencil
            stencil.setCentralNodeID(nod->getNodeID());
            ele_node_ids.push_back(nod->getNodeID());
            const std::set<size_t> &connected_nodes = topo.getConnectedNodes(i);
            for (std::set<size_t>::const_iterator itr=connected_nodes.begin(); itr!=connected_nodes.end(); ++itr) {
                stencil.addSurroundingNode(*itr);
                ele_node_ids.push_back(*itr);
            }
            ele_node_size_order.push_back(ele_node_ids.size());

            dofManager.mapEqsID(msh.getID(), ele_node_ids, local_dofmap);
            // previous and current results
            DiscreteLib::getLocalVector(dofManager, ele_node_ids, ele_node_size_order, *_u1, local_u_n1);
            DiscreteLib::getLocalVector(dofManager, ele_node_ids, ele_node_size_order, *_u0, local_u_n);
            // local assembly
            localEQS.create(local_dofmap.size());
            _transient_e_assembler->assembly(time, stencil, local_u_n1, local_u_n, localEQS);
            // update global
            eqs.addAsub(local_dofmap, *localEQS.getA());
            eqs.addRHSsub(local_dofmap, localEQS.getRHS());
        }
    }

private:
    IStencilWiseTransientLinearEQSLocalAssembler* _transient_e_assembler;
    const NumLib::TimeStep* _timestep;
    const std::vector<DiscreteLib::IDiscreteVector<double>*>* _u0;
    const std::vector<DiscreteLib::IDiscreteVector<double>*>* _u1;
};


}
