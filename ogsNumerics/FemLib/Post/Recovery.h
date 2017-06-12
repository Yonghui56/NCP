/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Recovery.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include "MeshLib/Topology/Topology.h"
#include "FemLib/Function/FemFunction.h"

namespace FemLib
{

/**
 * \brief A class to recover solution derivatives from one-order lower solutions
 */
class IFemRecovery {};

/**
 * \brief Weighted average method
 */
class FemRecoveryWeightedAverage : IFemRecovery
{
public:
    void compute(FemNodalFunctionScalar* u, FemNodalFunctionScalar* du)
    {
        const MeshLib::IMesh* msh = u->getMesh();

        std::vector<double> vec_dN_gp(msh->getNumberOfElements());
        for (size_t i=0; i<msh->getNumberOfElements(); i++) {
            MeshLib::IElement* e = msh->getElement(i);
            double *gp = 0; //gravity center of the element
            IFiniteElement* fe = u->getFiniteElement(e);
            double *local_ui = 0;
            fe->computeBasisFunctions(gp);
            MathLib::Matrix<double> *dN = fe->getGradBasisFunction();
            double dui = 0;
            //dN->multiply(local_ui, dui);
            vec_dN_gp[i] = dui;
        }

        std::vector<double> vec_du;
        MeshLib::TopologySequentialNodes2Elements node2ele(msh);
        for (size_t i=0; i<msh->getNumberOfNodes(); i++) {
            const std::vector<size_t> &eles = node2ele.getConnectedElements(i);
            for (size_t j=0; j<eles.size(); j++) {
                double ele_dui = vec_dN_gp[eles[j]];
                double r = .0; // dist(x_i, gp_j)
                vec_du[i] += ele_dui / r;
            }
            vec_du[i] /= eles.size();
        }

        du->setNodalValues(&vec_du[0]);
    }
};


class FemRecoveryLocalPatch : IFemRecovery
{
public:

};

}
