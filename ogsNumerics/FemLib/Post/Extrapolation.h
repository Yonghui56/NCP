/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Extrapolation.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <cmath>
#include <limits>

#include "MeshLib/Topology/TopologySequentialNodes2Elements.h"
#include "FemLib/Function/FemFunction.h"
#include "FemLib/Core/Integration/Integration.h"

namespace FemLib
{
/**
 * \brief Extrapolation methods for integration point values
 */
template<typename T_DIS_SYS, typename Tvalue>
class IFemExtrapolation
{
public:
    virtual ~IFemExtrapolation() {};
    virtual void extrapolate(TemplateFEMIntegrationPointFunction<T_DIS_SYS,Tvalue> &ele, TemplateFEMNodalFunction<T_DIS_SYS,Tvalue> &nod) = 0;
};

/**
 * \brief
 */
template<typename T_DIS_SYS, typename Tvalue>
class FemExtrapolationAverage : public IFemExtrapolation<T_DIS_SYS, Tvalue>
{
public:

    void extrapolate(TemplateFEMIntegrationPointFunction<T_DIS_SYS,Tvalue> &ele_var, TemplateFEMNodalFunction<T_DIS_SYS,Tvalue> &nod_var)
    {
        const MeshLib::IMesh* msh = ele_var.getMesh();
        IFeObjectContainer* feObjects = nod_var.getFeObjectContainer();
        MeshLib::TopologySequentialNodes2Elements node2eles(*msh);
        Tvalue v0 = ele_var.getIntegrationPointValues(0)[0];
        v0 *= .0;
        DiscreteLib::IDiscreteVector<Tvalue>* node_vec = nod_var.getDiscreteData();
        (*node_vec) = v0;

        for (size_t i=0; i<msh->getNumberOfElements(); i++) {
            MeshLib::IElement* e = msh->getElement(i);
            const typename TemplateFEMIntegrationPointFunction<T_DIS_SYS, Tvalue>::IntegrationPointVectorType &gp_values = ele_var.getIntegrationPointValues(i);
            std::vector<Tvalue> vec_gp_values(&gp_values[0], &gp_values[0]+gp_values.size());
            const size_t e_nnodes = e->getNumberOfNodes();
            std::vector<Tvalue> nodal_values(e_nnodes);
            IFiniteElement *fe = feObjects->getFeObject(*e);
            fe->extrapolate(vec_gp_values, nodal_values);

            for (size_t j=0; j<e_nnodes; j++) {
                const size_t nod_id = e->getNodeID(j);
                size_t n_conn_eles = node2eles.getConnectedElements(nod_id).size();
                (*node_vec)[nod_id] += nodal_values[j] / static_cast<double>(n_conn_eles);
            }
        }
    }
};


struct FEMExtrapolationMethod
{
    enum type {
        Average,
        INVALID
    };
};

template<typename T_DIS_SYS, typename Tvalue>
class FEMExtrapolationFactory
{
public:
    static IFemExtrapolation<T_DIS_SYS, Tvalue>* create(FEMExtrapolationMethod::type tp)
    {
        switch (tp) {
            case FEMExtrapolationMethod::Average:
                return new FemExtrapolationAverage<T_DIS_SYS,Tvalue>();
            default:
                return NULL;
        }
        return NULL;
    };
};


}
