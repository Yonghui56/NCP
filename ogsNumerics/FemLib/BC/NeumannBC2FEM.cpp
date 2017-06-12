/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NeumannBC2FEM.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "NeumannBC2FEM.h"

#include <map>

#include "GeoLib/Point.h"
#include "MeshLib/Tools/Tools.h"
#include "MathLib/DataType.h"

namespace FemLib
{

NeumannBC2FEM::NeumannBC2FEM(const MeshLib::IMesh &msh, const double &current_time, IFeObjectContainer &feObjects, const GeoLib::GeoObject &_geo, const NumLib::ITXFunction &_bc_func, std::vector<size_t> &_vec_nodes, std::vector<double> &_vec_values)
{
    // pickup nodes on geo
    MeshLib::findNodesOnGeometry(&msh, &_geo, &_vec_nodes);

    // distribute to RHS
    _vec_values.resize(_vec_nodes.size());
    switch (_geo.getGeoType()) {
        case  GeoLib::POINT:
        {
            // no need to integrate
            // get discrete values at nodes
            for (size_t i=0; i<_vec_nodes.size(); i++) {
                const GeoLib::Point* x = msh.getNodeCoordinatesRef(_vec_nodes[i]);
                NumLib::TXPosition pos(current_time, x->getData());
                _bc_func.eval(pos, _vec_values[i]);
            }
            break;
        }
        case GeoLib::POLYLINE:
        {
            // find edge elements on the geo
            std::vector<MeshLib::IElement*> vec_edge_eles;
            MeshLib::findBoundaryElementsOnGeometry(const_cast< MeshLib::IMesh*>(&msh), &_geo, &vec_edge_eles);
            // for each edge elements found
            std::map<size_t, double> map_nodeId2val;
            for (size_t i=0; i<vec_edge_eles.size(); i++) {
                MeshLib::IElement *e = vec_edge_eles[i];
                e->setCurrentOrder(msh.getCurrentOrder());
                const size_t edge_nnodes = e->getNumberOfNodes();
                // set values at nodes
                MathLib::LocalVector  nodal_val(edge_nnodes);
                for (size_t i_nod=0; i_nod<edge_nnodes; i_nod++) {
                    const GeoLib::Point* x = msh.getNodeCoordinatesRef(e->getNodeID(i_nod));
                    NumLib::TXPosition pos(current_time, x->getData());
                    double v = .0;
                    _bc_func.eval(pos, v);
                    nodal_val[i_nod] = v;
                }
                // compute integrals
                e->setCurrentOrder(msh.getCurrentOrder());
                IFiniteElement *fe_edge = feObjects.getFeObject(*e);
                fe_edge->getIntegrationMethod()->initialize(*e, 2);
                //IFiniteElement *fe_edge = _var->getFiniteElement(*e);
                MathLib::LocalVector  result(edge_nnodes);
                MathLib::LocalMatrix M = MathLib::LocalMatrix::Zero(edge_nnodes, edge_nnodes);
                IFemNumericalIntegration *q = fe_edge->getIntegrationMethod();
                double x_ref[3];
                MathLib::LocalMatrix fac(1,1);
                fac(0,0) = 1;
                for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
                    q->getSamplingPoint(j, x_ref);
                    fe_edge->computeBasisFunctions(x_ref);
                    fe_edge->integrateWxN(j, fac, M);
                }
                result = M * nodal_val;
                //M.axpy(1.0, &nodal_val[0], 0.0, &result[0]);
                // add into RHS values
                for (size_t k=0; k<edge_nnodes; k++)
                    map_nodeId2val[e->getNodeID(k)] += result[k];
            }
            for (size_t i=0; i<_vec_nodes.size(); i++) {
                _vec_values[i] = map_nodeId2val[_vec_nodes[i]];
            }
            break;
        }
        default:
            throw "Given GeoType is not supported yet in NeumannBC2FEM";
            break;
    }
}

}
