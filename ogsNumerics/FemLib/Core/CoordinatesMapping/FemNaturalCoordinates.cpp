/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemNaturalCoordinates.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "FemNaturalCoordinates.h"

#include <iostream>
#include <cassert>

namespace FemLib
{

/// initialize element
void FemNaturalCoordinates::initialize(MeshLib::IElement &ele)
{
    assert(ele.getMappedCoordinates()!=NULL);

    _ele = &ele;
    const size_t dim = _ele->getDimension();
    const size_t nnodes = _ele->getNumberOfNodes();
    _prop->shape_r->resize(1, nnodes);
    _prop->dshape_dr->resize(dim, nnodes);
    _prop->dshape_dx->resize(dim, nnodes);
    _prop->jacobian_dxdr->resize(dim, dim);
    _prop->inv_jacobian_drdx->resize(dim, dim);
    (*_prop->shape_r) *= 0;
    (*_prop->dshape_dr) *= 0;
    (*_prop->dshape_dx) *= 0;
    (*_prop->jacobian_dxdr) *= 0;
    (*_prop->inv_jacobian_drdx) *= 0;
    x.resize(dim, nnodes);
    MeshLib::IElementCoordinatesMapping *ele_local_coord = _ele->getMappedCoordinates();
	/*if (dim==1)
	{
		for (size_t i = 0; i<nnodes; i++)
			for (size_t j = 0; j<dim; j++)
				x(j, i) = ele_local_coord->getNodePoint(i)->getData()[j];
	}
	else
	{
		for (size_t i = 0; i < nnodes; i++){
			//for (size_t j = 0; j<dim; j++)
			x(0, i) = ele_local_coord->getNodePoint(i)->getData()[0];
			x(1, i) = ele_local_coord->getNodePoint(i)->getData()[2];
		}
	}*/
	for (size_t i = 0; i<nnodes; i++)
		for (size_t j = 0; j<dim; j++)
			x(j, i) = ele_local_coord->getNodePoint(i)->getData()[j];
    

}

/// compute mapping properties at the given location in natural coordinates
const CoordinateMappingProperty* FemNaturalCoordinates::compute(const double* natural_pt)
{
    //prepare
    const size_t dim = _ele->getDimension();
    const size_t nnodes = _ele->getNumberOfNodes();

    //shape, dshape
    MathLib::LocalMatrix *shape  = _prop->shape_r;
    MathLib::LocalMatrix *dshape_dr = _prop->dshape_dr;
    _shape->computeShapeFunction(natural_pt, &(*shape)(0,0));
    _shape->computeGradShapeFunction(natural_pt, &(*dshape_dr)(0,0));
//    _shape->computeShapeFunction(natural_pt, (double*)shape->getData());
//    _shape->computeGradShapeFunction(natural_pt, (double*)dshape_dr->getData());

    //jacobian: J=[dx/dr dy/dr // dx/ds dy/ds]
    MathLib::LocalMatrix *jac = _prop->jacobian_dxdr;
    (*jac) = MathLib::LocalMatrix::Zero(jac->rows(), jac->cols());
    for (size_t i_r=0; i_r<dim; i_r++) {
        for (size_t j_x=0; j_x<dim; j_x++) {
            for (size_t k=0; k<nnodes; k++) {
                (*jac)(i_r,j_x) += (*dshape_dr)(i_r, k) * x(j_x,k);
            }
        }
    }

    //determinant of J
    double det_j = jac->determinant();
    _prop->det_jacobian = det_j;
    if (det_j>.0) {
        //inverse of J
        MathLib::LocalMatrix *inv_jac = _prop->inv_jacobian_drdx;
        //jac->inverse(inv_jac);
        (*inv_jac) = jac->inverse();

        //dshape/dx = invJ * dNdr
        MathLib::LocalMatrix *dshape_dx = _prop->dshape_dx;
//        (*dshape_dx) = .0;
//        inv_jac->multiply(*dshape_dr, *dshape_dx);
        (*dshape_dx) = (*inv_jac) * (*dshape_dr);
    } else {
        std::cout << "***error: det_j is not positive. element id=" << _ele->getID() << std::endl;
    }

    return _prop;
};

/// compute physical coordinates at the given natural coordinates
/// \f[
///    \mathbf{x} = \mathbf{N(r)} * \mathbf{X}
/// \f]
///
/// @param natural_pt
/// @return physical_pt
void FemNaturalCoordinates::mapToPhysicalCoordinates(const double* natural_pt, double* physical_pt)
{
    const size_t dim = _ele->getDimension();
    const size_t nnodes = _ele->getNumberOfNodes();
    MathLib::LocalMatrix *shape  = _prop->shape_r;
    _shape->computeShapeFunction(natural_pt, &(*shape)(0,0));

    for (size_t i=0; i<dim; i++) {
        physical_pt[i] = .0;
        for (size_t j=0; j<nnodes; j++)
            physical_pt[i] += (*shape)(0,j) + x(i,j);
    }
}

/// compute physical coordinates at the given natural coordinates
/// \f[
///    \mathbf{x} = \mathbf{N(r)} * \mathbf{X}
/// \f]
///
/// @param natural_pt
/// @return physical_pt
void FemNaturalCoordinates::mapToPhysicalCoordinates(const CoordinateMappingProperty* prop, double* physical_pt)
{
    const size_t dim = _ele->getDimension();
    const size_t nnodes = _ele->getNumberOfNodes();
    MathLib::LocalMatrix *shape  = prop->shape_r;

    for (size_t i=0; i<dim; i++) {
        physical_pt[i] = .0;
        for (size_t j=0; j<nnodes; j++)
            physical_pt[i] += (*shape)(0,j) * x(i,j);
    }
}

/// compute natural coordinates at the given natural coordinates.
/// Assuming \f$ r=0 \f$ at \f$ x = \bar{x}_{avg} \f$, natural coordinates can be calculated as
/// \f[
///    \mathbf{r} = (\mathbf{J}^-1)^T * (\mathbf{x} - \bar{\mathbf{x}}_{avg})
/// \f]
///
/// @param physical_pt
/// @return natural_pt
void FemNaturalCoordinates::mapFromPhysicalCoordinates(const double* physical_pt, double* natural_pt)
{
    const size_t dim = _ele->getDimension();
    const size_t nnodes = _ele->getNumberOfNodes();
    MathLib::LocalMatrix *inv_jac = _prop->inv_jacobian_drdx;

    // calculate dx which is relative coordinates from element center
    std::vector<double> dx(dim, .0);
    // x_avg = sum_i {x_i} / n
    for (size_t i=0; i<dim; i++)
        for (size_t j=0; j<nnodes; j++)
            dx[i] += x(i,j);
    for (size_t i=0; i<dim; i++)
        dx[i] /= (double)nnodes;
    // dx = pt - x_avg
    for (size_t i=0; i<dim; i++)
        dx[i] = physical_pt[i] -dx[i];

    // r = invJ^T * dx
    for (size_t i=0; i<dim; i++) {
        natural_pt[i] = 0.0;
        for (size_t j=0; j<dim; j++)
            natural_pt[i] += (*inv_jac)(j * dim, i) * dx[j];
    }

}

} //end
