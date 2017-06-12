/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemLowerDimension.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Core/ElementCoordinatesMappingLocal.h"
#include "FemLib/Core/ShapeFunction/IFemShapeFunction.h"
#include "FemNaturalCoordinates.h"

namespace FemLib
{

/**
 * \brief Mapping element shapes to natural coordinates in element dimension
 *
 * Shape functions are not directional but derivatives of shape functions are
 * directional.
 * Any directional functions have to be transformed to appropriate coordinate
 * systems when given element dimension and a mesh dimension is not the same.
 *
 * In general, coordinate transformation can be done in the following way.
 * Consider that x is a point in global coordinates (x-y-z) and x' is that
 * in local coordinates (x'-y'). Coordinate transformation between the two
 * coordinates can be achieved with a rotational matrix R as
 *   x = R x'
 *
 * Hence, to transform derivatives of shape functions in natural coordinates to
 * those in global coordinates, which are directional,
 *   dN(x) = J(x)^-1 dN(x') = J(x)^-1 J(x')^-1 dN(r) = R^T (J(x')^-1 dN(r))
 */
class FemLowerDimension
: public FemNaturalCoordinates
{
public:

    /**
     *
     * @param shape
     * @param mesh_dim
     */
    FemLowerDimension(IFemShapeFunction *shape, size_t mesh_dim)
    : FemNaturalCoordinates(shape), _msh_dim(mesh_dim), _ele_dim(0)
    {
    }

    ///
    virtual ~FemLowerDimension()
    {
    };

    /// initialize element
    virtual void initialize(MeshLib::IElement &ele)
    {
        assert (ele.getDimension() < _msh_dim);

        _ele_dim = ele.getDimension();
        FemNaturalCoordinates::initialize(ele);
    }

    /// compute shape functions
    virtual const CoordinateMappingProperty* compute(const double* natural_pt)
    {
        // calculate mapping functions in natural coordinates
        CoordinateMappingProperty* prop = const_cast<CoordinateMappingProperty*>(FemNaturalCoordinates::compute(natural_pt));
        // modify derivatives of shape functions
        const MathLib::LocalMatrix &matR = getRotationMatrix();
        MathLib::LocalMatrix dshape_local = MathLib::LocalMatrix::Zero(matR.rows(), prop->dshape_dx->cols());
        dshape_local.topLeftCorner(prop->dshape_dx->rows(), prop->dshape_dx->cols()) = (*prop->dshape_dx);
        (*prop->dshape_dx) =  matR * dshape_local;
        //(*prop->dshape_dx) =  matR.transpose() * dshape_local;

        return prop;
    }

    /// map natural coordinates to physical coordinates
    virtual void mapToPhysicalCoordinates(const double* natural_pt, double* physical_pt)
    {
        MathLib::LocalVector x_local(_msh_dim);
        FemNaturalCoordinates::mapToPhysicalCoordinates(natural_pt, (double*)&x_local[0]);

        const MathLib::LocalMatrix &matR = getRotationMatrix();
        assert((unsigned)matR.rows() == _msh_dim);
        MathLib::LocalVector x_global(_msh_dim);
        x_global = matR * x_local;
        for (size_t i=0; i<_msh_dim; i++)
            physical_pt[i] = x_global[i];
    }

    /// map physical coordinates to natural coordinates
    virtual void mapFromPhysicalCoordinates(const double* physical_pt, double* natural_pt)
    {
        MathLib::LocalVector x_global(_msh_dim);
        for (size_t i=0; i<_msh_dim; i++)
            x_global[i] = physical_pt[i];
        const MathLib::LocalMatrix &matR = getRotationMatrix();
        assert(static_cast<size_t>(matR.rows()) == _msh_dim);
        MathLib::LocalVector x_local(_msh_dim);
        x_local = matR.transpose() * x_global;

        FemNaturalCoordinates::mapFromPhysicalCoordinates(x_local.data(), natural_pt);
    }

private:
    const MathLib::LocalMatrix& getRotationMatrix()
    {
        MeshLib::IElement* ele = FemNaturalCoordinates::getElement();
        MeshLib::ElementCoordinatesMappingLocal* ele_local_coord;
        ele_local_coord = (MeshLib::ElementCoordinatesMappingLocal*)ele->getMappedCoordinates();
        const MathLib::LocalMatrix &matR = ele_local_coord->getRotationMatrixToOriginal();
        return matR;
    }
private:
    size_t _msh_dim;
    size_t _ele_dim;
};

}
