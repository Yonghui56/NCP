/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemNaturalCoordinates.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Core/IElementCoordinatesMapping.h"
#include "FemLib/Core/ShapeFunction/IFemShapeFunction.h"
#include "IFemCoordinatesMapping.h"

namespace FemLib
{

/**
 * \brief Mapping element shapes to natural coordinates
 *
 * FemNaturalCoordinates mapping converts element shapes in physical coordinates
 * (x,y,z) to that in natural coordinates (r,s,t).
 * - Given physical coordinates should correspond to dimensions of the element, 
 *   i.e. (x,y) for triangles
 * - x(r,s,t) = N(r,s,t) x_i
 */
class FemNaturalCoordinates : public IFemCoordinatesMapping
{
public:
    ///
    explicit FemNaturalCoordinates(IFemShapeFunction *shape) 
    {
        _shape = shape;
        _prop = new CoordinateMappingProperty();
        _ele = NULL;
    }

    ///
    virtual ~FemNaturalCoordinates()
    {
        delete _shape;
        delete _prop;
    }

    /// initialize element
    virtual void initialize(MeshLib::IElement &ele);

    ///
    virtual const CoordinateMappingProperty* getProperties() const {return _prop;};

    /// compute mapping properties at the given location in natural coordinates
    virtual const CoordinateMappingProperty* compute(const double* natural_pt);

    /// compute physical coordinates at the given natural coordinates
    /// \f[
    ///    \mathbf{x} = \mathbf{N(r)} * \mathbf{X}
    /// \f]
    ///
    /// @param natural_pt
    /// @return physical_pt
    void mapToPhysicalCoordinates(const double* natural_pt, double* physical_pt);

    /// compute physical coordinates at the given natural coordinates
    /// \f[
    ///    \mathbf{x} = \mathbf{N(r)} * \mathbf{X}
    /// \f]
    ///
    /// @param natural_pt
    /// @return physical_pt
    void mapToPhysicalCoordinates(const CoordinateMappingProperty* prop, double* physical_pt);

    /// compute natural coordinates at the given natural coordinates.
    /// Assuming \f$ r=0 \f$ at \f$ x = \bar{x}_{avg} \f$, natural coordinates can be calculated as
    /// \f[
    ///    \mathbf{r} = (\mathbf{J}^-1)^T * (\mathbf{x} - \bar{\mathbf{x}}_{avg})
    /// \f]
    ///
    /// @param physical_pt
    /// @return natural_pt
    void mapFromPhysicalCoordinates(const double* physical_pt, double* natural_pt);

    ///
    MeshLib::IElement* getElement() const {return _ele;};


private:
    CoordinateMappingProperty* _prop;
    IFemShapeFunction* _shape;
    MeshLib::IElement* _ele;
    MathLib::LocalMatrix x;

};


}
