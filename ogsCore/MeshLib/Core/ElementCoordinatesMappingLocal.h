/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementCoordinatesMappingLocal
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include "Eigen"

#include "BaseLib/CodingTools.h"
#include "MathLib/DataType.h"
#include "GeoLib/Point.h"

#include "MeshLib/Core/IElementCoordinatesMapping.h"
#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Core/IElement.h"
#include "MeshLib/Core/CoordinateSystem.h"

namespace MeshLib
{

/**
 * \brief This classs mapps coordinates of element nodes to their intrinsic coordinates which is lower dimension
 *
 * For example, a triangle element in 2.5D is mapped into 2D.
 */
class ElementCoordinatesMappingLocal : public IElementCoordinatesMapping 
{
public:
    /**
     * 
     * \param msh                   Pointer to a mesh object
     * \param e                     Mesh element
     * \param org_coord_system      Original coordinate system
     */
    ElementCoordinatesMappingLocal(const IMesh* msh, IElement &e, const CoordinateSystem &org_coord_system);

    ///
    virtual ~ElementCoordinatesMappingLocal()
    {
    }

    /// return mapped coordinates of the node
    virtual GeoLib::Point* getNodePoint(size_t node_id) 
    {
        return &_point_vec[node_id];
    }

    /// return a rotation matrix converting to orinal coordinates
    const MathLib::LocalMatrix& getRotationMatrixToOriginal() const {return _matR2original;};

    /// return a rotation matrix converting to local coordinates
    const MathLib::LocalMatrix& getRotationMatrixToLocal() const {return _matR2local;};

    /// return a translation vector
    const GeoLib::Point& getTranslationVector() const {return _pt_translate;};

private:
    ///
    void translate(std::vector<GeoLib::Point> &point_vec);
    ///
    void flip(IElement &e, const CoordinateSystem &coordinate_system, std::vector<GeoLib::Point> &vec_pt);
    ///
    void rotate(IElement &e, const CoordinateSystem &coordinate_system, std::vector<GeoLib::Point> &vec_pt);
    // x=Rx' where x is original coordinates and x' is local coordinates
    void getRotationMatrixToOriginal(const IElement &e, const CoordinateSystem &coordinate_system, const std::vector<GeoLib::Point> &vec_pt);

private:
    const IMesh* _msh;
    std::vector<GeoLib::Point> _point_vec;
    GeoLib::Point _pt_translate;
    MathLib::LocalMatrix _matR2original;
    MathLib::LocalMatrix _matR2local;
    bool _is_R2orig_set;
};

}
