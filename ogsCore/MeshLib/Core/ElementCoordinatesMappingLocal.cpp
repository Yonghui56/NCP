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

#include "ElementCoordinatesMappingLocal.h"

#include <limits>
#include <cassert>
#include "MathLib/MathTools.h"

namespace MeshLib
{

ElementCoordinatesMappingLocal::ElementCoordinatesMappingLocal(const IMesh* msh, IElement &e, const CoordinateSystem &coordinate_system)
: _msh(msh), _pt_translate(0,0,0), _is_R2orig_set(false)
{
    assert (e.getDimension() <= coordinate_system.getDimension());

    // set initial coordinates
    for(size_t i = 0; i < e.getNumberOfNodes(); i++)
    {
        const GeoLib::Point *p = _msh->getNodeCoordinatesRef(e.getNodeID(i));
        _point_vec.push_back(GeoLib::Point((*p)[0], (*p)[1], (*p)[2]));
    }

    flip(e, coordinate_system, _point_vec);
    if (e.getDimension() < coordinate_system.getDimension()) {
        translate(_point_vec);
        rotate(e, coordinate_system, _point_vec);
    }
};

void ElementCoordinatesMappingLocal::translate(std::vector<GeoLib::Point> &vec_pt)
{
    _pt_translate =  vec_pt[0];
    for (size_t i=0; i<vec_pt.size(); ++i) {
        vec_pt[i] -= _pt_translate;
    }
}

///
void ElementCoordinatesMappingLocal::flip(IElement &ele, const CoordinateSystem &coordinate_system, std::vector<GeoLib::Point> &vec_pt)
{
    IElement* e = &ele;
    switch(coordinate_system.getType())
    {
    case CoordinateSystemType::Y:
        {
            assert(e->getDimension()==1);
            for(size_t i = 0; i < vec_pt.size(); i++)
            {
                double tmp_x = vec_pt[i][0];
                vec_pt[i][0] = vec_pt[i][1];
                vec_pt[i][1] = tmp_x;
            }
        }
        break;
    case CoordinateSystemType::Z:
        {
            assert(e->getDimension()==1);
            for(size_t i = 0; i < vec_pt.size(); i++)
            {
                double tmp_x = vec_pt[i][0];
                vec_pt[i][0] = vec_pt[i][2];
                vec_pt[i][2] = tmp_x;
            }
        }
        break;
    case CoordinateSystemType::XZ:
        {
            assert(e->getDimension()<3);
            for(size_t i = 0; i < vec_pt.size(); i++)
            {
                double tmp_y = vec_pt[i][1];
                vec_pt[i][1] = vec_pt[i][2];
                vec_pt[i][2] = tmp_y;
            }
        }
        break;
    case CoordinateSystemType::YZ:
        {
            assert(e->getDimension()<3);
            for(size_t i = 0; i < vec_pt.size(); i++)
            {
                double tmp_x = vec_pt[i][0];
                vec_pt[i][0] = vec_pt[i][1];
                vec_pt[i][1] = vec_pt[i][2];
                vec_pt[i][2] = tmp_x;
            }
        }
        break;
    default:
        break;
    }
}

///
void ElementCoordinatesMappingLocal::rotate(IElement &ele, const CoordinateSystem &coordinate_system, std::vector<GeoLib::Point> &vec_pt)
{
    const size_t global_dim = coordinate_system.getDimension();
    IElement* e = &ele;

    if (!_is_R2orig_set) {
        _matR2original = MathLib::LocalMatrix::Zero(global_dim, global_dim);
        getRotationMatrixToOriginal(*e, coordinate_system, vec_pt);
        _matR2local = _matR2original.transpose();
        _is_R2orig_set = true;
    }

    double const* const coords_node_0 (vec_pt[0].getData());
    MathLib::LocalVector dx = MathLib::LocalVector::Zero(global_dim);
    MathLib::LocalVector x_new = MathLib::LocalVector::Zero(3);
    for(size_t i = 0; i < e->getNumberOfNodes(); i++)
    {
        double const* const coords_node_i (vec_pt[i].getData());
        for (size_t j=0; j<global_dim; j++)
            dx[j] = (coords_node_i[j] - coords_node_0[j]);

        x_new.head(global_dim) = _matR2local * dx;
        _point_vec[i] = GeoLib::Point(x_new.data());
    }
};

// x=Rx' where x is original coordinates and x' is local coordinates
void ElementCoordinatesMappingLocal::getRotationMatrixToOriginal(const IElement &ele, const CoordinateSystem &coordinate_system, const std::vector<GeoLib::Point> &vec_pt)
{
    const size_t global_dim = coordinate_system.getDimension();
    double xx[3];
    double yy[3];
    double zz[3];
    const IElement* e = &ele;

    if (global_dim == e->getDimension()) {
        for (size_t i=0; i<global_dim; i++)
            _matR2original(i,i) = 1.;
    } else if (coordinate_system.getDimension() == 2 && e->getDimension() == 1) {
//    } else if (coordinate_system.getType() == CoordinateSystemType::XY && e->getDimension() == 1) {
        double const* const pnt0(vec_pt[0].getData());
        double const* const pnt1(vec_pt[1].getData());
        xx[0] = pnt1[0] - pnt0[0];
        xx[1] = pnt1[1] - pnt0[1];
        MathLib::normalizeVector(xx, 2);
        double cos_theta = xx[0];
        double sin_theta = xx[1];
        _matR2original(0,0) = _matR2original(1,1) = cos_theta;
        _matR2original(0,1) = - sin_theta;
        _matR2original(1,0) = sin_theta;
    } else if (coordinate_system.getType() == CoordinateSystemType::XYZ && e->getDimension() == 2) {
        // x"_vec
        //            xx[0] = nodes[1]->X() - nodes[0]->X();
        //            xx[1] = nodes[1]->Y() - nodes[0]->Y();
        //            xx[2] = nodes[1]->Z() - nodes[0]->Z();
        double const* const pnt0(vec_pt[0].getData());
        double const* const pnt1(vec_pt[1].getData());
        xx[0] = pnt1[0] - pnt0[0];
        xx[1] = pnt1[1] - pnt0[1];
        xx[2] = pnt1[2] - pnt0[2];
        MathLib::normalizeVector(xx, 3);
        // a vector on the plane
        //            yy[0] = nodes[2]->X() - nodes[1]->X();
        //            yy[1] = nodes[2]->Y() - nodes[1]->Y();
        //            yy[2] = nodes[2]->Z() - nodes[1]->Z();
        double const* const pnt2(vec_pt[2].getData());
        yy[0] = pnt2[0] - pnt1[0];
        yy[1] = pnt2[1] - pnt1[1];
        yy[2] = pnt2[2] - pnt1[2];
        // z"_vec. off plane
        MathLib::crossProd(xx, yy, zz);
        MathLib::normalizeVector(zz, 3);
        // y"_vec
        MathLib::crossProd(zz, xx, yy);
        MathLib::normalizeVector(yy, 3);

        for (size_t i=0; i<global_dim; ++i) {
            _matR2original(i, 0) = xx[i];
            _matR2original(i, 1) = yy[i];
            _matR2original(i, 2) = zz[i];
        }
    } else if (global_dim == 3 && e->getDimension() == 1) {
        // x"_vec
        double const* const pnt0(vec_pt[0].getData());
        double const* const pnt1(vec_pt[1].getData());
        xx[0] = pnt1[0] - pnt0[0];
        xx[1] = pnt1[1] - pnt0[1];
        xx[2] = pnt1[2] - pnt0[2];
        MathLib::normalizeVector(xx, 3);
        // an arbitrary vector
        for (size_t i = 0; i < 3; i++)
            yy[i] = 0.0;

        if (fabs(xx[0]) > 0.0 && fabs(xx[1]) + fabs(xx[2]) < std::numeric_limits<double>::epsilon()) {
            yy[2] = 1.0;
        } else if (fabs(xx[1]) > 0.0 && fabs(xx[0]) + fabs(xx[2]) < std::numeric_limits<double>::epsilon()) {
            yy[0] = 1.0;
        } else if (fabs(xx[2]) > 0.0 && fabs(xx[0]) + fabs(xx[1]) < std::numeric_limits<double>::epsilon()) {
            yy[1] = 1.0;
        } else {
            for (size_t i = 0; i < 3; i++) {
                if (fabs(xx[i]) > 0.0) {
                    yy[i] = -xx[i];
                    break;
                }
            }
        }
        // z"_vec
        MathLib::crossProd(xx, yy, zz);
        MathLib::normalizeVector(zz, 3);
        // y"_vec
        MathLib::crossProd(zz, xx, yy);
        MathLib::normalizeVector(yy, 3);

        for (size_t i=0; i<global_dim; ++i) {
            _matR2original(i, 0) = xx[i];
            _matR2original(i, 1) = yy[i];
            _matR2original(i, 2) = zz[i];
        }
    }

}

}
