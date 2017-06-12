/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TXFunctionLinear.h
 *
 * Created on 2012-11-14 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include "GeoLib/Surface.h"
#include "TXPosition.h"
#include "ITXFunction.h"

namespace NumLib
{

/**
 * \brief Constant value
 */
class TXFunctionLinear : public ITXFunction
{
public:
    TXFunctionLinear(const GeoLib::Surface* sfc, const std::vector<std::pair<size_t, double> > &list_pt_val)
    : _sfc(sfc), _list_pt_val(list_pt_val)
    {
        ITXFunction::isConst(true);
        ITXFunction::isTemporallyConst(true);
        ITXFunction::isSpatiallyConst(false);
    };

    virtual ~TXFunctionLinear() {};

#if 0
    virtual void eval(const TXPosition x, double &val) const
    {
        //----------------------------------------------------------------------
        // Interpolation of polygon values to nodes_on_sfc

        // find a triangle in which the given position locates
        std::vector<size_t> list_tri_vertex_id;
        findTriangle(x.getSpace(), list_tri_vertex_id);
        // get values at the vertex
        std::vector<double> list_tri_vertex_value;
        // interpolate within the triangle
        double unit[3];
        for (size_t l = 0; l < 2; l++)
            unit[l] /= Area1;
        double NTri[3];
        ShapeFunctionTri(NTri, unit);
        for (size_t l = 0; l < 3; l++)
            val += list_tri_vertex_value[l] * NTri[l];
    }
#endif

#if 0
    bool findTriangle(const double* pn, std::vector<size_t> &list_tri_vertex) const
    {
        double Area1, Area2;
        double Tol = 0.1;
        double gC[3], p1[3], p2[3], vn[3], unit[3], NTri[3];
        //
        GeoLib::Polyline* m_polyline = NULL;
        const GeoLib::Surface* m_surface = _sfc;

        std::vector<GeoLib::Polyline*>::iterator p =
                m_surface->polyline_of_surface_vector.begin();

        // nodes close to first polyline
        p = m_surface->polyline_of_surface_vector.begin();
        while (p != m_surface->polyline_of_surface_vector.end())
        {
            m_polyline = *p;
            // Gravity center of this polygon
            for (size_t i = 0; i < 3; i++)
                gC[i] = 0.0;
            vn[2] = 0.0;
            const size_t nPointsPly = m_polyline->point_vector.size();
            for (size_t i = 0; i < nPointsPly; i++)
            {
                gC[0] += m_polyline->point_vector[i]->x;
                gC[1] += m_polyline->point_vector[i]->y;
                gC[2] += m_polyline->point_vector[i]->z;
                vn[2] += m_polyline->point_vector[i]->getPropert();
            }
            for (size_t i = 0; i < 3; i++)
                gC[i] /= (double) nPointsPly;
            // BC value at center is an average of all point values of polygon
            vn[2] /= (double) nPointsPly;
            // Area of this polygon by the gravity center
            for (size_t i = 0; i < nPointsPly; i++)
            {
                p1[0] = m_polyline->point_vector[i]->x;
                p1[1] = m_polyline->point_vector[i]->y;
                p1[2] = m_polyline->point_vector[i]->z;
                size_t k = i + 1;
                if (i == nPointsPly - 1)
                    k = 0;
                p2[0] = m_polyline->point_vector[k]->x;
                p2[1] = m_polyline->point_vector[k]->y;
                p2[2] = m_polyline->point_vector[k]->z;
                vn[0] = m_polyline->point_vector[i]->getPropert();
                vn[1] = m_polyline->point_vector[k]->getPropert();

                Area1 = fabs(ComputeDetTri(p1, gC, p2));

                Area2 = 0.0;
                // Check if pn is in the triangle by points (p1, gC, p2)
                Area2 = fabs(ComputeDetTri(p2, gC, pn));
                unit[0] = fabs(ComputeDetTri(gC, p1, pn));
                unit[1] = fabs(ComputeDetTri(p1, p2, pn));
                Area2 += unit[0] + unit[1];
                if (fabs(Area1 - Area2) < Tol)
                {
                    return true;
                }
            }
            //
            p++;
        } // while
        return false;

    }
#endif

    virtual TXFunctionLinear* clone() const OGS_DECL_OVERRIDE
    {
        return new TXFunctionLinear(_sfc, _list_pt_val);
    }

private:
    const GeoLib::Surface* _sfc;
    std::vector<std::pair<size_t, double> > _list_pt_val;
};

} //end


