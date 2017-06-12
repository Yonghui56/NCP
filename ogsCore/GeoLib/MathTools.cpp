/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file MathTools.cpp
 *
 * Created on 2012-xx-xx by Thomas Fischer
 */

#include "MathTools.h"
#include "MathLib/MathTools.h"

namespace GeoLib
{

double sqrNrm2 (const GeoLib::Point* p0)
{
    return MathLib::scpr (p0->getData(), p0->getData(), 3);
}

double sqrDist (const GeoLib::Point* p0, const GeoLib::Point* p1)
{
    const double v[3] = {(*p1)[0] - (*p0)[0], (*p1)[1] - (*p0)[1], (*p1)[2] - (*p0)[2]};
    return MathLib::scpr (v, v, 3);
}

bool checkDistance(GeoLib::Point const &p0, GeoLib::Point const &p1, double squaredDistance)
{
    return (sqrDist(&p0, &p1) < squaredDistance);
}


double triangleArea(GeoLib::Point const &p0, GeoLib::Point const &p1, GeoLib::Point const &p2)
{
    double A = 0.5*(p0[0]*(p1[1]-p2[1])+p1[0]*(p2[1]-p0[1])+p2[0]*(p0[1]-p1[1]));
    return A;
}

}
