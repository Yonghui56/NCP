/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file AxisAlignedBoundingBox.h
 *
 * Created on 2010-04-22 by Thomas Fischer
 */

#ifndef AXISALIGNEDBOUNDINGBOX_H_
#define AXISALIGNEDBOUNDINGBOX_H_

#include <limits>
#include <vector>

#include "Point.h"

namespace GeoLib
{
/**
 *
 * \ingroup GeoLib
 *
 * \brief Class AABB is a bounding box around a given geometric entity
 * */
class AxisAlignedBoundingBox
{
public:
    /**
     * construction of object, initialization the axis aligned bounding box
     * */
    AxisAlignedBoundingBox ();

    /**
     * construction of object using vector of points
     * */
    AxisAlignedBoundingBox ( const std::vector<GeoLib::Point*> *points );

    void update (GeoLib::Point const & pnt);
    /**
     * update axis aligned bounding box
     */
    void update (double x, double y, double z);

    /**
     * update axis aligned bounding box
     */
    void update (const double *pnt)
    {
        update (pnt[0], pnt[1], pnt[2]);
    }

    /**
     * check if point is in the axis aligned bounding box
     * (employing containsPoint (double x, double y, double z))
     */
    bool containsPoint (GeoLib::Point const & pnt, double eps = std::numeric_limits<double>::epsilon()) const;

	/**
	 * wrapper for GeoLib::Point
	 */
    bool containsPoint (const double *pnt, double eps = std::numeric_limits<double>::epsilon()) const;

    /**
     * check if point described by its coordinates x, y, z is in
     * the axis aligned bounding box
     */
    bool containsPoint (double x, double y, double z, double eps = std::numeric_limits<double>::epsilon()) const;

    GeoLib::Point getMinPoint () const { return _min_pnt; }
    GeoLib::Point getMaxPoint () const { return _max_pnt; }

private:
    GeoLib::Point _min_pnt;
    GeoLib::Point _max_pnt;
};

} // end namespace

#endif /* AXISALIGNEDBOUNDINGBOX_H_ */
