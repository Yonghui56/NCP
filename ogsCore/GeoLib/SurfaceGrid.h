/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SurfaceGrid.h
 *
 * Created on 2012-02-09 by Thomas Fischer
 */

#ifndef SURFACEGRID_H_
#define SURFACEGRID_H_

#include <vector>

// GEOLIB
#include "AxisAlignedBoundingBox.h"

namespace GeoLib
{

// forward declaration
class Surface;
class Triangle;

class SurfaceGrid : public AxisAlignedBoundingBox
{
public:
	SurfaceGrid(Surface const*const sfc);
	virtual ~SurfaceGrid();

	bool isPntInSurface(const double* pnt, double eps = 0) const;

private:
	double _step_sizes[3];
	double _inverse_step_sizes[3];
	size_t _n_steps[3];
	std::vector<GeoLib::Triangle const*>* _triangles_in_grid_box;
};

} // end namespace GeoLib

#endif /* SURFACEGRID_H_ */
