/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LinearInterpolation.h
 *
 * Created on 2010-09-07 by Thomas Fischer
 */

#ifndef LINEARINTERPOLATION_H_
#define LINEARINTERPOLATION_H_

#include <vector>
#include <limits>

namespace MathLib {

/**
 * \brief Linear interpolation method
 */
class LinearInterpolation
{
public:
    /**
     *
     * @param supporting_points
     * @param values_at_supp_pnts
     */
    LinearInterpolation(const std::vector<double>& supporting_points, const std::vector<double>& values_at_supp_pnts);

    /**
     *
     * @param supporting_points
     * @param values_at_supp_pnts
     * @param points_to_interpolate
     * @param values_at_interpol_pnts
     */
    LinearInterpolation(const std::vector<double>& supporting_points, const std::vector<double>& values_at_supp_pnts, const std::vector<double>& points_to_interpolate, std::vector<double>& values_at_interpol_pnts);

    /**
     * Copy constructor
     * @param src
     */
    LinearInterpolation(const LinearInterpolation &src);

    /**
     *
     */
    virtual ~LinearInterpolation();

    /**
     *
     * @param pnt_to_interpolate
	 * @return the slope at pnt_to_interpolate
     * @return
     */
    double getValue ( double pnt_to_interpolate );
	double getSlope ( double pnt_to_interpolate );

private:
    const std::vector<double> _supporting_points;
    const std::vector<double> _values_at_supp_pnts;
};

} // end namespace MathLib

#endif /* LINEARINTERPOLATION_H_ */
