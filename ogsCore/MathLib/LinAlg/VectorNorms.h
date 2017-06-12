/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file VectorNorms.h
 *
 * Created on 2011-06-06 by Thomas Fischer
 */

#ifndef VECTORNORMS_H_
#define VECTORNORMS_H_

#include <cmath>
#include <algorithm>
#include "MathLib/MathTools.h"

namespace MathLib {

inline double normEuklid (double const * const vec, size_t n)
{
    return sqrt (scpr (vec, vec, n));
}

template<class T>
inline double norm_p(T &v, size_t n, int p)
{
    double s = .0;
    for (size_t i=0; i<n; i++)
        s += pow(fabs(v[i]),p);
    s = pow(s, 1./p);
    return s;
};

template<class T>
inline double norm1(T &v, size_t n)
{
    double s = .0;
    for (size_t i=0; i<n; i++)
        s += fabs(v[i]);
    return s;
};

template<class T>
inline double norm_max(T &v, size_t n)
{
    double u_max = .0;
    for (size_t i=0; i<n; i++)
        u_max = std::max(u_max, fabs(v[i]));
    return u_max;
};


} // end namespace MathLib

#endif /* VECTORNORMS_H_ */
