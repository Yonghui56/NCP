/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file MathTools.cpp
 *
 * Created on 2010-01-13 by Thomas Fischer
 */

#include "MathTools.h"

namespace MathLib {

#ifdef _OPENMP
double scpr(double const * const v, double const * const w, OPENMP_LOOP_TYPE n)
{
    long double res (v[0]*w[0]);
    OPENMP_LOOP_TYPE k;

    #pragma omp parallel for reduction (+:res)
    for (k = 1; k<n; k++) {
        res += v[k] * w[k];
    }

    return res;
}
#endif


void crossProd(const double u[3], const double v[3], double r[3])
{
    r[0] = u[1] * v[2] - u[2] * v[1];
    r[1] = u[2] * v[0] - u[0] * v[2];
    r[2] = u[0] * v[1] - u[1] * v[0];
}

void normalizeVector(const double* u, size_t n, double* r)
{
    double nrm(u[0] * u[0]);
    for(size_t i = 1; i < n; i++)
        nrm += u[i] * u[i];
    double sqrt_nrm (sqrt(nrm));
    for(size_t i = 0; i < n; i++)
        r[i] = u[i] / sqrt_nrm;
}

void normalizeVector(double* u, size_t n)
{
    normalizeVector(u, n, u);
}

double calcProjPntToLineAndDists(const double p[3], const double a[3],
        const double b[3], double &lambda, double &d0)
{
    // g (lambda) = a + lambda v, v = b-a
    double v[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    // orthogonal projection: (g(lambda)-p) * v = 0 => in order to compute lambda we define a help vector u
    double u[3] = {p[0] - a[0], p[1] - a[1], p[2] - a[2]};
    lambda = scpr (u, v, 3) / scpr (v, v, 3);

    // compute projected point
    double proj_pnt[3];
    for (size_t k(0); k<3; k++) proj_pnt[k] = a[k] + lambda * v[k];

    d0 = sqrt (sqrDist (proj_pnt, a));

    return sqrt (sqrDist (p, proj_pnt));
}

double sqrDist(const double* p0, const double* p1)
{
    const double v[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
    return scpr (v, v, 3);
}

float normalize(float min, float max, float val)
{
    return ((val-min)/static_cast<float>(max-min));
}

double getAngle (const double p0[3], const double p1[3], const double p2[3])
{
    const double v0[3] = {p0[0]-p1[0], p0[1]-p1[1], p0[2]-p1[2]};
    const double v1[3] = {p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]};

    // apply Cauchy Schwarz inequality
    return acos (scpr (v0,v1,3) / (sqrt(scpr(v0,v0,3)) * sqrt(scpr (v1,v1,3))));
}

double erf(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = std::abs(x);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return sign*y;
}

double erfc(double x)
{
    return 1.0 - erf(x);
}

} // namespace
