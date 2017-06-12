/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Plane.h
 *
 * Created on xxxx-xx-xx by Thomas Fisher
 */

#pragma once

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/MathTools.h"
#include "Point.h"

namespace GeoLib
{

/**
 * 
 */
class Plane
{
public:
    /**
     * create a rotation matrix from a given unit normal vector \f$ \mathbf{n}=(n_1,n_2,n_3)^T \f$ of the plane.
     * The local basis vectors \f$\mathbf{\bar{e}} \f$ can be calculated as
     * \f[
     *    \bar{e}_1 = \left ( \frac{n_3}{n_1}\frac{1}{\sqrt{1+(\frac{n_3}{n_1})^2}}, 0, -\frac{1}{\sqrt{1+(\frac{n_3}{n_1})^2}} \right )^T
     * \f]
     * \f[
     *    \bar{e}_2 = \frac{1}{D} \left ( n_2 \bar{e}_{13}, n_3 \bar{e}_{11} - n_1 \bar{e}_{13}, -n_2 \bar{e}_{11} \right )^T
     * \f]
     * \f[
     *    \bar{e}_3 = \left ( n_1, n_2, n_3 \right )^T
     * \f]
     * with
     * \f[
     *    D = \sqrt{ (n_2 \bar{e}_{13})^2 + (n_3 \bar{e}_{11} - n_1 \bar{e}_{13})^2 + (-n_2 \bar{e}_{11})^2} 
     * \f]
     */
    Plane(const Point &pt1, const Point &pt2, const Point &pt3)
    {
        calculateNormalVector(pt1, pt2, pt3);
    }

    Plane(double azimuth, double dip)
    {
        calculateNormalVector(azimuth, dip);
    }

    Plane(double *angle_x1x2)
    {
        calculateNormalVector(angle_x1x2);
    }

    const double* getNormalVector() const
    {
        return _n;
    }

    const void getLocalCoordinateSystem()
    {

    }

private:
    double _n[3];

    /**
     * a1 = p1-p0 
     * a2 = p2-p0 
     * n = a1 x a2 / |a1 x a2|
     */
    void calculateNormalVector(const Point &pt1, const Point &pt2, const Point &pt3)
    {
        double a1[3];
        double a2[3];
        for (size_t i=0; i<3; i++)
            a1[i] = pt2[i] - pt1[0];
        for (size_t i=0; i<3; i++)
            a2[i] = pt3[i] - pt1[0];
        MathLib::crossProd(a1,a2,_n);
        MathLib::normalizeVector(_n,3);
    }

    void calculateNormalVector(double azimuth, double dip)
    {
        _n[0] = sin(azimuth)*sin(dip);
        _n[1] = cos(azimuth)*sin(dip);
        _n[2] = cos(dip);
    };

    void calculateNormalVector(double *angle)
    {
        _n[0] = sin(angle[0])*cos(angle[1]);
        _n[1] = -sin(angle[1]);
        _n[2] = cos(angle[0])*cos(angle[1]);
    };

    void setRotationMatrix( double u[3] ) 
    {
        MathLib::Matrix<double> _matVec(3,3);
        const double ax = u[0];
        const double ay = u[1];
        const double az = u[2];

        double e1[3]={0.0},e2[3]={0.0},e3[3]={0.0};

        if (ax==0.0) {
            e1[0] = 1.0;
            e1[1] = 0.0;
            e1[2] = 0.0;

            e2[0] = 0.0;
            e2[1] = az;
            e2[2] = -ay;

            e3[0] = 0.0;
            e3[1] = ay;
            e3[2] = az;
        } else {
            double Disk = 1.0/sqrt(1.0+((az*az)/(ax*ax)));

            e1[0] = (az/ax)*Disk;
            e1[1] = 0.0;
            e1[2] = -Disk;

            double t1    = ay*e1[2];
            double t2    = az*e1[0]-ax*e1[2];
            double t3    = -ay*e1[0];

            Disk  = t1*t1;
            Disk += t2*t2;
            Disk += t3*t3;

            Disk = 1.0/sqrt(Disk);

            e2[0] = Disk*t1;
            e2[1] = Disk*t2;
            e2[2] = Disk*t3;

            e3[0] = u[0];
            e3[1] = u[1];
            e3[2] = u[2];
        }

        _matVec(0,0) = e1[0];
        _matVec(0,1) = e1[1];
        _matVec(0,2) = e1[2];

        _matVec(1,0) = e2[0];
        _matVec(1,1) = e2[1];
        _matVec(1,2) = e2[2];

        _matVec(2,0) = e3[0];
        _matVec(2,1) = e3[1];
        _matVec(2,2) = e3[2];
    }
};

}
