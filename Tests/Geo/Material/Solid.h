/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Solid.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

namespace Geo
{

struct Solid
{
    double density;
    double thermal_expansion;
    double poisson_ratio;
    double Youngs_modulus;
};

inline void calculateLameConstant(const double nv, const double E, double &Lambda, double &G, double &K)
{
    Lambda = E * nv / ((1. + nv) * (1. - 2. * nv));
    G = 0.5 * E / (1. + nv);
    K = (3.0 * Lambda + 2.0 * G) / 3.0;
}

template<class T_MATRIX>
inline void setElasticConsitutiveTensor(const size_t Dimension, double Lambda, double G, T_MATRIX &D_e)
{
    //D_e *= 0.0;
    D_e(0,0) = Lambda + 2 * G;
    D_e(0,1) = Lambda;
    D_e(0,2) = Lambda;

    D_e(1,0) = Lambda;
    D_e(1,1) = Lambda + 2 * G;
    D_e(1,2) = Lambda;

    D_e(2,0) = Lambda;
    D_e(2,1) = Lambda;
    D_e(2,2) = Lambda + 2 * G;

    D_e(3,3) = G;
    //Plane stress
    // plane stress, only for test
    //(*D_e)(0,0) = (1.0-Mu)*Lambda + 2 * G;
    //(*D_e)(0,1) = Lambda;

    //(*D_e)(1,0) = Lambda;
    // (*D_e)(1,1) = (1.0-Mu)*Lambda + 2 * G;
    // (*D_e)(3,3) = G;

    if(Dimension == 3)
    {
        D_e(4,4) = G;
        D_e(5,5) = G;
    }
}


}
