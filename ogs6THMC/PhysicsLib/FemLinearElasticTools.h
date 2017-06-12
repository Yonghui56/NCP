/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemLinearElasticToolss.h
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */


#pragma once

#include <string>
#include "MathLib/DataType.h"

template <class T_MATRIX>
inline void setNu_Matrix_byPoint(const size_t dim, const size_t nnodes, const T_MATRIX &N, T_MATRIX &matN)
{
    matN *= .0;
    for (size_t in=0; in<nnodes; in++) {
        const size_t offset = in*dim;
        for (size_t k=0; k<dim; k++)
            matN(k,offset+k) = N(in);
    }
}

template <class T_MATRIX>
inline void setNu_Matrix_byComponent(const size_t dim, const size_t nnodes, const T_MATRIX &N, T_MATRIX &matN)
{
    matN *= .0;
    for (size_t i_dim = 0; i_dim<dim; i_dim++) {
        const size_t offset = i_dim*nnodes;
        for (size_t in=0; in<nnodes; in++) {
            matN(i_dim, offset + in) = N(in);
        }
    }
}

template <class T_MATRIX>
inline void setB_Matrix_byPoint(const size_t dim, const size_t nnodes, const T_MATRIX &dN, T_MATRIX &matB)
{
    matB *= .0;
    if (dim==2) {
        for (size_t in=0; in<nnodes; in++) {
            const size_t offset = in*dim;
            const double dshape_dx = dN(0, in);
            const double dshape_dy = dN(1, in);
            matB(0,offset) = dshape_dx;
            matB(1,offset+1) = dshape_dy;
            matB(3,offset) = dshape_dy;
            matB(3,offset+1) = dshape_dx;
        }
    } else if (dim==3) {
        for (size_t in=0; in<nnodes; in++) {
            const size_t offset = in*dim;
            const double dshape_dx = dN(0, in);
            const double dshape_dy = dN(1, in);
            const double dshape_dz = dN(2, in);
            matB(0,offset) = dshape_dx;
            matB(1,offset+1) = dshape_dy;
            matB(2,offset+2) = dshape_dz;
            matB(3,offset) = dshape_dy;
            matB(3,offset+1) = dshape_dx;
            matB(4,offset) = dshape_dz;
            matB(4,offset+2) = dshape_dx;
            matB(5,offset+1) = dshape_dz;
            matB(5,offset+2) = dshape_dy;
        }
    }
}

template <class T_MATRIX>
inline void setB_Matrix_byComponent(const size_t dim, const size_t nnodes, const T_MATRIX &dN, T_MATRIX &matB)
{
    matB *= .0;

    if (dim==2) {
        for (size_t i_dim=0; i_dim<dim; i_dim++) {
            for (size_t in=0; in<nnodes; in++) {
                const double dshape_dx = dN(0, in);
                const double dshape_dy = dN(1, in);
                matB(0,in) = dshape_dx;
                matB(1,nnodes+in) = dshape_dy;
                matB(3,in) = dshape_dy;
                matB(3,nnodes+in) = dshape_dx;
            }
        }
    } else if (dim==3) {
        for (size_t i_dim=0; i_dim<dim; i_dim++) {
            for (size_t in=0; in<nnodes; in++) {
                const double dshape_dx = dN(0, in);
                const double dshape_dy = dN(1, in);
                const double dshape_dz = dN(2, in);
                matB(0,in) = dshape_dx;
                matB(1,nnodes+in) = dshape_dy;
                matB(2,nnodes*2+in) = dshape_dz;
                matB(3,in) = dshape_dy;
                matB(3,nnodes+in) = dshape_dx;
                matB(4,in) = dshape_dz;
                matB(4,nnodes*2+in) = dshape_dx;
                matB(5,nnodes+in) = dshape_dz;
                matB(5,nnodes*2+in) = dshape_dy;
            }
        }
    }
}

template <class T_MATRIX>
inline void setB_Matrix(const size_t dim, double dshape_dx, double dshape_dy, double dshape_dz, T_MATRIX &B_matrix)
{
    switch(dim)
    {
    case 2:
        // B_11, dN/dx
        (*B_matrix)(0,0) = dshape_dx;
        // B_12, 0.0
        (*B_matrix)(0,1) = 0.0;
        // B_21, 0.0
        (*B_matrix)(1,0) = 0.0;
        // B_22, dN/dy
        (*B_matrix)(1,1) = dshape_dy;
        // B_31, 0.0
        (*B_matrix)(2,0) = 0.0;
        // B_32, 0.0
        (*B_matrix)(2,1) = 0.0;
        // B_41, dN/dy
        (*B_matrix)(3,0) = dshape_dy;
        // B_42, dN/dx
        (*B_matrix)(3,1) = dshape_dx;

        break;
    case 3:
        // B_11, dN/dx
        (*B_matrix)(0,0) = dshape_dx;
        // B_22, dN/dy
        (*B_matrix)(1,1) = dshape_dy;
        // B_33, dN/dz
        (*B_matrix)(2,2) = dshape_dz;
        //
        // B_41, dN/dy
        (*B_matrix)(3,0) = dshape_dy;
        // B_42, dN/dx
        (*B_matrix)(3,1) = dshape_dx;
        //
        // B_51, dN/dz
        (*B_matrix)(4,0) = dshape_dz;
        // B_53, dN/dx
        (*B_matrix)(4,2) = dshape_dx;
        //
        // B_62, dN/dz
        (*B_matrix)(5,1) = dshape_dz;
        // B_63, dN/dy
        (*B_matrix)(5,2) = dshape_dy;

        break;
    }
}

template <class T_MATRIX>
inline void setB_Matrix4axisymmetry(const size_t dim, double radius, double dshape_dx, double dshape_dy, T_MATRIX &B_matrix)
{
    // B_11, dN/dx
    (*B_matrix)(0,0) = dshape_dx;
    // B_12, 0.0
    (*B_matrix)(0,1) = 0.0;
    // B_21, N/r
    (*B_matrix)(1,0) = dshape_dx / radius;
    // B_22, 0.0
    (*B_matrix)(1,1) = 0.0;
    // B_31, 0.0
    (*B_matrix)(2,0) = 0.0;
    // B_32, dN/dz
    (*B_matrix)(2,1) = dshape_dy;
    // B_41, dN/dy
    (*B_matrix)(3,0) = dshape_dy;
    // B_42, dN/dx
    (*B_matrix)(3,1) = dshape_dx;
}



