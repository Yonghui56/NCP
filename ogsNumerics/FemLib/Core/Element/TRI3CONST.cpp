/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TRI3CONST.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "TRI3CONST.h"

#include "GeoLib/Point.h"

namespace FemLib
{

void TRI3CONST::configure( MeshLib::IElement &e )
{
    _ele = static_cast<MeshLib::Triangle*>(&e);
    double nodes_x[3], nodes_y[3], nodes_z[3];
    // xyz
    for (size_t i=0; i<3; i++) {
        const GeoLib::Point *pt = getMesh()->getNodeCoordinatesRef(_ele->getNodeID(i));
        nodes_x[i] = pt->getData()[0];
        nodes_y[i] = pt->getData()[1];
        nodes_z[i] = pt->getData()[2];
    }
    // area
    A = GeoLib::triangleArea(getMesh()->getNodeCoordinates(_ele->getNodeID(0)),getMesh()->getNodeCoordinates(_ele->getNodeID(1)),getMesh()->getNodeCoordinates(_ele->getNodeID(2)));
    // set a,b,c
    a[0] = 0.5/A*(nodes_x[1]*nodes_y[2]-nodes_x[2]*nodes_y[1]);
    b[0] = 0.5/A*(nodes_y[1]-nodes_y[2]);
    c[0] = 0.5/A*(nodes_x[2]-nodes_x[1]);
    a[1] = 0.5/A*(nodes_x[2]*nodes_y[0]-nodes_x[0]*nodes_y[2]);
    b[1] = 0.5/A*(nodes_y[2]-nodes_y[0]);
    c[1] = 0.5/A*(nodes_x[0]-nodes_x[2]);
    a[2] = 0.5/A*(nodes_x[0]*nodes_y[1]-nodes_x[1]*nodes_y[0]);
    b[2] = 0.5/A*(nodes_y[0]-nodes_y[1]);
    c[2] = 0.5/A*(nodes_x[1]-nodes_x[0]);
    // gravity point
    for (size_t i=0; i<3; i++) x_cp[i] = 0;
    for (size_t i=0; i<3; i++) {
        x_cp[0] += nodes_x[i];
        x_cp[1] += nodes_y[i];
        x_cp[2] += nodes_z[i];
    }
    for (size_t i=0; i<3; i++) 
        x_cp[i] /= 3.0;
}

void TRI3CONST::computeBasisFunctions(const double* /*x*/)
{
    //computeBasisFunction(x, (double*)_shape.getData());
    //computeGradBasisFunction(x, _dshape);
}

void TRI3CONST::getRealCoordinates(double* x_real)
{
    for (size_t i=0; i<3; i++) 
        x_real[i] = x_cp[i];
};

MathLib::LocalMatrix* TRI3CONST::getBasisFunction()
{
    return &_shape;
}

MathLib::LocalMatrix* TRI3CONST::getGradBasisFunction()
{
    return &_dshape;
}

void TRI3CONST::computeBasisFunction(const double *x,  double *shape)
{
    for (size_t i=0; i<3; i++)
        shape[i] = a[i]+b[i]*x[0]+c[i]*x[1];
}

void TRI3CONST::computeGradBasisFunction(const double*,  MathLib::LocalMatrix &dshape)
{
    for (size_t i=0; i<3; i++) {
        dshape(0,i) = b[i];
        dshape(1,i) = c[i];
    }
}

double TRI3CONST::interpolate(double *x, double *nodal_values)
{
    double shape[3] = {};
    computeBasisFunction(x, shape);
    double v = 0;
    for (size_t i=0; i<3; i++)
        v += shape[i]*nodal_values[i];
    return v;
}

///// compute an matrix M = Int{W^T F N} dV
//void TRI3CONST::integrateWxN( MathLib::SpatialFunctionScalar* f, MathLib::LocalMatrix &mat)
//{
//    double v = .0;
//    f->eval(0, v);
//    mat(0,0) = 1.0;
//    mat(0,1) = 0.5;
//    mat(0,2) = 0.5;
//    mat(1,1) = 1.0;
//    mat(1,2) = 0.5;
//    mat(2,2) = 1.0;
//    mat *= v*A/6.0;
//    // make symmetric
//    for (size_t i=0; i<3; i++)
//        for (size_t j=0; j<i; j++)
//            mat(i,j) = mat(j,i);
//}
//
///// compute an matrix M = Int{W^T F dN} dV
//void TRI3CONST::integrateWxDN( MathLib::SpatialFunctionVector* f, MathLib::LocalMatrix &mat)
//{
//    MathLib::Vector v;
//    f->eval(0, v);
//    for (int i=0; i<3; i++)
//        for (int j=0; j<3; j++)
//            mat(i,j) = v[0]*b[j] + v[1]*c[j];
//    mat *= A/3.0;
//}
//
///// compute an matrix M = Int{dW^T F dN} dV
//void TRI3CONST::integrateDWxDN( MathLib::SpatialFunctionScalar *f, MathLib::LocalMatrix &mat)
//{
//    double v;
//    f->eval(0, v);
//    mat(0,0) = b[0]*b[0] + c[0]*c[0];
//    mat(0,1) = b[0]*b[1] + c[0]*c[1];
//    mat(0,2) = b[0]*b[2] + c[0]*c[2];
//    mat(1,1) = b[1]*b[1] + c[1]*c[1];
//    mat(1,2) = b[1]*b[2] + c[1]*c[2];
//    mat(2,2) = b[2]*b[2] + c[2]*c[2];
//    mat *= v*A;
//    // make symmetric
//    for (size_t i=0; i<3; i++)
//        for (size_t j=0; j<i; j++)
//            mat(i,j) = mat(j,i);
//}

/// compute an matrix M = Int{W^T F N} dV
void TRI3CONST::integrateWxN(size_t igp, MathLib::LocalMatrix &v, MathLib::LocalMatrix &mat)
{
    assert(igp==0);

    mat(0,0) = 1.0;
    mat(0,1) = 0.5;
    mat(0,2) = 0.5;
    mat(1,1) = 1.0;
    mat(1,2) = 0.5;
    mat(2,2) = 1.0;
    mat *= v(0,0)*A/6.0;
    // make symmetric
    for (size_t i=0; i<3; i++)
        for (size_t j=0; j<i; j++)
            mat(i,j) = mat(j,i);
}

/// compute an matrix M = Int{W^T F dN} dV
void TRI3CONST::integrateWxDN(size_t igp, MathLib::LocalMatrix &v, MathLib::LocalMatrix &mat)
{
    assert(igp==0);

    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            mat(i,j) = v(0,0)*b[j] + v(1,0)*c[j];
    mat *= A/3.0;
}

/// compute an matrix M = Int{dW^T F dN} dV
void TRI3CONST::integrateDWxDN(size_t igp, MathLib::LocalMatrix &v, MathLib::LocalMatrix &mat)
{
    assert(igp==0);

    mat(0,0) = b[0]*b[0] + c[0]*c[0];
    mat(0,1) = b[0]*b[1] + c[0]*c[1];
    mat(0,2) = b[0]*b[2] + c[0]*c[2];
    mat(1,1) = b[1]*b[1] + c[1]*c[1];
    mat(1,2) = b[1]*b[2] + c[1]*c[2];
    mat(2,2) = b[2]*b[2] + c[2]*c[2];
    mat *= v(0,0)*A;
    // make symmetric
    for (size_t i=0; i<3; i++)
        for (size_t j=0; j<i; j++)
            mat(i,j) = mat(j,i);
}

void TRI3CONST::extrapolate(const std::vector<MathLib::LocalVector> &gp_values, std::vector<MathLib::LocalVector> &nodal_values)
{
    // gp_values are all same
    MathLib::LocalVector v = gp_values[0];
    for (size_t i=0; i<nodal_values.size(); i++) {
        nodal_values[i] = v;
    }
}


} // end namespace

