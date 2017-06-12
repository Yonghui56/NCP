/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TRI3CONST.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <algorithm>
#include "MeshLib/Core/Element.h"
#include "TemplateFeBase.h"

namespace FemLib
{

/**
 * \brief TRI3CONST finite element class
 *
 * TRI3CONST is a finite element class for specific case of linear triangle elements with
 * element constant parameters. Efficient integration can be expected due to use of analytical
 * integration methods.
 * 
 * Remark:
 * - Given coefficients are assumed to be constant over a element.
 * - Axisymmetric model is not supported!
 */    
class TRI3CONST : public TemplateFeBase<3>
{
private:
    MeshLib::Triangle *_ele;
    double a[3], b[3], c[3];
    double A;
    FemIntegrationAnalytical _integration;
    MathLib::LocalMatrix _shape, _dshape;
    double x_cp[3];

    void computeBasisFunction(const double *x,  double *shape);
    void computeGradBasisFunction(const double *x,  MathLib::LocalMatrix &mat);
public:
    explicit TRI3CONST(MeshLib::IMesh* msh)
    : TemplateFeBase<3>(msh), _ele(NULL), A(.0), _shape(1,3), _dshape(2,3) {};

    /// initialize object for given mesh elements
    void configure( MeshLib::IElement &e );

    /// 
    void computeBasisFunctions(const double *x);
    /// compute real coordinates from the given position in reference coordinates
    virtual void getRealCoordinates(double* x_real);
    MathLib::LocalMatrix* getBasisFunction();
    MathLib::LocalMatrix* getGradBasisFunction();
    virtual double getDetJ() const {return 1.0;}


    /// make interpolation from nodal values
    double interpolate(double *x, double *nodal_values);

//    /// compute an matrix M = Int{W^T F N} dV
//    void integrateWxN(MathLib::SpatialFunctionScalar* f, LocalMatrix &mat);
//
//    /// compute an matrix M = Int{W^T F dN} dV
//    void integrateWxDN(MathLib::SpatialFunctionVector* f, LocalMatrix &mat);
//
//    /// compute an matrix M = Int{dW^T F dN} dV
//    void integrateDWxDN(MathLib::SpatialFunctionScalar *f, LocalMatrix &mat);

    /// compute an matrix M = Int{W^T F N} dV
    virtual void integrateWxN(size_t igp, MathLib::LocalMatrix &f, MathLib::LocalMatrix &mat);

    /// compute an matrix M = Int{W^T F dN} dV
    virtual void integrateWxDN(size_t igp, MathLib::LocalMatrix &f, MathLib::LocalMatrix &mat);

    /// compute an matrix M = Int{dW^T F dN} dV
    virtual void integrateDWxDN(size_t igp, MathLib::LocalMatrix &f, MathLib::LocalMatrix &mat);

    /// get the integration method
    IFemNumericalIntegration* getIntegrationMethod() const {return (IFemNumericalIntegration*)&_integration;};

    void extrapolate(const std::vector<MathLib::LocalVector> &gp_values, std::vector<MathLib::LocalVector> &nodal_values);
};

}
