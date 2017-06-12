/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TestMatrix.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include <gtest/gtest.h>

#include <Eigen>
#include "MathLib/LinAlg/Dense/Matrix.h"
//#include "MathLib/LinAlg/LinearEquations/LinearEquationsFactory.h"
#include "MathLib/LinAlg/LinearEquation/DenseLinearEquation.h"
//#include "MathLib/LinAlg/LinearEquations/LisInterface.h"
//#include "MathLib/Nonlinear/NewtonRaphson.h"
//#include "MathLib/Nonlinear/Picard.h"
//#include "MathLib/Nonlinear/BisectionMethod.h"
//#include "MathLib/Function/IFunction.h"
#include "TestUtil.h"

using namespace MathLib;


// --------------------------------------------------------------------------------
TEST(Math, Matrix_transposeAndMultiply1)
{
    MathLib::Matrix<double> matA(1,2);
    MathLib::Matrix<double> matB(1,2);
    MathLib::Matrix<double> matC(2,2);

    matA(0,0) = 1.0;
    matA(0,1) = 2.0;
    matB(0,0) = 3.0;
    matB(0,1) = 4.0;
    matC = .0;
    matA.transposeAndMultiply(matB, matC);

    ASSERT_EQ(matC(0,0), 3.0);
    ASSERT_EQ(matC(0,1), 4.0);
    ASSERT_EQ(matC(1,0), 6.0);
    ASSERT_EQ(matC(1,1), 8.0);
}

TEST(Math, Matrix_transposeAndMultiply2)
{
    MathLib::Matrix<double> matA(1,2);
    double vec[2] = {};
    MathLib::Matrix<double> matB(2,2);
    MathLib::Matrix<double> matC(2,2);

    matA(0,0) = 1.0;
    matA(0,1) = 2.0;
    vec[0] = 1.0;
    vec[1] = -1.0;
    matB(0,0) = 1.0;
    matB(0,1) = 4.0;
    matB(1,0) = 3.0;
    matB(1,1) = 2.0;
    matC = .0;
    matA.transposeAndMultiply(matB, vec, matC);

    ASSERT_EQ(-2., matC(0,0));
    ASSERT_EQ(2., matC(0,1));
    ASSERT_EQ(-4., matC(1,0));
    ASSERT_EQ(4., matC(1,1));
}


TEST(Math, MatrixAddSubMatrix)
{
    //Matrix<double> _m(4,4);
    DenseLinearEquation::LocalMatrix _m(4,4);
    _m(0,0) = 4.0; _m(0,1) = -1.0; _m(0,2) = -2.0; _m(0,3) = -1.0; 
    _m(1,1) = 4.0; _m(1,2) = -1.0; _m(1,3) = -2.0;
    _m(2,2) = 4.0; _m(2,3) = -1.0;
    _m(3,3) = 4.0;
    for (size_t i=0; i<4; i++)
        for (size_t j=0; j<i; j++) _m(i,j) = _m(j,i);

    DenseLinearEquation eqs;
    eqs.create(10);

    for (size_t i=0; i<3; i++) {
        std::vector<size_t> pos(4);
        for (size_t j=0; j<4; j++)
            pos[j] = i*3 + j;
        eqs.addAsub(pos, _m);
    }

    double expected_m[] = {
         4, -1, -2, -1,  0,  0,  0,  0,  0,  0,
        -1,  4, -1, -2,  0,  0,  0,  0,  0,  0,
        -2, -1,  4, -1,  0,  0,  0,  0,  0,  0,
        -1, -2, -1,  8, -1, -2, -1,  0,  0,  0,
         0,  0,  0, -1,  4, -1, -2,  0,  0,  0,
         0,  0,  0, -2, -1,  4, -1,  0,  0,  0,
         0,  0,  0, -1, -2, -1,  8, -1, -2, -1,
         0,  0,  0,  0,  0,  0, -1,  4, -1, -2,
         0,  0,  0,  0,  0,  0, -2, -1,  4, -1,
         0,  0,  0,  0,  0,  0, -1, -2, -1,  4
    };

    ASSERT_DOUBLE_ARRAY_EQ(expected_m, eqs.getA()->getData(), 100);
}
