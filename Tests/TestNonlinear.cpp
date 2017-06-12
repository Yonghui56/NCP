/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TestNonlinear.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include <gtest/gtest.h>

#include <valarray>

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/LinearEquation/DenseLinearEquation.h"
#include "MathLib/LinAlg/LinearEquation/SparseLinearEquation.h"
#ifdef USE_LIS
#include "MathLib/LinAlg/LinearEquation/LisLinearEquation.h"
#endif
#include "MathLib/Nonlinear/NewtonRaphson.h"
#include "MathLib/Nonlinear/Picard.h"
#include "MathLib/Nonlinear/BisectionMethod.h"
#include "MathLib/Nonlinear/NRIterationStepInitializerDummy.h"
#include "NumLib/Function/IFunction.h"
#include "TestUtil.h"

using namespace MathLib;
using namespace NumLib;

template<class T_F, class T_DF>
class TemplateFixedPointFunctonUsingC1Equation
{
    MathLib::DenseLinearEquation solver;
    std::valarray<double> b;
    std::valarray<double> tmp;
    T_F nr;
    T_DF nr_df;
public:
    TemplateFixedPointFunctonUsingC1Equation(size_t n)
    {
        solver.create(n);
        b.resize(n);
        tmp.resize(n);
    }
    // x = g(x) = x - f'(x)^-1 * f(x)
    // f'*g = f'*x - f
    void eval(std::valarray<double> &x, std::valarray<double> &x_new)
    {
        MathLib::Matrix<double> *j = solver.getA();
        solver.reset();
        b = 0;
        tmp = 0;

        nr_df.eval(x, *j);
        nr.eval(x, b);

        b*=-1;
        j->axpy(1.0, &x[0], .0, &tmp[0]);
        b+=tmp;
        for (size_t i=0; i<x.size(); i++)
            solver.setRHS(i, b[i]);

        solver.solve();

        double *u = solver.getX();
        for (size_t i=0; i<x.size(); i++)
            x_new[i] = u[i];
    }

};



// Example problems (one variable) ------------------------------------------------------------------
// f(x) = x*x -4 = 0
// x = 2,-2
class NL1_NR : public AbstractDefaultCloneFunction<double, double, NL1_NR>
{
public:
    // r = f(x)
    void eval(const double &x, double &r) { r = x*x-4.; }
};

class NL1_NR_D1 : public AbstractDefaultCloneFunction<double, double, NL1_NR_D1>
{
public:
    // j = f'(x)
    void eval(const double &x, double &j) { j = 2*x; }
};

class NL1_PC : public AbstractDefaultCloneFunction<double, double, NL1_PC>
{
public:
    // x = g(x) = x - 1/f'(x) * f(x)
    //void eval0(double &x, double &x_new) { x_new = 0.5*(2./x+x); }
    void eval(const double &x, double &x_new) { x_new = (x*x+4.)/(2.*x); }
};


// Example problems (2 variables) ------------------------------------------------------------------
// 3x-y=-2
// 2x^2-y=0
// (x,y) = (-1/2, 1/2) and (2, 8)
class NL2_NR
{
public:
    void eval(const std::valarray<double> &x, std::valarray<double> &r)
    {
        r[0] = 3*x[0]-x[1]+2.;
        r[1] = 2*x[0]*x[0]-x[1];
    }
};

class NL2_NR_D1
{
public:
    void eval(const std::valarray<double> &x, MathLib::Matrix<double> &j)
    {
        j(0,0) = 3.;
        j(0,1) = -1.0;
        j(1,0) = 4.*x[0];
        j(1,1) = -1.0;
    }
};

typedef TemplateFixedPointFunctonUsingC1Equation<NL2_NR,NL2_NR_D1> NL2_PC1;


class NL2_PC2
{
    MathLib::DenseLinearEquation solver;
    std::valarray<double> b;
public:
    NL2_PC2(size_t n)
    {
        solver.create(n);
        b.resize(n);
    }
    // if f=A*x+b, A*g=b
    void eval(std::valarray<double> &x, std::valarray<double> &x_new)
    {
        MathLib::Matrix<double> *j = solver.getA();
        solver.reset();
        b = 0;

        (*j)(0,0) = 3;
        (*j)(0,1) = -1.0;
        (*j)(1,0) = 2*x[0];
        (*j)(1,1) = -1;
        b[0] = -2;
        b[1] = 0.;
        //(*j)(0,0) = x[1];
        //(*j)(0,1) = -1.0;
        //(*j)(1,0) = .0;
        //(*j)(1,1) = x[1];
        //b[0] = .0;
        //b[1] = 1.;
        for (size_t i=0; i<x.size(); i++)
            solver.setRHS(i, b[i]);

        solver.solve();

        double *u = solver.getX();
        for (size_t i=0; i<x.size(); i++)
            x_new[i] = u[i];
    }

};


class NL3_NR
{
public:
  NL3_NR()
  {
      const size_t n = 2;
      MathLib::RowMajorSparsity sparse(n);
      sparse[0].insert(0);
      sparse[0].insert(1);
      sparse[1].insert(0);
      sparse[1].insert(1);
      _linear_solver.create(n, &sparse);
      _linear_solver.getOption().ls_method = MathLib::LIS_option::BiCGSTAB;
  }

  void eval(const std::valarray<double> &x, std::valarray<double> &r)
  {
        r[0] = 3*x[0]-x[1]+2.;
        r[1] = 2*x[0]*x[0]-x[1];
  }

  MathLib::LisLinearEquation* getLinearSolver() {return &_linear_solver;};


private:
  MathLib::LisLinearEquation _linear_solver;
};

class NL3_NR_D1
{
public:
    NL3_NR_D1(MathLib::LisLinearEquation* solver) : _linear_solver(solver)
    {
    }

    void eval(const std::valarray<double> &x, MathLib::CRSMatrix<double, signed> &j)
    {
        j.setValue(0, 0, 3.);
        j.setValue(0, 1, -1.0);
        j.setValue(1, 0, 4.*x[0]);
        j.setValue(1, 1,  -1.0);
    }

    MathLib::LisLinearEquation* getLinearSolver() {return _linear_solver;};


private:
  MathLib::LisLinearEquation* _linear_solver;
};

class NL4_Residual
{
public:
    void eval(const std::valarray<double> &x, std::valarray<double> &r)
    {
        double P = 1.;
        double R = 10.;
        double s = sqrt(2.);
        r[1-1]= (9*P*x[1-1])/4 + (9*x[2-1]*x[3-1])/(8*s) + (P*R*x[7-1])/s;
        r[2-1]= (81*P*x[2-1])/4 + (9*x[1-1]*x[3-1])/(8*s) + (P*R*x[8-1])/s;
        r[3-1]= (-9*x[1-1]*x[2-1])/(4*s) + 9*P*x[3-1] + s*P*R*x[9-1];
        r[4-1]= 36*P*x[4-1] + s*P*R*x[10-1];
        r[5-1]= -2*x[5-1] + (x[2-1]*x[7-1])/(2*s) + (x[1-1]*x[8-1])/(2*s) - (x[4-1]*x[9-1])/s + s*x[4-1]*x[9-1] - (x[3-1]*x[10-1])/s + s*x[3-1]*x[10-1];
        r[6-1]= -8*x[6-1] - (x[1-1]*x[7-1])/s - s*x[3-1]*x[9-1];
        r[7-1]= -(x[1-1]/s) - (x[2-1]*x[5-1])/(2*s) + (x[1-1]*x[6-1])/s - (3*x[7-1])/2.0 + (3*x[3-1]*x[8-1])/(4*s) + (3*x[2-1]*x[9-1])/(4*s);
        r[8-1]= -(x[2-1]/s) - (x[1-1]*x[5-1])/(2*s) - (3*x[3-1]*x[7-1])/(4*s) - (9*x[8-1])/2.0 - (3*x[1-1]*x[9-1])/(4*s);
        r[9-1]= -(s*x[3-1]) - (x[4-1]*x[5-1])/s + s*x[3-1]*x[6-1] - (3*x[2-1]*x[7-1])/(4*s) + (3*x[1-1]*x[8-1])/(4*s) - 3*x[9-1];
        r[10-1]= -(s*x[4-1]) - (x[3-1]*x[5-1])/s - 6*x[10-1];
    }
};

class NL4_Jacobian
{
public:
    void eval(const std::valarray<double> &x, MathLib::Matrix<double> &j)
    {
        double P = 1.;
        double R = 10.;
        double s = sqrt(2.);
        j = .0;
        j(1-1,1-1) = (9*P)/4.0;
        j(7-1,1-1) = -(1/s)+x[6-1]/s;
        j(1-1,2-1) = (9*x[3-1])/(8*s);
        j(7-1,2-1) = -x[5-1]/(2*s) + (3*x[9-1])/(4*s);
        j(1-1,3-1) = (9*x[2-1])/(8*s);
        j(7-1,3-1) = (3*x[8-1])/(4*s);
        j(1-1,7-1) = (P*R)/s;
        j(7-1,5-1) = -x[2-1]/(2*s);
        j(2-1,1-1) = (9*x[3-1])/(8*s);
        j(7-1,6-1) = x[1-1]/s;
        j(2-1,2-1) = (81*P)/4.0;
        j(7-1,7-1) = -1.5;
        j(2-1,3-1) = (9*x[1-1])/(8*s);
        j(7-1,8-1) = (3*x[3-1])/(4*s);
        j(2-1,8-1) = (P*R)/s;
        j(7-1,9-1) = (3*x[2-1])/(4*s);
        j(3-1,1-1) = (-9*x[2-1])/(4*s);
        j(8-1,1-1) = -x[5-1]/(2*s) - (3*x[9-1])/(4*s);
        j(3-1,2-1) = (-9*x[1-1])/(4*s);
        j(8-1,2-1) = -(1/s);
        j(3-1,3-1) = 9*P;
        j(8-1,3-1) = (-3*x[7-1])/(4*s);
        j(3-1,9-1) = s*P*R;
        j(8-1,5-1) = -x[1-1]/(2*s);
        j(4-1,4-1) = 36*P;
        j(8-1,7-1) = (-3*x[3-1])/(4*s);
        j(4-1,10-1)= s*P*R;
        j(8-1,8-1) = -4.5;
        j(5-1,1-1) = x[8-1]/(2*s);
        j(8-1,9-1) = (-3*x[1-1])/(4*s);
        j(5-1,2-1) = x[7-1]/(2*s);
        j(9-1,1-1) = (3*x[8-1])/(4*s);
        j(5-1,3-1) = -(x[10-1]/s) + s*x[10-1];
        j(9-1,2-1) = (-3*x[7-1])/(4*s);
        j(5-1,4-1) = -(x[9-1]/s) + s*x[9-1];
        j(9-1,3-1) = -s + s*x[6-1];
        j(5-1,5-1) = -2.0;
        j(9-1,4-1) = -(x[5-1]/s);
        j(5-1,7-1) = x[2-1]/(2*s);
        j(9-1,5-1) = -(x[4-1]/s);
        j(5-1,8-1) = x[1-1]/(2*s);
        j(9-1,6-1) = s*x[3-1];
        j(5-1,9-1) = -(x[4-1]/s) + s*x[4-1];
        j(9-1,7-1) = (-3*x[2-1])/(4*s);
        j(5-1,10-1)= -(x[3-1]/s) + s*x[3-1];
        j(9-1,8-1) = (3*x[1-1])/(4*s);
        j(6-1,1-1) = -(x[7-1]/s);
        j(9-1,9-1) = -3.0;
        j(6-1,3-1) = -(s*x[9-1]);
        j(10-1,3-1) = -(x[5-1]/s);
        j(6-1,6-1) = -8.0;
        j(10-1,4-1) = -s;
        j(6-1,7-1) = -(x[1-1]/s);
        j(10-1,5-1) = -(x[3-1]/s);
        j(6-1,9-1) = -(s*x[3-1]);
        j(10-1,10-1)= -6.0;
    }
};

//------------------------------------------------------------------

TEST(Math, NonlinearNR_double)
{
    NL1_NR f1;
    NL1_NR_D1 df1;
    double x0 = 6.0;
    double x = .0;
    NewtonRaphsonMethod nr;
    nr.solve(f1, df1, x0, x);
    ASSERT_NEAR(2.0, x, 1e-5);
}

TEST(Math, NonlinearNR_dense)
{
    typedef std::valarray<double> MyVector;

    NL2_NR f2;
    NL2_NR_D1 df2;
    MyVector x0(6.0, 2);
    MyVector x(0.0, 2);
    NewtonRaphsonMethod nr;
    nr.solve(f2,df2, x0, x);

    double my_expect[] = {2., 8.};
//    double my_expect[] = {1., 1.};
    ASSERT_DOUBLE_ARRAY_EQ(my_expect, x, 2, 1e-5);
}

TEST(Math, NonlinearNR_dense2)
{
    typedef std::valarray<double> MyVector;

    NL4_Residual f2;
    NL4_Jacobian df2;
    size_t n = 10;
    MyVector x0(1.0, n);
    MyVector x(0.0, n);
    NewtonRaphsonMethod nr;
    nr.solve(f2, df2, x0, x);

    double my_expect[] = {3.39935, 3.70074e-018, -1.42576e-017, 1.4903e-021, 4.35602e-018, 0.325, -1.08167, -5.61495e-018, 7.58394e-018, -3.79368e-021};
////    double my_expect[] = {1., 1.};
    ASSERT_DOUBLE_ARRAY_EQ(my_expect, x, 10, 1e-5);
}

TEST(Math, NonlinearNR_sparse)
{
    typedef std::valarray<double> MyVector;
    typedef MathLib::CRSMatrix<double, signed> MyMatrix;
    typedef NRCheckConvergence<MyVector,NRErrorAbsResMNormOrRelDxMNorm > MyConverge;
    typedef NewtonFunctionDXVector<NL3_NR_D1, MathLib::LisLinearEquation> MyDxFunction;

    NL3_NR f;
    NL3_NR_D1 df(f.getLinearSolver());
    MyVector x0(6.0, 2);
    MyVector x(0.0, 2);
    MyVector r(2), dx(2);
    //MyMatrix* j = f.getLinearSolver()->getA();
    MyDxFunction f_dx(df, *f.getLinearSolver());

    NewtonRaphsonMethod nr;
    nr.solve<NL3_NR,MyDxFunction,MyVector,MyConverge, NRIterationStepInitializerDummy>(f, f_dx, x0, x, r, dx);

    double my_expect[] = {2., 8.};
    ASSERT_DOUBLE_ARRAY_EQ(my_expect, x, 2, 1e-5);
}

TEST(Math, NonlinearPicard_double)
{
    NL1_PC f1;
    double x0 = 6.0;
    double x = .0;
    PicardMethod nr;
    nr.solve(f1, x0, x);
//    ASSERT_NEAR(sqrt(2.0), x, 1e-5);
    ASSERT_NEAR(2.0, x, 1e-5);
}

TEST(Math, NonlinearPicard_dense1)
{
    typedef std::valarray<double> MyVector;

    NL2_PC1 f2(2);
    MyVector x0(6.0, 2);
    MyVector x(0.0, 2);
    PicardMethod nr;
    nr.solve(f2, x0, x);

//    double my_expect[] = {1., 1.};
    double my_expect[] = {2., 8.};
    ASSERT_DOUBLE_ARRAY_EQ(my_expect, x, 2, 1e-5);
}

TEST(Math, NonlinearPicard_dense2)
{
    typedef std::valarray<double> MyVector;

    NL2_PC2 f2(2);
    MyVector x0(6.0, 2);
    MyVector x(0.0, 2);
    PicardMethod nr;
    nr.solve(f2, x0, x);

    double my_expect[] = {-0.5, 0.5};
    ASSERT_DOUBLE_ARRAY_EQ(my_expect, x, 2, 1e-5);
//    ASSERT_DOUBLE_ARRAY_EQ(my_expect, x, 2, 1e-5);
}

TEST(Math, NonlinearBS_double)
{
    NL1_NR f1;
    double x = .0;
    double a = .0;
    double c = 6.0;
    BisectionMethod nr;
    nr.solve(f1, a, c, x);
    ASSERT_NEAR(2.0, x, 1e-5);
}

