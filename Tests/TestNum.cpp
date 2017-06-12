/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TestNum.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include <gtest/gtest.h>

#include <vector>

#include "BaseLib/CodingTools.h"

#include "GeoLib/Point.h"

#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TimeStepping/TimeStepFunction.h"

#include "NumLib/Function/TXFunctionAnalytical.h"

//#include "NumLib/Output/Output.h"

#include "TestUtil.h"

TEST(Num, TimeStepping1)
{
    // define problems and solution strategy
    class TestProblem : public NumLib::ITransientSystem
    {
    private:
        double _t0, _tn, _dt, _v;
    public:
        TestProblem(double t0, double tn, double dt) : _t0(t0), _tn(tn), _dt(dt), _v(.0) {};
        virtual ~TestProblem() {};
        int solveTimeStep(const NumLib::TimeStep &/*time*/) { _v+=1.0; return 0; }
        double suggestNext(const NumLib::TimeStep &time_current)
        {
            if (time_current.getTime() < _t0) return _t0;
            else return (time_current.getTime()+_dt <= _tn) ? time_current.getTime()+_dt : time_current.getTime();
        }
        bool isAwake(const NumLib::TimeStep &time)  { return (time.getTime()>=_t0 && time.getTime()<=_tn); }
        void accept(const NumLib::TimeStep &/*time*/) {};
        double getValue() {return _v;};
    };
    TestProblem problem(.0, 100., 10.);

    // pass it to discretization systems
    NumLib::TimeSteppingController timeStepping;
    timeStepping.setTransientSystem(problem);
    // start time stepping
    timeStepping.setBeginning(.0);
    size_t n_steps = timeStepping.solve(100.);
    //check
    ASSERT_EQ(10U, n_steps);
    ASSERT_EQ(10., problem.getValue());
}

TEST(Num, TXFunctionAnalytical)
{
    std::string str_f1 = "100 + 1*x + 2*y + 3*z";
    NumLib::TXFunctionAnalytical tx_f1(str_f1);

    double v;
    {
        GeoLib::Point pt1(0., 0., 0.);
        tx_f1.eval(NumLib::TXPosition(pt1.getData()), v);
        ASSERT_NEAR(100., v, 1e-8);
    }

    {
        GeoLib::Point pt1(1., 1., 1.);
        tx_f1.eval(NumLib::TXPosition(pt1.getData()), v);
        ASSERT_NEAR(106., v, 1e-3);
    }

    std::string str_f2 = "-100+-1*x+-2*y+-3*z";
    NumLib::TXFunctionAnalytical tx_f2(str_f2);

    {
        GeoLib::Point pt1(0., 0., 0.);
        tx_f2.eval(NumLib::TXPosition(pt1.getData()), v);
        ASSERT_NEAR(-100., v, 1e-3);
    }

    {
        GeoLib::Point pt1(1., 1., 1.);
        tx_f2.eval(NumLib::TXPosition(pt1.getData()), v);
        ASSERT_NEAR(-106., v, 1e-3);
    }

}
