/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TestCoupling.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include <gtest/gtest.h>
#include <vector>

#include "BaseLib/Options.h"
#include "NumLib/Function/Function.h"
#include "NumLib/Coupling/MonolithicProblem.h"
#include "NumLib/Coupling/PartitionedProblem.h"
#include "NumLib/Coupling/CouplingStrucutreBuilder.h"
#include "NumLib/Coupling/Algorithm/PartitionedAlgorithm.h"
#include "NumLib/Coupling/Algorithm/TransientPartitionedAlgorithm.h"
#include "NumLib/Coupling/Algorithm/IConvergenceCheck.h"
#include "NumLib/TransientCoupling/AsyncPartitionedSystem.h"
#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"
#include "NumLib/TransientCoupling/TransientCouplingStructureBuilder.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"

using namespace MathLib;
using namespace NumLib;

typedef NumLib::FunctionConstant<double,double> MyFunction;

class MyConvergenceCheck : public IConvergenceCheck
{
public:
    virtual ~MyConvergenceCheck() {};

    virtual bool isConverged(UnnamedParameterSet& vars_prev, UnnamedParameterSet& vars_current, double eps, double &v_diff)
    {
        for (size_t i=0; i<vars_prev.size(); i++) {
            double v_prev = .0;
            vars_prev.get<MyFunction>(i)->eval(v_prev);
            double v_cur = .0;
            vars_current.get<MyFunction>(i)->eval(v_cur);
            v_diff = std::abs(v_cur - v_prev);
            if (v_diff>eps) {
                return false;
            }
        }
        return true;
    }
};

#if 0
// 2a + 2b + 3c = 15
// 3a + 5b + 2c = 19
// 5a + 3b + 3c = 20
// A. a=1, b=2, c=3
class StrongCouplingEQS1 : public IMonolithicProblem
{
public:
    enum InputParameters { b = 0, c = 1 };
    enum OutputParameters { a = 0 };
    size_t getNumberOfInputVarameters() const {return 2;};
    size_t getNumberOfOutputParameters() const {return 1;};
    StrongCouplingEQS1() {
        _vec_in_var.resize(getNumberOfInputVarameters(), 0);
        _vec_out_var.resize(getNumberOfOutputParameters(), 0);
    }
    int solve()
    {
        double vb = _vec_in_var[EQS1::b]->eval(.0);
        double vc = _vec_in_var[EQS1::c]->eval(.0);
        double va = 1./2.*(15. - 2.*vb - 3.*vc);
        if (_vec_out_var[EQS1::a]!=0)
            delete _vec_out_var[EQS1::a];
        _vec_out_var[EQS1::a] = new MathLib::FunctionConstant<double,double>(va);
        return 0;
    }
};

class StrongCouplingEQS2 : public IMonolithicProblem
{
public:
    enum InputParameters { a = 0, c = 1 };
    enum OutputParameters { b = 0 };
    size_t getNumberOfInputVarameters() const {return 2;};
    size_t getNumberOfOutputParameters() const {return 1;};
    StrongCouplingEQS2() {
        _vec_in_var.resize(getNumberOfInputVarameters(), 0);
        _vec_out_var.resize(getNumberOfOutputParameters(), 0);
    }
    int solve()
    {
        double va = _vec_in_var[EQS2::a]->eval(.0);
        double vc = _vec_in_var[EQS2::c]->eval(.0);
        double vb = 1./5.*(19.-3*va-2.*vc);
        if (_vec_out_var[EQS2::b]!=0)
            delete _vec_out_var[EQS2::b];
        _vec_out_var[EQS2::b] = new MathLib::FunctionConstant<double,double>(vb);
        return 0;
    }
};

class StrongCouplingEQS3 : public IMonolithicProblem
{
public:
    enum InputParameters { a = 0, b = 1 };
    enum OutputParameters { c = 0 };
    size_t getNumberOfInputVarameters() const {return 2;};
    size_t getNumberOfOutputParameters() const {return 1;};
    StrongCouplingEQS3() {
        _vec_in_var.resize(getNumberOfInputVarameters(), 0);
        _vec_out_var.resize(getNumberOfOutputParameters(), 0);
    }
    int solve()
    {
        double va = _vec_in_var[a]->eval(.0);
        double vb = _vec_in_var[b]->eval(.0);
        double vc = 1./3.*(20.-5*va-3.*vb);
        if (_vec_out_var[c]!=0)
            delete _vec_out_var[c];
        _vec_out_var[c] = new MathLib::FunctionConstant<double,double>(vc);
        return 0;
    }
private:
    double va, vb, vc;
};
#endif

// 2a + 2b + .3c = 6.9
// 3a + 5b + .2c = 13.6
// .5a + .3b + 3c = 10.1
// A. a=1, b=2, c=3
class WeakCouplingEQS1 : public TemplateSteadyMonolithicSystem
{
    enum In {b=0, c=1};
    enum Out {a=0};
public:
    WeakCouplingEQS1() 
    {
        resizeInputParameter(2);
        resizeOutputParameter(1);

        setOutput(a, new FunctionConstant<double,double>(.0));
    }
    MyConvergenceCheck _checker;
    virtual NumLib::IConvergenceCheck* getConvergenceChecker()
    {
        return &_checker;
    }

    int solve()
    {
        double vb, vc;
        getInput<MyFunction>(b)->eval(vb);
        getInput<MyFunction>(c)->eval(vc);
        double va = 1./2.*(6.9 - 2.*vb - 0.3*vc);
        setOutput(a, new FunctionConstant<double,double>(va));
        return 0;
    }
};

class WeakCouplingEQS2 :  public TemplateSteadyMonolithicSystem
{
    enum In {a=0, c=1};
    enum Out {b=0};
public:

    WeakCouplingEQS2() 
    {
        resizeInputParameter(2);
        resizeOutputParameter(1);
        setOutput(b, new FunctionConstant<double,double>(.0));
    }
    MyConvergenceCheck _checker;
    virtual NumLib::IConvergenceCheck* getConvergenceChecker()
    {
        return &_checker;
    }

    int solve()
    {
        double va, vc;
        getInput<MyFunction>(a)->eval(va);
        getInput<MyFunction>(c)->eval(vc);
        double vb = 1./5.*(13.6-3*va-0.2*vc);
        setOutput(b, new FunctionConstant<double,double>(vb));
        return 0;
    }
};

class WeakCouplingEQS3 : public TemplateSteadyMonolithicSystem
{
    enum In {a=0, b=1};
    enum Out {c=0};
public:

    WeakCouplingEQS3() 
    : va(.0), vb(.0), vc(.0)
    {
        resizeInputParameter(2);
        resizeOutputParameter(1);
        setOutput(c, new FunctionConstant<double,double>(.0));
    }
    MyConvergenceCheck _checker;
    virtual NumLib::IConvergenceCheck* getConvergenceChecker()
    {
        return &_checker;
    }

    int solve()
    {
        double va, vb;
        getInput<MyFunction>(a)->eval(va);
        getInput<MyFunction>(b)->eval(vb);
        double vc = 1./3.*(10.1-0.5*va-0.3*vb);
        setOutput(c, new FunctionConstant<double,double>(vc));
        return 0;
    }
private:
    double va, vb, vc;
};

// P(P(M1,M2),M3)
BaseLib::Options* defineOption4SteadyCoupling()
{
    BaseLib::Options* options = new BaseLib::Options();
    BaseLib::Options* coupling = options->addSubGroup("coupling");
    BaseLib::Options* P2 = coupling->addSubGroup("P");
    {
    //P2->addOption("name", "P2");
    P2->addOption("algorithm", "Jacobi");
    P2->addOption("convergence", "MyConvergenceCheck");
    P2->addOptionAsNum("max_itr", 100);
    P2->addOptionAsNum("epsilon", 1.e-4);
    std::vector<std::string> out_var;
//    out_var.push_back("a");
//    out_var.push_back("b");
//    out_var.push_back("c");
//    P2->addOptionAsArray("out", out_var);
    P2->addOption("out", "a");
    P2->addOption("out", "b");
    P2->addOption("out", "c");
    }
    BaseLib::Options* P2_sub = P2->addSubGroup("problems");
    BaseLib::Options* P1 = P2_sub->addSubGroup("P");
    {
    //P1->addOption("name", "P1");
    P1->addOption("algorithm", "Gauss");
    P1->addOption("convergence", "MyConvergenceCheck");
    P1->addOptionAsNum("max_itr", 100);
    P1->addOptionAsNum("epsilon", 1.e-4);
//    std::vector<std::string> out_var;
//    out_var.push_back("a");
//    out_var.push_back("b");
//    std::vector<std::string> in_var;
//    in_var.push_back("c");
//    P1->addOptionAsArray("out", out_var);
//    P1->addOptionAsArray("in", in_var);
    P1->addOption("out", "a");
    P1->addOption("out", "b");
    P1->addOption("in", "c");
    }
    BaseLib::Options* P1_sub = P1->addSubGroup("problems");
    BaseLib::Options* M1 = P1_sub->addSubGroup("M");
    {
    M1->addOption("type", "EQS1");
//    std::vector<std::string> out_var;
//    out_var.push_back("a");
//    std::vector<std::string> in_var;
//    in_var.push_back("b");
//    in_var.push_back("c");
//    M1->addOptionAsArray("out", out_var);
//    M1->addOptionAsArray("in", in_var);
    M1->addOption("out", "a");
    M1->addOption("in", "b");
    M1->addOption("in", "c");
    }
    BaseLib::Options* M2 = P1_sub->addSubGroup("M");
    {
    M2->addOption("type", "EQS2");
//    std::vector<std::string> out_var;
//    out_var.push_back("b");
//    std::vector<std::string> in_var;
//    in_var.push_back("a");
//    in_var.push_back("c");
//    M2->addOptionAsArray("out", out_var);
//    M2->addOptionAsArray("in", in_var);
    M2->addOption("out", "b");
    M2->addOption("in", "a");
    M2->addOption("in", "c");
    }
    BaseLib::Options* M3 = P2_sub->addSubGroup("M");
    {
    M3->addOption("type", "EQS3");
//    std::vector<std::string> out_var;
//    out_var.push_back("c");
//    std::vector<std::string> in_var;
//    in_var.push_back("a");
//    in_var.push_back("b");
//    M3->addOptionAsArray("out", out_var);
//    M3->addOptionAsArray("in", in_var);
    M3->addOption("out", "c");
    M3->addOption("in", "a");
    M3->addOption("in", "b");
    }

    return options;
};

class MyEQSFactory
{
public:
    TemplateSteadyMonolithicSystem* create(const std::string &eqs_name)
    {
        if (eqs_name.compare("EQS1")==0) {
            return new WeakCouplingEQS1();
        } else if (eqs_name.compare("EQS2")==0) {
                return new WeakCouplingEQS2();
        } else if (eqs_name.compare("EQS3")==0) {
            return new WeakCouplingEQS3();
        }
        return 0;
    };

};

TEST(Coupling, SteadyCouplingOption)
{
    BaseLib::Options* option = defineOption4SteadyCoupling();
    MyEQSFactory eqsFac;
    CouplingStrucutreBuilder cpl_builder;
    PartitionedProblem *coupled_sys = (PartitionedProblem*)cpl_builder.build(option, eqsFac);
    ASSERT_TRUE(coupled_sys->check());
    ASSERT_EQ(2u, coupled_sys->getNumberOfSubProblems());
    coupled_sys->solve();

    const double epsilon = 1.e-3;
    double v1, v2, v3;
    size_t outId_A = coupled_sys->getOutputParameterID("a");
    const MyFunction* f1 = (const MyFunction*)coupled_sys->getOutput(outId_A);
    f1->eval(v1);
    ((const MyFunction*)coupled_sys->getOutput(coupled_sys->getOutputParameterID("b")))->eval(v2);
    ((const MyFunction*)coupled_sys->getOutput(coupled_sys->getOutputParameterID("c")))->eval(v3);
    ASSERT_NEAR(1., v1, epsilon);
    ASSERT_NEAR(2., v2, epsilon);
    ASSERT_NEAR(3., v3, epsilon);

}

TEST(Coupling, SteadyCouplingCheck1)
{
    WeakCouplingEQS1 eqs1;
    eqs1.setOutputParameterName(0,"a");
    eqs1.setInputParameterName(0,"b");
    eqs1.setInputParameterName(1,"c");
    WeakCouplingEQS2 eqs2;
    eqs2.setOutputParameterName(0,"b");
    eqs2.setInputParameterName(0,"a");
    eqs2.setInputParameterName(1,"c");
    WeakCouplingEQS3 eqs3;
    eqs3.setOutputParameterName(0,"c");
    eqs3.setInputParameterName(0,"a");
    eqs3.setInputParameterName(1,"b");

    ASSERT_TRUE(eqs1.isValid());
    ASSERT_TRUE(eqs2.isValid());
    ASSERT_TRUE(eqs3.isValid());

    BlockJacobiMethod method(1.e-4, 100);
    MyConvergenceCheck checker;
    //method.setConvergenceCheck(checker);
    {
        // correct
        //BlockJacobiMethod<MyConvergenceCheck> method(1.e-4, 100);
        PartitionedProblem part1;
        part1.setAlgorithm(method);
        part1.resizeInputParameter(1);
        part1.resizeOutputParameter(2);
        part1.setOutputParameterName(0, "a");
        part1.setOutputParameterName(1, "b");
        part1.setInputParameterName(0, "c");
        part1.addProblem(eqs1);
        part1.addProblem(eqs2);
        part1.connectParameters();

        PartitionedProblem part2;
        part2.setAlgorithm(method);
        part2.resizeInputParameter(0);
        part2.resizeOutputParameter(3);
        part2.setOutputParameterName(0, "a");
        part2.setOutputParameterName(1, "b");
        part2.setOutputParameterName(2, "c");
        part2.addProblem(part1);
        part2.addProblem(eqs3);
        part2.connectParameters();

        ASSERT_TRUE(part2.check());
    }

    {
        // no source defined for c
        PartitionedProblem part1;
        part1.setAlgorithm(method);
        part1.resizeInputParameter(0);
        part1.resizeOutputParameter(2);
        part1.setOutputParameterName(0, "a");
        part1.setOutputParameterName(1, "b");
        //part1.setInputParameterName(0, "c");
        part1.addProblem(eqs1);
        part1.addProblem(eqs2);
        part1.connectParameters();

        ASSERT_FALSE(part1.check());
    }

    //{
    //    // no source defined for c
    //    BlockJacobiMethod<MyConvergenceCheck> method(1.e-4, 100);
    //    PartitionedProblem part1;
    //    part1.setAlgorithm(method);
    //    part1.resizeInputParameter(1);
    //    part1.resizeOutputParameter(2);
    //    part1.setOutputParameterName(0, "a");
    //    part1.setOutputParameterName(1, "b");
    //    part1.setInputParameterName(0, "c");
    //    part1.addProblem(eqs1);
    //    part1.addProblem(eqs2);
    //    part1.connectParameters();

    //    PartitionedProblem part2;
    //    part2.setAlgorithm(method);
    //    part2.resizeInputParameter(0);
    //    part2.resizeOutputParameter(3);
    //    part2.setOutputParameterName(0, "a");
    //    part2.setOutputParameterName(1, "b");
    //    part2.setOutputParameterName(2, "c");
    //    part2.addProblem(part1);
    //    part2.connectParameters();

    //    ASSERT_FALSE(part2.check());
    //}
}

void defineSteadyExample1(WeakCouplingEQS1 &eqs1, WeakCouplingEQS2 &eqs2, WeakCouplingEQS3 &eqs3, PartitionedProblem &part1, PartitionedProblem &part2, IPartitionedAlgorithm &algorithm1, IPartitionedAlgorithm &algorithm2)
{
    eqs1.setOutputParameterName(0,"a");
    eqs1.setInputParameterName(0,"b");
    eqs1.setInputParameterName(1,"c");
    eqs2.setOutputParameterName(0,"b");
    eqs2.setInputParameterName(0,"a");
    eqs2.setInputParameterName(1,"c");
    eqs3.setOutputParameterName(0,"c");
    eqs3.setInputParameterName(0,"a");
    eqs3.setInputParameterName(1,"b");

    part1.setAlgorithm(algorithm1);
    part1.resizeInputParameter(1);
    part1.resizeOutputParameter(2);
    part1.setOutputParameterName(0, "a");
    part1.setOutputParameterName(1, "b");
    part1.setInputParameterName(0, "c");
    part1.addProblem(eqs1);
    part1.addProblem(eqs2);
    part1.connectParameters();

    part2.setAlgorithm(algorithm2);
    part2.resizeInputParameter(0);
    part2.resizeOutputParameter(3);
    part2.setOutputParameterName(0, "a");
    part2.setOutputParameterName(1, "b");
    part2.setOutputParameterName(2, "c");
    part2.addProblem(part1);
    part2.addProblem(eqs3);
    part2.connectParameters();
}

TEST(Coupling, SteadyCouplingJacobi)
{
    WeakCouplingEQS1 eqs1;
    WeakCouplingEQS2 eqs2;
    WeakCouplingEQS3 eqs3;
    PartitionedProblem part1;
    PartitionedProblem part2;

    MyConvergenceCheck checker;
    BlockJacobiMethod method(1.e-4, 100);

    defineSteadyExample1(eqs1, eqs2, eqs3, part1, part2, method, method);
    ASSERT_TRUE(part2.check());

    part2.solve();

    const double epsilon = 1.e-3;
    double v1, v2, v3;
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(v3);
    ASSERT_NEAR(1., v1, epsilon);
    ASSERT_NEAR(2., v2, epsilon);
    ASSERT_NEAR(3., v3, epsilon);
}

TEST(Coupling, SteadyCouplingSeidel)
{
    WeakCouplingEQS1 eqs1;
    WeakCouplingEQS2 eqs2;
    WeakCouplingEQS3 eqs3;
    PartitionedProblem part1;
    PartitionedProblem part2;

    MyConvergenceCheck checker;
    BlockGaussSeidelMethod method(1.e-5, 100);

    defineSteadyExample1(eqs1, eqs2, eqs3, part1, part2, method, method);

    ASSERT_TRUE(part2.check());
    part2.solve();

    const double epsilon = 1.e-3;
    double v1, v2, v3;
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(v3);
    ASSERT_NEAR(1., v1, epsilon);
    ASSERT_NEAR(2., v2, epsilon);
    ASSERT_NEAR(3., v3, epsilon);
}

// 2a + 2*b + .3*c = 6.9*t 
// 3a + 5b + .2c = 13.6*t
// .5a + .3b + 3c = 10.1*t
// A. a=t, b=2*t, c=3*t
class TransientWeakCouplingEQS1 : public AbstractTransientMonolithicSystem
{
    double _dt;
    enum In {b=0, c=1};
    enum Out {a=0};
public:
    TransientWeakCouplingEQS1() : _dt(.0)
    {
        reset();
    };

    TransientWeakCouplingEQS1(double dt) : _dt(dt) 
    {
        reset();
    };
    MyConvergenceCheck _checker;
    virtual NumLib::IConvergenceCheck* getConvergenceChecker()
    {
        return &_checker;
    }

    void setTimeStepSize(double dt) {_dt = dt;};

    void reset()
    {
        resizeInputParameter(2);
        resizeOutputParameter(1);
        setOutput(a, new FunctionConstant<double,double>(.0));
    }

    int solveTimeStep(const TimeStep &ts) 
    {
        double t = ts.getTime();
        double vb, vc;
        getInput<MyFunction>(b)->eval(vb);
        getInput<MyFunction>(c)->eval(vc);
        double va = 1./2.*(6.9*t - 2.*vb - 0.3*vc);
        setOutput(a, new FunctionConstant<double,double>(va));
        return 0;
    }
    double suggestNext(const TimeStep &ts) {return ts.getTime()+_dt;};
    bool isAwake(const TimeStep &ts) 
    {
        double t = ts.getTime();
        double mod = t - static_cast<int>(t/_dt)*_dt;
        if (mod==.0) return true;
        return false;
    };
    void accept(const TimeStep &/*time*/) {};
};

class TransientWeakCouplingEQS2 :  public AbstractTransientMonolithicSystem
{
    double _dt;
    enum In {a=0, c=1};
    enum Out {b=0};
public:

    TransientWeakCouplingEQS2() : _dt(.0)
    {
        reset();
    };

    TransientWeakCouplingEQS2(double dt) : _dt(dt) 
    {
        reset();
    };
    MyConvergenceCheck _checker;
    virtual NumLib::IConvergenceCheck* getConvergenceChecker()
    {
        return &_checker;
    }

    void setTimeStepSize(double dt) {_dt = dt;};

    void reset()
    {
        resizeInputParameter(2);
        resizeOutputParameter(1);
        setOutput(b, new FunctionConstant<double,double>(.0));
    }

    int solveTimeStep(const TimeStep &ts) 
    {
        double t = ts.getTime();
        double va, vc;
        getInput<MyFunction>(a)->eval(va);
        getInput<MyFunction>(c)->eval(vc);
        double vb = 1./5.*(13.6*t-3*va-0.2*vc);
        setOutput(b, new FunctionConstant<double,double>(vb));
        return 0;
    }
    double suggestNext(const TimeStep &ts) {return (ts.getTime()+_dt);};
    bool isAwake(const TimeStep &ts) 
    {
        double t = ts.getTime();
        double mod = t - static_cast<int>(t/_dt)*_dt;
        if (mod==.0) return true;
        return false;
    };
    void accept(const TimeStep &/*time*/) {};
};

class TransientWeakCouplingEQS3 : public AbstractTransientMonolithicSystem
{
    double _dt;
    enum In {a=0, b=1};
    enum Out {c=0};
public:

    TransientWeakCouplingEQS3() : _dt(.0)
    {
        reset();
    };

    TransientWeakCouplingEQS3(double dt) : _dt(dt) 
    {
        reset();
    };
    MyConvergenceCheck _checker;
    virtual NumLib::IConvergenceCheck* getConvergenceChecker()
    {
        return &_checker;
    }

    void setTimeStepSize(double dt) {_dt = dt;};

    void reset()
    {
        resizeInputParameter(2);
        resizeOutputParameter(1);
        setOutput(c, new FunctionConstant<double,double>(.0));
    }

    int solveTimeStep(const TimeStep &ts) 
    {
        double t = ts.getTime();
        double va, vb;
        getInput<MyFunction>(a)->eval(va);
        getInput<MyFunction>(b)->eval(vb);
        double vc = 1./3.*(10.1*t-0.5*va-0.3*vb);
        setOutput(c, new FunctionConstant<double,double>(vc));
        return 0;
    }
    double suggestNext(const TimeStep &ts) {return (ts.getTime()+_dt);};
    bool isAwake(const TimeStep &ts) 
    {
        double t = ts.getTime();
        double mod = t - static_cast<int>(t/_dt)*_dt;
        if (mod==.0) return true;
        return false;
    };
    void accept(const TimeStep &/*time*/) {};
};

void defineTransientExample1(TransientWeakCouplingEQS1 &eqs1, TransientWeakCouplingEQS2 &eqs2, TransientWeakCouplingEQS3 &eqs3, AsyncPartitionedSystem &part1, AsyncPartitionedSystem &part2, ITransientPartitionedAlgorithm &algorithm1, ITransientPartitionedAlgorithm &algorithm2)
{
    eqs1.setOutputParameterName(0,"a");
    eqs1.setInputParameterName(0,"b");
    eqs1.setInputParameterName(1,"c");
    eqs2.setOutputParameterName(0,"b");
    eqs2.setInputParameterName(0,"a");
    eqs2.setInputParameterName(1,"c");
    eqs3.setOutputParameterName(0,"c");
    eqs3.setInputParameterName(0,"a");
    eqs3.setInputParameterName(1,"b");

    part1.setAlgorithm(algorithm1);
    part1.resizeInputParameter(1);
    part1.resizeOutputParameter(2);
    part1.setOutputParameterName(0, "a");
    part1.setOutputParameterName(1, "b");
    part1.setInputParameterName(0, "c");
    part1.addProblem(eqs1);
    part1.addProblem(eqs2);
    part1.connectParameters();

    part2.setAlgorithm(algorithm2);
    part2.resizeInputParameter(0);
    part2.resizeOutputParameter(3);
    part2.setOutputParameterName(0, "a");
    part2.setOutputParameterName(1, "b");
    part2.setOutputParameterName(2, "c");
    part2.addProblem(part1);
    part2.addProblem(eqs3);
    part2.connectParameters();
}

BaseLib::Options* defineOption4TransientCoupling()
{
    BaseLib::Options* options = new BaseLib::Options();
    BaseLib::Options* coupling = options->addSubGroup("coupling");
    BaseLib::Options* P2 = coupling->addSubGroup("P");
    {
    //P2->addOption("name", "P2");
    P2->addOption("algorithm", "Parallel");
    P2->addOption("convergence", "MyConvergenceCheck");
    P2->addOptionAsNum("max_itr", 100);
    P2->addOptionAsNum("epsilon", 1.e-4);
//    std::vector<std::string> out_var;
//    out_var.push_back("a");
//    out_var.push_back("b");
//    out_var.push_back("c");
//    P2->addOptionAsArray("out", out_var);
    P2->addOption("out", "a");
    P2->addOption("out", "b");
    P2->addOption("out", "c");
    }
    BaseLib::Options* P2_sub = P2->addSubGroup("problems");
    BaseLib::Options* P1 = P2_sub->addSubGroup("P");
    {
    //P1->addOption("name", "P1");
    P1->addOption("algorithm", "Serial");
    P1->addOption("convergence", "MyConvergenceCheck");
    P1->addOptionAsNum("max_itr", 100);
    P1->addOptionAsNum("epsilon", 1.e-4);
//    std::vector<std::string> out_var;
//    out_var.push_back("a");
//    out_var.push_back("b");
//    std::vector<std::string> in_var;
//    in_var.push_back("c");
//    P1->addOptionAsArray("out", out_var);
//    P1->addOptionAsArray("in", in_var);
    P1->addOption("out", "a");
    P1->addOption("out", "b");
    P1->addOption("in", "c");
    }
    BaseLib::Options* P1_sub = P1->addSubGroup("problems");
    BaseLib::Options* M1 = P1_sub->addSubGroup("M");
    {
    M1->addOption("type", "EQS1");
//    std::vector<std::string> out_var;
//    out_var.push_back("a");
//    std::vector<std::string> in_var;
//    in_var.push_back("b");
//    in_var.push_back("c");
//    M1->addOptionAsArray("out", out_var);
//    M1->addOptionAsArray("in", in_var);
    M1->addOption("out", "a");
    M1->addOption("in", "b");
    M1->addOption("in", "c");
    }
    BaseLib::Options* M2 = P1_sub->addSubGroup("M");
    {
    M2->addOption("type", "EQS2");
//    std::vector<std::string> out_var;
//    out_var.push_back("b");
//    std::vector<std::string> in_var;
//    in_var.push_back("a");
//    in_var.push_back("c");
//    M2->addOptionAsArray("out", out_var);
//    M2->addOptionAsArray("in", in_var);
    M2->addOption("out", "b");
    M2->addOption("in", "a");
    M2->addOption("in", "c");
    }
    BaseLib::Options* M3 = P2_sub->addSubGroup("M");
    {
    M3->addOption("type", "EQS3");
//    std::vector<std::string> out_var;
//    out_var.push_back("c");
//    std::vector<std::string> in_var;
//    in_var.push_back("a");
//    in_var.push_back("b");
//    M3->addOptionAsArray("out", out_var);
//    M3->addOptionAsArray("in", in_var);
    M3->addOption("out", "c");
    M3->addOption("in", "a");
    M3->addOption("in", "b");
    }

    return options;
};

class MyTransientEQSFactory
{
    double _dt1, _dt2, _dt3;
public:
    MyTransientEQSFactory(double dt1, double dt2, double dt3)
    {
        _dt1 = dt1;
        _dt2 = dt2;
        _dt3 = dt3;
    };

    AbstractTransientMonolithicSystem* create(const std::string &eqs_name)
    {
        if (eqs_name.compare("EQS1")==0) {
            return new TransientWeakCouplingEQS1(_dt1);
        } else if (eqs_name.compare("EQS2")==0) {
                return new TransientWeakCouplingEQS2(_dt2);
        } else if (eqs_name.compare("EQS3")==0) {
            return new TransientWeakCouplingEQS3(_dt3);
        }
        return 0;
    };
};

TEST(Coupling, TransientCouplingOption)
{
    BaseLib::Options* option = defineOption4TransientCoupling();
    MyTransientEQSFactory eqsFac(1.0, 1.0, 1.0);
    TransientCoulplingStrucutreBuilder cpl_builder;
    ITransientCoupledSystem *coupled_sys = cpl_builder.build(option, eqsFac);
    ASSERT_TRUE(coupled_sys->check());

    const double epsilon = 1.e-3;
    TimeSteppingController timestepping;
    timestepping.setTransientSystem(*coupled_sys);

    timestepping.setBeginning(.0);
    timestepping.solve(1.0);
    double v1, v2, v3;
    ((const MyFunction*)coupled_sys->getOutput(coupled_sys->getOutputParameterID("a")))->eval(v1);
    ((const MyFunction*)coupled_sys->getOutput(coupled_sys->getOutputParameterID("b")))->eval(v2);
    ((const MyFunction*)coupled_sys->getOutput(coupled_sys->getOutputParameterID("c")))->eval(v3);
    ASSERT_NEAR(1., v1, epsilon);
    ASSERT_NEAR(2., v2, epsilon);
    ASSERT_NEAR(3., v3, epsilon);

    timestepping.solve(2.0);
    ((const MyFunction*)coupled_sys->getOutput(coupled_sys->getOutputParameterID("a")))->eval(v1);
    ((const MyFunction*)coupled_sys->getOutput(coupled_sys->getOutputParameterID("b")))->eval(v2);
    ((const MyFunction*)coupled_sys->getOutput(coupled_sys->getOutputParameterID("c")))->eval(v3);
    ASSERT_NEAR(2., v1, epsilon);
    ASSERT_NEAR(4., v2, epsilon);
    ASSERT_NEAR(6., v3, epsilon);

}

TEST(Coupling, TransientCouplingParallelStaggered1)
{
    TransientWeakCouplingEQS1 eqs1(1.0);
    TransientWeakCouplingEQS2 eqs2(1.0);
    TransientWeakCouplingEQS3 eqs3(1.0);
    AsyncPartitionedSystem part1;
    AsyncPartitionedSystem part2;

    MyConvergenceCheck checker;
    ParallelStaggeredMethod method(1e-5, 100);
    defineTransientExample1(eqs1, eqs2, eqs3, part1, part2, method, method);

    TimeSteppingController timestepping;
    timestepping.setTransientSystem(part2);

    const double epsilon = 1.e-3;
    timestepping.setBeginning(.0);
    timestepping.solve(1.0);
    double v1, v2, v3;
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(1., v1, epsilon);
    ASSERT_NEAR(2., v2, epsilon);
    ASSERT_NEAR(3., v3, epsilon);

    timestepping.solve(2.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(2., v1, epsilon);
    ASSERT_NEAR(4., v2, epsilon);
    ASSERT_NEAR(6., v3, epsilon);
}

TEST(Coupling, TransientCouplingParallelStaggered2)
{
    TransientWeakCouplingEQS1 eqs1(1.0);
    TransientWeakCouplingEQS2 eqs2(2.0);
    TransientWeakCouplingEQS3 eqs3(4.0);
    AsyncPartitionedSystem part1;
    AsyncPartitionedSystem part2;

    MyConvergenceCheck checker;
    ParallelStaggeredMethod method(1e-5, 100);
    defineTransientExample1(eqs1, eqs2, eqs3, part1, part2, method, method);

    TimeSteppingController timestepping;
    timestepping.setTransientSystem(part2);

    const double epsilon = 1.e-3;
    timestepping.setBeginning(.0);
    timestepping.solve(1.0);
    double v1, v2, v3;
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(6.9/2., v1, epsilon);
    ASSERT_NEAR(0., v2, epsilon);
    ASSERT_NEAR(0., v3, epsilon);

    timestepping.solve(2.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(3.65, v1, epsilon);
    ASSERT_NEAR(3.25, v2, epsilon);
    ASSERT_NEAR(0., v3, epsilon);

    timestepping.solve(3.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(7.1, v1, epsilon);
    ASSERT_NEAR(3.25, v2, epsilon);
    ASSERT_NEAR(0., v3, epsilon);

    timestepping.solve(4.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(4., v1, epsilon);
    ASSERT_NEAR(8., v2, epsilon);
    ASSERT_NEAR(12., v3, epsilon);

    timestepping.solve(5.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(7.45, v1, epsilon);
    ASSERT_NEAR(8., v2, epsilon);
    ASSERT_NEAR(12., v3, epsilon);
}


TEST(Coupling, TransientCouplingParallelStaggered3)
{
    TransientWeakCouplingEQS1 eqs1(1.0);
    TransientWeakCouplingEQS2 eqs2(2.0);
    TransientWeakCouplingEQS3 eqs3(4.0);
    AsyncPartitionedSystem part1;
    AsyncPartitionedSystem part2;

    MyConvergenceCheck checker;
    ParallelStaggeredMethod method(1e-5, 1);
    defineTransientExample1(eqs1, eqs2, eqs3, part1, part2, method, method);

    TimeSteppingController timestepping;
    timestepping.setTransientSystem(part2);

    const double epsilon = 1.e-3;
    timestepping.setBeginning(.0);
    timestepping.solve(1.0);
    double v1, v2, v3;
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(3.45, v1, epsilon);
    ASSERT_NEAR(0., v2, epsilon);
    ASSERT_NEAR(0., v3, epsilon);

    eqs1.reset();
    eqs2.reset();
    eqs3.reset();
    timestepping.solve(2.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(6.9, v1, epsilon);
    ASSERT_NEAR(3.37, v2, epsilon);
    ASSERT_NEAR(0., v3, epsilon);

    eqs1.reset();
    eqs2.reset();
    eqs3.reset();
    timestepping.solve(3.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(6.98, v1, epsilon);
    ASSERT_NEAR(3.37, v2, epsilon);
    ASSERT_NEAR(0., v3, epsilon);

    eqs1.reset();
    eqs2.reset();
    eqs3.reset();
    timestepping.solve(4.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(10.43, v1, epsilon);
    ASSERT_NEAR(6.692, v2, epsilon);
    ASSERT_NEAR(11.96633, v3, epsilon);

    eqs1.reset();
    eqs2.reset();
    eqs3.reset();
    timestepping.solve(5.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(8.76305, v1, epsilon);
    ASSERT_NEAR(6.692, v2, epsilon);
    ASSERT_NEAR(11.96633, v3, epsilon);
}

TEST(Coupling, TransientCouplingSerialStaggered1)
{
    TransientWeakCouplingEQS1 eqs1(1.0);
    TransientWeakCouplingEQS2 eqs2(2.0);
    TransientWeakCouplingEQS3 eqs3(4.0);
    AsyncPartitionedSystem part1;
    AsyncPartitionedSystem part2;

    MyConvergenceCheck checker;
    SerialStaggeredMethod method(1e-5, 100);
    defineTransientExample1(eqs1, eqs2, eqs3, part1, part2, method, method);

    TimeSteppingController timestepping;
    timestepping.setTransientSystem(part2);

    const double epsilon = 1.e-3;
    timestepping.setBeginning(.0);
    timestepping.solve(1.0);
    double v1, v2, v3;
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(6.9/2., v1, epsilon);
    ASSERT_NEAR(0., v2, epsilon);
    ASSERT_NEAR(0., v3, epsilon);

    timestepping.solve(2.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(3.65, v1, epsilon);
    ASSERT_NEAR(3.25, v2, epsilon);
    ASSERT_NEAR(0., v3, epsilon);

    timestepping.solve(3.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(7.1, v1, epsilon);
    ASSERT_NEAR(3.25, v2, epsilon);
    ASSERT_NEAR(0., v3, epsilon);

    timestepping.solve(4.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(4, v1, epsilon);
    ASSERT_NEAR(8, v2, epsilon);
    ASSERT_NEAR(12., v3, epsilon);

    timestepping.solve(5.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(7.45, v1, epsilon);
    ASSERT_NEAR(8, v2, epsilon);
    ASSERT_NEAR(12., v3, epsilon);
}

TEST(Coupling, TransientCouplingSerialStaggered2)
{
    TransientWeakCouplingEQS1 eqs1(1.0);
    TransientWeakCouplingEQS2 eqs2(2.0);
    TransientWeakCouplingEQS3 eqs3(4.0);

    AsyncPartitionedSystem part1;
    AsyncPartitionedSystem part2;

    MyConvergenceCheck checker;
    SerialStaggeredMethod method(1e-5, 1);
    defineTransientExample1(eqs1, eqs2, eqs3, part1, part2, method, method);

    TimeSteppingController timestepping;
    timestepping.setTransientSystem(part2);

    const double epsilon = 1.e-3;
    timestepping.setBeginning(.0);
    timestepping.solve(1.0);
    double v1, v2, v3;
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(6.9/2., v1, epsilon);
    ASSERT_NEAR(0., v2, epsilon);
    ASSERT_NEAR(0., v3, epsilon);

    eqs1.reset();
    eqs2.reset();
    eqs3.reset();
    timestepping.solve(2.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(6.9, v1, epsilon);
    ASSERT_NEAR(1.3, v2, epsilon);
    ASSERT_NEAR(0., v3, epsilon);

    eqs1.reset();
    eqs2.reset();
    eqs3.reset();
    timestepping.solve(3.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(9.05, v1, epsilon);
    ASSERT_NEAR(1.3, v2, epsilon);
    ASSERT_NEAR(0., v3, epsilon);

    eqs1.reset();
    eqs2.reset();
    eqs3.reset();
    timestepping.solve(4.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0, v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0, v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0, v3);
    ASSERT_NEAR(12.5, v1, epsilon);
    ASSERT_NEAR(3.38, v2, epsilon);
    ASSERT_NEAR(11.04533, v3, epsilon);

    eqs1.reset();
    eqs2.reset();
    eqs3.reset();
    timestepping.solve(5.0);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("a"))->eval(0., v1);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("b"))->eval(0., v2);
    part2.getOutput<MyFunction>(part2.getOutputParameterID("c"))->eval(0., v3);
    ASSERT_NEAR(12.2132, v1, epsilon);
    ASSERT_NEAR(3.38, v2, epsilon);
    ASSERT_NEAR(11.04533, v3, epsilon);
}


