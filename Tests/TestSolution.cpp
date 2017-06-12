/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TestSolution.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include <gtest/gtest.h>

#include <vector>

#include "BaseLib/CodingTools.h"

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/LinearEquation/LisLinearEquation.h"
#include "MathLib/Nonlinear/NRIterationStepInitializerDummy.h"
#include "GeoLib/Rectangle.h"
#include "GeoLib/GeoDomain.h"

#include "MeshLib/Tools/MeshGenerator.h"
#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/Serial/DiscreteSystem.h"
#include "FemLib/Tools/FemElementObjectContainer.h"

#include "NumLib/Function/TXFunction.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerResidualLocalAssembler.h"
#include "NumLib/Nonlinear/TemplateDiscreteNonlinearSolver.h"
#include "NumLib/Nonlinear/DiscreteNRSolverWithStepInitFactory.h"

#include "SolutionLib/Fem/FemDirichletBC.h"
#include "SolutionLib/Fem/FemNeumannBC.h"
#include "SolutionLib/Fem/FemIVBVProblem.h"
#include "SolutionLib/Fem/SingleStepFEM.h"

#include "TestUtil.h"

using namespace GeoLib;
using namespace MathLib;
using namespace MeshLib;
using namespace FemLib;
using namespace NumLib;
using namespace SolutionLib;
using namespace DiscreteLib;

template <class T>
class GWTimeODEAssembler: public T
{
public:
    typedef MathLib::LocalVector LocalVectorType;
    typedef MathLib::LocalMatrix LocalMatrixType;

    GWTimeODEAssembler(FemLib::LagrangeFeObjectContainer &feObjects, NumLib::ITXFunction &mat)
    : _matK(&mat), _feObjects(&feObjects)
    {
    };

    virtual ~GWTimeODEAssembler() {};

protected:
    virtual void assembleODE(const NumLib::TimeStep &/*time*/, const MeshLib::IElement &e, const LocalVectorType &/*u1*/, const LocalVectorType &/*u0*/, LocalMatrixType &/*localM*/, LocalMatrixType &localK, LocalVectorType &/*localF*/)
    {
        IFiniteElement* fe = _feObjects->getFeObject(e);

        //localM = .0;
        //localF.resize(localF.size(), .0);
        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
            q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);
            LocalMatrixType k;
            _matK->eval(real_x, k);
            fe->integrateDWxDN(j, k, localK);
        }
    }
private:
    NumLib::ITXFunction* _matK;
    FemLib::LagrangeFeObjectContainer* _feObjects;
};

class GWJacobianAssembler : public NumLib::IElementWiseTransientJacobianLocalAssembler
{
private:
    NumLib::ITXFunction* _matK;
    FemLib::LagrangeFeObjectContainer* _feObjects;
    typedef MathLib::LocalMatrix LocalMatrixType;
    typedef MathLib::LocalVector LocalVectorType;
public:
    GWJacobianAssembler(FemLib::LagrangeFeObjectContainer &feObjects, NumLib::ITXFunction &mat)
    : _matK(&mat), _feObjects(&feObjects), test(.0)
    {
    };

    /// assemble a local Jacobian matrix for the given element
    /// @param time            time step
    /// @param e            element
    /// @param local_u_n1    guess of current time step value
    /// @param local_u_n    previous time step value
    /// @param local_J        local Jacobian
    virtual void assembly(const TimeStep &/*time*/,  const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &/*localDofManager*/, const LocalVectorType &local_u_n1, const LocalVectorType &/*local_u_n*/, LocalMatrixType &local_J)
    {
        IFiniteElement* fe = _feObjects->getFeObject(e);
        const size_t n = local_u_n1.size();

        //const double dt = time.getTimeStepSize();
        const double theta = 1.0;
        LocalVectorType localX0;
        LocalMatrixType localM = LocalMatrixType::Zero(n,n);
        LocalMatrixType localK = LocalMatrixType::Zero(n,n);
        LocalVectorType localF;
        //localF.resize(localF.size(), .0);
        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
            q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);
            MathLib::LocalMatrix k;
            _matK->eval(real_x, k);
            fe->integrateDWxDN(j, k, localK);
        }

        //localR = (1/dt * localM + theta * localK) * localX - (1/dt * localM + (1-theta) * localK) * localX0 - localF;
        LocalMatrixType tmpM(n, n);
        tmpM = localK;
        tmpM *= theta;
        //tmpM.axpy(1.0, &localX[0], 0.0, &localR[0]);

        local_J = localK;
    }

    void setTest(double v) {test = v;};
private:
    double test;
};

class NRIterationStepInitializerTest
{
public:
    /// Default setting
    NRIterationStepInitializerTest(
            GWTimeODEAssembler<ElementWiseTimeEulerResidualLocalAssembler>* assemblerR, 
            GWJacobianAssembler *assemblerJ) 
            : _assemblerR(assemblerR), _assemblerJ(assemblerJ) {};

    /**
     * do pre processing
     *
     * \param dx        solution increment
     * \param x_new     new solution
     */
    template<class T_X, class F_RESIDUALS, class F_DX>
    void pre_process(const T_X &/*dx*/, const T_X &/*x_new*/, F_RESIDUALS &/*f_residuals*/, F_DX &/*f_dx*/)
    {
        _assemblerJ->setTest(1.0);
    }

    /**
     * do post processing
     * @param
     * @param
     * @param
     * @param
     */
    template<class T_X, class F_RESIDUALS, class F_DX>
    void post_process(const T_X &/*dx*/, const T_X &/*x_new*/, F_RESIDUALS &/*f_residuals*/, F_DX &/*f_dx*/)
    {
        _assemblerJ->setTest(1.0);
    }

private:
    GWTimeODEAssembler<ElementWiseTimeEulerResidualLocalAssembler>* _assemblerR;
    GWJacobianAssembler* _assemblerJ;
};

template <
    class T_LINEAR_SOLVER
    >
class GWFemTestSystem : public NumLib::ITransientSystem
{
    typedef TemplateFemEquation
            <
                DiscreteSystem,
                T_LINEAR_SOLVER,
                GWTimeODEAssembler<ElementWiseTimeEulerEQSLocalAssembler>,
                GWTimeODEAssembler<ElementWiseTimeEulerResidualLocalAssembler>,
                GWJacobianAssembler
            > GWFemEquation;

    typedef FemIVBVProblem<DiscreteSystem,GWFemEquation> GWFemProblem;

    typedef typename GWFemEquation::LinearEQSType MyLinearFunction;
    typedef typename GWFemEquation::ResidualEQSType MyResidualFunction;
    typedef typename GWFemEquation::DxEQSType MyDxFunction;

    typedef TemplateDiscreteNonlinearSolver
            <
                DiscreteSystem,
                MyLinearFunction, 
                MyResidualFunction, 
                MyDxFunction
            > MyNonlinearFunction;

    typedef SingleStepFEM
            <
                GWFemProblem,
                T_LINEAR_SOLVER,
                NumLib::DiscreteNRSolverWithStepInitFactory
                <
                    NRIterationStepInitializerTest
                >
            > SolutionForHead;

public:
    GWFemTestSystem()
    {
        BaseLib::zeroObject(_rec, _head, _problem, _solHead, _feObjects);
    }

    virtual ~GWFemTestSystem()
    {
        BaseLib::releaseObject(_rec, _head, _problem);
    }

    double suggestNext(const TimeStep &time_current)
    {
        return _solHead->suggestNext(time_current);
    }
    bool isAwake(const TimeStep &time)
    {
        return _solHead->isAwake(time);
    }

    int solveTimeStep(const TimeStep &time)
    {
        _solHead->solveTimeStep(time);
        _head = _solHead->getCurrentSolution(0);
        return 0;
    }

    void accept(const TimeStep &time) 
    {
        _solHead->accept(time);
    };

    //#Define a problem
    void define(DiscreteSystem &dis, NumLib::ITXFunction &K, BaseLib::Options &option)
    {
        MeshLib::IMesh *msh = dis.getMesh();
        //size_t nnodes = msh->getNumberOfNodes();
        _feObjects = new LagrangeFeObjectContainer(msh);
        //equations
        typename GWFemEquation::LinearAssemblerType* local_linear = new typename GWFemEquation::LinearAssemblerType(*_feObjects, K);
        typename GWFemEquation::ResidualAssemblerType* local_r = new typename GWFemEquation::ResidualAssemblerType(*_feObjects, K);
        typename GWFemEquation::JacobianAssemblerType* local_J = new typename GWFemEquation::JacobianAssemblerType(*_feObjects, K);
        //IVBV problem
        _problem = new GWFemProblem(&dis);
        GWFemEquation* eqs = _problem->createEquation();
        eqs->initialize(local_linear, local_r, local_J);
        // var
        typename GWFemProblem::MyVariable* _head = _problem->addVariable("head");
        //IC
        SolutionLib::FemIC* var_ic = new SolutionLib::FemIC(dis.getMesh());
        var_ic->addDistribution(new GeoLib::GeoDomain(), new  NumLib::TXFunctionConstant(.0));
        _head->setIC(var_ic);
//        typename GWFemProblem::MyVariable::MyNodalFunctionScalar* h0 = new typename GWFemProblem::MyVariable::MyNodalFunctionScalar();
//        h0->initialize(dis, PolynomialOrder::Linear, .0);
//        _head->setIC(h0);
        //BC
        _rec = new GeoLib::Rectangle(Point(0.0, 0.0, 0.0),  Point(2.0, 2.0, 0.0));
        const GeoLib::Polyline &poly_left = _rec->getLeft();
        const GeoLib::Polyline &poly_right = _rec->getRight();
        //_head->setIC();
        _head->addDirichletBC(new SolutionLib::FemDirichletBC(msh, &poly_right,  new NumLib::TXFunctionConstant(.0)));
        _head->addNeumannBC(new SolutionLib::FemNeumannBC(msh, _feObjects, &poly_left, new NumLib::TXFunctionConstant(-1e-5)));
        //transient
        TimeStepFunctionConstant tim(.0, 100.0, 10.0);
        _problem->setTimeSteppingFunction(tim);
        //solution algorithm
        //MathLib::NRIterationStepInitializerDummy* nl_step_init(new MathLib::NRIterationStepInitializerDummy());
        NRIterationStepInitializerTest* nl_step_init(new NRIterationStepInitializerTest(local_r, local_J));
        DiscreteNRSolverWithStepInitFactory<NRIterationStepInitializerTest>* nl_factory(new DiscreteNRSolverWithStepInitFactory<NRIterationStepInitializerTest>(nl_step_init));
        _solHead = new SolutionForHead(&dis, _problem, nl_factory);
        //_solHead->getTimeODEAssembler()->setTheta(1.0);
        typename SolutionForHead::LinearSolverType* linear_solver = _solHead->getLinearEquationSolver();
        linear_solver->setOption(option);
        typename SolutionForHead::NonlinearSolverType* nonlinear_solver = _solHead->getNonlinearSolver();
        nonlinear_solver->setOption(option);


        //vel = new FEMIntegrationPointFunctionVector2d(msh);
    }

    FemLib::FemNodalFunctionScalar<DiscreteSystem>::type* getCurrentHead()
    {
        return _head;
    }

private:
    GWFemProblem* _problem;
    SolutionForHead* _solHead; 
    GeoLib::Rectangle *_rec;
    FemLib::FemNodalFunctionScalar<DiscreteSystem>::type *_head;
    LagrangeFeObjectContainer* _feObjects;

    DISALLOW_COPY_AND_ASSIGN(GWFemTestSystem);
};





TEST(Solution, Fem1_Linear)
{
    // create a discrete system
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    DiscreteSystem dis(msh);
    // mat
    NumLib::TXFunctionConstant K(1.e-11);
    // options
    BaseLib::Options options;
    BaseLib::Options* op_lis = options.addSubGroup("LinearSolver");
    op_lis->addOption("solver_type", "CG");
    op_lis->addOption("precon_type", "NONE");
    op_lis->addOptionAsNum("error_tolerance", 1e-10);
    op_lis->addOptionAsNum("max_iteration_step", 500);
    // define problems
    GWFemTestSystem<MathLib::LisLinearEquation> gwProblem;
    gwProblem.define(dis, K, options);

    gwProblem.solveTimeStep(TimeStep(1.0, 1.0));

    std::vector<double> expected(9);
    for (size_t i=0; i<9; i++) {
        if (i%3==0) expected[i] = 2.e+6;
        if (i%3==1) expected[i] = 1.e+6;
        if (i%3==2) expected[i] = 0.e+6;
    }

    FemLib::FemNodalFunctionScalar<DiscreteSystem>::type* h = gwProblem.getCurrentHead();

    ASSERT_DOUBLE_ARRAY_EQ(&expected[0], &(*h->getDiscreteData())[0], h->getNumberOfNodes());

    BaseLib::releaseObject(msh);
}

TEST(Solution, Fem1_Picard)
{
    // create a discrete system
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    DiscreteSystem dis(msh);
    // mat
    NumLib::TXFunctionConstant K(1.e-11);
    // options
    BaseLib::Options options;
    BaseLib::Options* op_lis = options.addSubGroup("LinearSolver");
    op_lis->addOption("solver_type", "CG");
    op_lis->addOption("precon_type", "NONE");
    op_lis->addOptionAsNum("error_tolerance", 1e-10);
    op_lis->addOptionAsNum("max_iteration_step", 500);
    BaseLib::Options* op_nl = options.addSubGroup("Nonlinear");
    op_nl->addOption("solver_type", "Picard");
    op_nl->addOptionAsNum("error_tolerance", 1e-6);
    op_nl->addOptionAsNum("max_iteration_step", 500);
    // define problems
    GWFemTestSystem<MathLib::LisLinearEquation> gwProblem;
    gwProblem.define(dis, K, options);

    gwProblem.solveTimeStep(TimeStep(1.0, 1.0));

    std::vector<double> expected(9);
    for (size_t i=0; i<9; i++) {
        if (i%3==0) expected[i] = 2.e+6;
        if (i%3==1) expected[i] = 1.e+6;
        if (i%3==2) expected[i] = 0.e+6;
    }

    FemLib::FemNodalFunctionScalar<DiscreteSystem>::type* h = gwProblem.getCurrentHead();

    ASSERT_DOUBLE_ARRAY_EQ(&expected[0], &(*h->getDiscreteData())[0], h->getNumberOfNodes());

    BaseLib::releaseObject(msh);
}

TEST(Solution, Fem1_Newton)
{
    // create a discrete system
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    DiscreteSystem dis(msh);
    // mat
    NumLib::TXFunctionConstant K(1.e-11);
    // options
    BaseLib::Options options;
    BaseLib::Options* op_lis = options.addSubGroup("LinearSolver");
    op_lis->addOption("solver_type", "CG");
    op_lis->addOption("precon_type", "NONE");
    op_lis->addOptionAsNum("error_tolerance", 1e-10);
    op_lis->addOptionAsNum("max_iteration_step", 500);
    BaseLib::Options* op_nl = options.addSubGroup("Nonlinear");
    op_nl->addOption("solver_type", "Newton");
    op_nl->addOptionAsNum("error_tolerance", 1e-6);
    op_nl->addOptionAsNum("max_iteration_step", 10);
    // define problems
    GWFemTestSystem<MathLib::LisLinearEquation> gwProblem;
    gwProblem.define(dis, K, options);

    gwProblem.solveTimeStep(TimeStep(1.0, 1.0));

    std::vector<double> expected(9);
    for (size_t i=0; i<9; i++) {
        if (i%3==0) expected[i] = 2.e+6;
        if (i%3==1) expected[i] = 1.e+6;
        if (i%3==2) expected[i] = 0.e+6;
    }

    FemLib::FemNodalFunctionScalar<DiscreteSystem>::type* h = gwProblem.getCurrentHead();

    ASSERT_DOUBLE_ARRAY_EQ(&expected[0], &(*h->getDiscreteData())[0], h->getNumberOfNodes());

    BaseLib::releaseObject(msh);
}

TEST(Solution, Fem2)
{
    // create a discrete system
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    DiscreteSystem dis(msh);
    // mat
    NumLib::TXFunctionConstant K(1.e-11);
    // options
    BaseLib::Options options;
    BaseLib::Options* op_lis = options.addSubGroup("LinearSolver");
    op_lis->addOption("solver_type", "CG");
    op_lis->addOption("precon_type", "NONE");
    op_lis->addOptionAsNum("error_tolerance", 1e-10);
    op_lis->addOptionAsNum("max_iteration_step", 500);
    // define problems
    GWFemTestSystem<MathLib::LisLinearEquation> gwProblem;
    gwProblem.define(dis, K, options);

    // start time stepping
    TimeSteppingController timeStepping;
    timeStepping.setTransientSystem(gwProblem);
    timeStepping.setBeginning(.0);
    timeStepping.solve(30.);


    std::vector<double> expected(9);
    for (size_t i=0; i<9; i++) {
        if (i%3==0) expected[i] = 2.e+6;
        if (i%3==1) expected[i] = 1.e+6;
        if (i%3==2) expected[i] = 0.e+6;
    }

    FemLib::FemNodalFunctionScalar<DiscreteSystem>::type* h = gwProblem.getCurrentHead();

    ASSERT_DOUBLE_ARRAY_EQ(&expected[0], &(*h->getDiscreteData())[0], h->getNumberOfNodes());


    BaseLib::releaseObject(msh);
}

