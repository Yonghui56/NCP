/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Test.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include <cmath>
#include <gtest/gtest.h>

#include "BaseLib/CodingTools.h"

#include "MathLib/Vector.h"
#include "MathLib/LinAlg/LinearEquation/LisLinearEquation.h"
#include "GeoLib/Line.h"
#include "GeoLib/GeoDomain.h"
#include "MeshLib/Tools/MeshGenerator.h"
#include "DiscreteLib/Serial/DiscreteSystem.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TransientCoupling/AsyncPartitionedSystem.h"
#include "NumLib/Nonlinear/TemplateDiscreteNonlinearSolver.h"
#include "NumLib/Coupling/Algorithm/TransientPartitionedAlgorithm.h"
#include "NumLib/Function/NormOfDiscreteDataFunction.h"
#include "SolutionLib/Fem/FemNeumannBC.h"
#include "Tests/Geo/Model/Head.h"
#include "Tests/Geo/Model/Velocity.h"
#include "Tests/Geo/Model/Concentration.h"
#include "Tests/Geo/Model/Displacement.h"
#include "Tests/Geo/Model/StressStrain.h"

#include "Tests/ExactSolution/OgataBank.h"

#include "TestUtil.h"

typedef FemLib::FEMIntegrationPointFunctionVector<DiscreteSystem>::type MyIntegrationPointFunctionVector;
typedef FemLib::FemNodalFunctionScalar<DiscreteSystem>::type MyNodalFunctionScalar;
//typedef NumLib::ITXDiscreteFunction<double> MyNodalFunctionScalar;
//typedef NumLib::ITXDiscreteFunction<MathLib::TemplateVectorX<MathLib::LocalVector> > MyIntegrationPointFunctionVector;

//class DiscreteDataConvergenceCheck : public IConvergenceCheck
//{
//public:
//    explicit DiscreteDataConvergenceCheck(DiscreteLib::DiscreteSystem *dis)
//    {
//
//    }
//
//    bool isConverged(UnnamedParameterSet& vars_prev, UnnamedParameterSet& vars_current, double eps, double &v_diff)
//    {
//        for (size_t i=0; i<vars_prev.size(); i++) {
//#if 1
//            if (vars_prev.getName(i).compare("h")==0 || vars_prev.getName(i).compare("c")==0) {
//                const MyNodalFunctionScalar* f_fem_prev = vars_prev.get<MyNodalFunctionScalar>(i);
//                const MyNodalFunctionScalar* f_fem_cur = vars_current.get<MyNodalFunctionScalar>(i);
//                //v_diff = f_fem_cur->norm_diff(*f_fem_prev);
//                NumLib::NormOfDiscreteDataFunction<double> _norm;
//                v_diff = _norm(*f_fem_prev, *f_fem_cur);
//            } else if (vars_prev.getName(i).compare("v")==0) {
//                const MyIntegrationPointFunctionVector* f_fem_prev = vars_prev.get<MyIntegrationPointFunctionVector>(i);
//                const MyIntegrationPointFunctionVector* f_fem_cur = vars_current.get<MyIntegrationPointFunctionVector>(i);
//                //v_diff = f_fem_cur->norm_diff(*f_fem_prev);
//                NumLib::NormOfDiscreteDataFunction<MathLib::TemplateVectorX<MathLib::LocalVector> > _norm;
//                v_diff = _norm(*f_fem_prev, *f_fem_cur);
//            }
//#endif
//            if (v_diff>eps) {
//                return false;
//            }
//        }
//        return true;
//    }
//};

typedef Geo::FunctionHead<LisLinearEquation> MyFunctionHead;
typedef Geo::FunctionVelocity MyFunctionVelocity;
typedef Geo::FunctionConcentration<LisLinearEquation> MyFunctionConcentration;
typedef Geo::FunctionDisplacement<LisLinearEquation> MyFunctionDisplacement;
typedef Geo::FunctionStressStrain MyFunctionStressStrain;

namespace Geo
{
typedef MyFunctionHead::GWFemProblem GWFemProblem;
typedef MyFunctionHead::GWFemEquation GWFemEquation;
typedef MyFunctionConcentration::MassFemProblem MassFemProblem;
typedef MyFunctionConcentration::MassFemEquation MassFemEquation;
typedef MyFunctionDisplacement::FemLinearElasticProblem FemLinearElasticProblem;
typedef MyFunctionDisplacement::FemLinearElasticEquation FemLinearElasticEquation;
}

Geo::GWFemProblem* createGWProblem(DiscreteSystem &dis, Geo::PorousMedia &pm)
{
    LagrangeFeObjectContainer* _feObjects = new LagrangeFeObjectContainer(dis.getMesh());
    //equations
    Geo::GWFemEquation::LinearAssemblerType* linear_assembler = new Geo::GWFemEquation::LinearAssemblerType(*_feObjects, pm);
    Geo::GWFemEquation::ResidualAssemblerType* r_assembler = new Geo::GWFemEquation::ResidualAssemblerType(*_feObjects, pm);
    Geo::GWFemEquation::JacobianAssemblerType* j_eqs = new Geo::GWFemEquation::JacobianAssemblerType(*_feObjects, pm);
    //IVBV problem
    Geo::GWFemProblem* _problem = new Geo::GWFemProblem(&dis);
    Geo::GWFemEquation* eqs = _problem->createEquation();
    eqs->initialize(linear_assembler, r_assembler, j_eqs);
    return _problem;
}


Geo::GWFemProblem* defineGWProblem(DiscreteSystem &dis, GeoLib::Rectangle &_rec, Geo::PorousMedia &pm, FemLib::LagrangeFeObjectContainer* feObjects)
{
    Geo::GWFemProblem* _problem = createGWProblem(dis, pm);
    // var
    Geo::GWFemProblem::MyVariable* head = _problem->addVariable("head");
    // IC
//    MyNodalFunctionScalar* h0 = new MyNodalFunctionScalar();
//    h0->initialize(dis, PolynomialOrder::Linear, 0);
//    h0->setFeObjectContainer(feObjects);
    SolutionLib::FemIC* var_ic = new SolutionLib::FemIC(dis.getMesh());
    var_ic->addDistribution(new GeoLib::GeoDomain(), new  NumLib::TXFunctionConstant(.0));
    head->setIC(var_ic);
    //BC
    const GeoLib::Polyline &poly_left = _rec.getLeft();
    const GeoLib::Polyline &poly_right = _rec.getRight();
    head->addDirichletBC(new FemDirichletBC(dis.getMesh(), &poly_right, new NumLib::TXFunctionConstant(.0)));
    head->addNeumannBC(new FemNeumannBC(dis.getMesh(), feObjects, &poly_left, new NumLib::TXFunctionConstant(-1e-5)));

    return _problem;
}

Geo::GWFemProblem* defineGWProblem1D(DiscreteSystem &dis, GeoLib::Line &line, Geo::PorousMedia &pm, FemLib::LagrangeFeObjectContainer* feObjects, double p_out, double q_in)
{
    Geo::GWFemProblem* _problem = createGWProblem(dis, pm);
    // var
    Geo::GWFemProblem::MyVariable* head = _problem->addVariable("head");
    // IC
//    MyNodalFunctionScalar* h0 = new MyNodalFunctionScalar();
//    h0->initialize(dis, PolynomialOrder::Linear, 0);
//    h0->setFeObjectContainer(feObjects);
//    head->setIC(h0);
    SolutionLib::FemIC* var_ic = new SolutionLib::FemIC(dis.getMesh());
    var_ic->addDistribution(new GeoLib::GeoDomain(), new  NumLib::TXFunctionConstant(.0));
    head->setIC(var_ic);
    //BC
    head->addDirichletBC(new FemDirichletBC(dis.getMesh(), &line.getPoint2(), new NumLib::TXFunctionConstant(p_out)));
    head->addNeumannBC(new FemNeumannBC(dis.getMesh(), feObjects, &line.getPoint1(), new NumLib::TXFunctionConstant(-q_in)));

    return _problem;
}

Geo::MassFemProblem* defineMassTransportProblem(DiscreteSystem &dis, GeoLib::Rectangle &_rec, Geo::PorousMedia &pm, Geo::Compound &comp, FemLib::LagrangeFeObjectContainer* /*feObjects*/)
{
    LagrangeFeObjectContainer* _feObjects = new LagrangeFeObjectContainer(dis.getMesh());
    //equations
    Geo::MassFemEquation::LinearAssemblerType* linear_assembler = new Geo::MassFemEquation::LinearAssemblerType(*_feObjects, pm, comp);
    Geo::MassFemEquation::ResidualAssemblerType* r_assembler = new Geo::MassFemEquation::ResidualAssemblerType(*_feObjects, pm, comp);
    Geo::MassFemEquation::JacobianAssemblerType* j_eqs = new Geo::MassFemEquation::JacobianAssemblerType(*_feObjects, pm, comp);
    //IVBV problem
    Geo::MassFemProblem* _problem = new Geo::MassFemProblem(&dis);
    Geo::MassFemEquation* eqs = _problem->createEquation();
    eqs->initialize(linear_assembler, r_assembler, j_eqs);
    // var
    Geo::MassFemProblem::MyVariable* c = _problem->addVariable("c");
    // IC
//    MyNodalFunctionScalar* c0 = new MyNodalFunctionScalar();
//    c0->initialize(dis, PolynomialOrder::Linear, 0);
//    c0->setFeObjectContainer(feObjects);
//    c->setIC(c0);
    SolutionLib::FemIC* var_ic = new SolutionLib::FemIC(dis.getMesh());
    var_ic->addDistribution(new GeoLib::GeoDomain(), new  NumLib::TXFunctionConstant(.0));
    c->setIC(var_ic);
    //BC
    const GeoLib::Polyline &poly_left = _rec.getLeft();
    c->addDirichletBC(new FemDirichletBC(dis.getMesh(), &poly_left, new NumLib::TXFunctionConstant(1.0)));

    return _problem;
}

Geo::FemLinearElasticProblem* defineLinearElasticProblem(DiscreteSystem &dis, GeoLib::Rectangle &_rec, Geo::PorousMedia &pm, FemLib::LagrangeFeObjectContainer* /*feObjects*/)
{
    LagrangeFeObjectContainer* _feObjects = new LagrangeFeObjectContainer(dis.getMesh());
    //equations
    Geo::FemLinearElasticEquation::LinearAssemblerType* linear_assembler = new Geo::FemLinearElasticEquation::LinearAssemblerType(*_feObjects, pm);
    Geo::FemLinearElasticEquation::ResidualAssemblerType* r_assembler = new Geo::FemLinearElasticEquation::ResidualAssemblerType(*_feObjects, pm);
    Geo::FemLinearElasticEquation::JacobianAssemblerType* j_eqs = new Geo::FemLinearElasticEquation::JacobianAssemblerType(*_feObjects, pm);
    //IVBV problem
    Geo::FemLinearElasticProblem* _problem = new Geo::FemLinearElasticProblem(&dis);
    Geo::FemLinearElasticEquation* eqs = _problem->createEquation();
    eqs->initialize(linear_assembler, r_assembler, j_eqs);
    // var
    Geo::FemLinearElasticProblem::MyVariable* u_x = _problem->addVariable("u_x");
    Geo::FemLinearElasticProblem::MyVariable* u_y = _problem->addVariable("u_y");
    // IC
//    MyNodalFunctionScalar* u0 = new MyNodalFunctionScalar();
//    u0->initialize(dis, PolynomialOrder::Linear, 0);
//    u0->setFeObjectContainer(feObjects);
//    u_x->setIC(u0);
//    u_y->setIC(u0);
    SolutionLib::FemIC* var_ic = new SolutionLib::FemIC(dis.getMesh());
    var_ic->addDistribution(new GeoLib::GeoDomain(), new  NumLib::TXFunctionConstant(.0));
    u_x->setIC(var_ic);
    u_y->setIC(var_ic);
    //BC
    const GeoLib::Polyline &poly_bottom = _rec.getBottom();
    u_y->addDirichletBC(new FemDirichletBC(dis.getMesh(), &poly_bottom, new NumLib::TXFunctionConstant(.0)));
    u_y->addNeumannBC(new SolutionLib::FemNeumannBC(dis.getMesh(), _feObjects, &_rec.getTop(), new NumLib::TXFunctionConstant((-1e+6)*(-1.))));


    return _problem;
}

static void getGWExpectedHead(std::vector<double> &expected)
{
    expected.resize(9);
    for (size_t i=0; i<9; i++) {
        if (i%3==0) expected[i] = 2.e+6;
        if (i%3==1) expected[i] = 1.e+6;
        if (i%3==2) expected[i] = 0.e+6;
    }
}



TEST(Solution, CouplingFem2D)
{
    try {
        MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
        GeoLib::Rectangle* _rec = new GeoLib::Rectangle(Point(0.0, 0.0, 0.0),  Point(2.0, 2.0, 0.0));
        Geo::PorousMedia pm;
        pm.hydraulic_conductivity = new NumLib::TXFunctionConstant(1.e-11);
        pm.porosity = new NumLib::TXFunctionConstant(0.2);
        DiscreteSystem dis(msh);
        FemLib::LagrangeFeObjectContainer feObjects(msh);
        Geo::GWFemProblem* pGW = defineGWProblem(dis, *_rec, pm, &feObjects);
        TimeStepFunctionConstant tim(.0, 100.0, 10.0);
        pGW->setTimeSteppingFunction(tim);
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

        MyFunctionHead f_head;
        f_head.define(&dis, pGW, options);
        f_head.setOutputParameterName(0, "h");
        MyFunctionVelocity f_vel;
        f_vel.define(dis, pm);
        f_vel.setInputParameterName(0, "h");
        f_vel.setOutputParameterName(0, "v");

        SerialStaggeredMethod method(1e-5, 100);
        AsyncPartitionedSystem apart1;
        apart1.setAlgorithm(method);
        apart1.resizeOutputParameter(2);
        apart1.setOutputParameterName(0, "h");
        apart1.setOutputParameterName(1, "v");
        apart1.addProblem(f_head);
        apart1.addProblem(f_vel);
        apart1.connectParameters();

        TimeSteppingController timestepping;
        timestepping.setTransientSystem(apart1);

        //const double epsilon = 1.e-3;
        timestepping.setBeginning(.0);
        timestepping.solve(tim.getEnd());

        const MyNodalFunctionScalar* r_f_head = apart1.getOutput<MyNodalFunctionScalar>(apart1.getOutputParameterID("h"));
        const MyIntegrationPointFunctionVector* r_f_v = apart1.getOutput<MyIntegrationPointFunctionVector>(apart1.getOutputParameterID("v"));
        const IDiscreteVector<double>* vec_h = r_f_head->getDiscreteData();
        //const FEMIntegrationPointFunctionVector2d::DiscreteVectorType* vec_v = r_f_v->getNodalValues();

        r_f_head->printout();
        r_f_v->printout();

        std::vector<double> expected;
        getGWExpectedHead(expected);

        ASSERT_DOUBLE_ARRAY_EQ(&expected[0], &(*vec_h)[0], vec_h->size());

    } catch (const char* e) {
        std::cout << "***Exception caught! " << e << std::endl;
    }

}

TEST(Solution, line)
{
    try {
        const double len = 2.0;
        const size_t div = 2;
        const double h = len / div;
        const double p_out = .0;
        const double q_in = 1e-5;
        const double k = 1e-11;
        const double poro = 0.2;
        MeshLib::IMesh *msh = MeshGenerator::generateLineMesh(len, div, .0, .0, .0);
        GeoLib::Line* line = new GeoLib::Line(Point(0.0, 0.0, 0.0),  Point(len, 0.0, 0.0));
        Geo::PorousMedia pm;
        pm.hydraulic_conductivity = new NumLib::TXFunctionConstant(k);
        pm.porosity = new NumLib::TXFunctionConstant(poro);
        DiscreteSystem dis(msh);
        FemLib::LagrangeFeObjectContainer feObjects(msh);
        Geo::GWFemProblem* pGW = defineGWProblem1D(dis, *line, pm, &feObjects, p_out, q_in);
        TimeStepFunctionConstant tim(.0, 10.0, 10.0);
        pGW->setTimeSteppingFunction(tim);
        // options
        BaseLib::Options options;
        BaseLib::Options* op_lis = options.addSubGroup("LinearSolver");
        op_lis->addOption("solver_type", "CG");
        op_lis->addOption("precon_type", "NONE");
        op_lis->addOptionAsNum("error_tolerance", 1e-10);
        op_lis->addOptionAsNum("max_iteration_step", 500);

        MyFunctionHead f_head;
        f_head.define(&dis, pGW, options);
        f_head.setOutputParameterName(0, "h");
        MyFunctionVelocity f_vel;
        f_vel.define(dis, pm);
        f_vel.setInputParameterName(0, "h");
        f_vel.setOutputParameterName(0, "v");

        SerialStaggeredMethod method(1e-5, 100);
        AsyncPartitionedSystem apart1;
        apart1.setAlgorithm(method);
        apart1.resizeOutputParameter(2);
        apart1.setOutputParameterName(0, "h");
        apart1.setOutputParameterName(1, "v");
        apart1.addProblem(f_head);
        apart1.addProblem(f_vel);
        apart1.connectParameters();

        TimeSteppingController timestepping;
        timestepping.setTransientSystem(apart1);

        //const double epsilon = 1.e-3;
        timestepping.setBeginning(.0);
        timestepping.solve(tim.getEnd());

        const MyNodalFunctionScalar* r_f_head = apart1.getOutput<MyNodalFunctionScalar>(apart1.getOutputParameterID("h"));
        const MyIntegrationPointFunctionVector* r_f_v = apart1.getOutput<MyIntegrationPointFunctionVector>(apart1.getOutputParameterID("v"));
        const IDiscreteVector<double>* vec_h = r_f_head->getDiscreteData();
        //const FEMIntegrationPointFunctionVector2d::DiscreteVectorType* vec_v = r_f_v->getNodalValues();

        r_f_head->printout();
        r_f_v->printout();

        std::vector<double> expected;
        expected.resize(div+1);
        //const double p_left = 2.e+6;
        const double p_left = p_out + q_in * len / k;
        const double p_right = p_out;
        const double x_len = len;
        for (size_t i=0; i<expected.size(); i++) {
            double x = i*h;
            expected[i] = (p_right-p_left) / x_len * x + p_left;
        }

        ASSERT_DOUBLE_ARRAY_EQ(&expected[0], &(*vec_h)[0], vec_h->size());

    } catch (const char* e) {
        std::cout << "***Exception caught! " << e << std::endl;
    }

}

TEST(Solution, CouplingFem2)
{
    // problem definition
    const size_t div = 20;
    const double poro = 1.0;
    const double mol_diff = 1e-6;
    const double darcy_vel = 1e-5;

    try {
        //space
        MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, div, .0, .0, .0);
        GeoLib::Rectangle* _rec = new GeoLib::Rectangle(Point(0.0, 0.0, 0.0),  Point(2.0, 2.0, 0.0));
        //time
        //TimeStepFunctionConstant tim(.0, 1e+3, 1e+3);
        TimeStepFunctionConstant tim(.0, 1e+4, 1e+3);
        //material
        Geo::PorousMedia pm;
        pm.hydraulic_conductivity = new NumLib::TXFunctionConstant(1.e-11);
        pm.porosity = new NumLib::TXFunctionConstant(poro);
        Geo::Compound tracer;
        tracer.molecular_diffusion = new NumLib::TXFunctionConstant(mol_diff);
        //problems
        DiscreteSystem dis(msh);
        FemLib::LagrangeFeObjectContainer feObjects(msh);
        Geo::GWFemProblem* pGW = defineGWProblem(dis, *_rec, pm, &feObjects);
        Geo::MassFemProblem* pMass = defineMassTransportProblem(dis, *_rec, pm, tracer, &feObjects);
        pGW->setTimeSteppingFunction(tim);
        pMass->setTimeSteppingFunction(tim);
        //options
        BaseLib::Options optionsGW;
        BaseLib::Options* op_lis = optionsGW.addSubGroup("LinearSolver");
        op_lis->addOption("solver_type", "CG");
        op_lis->addOption("precon_type", "NONE");
        op_lis->addOptionAsNum("error_tolerance", 1e-10);
        op_lis->addOptionAsNum("max_iteration_step", 1000);
        BaseLib::Options* op_nl = optionsGW.addSubGroup("Nonlinear");
        op_nl->addOption("solver_type", "Picard");
        op_nl->addOptionAsNum("error_tolerance", 1e-6);
        op_nl->addOptionAsNum("max_iteration_step", 500);
        BaseLib::Options optionsMT;
        op_lis = optionsMT.addSubGroup("LinearSolver");
        op_lis->addOption("solver_type", "BICG");
        op_lis->addOption("precon_type", "NONE");
        op_lis->addOptionAsNum("error_tolerance", 1e-10);
        op_lis->addOptionAsNum("max_iteration_step", 1000);
        BaseLib::Options* op_nl2 = optionsMT.addSubGroup("Nonlinear");
        op_nl2->addOption("solver_type", "Picard");
        op_nl2->addOptionAsNum("error_tolerance", 1e-6);
        op_nl2->addOptionAsNum("max_iteration_step", 500);
        // unknowns
        MyFunctionHead f_head;
        f_head.define(&dis, pGW, optionsGW);
        f_head.setOutputParameterName(0, "h");
        MyFunctionVelocity f_vel;
        f_vel.define(dis, pm);
        f_vel.setInputParameterName(0, "h");
        f_vel.setOutputParameterName(0, "v");
        MyFunctionConcentration f_c;
        f_c.define(&dis, pMass, optionsMT);
        f_c.setInputParameterName(0, "v");
        f_c.setOutputParameterName(0, "c");


        SerialStaggeredMethod method(1e-5, 100);
        AsyncPartitionedSystem apart1;
        apart1.setAlgorithm(method);
        apart1.resizeOutputParameter(3);
        apart1.setOutputParameterName(0, "h");
        apart1.setOutputParameterName(1, "v");
        apart1.setOutputParameterName(2, "c");
        apart1.addProblem(f_head);
        apart1.addProblem(f_vel);
        apart1.addProblem(f_c);
        apart1.connectParameters();

        TimeSteppingController timestepping;
        timestepping.setTransientSystem(apart1);

        //const double epsilon = 1.e-3;
        timestepping.setBeginning(.0);
        timestepping.solve(tim.getEnd());

        const MyNodalFunctionScalar* r_f_head = apart1.getOutput<MyNodalFunctionScalar>(apart1.getOutputParameterID("h"));
        //const MyIntegrationPointFunctionVector* r_f_v = apart1.getOutput<MyIntegrationPointFunctionVector>(apart1.getOutputParameterID("v"));
        const MyNodalFunctionScalar* r_f_c = apart1.getOutput<MyNodalFunctionScalar>(apart1.getOutputParameterID("c"));
        const IDiscreteVector<double>* vec_h = r_f_head->getDiscreteData();
        //const FEMIntegrationPointFunctionVector2d::DiscreteVectorType* vec_v = r_f_v->getNodalValues();
        const IDiscreteVector<double>* vec_c = r_f_c->getDiscreteData();

        //r_f_head->printout();
        //r_f_v->printout();
//#undef OUTPUT_C
#define OUTPUT_C
        std::vector<double> expectedHead(msh->getNumberOfNodes());
        const double p_left = 2.e+6;
        const double p_right = .0;
        const double x_len = 2.0;
#ifdef OUTPUT_C
        std::cout << std::endl << "expected p=";
#endif
        for (size_t i=0; i<expectedHead.size(); i++) {
            double x = msh->getNodeCoordinatesRef(i)->getData()[0];
            expectedHead[i] = (p_right-p_left) / x_len * x + p_left;
#ifdef OUTPUT_C
            std::cout << expectedHead[i] << " ";
#endif
        }
//        getGWExpectedHead(expectedHead);
        ASSERT_DOUBLE_ARRAY_EQ(&expectedHead[0], &(*vec_h)[0], expectedHead.size(), 1e-5);


        std::vector<double> expectedC(msh->getNumberOfNodes());

#ifdef OUTPUT_C
        std::cout << std::endl << "simulated C:"<< std::endl;
        r_f_c->printout();
        std::cout << "expected C=";
#endif
        for (size_t i=0; i<msh->getNumberOfNodes(); i++) {
            double x = msh->getNodeCoordinatesRef(i)->getData()[0];
            expectedC[i] = analyticalOgataBank(x, tim.getEnd(), darcy_vel/poro, mol_diff);
#ifdef OUTPUT_C
            std::cout << expectedC[i] << " ";
#endif
        }
#ifdef OUTPUT_C
        std::cout << std::endl;
#endif
        ASSERT_DOUBLE_ARRAY_EQ(&expectedC[0], &(*vec_c)[0], expectedC.size(), 5e-2);


    } catch (const char* e) {
        std::cout << "***Exception caught! " << e << std::endl;
    }

}

#if 1
TEST(Solution, LinearElastic2D)
{
    try {
        MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
        GeoLib::Rectangle* _rec = new GeoLib::Rectangle(Point(0.0, 0.0, 0.0),  Point(2.0, 2.0, 0.0));
        Geo::PorousMedia pm;
        pm.hydraulic_conductivity = new NumLib::TXFunctionConstant(1.e-11);
        pm.porosity = new NumLib::TXFunctionConstant(0.2);
        Geo::Solid solid;
        solid.density = 3e+3;
        solid.poisson_ratio = 0.2;
        solid.Youngs_modulus = 10e+9;
        pm.solidphase = &solid;
        DiscreteSystem dis(msh);
        FemLib::LagrangeFeObjectContainer feObjects(msh);
        Geo::FemLinearElasticProblem* pDe = defineLinearElasticProblem(dis, *_rec, pm, &feObjects);
        TimeStepFunctionConstant tim(.0, 10.0, 10.0);
        pDe->setTimeSteppingFunction(tim);
        // options
        BaseLib::Options options;
        BaseLib::Options* op_lis = options.addSubGroup("LinearSolver");
        op_lis->addOption("solver_type", "CG");
        op_lis->addOption("precon_type", "NONE");
        op_lis->addOptionAsNum("error_tolerance", 1e-10);
        op_lis->addOptionAsNum("max_iteration_step", 500);
        BaseLib::Options* op_nl = options.addSubGroup("Nonlinear");
        op_nl->addOption("solver_type", "Linear");
        op_nl->addOptionAsNum("error_tolerance", 1e-6);
        op_nl->addOptionAsNum("max_iteration_step", 500);

        MyFunctionDisplacement f_u;
        f_u.define(&dis, pDe, &pm, options);
        f_u.setOutputParameterName(0, "u_x");
        f_u.setOutputParameterName(1, "u_y");
        f_u.setOutputParameterName(2, "Strain");
        f_u.setOutputParameterName(3, "Stress");
//        MyFunctionStressStrain f_sigma;
//        f_sigma.setOutputParameterName(0, "strain_xx");

        SerialStaggeredMethod method(1e-5, 100);
        AsyncPartitionedSystem apart1;
        apart1.setAlgorithm(method);
        apart1.resizeOutputParameter(4);
        apart1.setOutputParameterName(0, "u_x");
        apart1.setOutputParameterName(1, "u_y");
        apart1.setOutputParameterName(2, "Strain");
        apart1.setOutputParameterName(3, "Stress");
        apart1.addProblem(f_u);
        apart1.connectParameters();

        TimeSteppingController timestepping;
        timestepping.setTransientSystem(apart1);

        //const double epsilon = 1.e-3;
        timestepping.setBeginning(.0);
        timestepping.solve(tim.getEnd());

//        const MyNodalFunctionScalar* r_f_ux = apart1.getOutput<MyNodalFunctionScalar>(apart1.getOutputParameterID("u_x"));
//        const MyNodalFunctionScalar* r_f_uy = apart1.getOutput<MyNodalFunctionScalar>(apart1.getOutputParameterID("u_y"));
        const MyIntegrationPointFunctionVector* r_f_strain = apart1.getOutput<MyIntegrationPointFunctionVector>(apart1.getOutputParameterID("Strain"));
        const MyIntegrationPointFunctionVector* r_f_stress = apart1.getOutput<MyIntegrationPointFunctionVector>(apart1.getOutputParameterID("Stress"));
//        const IDiscreteVector<double>* vec_r_f_ux = r_f_ux->getDiscreteData();
//        const IDiscreteVector<double>* vec_r_f_uy = r_f_uy->getDiscreteData();
        const MyIntegrationPointFunctionVector::MyDiscreteVector* vec_strain = r_f_strain->getDiscreteData();
        const MyIntegrationPointFunctionVector::MyDiscreteVector* vec_stress = r_f_stress->getDiscreteData();

//        r_f_ux->printout();
//        r_f_uy->printout();
//        r_f_strain->printout();
//        r_f_stress->printout();

        const MathLib::LocalVector &strain1 = (*vec_strain)[0][0];
        const MathLib::LocalVector &stress1 = (*vec_stress)[0][0];
        double E = solid.Youngs_modulus;
        double nu = solid.poisson_ratio;
        double sx = .0; //E/((1.+nu)*(1-2*nu))*((1-nu)*ex+nu*ey);
        double sy = -1e+6; //E/((1.-nu)*(1-2*nu))*(nu*ex+(1-nu)*ey);
        double sxy = .0; //0.5*E/(1.+nu)*exy;
        double ex = (1+nu)/E*((1-nu)*sx-nu*sy);
        double ey = (1+nu)/E*((1-nu)*sy-nu*sx);
        double exy = 2*(1+nu)/E*sxy;
        double sz = nu*E/((1.+nu)*(1-2*nu))*(ex+ey);
        double epsilon = 1e-6;
        ASSERT_NEAR(ex, strain1(0), epsilon);
        ASSERT_NEAR(ey, strain1(1), epsilon);
        ASSERT_NEAR(.0, strain1(2), epsilon);
        ASSERT_NEAR(exy, strain1(3), epsilon);
        epsilon = 1;
        ASSERT_NEAR(sx, stress1(0), epsilon);
        ASSERT_NEAR(sy, stress1(1), epsilon);
        ASSERT_NEAR(sz, stress1(2), epsilon);
        ASSERT_NEAR(sxy, stress1(3), epsilon);

    } catch (const char* e) {
        std::cout << "***Exception caught! " << e << std::endl;
    }

}
#endif

