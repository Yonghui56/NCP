/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TestFEM.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include <gtest/gtest.h>

#include "MathLib/Vector.h"
#include "MathLib/Interpolation/LinearInterpolation.h"
//#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"
#include "NumLib/Function/Function.h"
#include "MathLib/LinAlg/LinearEquation/DenseLinearEquation.h"

#include "GeoLib/Rectangle.h"

#include "MeshLib/Tools/MeshGenerator.h"

#include "DiscreteLib/Serial/DiscreteSystem.h"

#include "FemLib/Function/FemFunction.h"
#include "FemLib/BC/DirichletBC2FEM.h"
#include "FemLib/BC/NeumannBC2FEM.h"
#include "FemLib/Post/Extrapolation.h"

#include "TestExamples.h"
#include "TestUtil.h"

#include <vector>
#include <memory>

using namespace FemLib;
using namespace GeoLib;
using namespace MeshLib;
using namespace DiscreteLib;

typedef MathLib::Matrix<double> GlobalMatrixType;
typedef std::vector<double> GlobalVectorType;
typedef FemLib::FemNodalFunctionScalar<DiscreteLib::DiscreteSystem>::type MyNodalFunctionScalar;
typedef FemLib::FemNodalFunctionVector<DiscreteLib::DiscreteSystem>::type MyNodalFunctionVector;
typedef FemLib::FEMIntegrationPointFunctionVector<DiscreteLib::DiscreteSystem>::type MyIntegrationPointFunctionVector;

//static void outputLinearEQS(MathLib::Matrix<double> &globalA, std::vector<double> &globalRHS)
//{
//    std::cout << "A=" << std::endl;
//    globalA.write(std::cout);
//    std::cout << "x=" << std::endl;
//    for (size_t i=0; i<globalRHS.size(); i++)
//        std::cout << globalRHS[i] << " ";
//    std::cout << std::endl;
//}
//
//static void outputLinearEQS(MathLib::Matrix<double> &globalA, double* globalRHS)
//{
//    std::cout << "A=" << std::endl;
//    globalA.write(std::cout);
//    std::cout << "x=" << std::endl;
//    for (size_t i=0; i<globalA.getNRows(); i++)
//        std::cout << globalRHS[i] << " ";
//    std::cout << std::endl;
//}




static void getGWExpectedHead(std::vector<double> &expected)
{
    expected.resize(9);
    for (size_t i=0; i<9; i++) {
        if (i%3==0) expected[i] = 2.e+6;
        if (i%3==1) expected[i] = 1.e+6;
        if (i%3==2) expected[i] = 0.e+6;
    }
}

TEST(FEM, testUnstructuredMesh)
{
    GWFemTest gw;
    MeshLib::IMesh *msh = MeshGenerator::generateRegularQuadMesh(2.0, 2, .0, .0, .0);
    gw.define(msh);
    //#Solve
    GWFemTest::calculateHead(gw);

    DiscreteLib::IDiscreteVector<double>* h = gw.head->getDiscreteData();
    std::vector<double> expected;
    getGWExpectedHead(expected);

    ASSERT_DOUBLE_ARRAY_EQ(&expected[0], &(*h)[0], gw.head->getNumberOfNodes());
}

TEST(FEM, testStructuredMesh)
{
    GWFemTest gw;
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    gw.define(msh);
    //#Solve
    GWFemTest::calculateHead(gw);

    DiscreteLib::IDiscreteVector<double>* h = gw.head->getDiscreteData();
    std::vector<double> expected;
    getGWExpectedHead(expected);

    ASSERT_DOUBLE_ARRAY_EQ(&expected[0], &(*h)[0], gw.head->getNumberOfNodes());
}



TEST(FEM, ExtrapolateAverage1)
{
    GWFemTest gw;
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    gw.define(msh);
    GWFemTest::calculateHead(gw);
    GWFemTest::calculateVelocity(gw);

    MyNodalFunctionVector nodal_vel;
    nodal_vel.initialize(*gw.dis, PolynomialOrder::Linear);
    nodal_vel.setFeObjectContainer(gw._feObjects);
    FemExtrapolationAverage<DiscreteLib::DiscreteSystem, MathLib::LocalVector> extrapo;
    extrapo.extrapolate(*gw.vel, nodal_vel);

    DiscreteLib::IDiscreteVector<MathLib::LocalVector> *v = nodal_vel.getDiscreteData();
    MathLib::LocalVector expected(2);
    expected[0] = 1.e-5;
    expected[1] = .0;
    ASSERT_DOUBLE_ARRAY_EQ(expected, (*v)[0], expected.size());
}

class MyFunction : public NumLib::ITXFunction
{
public:
    virtual ~MyFunction() {};
    virtual void eval(const double* x, double &v) const
    {
        if (x[0]<1.0) v = 1e-11;
        else v = 2e-11;
    };

    MyFunction* clone() const {return 0;};

};

TEST(FEM, ExtrapolateAverage2)
{
    GWFemTest gw;
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    gw.define(msh);
    delete gw._K;
    gw._K = new MyFunction();
    NumLib::TXFunctionConstant f_bc2(2.e+6);
    std::vector<size_t> bc2_nodes;
    std::vector<double> bc2_vals;
    FemLib::DirichletBC2FEM bc2(*msh, gw.rec->getLeft(), f_bc2, bc2_nodes, bc2_vals);
    gw.vec_bc1_nodes.insert(gw.vec_bc1_nodes.end(), bc2_nodes.begin(), bc2_nodes.end());
    gw.vec_bc1_vals.insert(gw.vec_bc1_vals.end(), bc2_vals.begin(), bc2_vals.end());
    gw.vec_bc2_nodes.clear();
    gw.vec_bc2_vals.clear();

    GWFemTest::calculateHead(gw);
    GWFemTest::calculateVelocity(gw);

    MyNodalFunctionVector nodal_vel;
    nodal_vel.initialize(*gw.dis, PolynomialOrder::Linear);
    nodal_vel.setFeObjectContainer(gw._feObjects);
    FemExtrapolationAverage<DiscreteLib::DiscreteSystem, MathLib::LocalVector> extrapo;
    extrapo.extrapolate(*gw.vel, nodal_vel);

    std::vector<double> exH(9);
    for (size_t i=0; i<9; i++) {
        if (i%3==0) exH[i] = 2.e+6;
        if (i%3==1) exH[i] = 2./3.*1.e+6;
        if (i%3==2) exH[i] = 0.e+6;
    }

    DiscreteLib::IDiscreteVector<double> *h = gw.head->getDiscreteData();
    ASSERT_DOUBLE_ARRAY_EQ(&exH[0], &(*h)[0], gw.head->getNumberOfNodes(), 10);

    DiscreteLib::IDiscreteVector<MathLib::LocalVector> *v = nodal_vel.getDiscreteData();
    MathLib::LocalVector expected(2);
    expected[0] = 4./3.*1.e-5;
    expected[1] = .0;
    ASSERT_DOUBLE_ARRAY_EQ(expected, (*v)[0], expected.size());
}

TEST(FEM, LIE_LINE_IN_2D)
{
    //--------------------------------------------------------------------------
    // make 2d quad + 1d line
    MeshLib::UnstructuredMesh *msh = MeshGenerator::generateRegularQuadMesh(1.0, 1, .0, .0, .0);
    assert (*msh->getNodeCoordinatesRef(1) == GeoLib::Point(1.0, .0, .0));
    assert (*msh->getNodeCoordinatesRef(3) == GeoLib::Point(1.0, 1.0, .0));
    // horizontal line
    MeshLib::Line* e_line1 = new MeshLib::Line(0);
    e_line1->setNodeID(0, 0);
    e_line1->setNodeID(1, 1);
    msh->addElement(e_line1);
    // vertical line
    MeshLib::Line* e_line2 = new MeshLib::Line(1);
    e_line2->setNodeID(0, 1);
    e_line2->setNodeID(1, 3);
    msh->addElement(e_line2);
    // diagonal line
    MeshLib::Line* e_line3 = new MeshLib::Line(2);
    e_line3->setNodeID(0, 0);
    e_line3->setNodeID(1, 3);
    msh->addElement(e_line3);

    FemLib::FemNaturalCoordinates local_fe_map(new FemLib::FemShapeLine2());

    //--------------------------------------------------------------------------
    // check horizontal Line
    FemLib::LINE2 fe_line(msh);
    {
        fe_line.configure(*e_line1);
        double natural_x[2] = {.0, .0};
        fe_line.computeBasisFunctions(natural_x);
        MathLib::LocalMatrix &N = *fe_line.getBasisFunction();
        MathLib::LocalMatrix &dN = *fe_line.getGradBasisFunction();
        double det_j = fe_line.getDetJ();

//        std::cout << "N=\n" << N << std::endl;
//        std::cout << "dN=\n" << dN << std::endl;
//        std::cout << "det J=" << det_j << std::endl;
//        MeshLib::ElementCoordinatesMappingLocal* e_map = (MeshLib::ElementCoordinatesMappingLocal*)e_line1->getMappedCoordinates();
//        std::cout << "R=\n" << e_map->getRotationMatrixToOriginal() << std::endl;
//        std::cout << "R^T=\n" << e_map->getRotationMatrixToLocal() << std::endl;
//        std::cout << "x'0="; e_map->getNodePoint(0)->write(std::cout); std::cout << std::endl;
//        std::cout << "x'1="; e_map->getNodePoint(1)->write(std::cout); std::cout << std::endl;

        ASSERT_EQ(1, N.rows());
        ASSERT_EQ(2, N.cols());
        ASSERT_EQ(2, dN.rows());
        ASSERT_EQ(2, dN.cols());

        local_fe_map.initialize(*e_line1);
        //const FemLib::CoordinateMappingProperty* prop = local_fe_map.compute(natural_x);

        ASSERT_EQ(0.5, N(0,0));
        ASSERT_EQ(0.5, N(0,1));
        ASSERT_EQ(-1.0, dN(0,0));
        ASSERT_EQ(1.0, dN(0,1));
        ASSERT_EQ(0.0, dN(1,0));
        ASSERT_EQ(0.0, dN(1,1));
        ASSERT_EQ(0.5, det_j);
    }
    //--------------------------------------------------------------------------
    // check vertical Line
    {
        fe_line.configure(*e_line2);
        double natural_x[2] = {.0, .0};
        fe_line.computeBasisFunctions(natural_x);
        MathLib::LocalMatrix &N = *fe_line.getBasisFunction();
        MathLib::LocalMatrix &dN = *fe_line.getGradBasisFunction();
        double det_j = fe_line.getDetJ();

//        std::cout << "N=\n" << N << std::endl;
//        std::cout << "dN=\n" << dN << std::endl;
//        std::cout << "det J=" << det_j << std::endl;
//        MeshLib::ElementCoordinatesMappingLocal* e_map = (MeshLib::ElementCoordinatesMappingLocal*)e_line2->getMappedCoordinates();
//        std::cout << "R=\n" << e_map->getRotationMatrixToOriginal() << std::endl;
//        std::cout << "R^T=\n" << e_map->getRotationMatrixToLocal() << std::endl;
//        std::cout << "x'0="; e_map->getNodePoint(0)->write(std::cout); std::cout << std::endl;
//        std::cout << "x'1="; e_map->getNodePoint(1)->write(std::cout); std::cout << std::endl;

        ASSERT_EQ(0.5, N(0,0));
        ASSERT_EQ(0.5, N(0,1));
        ASSERT_EQ(0., dN(0,0));
        ASSERT_EQ(0., dN(0,1));
        ASSERT_EQ(-1., dN(1,0));
        ASSERT_EQ(1., dN(1,1));
//        ASSERT_EQ(1., dN(1,0));
//        ASSERT_EQ(-1., dN(1,1));
        ASSERT_EQ(.5, det_j);
    }
    //--------------------------------------------------------------------------
    // check diagonal Line
    fe_line.configure(*e_line3);
    {
        double natural_x[2] = {.0, .0};
        fe_line.computeBasisFunctions(natural_x);
        MathLib::LocalMatrix &N = *fe_line.getBasisFunction();
        MathLib::LocalMatrix &dN = *fe_line.getGradBasisFunction();
        double det_j = fe_line.getDetJ();

//        std::cout << "N=\n" << N << std::endl;
//        std::cout << "dN=\n" << dN << std::endl;
//        std::cout << "det J=" << det_j << std::endl;
//        MeshLib::ElementCoordinatesMappingLocal* e_map = (MeshLib::ElementCoordinatesMappingLocal*)e_line3->getMappedCoordinates();
//        std::cout << "R=\n" << e_map->getRotationMatrixToOriginal() << std::endl;
//        std::cout << "R^T=\n" << e_map->getRotationMatrixToLocal() << std::endl;
//        std::cout << "x'0="; e_map->getNodePoint(0)->write(std::cout); std::cout << std::endl;
//        std::cout << "x'1="; e_map->getNodePoint(1)->write(std::cout); std::cout << std::endl;

        double epsilon = 1e-3;
        ASSERT_NEAR(0.5, N(0,0), epsilon);
        ASSERT_NEAR(0.5, N(0,1), epsilon);
        ASSERT_NEAR(-0.5, dN(0,0), epsilon);
        ASSERT_NEAR(0.5, dN(0,1), epsilon);
        ASSERT_NEAR(-0.5, dN(1,0), epsilon);
        ASSERT_NEAR(0.5, dN(1,1), epsilon);
//        ASSERT_NEAR(0.5, dN(1,0), epsilon);
//        ASSERT_NEAR(-0.5, dN(1,1), epsilon);
        ASSERT_NEAR(.7071, det_j, epsilon);
    }

}
