/* File : example.c */

#include "example.h"

#include "MathLib/Coupling/Algorithm/IConvergenceCheck.h"
#include "MathLib/Coupling/Algorithm/SerialStaggeredMethod.h"
#include "MathLib/LinAlg/LinearEquations/LisInterface.h"
#include "GeoLib/Shape/Line.h"
#include "MeshLib/Tools/MeshGenerator.h"
#include "NumLib/TransientCoupling/AsyncPartitionedSystem.h"
#include "Tests/TestExamples.h"
#include "Tests/Geo/Model/Head.h"

#define M_PI 3.14159265358979323846

/* Move the shape to a new location */
void Shape::move(double dx, double dy) {
  x += dx;
  y += dy;
}

int Shape::nshapes = 0;

double Circle::area(void) {
  return M_PI*radius*radius;
}

double Circle::perimeter(void) {
  return 2*M_PI*radius;
}

double Square::area(void) {
  return width*width;
}

double Square::perimeter(void) {
  return 4*width;
}

Geo::GWFemProblem* defineGWProblem1D(DiscreteSystem &dis, GeoLib::Line &line, Geo::PorousMedia &pm)
{
    LagrangeFeObjectContainer* _feObjects = new LagrangeFeObjectContainer(*dis.getMesh());
    //equations
    Geo::GWFemProblem::LinearAssemblerType* linear_assembler = new Geo::GWFemProblem::LinearAssemblerType(*_feObjects, pm);
    Geo::GWFemProblem::ResidualAssemblerType* r_assembler = new Geo::GWFemProblem::ResidualAssemblerType(*_feObjects, pm);
    Geo::GWFemProblem::JacobianAssemblerType* j_eqs = new Geo::GWFemProblem::JacobianAssemblerType(*_feObjects, pm);
    //IVBV problem
    Geo::GWFemProblem* _problem = new Geo::GWFemProblem(dis, *dis.getMesh(), linear_assembler, r_assembler, j_eqs);
    //BC
    size_t headId = _problem->createField(PolynomialOrder::Linear);
    FemNodalFunctionScalar* _head = _problem->getField(headId);
    _problem->setIC(headId, *_head);
    MathLib::SpatialFunctionConstant<double> f1(.0);
    _problem->addDirichletBC(headId, *line.getPoint2(), false, f1);
    MathLib::SpatialFunctionConstant<double> f2(-1e-5);
    _problem->addNeumannBC(headId, *line.getPoint1(), false, f2);

    return _problem;
}

class FemFunctionConvergenceCheck : public MathLib::IConvergenceCheck
{
public:
    bool isConverged(UnnamedParameterSet& vars_prev, UnnamedParameterSet& vars_current, double eps, double &v_diff)
    {
        for (size_t i=0; i<vars_prev.size(); i++) {
#if 1
            if (vars_prev.getName(i).compare("h")==0 || vars_prev.getName(i).compare("c")==0) {
                const FemNodalFunctionScalar* f_fem_prev = vars_prev.get<FemNodalFunctionScalar>(i);
                const FemNodalFunctionScalar* f_fem_cur = vars_current.get<FemNodalFunctionScalar>(i);
                v_diff = f_fem_cur->norm_diff(*f_fem_prev);
            } else if (vars_prev.getName(i).compare("v")==0) {
                const FEMIntegrationPointFunctionVector2d* f_fem_prev = vars_prev.get<FEMIntegrationPointFunctionVector2d>(i);
                const FEMIntegrationPointFunctionVector2d* f_fem_cur = vars_current.get<FEMIntegrationPointFunctionVector2d>(i);
                v_diff = f_fem_cur->norm_diff(*f_fem_prev);
            }
#endif
            if (v_diff>eps) {
                return false;
            }
        }
        return true;
    }
};

void gw1(double len, size_t div, std::vector<double> *x, std::vector<double> *results)
{
    MeshLib::IMesh *msh = MeshGenerator::generateLineMesh(len, div, .0, .0, .0);
#if 1
    GeoLib::Line* line = new GeoLib::Line(Point(0.0, 0.0, 0.0),  Point(len, 0.0, 0.0));
    Geo::PorousMedia pm;
    pm.hydraulic_conductivity = new MathLib::SpatialFunctionConstant<double>(1.e-11);
    pm.porosity = new MathLib::SpatialFunctionConstant<double>(0.2);
    DiscreteSystem dis(*msh);
    Geo::GWFemProblem* pGW = defineGWProblem1D(dis, *line, pm);
    TimeStepFunctionConstant tim(.0, 1e+4, 1e+3);
    pGW->setTimeSteppingFunction(tim);

    Base::Options options;
    Base::Options* op_lis = options.addSubGroup("LinearSolver");
    op_lis->addOption("solver_type", "CG");
    op_lis->addOption("precon_type", "NONE");
    op_lis->addOptionAsNum("error_tolerance", 1e-10);
    op_lis->addOptionAsNum("max_iteration_step", 500);

    Geo::FunctionHead<MathLib::DenseLinearEquation> f_head;
//    Geo::FunctionHead<MathLib::CRSLisSolver> f_head;
    f_head.define(&dis, pGW, options);
    f_head.setOutputParameterName(0, "h");

    FemFunctionConvergenceCheck checker;
    MathLib::SerialStaggeredMethod method(checker, 1e-5, 100);
    NumLib::AsyncPartitionedSystem apart1;
    apart1.setAlgorithm(method);
    apart1.resizeOutputParameter(1);
    apart1.setOutputParameterName(0, "h");
    apart1.addProblem(f_head);
    apart1.connectParameters();

    TimeSteppingController timestepping;
    timestepping.addTransientSystem(apart1);

    //const double epsilon = 1.e-3;
    timestepping.setBeginning(.0);
    timestepping.solve(1.0);

    const FemLib::FemNodalFunctionScalar* r_f_head = apart1.getOutput<FemLib::FemNodalFunctionScalar>(apart1.getOutputParameterID("h"));
    const DiscreteLib::DiscreteVector<double>* h = r_f_head->getNodalValues();
#endif
    x->resize(msh->getNumberOfNodes());
    for (size_t i=0; i<x->size(); i++)
        (*x)[i] = msh->getNodeCoordinatesRef(i)->getData()[0];

#if 1
    results->resize(h->size());
    for (size_t i=0; i<results->size(); i++)
        (*results)[i] = (*h)[i];
#else
    results->resize(x->size(), 1.0);
#endif
}

void range(std::vector<int> *rangevec)
{
    int i;
    for (i=0; i< rangevec->size(); i++)
        (*rangevec)[i] = i;
}
