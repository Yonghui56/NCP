

#include "logog.hpp"

#include "MeshLib/Tools/MeshGenerator.h"
#include "FemLib/Function/FemNodalFunction.h"
#include "NumLib/Function/TXFunctionBuilder.h"
#include "OutputIO/OutputBuilder.h"
#include "OutputIO/OutputTimingBuilder.h"
#include "SolutionLib/Fem/FemSourceTerm.h"
#include "SolutionLib/Fem/FemVariable.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Solid.h"
#include "PhysicsLib/FemLinearElasticTools.h"
#include "Ogs6FemData.h"
#include "FemVariableBuilder.h"

template <class T1, class T2>
typename FunctionDisplacementPressure<T1,T2>::MyVariable* FunctionDisplacementPressure<T1,T2>::getDisplacementComponent(MyVariable *u_x, MyVariable* u_y, MyVariable* u_z, const std::string &var_name)
{
    if (var_name.find("_X")!=std::string::npos) {
        return u_x;
    } else if (var_name.find("_Y")!=std::string::npos) {
        return u_y;
    } else {
        return u_z;
    }
}

template <class T1, class T2>
size_t FunctionDisplacementPressure<T1,T2>::getDisplacementComponentIndex(const std::string &var_name) const
{
    if (var_name.find("_X")!=std::string::npos) {
        return 0;
    } else if (var_name.find("_Y")!=std::string::npos) {
        return 1;
    } else {
        return 2;
    }
}

template <class T1, class T2>
bool FunctionDisplacementPressure<T1,T2>::initialize(const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    size_t msh_id = option.getOptionAsNum<size_t>("MeshID");
    size_t time_id = option.getOptionAsNum<size_t>("TimeGroupID");

    //--------------------------------------------------------------------------
    // set up mesh and FE objects
    //--------------------------------------------------------------------------
    MeshLib::IMesh* msh = femData->list_mesh[msh_id];
    INFO("->generating higher order mesh...");
    MeshLib::MeshGenerator::generateHigherOrderUnstrucuredMesh(*(MeshLib::UnstructuredMesh*)msh, 2);
    INFO("* mesh id %d: order=%d, nodes=%d, elements=%d", msh_id, msh->getMaxiumOrder(), msh->getNumberOfNodes(msh->getMaxiumOrder()), msh->getNumberOfElements());
    MyDiscreteSystem* dis = 0;
    dis = DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(msh);
    _feObjects = new FemLib::LagrangeFeObjectContainer(msh);

    //--------------------------------------------------------------------------
    // set up problem
    //--------------------------------------------------------------------------
    _problem = new MyProblemType(dis);
    NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];
    _problem->setTimeSteppingFunction(*tim);

    //--------------------------------------------------------------------------
    // set up variables
    //--------------------------------------------------------------------------
    // definitions
    MyVariable* u_x = _problem->addVariable("u_x", FemLib::PolynomialOrder::Quadratic);
    MyVariable* u_y = _problem->addVariable("u_y", FemLib::PolynomialOrder::Quadratic);
    MyVariable* p = _problem->addVariable("p", FemLib::PolynomialOrder::Linear);
    _var_p_id = p->getID();
    FemVariableBuilder var_builder;
    var_builder.doit(this->getOutputParameterName(Displacement)+"_X1", option, msh, femData->geo, femData->geo_unique_name, _feObjects, u_x);
    var_builder.doit(this->getOutputParameterName(Displacement)+"_Y1", option, msh, femData->geo, femData->geo_unique_name, _feObjects, u_y);
    var_builder.doit(this->getOutputParameterName(Pressure), option, msh, femData->geo, femData->geo_unique_name, _feObjects, p);

    //--------------------------------------------------------------------------
    // set up equations
    //--------------------------------------------------------------------------
    MyEquationType* eqs = _problem->createEquation();
    std::vector<size_t> vec_orders;
    for (size_t i=0; i<_problem->getNumberOfVariables(); i++)
        vec_orders.push_back(_problem->getVariable(i)->getCurrentOrder());

    MyLinearAssemblerType* linear_assembler = new MyLinearAssemblerType(_feObjects, _problem->getNumberOfVariables(), vec_orders);
    MyResidualAssemblerType* r_assembler = new MyResidualAssemblerType(_feObjects, _problem->getNumberOfVariables(), vec_orders);
    MyJacobianAssemblerType* j_eqs = new MyJacobianAssemblerType(_feObjects, _problem->getNumberOfVariables(), vec_orders);
    eqs->initialize(linear_assembler, r_assembler, j_eqs);


    //--------------------------------------------------------------------------
    // set up this solution algorithm
    //--------------------------------------------------------------------------
    _solution = new MySolutionType(dis, _problem);
    typename MySolutionType::LinearSolverType* linear_solver = _solution->getLinearEquationSolver();
    const BaseLib::Options* optNum = option.getSubGroup("Numerics");
    linear_solver->setOption(*optNum);
    _solution->getNonlinearSolver()->setOption(*optNum);
    _solution->getDofEquationIdTable()->setNumberingType(DiscreteLib::DofNumberingType::BY_POINT);

    //--------------------------------------------------------------------------
    // set up data for output
    //--------------------------------------------------------------------------
    // create u variable which is vector
    MathLib::LocalVector tmp_u0(3);
    tmp_u0 *= .0;
    _displacement = new MyNodalFunctionVector();
    _displacement->initialize(*dis, u_x->getCurrentOrder(), tmp_u0);
    for (size_t i=0; i<_displacement->getNumberOfNodes(); i++) {
        _displacement->getValue(i)(0) = _solution->getCurrentSolution(0)->getValue(i);
        _displacement->getValue(i)(1) = _solution->getCurrentSolution(1)->getValue(i);
    }
    for (size_t i=0; i<2; i++) {
        _vec_u_components.push_back(new NodalPointScalarWrapper(_displacement, i));
    }

    // set initial output
    OutputVariableInfo var(this->getOutputParameterName(Displacement), msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 2, _displacement);
    femData->outController.setOutput(var.name, var);
    for (size_t i=0; i<_vec_u_components.size(); i++) {
        OutputVariableInfo var1(this->getOutputParameterName(Displacement) + getDisplacementComponentPostfix(i), msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _vec_u_components[i]);
        femData->outController.setOutput(var1.name, var1);
    }
    OutputVariableInfo outP(this->getOutputParameterName(Pressure), msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(_var_p_id));
    femData->outController.setOutput(outP.name, outP);

    //--------------------------------------------------------------------------
    // set initial output
    //--------------------------------------------------------------------------
    // initial output parameter
    this->setOutput(Displacement, _displacement);
    this->setOutput(Pressure, _solution->getCurrentSolution(_var_p_id));


    return true;
}

template <class T1, class T2>
void FunctionDisplacementPressure<T1,T2>::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{
    for (size_t i=0; i<_displacement->getNumberOfNodes(); i++) {
        _displacement->getValue(i)(0) = _solution->getCurrentSolution(0)->getValue(i);
        _displacement->getValue(i)(1) = _solution->getCurrentSolution(1)->getValue(i);
    }
    setOutput(Displacement, _displacement);
    this->setOutput(Pressure, _solution->getCurrentSolution(_var_p_id));

    //calculateStressStrain();
}

template <class T1, class T2>
void FunctionDisplacementPressure<T1,T2>::output(const NumLib::TimeStep &/*time*/)
{
    //update data for output
    const size_t msh_id = _problem->getDiscreteSystem()->getMesh()->getID();
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    OutputVariableInfo var(this->getOutputParameterName(Displacement), msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 2, _displacement);
    femData->outController.setOutput(var.name, var);
    for (size_t i=0; i<_vec_u_components.size(); i++) {
        OutputVariableInfo var1(this->getOutputParameterName(Displacement) + getDisplacementComponentPostfix(i), msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _vec_u_components[i]);
        femData->outController.setOutput(var1.name, var1);
    }
    OutputVariableInfo outP(this->getOutputParameterName(Pressure), msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(_var_p_id));
    femData->outController.setOutput(outP.name, outP);
};

