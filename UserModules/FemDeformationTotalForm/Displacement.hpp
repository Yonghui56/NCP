
//#include "Displacement.h"

#include "logog.hpp"

#include "FemLib/Function/FemNodalFunction.h"
#include "NumLib/Function/TXFunctionBuilder.h"
#include "OutputIO/OutputBuilder.h"
#include "OutputIO/OutputTimingBuilder.h"
#include "SolutionLib/Fem/FemSourceTerm.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Solid.h"
#include "PhysicsLib/DeformationTools.h"
#include "Ogs6FemData.h"
#include "FemVariableBuilder.h"

template <class T1, class T2>
typename FunctionDisplacement<T1,T2>::MyVariable* FunctionDisplacement<T1,T2>::getDisplacementComponent(MyVariable *u_x, MyVariable* u_y, MyVariable* u_z, const std::string &var_name)
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
bool FunctionDisplacement<T1,T2>::initialize(const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    size_t msh_id = option.getOptionAsNum<size_t>("MeshID");
    size_t time_id = option.getOptionAsNum<size_t>("TimeGroupID");
    NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];

    //mesh and FE objects
    MeshLib::IMesh* msh = femData->list_mesh[msh_id];
    MyDiscreteSystem* dis = 0;
    dis = DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(msh);
    _feObjects = new FemLib::LagrangeFeObjectContainer(msh);

    // equations
    MyLinearAssemblerType* linear_assembler = new MyLinearAssemblerType(_feObjects);
    MyResidualAssemblerType* r_assembler = new MyResidualAssemblerType(_feObjects);
    MyJacobianAssemblerType* j_eqs = new MyJacobianAssemblerType(_feObjects);

    // set up problem
    _problem = new MyProblemType(dis);
    MyEquationType* eqs = _problem->createEquation();
    eqs->initialize(linear_assembler, r_assembler, j_eqs);
    _problem->setTimeSteppingFunction(*tim);
    // set up variable
    MyVariable* u_x = _problem->addVariable("u_x");
    MyVariable* u_y = _problem->addVariable("u_y");
    FemVariableBuilder var_builder;
    var_builder.doit(this->getOutputParameterName(Displacement)+"_X", option, msh, femData->geo, femData->geo_unique_name, _feObjects, u_x);
    var_builder.doit(this->getOutputParameterName(Displacement)+"_Y", option, msh, femData->geo, femData->geo_unique_name, _feObjects, u_y);

    // set up solution
    _solution = new MySolutionType(dis, _problem);
    typename MySolutionType::LinearSolverType* linear_solver = _solution->getLinearEquationSolver();
    const BaseLib::Options* optNum = option.getSubGroup("Numerics");
    linear_solver->setOption(*optNum);
    _solution->getNonlinearSolver()->setOption(*optNum);
    DiscreteLib::DofEquationIdTable* dofMapping = _solution->getDofEquationIdTable();
    dofMapping->setNumberingType(DiscreteLib::DofNumberingType::BY_POINT);
    dofMapping->setLocalNumberingType(DiscreteLib::DofNumberingType::BY_POINT);

    // create u variable which is vector
    MathLib::LocalVector tmp_u0(3);
    tmp_u0 *= .0;
    _displacement = new MyNodalFunctionVector();
    _displacement->initialize(*dis, FemLib::PolynomialOrder::Linear, tmp_u0);
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

    // initial output parameter
    this->setOutput(Displacement, _displacement);


    return true;
}

template <class T1, class T2>
void FunctionDisplacement<T1,T2>::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{
    for (size_t i=0; i<_displacement->getNumberOfNodes(); i++) {
        _displacement->getValue(i)(0) = _solution->getCurrentSolution(0)->getValue(i);
        _displacement->getValue(i)(1) = _solution->getCurrentSolution(1)->getValue(i);
    }
    setOutput(Displacement, _displacement);

    //calculateStressStrain();
}

template <class T1, class T2>
void FunctionDisplacement<T1,T2>::output(const NumLib::TimeStep &/*time*/)
{
    const size_t msh_id = _problem->getDiscreteSystem()->getMesh()->getID();
    //update data for output
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    OutputVariableInfo var(this->getOutputParameterName(Displacement),  msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 2, _displacement);
    femData->outController.setOutput(var.name, var);
    for (size_t i=0; i<_vec_u_components.size(); i++) {
        OutputVariableInfo var1(this->getOutputParameterName(Displacement) + getDisplacementComponentPostfix(i), msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _vec_u_components[i]);
        femData->outController.setOutput(var1.name, var1);
    }
};

