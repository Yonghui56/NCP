/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SingleStepFDM.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "BaseLib/CodingTools.h"
#include "MathLib/LinAlg/LinearEquation/ILinearEquation.h"
#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/Serial/DiscreteSystem.h"
#include "DiscreteLib/Serial/DiscreteLinearEquation.h"
#include "DiscreteLib/Utils/SparsityBuilderFromNodeConnectivity.h"
#include "NumLib/Nonlinear/TemplateDiscreteNonlinearSolver.h"

#include "SolutionLib/Core/AbstractTimeSteppingAlgorithm.h"

#include "FdmFunction.h"
#include "BoundaryConditions.h"
#include "TemplateTransientLinearFDMFunction.h"

namespace FdmLib
{

/**
 * \brief Solution algorithm for linear transient problems using FEM with single time stepping method
 *
 * @tparam T_USER_FEM_PROBLEM      FEM problem class
 * @tparam T_LINEAR_SOLVER         Linear equation solver class
 */
template <
    class T_USER_PROBLEM,
    class T_LINEAR_SOLVER
    >
class SingleStepFDM
    : public SolutionLib::AbstractTimeSteppingAlgorithm
{
public:
    typedef T_USER_PROBLEM UserFemProblem;
    typedef TemplateTransientLinearFDMFunction<UserFemProblem, typename UserFemProblem::LinearAssemblerType> UserLinearFemFunction;
    typedef T_LINEAR_SOLVER LinearSolverType;

    /// constructor
    /// - initialize solution vectors
    /// - set up DoFs and equation index
    /// - prepare linear equation and solver
    /// - prepare linear functions
    SingleStepFDM(DiscreteLib::DiscreteSystem* dis, UserFemProblem* problem)
        : AbstractTimeSteppingAlgorithm(*problem->getTimeSteppingFunction()), 
          _problem(problem), _discrete_system(dis)
    {
        //this->caterogorizeGridPoint();

        const size_t n_var = problem->getNumberOfVariables();
        // create dof map
        for (size_t i=0; i<n_var; i++) {
            _dofManager.addVariableDoFs(dis->getMesh()->getID(), 0, problem->getField(i)->getNumberOfNodes());
        }
        _dofManager.construct();
        // initialize solution vectors
        _u_n1.resize(n_var, 0);
        _u_n1[0] = (FdmLib::FdmFunctionScalar*) problem->getIC(0)->clone();
        const size_t n_dofs = _dofManager.getTotalNumberOfActiveDoFs();
        _vec_n0 = dis->createVector<double>(n_dofs);
        _vec_n1 = dis->createVector<double>(n_dofs);
        _vec_n1_0 = dis->createVector<double>(n_dofs);
        _vec_st = dis->createVector<double>(n_dofs);
        FdmLib::FdmFunctionScalar *f_ic = (FdmLib::FdmFunctionScalar*) problem->getIC(0); //TODO one var
        *_vec_n1 = *f_ic->getNodalValues();
        // create linear equation systems
        _linear_solver = new LinearSolverType();
        _linear_eqs = _discrete_system->createLinearEquation<LinearSolverType, DiscreteLib::SparsityBuilderFromNodeConnectivity>(_linear_solver, &_dofManager);
        // setup functions
        _f_linear = new UserLinearFemFunction(problem, problem->getLinearAssembler(), _linear_eqs);
    };

    ///
    virtual ~SingleStepFDM()
    {
        _discrete_system->deleteLinearEquation(_linear_eqs);
        _discrete_system->deleteVector(_vec_n0);
        _discrete_system->deleteVector(_vec_n1);
        _discrete_system->deleteVector(_vec_n1_0);
        _discrete_system->deleteVector(_vec_st);
        BaseLib::releaseObject(_linear_solver);
        BaseLib::releaseObject(_f_linear);
    }

    /// solve 
    int solveTimeStep(const NumLib::TimeStep &t_n1)
    {
        // time step
        double dt = t_n1.getTime() - AbstractTimeSteppingAlgorithm::getTimeStepFunction()->getPrevious();
        NumLib::TimeStep this_t_n1;
        this_t_n1.assign(t_n1);
        this_t_n1.setTimeStepSize(dt);

        const size_t msh_id = _discrete_system->getMesh()->getID();

        // bc1
        std::vector<size_t> list_bc1_eqs_id;
        std::vector<double> list_bc1_val;
        for (size_t i=0; i<_problem->getNumberOfDirichletBC(); i++) {
            FdmLib::FdmDirichletBC<double> *bc1 = _problem->getFdmDirichletBC(i);
            bc1->setup();
            size_t varId = 0; //TODO var id
            std::vector<size_t> &list_bc_nodes = bc1->getListOfBCNodes();
            std::vector<double> &list_bc_values = bc1->getListOfBCValues();
            convertToEqsValues(_dofManager, varId, msh_id, list_bc_nodes, list_bc_values, list_bc1_eqs_id, list_bc1_val);
        }

        // st
        std::vector<size_t> list_st_eqs_id;
        std::vector<double> list_st_val;
        for (size_t i=0; i<_problem->getNumberOfNeumannBC(); i++) {
            FdmLib::FdmNeumannBC<double, double> *bc2 = _problem->getFdmNeumannBC(i);
            bc2->setup();
            std::vector<size_t> &list_bc_nodes = bc2->getListOfBCNodes();
            std::vector<double> &list_bc_values = bc2->getListOfBCValues();

            size_t varid = 0; //TODO var id
            const size_t msh_id = _problem->getMesh()->getID();
            convertToEqsValues(_dofManager, varid, msh_id, list_bc_nodes, list_bc_values, list_st_eqs_id, list_st_val);
        }
        (*_vec_st) = .0;
        for (size_t i=0; i<list_st_eqs_id.size(); i++) {
            (*_vec_st)[list_st_eqs_id[i]] = list_st_val[i];
        }

        // setup functions
        _f_linear->reset(&t_n1, _vec_n0);

        // initial guess
        *_vec_n1_0 = *_vec_n0;
        for (size_t i=0; i<list_bc1_eqs_id.size(); i++) {
            (*_vec_n1_0)[list_bc1_eqs_id[i]] = list_bc1_val[i];
        }
        
        // solve
        _f_linear->eval(*_vec_n1_0, *_vec_n1);

        _u_n1[0]->setNodalValues(*_vec_n1);

        return 0;
    }

    /// get the current solution
    FdmLib::FdmFunctionScalar* getCurrentSolution(int var_id)
    {
        return _u_n1[var_id];
    }

    ///
    UserFemProblem* getProblem() {return _problem;};

    /// get a linear solver
    LinearSolverType* getLinearEquationSolver() { return _linear_solver; };

    ///
    virtual void accept(const NumLib::TimeStep &t)
    {
        AbstractTimeSteppingAlgorithm::accept(t);
        *_vec_n0 = *_vec_n1; //copy current value to previous value
    };

private:
    DISALLOW_COPY_AND_ASSIGN(SingleStepFDM);

    void convertToEqsValues(const DiscreteLib::DofEquationIdTable &eqs_map, size_t var_id, size_t msh_id, const std::vector<size_t> &list_node_id, const std::vector<double> &list_node_values, std::vector<size_t> &list_eqs_id, std::vector<double> &list_eqs_val)
    {
        for (size_t j=0; j<list_node_id.size(); j++) {
            size_t pt_id = list_node_id[j];
            if (eqs_map.isActiveDoF(var_id, msh_id, pt_id)) {
                size_t eqs_id = eqs_map.mapEqsID(var_id, msh_id, pt_id);
                double bc_val = list_node_values[j];
                list_eqs_id.push_back(eqs_id);
                list_eqs_val.push_back(bc_val);
            }
        }
    }

private:
    DiscreteLib::DofEquationIdTable _dofManager;
    UserFemProblem* _problem;
    LinearSolverType* _linear_solver;
    DiscreteLib::DiscreteSystem *_discrete_system;
    DiscreteLib::IDiscreteLinearEquation* _linear_eqs;
    std::vector<FdmLib::FdmFunctionScalar*> _u_n1;
    UserLinearFemFunction* _f_linear;
    MyFemVector *_vec_n0;
    MyFemVector *_vec_n1_0;
    MyFemVector *_vec_n1;
    MyFemVector *_vec_st;
};


}
