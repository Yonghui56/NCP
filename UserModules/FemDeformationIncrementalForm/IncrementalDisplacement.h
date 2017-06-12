/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IncrementalDisplacement.h
 *
 * Created on 2012-11-29 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <string>

#include "BaseLib/CodingTools.h"
#include "FemLib/Function/FemIntegrationPointFunction.h"
#include "NumLib/Function/TXVectorFunctionAsColumnData.h"
#include "SolutionLib/Fem/TemplateFemEquation.h"
#include "SolutionLib/Fem/FemIVBVProblem.h"
#include "SolutionLib/Fem/SingleStepFEM.h"

#include "ProcessLib/AbstractTransientProcess.h"

#include "FemElastoPlasticLinearLocalAssembler.h"
#include "FemElastoPlasticResidualLocalAssembler.h"
#include "FemElastoPlasticJacobianLocalAssembler.h"

/**
 *
 */
template<class T_DISCRETE_SYSTEM, class T_LINEAR_SOLVER>
class FunctionIncrementalDisplacement:
        public ProcessLib::AbstractTransientProcess
{
public:
    enum Out
    {
        Displacement = 0, Strain = 1, Stress = 2
    };

    typedef T_DISCRETE_SYSTEM MyDiscreteSystem;
    typedef T_LINEAR_SOLVER MyLinearSolver;
    // local assembler
    typedef FemElastoPlasticLinearLocalAssembler MyLinearAssemblerType;
    typedef FemElastoPlasticResidualLocalAssembler MyResidualAssemblerType;
    typedef FemElastoPlasticJacobianLocalAssembler MyJacobianAssemblerType;
    // Equation definition
    typedef SolutionLib::TemplateFemEquation<MyDiscreteSystem, MyLinearSolver,
            MyLinearAssemblerType, MyResidualAssemblerType,
            MyJacobianAssemblerType> MyEquationType;
    // IVBV problem definition
    typedef SolutionLib::FemIVBVProblem<MyDiscreteSystem, MyEquationType> MyProblemType;
    // Solution algorithm definition
    typedef SolutionLib::SingleStepFEM<MyProblemType, MyLinearSolver> MySolutionType;

    typedef typename FemLib::FemNodalFunctionScalar<MyDiscreteSystem>::type MyNodalFunctionScalar;
    typedef typename FemLib::FemNodalFunctionVector<MyDiscreteSystem>::type MyNodalFunctionVector;
    typedef typename NumLib::TXVectorFunctionAsColumnData<MyNodalFunctionVector> NodalPointScalarWrapper;
    typedef typename FemLib::FEMIntegrationPointFunctionVector<MyDiscreteSystem>::type MyIntegrationPointFunctionVector;
    typedef typename NumLib::TXVectorFunctionAsColumnData<
            MyIntegrationPointFunctionVector> IntegrationPointScalarWrapper;
    typedef typename MyProblemType::MyVariable MyVariable;

    ///
    FunctionIncrementalDisplacement()
            : ProcessLib::AbstractTransientProcess("INCREMENTAL_DEFORMATION", 0,
                    3), _problem(nullptr), _solution(nullptr), _feObjects(
                    nullptr), _delta_displacement(nullptr), _previous_displacement(
                    nullptr), _current_displacement(nullptr), _previous_strain(
                    nullptr), _current_strain(nullptr), _previous_stress(
                    nullptr), _current_stress(nullptr)
    {
        // set default parameter name
        this->setOutputParameterName(Displacement, "Displacement");
        this->setOutputParameterName(Strain, "Strain");
        this->setOutputParameterName(Stress, "Stress");
    }
    ;

    ///
    virtual ~FunctionIncrementalDisplacement()
    {
        BaseLib::releaseObject(_problem, _solution, _feObjects);
        BaseLib::releaseObject(_delta_displacement, _previous_displacement,
                _current_displacement);
        BaseLib::releaseObject(_previous_strain, _current_strain,
                _previous_stress, _current_stress);
    }

    /// initialize this process
    virtual bool initialize(const BaseLib::Options &option);

    /// finalize but nothing to do here
    virtual void finalize()
    {
    }

    ///
    virtual NumLib::IConvergenceCheck* getConvergenceChecker()
    {
        return &_checker;
    }

protected:
    virtual void initializeTimeStep(const NumLib::TimeStep &time);
    virtual void postSolutionAlgorithm(const NumLib::TimeStep &time);
    virtual void postTimeStep(const NumLib::TimeStep &time);
    virtual void updateOutputParameter(const NumLib::TimeStep &time);

    virtual MySolutionType* getSolution()
    {
        return _solution;
    }

    virtual void output(const NumLib::TimeStep &time);

private:

    DISALLOW_COPY_AND_ASSIGN(FunctionIncrementalDisplacement);

private:
    MyProblemType* _problem;
    MySolutionType* _solution;
    FemLib::LagrangeFeObjectContainer* _feObjects;
    NumLib::DiscreteDataConvergenceCheck _checker;
    MyNodalFunctionVector* _delta_displacement;
    MyNodalFunctionVector* _previous_displacement;
    MyNodalFunctionVector* _current_displacement;
    MyIntegrationPointFunctionVector* _previous_strain;
    MyIntegrationPointFunctionVector* _current_strain;
    MyIntegrationPointFunctionVector* _previous_stress;
    MyIntegrationPointFunctionVector* _current_stress;
    std::vector<NodalPointScalarWrapper*> _vec_u_components;
    std::vector<IntegrationPointScalarWrapper*> _vec_strain_components;
    std::vector<IntegrationPointScalarWrapper*> _vec_stress_components;

};

#include "IncrementalDisplacement.tpp"

