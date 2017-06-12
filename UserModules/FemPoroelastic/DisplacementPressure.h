/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DisplacementPressure.h
 *
 * Created on 2012-08-21 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include "BaseLib/CodingTools.h"
#include "FemLib/Function/FemIntegrationPointFunction.h"
#include "NumLib/Function/TXVectorFunctionAsColumnData.h"
#include "SolutionLib/Fem/TemplateFemEquation.h"
#include "SolutionLib/Fem/FemIVBVProblem.h"
#include "SolutionLib/Fem/SingleStepFEM.h"

#include "ProcessLib/AbstractTransientProcess.h"

#include "FemPoroelasticLinearEQSLocalAssembler.h"
#include "FemPoroelasticResidualLocalAssembler.h"
#include "FemPoroelasticJacobianLocalAssembler.h"

/**
 * \brief Monolithic poro-elasticity calculator using FEM
 */
template <class T_DISCRETE_SYSTEM, class T_LINEAR_SOLVER>
class FunctionDisplacementPressure
: public ProcessLib::AbstractTransientProcess
{
public:
    enum Out { Displacement=0, Pressure=1 };

    typedef T_DISCRETE_SYSTEM MyDiscreteSystem;
    typedef T_LINEAR_SOLVER MyLinearSolver;
    // local assembler
    typedef FemPoroelasticLinearEQSLocalAssembler MyLinearAssemblerType;
    typedef FemPoroelasticResidualLocalAssembler MyResidualAssemblerType;
    typedef FemPoroelasticJacobianLocalAssembler MyJacobianAssemblerType;
    // Equation definition
    typedef SolutionLib::TemplateFemEquation<
            MyDiscreteSystem,
            MyLinearSolver,
            MyLinearAssemblerType,
            MyResidualAssemblerType,
            MyJacobianAssemblerType
            >
    MyEquationType;
    // IVBV problem definition
    typedef SolutionLib::FemIVBVProblem<
            MyDiscreteSystem,
            MyEquationType
            > MyProblemType;
    // Solution algorithm definition
    typedef SolutionLib::SingleStepFEM
            <
                MyProblemType,
                MyLinearSolver
            > MySolutionType;

    typedef typename FemLib::FemNodalFunctionScalar<MyDiscreteSystem>::type MyNodalFunctionScalar;
    typedef typename FemLib::FemNodalFunctionVector<MyDiscreteSystem>::type MyNodalFunctionVector;
    typedef typename NumLib::TXVectorFunctionAsColumnData<MyNodalFunctionVector> NodalPointScalarWrapper;
    typedef typename MyProblemType::MyVariable MyVariable;

    ///
    FunctionDisplacementPressure()
    : ProcessLib::AbstractTransientProcess("DEFORMATION_FLOW", 0, 2),
      _problem(0), _solution(0), _feObjects(0), _displacement(0), _var_p_id(0)//, _strain(0), _stress(0)
    {
        // set default parameter name
        this->setOutputParameterName(Displacement, "Displacement");
        this->setOutputParameterName(Pressure, "Pressure");
    };

    ///
    virtual ~FunctionDisplacementPressure()
    {
        BaseLib::releaseObject(_problem, _solution, _feObjects, _displacement /*, _strain, _stress*/);
    }

    /// initialize this process
    virtual bool initialize(const BaseLib::Options &option);

    /// finalize but nothing to do here
    virtual void finalize() {};

    ///
    virtual NumLib::IConvergenceCheck* getConvergenceChecker()
    {
        return &_checker;
    }

protected:
    virtual void updateOutputParameter(const NumLib::TimeStep &time);

    virtual MySolutionType* getSolution() {return _solution;};

    virtual void output(const NumLib::TimeStep &time);

private:
    MyVariable* getDisplacementComponent(MyVariable *u_x, MyVariable* u_y, MyVariable* u_z, const std::string &var_name);
    size_t getDisplacementComponentIndex(const std::string &var_name) const;
    //void calculateStressStrain();

    DISALLOW_COPY_AND_ASSIGN(FunctionDisplacementPressure);

private:
    MyProblemType* _problem;
    MySolutionType* _solution;
    FemLib::LagrangeFeObjectContainer* _feObjects;
    NumLib::DiscreteDataConvergenceCheck _checker;
    MyNodalFunctionVector* _displacement;
    std::vector<NodalPointScalarWrapper*> _vec_u_components;
    size_t _var_p_id;
//    FemLib::FEMIntegrationPointFunctionVector* _strain;
//    FemLib::FEMIntegrationPointFunctionVector* _stress;

};

#include "DisplacementPressure.hpp"

