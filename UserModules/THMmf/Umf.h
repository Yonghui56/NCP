/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Umf.h
 *
 * Created on 2012-11-15 by Norihiro Watanabe
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

#include "DeformationInPorousMediaLinearEQSLocalAssembler.h"

namespace THMmf
{

/**
 * \brief Monolithic poro-elasticity calculator using FEM
 */
template <class T_DISCRETE_SYSTEM, class T_LINEAR_SOLVER>
class Umf
: public ProcessLib::AbstractTransientProcess
{
public:
    enum Out
    {
        Displacement=0
    };

    typedef T_DISCRETE_SYSTEM MyDiscreteSystem;
    typedef T_LINEAR_SOLVER MyLinearSolver;
    // local assembler
    typedef DeformationInPorousMediaLinearEQSLocalAssembler MyLinearAssemblerTypeForPorousMedia;
    typedef NumLib::ElementWiseDimDependentTransientLinearEQSLocalAssembler MyLinearAssemblerType;
    // Equation definition
    typedef SolutionLib::TemplateFemEquation<
            MyDiscreteSystem,
            MyLinearSolver,
            MyLinearAssemblerType
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
    Umf()
    : ProcessLib::AbstractTransientProcess("DEFORMATION", 0, 1),
      _problem(0), _solution(0), _feObjects(0), _displacement(0)
    {
        // set default parameter name
        this->setOutputParameterName(Displacement, "Displacement");
    };

    ///
    virtual ~Umf()
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

    DISALLOW_COPY_AND_ASSIGN(Umf);

private:
    MyProblemType* _problem;
    MySolutionType* _solution;
    FemLib::IFeObjectContainer* _feObjects;
    NumLib::DiscreteDataConvergenceCheck _checker;
    MyNodalFunctionVector* _displacement;
    std::vector<NodalPointScalarWrapper*> _vec_u_components;

};

}

#include "Umf.tpp"

