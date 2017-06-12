/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file StressStrain.h
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "BaseLib/CodingTools.h"
#include "BaseLib/Options.h"
#include "DiscreteLib/Serial/DiscreteSystem.h"
#include "FemLib/Function/FemFunction.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "NumLib/Function/DiscreteDataConvergenceCheck.h"
#include "ProcessLib/AbstractTimeIndependentProcess.h"

/**
 * \brief Stress, strain evaluator based on total displacements
 */
template <class T_DISCRETE_SYSTEM>
class FunctionElementStressStrain
    : public ProcessLib::AbstractTimeIndependentProcess
{
public:
    enum In { Displacement=0 };
    enum Out { Strain=0, Stress=1};

    typedef T_DISCRETE_SYSTEM MyDiscreteSystem;
    typedef typename FemLib::FemNodalFunctionVector<MyDiscreteSystem>::type MyNodalFunctionVector;
    typedef typename FemLib::FEMIntegrationPointFunctionVector<MyDiscreteSystem>::type MyIntegrationPointFunctionVector;
    typedef typename NumLib::TXVectorFunctionAsColumnData<MyIntegrationPointFunctionVector> IntegrationPointScalarWrapper;
    typedef typename NumLib::TXVectorFunctionAsColumnData<MyNodalFunctionVector> NodalPointScalarWrapper;

    FunctionElementStressStrain()
    : ProcessLib::AbstractTimeIndependentProcess("ELEMENT_STRESS_STRAIN", 1, 2), _dis(0), _strain(0), _stress(0), _feObjects(0)
    {
        this->setInputParameterName(Displacement, "Displacement");
        this->setOutputParameterName(Strain, "Strain");
        this->setOutputParameterName(Stress, "Stress");
    };

    virtual ~FunctionElementStressStrain()
    {
        BaseLib::releaseObject(_strain, _stress, _feObjects);
        BaseLib::releaseObjectsInStdVector(_vec_strain_components);
        BaseLib::releaseObjectsInStdVector(_vec_stress_components);
    };

    virtual bool initialize(const BaseLib::Options &op);

    virtual void finalize() {};

    int solveTimeStep(const NumLib::TimeStep &/*time*/);

    void finalizeTimeStep(const NumLib::TimeStep &/*time*/);

    ///
    virtual NumLib::IConvergenceCheck* getConvergenceChecker() { return &_checker; };

private:
    size_t getNumberOfStrainComponents() const;

    DISALLOW_COPY_AND_ASSIGN(FunctionElementStressStrain);

private:
    MyDiscreteSystem* _dis;
    MyIntegrationPointFunctionVector* _strain;
    MyIntegrationPointFunctionVector* _stress;
    std::vector<IntegrationPointScalarWrapper*> _vec_strain_components;
    std::vector<IntegrationPointScalarWrapper*> _vec_stress_components;
    NumLib::DiscreteDataConvergenceCheck _checker;
    FemLib::LagrangeFeObjectContainer* _feObjects;
};

#include "ElementStressStrain.hpp"
