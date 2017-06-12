/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NodalStressStrain.h
 *
 * Created on 2012-08-14 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "BaseLib/CodingTools.h"
#include "BaseLib/Options.h"
#include "DiscreteLib/Serial/DiscreteSystem.h"
#include "FemLib/Function/FemFunction.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "NumLib/Function/TXVectorFunctionAsColumnData.h"
#include "NumLib/Function/DiscreteDataConvergenceCheck.h"
#include "ProcessLib/AbstractTimeIndependentProcess.h"

/**
 * \brief Stress, strain evaluator based on total displacements
 */
template <class T_DISCRETE_SYSTEM>
class FunctionNodalStressStrain
    : public ProcessLib::AbstractTimeIndependentProcess
{
public:
    enum In { GpStrain=0, GpStress=1};
    enum Out { NodStrain=0, NodStress=1};

    typedef T_DISCRETE_SYSTEM MyDiscreteSystem;
    typedef typename FemLib::FemNodalFunctionVector<MyDiscreteSystem>::type MyNodalFunctionVector;
    typedef typename FemLib::FEMIntegrationPointFunctionVector<MyDiscreteSystem>::type MyIntegrationPointFunctionVector;
    typedef typename NumLib::TXVectorFunctionAsColumnData<MyIntegrationPointFunctionVector> IntegrationPointScalarWrapper;
    typedef typename NumLib::TXVectorFunctionAsColumnData<MyNodalFunctionVector> NodalPointScalarWrapper;

    FunctionNodalStressStrain()
    : ProcessLib::AbstractTimeIndependentProcess("NODAL_STRESS_STRAIN",2,2), _dis(0), _nodal_strain(0), _nodal_stress(0)
    {
        this->setInputParameterName(GpStrain, "GpStrain");
        this->setInputParameterName(GpStress, "GpStress");
        this->setOutputParameterName(NodStrain, "NodStrain");
        this->setOutputParameterName(NodStress, "NodStress");
    };

    virtual ~FunctionNodalStressStrain()
    {
        BaseLib::releaseObject(_nodal_strain, _nodal_stress, _feObjects);
        BaseLib::releaseObjectsInStdVector(_vec_nodal_strain_components);
        BaseLib::releaseObjectsInStdVector(_vec_nodal_stress_components);
    };

    virtual bool initialize(const BaseLib::Options &op);

    virtual void finalize() {};

    int solveTimeStep(const NumLib::TimeStep &/*time*/);

    void finalizeTimeStep(const NumLib::TimeStep &/*time*/);

    ///
    virtual NumLib::IConvergenceCheck* getConvergenceChecker() { return &_checker; };

private:
    size_t getNumberOfStrainComponents() const;

    DISALLOW_COPY_AND_ASSIGN(FunctionNodalStressStrain);

private:
    MyDiscreteSystem* _dis;
    MyNodalFunctionVector* _nodal_strain;
    MyNodalFunctionVector* _nodal_stress;
    std::vector<NodalPointScalarWrapper*> _vec_nodal_strain_components;
    std::vector<NodalPointScalarWrapper*> _vec_nodal_stress_components;
    FemLib::LagrangeFeObjectContainer* _feObjects;
    NumLib::DiscreteDataConvergenceCheck _checker;
};

#include "NodalStressStrain.hpp"
