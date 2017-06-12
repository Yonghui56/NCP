/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file XFEM_EXAMPLE_CRACK1.h
 *
 * Created on 2012-09-20 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/Serial/DiscreteSystem.h"
#include "NumLib/Function/DiscreteDataConvergenceCheck.h"
#include "FemLib/Function/FemNodalFunction.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "ProcessLib/AbstractTimeIndependentProcess.h"

namespace xfem
{

/**
 * \brief XFEM example: 2D-XFEM with sign-enrichment along crack
 *
 * 2D-XFEM problem with sign-enrichment along crack (NO branch-enrichment at crack-tip)
 * at y=0 in a square domain (x [-1,1], y [-1, 1]).
 * This example is created based on a MATLAB script (XFEM2dCrack_SignEnr.m)
 * provided by Thomas-Peter Fries, RWTH Aachen University. His license statement
 * is included in License.txt.
 */
class FunctionXFEM_EXAMPLE_CRACK1
: public ProcessLib::AbstractTimeIndependentProcess
{
public:

    typedef DiscreteLib::DiscreteSystem MyDiscreteSystem;
    typedef FemLib::FemNodalFunctionVector<MyDiscreteSystem>::type MyNodalFunctionVector;

    ///
    FunctionXFEM_EXAMPLE_CRACK1()
    : ProcessLib::AbstractTimeIndependentProcess("XFEM_EXAMPLE_CRACK1", 0, 0),
      _feObjects(NULL), _displacement(NULL), _exact_displacement(NULL), _msh(NULL), _dis(NULL)
    {
        // set default parameter name
    };

    ///
    virtual ~FunctionXFEM_EXAMPLE_CRACK1()
    {
        BaseLib::releaseObject(_feObjects, _displacement , _exact_displacement);
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

    int solveTimeStep(const NumLib::TimeStep &/*time*/);

    void finalizeTimeStep(const NumLib::TimeStep &/*time*/);

protected:
    virtual void updateOutputParameter(const NumLib::TimeStep &/*time*/) {};

private:
    DISALLOW_COPY_AND_ASSIGN(FunctionXFEM_EXAMPLE_CRACK1);

private:
    FemLib::LagrangeFeObjectContainer* _feObjects;
    NumLib::DiscreteDataConvergenceCheck _checker;
    MyNodalFunctionVector* _displacement;
    MyNodalFunctionVector* _exact_displacement;
    MeshLib::IMesh* _msh;
    MyDiscreteSystem* _dis;

};

}

