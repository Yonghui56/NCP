/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemEquation.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "NumLib/TransientAssembler/DummyElementWiseTransientLinearEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/DummyElementWiseTransientJacobianLocalAssembler.h"
#include "NumLib/TransientAssembler/DummyElementWiseTransientResidualLocalAssembler.h"
//#include "FemLinearEQS.h"
#include "FemResidualEQS_GIAReduct.h"
#include "FemDxEQS_GIAReduct.h"

/**
 * \brief FEM equations
 *
 * - Linear equation: Ax = b
 * - Residual equation: r = Ax - b
 * - Dx equation for Newton: J dx = -r
 *
 * \tparam T_LOCAL_ASSEMBLER_LINEAR     Local assembler for linear equation
 * \tparam T_LOCAL_ASSEMBLER_LINEAR     Local assembler for residual vector
 * \tparam T_LOCAL_ASSEMBLER_LINEAR     Local assembler for Jacobian matrix
 */
template
    <
    class T_DIS_SYS,
    class T_USER_FUNCTION_DATA,
    class T_LINEAR_SOLVER,
    class T_LOCAL_ASSEMBLER_LINEAR = NumLib::DummyElementWiseTransientLinearEQSLocalAssembler,
    class T_LOCAL_ASSEMBLER_RESIDUAL = NumLib::DummyElementWiseTransientResidualLocalAssembler
//    class T_LOCAL_ASSEMBLER_JACOBIAN = NumLib::DummyElementWiseTransientJacobianLocalAssembler
    >
class TemplateFemEquation_GIAReduct
{
public:
    typedef T_DIS_SYS MyDiscreteSystem;
    typedef T_LINEAR_SOLVER LinearSolverType;
    typedef T_LOCAL_ASSEMBLER_LINEAR LinearAssemblerType;
    typedef SolutionLib::TemplateTransientLinearFEMFunction<MyDiscreteSystem,LinearSolverType,LinearAssemblerType> LinearEQSType;
    typedef TemplateTransientResidualFEMFunction_GIA_Reduct<MyDiscreteSystem,T_USER_FUNCTION_DATA> ResidualEQSType;
    typedef TemplateTransientDxFEMFunction_GIA_Reduct<MyDiscreteSystem,LinearSolverType,T_USER_FUNCTION_DATA> DxEQSType;

    TemplateFemEquation_GIAReduct()
    {};

    ///
    ~TemplateFemEquation_GIAReduct()
    {
    }

private:

};


