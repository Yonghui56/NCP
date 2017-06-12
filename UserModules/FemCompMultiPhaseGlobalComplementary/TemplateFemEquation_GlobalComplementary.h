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
//#include "NumLib/TransientAssembler/DummyElementWiseTransientLinearEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/DummyElementWiseTransientJacobianLocalAssembler.h"
#include "NumLib/TransientAssembler/DummyElementWiseTransientResidualLocalAssembler.h"
//#include "FemLinearEQS.h"
#include "FemResidualEQS_GlobalComplementary.h"
#include "FemDxEQS_GlobalComplementary.h"

#include "SolutionLib/DataType.h"
#include "SolutionLib/Core/AbstractTimeSteppingAlgorithm.h"
#include "SolutionLib/Fem/TemplateFemEquation.h"
#include "SolutionLib/Fem/FemDirichletBC.h"
#include "SolutionLib/Fem/FemNeumannBC.h"
#include "SolutionLib/Fem/FemLinearEQS.h"


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
	class T_LINEAR_SOLVER,
	class T_LOCAL_ASSEMBLER_LINEAR,
	class T_FUNCTION_DATA,
	class T_LOCAL_ASSEMBLER_RESIDUAL = NumLib::DummyElementWiseTransientResidualLocalAssembler,
	class T_LOCAL_ASSEMBLER_JACOBIAN = NumLib::DummyElementWiseTransientJacobianLocalAssembler
>
class TemplateFemEquation_GlobalComplementary
{
public:
    typedef T_DIS_SYS MyDiscreteSystem;
    typedef T_LINEAR_SOLVER LinearSolverType;
    typedef T_LOCAL_ASSEMBLER_LINEAR LinearAssemblerType;
	typedef T_LOCAL_ASSEMBLER_RESIDUAL ResidualAssemblerType;
	typedef T_LOCAL_ASSEMBLER_JACOBIAN JacobianAssemblerType;
    typedef SolutionLib::TemplateTransientLinearFEMFunction<MyDiscreteSystem,LinearSolverType,LinearAssemblerType> LinearEQSType;
	typedef TemplateTransientResidualFEMFunction_Global_Complementary<MyDiscreteSystem, T_FUNCTION_DATA, ResidualAssemblerType> ResidualEQSType;
	typedef TemplateTransientDxFEMFunction_Global_Complementary<MyDiscreteSystem, LinearSolverType, T_FUNCTION_DATA, JacobianAssemblerType> DxEQSType;

	TemplateFemEquation_GlobalComplementary():
		_linear_assembler(NULL),
		_residual_assembler(NULL),
		_jacobian_assembler(NULL)
    {};

    ///
	~TemplateFemEquation_GlobalComplementary()
    {
		BaseLib::releaseObject(_residual_assembler, _jacobian_assembler);
    }
	void initialize(
		LinearAssemblerType *linear_assembly,
		ResidualAssemblerType *residual_assembly,
		JacobianAssemblerType *jacobian_assembly
		)
	{
		_linear_assembler = linear_assembly;
		_residual_assembler = residual_assembly;
		_jacobian_assembler = jacobian_assembly;
	}

	void initialize(
		LinearAssemblerType *linear_assembly
		)
	{
		_linear_assembler = linear_assembly;
	}


	///
	LinearAssemblerType* getLinearAssembler() const { return _linear_assembler; }

	///
	ResidualAssemblerType* getResidualAssembler() const { return _residual_assembler; }

	///
	JacobianAssemblerType* getJacobianAssembler() const { return _jacobian_assembler; }

private:
	LinearAssemblerType* _linear_assembler;
	ResidualAssemblerType* _residual_assembler;
	JacobianAssemblerType* _jacobian_assembler;
};


