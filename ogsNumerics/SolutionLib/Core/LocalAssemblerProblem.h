/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LocalAssemblerProblem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once


namespace SolutionLib
{

template
    <
    class T_LOCAL_ASSEMBLER_LINEAR,
    class T_LOCAL_ASSEMBLER_RESIDUAL,
    class T_LOCAL_ASSEMBLER_JACOBIAN
    >
class LocalAssemblerProblem
{
public:
    typedef T_LOCAL_ASSEMBLER_LINEAR LinearAssemblerType;
    typedef T_LOCAL_ASSEMBLER_RESIDUAL ResidualAssemblerType;
    typedef T_LOCAL_ASSEMBLER_JACOBIAN JacobianAssemblerType;

    LocalAssemblerProblem(
            LinearAssemblerType *linear_assembly,
            ResidualAssemblerType *residual_assembly,
            JacobianAssemblerType *jacobian_assembly
    ) :
    _linear_assembler(linear_assembly),
    _residual_assembler(residual_assembly),
    _jacobian_assembler(jacobian_assembly)
    {};

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


}
