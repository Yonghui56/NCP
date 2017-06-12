/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SequentialElementWiseLinearEquationAssembler.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "DiscreteLib/Core/IDiscreteLinearEquationAssembler.h"
#include "DiscreteLib/Core/IElemenetWiseLinearEquationLocalAssembler.h"
#include "DiscreteLib/Utils/DofEquationIdTable.h"
//#include "ogs6THMC/ogs6THMC/UserModules/FemCompMultiPhaseGlobalComplementary/FunctionCMP_GlobalComplementaryForm.h"

namespace MeshLib 
{
class IMesh;
}

namespace DiscreteLib
{

/**
 * \brief Element-based linear equation assembler
 */
template <class T_UPDATER, class T_SOLVER>
class SequentialElementWiseLinearEquationAssembler : public IDiscreteLinearEquationAssembler
{
public:
    typedef T_UPDATER UpdaterType;
    typedef T_SOLVER SolverType;

    /**
     *
     * @param updater
     */
    explicit SequentialElementWiseLinearEquationAssembler(UpdaterType* updater) : _e_assembler(updater) {};

    ///
    virtual ~SequentialElementWiseLinearEquationAssembler(){};

    /// Conduct the element by element assembly procedure
    ///
    /// @param msh Mesh
    /// @param list_dofId List of Dof IDs used in this problem
    /// @param eqs Linear equation solver
    void assembly(const MeshLib::IMesh &msh, const DofEquationIdTable &dofEquationIdTable, SolverType &eqs);

    virtual void assembly(const MeshLib::IMesh &msh, const DofEquationIdTable &dofEquationIdTable, MathLib::ILinearEquation &eqs)
    {
        assembly(msh, dofEquationIdTable, *((SolverType*)&eqs));
    }

private:
    UpdaterType* _e_assembler;
};

template <class T1, class T2>
void SequentialElementWiseLinearEquationAssembler<T1, T2>::assembly(const MeshLib::IMesh &msh, const DofEquationIdTable &dofEquationIdTable, SolverType &eqs)
{
    const size_t n_ele = msh.getNumberOfElements();
    for (size_t i=0; i<n_ele; i++) {

    	MeshLib::IElement *e = msh.getElement(i);

        _e_assembler->update(*e, dofEquationIdTable, eqs);

    }
};

}
