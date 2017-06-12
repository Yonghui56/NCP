/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file OMPElementWiseLinearEquationAssembler.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#ifdef _OPENMP
#include <omp.h>
#endif

#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/Core/IDiscreteLinearEquationAssembler.h"

namespace DiscreteLib
{

/**
 * \brief Element-based discrete system assembler classes
 */
template <class T_UPDATER, class T_SOLVER>
class OMPElementWiseLinearEquationAssembler : public IDiscreteLinearEquationAssembler
{
public:
    typedef T_UPDATER UpdaterType;
    typedef T_SOLVER SolverType;

    ///
    explicit OMPElementWiseLinearEquationAssembler(UpdaterType* a) : _e_assembler(a) {};

    virtual ~OMPElementWiseLinearEquationAssembler() {};


    /// Conduct the element by element assembly procedure
    ///
    /// @param msh Mesh
    /// @param dofManager Dof map manager
    /// @param list_dofId List of Dof IDs used in this problem
    /// @param eqs Linear equation solver
    void assembly(const MeshLib::IMesh &msh, SolverType &eqs);

    virtual void assembly(const MeshLib::IMesh &msh, MathLib::ILinearEquation &eqs)
    {
        assembly(msh, *((SolverType*)&eqs));
    }

private:
    UpdaterType* _e_assembler;
};

template <class T1, class T2>
void OMPElementWiseLinearEquationAssembler<T1,T2>::assembly(const MeshLib::IMesh &msh, SolverType &eqs)
{
    const size_t n_ele = msh.getNumberOfElements();
    UpdaterType assembler(*_e_assembler);

#ifdef _OPENMP
    #pragma omp parallel for default(none), shared(std::cout, msh, eqs), firstprivate(assembler)
#endif
    for (size_t i=0; i<n_ele; i++) {
        MeshLib::IElement *e = msh.getElement(i);
        assembler.update(*e, eqs);
    }
};

//#ifdef _OPENMP
//        int thread_id = omp_get_thread_num();
//        std::cout << "thread " << thread_id << ": e_id = " << e->getID() << std::endl;
//#endif

}
