/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file OMPElementWiseVectorAssembler.h
 *
 * Created on 2012-08-20 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "MeshLib/Core/IMesh.h"
#include "MathLib/DataType.h"
#include "DiscreteLib/Core/IDiscreteVectorAssembler.h"

namespace MeshLib 
{
class IMesh;
}

namespace DiscreteLib
{

/**
 * \brief Element-based discrete vector assembler classes
 */
template <class T_VALUE, class T_UPDATER>
class OMPElementWiseVectorAssembler : public IDiscreteVectorAssembler<T_VALUE>
{
public:
    typedef typename IDiscreteVectorAssembler<T_VALUE>::VectorType GlobalVectorType;
    typedef T_UPDATER UpdaterType;

    ///
    explicit OMPElementWiseVectorAssembler(UpdaterType* a) : _e_assembler(a) {};
    ///
    virtual ~OMPElementWiseVectorAssembler() {};

    /// Conduct the element by element assembly procedure
    ///
    /// @param msh Mesh
    /// @param dofManager Dof map manager
    /// @param vec Discrete vector
    void assembly(const MeshLib::IMesh &msh, GlobalVectorType &globalVec);

private:
    UpdaterType* _e_assembler;
};


template <class T1, class T2>
void OMPElementWiseVectorAssembler<T1,T2>::assembly(const MeshLib::IMesh &msh, GlobalVectorType &globalVec)
{
    const size_t n_ele = msh.getNumberOfElements();
    UpdaterType assembler(*_e_assembler);

#ifdef _OPENMP
    #pragma omp parallel for default(none), shared(msh, globalVec), firstprivate(assembler)
#endif
    for (size_t i=0; i<n_ele; i++) {
        MeshLib::IElement *e = msh.getElement(i);
        assembler.update(*e, globalVec);
    }
};

}
