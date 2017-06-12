/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file OMPDiscreteSystem
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#ifdef _OPENMP
#include <omp.h>
#endif

#include "logog.hpp"

#include "BaseLib/CodingTools.h"

#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "DiscreteLib/Core/DiscreteEnums.h"
#include "DiscreteLib/Serial/DiscreteLinearEquation.h"
#include "DiscreteLib/Serial/DiscreteVector.h"
#include "OMPElementWiseLinearEquationAssembler.h"
#include "OMPElementWiseVectorAssembler.h"

namespace MeshLib
{
class IMesh;
}

namespace DiscreteLib
{

class DofEquationIdTable;

/**
 * \brief Discrete system using OpenMP parallelization for assembly
 *
 */
class OMPDiscreteSystem : public IDiscreteSystem
{
public:
    template<typename T_LINEAR_SOLVER, typename T_SPARSITY_BUILDER>
    struct MyLinearEquation
    {
        typedef DiscreteLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER> type;
    };
    template<typename T>
    struct MyVector
    {
        typedef DiscreteVector<T> type;
    };
    template <class T_UPDATER, typename T_LINEAR_SOLVER>
    struct MyLinearEquationAssembler
    {
        typedef OMPElementWiseLinearEquationAssembler<T_UPDATER, T_LINEAR_SOLVER> type;
    };
    template <typename T_VALUE, typename T_UPDATER>
    struct MyVectorAssembler
    {
        typedef OMPElementWiseVectorAssembler<T_VALUE, T_UPDATER> type;
    };

    /**
     * 
     * @return
     */
    static DiscreteSystemType::type getSystemType() { return DiscreteSystemType::OpenMP;};

    /// constructor
    /// \param msh  a mesh to represent discrete systems by nodes or elements or edges
    explicit OMPDiscreteSystem(const MeshLib::IMesh* msh) : _msh(msh)
    {
#ifdef _OPENMP
        INFO("OMPDiscreteSystem() is created. The number of processors is %d", omp_get_num_procs());
#else
        WARN("***Warning: OMPDiscreteSystem() is used however OpenMP was disabled when compiling.");
#endif
    };

    ///
    virtual ~OMPDiscreteSystem()
    {
        BaseLib::releaseObject(_msh);
    };

    /// return this mesh
    virtual MeshLib::IMesh* getMesh() const {return (MeshLib::IMesh*)_msh;};

    /// create a new linear equation
    /// @tparam T_LINEAR_EQUATION
    /// @tparam T_LINEAR_SOLVER
    /// @tparam T_SPARSITY_BUILDER
    /// @param linear_solver         Linear solver
    /// @param dofManager            Equation index table
    template <class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
    typename MyLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>::type* createLinearEquation(T_LINEAR_SOLVER* linear_solver, DofEquationIdTable* dofManager)
    {
        return MyLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>::type::createInstance(*this, _msh, linear_solver, dofManager);
    }

    /// create a new vector
    /// @param n    vector length
    /// @return vector object
    template<typename T>
    typename MyVector<T>::type* createVector(size_t n)
    {
        return MyVector<T>::type::createInstance(*this, n);
    };

private:
    DISALLOW_COPY_AND_ASSIGN(OMPDiscreteSystem);

protected:
    const MeshLib::IMesh* _msh;
};

} //end
