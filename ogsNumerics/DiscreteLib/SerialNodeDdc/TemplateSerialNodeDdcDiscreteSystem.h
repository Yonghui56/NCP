/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NodeDDCSerialSharedDiscreteSystem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "DiscreteLib/Core/DiscreteEnums.h"
#include "DiscreteLib/DDC/DecomposedDomain.h"
#include "DiscreteLib/DDC/DecomposedMesh.h"
#include "DiscreteLib/Serial/SequentialElementWiseLinearEquationAssembler.h"
#include "DiscreteLib/Serial/SequentialElementWiseVectorAssembler.h"
#include "DecomposedVector.h"

namespace DiscreteLib
{

/**
 * \brief template class for serial node-based decomposition system
 * 
 * \tparam T_EQS Linear equation class
 */
template <template <typename, typename> class T_EQS, DiscreteSystemType::type T_TYPE>
class TemplateSerialNodeDdcDiscreteSystem : public IDiscreteSystem
{
public:
    template<typename T_LINEAR_SOLVER, typename T_SPARSITY_BUILDER>
    struct MyLinearEquation
    {
        typedef T_EQS<T_LINEAR_SOLVER, T_SPARSITY_BUILDER> type;
    };
    template<typename T>
    struct MyVector
    {
        typedef DecomposedVector<T> type;
    };

    template <typename T_UPDATER, typename T_LINEAR_SOLVER>
    struct MyLinearEquationAssembler
    {
        typedef SequentialElementWiseLinearEquationAssembler<T_UPDATER, T_LINEAR_SOLVER> type;
    };

    template <typename T_VALUE, typename T_UPDATER>
    struct MyVectorAssembler
    {
        typedef SequentialElementWiseVectorAssembler<T_VALUE, T_UPDATER> type;
    };
    
    ///
    static DiscreteSystemType::type getSystemType() { return T_TYPE;};

    /**
     * 
     * @param msh
     */
    explicit TemplateSerialNodeDdcDiscreteSystem(const MeshLib::IMesh *msh)
    : _ddc_global(NULL), _msh(msh)
    { };
    
    /// 
    virtual ~TemplateSerialNodeDdcDiscreteSystem()
    {
        BaseLib::releaseObject(_msh);
    };

    /**
     * 
     * @param ddc_global
     */
    void initialize(DecomposedDomain* ddc_global)
    {
        _ddc_global = ddc_global;
        //_msh = new DecomposedMesh(ddc_global->getID(), ddc_global);
    }
    
    /**
     * 
     * @return
     */
    virtual MeshLib::IMesh* getMesh() const { return (MeshLib::IMesh*)_msh; };

    /// create a new linear equation
    /// \tparam T_LINEAR_SOLVER     Linear solver class
    /// \tparam T_SPARSITY_BUILDER  Sparse pattern builder class
    /// \param linear_solver
    /// \param dofManager
    template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
    typename MyLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>::type* createLinearEquation(T_LINEAR_SOLVER* linear_solver, DofEquationIdTable* dofManager)
    {
        return MyLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>::type::createInstance(*this, _ddc_global, linear_solver, dofManager);
    }

    //// create a vector
    /// \tparam T   data type
    /// \param n    vector size
    template<typename T>
    typename MyVector<T>::type* createVector(const size_t &n)
    {
        std::vector<size_t> list_range_begin;
        for (size_t i=0; i<_ddc_global->getNumberOfSubDomains(); i++) {
            size_t cnt = (size_t)((double)n / (double)_ddc_global->getNumberOfSubDomains() * (double)i);
            list_range_begin.push_back(cnt);
        }

        return MyVector<T>::type::createInstance(*this, n, list_range_begin);
    };

private:
    DISALLOW_COPY_AND_ASSIGN(TemplateSerialNodeDdcDiscreteSystem);

private:
    DecomposedDomain* _ddc_global;
    const MeshLib::IMesh* _msh;
};

} //end
