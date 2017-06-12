/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file AbstractMeshBasedDiscreteLinearEquation.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"

#include "MathLib/LinAlg/Sparse/Sparsity.h"

#include "MeshLib/Core/IMesh.h"

#include "DiscreteLib/Core/IDiscreteLinearEquation.h"

namespace DiscreteLib
{

/**
 * \brief Abstract class for mesh-based discrete linear equations
 * 
 * - Mesh
 * - DoF manager
 * - Sparse pattern
 */
class AbstractMeshBasedDiscreteLinearEquation : public IDiscreteLinearEquation
{
public:
    /// \param  msh
    /// \param  dofManager
    AbstractMeshBasedDiscreteLinearEquation(const MeshLib::IMesh* msh, DofEquationIdTable* dofManager)
    : _msh(msh), _dofManager(dofManager), _sparsity(0)
    {
    }

    virtual ~AbstractMeshBasedDiscreteLinearEquation()
    {
        //Base::releaseObject(_dofManager);
        BaseLib::releaseObject(_sparsity);
    }

    /// return a mesh
    const MeshLib::IMesh* getMesh() const { return _msh; }

    /// return DoF map table
    virtual DofEquationIdTable* getDofMapManger() const { return _dofManager; }

    /// return a sparse pattern of equation matrix
    MathLib::RowMajorSparsity* getSparsity() const { return _sparsity; };

protected:
    void setSparsity(MathLib::RowMajorSparsity* sp) { _sparsity = sp; }

private:
    DISALLOW_COPY_AND_ASSIGN(AbstractMeshBasedDiscreteLinearEquation);

private:
    const MeshLib::IMesh* _msh;
    DofEquationIdTable* _dofManager;
    MathLib::RowMajorSparsity* _sparsity;
};


} //end
