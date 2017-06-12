/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IDiscreteVectorAssembler.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "IDiscreteVector.h"

namespace MeshLib 
{
class IMesh;
}

namespace DiscreteLib
{
class DofEquationIdTable;

/**
 * \brief Interface of discrete vector assembler classes
 */
template <class T>
class IDiscreteVectorAssembler
{
public:
    typedef IDiscreteVector<T> VectorType;

    virtual ~IDiscreteVectorAssembler() {};

    /// Conduct the element by element assembly procedure
    ///
    /// @param msh Mesh
    /// @param vec Discrete vector
    virtual void assembly( const MeshLib::IMesh &msh, const DofEquationIdTable &dofEquationIdTable, IDiscreteVector<T> &vec) = 0;
};

}
