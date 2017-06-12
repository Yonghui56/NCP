/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IVBVProblem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "GeoLib/GeoObject.h"
#include "NumLib/Function/TXFunction.h"

namespace SolutionLib
{

class Variable
{

};


/**
 * \brief Initial value - boundary value (IVBV) problems
 *
 * IVBV problem has the following elements
 * - variables
 * - equations
 * - IC
 * - BC
 * - domain (space & time)?
 */
class IVBVProblem
{
public:
    ///
    virtual ~IVBVProblem() {};

    /// get the number of variables
    virtual size_t getNumberOfVariables() const = 0;


    /// set initial condition
    virtual void setIC(size_t var_type, NumLib::ITXFunction &ic) = 0;

    /// get initial condition
    virtual NumLib::ITXFunction* getIC(size_t var_type) const = 0;


    /// add a Dirichlet boundary condition
    virtual void addDirichletBC(size_t var_type,  GeoLib::GeoObject &geo, NumLib::ITXFunction &bc1) = 0;

    /// get the number of Dirichlet BCs
    virtual size_t getNumberOfDirichletBC(size_t var_type) const = 0;

    /// get the Dirichlet BC
    virtual NumLib::ITXFunction* getDirichletBC(size_t var_type, size_t bc_id) const = 0;


    /// add a Neumann boundary condition
    virtual void addNeumannBC(size_t var_type,  GeoLib::GeoObject &geo, NumLib::ITXFunction &bc2) = 0;

    /// get the number of Neumann BCs
    virtual size_t getNumberOfNeumannBC(size_t var_type) const = 0;

    /// get the Neumann BC
    virtual NumLib::ITXFunction* getNeumannBC(size_t var_type, size_t bc_id) const = 0;
};

}
