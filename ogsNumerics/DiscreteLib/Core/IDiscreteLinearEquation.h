/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IDiscreteLinearEquation.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "IDiscreteObject.h"
#include "IDiscreteVector.h"

namespace DiscreteLib
{

// forward declaration
class DofEquationIdTable;
class IDiscreteLinearEquationAssembler;

/** 
 * \brief Interface for discrete linear equations
 */
class IDiscreteLinearEquation : public IDiscreteObject
{
public:
    typedef IDiscreteVector<double> GlobalVectorType;

    ///
    virtual ~IDiscreteLinearEquation() {};
    /// 
    virtual void initialize() = 0;
    /// construct 
    virtual void construct(IDiscreteLinearEquationAssembler& assemler) = 0;
    /// solve
    virtual void solve() = 0;
    /// get solution
    virtual void getGlobalX(std::vector<double> &x) = 0;
    ///
    virtual double* getLocalX() = 0;
    ///
    virtual void getX(GlobalVectorType &v) = 0;
    ///
    virtual void setX(const GlobalVectorType &v) = 0;
    /// get a Dof map manager
    virtual DofEquationIdTable* getDofMapManger() const = 0;
    /// set prescribed dof
    virtual void setPrescribedDoF(size_t dofId, std::vector<size_t> &list_discrete_pt_id, std::vector<double> &list_prescribed_values) = 0;
    /// set additional RHS values
    virtual void addRHS(size_t dofId, std::vector<size_t> &list_discrete_pt_id, std::vector<double> &list_rhs_values, double fkt) = 0;
    /// set additional RHS values
    virtual void addRHS(const GlobalVectorType &v, double fkt) = 0;
};


}
