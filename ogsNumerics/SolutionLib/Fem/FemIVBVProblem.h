/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemIVBVProblem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <string>

#include "BaseLib/CodingTools.h"
#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "SolutionLib/Core/TimeSteppingProblem.h"
#include "FemVariable.h"

namespace SolutionLib
{

/**
 * \brief IVBV problems for FEM
 *
 * This class contains
 * - Variables
 * - Equations
 * - Reference to discrete system
 *
 * \tparam T_FEM_EQUATION   FEM equation
 */
template < class T_DIS_SYS, class T_FEM_EQUATION >
class FemIVBVProblem : public TimeSteppingProblem
{
public:
    typedef T_DIS_SYS MyDiscreteSystem;
    typedef FemVariable MyVariable;
    typedef T_FEM_EQUATION EquationType;

    ///
    explicit FemIVBVProblem(MyDiscreteSystem* dis)
     : _discrete_system(dis), _eqs(NULL)
    {
    }

    ///
    virtual ~FemIVBVProblem()
    {
        BaseLib::releaseObjectsInStdVector(_variables);
        BaseLib::releaseObject(_eqs);
    }

    /// get this discrete system
    MyDiscreteSystem* getDiscreteSystem() {return _discrete_system;};

    /// create FE approximation field
    MyVariable* addVariable(const std::string &name, FemLib::PolynomialOrder::type initial_order = FemLib::PolynomialOrder::Linear)
    {
        MyVariable* var = new MyVariable(_variables.size(), name, initial_order);
        FemLib::LagrangeFeObjectContainer* feContainer = new FemLib::LagrangeFeObjectContainer(_discrete_system->getMesh());
        var->setFeObjectContainer(feContainer);
        _variables.push_back(var);
        return var;
    }

    /// get a variable
    MyVariable* getVariable(size_t i) const { return _variables[i]; }

    /// get the number of variables
    size_t getNumberOfVariables() const { return _variables.size(); }

    /// create an equation
    EquationType* createEquation()
    {
        if (_eqs == NULL)
            _eqs = new EquationType();
        return _eqs;
    };

    /// get this equation
    EquationType* getEquation() const {return _eqs;};

private:
    DISALLOW_COPY_AND_ASSIGN(FemIVBVProblem);

private:
    MyDiscreteSystem* _discrete_system;
    std::vector<MyVariable*> _variables;
    EquationType* _eqs;
};

} //end
