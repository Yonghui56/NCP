/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemReactOPSsolution.h
 *
 * Created on 2013-03-28 by Haibing Shao
 */

#ifndef FEM_REACT_OPS_SOLUTION_H
#define FEM_REACT_OPS_SOLUTION_H

#include <vector>
#include <string>

#include "BaseLib/CodingTools.h"
#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "SolutionLib/Core/TimeSteppingProblem.h"
#include "SolutionLib/Fem/FemVariable.h"
#include "ChemLib/chemEqReactSysActivity.h"

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
template < class T_DIS_SYS >
class FemReactOPSsolution : public TimeSteppingProblem
{
public:
    typedef T_DIS_SYS MyDiscreteSystem;
    typedef FemVariable MyVariable;

    /**
      * constructor
      */
    explicit FemReactOPSsolution(MyDiscreteSystem* dis, ogsChem::chemEqReactSysActivity* mychemEqReactSys  )
     : _discrete_system(dis), _EqReactSys(mychemEqReactSys)
    {
    }

    /**
      * destructor
      */
    virtual ~FemReactOPSsolution()
    {
        BaseLib::releaseObjectsInStdVector(_variables);
        _EqReactSys = NULL;
        _discrete_system = NULL;
    }

    /**
      * get this discrete system
      */
    MyDiscreteSystem* getDiscreteSystem() {return _discrete_system;};

    /**
      * create FE approximation field
      */
    MyVariable* addVariable(const std::string &name)
    {
        _variables.push_back(new MyVariable(_variables.size(), name));
        return _variables.back();
    }

    /**
      * get a variable
      */
    MyVariable* getVariable(size_t i) const { return _variables[i]; }

    /**
      * get the number of variables
      */
    size_t getNumberOfVariables() const { return _variables.size(); }

	/**
      * get the pointer to the equilibrium reaction system
      */
    ogsChem::chemEqReactSysActivity* getEqReactSys() const { return _EqReactSys; }

private:
    DISALLOW_COPY_AND_ASSIGN(FemReactOPSsolution);

private:
    /**
      * discretization 
      */ 
    MyDiscreteSystem* _discrete_system;

    /**
      * vector of variables
      */ 
    std::vector<MyVariable*> _variables;

    /**
      * pointer to equilibrium reaction system
      */ 
	ogsChem::chemEqReactSysActivity* _EqReactSys; // TO BE CHANGED. 

};

} //end

#endif  // end of ifndef
