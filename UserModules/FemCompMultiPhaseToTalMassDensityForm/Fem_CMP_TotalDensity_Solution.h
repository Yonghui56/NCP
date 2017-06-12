/**
* Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.com/LICENSE.txt
*
*
* \file Fem_CMP_TotalDensity_Solution.h
*
* Created on 2015-04-20 by Yonghui HUANG
*/

#ifndef FEM_CMP_TOTALDENSITY_SOLUTION_H
#define FEM_CMP_TOTALDENSITY_SOLUTION_H

#include <vector>
#include <string>

#include "BaseLib/CodingTools.h"
#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "SolutionLib/Core/TimeSteppingProblem.h"
#include "SolutionLib/Fem/FemVariable.h"

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
	class Fem_CMP_TotalDensity_Solution : public TimeSteppingProblem
	{
	public:
		typedef T_DIS_SYS MyDiscreteSystem;
		typedef FemVariable MyVariable;

		/**
		  * constructor
		  */
		explicit Fem_CMP_TotalDensity_Solution(MyDiscreteSystem* dis)
			: _discrete_system(dis)
		{
		}

		/**
		  * destructor
		  */
		virtual ~Fem_CMP_TotalDensity_Solution()
		{
			BaseLib::releaseObjectsInStdVector(_variables);
			_discrete_system = NULL;
		}

		/**
		  * get this discrete system
		  */
		MyDiscreteSystem* getDiscreteSystem() { return _discrete_system; };

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

	private:
		DISALLOW_COPY_AND_ASSIGN(Fem_CMP_TotalDensity_Solution);

	private:
		/**
		  * discretization
		  */
		MyDiscreteSystem* _discrete_system;

		/**
		  * vector of variables
		  */
		std::vector<MyVariable*> _variables;

	};

} //end

#endif  // end of ifndef
