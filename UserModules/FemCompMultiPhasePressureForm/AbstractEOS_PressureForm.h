/**
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file AbstractEOS.h
 *
 * Created on 2014-05-28 by Haibing Shao
 */
 
#ifndef ABSTRACTEOS_PRESSUREFORM_H 
#define ABSTRACTEOS_PRESSUREFORM_H 

#include "ChemLib/chemconst.h"

class AbstractEOS_PressureForm {
public:
	/**
	  * constructor, input is the number of unknowns
	  */
	AbstractEOS_PressureForm(std::size_t n_unknowns)
		:N(n_unknowns)
	{};

	/**
	  * destructor will be overriden depending on the real EOS class. 
	  */
	virtual ~AbstractEOS_PressureForm() {};

	/**
	  * set the environmental condition of P, T, X ... etc.
	  * these values are stored in the env_condition vector. 
	  */
	virtual void set_env_condition(ogsChem::LocalVector & env_conditon) = 0;

	//virtual void eval(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalVector & res) = 0;

	/**
	  * calculate the Jacobian matrix based on the values of unknowns given
	  */
	virtual void calc_Jacobian(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalMatrix & J) = 0;
	/**
	* calculate REST SECONDARY VARIABLES BASED ON THE SATURATION AND p_l
	*/
	//virtual void calc_Rest_SecVar(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalVector & vec_rest_var) = 0;
	/**
	  * get the number of unknowns and number of equations
	  */
	std::size_t get_n() { return N; };

private: 

	/**
	  * number of unknowns and governing equtions
	  */
	const std::size_t N; 


}; 

#endif