/**
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file AbstractEOS_TotalDensityForm.h
 *
 * Created on 2015-04-28 by Yonghui Huang
 */
 
#ifndef ABSTRACTEOS_GLOBALCOMPLEMENTARYFORM_H 
#define ABSTRACTEOS_GLOBALCOMPLEMENTARYFORM_H 

#include "ChemLib/chemconst.h"

class AbstractEOS_GlobalComplementaryForm {
public:
	/**
	  * constructor, input is the number of unknowns
	  */
	AbstractEOS_GlobalComplementaryForm(std::size_t n_unknowns)
		:N(n_unknowns)
	{};

	/**
	  * destructor will be overriden depending on the real EOS class. 
	  */
	virtual ~AbstractEOS_GlobalComplementaryForm() {};

	/**
	  * set the environmental condition of P, T, X ... etc.
	  * these values are stored in the env_condition vector. 
	  */
	virtual void set_env_condition(ogsChem::LocalVector & env_conditon) = 0;

	virtual double getPcbySg(double Sg) = 0;
	virtual double Deriv_dPCdS(double Sg) = 0; 
		virtual double get_RHO_L_H(double Sg) = 0;
		virtual double get_RHO_G_H(double Sg) = 0;
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