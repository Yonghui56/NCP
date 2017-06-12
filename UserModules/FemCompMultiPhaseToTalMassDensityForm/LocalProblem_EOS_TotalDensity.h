/**
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LocalProblem_EOS.h
 *
 * Created on 2014-05-14 by Yonghui HUANG & Haibing Shao
 */

#ifndef LOCAL_PROBLEM_EOS_TOTALDENSITY_H
#define LOCAL_PROBLEM_EOS_TOTALDENSITY_H

#include "ChemLib/chemconst.h"
#include "AbstractEOS_TotalDensityForm.h"

/**
  * TODO: describe this class
  */
class LocalProblem_EOS_TotalDensity
{
public:
	/**
      * constructor of the class
      */
	LocalProblem_EOS_TotalDensity();
	
	/**
      * destructor of the class
      */
	~LocalProblem_EOS_TotalDensity(void);

	/**
	  * solve the EOS
	  * In our case, we have 2 inputs, P and X
	  * and 8 outputs, which are the secondary variables. 
	  * the initial guess will be embedded into the Output vector.
	  */
	void solve(ogsChem::LocalVector & Input, ogsChem::LocalVector & Output );
	//void deriv_solve(ogsChem::LocalVector & Input, ogsChem::LocalVector & Output)
	/**
	  * return the number of secondary variables
	  */
	std::size_t get_n() { return N; };

	/**
	  * TODO: describe the parameter
	  */
	const double eps = 1.0e-7;

	/**
	  * TODO: describe the parameter
	  */
	const double Iter_tol=200;

	/**
	  * TODO: describe the parameter
	  */
	const double theta = 0.7;

	/**
	  * tolorence of newton iteration
	  */
	const double tot = 1.0e-12; //a small tolorence is necessary for the iteration
	
	/**
	  * parameter to calc the N_G
	  */
	const double s= 1.0e-4;

private:
	/**
      * parameter to define the number of secondary variables
	  */
	std::size_t N;

	/**
	  * mean pressure ---primary variable
	  */
	double P_L;

	/**
	  * Total hydrogen mass density ---primary variable
	  */
	double X_L;

	/**
	  *  initial guess vector of the secondary variables 
	  */
	ogsChem::LocalVector U_ini;

	/**
	  * internal output vector
	  */
	ogsChem::LocalVector U_cur;

	/**
	  * Internal Residual vector
	  */
	ogsChem::LocalVector Res;

	/**
	  * Internal matrix of Jacobian 
	  */
	ogsChem::LocalMatrix Matrix_Jacobian;

	/**
	  * pointer the abstract EOS class
	  */
	AbstractEOS_TotalDensityForm * _EOS;

	/**
	  * minimization solve function. 
	  * it is capable of handling rank(J) < n. 
	  */
	void solve_minimization(ogsChem::LocalMatrix & J,
		ogsChem::LocalVector & b,
		ogsChem::LocalVector & dx);

	/**
	  * Newton iteration with line search
	  */
	void solve_LocalProblem_Newton_LineSearch(std::size_t & flag);
};

#endif
