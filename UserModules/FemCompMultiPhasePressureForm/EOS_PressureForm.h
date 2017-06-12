/**
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file EOS_H2_H2O_ISOT.h
 *
 * Created on 2014-05-12 by Yonghui HUANG
 */

#pragma once


#include "math.h"
#include <algorithm>
#include "AbstractEOS_PressureForm.h"
#include <cmath>

/**
 *  EOS formulation for calculating each secondary variable
 *  PC PG PL N_L N_G Kr_g Kr_l X_m X_M Sg
 */
class EOS_PressureForm : public AbstractEOS_PressureForm
{
public:
    
	/**
	  * Constructor
	  */
	EOS_PressureForm() : AbstractEOS_PressureForm(1)
    {
    };

	/**
	  * DESTRUCTOR
	  */
	virtual ~EOS_PressureForm()
    {

    };

	/**
	  * realization of the eval function. 
	  * this function will evaluate the resiudal of the governing equation, 
	  * based on the unknown values given. 
	  */

	/**
	  * realization of the calc_Jacobian function.
	  * this function will evaluate the Jacobian matrix,
	  * based on the unknown values given.
	  */
	virtual void calc_Jacobian(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalMatrix & J)
	{
		assert(J.cols() == this->get_n() && J.rows() == this->get_n() );

		
	};

	virtual void set_env_condition(ogsChem::LocalVector & env_conditon)
	{
		// mean pressure
		P_L = env_conditon(0);
		// molar fraction of the light component
		X_L = env_conditon(1);
	}; 

	
	virtual double getS_bar(double Sg, double res_S_l = 0.40, double res_S_g = 0.0, double es=0.59)
	{
		double S_bar = 0.0;
		
		S_bar = res_S_g + (1 - es)*(Sg - res_S_g) + 0.5*es*(1 - res_S_g - res_S_l);
		return S_bar;
	}

	virtual double getSat_byPC( double n = 1.49, double m = 1.0 - 1.0 / 1.49)//for MoMaS benchmark2 n = 1.49
	{
		double P_c = 0.0;
		double Sg = 0.0;
		P_c = (X_L / C_h) - P_L;
		if (P_c > 0){
			Sg = 1 - (((1 - S_gr - S_lr) / pow((pow((P_c / P_r), n) + 1), m)) + S_lr);
		}
		else
			Sg = 0.0;//just for test
		return Sg;

	}

	virtual double get_diriv_PC(double n = 1.49, double m = 1.0 - 1.0 / 1.49)//for MoMaS benchmark2 n = 1.49
	{
		double dPC(0.0);
		double P_c = 0.0;
		P_c = (X_L / C_h) - P_L;
		dPC = m*n*(1 - S_gr - S_lr)*(1 / P_r)*(pow((P_c / P_r), (n - 1)))*pow(((pow((P_c / P_r), n)) + 1), (-m - 1));
		if (P_c<=0)
		{
			dPC = 0.0;//just for test
		}
		return dPC;

	}

private:

	/**
	  * parameter for calculate the molar density of gas phase
      */
	const double RT=1.0/37.0;

	const double Hen = 7.65e-6; //Henry constant

	const double C_h = 1.530e-8;// Here I define a constant value C_h which is equal to Hen*M_G(Molar mass of gas)
	const double eps = 1e-7;

	const double S_gr = 0.0;//residual saturation of the gas phase
	const double S_lr = 0.4;//residual saturation of the liquid phase
	const double P_r = 2e+6;//pa

	double P_L;
	
	double X_L;
};


