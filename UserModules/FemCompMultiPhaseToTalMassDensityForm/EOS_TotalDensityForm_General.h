/**
* Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.com)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.com/LICENSE.txt
*
*
* \file EOS_H2_H2O_ISOT_3x_Complementary.h
*
* Created on 2014-05-30 by Haibing SHAO
*/

#pragma once


#include "math.h"
#include <algorithm>
#include "AbstractEOS_PressureForm.h"
#include <cmath>


/**
*  EOS formulation for calculating each secondary variable
*  PC PG PL N_L N_G Kr_g Kr_l X_m X_M Sg, with only 3 unknows! 
*/
class EOS_PressureForm_General : public AbstractEOS_PressureForm
{
public:

	/**
	  * Constructor
	  */
	EOS_PressureForm_General() : AbstractEOS_PressureForm(1)
	{
	};

	/**
	  * DESTRUCTOR
	  */
	virtual ~EOS_PressureForm_General()
	{

	};

	/**
	  * realization of the eval function.
	  * this function will evaluate the resiudal of the governing equation,
	  * based on the unknown values given.
	  */
	virtual void set_env_condition(ogsChem::LocalVector & env_conditon)
	{
		// mean pressure
		P_L = env_conditon(0);
		// molar fraction of the light component
		X_L = env_conditon(1);
	};

	virtual double getSat_byPC(double n = 1.49, double m = 1.0 - 1.0 / 1.49)//for MoMaS benchmark2 n = 1.49
	{
		double P_c = 0.0;
		double Sg = 0.0;
		double PG = 0.0;
		PG = X_L / C_h + P_vapor*(rho_l_std / (rho_l_std + (M_L / M_G)*X_L));
		P_c = PG - P_L;
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
		double PG = 0.0;
		PG = X_L / C_h + P_vapor*(rho_l_std / (rho_l_std + (M_L / M_G)*X_L));
		P_c = PG - P_L;
		dPC = m*n*(1 - S_gr - S_lr)*(1 / P_r)*(pow((P_c / P_r), (n - 1)))*pow(((pow((P_c / P_r), n)) + 1), (-m - 1));
		if (P_c <0)
		{
			dPC = 0.0;//just for test
		}
		return dPC;

	}

private:
	const double S_gr = 0.0;
	const double S_lr = 0.4;
	const double P_r = 2e+6;
	 const double C_h = 1.530e-8;// Here I define a constant value C_h which is equal to Hen*M_G(Molar mass of gas)
	 const double P_vapor =1.0 ;
	 const double rho_l_std = 1000.0;
	 const double M_G;
	 const double M_L;
		
	double P_L;

	double X_L;
};


