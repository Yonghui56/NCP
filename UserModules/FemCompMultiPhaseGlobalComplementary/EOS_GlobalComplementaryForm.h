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
#include "AbstractEOS_GlobalComplementaryForm.h"
#include <cmath>

/**
 *  EOS formulation for calculating each secondary variable
 *  PC PG PL N_L N_G Kr_g Kr_l X_m X_M Sg
 */
class EOS_GlobalComplementaryForm : public AbstractEOS_GlobalComplementaryForm
{
public:
    
	/**
	  * Constructor
	  */
	EOS_GlobalComplementaryForm() : AbstractEOS_GlobalComplementaryForm(3)
    {
    };

	/**
	  * DESTRUCTOR
	  */
	virtual ~EOS_GlobalComplementaryForm()
    {

    };

	


	virtual void set_env_condition(ogsChem::LocalVector & env_condition)
	{
		// mean pressure
		P_L = env_condition(0);
		// molar fraction of the light component
		X_L = env_condition(1);
		//SATURATION of the gas phase
		S_G = env_condition(2);

	}; 
	/**
	* realization of the calc the residual of saturation
	* EOS--Function 1
	* X=Sl*rho_L_h+Sg*rho_G_h
	*/
	

	virtual double getS_bar(double Sg, double xi=1e-5)
	{
		double S_bar = 0.0;
		
		S_bar = S_gr + (1 - xi)*(Sg - S_gr) + 0.5*xi*(1 - S_gr - S_lr);
		return S_bar;
	}
	/**
	*  van Genuchten capillary pressure-saturation Model
	*/
	virtual double getPc_vG_Sg(double Sg)
	{
		double Pc_vG = 0.0;
		double S_le = 0.0;
	    //effective saturation
		S_le =  (1 - Sg - S_lr) / (1 - S_gr - S_lr);
		//Pc_vG = P_r*(S_le. ^ (-1 / m) - 1). ^ (1 / n);
		Pc_vG = P_r*pow(pow(S_le, -1 / m) - 1, 1 / n);
		return Pc_vG;
	}
	/**
	* regularized van Genuchten capillary pressure-saturation Model
	*/
	virtual double getPc_bar_vG_Sg(double Sg, double xi=1e-5)
	{
		double Pc_bar_vG = 0.0;
		double S_bar;
		S_bar = getS_bar(Sg);
		Pc_bar_vG = getPc_vG_Sg(S_bar) - getPc_vG_Sg(S_gr + (1 - S_gr - S_lr)*xi / 2);
		return Pc_bar_vG;
	}
	/**
	* derivative dPCdS based on standard van Genuchten capillary pressure-saturation Model
	*/
	virtual double get_dPCdS_vG(double Sg)
	{
		double dPcdSg = 0.0;
		double S_le = (1 - Sg - S_lr) / (1 - S_gr - S_lr);
		dPcdSg = P_r*(1 / (m*n))*(1 / (1 - S_lr - S_gr))*pow(pow(S_le, (-1 / m)) - 1, (1 / n) - 1)*pow(S_le, (-1 / m)) / S_le;
		return dPcdSg;
	}
	/**
	* derivative dPCdS based on regularized van Genuchten capillary pressure-saturation Model
	*/
	virtual double get_dPCdS_vG_bar(double Sg, double xi=1e-5)
	{
		double dPCdS(0.0);
		double S_bar = 0.0;
		S_bar = getS_bar(Sg);
		dPCdS = get_dPCdS_vG(S_bar)*(1 - xi);
		return dPCdS;
	}

	/**
	* Derivation of PC (capillary pressure) in terms of saturation based on van Genuchten model
	*/

	virtual double getPcbySg(double Sg)//for waste package case1 n=2.0
	{
		double Pc = 0.0;
		if (Sg >= S_gr && Sg <= 1 - S_lr)
		{
			Pc = getPc_bar_vG_Sg(Sg);
		}
		else if(Sg < S_gr)
		{
			Pc = getPc_bar_vG_Sg(S_gr) + get_dPCdS_vG_bar(S_gr)*(Sg - S_gr);
		}
		else
		{
			Pc = getPc_bar_vG_Sg(1 - S_lr) + get_dPCdS_vG_bar(1 - S_lr)*(Sg - 1 + S_lr);
		}
		/*double Pc = 0.0;
		double EffectSat_l_1 = 0.0;

		EffectSat_l_1 = getEffectSat_lbySg(Sg);
		// PC = PC_0*(Se_l ^ (-1 / m) - 1) ^ (1 / n);
		Pc = P_r*pow((pow(EffectSat_l_1, -1 / m) - 1), 1 / n);
		if (Sg < 0.0)
			Pc = 0.0;*/
		return Pc;
	}
	/**
	* Derivation of PC (capillary pressure) in terms of saturation based on van Genuchten model
	*/
	virtual double Deriv_dPCdS(double Sg)
	{
		double dPcdSg(0.0);
		if (Sg >= S_gr && Sg <= 1 - S_lr)
		{
			dPcdSg = get_dPCdS_vG_bar(Sg);
		}
		else if (Sg < S_gr)
		{
			dPcdSg = get_dPCdS_vG_bar(S_gr);
		}
		else
		{
			dPcdSg = get_dPCdS_vG_bar(1 - S_lr);
		}
		/*double S_le(0.0);
		S_le = getEffectSat_lbySg(Sg);
		double dPcdSg(0.0);

		dPcdSg = P_r*(1 / (m*n))*(1 / (1 - S_lr - S_gr))*pow(pow(S_le, (-1 / m)) - 1, (1 / n) - 1)*pow(S_le, (-1 / m)) / S_le;

		if (Sg<=0)
		{
			dPcdSg = 0.0;
		}*/
		return dPcdSg;
	}
	/**
	* calculate the effective saturation of Liquid phase by Gas saturation
	*/
	virtual double getEffectSat_lbySg(double Sg)
	{
		double EffectSat_l = 0.0;
		// this->res_saturation_g->eval(Sg, res_S_g);

		/*if (Sg < S_gr)
			EffectSat_l = (1 - Sg - S_lr) / (1 - S_gr - S_lr);
		else if (Sg>(1 - S_lr))
			EffectSat_l = 0.0;
		else*/
		EffectSat_l = (1 - Sg - S_lr) / (1 - S_gr - S_lr);

		return EffectSat_l;
	}

	
	/**
	* based on primary variable mean_pressure to get the value of PG ---Gas phase pressure
	*/
	virtual double getPG(double Sg)
	{
		double PG = 0.0;
		double PC = 0.0;
		PC= getPcbySg(Sg);
		PG = P_L + PC;
		return PG;
	}


	/**
	* based on Sg  get the value of dissoved hydrogen mass density in the LIQUID phase
	*/
	virtual double get_RHO_L_H(double sg)
	{
		double rho_L_h(0.0);
		double PG(0.0);
		PG = getPG(sg);

		//rho_L_h = std::min(C_h*PG, X_L*M_G*rho_L_std / M_L);//Henry law
		rho_L_h = X_L*M_G*rho_L_std / M_L;
		return rho_L_h;
	}
	/**
	* based on Sg  get the value of hydrogen mass density in the GAS phase
	*/
	virtual double get_RHO_G_H(double sg)
	{
		double rho_G_h(0.0);
		double PG(0.0);
		PG = getPG(sg);
		rho_G_h =C_v*PG ;//ideal gas law
		return rho_G_h;
	}


	/**
	* water vapor pressure of pure water
	* Function of Temperature
	* comes from Amaziane's paper
	*/
	virtual double get_Vapor_Pressure(double T)
	{
		// Here unit of T is Celsius;
		double P_vapor(0.0);
		double tmp(0.0);
		tmp = 2.786 + 0.031514*T - 1.2373e-4*pow(T, 2) + 4.2267e-7*pow(T, 3) - 8.1308e-10*pow(T, 4);
		P_vapor = pow(10, tmp);
		return P_vapor;
	}

private:

	/**
	  * parameter for calculate the molar density of gas phase
      */
	const double RT=1.0/37.0;
	const double T0 = 303;// [K]
	const double Hen = 7.65e-6; //Henry constant
	//const double P_vapor = 0.0;
	const double C_h = 1.530e-8;// Here I define a constant value C_h which is equal to Hen*M_G(Molar mass of gas) Hen* M_G
	const double C_v = 7.938638137741564e-07;// 7.939211048841233e-07;//M_G/RT
	const double C_w = 3.969319068870782e-06;//M_L/RT
	const double eps = 1e-12;
	const double rho_L_std = 1000;
	const double M_G = 0.002;
	const double M_L = 0.01;
	/**
	* parameters for van-Genuchten capillary pressure-saturation model
	*/
	const double S_gr = 0.0;//residual saturation of the gas phase
	const double S_lr = 0.4;//residual saturation of the liquid phase
	const double P_r = 2e+6;//pa
	const double n = 1.49;
	const double m = 1 - 1 / n;

	double P_L;
	
	double X_L;
	double S_G;
	double _alpha;
	double _beta;
};


