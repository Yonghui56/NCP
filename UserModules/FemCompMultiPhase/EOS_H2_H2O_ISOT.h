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
#include "AbstractEOS.h"


/**
 *  EOS formulation for calculating each secondary variable
 *  PC PG PL N_L N_G Kr_g Kr_l X_m X_M Sg
 */
class EOS_H2_H2O_ISOT : public AbstractEOS
{
public:
    
	/**
	  * Constructor
	  */
	EOS_H2_H2O_ISOT() : AbstractEOS(8), RT(1.0 / 39.7)
    {
    };

	/**
	  * DESTRUCTOR
	  */
	virtual ~EOS_H2_H2O_ISOT()
    {

    };

	/**
	  * realization of the eval function. 
	  * this function will evaluate the resiudal of the governing equation, 
	  * based on the unknown values given. 
	  */
	virtual void eval(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalVector & res)
	{
		double PC = vec_unknowns(0);
		double PL = vec_unknowns(1);
		double PG = vec_unknowns(2);
		double X_m = vec_unknowns(3);
		double X_M = vec_unknowns(4);
		double N_G = vec_unknowns(5);
		double N_L = vec_unknowns(6);
		double Sg = vec_unknowns(7);
		//Res = ogsChem::LocalVector::Zero(N);
		res(0) = PC - this->getPCbySg_Mod(Sg);//getPcbySg(Sg);
		res(1) = PL - this->getPL(P_L, Sg, PC);
		res(2) = PG - this->getPG(P_L, Sg, PC);
		res(3) = X_m - this->getX_Comp1_Liq(P_L, Sg, PG);
		res(4) = X_M - this->getX_Comp1_Gas(P_L, Sg, PG);
		res(5) = N_G - this->getN_G(P_L, Sg, PG);
		res(6) = N_L - this->getN_L(P_L, PL, X_m);
		res(7) = Sg - this->getSg(X_L, N_G, N_L, X_m, X_M);

	};

	/**
	  * realization of the calc_Jacobian function.
	  * this function will evaluate the Jacobian matrix,
	  * based on the unknown values given.
	  */
	virtual void calc_Jacobian(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalMatrix & J)
	{
		assert(J.cols() == this->get_n() && J.rows() == this->get_n() );

		double PC = vec_unknowns(0);
		double PL = vec_unknowns(1);
		double PG = vec_unknowns(2);
		double X_m = vec_unknowns(3);
		double X_M = vec_unknowns(4);
		double N_G = vec_unknowns(5);
		double N_L = vec_unknowns(6);
		double Sg = vec_unknowns(7);

		// Matrix_Jacobian
		//Matrix_Jacobian = ogsChem::LocalMatrix::Zero(N, N);

		J.setZero(); 
		J(0, 7) = -(this->getPCbySg_Mod(Sg + eps) - this->getPCbySg_Mod(Sg)) / eps;

		J(1, 7) = -(this->getPL(P_L, Sg + eps, PC) - this->getPL(P_L, Sg, PC)) / eps;
		J(1, 0) = this->getweightedFunc(Sg);

		J(2, 0) = -(1 - this->getweightedFunc(Sg));
		J(2, 7) = -(this->getPG(P_L, Sg + eps, PC) - this->getPG(P_L, Sg, PC)) / eps;

		J(3, 2) = -(this->getX_Comp1_Liq(P_L, Sg, PG + eps) - this->getX_Comp1_Liq(P_L, Sg, PG)) / eps;
		J(4, 2) = -(this->getX_Comp1_Gas(P_L, Sg, PG + eps) - this->getX_Comp1_Gas(P_L, Sg, PG)) / eps;
		J(5, 2) = -1 / RT;

		J(6, 1) = -(this->getN_L(P_L, PL + eps, X_m) - this->getN_L(P_L, PL, X_m)) / eps;
		J(6, 3) = -(this->getN_L(P_L, PL, X_m + eps) - this->getN_L(P_L, PL, X_m)) / eps;

		J(7, 3) = -(this->getSg(X_L, N_G, N_L, X_m + eps, X_M) - this->getSg(X_L, N_G, N_L, X_m, X_M)) / eps;
		J(7, 4) = -(this->getSg(X_L, N_G, N_L, X_m, X_M + eps) - this->getSg(X_L, N_G, N_L, X_m, X_M)) / eps;
		J(7, 5) = -(this->getSg(X_L, N_G + eps, N_L, X_m, X_M) - this->getSg(X_L, N_G, N_L, X_m, X_M)) / eps;
		J(7, 6) = -(this->getSg(X_L, N_G, N_L + eps, X_m, X_M) - this->getSg(X_L, N_G, N_L, X_m, X_M)) / eps;
	};

	virtual void set_env_condition(ogsChem::LocalVector & env_conditon)
	{
		// mean pressure
		P_L = env_conditon(0);
		// molar fraction of the light component
		X_L = env_conditon(1);
	}; 

	/**
	  * return the capillary pressure by saturation
	  */
	virtual double getPcbySg(double Sg, double PC_0 = 20, double n = 1.49, double m = 1.0 - 1.0 / 1.49)
	{
		double Pc = 0.0;
		double EffectSat_l_1 = 0.0;

		EffectSat_l_1 = getEffectSat_lbySg(Sg);
		// PC = PC_0*(Se_l ^ (-1 / m) - 1) ^ (1 / n);
		Pc = PC_0*pow((pow(EffectSat_l_1, -1 / m) - 1), 1 / n);
		return Pc;
	}
	virtual double getPCbySg_Mod(double Sg, double res_S_l = 0.10, double res_S_g = 0.0)
	{
		double S_bar = 0.0;
		double Pc = 0.0;
		double S_lr_bar = 0.0;
		double S_lr_eps_bar = 0.0;
		if ( Sg<1-res_S_l && Sg>res_S_g)
		{
			S_bar = getS_bar(Sg);
			Pc = getPcbySg(S_bar);

		}
		else if (Sg == res_S_g)
		{
			Pc = 0.0;
		}
		else if (Sg < res_S_g)
		{
			S_bar = getS_bar(0.0);
			Pc = getPcbySg(S_bar);
		}
		else if (Sg>res_S_l && Sg<=1)
		{
			S_lr_bar = getS_bar(1 - res_S_l);
			S_lr_eps_bar = getS_bar(1 - res_S_l+eps);
			Pc = getPcbySg(S_lr_bar) + ((getPcbySg(S_lr_eps_bar) - getPcbySg(S_lr_bar)) / eps) *(Sg - 1 + res_S_l);
		}
		else if (Sg > 1)
		{
			S_lr_bar = getS_bar(1 - res_S_l);
			S_lr_eps_bar = getS_bar(1 - res_S_l + eps);
			Pc = getPcbySg(S_lr_bar) + ((getPcbySg(S_lr_eps_bar) - getPcbySg(S_lr_bar)) / eps) *(res_S_l);
		}
		return Pc;
	}
	virtual double getS_bar(double Sg, double res_S_l = 0.10, double res_S_g = 0.0, double es=0.59)
	{
		double S_bar = 0.0;
		
		S_bar = res_S_g + (1 - es)*(Sg - res_S_g) + 0.5*es*(1 - res_S_g - res_S_l);
		return S_bar;
	}

	/**
	  * based on primary variable mean_pressure to get the value of PG ---Gas phase pressure
	  */
	virtual double getPG(double mean_P, double Sg, double PC)
	{
		double PG = 0.0;
		//double Pc = getPcbySg(Sg);
		double gamma = getweightedFunc(Sg);
		PG = mean_P + (1 - gamma)*PC;
		return PG;
	}

	/**
	  * based on primary variables to mean_pressure get the value of PL ---Liquid phase pressure
	  */ 
	virtual double getPL(double mean_P, double Sg, double PC)
	{
		double PL = 0.0;
		//double Pc = getPcbySg(Sg);
		double gamma = getweightedFunc(Sg);
		PL = mean_P - gamma*PC;
		return PL;
	}

	/** 
	  * weighted function for mean pressure based on  gas/liquid pressure function of Saturation 
	  * extend S to R because of complemntaty condition
	  */ 
	virtual double getweightedFunc(double Sg)
	{
		double gamma;
		if (Sg >= 1)
			gamma = 1;
		else if (Sg <= 0)
			gamma = 0;
		else
			gamma = Sg*(2 - Sg);
		return gamma; 
	}


	/**
	  * calculate the effective saturation of Liquid phase by Gas saturation
	  */
	virtual double getEffectSat_lbySg(double Sg, double res_S_l=0.10, double res_S_g=0.0)
	{
		double EffectSat_l = 0.0;
		// this->res_saturation_g->eval(Sg, res_S_g);

		if (Sg < res_S_g)
			EffectSat_l = (1 - Sg - res_S_l) / (1 - res_S_g - res_S_l);
		else if (Sg>(1 - res_S_l))
			EffectSat_l = 0.0;
		else
			EffectSat_l = (1 - Sg - res_S_l) / (1 - res_S_g - res_S_l);

		return EffectSat_l;
	}

	/**
	  * calculate the effective saturation of gas phase by Gas saturation
	  */
	virtual double getEffectSat_gbySg(double Sg, double res_S_l, double res_S_g)
	{
		double EffectSat_g = 0.0;
		// this->res_saturation_g->eval(Sg, res_S_g);

		EffectSat_g = (Sg - res_S_g) / (1 - res_S_g - res_S_l);

		return EffectSat_g;
	}

	/**
	  * get the relative permeability based on Sg
	  */
	virtual double getKr_L_bySg(double Sg/*gas phase pressure*/, double m, double res_S_l, double res_S_g)
	{
		double Kr_L = 0.0;
		double EffectSat_l = 0.0;
		//double EffectSat_g = 0.0;
		EffectSat_l = getEffectSat_lbySg(Sg, 0.4, 0.0);
		//EffectSat_g = getEffectSat_gbySg(Sg, 0.4, 0.0);
		Kr_L = sqrt(EffectSat_l)*pow(1 - pow(1 - pow(EffectSat_l, 1 / m), m), 2);
		return Kr_L;
	}

	/**
	  * get the relative permeability of gas phase based on Sg
	  */
	virtual double getKr_g_bySg(double Sg/*gas phase pressure*/, double m, double res_S_l, double res_S_g)
	{
		double Kr_G = 0.0;
		double EffectSat_g = 0.0;

		EffectSat_g = getEffectSat_gbySg(Sg, res_S_l, res_S_g);
		Kr_G = sqrt(EffectSat_g)*pow(1 - pow(1 - EffectSat_g, 1 / m), 2 * m);

		return Kr_G;
	}

	/**
	  * get the molar fraction of the light component in liquid phase
	  */
	virtual double getX_Comp1_Liq(double mean_P, double Sg, double PG, double X_m_crit=0.80, double P_crit=2000.0)
	{
		double X_m = 0.0;
		//double PG = 0.0;
		double K = 0;

		//PG = getPG(mean_P, Sg);
		K = X_m_crit*((PG) / (PG + P_crit));
		X_m = std::max(0.0, K);/////? problem?
		return X_m;
	}
	
	/**
	  *   get the molar fraction of the light component in gas phase
	  */
	virtual double getX_Comp1_Gas(double mean_P, double Sg, double PG, double X_M_crit=0.85, double P_crit=2000.0)
	{
		double X_M = 0.0;
		//double PG = 0.0;
		double K = 0;
		K = X_M_crit*sqrt(PG / P_crit);
		X_M = std::min(X_M_crit, K);
		return X_M;
	}

	/**
	  * Calculate the molar density in gas phase and liquid phase, respectively
	  */
	virtual double getN_G(double mean_P, double Sg, double PG, double A=39.7)
	{
		double N_G = 0.0;
		//double PG = 0.0;
		//PG = getPG(mean_P, Sg);
		N_G = A*PG;
		return N_G;
	}

	/**
	  * intermedia function for the calculation of molar density in liquid phase
	  */
	virtual double getB_L(double mean_P, double PL, double s=1.0e-4)
	{
		double B_L = 0.0;
		//double PL = 0.0;
		//PL = getPL(mean_P, Sg);
		B_L = std::min(1.0, 1.0 / (1.0 + s*PL));
		return B_L;
	}

	/**
	  * get the molar density of liquid phase 
	  */
	virtual double getN_L(double mean_P, double PL, double X_m, double X_L_std =0.9995, double N_L_std=85000.0, double s=1e-4)
	{
		double N_L = 0.0;
		double B_L = 0.0;

		B_L = getB_L(mean_P, PL);
		double K1 = 0.0;
		double K2 = 0.0;
		K1 = N_L_std / B_L / (1 - X_L_std);
		K2 = N_L_std / B_L / (1 - X_m);
		N_L = std::min(K1, K2);
		return N_L;
	}

	/**
	  * get Saturation of gas phase based on molar fraction-X, X_m, X_M
	  */
	virtual double getSg(double Molar_Fraction, double N_G, double N_L, double X_m, double X_M)
	{
		double Sg = 0.0;

		Sg = (N_L)*(Molar_Fraction - X_m) / (N_G*(X_M - Molar_Fraction) + N_L*(Molar_Fraction - X_m));
		return Sg;
	}

private:

	/**
	  * parameter for calculate the molar density of gas phase
      */
	const double RT;

	const double eps = 1e-6;

	double P_L;
	
	double X_L;
};


