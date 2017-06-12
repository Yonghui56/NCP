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
#include "AbstractEOS.h"
#include <cmath>


/**
*  EOS formulation for calculating each secondary variable
*  PC PG PL N_L N_G Kr_g Kr_l X_m X_M Sg, with only 3 unknows! 
*/
class EOS_H2_H2O_ISOT_3x_Complementary : public AbstractEOS
{
public:

	/**
	  * Constructor
	  */
	EOS_H2_H2O_ISOT_3x_Complementary() : AbstractEOS(3)
	{
	};

	/**
	  * DESTRUCTOR
	  */
	virtual ~EOS_H2_H2O_ISOT_3x_Complementary()
	{

	};

	/**
	  * realization of the eval function.
	  * this function will evaluate the resiudal of the governing equation,
	  * based on the unknown values given.
	  */
	virtual void eval(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalVector & res)
	{
		// getting unknowns
		double Sg = vec_unknowns(0);
		double X_m = vec_unknowns(1);
		double X_M = vec_unknowns(2);

		// calculating residual
		res(0) = Calc_mole_fraction_X_m(P_L, Sg, X_m);
		res(1) = Calc_mole_fraction_X_M(P_L, Sg, X_M);
		res(2) = Calc_Saturation(P_L, X_L, Sg, X_m, X_M);
	};

	/**
	  * realization of the calc_Jacobian function.
	  * this function will evaluate the Jacobian matrix,
	  * based on the unknown values given.
	  */
	virtual void calc_Jacobian(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalMatrix & J)
	{
		assert(J.cols() == this->get_n() && J.rows() == this->get_n());

		// getting unknowns
		double Sg = vec_unknowns(0);
		double X_m = vec_unknowns(1);
		double X_M = vec_unknowns(2);

		// evaluate J
		J.setZero();
		J(0, 0) = (Calc_mole_fraction_X_m(P_L, Sg + eps, X_m) - Calc_mole_fraction_X_m(P_L, Sg, X_m)) / eps;
		J(0, 1) = (Calc_mole_fraction_X_m(P_L, Sg, X_m + eps) - Calc_mole_fraction_X_m(P_L, Sg, X_m)) / eps;
		J(0, 2) = 0.0;
		J(1, 0) = (Calc_mole_fraction_X_M(P_L, Sg + eps, X_M) - Calc_mole_fraction_X_M(P_L, Sg, X_M)) / eps;
		J(1, 1) = 0.0;
		J(1, 2) = (Calc_mole_fraction_X_M(P_L, Sg, X_M + eps) - Calc_mole_fraction_X_M(P_L, Sg, X_M)) / eps;
		J(2, 0) = (Calc_Saturation(P_L, X_L, Sg + eps, X_m, X_M) - Calc_Saturation(P_L, X_L, Sg, X_m, X_M)) / eps;
		J(2, 1) = (Calc_Saturation(P_L, X_L, Sg, X_m + eps, X_M) - Calc_Saturation(P_L, X_L, Sg, X_m, X_M)) / eps;
		J(2, 2) = (Calc_Saturation(P_L, X_L, Sg, X_m, X_M + eps) - Calc_Saturation(P_L, X_L, Sg, X_m, X_M)) / eps;
		
		
		/*assert(J.cols() == this->get_n() && J.rows() == this->get_n());

		double PG(0.0);
		PG = getPG(P_L,Sg);
		double NG(0.0);
		NG = getN_G(P_L, Sg);
		double NL(0.0);
		NL = getN_L(P_L, Sg);
		double F1 = Sg;
		double G1 = getX_Comp1_Liq(P_L, Sg) - X_m;
		double F2 = 1 - Sg;
		double G2 = X_M - getX_Comp1_Gas(P_L, Sg);
		// evaluate J
		J.setZero();
		J(2, 0) = ((NG*(X_L - X_M) - NL*(X_L - X_m))*(Sg*NG + (1 - Sg)*NL) - (Sg*NG*(X_L - X_M) + (1 - Sg)*NL*(X_L-X_m))*(NG - NL)) / pow((Sg*NG + (1 - Sg)*NL), 2);
		J(2, 1) = -(1 - Sg)*NL / (Sg*NG + (1 - Sg)*NL);
		J(2, 2) = Sg*NG / (Sg*NG + (1 - Sg)*NL);

		if (F1 <= G1) {
			J(0, 0) = 1.0;
			J(0, 1) = 0.0;
			J(0, 2) = 0.0;
		}
		else{
			
			J(0, 0) = (getX_Comp1_Liq(P_L, Sg + eps) - getX_Comp1_Liq(P_L, Sg - eps)) / 2 / eps;// X_m_crit *P_crit*PC_0 / pow(PG + P_crit, 2);
			J(0, 1) = -1.0;
			J(0, 2) = 0.0;
		}
		if (F2 <= G2) {
			J(1, 0) = -1.0;
			J(1, 1) = 0.0;
			J(1, 2) = 0.0;
		}
		else{
			
			J(1, 0) = (getX_Comp1_Gas(P_L, Sg + eps) - getX_Comp1_Gas(P_L, Sg-eps))/2/eps;//; -0.5*X_M_crit*pow(P_crit, -1 / 2)*pow(PG, -1 / 2)*PC_0;
			J(1, 1) = 0.0;
			J(1, 2) = 1.0;
			

		}
		*/
	};
	/**
	* Complementary condition 1
	* for calculating molar fraction of light component in the liquid phase
	*/
	virtual double Calc_mole_fraction_X_m(double mean_P, double S, double X_m)
	{
		double PC_ = 0.0;
		double Res_X_m = 0.0;
		double PG_ = 0.0;
		double PL_ = 0.0;
		double gamma_ = 0.0;
		double X_m_P_S = 0.0;
		PC_ = getPCbySg_Mod(S);
		gamma_ = getweightedFunc(S);
		PG_ = mean_P + (1 - gamma_)*PC_;
		PL_ = mean_P - gamma_*PC_;
		X_m_P_S = getX_Comp1_Liq(mean_P, S);
		Res_X_m = std::min(S, X_m_P_S - X_m);
		return Res_X_m;
	}
	/**
	* Complementary condition 2
	* for calculating molar fraction of light component in the gas phase
	*/
	virtual double Calc_mole_fraction_X_M(double mean_P, double S, double X_M)
	{
		double PC_ = 0.0;
		double Res_X_M = 0.0;
		double PG_ = 0.0;
		double PL_ = 0.0;
		double gamma_ = 0.0;
		double X_M_P_S = 0.0;
		PC_ = getPCbySg_Mod(S);
		gamma_ = getweightedFunc(S);
		PG_ = mean_P + (1 - gamma_)*PC_;
		PL_ = mean_P - gamma_*PC_;
		X_M_P_S = getX_Comp1_Gas(mean_P, S);
		Res_X_M = std::min(1 - S, X_M - X_M_P_S);
		return Res_X_M;
	}
	/**
	* Complementary condition 3
	* for calculating the saturation
	*/
	virtual double Calc_Saturation(double mean_P, double mean_X, double S, double X_m, double X_M)
	{
		double PC_ = 0.0;
		double gamma_ = 0.0;
		double PG_ = 0.0;
		double PL_ = 0.0;
		double N_G_ = 0.0;
		double N_L_ = 0.0;
		double Res_S = 0.0;
		PC_ = getPCbySg_Mod(S);
		gamma_ = getweightedFunc(S);
		PG_ = mean_P + (1 - gamma_)*PC_;
		PL_ = mean_P - gamma_*PC_;
		N_G_ = getN_G(mean_P, S);
		N_L_ = getN_L(mean_P, S);
		Res_S = (S*N_G_*(mean_X - X_M) + (1 - S)*N_L_*(mean_X - X_m)) / (S*N_G_ + (1 - S)*N_L_);
		return Res_S;

	}
	virtual void calc_Rest_SecVar(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalVector & vec_rest_var)
	{
		vec_rest_var = ogsChem::LocalVector::Zero(5);
		double Sg = vec_unknowns(0);
		double X_m = vec_unknowns(1);
		double X_M = vec_unknowns(2);
		vec_rest_var(0) = getPCbySg_Mod(Sg);
		vec_rest_var(1) = getPL(P_L, Sg);
		vec_rest_var(2) = getPG(P_L, Sg);
		vec_rest_var(3) = getN_L(P_L, Sg);
		vec_rest_var(4) = getN_G(P_L, Sg);
	}

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
		double PC_0 = 20;
		double Pc = 0;
		Pc = PC_0*Sg;
		return Pc;
	}
	/**
	* calculate saturaion based on the value of capillary pressure 
	* explicit calculation
	*/
	virtual double getSat_byPC(double mean_P, double mean_X, double n = 1.49, double m = 1.0 - 1.0 / 1.49)
	{
		double P_c = 0.0;
		double Sg = 0.0;
		P_c = (mean_X / C_h) - mean_P;
		if (P_c >= 0){
			Sg = 1 - (((1 - S_gr - S_lr) / pow((pow((P_c / P_r), n) + 1), m)) + S_lr);
		}
		else
			Sg = 0.0;
		return Sg;

	}

	virtual double get_diriv_PC(double mean_P, double mean_X, double n = 1.49, double m = 1.0 - 1.0 / 1.49)
	{
		double dPC(0.0);
		double P_c = 0.0;
		P_c = (mean_X / C_h) - mean_P;
		dPC = m*n*(1 - S_gr - S_lr)*(1 / P_r)*(pow((P_c / P_r) , (n - 1)))*pow(((pow((P_c / P_r) , n)) + 1) ,(-m - 1));
		return dPC;

	}
	virtual double getS_bar(double Sg, double res_S_l = 0.10, double res_S_g = 0.0, double es = 0.59)
	{
		double S_bar = 0.0;

		S_bar = res_S_g + (1 - es)*(Sg - res_S_g) + 0.5*es*(1 - res_S_g - res_S_l);
		return S_bar;
	}

	/**
	  * based on primary variable mean_pressure to get the value of PG ---Gas phase pressure
	  */
	/*virtual double getPG(double mean_P, double Sg, double PC)
	{
		double PG = 0.0;
		//double Pc = getPcbySg(Sg);
		double gamma = getweightedFunc(Sg);
		PG = mean_P + (1 - gamma)*PC;
		return PG;
	}*/
	virtual double getPG(double mean_P, double Sg)//modified on 2,7,2014
	{
		double PG = 0.0;
		double PC = getPCbySg_Mod(Sg);
		double gamma = getweightedFunc(Sg);
		PG = mean_P + (1 - gamma)*PC;
		return PG;
	}

	/**
	  * based on primary variables to mean_pressure get the value of PL ---Liquid phase pressure
	  */
	/*virtual double getPL(double mean_P, double Sg, double PC)
	{
		double PL = 0.0;
		//double Pc = getPcbySg(Sg);
		double gamma = getweightedFunc(Sg);
		PL = mean_P - gamma*PC;
		return PL;
	}*/
	virtual double getPL(double mean_P, double Sg)//modified on 2,7,2014
	{
	double PL = 0.0;
	double PC = getPCbySg_Mod(Sg);
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
	* 
	* 
	*/
	virtual double get_deriv_weight_Func(double Sg)
	{
		double dgamma_dSg(0.0);
		dgamma_dSg = 2 - 2 * Sg;
		/*if (Sg >= 1)
			dgamma_dSg = 0;
		else if (Sg <= 0)
			dgamma_dSg = 0;*/
		return dgamma_dSg;
	}


	/**
	  * calculate the effective saturation of Liquid phase by Gas saturation
	  */
	virtual double getEffectSat_lbySg(double Sg, double res_S_l = 0.10, double res_S_g = 0.0)
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
	virtual double getKr_L_bySg(double Sg/*gas phase saturation*/, double m = 1.0 - 1.0 / 1.49, double res_S_l = 0.10, double res_S_g = 0.0)
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
	virtual double getKr_g_bySg(double Sg/*gas phase saturation*/, double m = 1.0 - 1.0 / 1.49, double res_S_l = 0.10, double res_S_g = 0.0)
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
	/*virtual double getX_Comp1_Liq(double mean_P, double Sg, double PG, double X_m_crit = 0.80, double P_crit = 2000.0)
	{
		double X_m = 0.0;
		//double PG = 0.0;
		double K = 0;

		//PG = getPG(mean_P, Sg);
		K = X_m_crit*((PG) / (PG + P_crit));
		X_m = std::max(0.0, K);/////? problem?
		return X_m;
	}*/
	/**
	* Calculate X_m
	*/
	virtual double getX_Comp1_Liq(double mean_P, double Sg,  double X_m_crit = 0.80, double P_crit = 2000.0)//modified on 2,7,2014
	{
		double X_m = 0.0;
		double PG = 0.0;
		double K = 0;

		PG = getPG(mean_P, Sg);
		K = X_m_crit*((PG) / (PG + P_crit));
		X_m = std::max(0.0, K);
		return X_m;
	}
	virtual double getX_Comp1_Gas(double mean_P, double Sg,  double X_M_crit = 0.85, double P_crit = 2000.0)//modified on 2,7,2014
	{
		double X_M = 0.0;
		double PG = 0.0;
		double K = 0;
		PG = getPG(mean_P, Sg);
		K = X_M_crit*sqrt(PG / P_crit);
		X_M = std::min(X_M_crit, K);
		return X_M;
	}

	virtual double getN_G(double mean_P, double Sg, double A = 39.7)//modified on 2,7,2014
	{
		double N_G = 0.0;
		double PG = 0.0;
		PG = getPG(mean_P, Sg);
		N_G = A*PG;
		return N_G;
	}
	/**
	  * intermedia function for the calculation of molar density in liquid phase
	  */
	virtual double getB_L(double mean_P, double PL, double s = 1.0e-4)
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
	virtual double getN_L(double mean_P, double Sg, double X_L_std = 0.9995, double N_L_std = 85000.0)// modified
	{
		double N_L = 0.0;
		double B_L = 0.0;
		double X_m = 0.0;
		double PL = 0.0;
		X_m = getX_Comp1_Liq(mean_P, Sg);
		PL = getPL(mean_P, Sg);
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
	
	/**
	* Calculate the derivatives using the analytical way 
	* Add on 2,7,2014
	*/
	/**
	* 
	*/
	virtual double Calculate_dSgdP(double P, double X, double S, double PL, double PG, double T = 303.0)
	{
		double dSgdP(0.0);
		double GAMMA = getweightedFunc(S);
		double BL = 1 / (1 + s*PL);
		double A = (1E+5) / R / T;
		dSgdP = ((0.5*A*PG*X_M_crit*S) / (P_crit*sqrt((P + (PC_0 - 1.*GAMMA*PC_0)*S) / P_crit)) + (N_L_std*P_crit* X_m_crit*(1. - 1.*X + (-1. + 1.*X)*S)) / (BL*pow((-1.*P_crit + P*(-1. + 1.*X_m_crit) + \
			PC_0*(-1. + GAMMA*(1. - 1.*X_m_crit) + 1.*X_m_crit)*S), 2))) / ((0.5*A*(-1. + GAMMA)*PC_0*PG*X_M_crit*S) / (P_crit*sqrt((P + (PC_0 - 1.*GAMMA*PC_0)*S) / \
			P_crit)) + ((-1. + GAMMA)*N_L_std* PC_0*(-1. + X_m_crit)*(-1. + S)*((P + P_crit)*X - 1.*P*X_m_crit - 1.*(-1. + GAMMA)*PC_0*(X - 1.*X_m_crit)*S)) / \
			(BL*pow((P + P_crit - 1.*P*X_m_crit + (-1. + GAMMA)*PC_0*(-1. + X_m_crit)*S), 2)) - \
			(1.*(-1. + GAMMA)*N_L_std* PC_0*(X - 1.*X_m_crit)*(-1. + S)) / (BL*(-1.*P_crit + \
			P*(-1. + X_m_crit) + PC_0*(-1. + GAMMA + X_m_crit - 1.*GAMMA*X_m_crit)*S)) + (N_L_std*((P + P_crit)*X - 1.*P*X_m_crit - 1.*(-1. + GAMMA)*PC_0*(X - 1.*X_m_crit)*S)) / \
			(BL*(-1.*P_crit + P*(-1. + X_m_crit) + PC_0*(-1. + GAMMA + X_m_crit - 1.*GAMMA*X_m_crit)*S)) + A*PG*(X - 1.*X_M_crit*sqrt((P + (PC_0 - 1.*GAMMA*PC_0)*S) / P_crit)));
		return dSgdP;

	}
	virtual double Calculate_dSgdX(double P, double X, double S, double PL, double PG,  double T=303.0)
	{
		double dSgdX(0.0);
		double GAMMA = getweightedFunc(S);
		double A = (1e+5) / R / T;
		double BL = 1 / (1 + s*PL);
		dSgdX=(N_L_std*(-1.*P - 1.*P_crit) + S*(N_L_std*(1.*P - 1.*PC_0 + 1.*GAMMA*PC_0 + 1.*P_crit) + A*BL*PG*(-1.*P - 1.*P_crit + 1.*P*X_m_crit) + PC_0*((1. - 1.*GAMMA)*N_L_std + A*BL*PG*(-1. + 1.*GAMMA + 1.*X_m_crit - 1.*GAMMA*X_m_crit))*\
			S)) / (BL*(1.*P_crit + P*(1. - 1.*X_m_crit) + PC_0*(1. - 1.*X_m_crit + GAMMA*(-1. + 1.*X_m_crit))*S)*((0.5*A*(-1. + GAMMA)*PC_0*PG*X_M_crit*S) / (P_crit*sqrt((P + (PC_0 - 1.*GAMMA*PC_0)*S) / \
			P_crit)) + ((-1. + GAMMA)*N_L_std* PC_0*(-1. + X_m_crit)*(-1. + S)*((P + P_crit)*X - 1.*P*X_m_crit - 1.*(-1. + GAMMA)*PC_0*(X - 1.*X_m_crit)*S)) / (BL*pow((P + P_crit - 1.*P*X_m_crit + (-1. + GAMMA)*PC_0*(-1. + X_m_crit)* S), 2)) - \
			(1.*(-1. + GAMMA)*N_L_std* PC_0*(X - 1.*X_m_crit)*(-1. + S)) / (BL*(-1.*P_crit + P*(-1. + X_m_crit) + PC_0*(-1. + GAMMA + X_m_crit - 1.*GAMMA*X_m_crit)*S)) + \
			(N_L_std*((P + P_crit)*X - 1.*P*X_m_crit - 1.*(-1. + GAMMA)*PC_0*(X - 1.*X_m_crit)*S)) / (BL*(-1.*P_crit + P*(-1. + X_m_crit) + PC_0*(-1. + GAMMA + X_m_crit - 1.*GAMMA*X_m_crit)*S)) + A*PG*(X - 1.*X_M_crit*sqrt((P + (PC_0 - 1.*GAMMA*PC_0)*S) / P_crit))));
		return dSgdX;


	}
	/*
	* return the value of dPCdP(derivative of capillary pressure based on pressure)
	*/
	virtual double Calculate_dPCdP(double dSgdP)
	{
		double dPCdP(0.0);
		dPCdP = PC_0*dSgdP;
		return dPCdP;
	}
	/*
	* return the value of dPCdX(derivative of capillary pressure based on X)
	*/
	virtual double Calculate_dPCdX(double dSgdX)
	{
		double dPCdX(0.0);
		dPCdX = PC_0*dSgdX;
		return dPCdX;
	}
	/*
	* return derivative of dPGdX
	*/
	virtual double Calculate_dPGdX(double S, double PC, double dSgdX,  double dPCdX)
	{
		double dPGdX(0.0);
		double GAMMA(0.0);
		double dGammadS(0.0);
		dGammadS = get_deriv_weight_Func(S);
		GAMMA = getweightedFunc(S);
		dPGdX = - PC*dGammadS*dSgdX + (1 - GAMMA)*dPCdX;
		return dPGdX;
	}

	virtual double Calculate_dPGdP(double S, double PC, double dSgdP,  double dPCdP)
	{
		double dPGdP(0.0);
		double GAMMA(0.0);
		double dGammadS(0.0);
		GAMMA = getweightedFunc(S);
		dGammadS = get_deriv_weight_Func(S);
		
		dPGdP = 1 - (PC*dGammadS*dSgdP) + (1 - GAMMA)*dPCdP;
		return dPGdP;
	}

	/*
	*return derivative of dPL
	*/
	virtual double Calculate_dPLdP(double S, double PC, double dSgdP, double dPCdP)
	{
		double dPLdP(0.0);
		double GAMMA(0.0);
		double dGammadS(0.0);
		GAMMA = getweightedFunc(S);
		dGammadS = get_deriv_weight_Func(S);
		dPLdP = 1 - (GAMMA*dPCdP + PC*dGammadS*dSgdP);
		return dPLdP;
	}

	virtual double Calculate_dPLdX(double S, double PC, double dSgdX, double dPCdX)
	{
		double dPLdX(0.0);
		double GAMMA(0.0);
		double dGammadS(0.0);
		GAMMA = getweightedFunc(S);
		dGammadS = get_deriv_weight_Func(S);

		dPLdX = - (PC*dGammadS*dSgdX + GAMMA*dPCdX);
		return dPLdX;
	}
	/*
	* return the value of dX_MdP(derivative of light component molar fracrion based on pressure )
	*/
	virtual double Calculate_dXMdP(double PG, double dPGdP)
	{
		double dX_MdP(0.0);
		
		dX_MdP = 0.5*X_M_crit*(1 / sqrt(P_crit))*(1 / sqrt(PG))*dPGdP;

		return dX_MdP;
	}
	/*
	* return the value of dX_MdX(derivative of light component molar fracrion based on total molar fraction )
	*/

	virtual double Calculate_dXMdX(double PG, double dPGdX)
	{
		double dX_MdX(0.0);
		dX_MdX = 0.5*X_M_crit*(1 / sqrt(P_crit))*(1 / sqrt(PG))*dPGdX;
		return dX_MdX;
	}
	virtual double Calculate_dX_mdX(double PG, double dPGdX)
	{
		double X_mdX(0.0);

		X_mdX = (X_m_crit*P_crit / pow((PG + P_crit), 2))*dPGdX;
		
		return X_mdX;
	}
	virtual double Calculate_dX_mdP(double PG, double dPGdP)
	{
		double X_m_dP(0.0);

		X_m_dP = (X_m_crit*P_crit / pow((PG + P_crit), 2))*dPGdP;
		return X_m_dP;
	}

	virtual double Calculate_dN_GdP(double dPGdP,double T=303.0)
	{
		double dN_GdP(0.0);
		double RT;
		RT = (1e+5) / R / T;
		dN_GdP = RT*dPGdP;
		return dN_GdP;
	}

	virtual double Calculate_dN_GdX(double dPGdX, double T = 303.0)
	{

		double dN_GdX(0.0);
		double RT;
		RT = (1e+5) / R / T;
		dN_GdX = RT*dPGdX;
		return dN_GdX;
	}
	virtual double Calculate_dN_LdP(double X_m, double PL, double dX_mdP, double dPLdP)
	{
		double dN_LdP(0.0);
		double dBLdP(0.0);
		double BL(0.0);
			BL = 1 / (1 + s*PL);
		dBLdP = -s*dPLdP / pow((1 + s*PL), 2);
		dN_LdP=N_L_std*((-dBLdP / (1 - X_m) / BL / BL) + dX_mdP / BL / pow((1 - X_m) , 2));
		return dN_LdP;
	}
	virtual double Calculate_dN_LdX(double X_m, double PL, double dX_mdX,  double dPLdX)
	{
		double dN_LdX(0.0);
		double dBLdX(0.0);
		double BL(0.0);
		BL = 1 / (1 + s*PL);
		dBLdX = -s*dPLdX / pow((1 + s*PL), 2);
		dN_LdX=	N_L_std*((-dBLdX / (1 - X_m) / BL / BL) + dX_mdX / BL / pow((1 - X_m) , 2));
		return dN_LdX;
	}

	

	/**
	* Calculate the gradient of 
	*/

	/**
	* Define three character function for calculate the gradient
	*/
	std::size_t Characteristic_Function_xi(double Sg)
	{
		size_t xi(0);
		if (std::abs(Sg - 1) < 1e-10)
			xi = 0;// here the gas saturation is close to 1 and it is gas phase
		else if (abs(Sg) < 1e-10)
			xi = 0;//here the 
		else
			xi = 1;
		return xi;
	}
	std::size_t Characteristic_Function_alpha(double Sg)
	{
		size_t alpha(0);
		if (std::abs(Sg) < 1e-5)
			alpha = 1;
		else
			alpha = 0;
		return alpha;
	}
	std::size_t Characteristic_Function_Beta(double Sg)
	{
		size_t Beta(0);
		if (std::abs(Sg - 1) < 1e-5)
			Beta = 1;
		else
			Beta = 0;
		return Beta;
	}

private:
	/* some const parameters*/
	const double eps = 1e-7;
	const double R = 8.314;
	const double X_m_crit = 0.8;
	const double X_M_crit = 0.85;
	const double PC_0 = 20.0;
	const double P_crit = 100 * PC_0;
	const double N_L_std = 85000.0;
	const double s = 1e-4;
	const double C_h = 1.53e-8;
	const double S_gr = 0.0;
	const double S_lr = 0.4;
	const double P_r = 2e+6;
	double P_L;

	double X_L;
};


