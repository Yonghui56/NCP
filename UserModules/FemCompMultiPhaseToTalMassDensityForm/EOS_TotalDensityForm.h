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
#include "AbstractEOS_TotalDensityForm.h"
#include <cmath>

/**
 *  EOS formulation for calculating each secondary variable
 *  PC PG PL N_L N_G Kr_g Kr_l X_m X_M Sg
 */
class EOS_TotalDensityForm : public AbstractEOS_TotalDensityForm
{
public:
    
	/**
	  * Constructor
	  */
	EOS_TotalDensityForm() : AbstractEOS_TotalDensityForm(3)
    {
    };

	/**
	  * DESTRUCTOR
	  */
	virtual ~EOS_TotalDensityForm()
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
		double rho_L_h = vec_unknowns(1);
		double rho_G_h = vec_unknowns(2);

		// calculating residual
		res(0) = Calc_Res_Sg(Sg, rho_L_h, rho_G_h);
		res(1) = Calc_Res_rho_L_h(Sg, rho_L_h);
		res(2) = Calc_Res_rho_G_h(Sg, rho_G_h);
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
		double rho_L_h = vec_unknowns(1);
		double rho_G_h = vec_unknowns(2);
		double PG(0.0);
		PG = getPG(Sg);
		double PG_h = Function_PG_h(PG);
		double F1 = Sg;
		double G1 = std::min(C_h*PG_h, X_L) - rho_L_h;
		double F2 = 1 - Sg;
		double G2 = rho_G_h - std::max(C_v*PG_h, X_L);
		// evaluate J
		J.setZero();
		J(0, 0) = rho_L_h - rho_G_h;//-----dF(1)/dSg
		J(0, 1) = -(1 - Sg);
		J(0, 2) = -Sg;

		if (F1 <= G1) {
			J(1, 0) = 1.0;
			J(1, 1) = 0.0;
			J(1, 2) = 0.0;
		}
		else{
			if (C_h*PG_h <= X_L)
			{
				J(1, 0) = C_h*Deriv_dPGH_dPG(Sg)*Deriv_dPGdSg(Sg);
				J(1, 1) = -1.0;
				J(1, 2) = 0.0;
			}
			else{
				J(1, 0) = 0.0;
				J(1, 1) = -1.0;
				J(1, 2) = 0.0;
			}

		}
		if (F2 <= G2) {
			J(2, 0) = -1.0;
			J(2, 1) = 0.0;
			J(2, 2) = 0.0;
		}
		else{
			if (C_v*PG_h >= X_L)
			{
				J(2, 0) = -C_v*Deriv_dPGH_dPG(Sg)*Deriv_dPGdSg(Sg);
				J(2, 1) = 0.0;
				J(2, 2) = 1.0;
			}
			else{
				J(2, 0) = 0.0;
				J(2, 1) = 0.0;
				J(2, 2) = 1.0;
			}

		}
	};
	
	/*
	virtual void calc_Jacobian(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalMatrix & J)
	{
		assert(J.cols() == this->get_n() && J.rows() == this->get_n());

		// getting unknowns
		double Sg = vec_unknowns(0);
		double rho_L_h = vec_unknowns(1);
		double rho_G_h = vec_unknowns(2);

		// evaluate J
		J.setZero();
		J(0, 0) = (Calc_Res_Sg(Sg + eps, rho_L_h, rho_G_h) - Calc_Res_Sg(Sg - eps, rho_L_h, rho_G_h)) / eps / 2;
		J(0, 1) = (Calc_Res_Sg(Sg, rho_L_h + eps, rho_G_h) - Calc_Res_Sg(Sg, rho_L_h - eps, rho_G_h)) / eps / 2;
		J(0, 2) = (Calc_Res_Sg(Sg, rho_L_h, rho_G_h + eps) - Calc_Res_Sg(Sg, rho_L_h, rho_G_h - eps)) / eps / 2;
		J(1, 0) = (Calc_Res_rho_L_h(Sg + eps, rho_L_h) - Calc_Res_rho_L_h(Sg - eps, rho_L_h)) / eps / 2;
		J(1, 1) = (Calc_Res_rho_L_h(Sg, rho_L_h + eps) - Calc_Res_rho_L_h(Sg, rho_L_h - eps)) / eps / 2;
		J(1, 2) = 0.0;
		J(2, 0) = (Calc_Res_rho_G_h(Sg + eps, rho_G_h) - Calc_Res_rho_G_h(Sg - eps, rho_G_h)) / eps / 2;
		J(2, 1) = 0.0;
		J(2, 2) = (Calc_Res_rho_G_h(Sg, rho_G_h + eps) - Calc_Res_rho_G_h(Sg, rho_G_h - eps)) / eps / 2;
	};
	*/
	/**
	* based on 
	*/
	virtual void calc_Rest_SecVar(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalVector & vec_rest_var)
	{
		vec_rest_var = ogsChem::LocalVector::Zero(4);
		double Sg = vec_unknowns(0);
		double rho_L_h = vec_unknowns(1);
		double rho_G_h = vec_unknowns(2);
		vec_rest_var(0) = getPcbySg(Sg);//PC
		vec_rest_var(1) = getPL(Sg);//PL
		vec_rest_var(2) = getPG(Sg);//PG
		vec_rest_var(3) = get_RHO_G_W(Sg, rho_G_h);
	};
	/**
	  * realization of the calc_Jacobian function.
	  * this function will evaluate the Jacobian matrix,
	  * based on the unknown values given.
	  */


	virtual void set_env_condition(ogsChem::LocalVector & env_conditon)
	{
		// mean pressure
		P_L = env_conditon(0);
		// molar fraction of the light component
		X_L = env_conditon(1);
	}; 
	/**
	* realization of the calc the residual of saturation
	* EOS--Function 1
	* X=Sl*rho_L_h+Sg*rho_G_h
	*/
	virtual double Calc_Res_Sg(double Sg, double rho_L_h, double rho_G_h)
	{
		double Res_Sg(0.0);
		Res_Sg = X_L - ((1 - Sg)*rho_L_h + Sg*rho_G_h);
		return Res_Sg;
	}
	/**
	* realization of the calc the residual of rho_L_h
	* EOS--Function 2
	* rho_L_h=min(C_h*p_g^h,X_L)
	*/
	virtual double Calc_Res_rho_L_h(double Sg, double rho_L_h)
	{
		double Res_rho_L_h(0.0);
		Res_rho_L_h = std::min(Sg, -rho_L_h + get_RHO_L_H(Sg));
		return Res_rho_L_h;
	}
	/**
	* realization of the calc the residual of rho_L_h
	* EOS--Function 3
	* rho_G_h=MAX(C_V*p_g^h,X_L)
	*/
	virtual double Calc_Res_rho_G_h(double Sg, double rho_G_h)
	{
		double Res_rho_G_h(0.0);
		Res_rho_G_h = std::min((1-Sg), rho_G_h - get_RHO_G_H(Sg));
		return Res_rho_G_h;
	}

	virtual double getS_bar(double Sg)
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
	virtual double getPc_bar_vG_Sg(double Sg)
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
	virtual double get_dPCdS_vG_bar(double Sg)
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
		/*double Pc = 0.0;
		if (Sg >= S_gr && Sg <= 1 - S_lr)
		{
			double S_le = (1 - Sg - S_lr) / (1 - S_gr - S_lr);
			Pc = P_r*pow(pow(S_le, -1 / m) - 1, 1 / n);
		}
		else if(Sg < S_gr)
		{
			Pc = 0.0;
		}
		else
		{
			Pc = getPc_bar_vG_Sg(1 - S_lr) + get_dPCdS_vG_bar(1 - S_lr)*(Sg - 1 + S_lr);
		}
		return Pc;*/
		double Pc = 0.0;
		if (Sg >= S_gr && Sg <= 1 - S_lr)
		{
			Pc = getPc_bar_vG_Sg(Sg);
		}
		else if (Sg < S_gr)
		{
			Pc = getPc_bar_vG_Sg(S_gr) + get_dPCdS_vG_bar(S_gr)*(Sg - S_gr);
		}
		else
		{
			Pc = getPc_bar_vG_Sg(1 - S_lr) + get_dPCdS_vG_bar(1 - S_lr)*(Sg - 1 + S_lr);
		}
		return Pc;
	}
	/**
	* Derivation of PC (capillary pressure) in terms of saturation based on van Genuchten model
	*/
	virtual double Deriv_dPCdS(double Sg)
	{
		/*double dPcdSg(0.0);
		if (Sg >= S_gr && Sg <= 1 - S_lr)
		{
			//dPcdSg = get_dPCdS_vG_bar(Sg);
			double S_le = (1 - Sg - S_lr) / (1 - S_gr - S_lr);
			dPcdSg=P_r*(1 / (m*n))*(1 / (1 - S_lr - S_gr))*pow(pow(S_le, (-1 / m)) - 1, (1 / n) - 1)*pow(S_le, (-1 / m)) / S_le;
		}
		else if (Sg < S_gr)
		{
			dPcdSg = 0.0;
		}
		else
		{
			dPcdSg = get_dPCdS_vG_bar(1 - S_lr);
		}
		return dPcdSg;
		*/
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
		return dPcdSg;
	}
	/**
	* calculate the effective saturation of Liquid phase by Gas saturation
	*/
	virtual double getEffectSat_lbySg(double Sg)
	{
		double EffectSat_l = 0.0;
		
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
		double gamma = getweightedFunc(Sg);
		PG = P_L + PC;// (1 - gamma)*PC;
		return PG;
	}

	/**
	* based on primary variables to mean_pressure get the value of PL ---Liquid phase pressure
	*/
	virtual double getPL(double Sg)
	{
		double PL = 0.0;
		double PC = 0.0;
		PC = getPcbySg(Sg);
		double gamma = getweightedFunc(Sg);
		PL = P_L - gamma*PC;
		return PL;
	}
	/**
	* weighted function 
	*/
	virtual double getweightedFunc(double Sg)
	{
		double gamma=0.0;
		
		/*if (Sg >= 1)
			gamma = 1;
		else if (Sg <= 0)
			gamma = 0;
		else
			gamma = Sg*(2 - Sg);
		*/
		return gamma;
	}
	/**
	* Derivation of weighted function based on Sg
	*/
	virtual double get_deriv_weight_Func(double Sg)
	{
		double dgamma_dSg(0.0);
		/*
		dgamma_dSg = 2 - 2 * Sg;
		if (Sg >= 1)
			dgamma_dSg = 0;
		else if (Sg <= 0)
			dgamma_dSg = 0;
			*/
		return dgamma_dSg;
	}

	/**
	* based on Sg  get the value of dissoved hydrogen mass density in the LIQUID phase
	*/
	virtual double get_RHO_L_H(double sg)
	{
		double rho_L_h(0.0);
		double PG(0.0);
		PG = getPG(sg);
		//double PG_H(0.0);
		//PG_H = Function_PG_h(PG);
		rho_L_h = std::min(C_h*PG, X_L);//Henry law
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
		//double PG_H(0.0);
		//PG_H = Function_PG_h(PG);
		rho_G_h = C_v*PG;// std::max(C_v*PG_H, X_L);//ideal gas law
		return rho_G_h;
	}
	/**
	* based on Sg  get the value of water vapor mass density in the GAS phase
	*/
	virtual double get_RHO_G_W(double sg, double rho_G_h)
	{
		double rho_G_w(0.0);
		/*double PG(0.0);
		PG = getPG(sg); 
		//rho_G_w = C_w*PG - C_w*rho_G_h / C_v;
		rho_G_w = C_w*PG - M_L*rho_G_h / M_G;*/
		return rho_G_w;
	}
	/**
	* based on pg  get the value of P_G^h ---PARTIAL GAS PRESSURE OF HYDROGEN
	*/

	virtual double Function_PG_h(double PG)
	{
		double PG_H(0.0);
		double P_vapor(0.0);
		P_vapor = get_Vapor_Pressure(T0-273);
		//PG_H = (1 / (2 * C_h*M_L))* (C_h*M_L*PG - M_G* rho_L_std + sqrt(pow((-C_h*M_L*PG + M_G*rho_L_std), 2) - 4 * C_h*M_L*(-M_G*PG*rho_L_std + M_G* P_vapor* rho_L_std)));
		PG_H = PG;
		return PG_H;

	}


	/**
	* derivative dP_G^HdPG
	*/
	virtual double Deriv_dPGH_dPG(double Sg)
	{
		double PG_h(0.0);
		double PG(0.0);
		double dPGH_dPG(0.0);
		PG = getPG(Sg);
		PG_h = Function_PG_h(PG);
		double P_vapor(0.0);
		P_vapor = get_Vapor_Pressure(T0-273);
		double A = 0.0;
		A = M_L / M_G;
		double B = 0.0;
		B = P_vapor*rho_L_std*A*C_h / (rho_L_std + A*C_h*PG_h);
		dPGH_dPG = 1 / (1 - B);
		return dPGH_dPG;
	}
	/**
	* Intermedia variables Alpha and Beta
	*/

	virtual double Function_Alpha(double Sg)
	{
		double dPGH_dPG(0.0);
		dPGH_dPG = Deriv_dPGH_dPG(Sg);
		double Alpha(0.0);
		Alpha = ((C_v - C_h)*Sg + C_h)*dPGH_dPG;
		return Alpha;
	}
	virtual double Function_Beta(double Sg)
	{
		double PGh(0.0);
		double PG(0.0);
		double Beta(0.0);
		PG = getPG(Sg);
		PGh = Function_PG_h(PG);
		Beta = (C_v - C_h)*PGh;
		return Beta;
	}

	/**
	* Derivation of PG based on saturation
	*/
	virtual double Deriv_dPGdSg(double Sg, double n = 1.49, double m = 1.0 - 1.0 / 1.49)
	{
		double dPcdSg(0.0);
		dPcdSg = Deriv_dPCdS(Sg);
		double PC(0.0);
		PC = getPcbySg(Sg);
		double Gamma(0.0);
		Gamma = getweightedFunc(Sg);
		double dPGdSg(0.0);
		dPGdSg = dPcdSg*(1 - Gamma) - PC*get_deriv_weight_Func(Sg);
		//if (Sg < 1e-14)
			//dPGdSg = 0.0;
		return dPGdSg;
	}
	
	virtual double Deriv_dSgdP(double Sg)
	{
		
		double alpha(0.0);
		alpha = Function_Alpha(Sg);
		double beta(0.0);
		beta = Function_Beta(Sg);
		double dPG_dSg(0.0);
		dPG_dSg = Deriv_dPGdSg(Sg);
		double dSgdP(0.0);
		dSgdP = -alpha/ (beta + alpha*dPG_dSg);
		return dSgdP;
	}
	virtual double Deriv_dSgdX(double Sg)
	{
		double alpha(0.0);
		alpha = Function_Alpha(Sg);
		double beta(0.0);
		beta = Function_Beta(Sg);
		double dPG_dSg(0.0);
		dPG_dSg = Deriv_dPGdSg(Sg);
		double dSgdX(0.0);
		dSgdX = 1/ (beta + alpha*dPG_dSg);
		return dSgdX;
	}
	/**
	* The intermedia function
	* Omega=beta(p,Sg)*Func_characteristic(P,X)/(beta(P,Sg)+alpha(P,Sg)*dPgdSg(P,Sg))
	* M=Func_characteristic(P,X)*dPgdSg(P,Sg)/(beta(P,Sg)+alpha(P,Sg)*dPgdSg(P,Sg))
	*/
	virtual double Func_Omega(double Sg)
	{
		double alpha(0.0);
		alpha = Function_Alpha(Sg);
		double beta(0.0);
		beta = Function_Beta(Sg);
		double Omega(0.0);
		double dPG_dSg(0.0);
		dPG_dSg = Deriv_dPGdSg(Sg);
		Omega = beta/ (beta + alpha*dPG_dSg);
		return Omega;

	}
	virtual double Func_InterM(double Sg)
	{
		double M(0.0);
		double dPG_dSg(0.0);
		dPG_dSg = Deriv_dPGdSg(Sg);
		double alpha(0.0);
		alpha = Function_Alpha(Sg);
		double beta(0.0);
		beta = Function_Beta(Sg);
		M = 1 / (beta / dPG_dSg + alpha);
		return M;

	}
	/**
	* Deriv rho_L_h based on mean pressure
	*/
	virtual double Deriv_drhoLh_dP(double Sg, double Omega)
	{
		double drho_L_h_dP(0.0);
		double dPG_H_dPG(0.0);
		dPG_H_dPG = Deriv_dPGH_dPG(Sg);
		drho_L_h_dP = C_h*dPG_H_dPG*Omega;
		return drho_L_h_dP;
	}

	/**
	* Deriv rho_L_h based on total mass density
	*/
	virtual double Deriv_drhoLh_dX(double Sg, double M)
	{
		double drho_L_h_dX(0.0);
		double dPG_H_dPG(0.0);
		dPG_H_dPG = Deriv_dPGH_dPG(Sg);
		drho_L_h_dX = C_h*dPG_H_dPG*M;
		return drho_L_h_dX;
	}
	/**
	* Deriv rho_G_h based on P mean pressure
	*/
	virtual double Deriv_drhoGh_dP(double Sg, double Omega)
	{
		double drho_G_h_dP(0.0);
		double dPG_H_dPG(0.0);
		dPG_H_dPG = Deriv_dPGH_dPG(Sg);
		drho_G_h_dP = C_v*dPG_H_dPG*Omega;
		return drho_G_h_dP;
	}

	/**
	* Deriv rho_G_h based on X --total mass density
	*/
	virtual double Deriv_drhoGh_dX(double Sg, double M )
	{
		double drho_G_h_dX(0.0);
		double dPG_H_dPG(0.0);
		dPG_H_dPG = Deriv_dPGH_dPG(Sg);
		drho_G_h_dX = C_v*dPG_H_dPG*M;
		return drho_G_h_dX;
	}

	/**
	* The characterastic function
	* if return 1 indicates the two phase region
	* if return 0 indicates the single phase region
	*/
	virtual std::size_t Func_Characteristic(double mean_P, double total_X)
	{
		std::size_t Func_C(0);
		double PG_h_max(0.0);
		double PG_h_min(0.0);
		double PG(0.0);
		PG_h_min = mean_P + getPcbySg(0.0);
		PG_h_max = mean_P + getPcbySg(0.60);
		if (total_X<C_v*PG_h_max && total_X> C_h*PG_h_min)
		{
			Func_C = 1;
		}
		else
		{
			Func_C = 0;
		}
		return Func_C;
	}
	virtual std::size_t Func_Characteristic_1(double X)
	{
		std::size_t Func_C(0);
		double PG_h_max(0.0);
		double PG_h_min(0.0);
		double PG(0.0);
		PG_h_min = Function_PG_h(getPG(0));
		PG_h_max = Function_PG_h(getPG(0.6));
		if (X<C_v*PG_h_max && X> C_h*PG_h_min)
		{
			Func_C = 1;
		}
		else
		{
			Func_C = 0;
		}
		return Func_C;
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
		/*double tmp(0.0);
		tmp = 2.786 + 0.031514*T - 1.2373e-4*pow(T, 2) + 4.2267e-7*pow(T, 3) - 8.1308e-10*pow(T, 4);
		P_vapor = pow(10, tmp);*/
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

	double _alpha;
	double _beta;

	double xi = 1e-5;
};


