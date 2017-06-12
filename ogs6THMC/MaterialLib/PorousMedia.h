/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file PorousMedia.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "IMedium.h"
#include "NumLib/Function/FunctionLinear.h"
#include "logog.hpp"
#include "math.h"

namespace NumLib
{
class ITXFunction;
class FunctionLinear1D;
}

namespace MaterialLib
{

struct PorousMedia : public IMedium
{
    NumLib::ITXFunction* hydraulic_conductivity;
    NumLib::ITXFunction* permeability;
    NumLib::ITXFunction* porosity;
    NumLib::ITXFunction* storage;
    NumLib::ITXFunction* geo_area;
    NumLib::ITXFunction* dispersivity_long;
    NumLib::ITXFunction* dispersivity_trans; 
	NumLib::ITXFunction* res_saturation; 
	NumLib::ITXFunction* res_saturation_g;
	NumLib::ITXFunction* res_saturation_l;
	/**
	*PARAMETER DEFINED FOR MULTIPHASE
	*g/L means gas phase/Liquid phase respectively
	*/
	NumLib::ITXFunction* Diffusion_Coefficient_G;//
	NumLib::ITXFunction* Diffusion_Coefficient_L;
	NumLib::ITXFunction* Viscosity_G;
	NumLib::ITXFunction* Viscosity_L;

    NumLib::ITXFunction* max_saturation;      
    NumLib::ITXFunction* exp_saturation;      
    NumLib::FunctionLinear1D* capp_sat_curve; 
    std::vector<size_t> perm_saturation_model;    
    std::vector<NumLib::FunctionLinear1D*> perm_saturation_curve; 
    size_t capp_sat_model;          
	double minimum_relative_permeability;  

	double Pb_van_Genuchten;
	double S_lr;
	double S_gr;
	double n_van_Genuchten;

	double P_entry_Brook;
	double Lambda_Brook;
    PorousMedia()
    {
		BaseLib::zeroObject(
			hydraulic_conductivity,
			permeability,
			porosity,
			storage,
			geo_area,
			res_saturation, //?
			max_saturation, //?
			exp_saturation
                ); 
		BaseLib::zeroObject(
			    res_saturation_g,
			    res_saturation_l
			    );
        BaseLib::zeroObject(
                capp_sat_curve,
                dispersivity_long,
                dispersivity_trans
                );
		
        capp_sat_model = 0;
        minimum_relative_permeability = .0;
		Pb_van_Genuchten=0.0;
		n_van_Genuchten=0.0;
		P_entry_Brook=0.0;
		Lambda_Brook=0.0;
		S_gr=0.0;
		S_lr=0.0;
    }

    virtual ~PorousMedia()
    {
        BaseLib::releaseObject(
                hydraulic_conductivity,
                permeability,
                porosity,
                storage,
                geo_area,
                res_saturation,
				max_saturation,
                exp_saturation
                );
		BaseLib::releaseObject(
			res_saturation_g,
			res_saturation_l
			);
        perm_saturation_model.clear();
        BaseLib::releaseObjectsInStdVector(perm_saturation_curve);
        BaseLib::releaseObject(capp_sat_curve);
		BaseLib::releaseObject(
			Diffusion_Coefficient_G,
			Diffusion_Coefficient_L,
			Viscosity_G,
			Viscosity_L
			);
    }

    virtual MediumType::type getMediumType() const {return MediumType::PorousMedium;};

	virtual double getKrelbySw(const double Sw/*wetting saturation*/, size_t idx_phase)
    {
        double kr = 0.0; 
        //bool phase_shift = false;
        size_t model; 
        model = this->perm_saturation_model[idx_phase];
		minimum_relative_permeability =  1.0e-9; // Could be set as default at a beter place
		
        switch(model)
        {
	     default: // Error occurs! 
		    ERR("ERROR in PermeabilitySaturationFunction(): Unrecognized relative permeability method.\n");
            exit(0);
	        break;
	    case 0: // CURVE
			this->perm_saturation_curve[idx_phase]->eval( Sw, kr ); 
            if( kr < minimum_relative_permeability )
	        kr = minimum_relative_permeability;
            break;
        }
        return kr; 
    }

    /**
    * return the water saturation value by capillary pressure
    */
    virtual double getSwbyPc(double Pc)
    {	
        double Sw = 0.0;

        switch ( this->capp_sat_model )
        {
        case 0:  // curve value
            // get the value from curve. 
            capp_sat_curve->eval( Pc, Sw );
            break;
        default: 
            ERR("No valid capilary pressure vs water saturation model! "); 
            exit(0); 
            break;
        } // end of switch capp_sat_model
        return Sw; 
    }


	virtual double getdSwdPc (double Pc)
    {
        double dSwdPc = 0.0;

        switch ( this->capp_sat_model )
        {
        case 0:  // curve value
            // get the value from curve. 
            capp_sat_curve->eval_slope( Pc, dSwdPc );
            break;
        default: 
            ERR("Error in getSwbyPc: No valid capilary pressure vs water saturation model! "); 
            exit(0); 
            break;
        } // end of switch capp_sat_model
        return dSwdPc; 
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
	virtual double getEffectSat_lbySg(double Sg)
	{
		double EffectSat_l = 0.0;
		// this->res_saturation_g->eval(Sg, res_S_g);

		if (Sg < S_gr)
			EffectSat_l = (1 - Sg - S_lr) / (1 - S_gr - S_lr);
		else if (Sg>(1 - S_lr))
			EffectSat_l = (1 - Sg - S_lr) / (1 - S_gr - S_lr);
		else
			EffectSat_l = (1 - Sg - S_lr) / (1 - S_gr - S_lr);

		return EffectSat_l;
	}

	/**
	* calculate the effective saturation of gas phase by Gas saturation
	*/
	virtual double getEffectSat_gbySg(double Sg)
	{
		double EffectSat_g = 0.0;
		// this->res_saturation_g->eval(Sg, res_S_g);

		EffectSat_g = (Sg - S_gr) / (1 - S_gr - S_lr);

		return EffectSat_g;
	}

	/**
	* get the relative permeability based on Sg
	*n = 1.49 is for MoMaS Benchmark 2
	* n=1.54 is for  MoMaS Benchmark 3
	*/
	virtual double getKr_L_bySg(double Sg/*gas phase saturation*/)//double m = 1.0 - 1.0 / 1.49, double res_S_l = 0.40, double res_S_g = 0.0 
	{
		double Kr_L = 0.0;
		double EffectSat_l = 0.0;
		double m = 0.0;
		EffectSat_l = getEffectSat_lbySg(Sg);//for MoMaS benchmark 2

		switch (capp_sat_model)
		{
		case 4:
			m = 1.0 - 1.0 / n_van_Genuchten;
			Kr_L = sqrt(EffectSat_l)*pow(1 - pow(1 - pow(EffectSat_l, 1 / m), m), 2);
			if (Sg > 1)
				Kr_L = 0;
			else if (Sg < 0)
				Kr_L = 1;
			break;
		case 6:
			// krw = pow(se,3.0+2.0/m)
			// Se  = (sl - slr) / (slm - slr)
			Kr_L = pow(EffectSat_l, 3.0 + 2.0 / Lambda_Brook);
			//Kr_L = pow(1 - Sg, 2);
			if (Sg > 1)
				Kr_L = 0;
			else if (Sg < 0)
				Kr_L = 1;
			break;
		default:
			ERR("Error in getKr_L_bySg: No valid relative permeability vs saturation model! ");
			break;
		}
		//if (Kr_L < 1e-9)
			//Kr_L = 1e-9;
		
		return Kr_L;
	}

	/**
	* get the relative permeability of gas phase based on Sg
	*/
	virtual double getKr_g_bySg(double Sg/*gas phase saturation*/)
	{
		double Kr_G = 0.0;

		double EffectSat_g = 0.0;
		double EffectSat_l = 0.0;
		EffectSat_g = getEffectSat_gbySg(Sg);
		EffectSat_l = getEffectSat_lbySg(Sg);
		double m = 0.0;
		switch (capp_sat_model)
		{
		case 4:// index 4 stands for the van Genuchten Model
			m = 1.0 - 1.0 / n_van_Genuchten;
			Kr_G = sqrt(EffectSat_g)*pow(1 - pow(1 - EffectSat_g, 1 / m), 2 * m);
			if (Sg > 1)
				Kr_G = 1;
			else if (Sg < 0)
				Kr_G = 0;
			break;
		case 6: //index 6 stands for the Brooks corey  model
		// krg = pow(1.0-se,2)*(1.0-pow(se,1.0+2.0/m))
		// Se  = (1-Sg - S_lr) / (1-S_lr - S_gr)
			Kr_G = pow(1.0 - EffectSat_l, 2)*(1.0 - pow(EffectSat_l, 1.0 + 2.0 / Lambda_Brook));
			//Kr_G = pow(Sg, 2);
			if (Sg > 1)
				Kr_G = 1;
			else if (Sg < 0)
				Kr_G = 0;
			break;
		default:
			ERR("Error in getKr_g_bySg: No valid relative permeability vs saturation model! ");
			break;
		}
		return Kr_G;
	}
	virtual double getSat_byPC(double P_c)//for MoMaS benchmark2 n = 1.49
	{
		
		double Sg = 0.0;
		double m = 0.0;
		switch (capp_sat_model)
		{
		case 4:// index 4 stands for the van Genuchten Model

			m = 1.0 - 1.0 / n_van_Genuchten;
			if (P_c > 0){
				Sg = 1 - (((1 - S_gr - S_lr) / pow((pow((P_c / Pb_van_Genuchten), n_van_Genuchten) + 1), m)) + S_lr);
			}
			else
				Sg = 0.0;//just for test
			break;
		case 6: // index 6 stands for the Brooks Corey Model
			if (P_c<P_entry_Brook)
			{
				Sg = (P_c - P_entry_Brook)*Lambda_Brook / P_entry_Brook;
				//Sg = 1 - (S_lr + (1 - S_lr - S_gr)*pow(P_entry_Brook / P_c, Lambda_Brook));
			}
			else if (P_c > 4 * P_entry_Brook)
			{
				Sg = 1 - 1 / pow(4, Lambda_Brook) + Lambda_Brook*(P_c - 4 * P_entry_Brook) / P_entry_Brook / pow(4, 1 + Lambda_Brook);
			}

			else
			{
				Sg = 1 - (S_lr + (1 - S_lr - S_gr)*pow(P_entry_Brook / P_c, Lambda_Brook));// 1 - pow(P_entry_Brook / P_c, Lambda_Brook);
			}
			break;
		default:
			ERR("Error in getSat_byPC: No valid capillary pressure vs saturation model! ");
			break;
		}
		//if (Sg < 1e-5)
		//Sg = 1e-5;
		return Sg;

	}

	virtual double get_diriv_PC(double P_c)//for MoMaS benchmark2 n = 1.49
	{
		double dPC(0.0);
		double m = 0.0;
		switch (capp_sat_model)
		{
		case 4:
			m = 1.0 - 1.0 / n_van_Genuchten;
			dPC = m*n_van_Genuchten*(1 - S_gr - S_lr)*(1 / Pb_van_Genuchten)*(pow((P_c / Pb_van_Genuchten), (n_van_Genuchten - 1)))*pow(((pow((P_c / Pb_van_Genuchten), n_van_Genuchten)) + 1), (-m - 1));
			if (P_c <= 0)
			{
				dPC = 0.0;//just for test
			}
			break;
		case 6:
			//dPC = (1 - S_lr - S_gr)*Lambda_Brook*pow(P_entry_Brook, Lambda_Brook)*pow(P_c, -Lambda_Brook - 1);
			if (P_c <= P_entry_Brook)
			{
				dPC = Lambda_Brook / P_entry_Brook;
				//Sg = 1 - (S_lr + (1 - S_lr - S_gr)*pow(P_entry_Brook / P_c, -Lambda_Brook));
			}
			else if (P_c > 4 * P_entry_Brook)
			{
				dPC = Lambda_Brook / P_entry_Brook / pow(4, 1 + Lambda_Brook);
			}
			else
			{
				dPC = (1 - S_lr - S_gr)* Lambda_Brook*pow(P_entry_Brook, Lambda_Brook)*pow(P_c, -Lambda_Brook - 1);

			}
			break;
		default:
			ERR("Error in get_diriv_PC: No valid inverse capillary pressure vs saturation model! ");
			break;

		}
		//if (dPC < 1e-5)
		//dPC = 1e-5;
		return dPC;

	}

	virtual double getPc_bySat(double Sg)
	{

		double P_c = 0.0;
		double m = 0.0;
		switch (capp_sat_model)
		{
		case 4:
			//m = 1.0 - 1.0 / n_van_Genuchten;
			if (Sg >= S_gr && Sg <= 1 - S_lr)
			{
				P_c = getPc_bar_vG_Sg(Sg);
			}
			else if (Sg < S_gr)
			{
				P_c = getPc_bar_vG_Sg(S_gr) + get_dPCdS_vG_bar(S_gr)*(Sg - S_gr);
			}
			else
			{
				P_c = getPc_bar_vG_Sg(1 - S_lr) + get_dPCdS_vG_bar(1 - S_lr)*(Sg - 1 + S_lr);
			}
			break;
		case 6:
			//dPC = (1 - S_lr - S_gr)*Lambda_Brook*pow(P_entry_Brook, Lambda_Brook)*pow(P_c, -Lambda_Brook - 1);
			if (Sg < S_gr)
			{
				P_c = ((Sg / (1 - S_lr))*P_entry_Brook / Lambda_Brook) + P_entry_Brook;
				//Sg = 1 - (S_lr + (1 - S_lr - S_gr)*pow(P_entry_Brook / P_c, -Lambda_Brook));
			}
			else if (Sg>S_lr)
			{
				P_c = pow(4, 1 + Lambda_Brook)*P_entry_Brook*((1 / pow(4, Lambda_Brook)) - (1 - Sg - S_lr) / (1 - S_lr)) / Lambda_Brook + 4 * P_entry_Brook;
			}
			else
			{
				P_c = P_entry_Brook*pow((1 - Sg - S_lr) / (1 - S_lr), -1 / Lambda_Brook);

			}
			break;
		default:
			ERR("Error in getPc_bySat: No valid  capillary pressure vs saturation model! ");
			break;

		}
		//if (dPC < 1e-5)
		//dPC = 1e-5;
		return P_c;
	}
	virtual double Deriv_dPCdS(double Sg)
	{
		double dPcdSg(0.0);
		switch (capp_sat_model)
		{
		case 4:
			//m = 1.0 - 1.0 / n_van_Genuchten;
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
			break;
		case 6:
			//dPC = (1 - S_lr - S_gr)*Lambda_Brook*pow(P_entry_Brook, Lambda_Brook)*pow(P_c, -Lambda_Brook - 1);
			if (Sg < S_gr)
			{
				dPcdSg = P_entry_Brook / Lambda_Brook / (1 - S_lr);// ((Sg / 1 - S_lr)*P_entry_Brook / Lambda_Brook) + P_entry_Brook;
				//Sg = 1 - (S_lr + (1 - S_lr - S_gr)*pow(P_entry_Brook / P_c, -Lambda_Brook));
			}
			else if (Sg>S_lr)
			{
				dPcdSg = pow(4, 1 + Lambda_Brook)*P_entry_Brook / Lambda_Brook / (1 - S_lr);
			}
			else
			{
				dPcdSg = (1 / Lambda_Brook)*P_entry_Brook*pow((1 - Sg - S_lr) / (1 - S_lr), -1 - 1 / Lambda_Brook) / (1 - S_lr);

			}
			break;
		default:
			ERR("Error in getPc_bySat: No valid  capillary pressure vs saturation model! ");
			break;

		}

		return dPcdSg;
	}

	/// modified van Genuchten model
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
		double m = 1 - 1 / n_van_Genuchten;
		//effective saturation
		S_le = (1 - Sg - S_lr) / (1 - S_gr - S_lr);
		//Pc_vG = P_r*(S_le. ^ (-1 / m) - 1). ^ (1 / n);
		Pc_vG = Pb_van_Genuchten*pow(pow(S_le, -1 / m) - 1, 1 / n_van_Genuchten);
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
		double m = 0.0;
		m = 1 - 1 / n_van_Genuchten;
		double S_le = (1 - Sg - S_lr) / (1 - S_gr - S_lr);
		dPcdSg = Pb_van_Genuchten*(1 / (m*n_van_Genuchten))*(1 / (1 - S_lr - S_gr))*pow(pow(S_le, (-1 / m)) - 1, (1 / n_van_Genuchten) - 1)*pow(S_le, (-1 / m)) / S_le;
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
	* Define three character function for calculate the gradient
	*/
	std::size_t Characteristic_Function_xi(double Sg)
	{
		size_t xi(0);
		if (std::abs(Sg)>1e-14 && std::abs(Sg-1)<1e-14)
			xi = 1;
		else
			xi = 0;
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
		/**
		* parameters for van-Genuchten capillary pressure-saturation model
		*/
		/*
		const double S_gr = 0.0;//residual saturation of the gas phase
		const double S_lr = 0.4;//residual saturation of the liquid phase
		const double P_r = 2e+6;//pa
		const double n = 1.49;
		const double m = 1 - 1 / n;*/
		double xi=1e-5;
};

} //end
 
