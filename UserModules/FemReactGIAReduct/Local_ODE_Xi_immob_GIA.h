/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Local_ODE_Xi_immob.h
 *
 * Created on 2013-17-09 by Reza Zolfaghari and Haibing Shao
 */
 
#ifndef LOCAL_ODE_XI_IMMOB_GIA_H
#define LOCAL_ODE_XI_IMMOB_GIA_H

#include "ChemLib/chemReductionGIA.h"

class Local_ODE_Xi_immob_GIA
{
public: 

	/**
      * constructor
      */ 
    Local_ODE_Xi_immob_GIA(ogsChem::chemReductionGIA* reductionGIA)
		: _reductionGIA(reductionGIA), _n_eta( _reductionGIA->get_n_eta()), _n_etabar(_reductionGIA->get_n_eta_bar()),
		  _n_xi_global(_reductionGIA->get_n_xi_global()), _n_xi_local(_reductionGIA->get_n_xi_local()),_n_xi_Kin_bar(_reductionGIA->get_n_xi_Kin_bar()), _J_tot_kin(_reductionGIA->get_n_xi_Kin_total())
    {
		_vec_eta_mob      = MathLib::LocalVector::Zero(_n_eta);
	    _vec_eta_immob    = MathLib::LocalVector::Zero(_n_etabar);
	    _vec_xi_global    = MathLib::LocalVector::Zero(_n_xi_global);
	    _vec_xi_local     = MathLib::LocalVector::Zero(_n_xi_local);
		_vec_dxi_immob_dt = MathLib::LocalVector::Zero( _J_tot_kin );
		_vec_dxi_immob_dt_new = MathLib::LocalVector::Zero( _n_xi_Kin_bar );
		_vec_loc_conc     = MathLib::LocalVector::Zero(_reductionGIA->get_n_Comp());
	}
	 
    /**
      * destructor
      */ 
	~Local_ODE_Xi_immob_GIA()
    {
        _reductionGIA = NULL;
    }

    /**
      * update the eta and xi values
      */ 
	void update_eta_xi( MathLib::LocalVector & vec_eta_mob, 
		                MathLib::LocalVector & vec_eta_immob, 
						MathLib::LocalVector & vec_xi_global,
						MathLib::LocalVector & vec_xi_local)
	{
		_vec_eta_mob   = vec_eta_mob; 
		_vec_eta_immob = vec_eta_immob; 
	    _vec_xi_global = vec_xi_global;
	    _vec_xi_local  = vec_xi_local;
	}

    /**
      * evaluate the change of xi_immob over time
      */ 
	MathLib::LocalVector operator() (double /*time*/, MathLib::LocalVector vec_xi_kin_bar )
	{

		const double theta_water_content(1.0);  // HS: testing, will be removed. 
		MathLib::LocalMatrix mat_A2kin = _reductionGIA->get_matrix_A2kin();

		_vec_xi_local.tail(_n_xi_Kin_bar) = vec_xi_kin_bar;
		_reductionGIA->EtaXi2Conc(_vec_eta_mob, _vec_eta_immob, _vec_xi_global, _vec_xi_local, _vec_loc_conc); 
		this->_reductionGIA->Calc_Kin_Rate(_vec_loc_conc, _vec_dxi_immob_dt);
		_vec_dxi_immob_dt_new  = (theta_water_content * mat_A2kin) * _vec_dxi_immob_dt;
		return _vec_dxi_immob_dt_new;
	}

private: 

    /**
      * pointer to reduction scheme class
      */ 
	ogsChem::chemReductionGIA* _reductionGIA;
	
    /**
      * local vector of eta mobile
      */ 
    MathLib::LocalVector _vec_eta_mob; 

    /**
      * local vector of eta immobile or bar
      */ 
    MathLib::LocalVector _vec_eta_immob; 
    
    /**
      * local vector of xi global
      */ 
    MathLib::LocalVector _vec_xi_global;

    /**
      * local vector of xi local
      */ 
	MathLib::LocalVector _vec_xi_local;
	
	/**
	  * local vector of concentration values
	  */
	MathLib::LocalVector _vec_loc_conc; 

    /**
      * local vector of dxi_immob_dt
      */ 
	MathLib::LocalVector _vec_dxi_immob_dt, _vec_dxi_immob_dt_new;

	std::size_t _n_eta, _n_etabar, _n_xi_global, _n_xi_local,_n_xi_Kin_bar, _J_tot_kin;
 
}; 
 
 #endif
