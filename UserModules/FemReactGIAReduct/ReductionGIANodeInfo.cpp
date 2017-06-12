/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SingleStepKinReduction.h
 *
 * Created on 24.05.2013 by Reza Zolfaghari and Haibing Shao
 */
 
#include "ReductionGIANodeInfo.h"
 
ReductionGIANodeInfo::ReductionGIANodeInfo(size_t node_id,
	                                       size_t n_comp, 
	                                       size_t n_eta,
										   size_t n_eta_bar,
										   size_t n_xi_global,
										   size_t n_xi_local,
										   ogsChem::chemReductionGIA* ReductionGIA)
    : _node_id(node_id), _n_comp(n_comp), _n_eta(n_eta),
	  _n_eta_bar(n_eta_bar), _n_xi_global(n_xi_global), _n_xi_local(n_xi_local), _ReductionGIA(ReductionGIA)
{
	_Comp_Conc = MathLib::LocalVector::Zero( _n_comp );
	loc_eta  	   = MathLib::LocalVector::Zero( _n_eta );   //changing local variable names to avoid the confusion with global variables.
	loc_eta_bar   = MathLib::LocalVector::Zero( _n_eta_bar );
	loc_xi_global = MathLib::LocalVector::Zero( _n_xi_global );
	loc_xi_local  = MathLib::LocalVector::Zero( _n_xi_local );
}

ReductionGIANodeInfo::~ReductionGIANodeInfo()
{
    _ReductionGIA = NULL;
}

void ReductionGIANodeInfo::set_comp_conc( size_t comp_idx, double val )
{
	_Comp_Conc(comp_idx) = val; 
}

void ReductionGIANodeInfo::transform(void)
{
	_ReductionGIA->Conc2EtaXi( _Comp_Conc, loc_eta, loc_eta_bar, loc_xi_global, loc_xi_local );
}

