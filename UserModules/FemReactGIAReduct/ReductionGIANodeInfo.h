/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SingleStepKinReduction.h
 *
 * Created on 24.05.2013 by Reza Zolfaghari & Haibing Shao
 */

#ifndef REDUCTION_GIA_NODE_INFO_H
#define REDUCTION_GIA_NODE_INFO_H

#include "ChemLib/chemReductionGIA.h"

class ReductionGIANodeInfo
{
public: 
    /**
      * constructor
      */ 
	ReductionGIANodeInfo(size_t node_id,
		                 size_t n_comp, 
	                     size_t n_eta,
						 size_t n_eta_bar,
						 size_t n_xi_global,
						 size_t n_xi_local,
						 ogsChem::chemReductionGIA* ReductionGIA);
    /**
      * destructor
      */ 
	~ReductionGIANodeInfo();

    /**
      * set the value of particular component
      */ 
	void set_comp_conc( size_t comp_idx, double val); 

    /**
      * return the node ID
      */ 
	size_t get_node_id(void) { return _node_id; }

    /**
      * return the value in eta mobile vector
      */ 
	double get_eta_value( size_t eta_idx ) { return loc_eta(eta_idx); }

    /**
      * return the value in eta immobile vector
      */ 
  //  double get_eta_bar_value( size_t eta_bar_idx ) { return loc_eta_bar(eta_bar_idx); }  //never called

    /**
      * return the value in xi global vector
      */ 
	double get_xi_global_value( size_t xi_global_idx ) { return loc_xi_global(xi_global_idx); }

    /**
      * return the value in xi local vector
      */ 
//	double get_xi_local_value( size_t xi_local_idx ) { return loc_xi_local(xi_local_idx); }  //never called

    /**
      * transform from concentrations to eta and xi values
      */ 
	void transform(void); 

private: 
    /**
      * node ID
      */ 
	const size_t _node_id; 

    /**
      * number of components
      */ 
	const size_t _n_comp; 

    /**
      * number of eta_mob
      */ 
	const size_t _n_eta;

    /**
      * number of eta_immob
      */ 
	const size_t _n_eta_bar;

    /**
      * number of xi_mob
      */
	const size_t _n_xi_global;

    /**
      * number of xi_immob
      */ 
	const size_t _n_xi_local;

    /**
      * component concentrations vector
      */ 
	MathLib::LocalVector _Comp_Conc; 

    /**
      * eta vector
      */ 
	MathLib::LocalVector loc_eta;

    /**
      * eta_bar vector
      */ 
	MathLib::LocalVector loc_eta_bar;

    /**
      * xi_global vector
      */ 
	MathLib::LocalVector loc_xi_global;

    /**
      * xi_local vector
      */ 
	MathLib::LocalVector loc_xi_local;

    /**
      * pointer to the ReductionGIA class
      */ 
	ogsChem::chemReductionGIA* _ReductionGIA;
}; 

#endif  // end of ifndef
