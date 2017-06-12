/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ConcNodeInfo.h
 *
 * Created on 2013-03-31 by Haibing Shao
 */

#ifndef CONC_NODE_INFO_H
#define CONC_NODE_INFO_H

#include "ChemLib/chemEqReactSysActivity.h"

class ConcNodeInfo
{
public: 
    /**
      * constructor
      */ 
	ConcNodeInfo(size_t node_id, 
                 size_t n_comp, 
                 ogsChem::chemEqReactSysActivity* EqReactSys); 
    /**
      * destructor
      */ 
	~ConcNodeInfo(); 

    /**
      * set the value of particular component
      */ 
	void set_comp_conc( size_t comp_idx, double val);

    /**
      * get the value of particular component
      */ 
    double get_comp_conc( size_t comp_idx ) { return _Comp_Conc(comp_idx); };

    /**
      * return the node ID
      */ 
	size_t get_node_id(void) { return _node_id; }

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
      * component concentrations vector
      */ 
	MathLib::LocalVector _Comp_Conc; 

    /**
      * pointer to the equilibrium reaction class
      */ 
	ogsChem::chemEqReactSysActivity* _EqReactSys;
}; 

#endif  // end of ifndef