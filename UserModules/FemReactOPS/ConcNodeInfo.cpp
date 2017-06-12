/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ConcNodeInfo.cpp
 *
 * Created on 2013-03-31 by Haibing Shao
 */
 
#include "ConcNodeInfo.h"
 
ConcNodeInfo::ConcNodeInfo(size_t node_id, 
                           size_t n_comp,
                           ogsChem::chemEqReactSysActivity* EqReactSys)
    : _node_id(node_id), _n_comp(n_comp), _EqReactSys(EqReactSys)
{
	_Comp_Conc = MathLib::LocalVector::Zero( _n_comp ); 
}

ConcNodeInfo::~ConcNodeInfo()
{
    _EqReactSys = NULL;
}

void ConcNodeInfo::set_comp_conc( size_t comp_idx, double val )
{
	_Comp_Conc(comp_idx) = val; 
}

