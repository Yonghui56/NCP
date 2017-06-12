/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemReactionEq.h
 *
 * Created on 2013-03-18 by Haibing Shao & Reza Zolfaghari
 */
 
#include "chemReactionEq.h"
#include "logog.hpp"

namespace ogsChem
{

chemReactionEq::chemReactionEq()
{
	// set initial values
	_eq_const_K = 1.0000001;
	_ln_K = log(_eq_const_K);  

}

chemReactionEq::~chemReactionEq(void)
{
}

void chemReactionEq::readReactionStr(std::string & /*reaction_str*/)
{
	// TODO
	std::string tmp_str;
}


void chemReactionEq::readReactionKRC(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & list_chemComp, 
                                      ogs5::CKinReact* KRC_reaction)
{
	// get the list of components 
	_vecCompNames = KRC_reaction->reactionpartner; 

	// find the components based on component names
	size_t i; 
	std::string str_comp; 
	BaseLib::OrderedMap<std::string, ogsChem::ChemComp*>::iterator comp_iterator; 
	for (i=0; i<_vecCompNames.size(); i++)
	{
		str_comp = _vecCompNames[i];
		comp_iterator = list_chemComp.find(str_comp);

		if ( comp_iterator != list_chemComp.end() )
		{
			ogsChem::ChemComp* component = comp_iterator->second; 
			this->_vecComponents.push_back(component);		
		}  // end of if	
	}  // end of for i

    // copy the stoichiometric vector 
	_vecStoi = KRC_reaction->stochmet;

	// read the rate parameters
    if ( KRC_reaction->getType().find("EQ_REACT") != std::string::npos )
	{   // this is an equilibrium reaction
        if ( KRC_reaction->getType().find("SORP_EQ_REACT") == 0 )
            this->_eqReactType = ogsChem::SORP_EQ_REACT; 
        else if ( KRC_reaction->getType().find("MIN_EQ_REACT") == 0 )
            this->_eqReactType = ogsChem::MIN_EQ_REACT;
        else if ( KRC_reaction->getType().find("MOB_EQ_REACT") == 0 )
            this->_eqReactType = ogsChem::MOB_EQ_REACT; 

        // read the equilibrium constant K
        this->_eq_const_K = KRC_reaction->eq_const_k; 
        // convert to logK and save it
        this->_ln_K = log(KRC_reaction->eq_const_k);
	}  // end of if KRC_reaction

}


}  // end of namespace
