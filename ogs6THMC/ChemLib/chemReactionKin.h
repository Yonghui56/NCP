/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemReaction.h
 *
 * Created on 2012-09-05 by Haibing Shao
 */
#include "chemReaction.h"
#include "ogsFileIO/FemIO/ogs5/rf_kinreact.h"
#include "BaseLib/OrderedMap.h" 
#include "muParser.h"

#ifndef CHEM_REACTION_KIN_H
#define CHEM_REACTION_KIN_H

namespace ogsChem
{

class chemReactionKin : public chemReaction
{
public: 
	/**
      * constructor and destructor
      */
	chemReactionKin();
	virtual ~chemReactionKin(void);

	/**
      * force the program to calculate the reaction rate
      */
	void calcReactionRate(ogsChem::LocalVector & vec_Comp_Conc); 

	/**
      * calculate the reaction rate with monod kinetics
      */
	double calcReactionRateMonod(ogsChem::LocalVector & vec_Comp_Conc); 

    /**
      * calculate the reaction rate with sum of monod kinetic rates
      */
	double calcReactionRateMonodSum(ogsChem::LocalVector & vec_Comp_Conc); 

    /**
      * calculate the reaction rate with user defined rate expression
      */
    double calcReactionRateUserExp(ogsChem::LocalVector & vec_Comp_Conc);

	/**
      * return the rate of current reaction. 
      */
	double getRate(void) { return _rate; }

	/**
      * override reading str function
      */
	void readReactionStr(std::string & reaction_str); 

	/**
      * reading one reaction from KRC data structure
      */
	void readReactionKRC(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & list_chemComp, 
                         ogs5::CKinReact* KRC_reaction); 

private:
	/**
      * type of kinetic reaction
      */
	KinReactType _kinReactType; 
	
	/**
      * the rate of reaction
      */
	double _rate; 

	/**
      * the rate constant
      */
	double _rate_constant; 

	/**
      * the rate constant order
      */
	double _rate_constant_order; 

    /**
      * index of the bacteria
      */
	size_t _idx_bacteria; 

    /**
      * the decay rate
      */
    double _decay_rate; 

	/**
	  * string defined by the user regarding the kinetic rate expression
	  */
	std::string _user_rate_Exp; 

	/**
	  * parser library handler
	  */
	mu::Parser* _userExp_parser;

    /**
      * flag whether the parser library has bee successfully initialized
      */
    bool _flag_parser_initialized;

	/**
      * a vector of monod rate components
      */
	std::vector<size_t> _vec_Monod_Comps_Idx; 

    /**
      * a vector of monod components' concentrations
      */
	std::vector<double> _vec_Monod_Comps_Conc; 

	/**
      * a vector of monod components' order
      */
	std::vector<double> _vec_Monod_Comps_order; 

    /**
      * a vector of inhibition rate components
      */
	std::vector<size_t> _vec_Inhibition_Comps_Idx; 

    /**
      * a vector of monod components' order
      */
	std::vector<double> _vec_Monod_Rate_Constants; 

    /**
      * a vector of inhibition components' concentrations
      */
	std::vector<double> _vec_Inhibition_Comps_Conc; 

    /**
      * a vector of inhibition components' order
      */
	std::vector<double> _vec_Inhibition_Comps_order; 

    /**
      * a vector of local concentration values.
      * this is kept here so that the user defined reactions will
      * find the corresponding component concentration at the same location. 
      */
    ogsChem::LocalVector  _vec_loc_Comp_Conc; 
};

}  // end of namespace

#endif
