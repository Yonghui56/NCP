/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemconst.h
 *
 * Created on 2012-08-23 by Haibing Shao
 */

#ifndef CHEMCONST_H
#define CHEMCONST_H

#include <string>
#include <vector>
#include "MathLib/DataType.h"

namespace ogsChem
{
    // chemical and physical constants used in namespace ChemSys
    /**
      * Kelvin for zero degree Celcius
      */
    const double TK0 = 273.15;

    /**
      * Ideal gas constant, unit in J/mol/K
      */
    const double IDEAL_GAS_CONST = 8.3144621;

    /**
      * constant ln(10)
      */
    const double LN10 = 2.3025850929940456840179914546844;

    // enum structures
    // type of component
    enum Comp_Type {
         BASIS_COMP,      /* basis component                */
         AQ_PHASE_COMP,   /* aquous solution                */
         GAS_PHASE_COMP,  /* ideal gas mixture              */
         MIN_PHASE_COMP,  /* singal component mineral phase */
         SORPTION_COMP,   /* sorption component             */ 
         SS_PHASE_COMP,   /* solid solution component       */
		 KIN_COMP         /* kinetic component              */
    };

	enum Comp_Mobility {
		MOBILE,           /* mobile components*/
		IMMOBILE          /* immobile components*/
	}; 

    // type of chemical reaction
    enum ReactionType {
        EQ_LOGK,     /* - equilibrium reaction with constant logK */
        EQ_LOGK_T,   /* - equilibrium reaction with T dependent logK */
        KIN_SI_RATE, /* - kinetic reaction with saturation index based rate */
        KIN_ARB_RATE /* - kinetic reaction with arbitarily defined rate expression. */
    };

    // type of equilibrium reaction
    enum EqReactType {
        MOB_EQ_REACT,   /* - equilibrium reaction involves only mobile components */
        SORP_EQ_REACT,  /* - equilibrium reaction involves a sorption reaction */
        MIN_EQ_REACT    /* - equilibrium reaction involves a mineral reaction */
    };

	// type of kinetic reactions
	enum KinReactType {
		Monod,             /* - monod reaction */ 
		// DoubleMonodDecay,  /* - double monod reaction with decay term */ 
        MonodSum,          /* - sum of monod rate terms */
		UserExp,           /* - user defined rate expression */
		NoType             /* - uninitialized type */ 
	}; 

    // different thermodynamic databases
    enum DB_SOURCE {
        DB_UNKNOWN,
        THERMODDEM,   // THERMODDEM database from BRGM http://thermoddem.brgm.fr/
        NAGRA_PSI_01  // NAGRA-PSI database 01/01      http://les.web.psi.ch/TDBbook/index.htm
    };

    // format of the thermodynamic database file
    enum DB_FORMAT {
        PQC,          // PhreeQC format
        GRE_XML       // XML database format of this program
    };

    // the activity model included
    enum ACTIVITY_MODEL {
        ACT_MOD_AQ_DBH,  // Debyle-Huekel model for aqueous phase
        ACT_MOD_AQ_DVS,  // Davis model for aqueous phase
        ACT_MOD_UNITY    // activity coefficients are all ones
    };

	// type definitions
	typedef MathLib::LocalMatrix LocalMatrix;
	typedef MathLib::LocalVector LocalVector;
}

#endif // CHEMCONST_H
