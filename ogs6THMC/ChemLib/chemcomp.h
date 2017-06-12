/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemcomp.h
 *
 * Created on 2012-08-23 by Haibing Shao
 */

#ifndef CHEMCOMP_H
#define CHEMCOMP_H

#include "chemconst.h"
#include "ogs6THMC/MaterialLib/Compound.h"
#include "BaseLib/FileTools.h"

namespace ogsChem
{

// decalre ChemPhase class here to avoid cross inclusion.
class ChemPhase;
class ChemReaction;

class ChemComp
{
public:
    // constructor and destructor
    /**
      * constructor of the class
      */
    ChemComp(void);
	/**
      * constructor of the class, using the Compound class
      */
	ChemComp(MaterialLib::Compound* mCompound);

    /**
      * distructor of the class
      */
    ~ChemComp(void);

    /**
      * getter for mole mass in gram/mol.
      */
    double get_mole_mass(void) const;

    /**
      * getter for mole volume in m3.
      */
    double get_mole_volume(void) const;

    /**
      * getter for charge in the unit of equivalent.
      */
    double get_charge(void) const;

    double get_neg_gamma_1(void)
    {return neg_gamma_1;}

    double get_neg_gamma_2(void)
    {return neg_gamma_2;}

    /**
      * getter for chemical component name.
      */
    std::string& get_name(void)
    {return _name;}

    void set_name(std::string new_name)
    {_name = new_name;}

    void set_composition(std::string new_composition)
    {_composition = new_composition;}

	void set_mobility(Comp_Mobility new_mobility)
	{_Mobility = new_mobility; }

    void set_charge(double charge)
    {_charge = charge; }

    std::string get_composition(void)
    {return _composition;}

    void set_compTyps(Comp_Type new_type)
    {_mCompType = new_type;}

    void set_neg_gamma(double n_gamma_1, double n_gamma_2)
    {neg_gamma_1 = n_gamma_1; neg_gamma_2 = n_gamma_2;}

    void set_pReaction(ChemReaction* new_Reaction)
    {_mReaction = new_Reaction;}

    ChemReaction* get_pReaction(void)
    {return _mReaction;}

	size_t getIndex(void) {return _index;};
	
	void   setIndex(size_t new_index) {_index = new_index;};

	Comp_Mobility getMobility(void)
	{return _Mobility; }

    Comp_Type getCompType(void)
    {return _mCompType; }

    double getCharge(void)
    {return _charge; }

private:
    /**
      * The name of the chemical component.
      */
    std::string _name;

    /**
      * The molecular composition of the component.
      */
    std::string _composition;

    /**
      * Molar mass of the component, unit in g/mol.
      */
    double _mole_mass;

    /**
      * Molar volume of the component, unit in m3/mol.
      */
    double _mole_volume;

    /**
      * charge of the component, unit in eq/mol,
      * e.g. calsium cation Ca2+ has a charge of +2.0 eq/mol
      */
    double _charge;

    /**
      * Gibbs formation energy of the component.
      */
    double _G0;

    /**
      * pointer of phase this component belongs to.
      */
    ChemPhase* _mPhase;

    /**
      * pointer of phase this component belongs to.
      */
    ChemReaction* _mReaction;

    /**
      * type of phase this component belongs to.
      */
    Comp_Type _mCompType;

	/**
      * whether this component is mobile, sorption or mineral component
      */
	Comp_Mobility _Mobility; 

    /**
      * negative gamma values
      */
    double neg_gamma_1, neg_gamma_2;

	/**
      * index of the component
      */
    size_t _index;

};

} // end of namespace

#endif // CHEMCOMP_H
