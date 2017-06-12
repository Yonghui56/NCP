/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemelem.h
 *
 * Created on 2012-08-23 by Haibing Shao
 */
#ifndef CHEMELEM_H
#define CHEMELEM_H

#include "chemconst.h"

namespace ogsChem
{
class ChemElem
{
public:
    /**
      * constructor of the class
      */
    ChemElem();

    /**
      * distructor of the class
      */
    ~ChemElem();

    /**
      * getter for element name.
      */
    std::string& get_name(void)
    {return _name;}

    /**
      * setter for element name.
      */
    void set_name(std::string new_name)
    {_name = new_name;}

    /**
      * getter for alkalinity.
      */
    double get_alk(void)
    {return _alkalinity;}

    /**
      * setter for alkalinity.
      */
    void set_alk(double new_alk)
    {_alkalinity = new_alk;}

    /**
      * getter for mole mass.
      */
    double get_mole_mass(void)
    {return _mole_mass;}

    /**
      * setter for mole mass.
      */
    void set_mole_mass(double new_mole_mass)
    {_mole_mass = new_mole_mass;}

    /**
      * getter for charge.
      */
    double get_charge(void)
    {return _charge;}

    /**
      * setter for charge.
      */
    void set_charge(double new_charge)
    {_charge = new_charge;}

    /**
      * getter for atm number.
      */
    size_t get_atom_num(void)
    {return _atm_num;}

    /**
      * setter for atm number.
      */
    void set_atom_num(size_t new_atm_num)
    {_atm_num = new_atm_num;}

private:
    /**
      * The name of the chemical element.
      */
    std::string _name;

    /**
      * The alkalinity of the element.
      */
    double _alkalinity;

    /**
      * Molar mass of the element, unit in g/mol.
      */
    double _mole_mass;

    /**
      * charge of the component, unit in eq/mol,
      * e.g. calsium cation Ca2+ has a charge of +2.0 eq/mol
      */
    double _charge;

    /**
      * atomic number of the element
      * e.g. the atomic number of element carbon is 6.
      */
    size_t _atm_num;

};
}

#endif // CHEMELEM_H
