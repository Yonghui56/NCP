/**
* Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.com/LICENSE.txt
*
*
* \file Local_ODE_Xi_immob_GIA.h
*
* Created on 2014-Jan-19 by Haibing Shao
*/

#ifndef LOCAL_ODE_KIN_REACT_H
#define LOCAL_ODE_KIN_REACT_H

#include "ChemLib/chemcomp.h"
#include "ChemLib/chemReactionKin.h"

class Local_ODE_KinReact
{
public:

    /**
    * constructor
    */
    Local_ODE_KinReact(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp, 
                       std::vector<ogsChem::chemReactionKin*> & list_kin_reactions, 
                       ogsChem::LocalMatrix & mat_Stoi_kin)
                       : _list_kin_reactions(list_kin_reactions)
    {
        _I = map_chemComp.size(); 
        _J_kin = _list_kin_reactions.size();
        
        _vec_conc_rates      = ogsChem::LocalVector::Zero(_I); 
        _vec_kin_react_rates = ogsChem::LocalVector::Zero(_J_kin);

        _Stoi_kin = mat_Stoi_kin;

    }

    /**
      * destructor
      */
    ~Local_ODE_KinReact()
    {
    }

    /**
      * evaluate the change of xi_immob over time
      */
    ogsChem::LocalVector operator() (double dt, MathLib::LocalVector & vec_Comp_Conc)
    {
        for (std::size_t i = 0; i < _J_kin; i++)
        {
            _list_kin_reactions[i]->calcReactionRate(vec_Comp_Conc);
            _vec_kin_react_rates(i) = _list_kin_reactions[i]->getRate();
        }
        _vec_conc_rates = _Stoi_kin * _vec_kin_react_rates; 
        return _vec_conc_rates;
    }

private:

    /**
      * vector storing the kinetic reaction rates
      */
    ogsChem::LocalVector _vec_conc_rates; 

    /**
      * vector storing the kinetic reaction rates
      */
    ogsChem::LocalVector _vec_kin_react_rates;

    /**
      * stoichiometric matrix of kinetic reactions 
      */
    ogsChem::LocalMatrix _Stoi_kin; 

    /**
      * address of kinetic reactions
      */
    std::vector<ogsChem::chemReactionKin*> & _list_kin_reactions; 

    /**
      * number of kinetic reactions
      */
    std::size_t _J_kin, _I; 
};

#endif
