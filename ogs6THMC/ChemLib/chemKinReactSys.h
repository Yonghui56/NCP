/**
* Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.com)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.com/LICENSE.txt
*
*
* \file chemKinReactSys.h
*
* Created on 2014-01-19 by Haibing Shao
*/
#ifndef CHEM_KIN_REACT_SYS_H
#define CHEM_KIN_REACT_SYS_H

#include "chemconst.h"
#include "chemcomp.h"
#include "chemReactionKin.h"
#include "BaseLib/OrderedMap.h"
#include "chemActivityModelUnity.h"
#include "UserModules/FemReactOPS/Local_ODE_KinReact.h"
#include "UserModules/FemReactGIAReduct/StepperBulischStoer.h"

namespace ogsChem
{

class chemKinReactSys
{
public:
    /**
    * constructor of the class
    */
    chemKinReactSys(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp,
                    std::vector<ogsChem::chemReactionKin*>               & list_kin_reactions,
                    ogsChem::chemActivityModelAbstract                  *a)
                    : _list_kin_reactions(list_kin_reactions), 
                      _I(map_chemComp.size()), _J_kin(list_kin_reactions.size()),
                      _activity_model(a)
    {
        double t0 = 0.0;

		countComp(map_chemComp); 

        buildStoi(map_chemComp, list_kin_reactions);

        _p_local_ODE = new Local_ODE_KinReact(map_chemComp, list_kin_reactions, _matStoi_kin);

        _vec_Comp_Conc = ogsChem::LocalVector::Zero(_I); 

        _vec_conc_rate = ogsChem::LocalVector::Zero(_J_kin);

        _sbs = new MathLib::StepperBulischStoer<Local_ODE_KinReact>(_vec_Comp_Conc, _vec_conc_rate,
                                                                    t0, 1.0e-12, 1.0e-12, true);

        if (_I > 0 && _J_kin > 0 && _p_local_ODE && _sbs)
            isInitialized = true; 

    };

    /**
      * destructor of the class
      */
    ~chemKinReactSys(void)
    {
        BaseLib::releaseObject(_p_local_ODE);
        BaseLib::releaseObject(_sbs);
    };

    /**
      * whether the reduction scheme has been initialized
      */
    bool IsInitialized(void) { return isInitialized; };

    /**
      * solve the kinetic reaction system
      * using the StepperBulischStoer ODE solver. 
      * input values are: 
      * the vector of concentration values
      * the convergence tolerance,
      * the maximum number of iterations.
      * and the status of convergence
      */
    void solve_KinSys(LocalVector & vec_conc,
                      size_t & result,
                      size_t & node_idx,
                      double iter_tol,
                      double rel_tol,
                      double dt, 
                      double max_iter)
    {
        // getting the initial rate evaluation. 
        if (_vec_conc_rate.rows() > 0)
            _vec_conc_rate = (*_p_local_ODE)(dt, vec_conc);

        //_sbs->set_y(Xi_Kin_bar);
        _sbs->set_y(vec_conc);
        _sbs->set_dydx(_vec_conc_rate);

        // solve the local ODE problem for xi kin bar
        _sbs->step(dt, _p_local_ODE);
        //vec_Xi_Kin_bar_new = _sbs->get_y();
        vec_conc = _sbs->get_y();
		// successfully run. 
		result = 0; 
    }; 


	std::size_t get_n_Kin_React(void) { return _J_kin;  };

	std::size_t get_n_Comp(void) { return _I; }; 

	std::size_t get_n_Comp_mob(void) { return _I_mob + _I_sec_mob;  }

private:

	void countComp(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp)
	{
		_I_mob = 0;
		_I_sec_mob = 0;
		_I_sec_sorp = 0;
		_I_sec_min = 0;
		_I_kin_mob = 0; 
		_I_kin_immob = 0; 

		_I = map_chemComp.size(); 
		BaseLib::OrderedMap<std::string, ogsChem::ChemComp*>::iterator it;
		for (it = map_chemComp.begin(); it != map_chemComp.end(); it++)
		{
			switch (it->second->getCompType())
			{
			case ogsChem::AQ_PHASE_COMP:
				_I_mob++;
				break;
			case ogsChem::SORPTION_COMP:
				_I_sec_sorp++;
				break;
			case ogsChem::MIN_PHASE_COMP:
				_I_sec_min++;
				break;
			case ogsChem::KIN_COMP: 
				if (it->second->getMobility() == ogsChem::MOBILE)
					_I_kin_mob++;
				else
					_I_kin_immob++; 
				break; 
			default:
				_I_sec_min++;
				break;
			}
		}

		assert(_I == _I_mob + _I_sec_sorp + _I_sec_min + _I_kin_mob + _I_kin_immob );
	};

    void buildStoi(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp, 
                   std::vector<ogsChem::chemReactionKin*>               & list_kin_reactions)
    {
        size_t i, j, tmp_idx;
        double tmp_stoi;
        std::string tmp_str;
        BaseLib::OrderedMap<std::string, ogsChem::ChemComp*>::iterator tmp_Comp;

        // size info
        _I = map_chemComp.size(); 
        _J_kin = list_kin_reactions.size();

        // creat the memory for Stoi matrix
        _matStoi_kin = LocalMatrix::Zero(_I, _J_kin);

        // based on the reactions, fill in the stoi matrix
        // loop over all kinetic reactions
        for (j = 0; j < list_kin_reactions.size(); j++)
        {	// for each reaction
            // find each participating components
            for (i = 0; i < list_kin_reactions[j]->get_vecCompNames().size(); i++){
                tmp_str = list_kin_reactions[j]->get_vecCompNames()[i];
                tmp_Comp = map_chemComp.find(tmp_str);
				// TODO: 
				// at this location, it is very likely that the user will type 
				// in some illegal text in the equation. Then it will not be found 
				// in the map of chemical components. 
				// here, we shall use the std::map in the future and catch the exception

                tmp_idx = tmp_Comp->second->getIndex();

				// for kinetic reactions, we do not differentiate 
				// basis species. Therefore treat as normal reactions
                tmp_stoi = list_kin_reactions[j]->get_vecStoi()[i];
                // and put them into Stoi matrix
                _matStoi_kin(tmp_idx, j) = tmp_stoi;

            }  // end of for i
        }  // end of for j

#ifdef _DEBUG
        // debugging--------------------------
        // std::cout << "Stoichiometric Matrix S: " << std::endl;
        // std::cout << _matStoi_kin << std::endl;
        // end of debugging-------------------
#endif
    }


    /**
      * private flag indicating initialization
      */
    bool isInitialized;

    /**
      * a list of all kinetic reactions
      */
    std::vector<ogsChem::chemReactionKin*> & _list_kin_reactions;

    /**
      * pointer to the activity model
      */
    ogsChem::chemActivityModelAbstract *_activity_model;

    /**
      * _I is the number of components and _J_kin is the number of kinetic reactions
      */
    size_t _I, _I_mob, _J_kin;

	size_t _I_sec_mob, _I_sec_sorp, _I_sec_min, _I_kin_mob, _I_kin_immob; 

    /**
      * pointer to the Local_ODE_KinReact
      */
    Local_ODE_KinReact* _p_local_ODE; 

    /**
      * pointer to the ODE solver StepperBulischStoer
      */
    MathLib::StepperBulischStoer<Local_ODE_KinReact>* _sbs;

    ogsChem::LocalVector _vec_Comp_Conc; 

    ogsChem::LocalVector _vec_conc_rate;

    ogsChem::LocalMatrix _matStoi_kin; 

};

}

#endif