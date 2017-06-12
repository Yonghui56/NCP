/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemReductionGIA.h
 *
 * Created on 2013-05-02 by Reza Zolfaghari and Haibing Shao
 */
#ifndef CHEMREDUCTIONGIA_H
#define CHEMREDUCTIONGIA_H

#include "chemconst.h"
#include "chemReactionKin.h"
#include "chemReactionEq.h"
#include "chemcomp.h"
#include "BaseLib/OrderedMap.h"
#include "chemActivityModelUnity.h"

namespace ogsChem
{

class chemReductionGIA
{
public:
	/**
      * constructor of the class
      */
	chemReductionGIA(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & list_chemComp,
		             std::vector<ogsChem::chemReactionEq*>                & list_eq_reactions,
		             std::vector<ogsChem::chemReactionKin*>               & list_kin_reactions,
		             ogsChem::chemActivityModelAbstract                  *a );
	
	/**
      * destructor of the class
      */
	~chemReductionGIA(void){};

    /**
      * convert concentration vector to eta and xi vector

	void Conc2EtaXi(ogsChem::LocalVector &local_conc,
					ogsChem::LocalVector &local_eta,
					ogsChem::LocalVector &local_eta_bar,
					ogsChem::LocalVector &local_xi_Mob,
					ogsChem::LocalVector &local_xi_Sorp_tilde,
					ogsChem::LocalVector &local_xi_Sorp_bar,
					ogsChem::LocalVector &local_xi_Min_tilde,
					ogsChem::LocalVector &local_xi_Min_bar,
					ogsChem::LocalVector &local_xi_Kin,
					ogsChem::LocalVector &local_xi_Kin_bar);
*/
	void Conc2EtaXi(ogsChem::LocalVector &local_conc,
					ogsChem::LocalVector &local_eta,
					ogsChem::LocalVector &local_eta_bar,
				    ogsChem::LocalVector &local_xi_global,
					ogsChem::LocalVector &local_xi_local);

	/**
      * convert eta and xi vectors to concentration vector
 */
	void EtaXi2Conc(ogsChem::LocalVector &local_eta,
					ogsChem::LocalVector &local_eta_bar,
					ogsChem::LocalVector &local_xi_global,
					ogsChem::LocalVector &local_xi_local,
					ogsChem::LocalVector &local_conc);

	void EtaXi2Conc_JH(ogsChem::LocalVector &local_eta,
					ogsChem::LocalVector &local_eta_bar,
					ogsChem::LocalVector &local_xi_Mob,
					ogsChem::LocalVector &local_xi_Sorp_tilde,
					ogsChem::LocalVector &local_xi_Sorp_bar,
					ogsChem::LocalVector &local_xi_Min_tilde,
					ogsChem::LocalVector &local_xi_Min_bar,
					ogsChem::LocalVector &local_xi_Kin,
					ogsChem::LocalVector &local_xi_Kin_bar,
					ogsChem::LocalVector &local_conc);

	void EtaXi2Conc_JH_NOCUTOFF(ogsChem::LocalVector &local_eta,
					ogsChem::LocalVector &local_eta_bar,
					ogsChem::LocalVector &local_xi_Mob,
					ogsChem::LocalVector &local_xi_Sorp_tilde,
					ogsChem::LocalVector &local_xi_Sorp_bar,
					ogsChem::LocalVector &local_xi_Min_tilde,
					ogsChem::LocalVector &local_xi_Min_bar,
					ogsChem::LocalVector &local_xi_Kin,
					ogsChem::LocalVector &local_xi_Kin_bar,
					ogsChem::LocalVector &local_conc);

    /**
      * whether the reduction scheme has been initialized
      */
	bool IsInitialized(void) {return isInitialized;}; 

	/**
	 * calculate the reaction rates
	 */
	void Calc_Kin_Rate(ogsChem::LocalVector &vec_conc_linear,
					   ogsChem::LocalVector &vec_rates);

	void Calc_Kin_Rate_temp(ogsChem::LocalVector &local_xi_Mob,
										 ogsChem::LocalVector &local_xi_Sorp,
				                         ogsChem::LocalVector &local_xi_Sorp_tilde,
										 ogsChem::LocalVector &local_xi_Sorp_bar,
										 ogsChem::LocalVector &local_xi_Min,
										 ogsChem::LocalVector &local_xi_Min_tilde,
										 ogsChem::LocalVector &local_xi_Min_bar,
										 ogsChem::LocalVector &local_xi_Kin,
									     ogsChem::LocalVector &local_xi_Kin_bar,
										 ogsChem::LocalVector &local_eta,
										 ogsChem::LocalVector &local_eta_bar,
										 ogsChem::LocalVector &vec_rates);

	/**
      * get the number of components
      */
	std::vector<ogsChem::chemReactionKin*> & get_list_kin_reactions(void) {return _list_kin_reactions; };

	/**
      * get the number of components
      */
	size_t get_n_Comp(void) {return _I_tot; };

	/**
      * get the number of mobile components
      */
	size_t get_n_Comp_mob(void) {return _I_mob; };

	/**
      * get the number of immobile sorbed components
      */
	size_t get_n_Comp_sorb(void) {return _I_sorp; };
	
	/**
      * get the number of immobile mieral components
      */
	size_t get_n_Comp_min(void) {return _I_min; };

	/**
	  * get the number of kinetic components
	  */
	size_t get_n_Comp_kin(void) {return _I_kin; };

    /**
      * get the number of mobile reactions
      */
    size_t get_J_mob(void) { return _Jmob; };

    /**
      * get the number of sorption reactions
      */
    size_t get_J_sorp(void) { return _Jsorp; };

    /**
      * get the number of mineral reactions
      */
    size_t get_J_min(void) { return _Jmin; };

	/**
      * get the length of eta_mob
      */
	size_t get_n_eta(void) {return _n_eta;};

	/**
      * get the length of eta_immob
      */
	size_t get_n_eta_bar(void) {return _n_eta_bar;};

	/**
      * get the length of eta_global
      */
	size_t get_n_xi_global(void) {return _n_xi_global;};

	/**
      * get the length of eta_local
      */
	size_t get_n_xi_local(void) {return _n_xi_local;};


	/**
      * get the length of eta_mob
      */
	size_t get_n_xi_mobile(void) {return _n_xi_mobile;};

	/**
      * get the length of eta_mob
      */
	size_t get_n_xi_immobile(void) {return _n_xi_immobile;};

	/**
      * get the length of xi_Mob = Jmob
      */
	size_t get_n_xi_Mob(void) {return _n_xi_Mob; };

	/**
      * get the length of xi_Sorp_tilde = Jsorpli
      */
	size_t get_n_xi_Sorp(void) {return _n_xi_Sorp; };

	/**
      * get the length of xi_Min_bar = Jmin
      */
	size_t get_n_xi_Min(void) {return _n_xi_Min; };


	/**
      * get the length of xi_Sorp_tilde = Jsorpli
      */
	size_t get_n_xi_Sorp_tilde(void) {return _n_xi_Sorp_tilde; };

	/**
      * get the length of xi_Sorp_bar = Jsorp
      */
	size_t get_n_xi_Sorp_bar(void) {return _n_xi_Sorp_bar; };

	/**
      * get the length of xi_Sorp_bar_li = Jsorp_li
      */
	size_t get_n_xi_Sorp_bar_li(void) {return _n_xi_Sorp_bar_li; };

	/**
      * get the length of xi_Sorp_bar_ld = Jsorp_ld
      */
	size_t get_n_xi_Sorp_bar_ld(void) {return _n_xi_Sorp_bar_ld; };

	/**
      * get the length of xi_Min_tilde = Jmin
      */
	size_t get_n_xi_Min_tilde(void) {return _n_xi_Min_tilde; };

	/**
      * get the length of xi_Min_bar = Jmin
      */
	size_t get_n_xi_Min_bar(void) {return _n_xi_Min_bar; };
	/**
      * get the length of xi_ki = j1kin
      */
	size_t get_n_xi_Kin(void) {return _n_xi_Kin; };

	/**
      * get the length of xi_Kin_bar = j2kin
      */
	size_t get_n_xi_Kin_bar(void) {return _n_xi_Kin_bar; };

	/**
      * get the length of xi_ki total
      */
	size_t get_n_xi_Kin_total(void) {return _J_tot_kin; };

	/**
	  * get _J_2_kin_ast
	  */
	size_t get_J_2_kin_ast(void) { return _J_2_kin_ast; };

	/**
      * get the _mat_c_mob_2_xi_mob
      */
	LocalMatrix get_matrix_C2Xi(void) {return _mat_c_mob_2_xi_mob; };

	/**
      * get the _mat_c_immob_2_xi_immob
      */
	LocalMatrix get_matrix_Cbar2XiBar(void) {return _mat_c_immob_2_xi_immob; };

	/**
      * get the _mat_c_mob_2_eta
      */
	LocalMatrix get_matrix_C2Eta(void) {return _mat_c_mob_2_eta_mob; };

	/**
      * get the _mat_c_mob_2_eta bar
      */
	LocalMatrix get_matrix_C2EtaBar(void) {return _mat_c_immob_2_eta_immob; };

	/**
      * get the _mat_Ald matrix
      */
	LocalMatrix get_matrix_Ald(void) {return _mat_Ald; };

	/**
      * get the _mat_A1sorp matrix
      */
	LocalMatrix get_matrix_A1sorp(void) {return _mat_A1sorp; };

	/**
      * get the _mat_A2sorpli matrix
      */
	LocalMatrix get_matrix_A2sorpli(void) {return _mat_A2sorpli; };

	/**
      * get the _mat_A2sorpld matrix
      */
	LocalMatrix get_matrix_A2sorpld(void) {return _mat_A2sorpld; };

	/**
      * get the _mat_A1min matrix
      */
	LocalMatrix get_matrix_A1min(void) {return _mat_A1min; };

	/**
      * get the _mat_A2kin matrix
      */
	LocalMatrix get_matrix_A2kin(void) {return _mat_A2kin; };

	/**
      * get the _mat_A1kin matrix
      */
	LocalMatrix get_matrix_A1kin(void) {return _mat_A1kin; };

	/**
      * get the _mat_S1min matrix
      */
	LocalMatrix get_matrix_S1min(void) {return _mat_S1min; };

	/**
      * get the _mat_S1mob matrix
      */
	LocalMatrix get_matrix_S1mob(void) {return _mat_S1mob; };

	/**
      * get the _mat_S1sorp matrix
      */
	LocalMatrix get_matrix_S1sorp(void) {return _mat_S1sorp; };

	/**
      * get the _mat_S1sorpli matrix
      */
	LocalMatrix get_matrix_S1sorpli(void) {return _mat_S1sorp_li; };

	/**
      * get the _mat_S1kin_ast matrix
      */
	LocalMatrix get_matrix_S1kin_ast(void) {return _mat_S1kin_ast; };

	/**
      * get the _mat_S2sorp matrix
      */
	LocalMatrix get_matrix_S2sorp(void) {return _mat_S2sorp; };


	/**
      * get the _matrix_Ssorp matrix
      */
	LocalMatrix get_matrix_Ssorp(void) {return _matrix_Ssorp; };

	/**
      * get the log k of equilibrium mobile reactions
      */
	LocalMatrix get_logk_mob(void) {return _logk_mob; };
	/**
      * get the log k of equilibrium sorbed reactions
      */
	LocalMatrix get_logk_sorp(void) {return _logk_sorp; };
	/**
      * get the log k of equilibrium mineral reactions
      */
	LocalMatrix get_logk_min(void) {return _logk_min; };

	/**
      * get _mat_S1_ast matrix
      */
	LocalMatrix get_mat_S1_ast(void) {return _mat_S1_ast; };
	/**
      * get _mat_S2_ast matrix
      */
	LocalMatrix get_mat_S2_ast(void) {return _mat_S2_ast; };
	/**
      * get _mat_S1_orth
      */
	LocalMatrix get_mat_S1_orth(void) {return _mat_S1_orth; };
	/**
      * get _mat_S2_orth
      */
	LocalMatrix get_mat_S2_orth(void) {return _mat_S2_orth; };

	/**
      * get _activity_model
      */
	ogsChem::chemActivityModelAbstract * get_activity_model(void) {return _activity_model; };



private:
	/**
      * private flag indicating initialization
      */
	bool isInitialized; 


	/**
      * a list of all equilibrium reactions
      */
	std::vector<ogsChem::chemReactionEq*> & _list_eq_reactions;

	/**
      * a list of all kinetic / equilibrium reactions
      */
	std::vector<ogsChem::chemReactionKin*> & _list_kin_reactions;

    /**
      * pointer to the activity model
      */
    ogsChem::chemActivityModelAbstract *_activity_model;

    /**
      * stoichiometric matrix S
      */	
	LocalMatrix _matStoi; 

	/**
      * sub-matrixes of S, 
      */	
	LocalMatrix _mat_S1, _mat_S1mob,_mat_S2mob, _mat_Ssorp,_mat_Ssorp_ast,_mat_Ssorp_li,_mat_Ssorp_ld, _mat_S1sorp_li, _mat_S1sorp_ld, _mat_S1sorp, _mat_S1min, _mat_S2min, _mat_S2, _mat_S2sorp, _mat_S1kin, _mat_S2kin,
	            _mat_Skin,_mat_Skin_li, _mat_S2sorp_li,_mat_S2sorp_ld,_mat_S1_preserve, _mat_S2_preserve, _matrix_Ssorp;

//	/**
//      * LocalVector sub-vectors of xi and eta,    ?
//      */
//	LocalVector local_xi_Mob, local_xi_Sorp,local_xi_Sorp_li,local_xi_Sorp_ld,local_xi_Sorp_tilde,local_xi_Sorp_bar,local_xi_Min,local_xi_Min_tilde,local_xi_Min_bar,local_xi_Kin,local_xi_Kin_bar,
//	            local_xi,local_xi_bar;

	LocalVector local_conc,local_eta,local_eta_bar;


	/**
      * Max linearly independent columns (astride/ast) of the matrixes
      */
	LocalMatrix _mat_S1kin_ast, _mat_S2kin_ast, _mat_S1_ast, _mat_S2_ast;

	/**
      * orthogonal complementary of the matrixes
      */
	LocalMatrix _mat_S1_orth, _mat_S2_orth, _mat_S2_orth_tilde;

	/**
      * matrix for the calculation of right hand side rate transformation
      */
	LocalMatrix _mat_A1,_mat_A2,_mat_A1sorp, _mat_A1min, _mat_A1kin, _mat_A2kin, _mat_Ald, _mat_A2sorp, _mat_A2sorpli, _mat_A2sorpld;

	/**
      * transformation matrix for the calculation between eta, xi and concentrations
      */
	LocalMatrix _mat_c_mob_2_eta_mob,
		        _mat_c_immob_2_eta_immob,   
	            _mat_c_mob_2_xi_mob, 
	            m,
				_mat_c_immob_2_xi_immob;

	/**
	 *
	 */
	LocalVector _logk_mob, _logk_sorp, _logk_min;
	/**
      * the size of eta and xi vector
      */
	size_t _n_eta,_n_eta_bar, _n_xi_mobile, _n_xi_immobile, _n_xi_local,  _n_xi_global, _n_xi_Mob, _n_xi_Sorp_tilde,_n_xi_Sorp,
	       _n_xi_Sorp_bar,_n_xi_Sorp_bar_li,_n_xi_Sorp_bar_ld,_n_xi_Min, _n_xi_Min_tilde, _n_xi_Min_bar, _n_xi_Kin, _n_xi_Kin_bar;

	/**
      * _I_tot is the number of all components and _J_tot is the number of all reactions in the system
      */
	//size_t _I_tot, _J_tot;

	/**
      * _I_tot is the number of all components and _J_tot is the number of all reactions in the system
      */
	//size_t 	Linear_Indep, Linear_Dep;
	/**
      * _I_tot is the number of all components and _J_tot is the number of all reactions in the system
      */
	size_t _J_tot_eq, _J_tot_kin,_J_tot_li,_J_tot;

	/**
      * number of mobile, sorption and mineral and kinetic reactions
      */
	size_t _Jeq_li, _Jmob, _Jsorp, _Jsorp_li, _Jsorp_ld, _Jmin, _J_1_kin_ast, _J_2_kin_ast;

	/**
      * number of mobile and immobile components and species
      */
	size_t _I_bar, _I_NMin_bar, _I_tot,_I_mob, _I_sorp,_I_min, _I_kin, _I_immob, _I_immob_min, _I_immob_nonMin;

    /**
     * will be moved to the main function
      * reading in the logK values of each equilibrium reactions
      */
    void read_logK(std::vector<ogsChem::chemReactionEq*>                & list_eq_reactions);


	/**
      * construct stoichiometric matrix out of list of components and reactions
      */
	void buildStoi(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & list_chemComp, 
		           std::vector<ogsChem::chemReactionEq*>                & list_eq_reactions,
		           std::vector<ogsChem::chemReactionKin*>               & list_kin_reactions);

	/**
      * calculate the intemediate parameters for the reduction scheme
      */
	void update_reductionScheme(void); 

	/**
      * return the orthogonal complement of the given matrix
      */
	LocalMatrix orthcomp( LocalMatrix & inMat ); 

	/**
      * count how many mobile and immobile components
      */
	void countComp(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp); 

};

}

#endif
