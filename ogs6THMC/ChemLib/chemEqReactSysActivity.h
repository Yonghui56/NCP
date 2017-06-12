/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemEqReactSys.h
 *
 * Created on 2013-09-18 by Haibing Shao
 */
#ifndef CHEM_EQ_REACT_SYS_ACTIVITY_H
#define CHEM_EQ_REACT_SYS_ACTIVITY_H

#include "logog.hpp"
#include "chemconst.h"
#include "chemcomp.h"
#include "chemReactionEq.h"
#include "BaseLib/OrderedMap.h"
#include "chemActivityModelUnity.h"

namespace ogsChem
{

/**
  * The main difference btw. this class and the class chemEqReactSys is, 
  * we us [logC; c_min_bar] as the unknown vector. 
  * The size of this vector is equal to the number of components. 
  * by keeping the log, we will keep concentration values to positive
  */
class chemEqReactSysActivity
{
public:
    /**
      * constructor of the class
      */
	chemEqReactSysActivity(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp, 
                           std::vector<ogsChem::chemReactionEq*>                & list_eq_reactions, 
                           ogsChem::chemActivityModelAbstract                  *a )
      : _list_eq_reactions(list_eq_reactions), _I(0), _J(0), _activity_model(a)
    {
	    // by default, the class is not yet initialized
	    isInitialized = false; 

	    if ( map_chemComp.size() > 0 && _list_eq_reactions.size() > 0 )
	    {   // if there are reactions. 
		    // cout how many mobile and how many immobile components
		    countComp(map_chemComp); 
            // count how many which type of reactions we have
            countReactions( map_chemComp, _list_eq_reactions ); 
		    // make stoichiometric matrix
		    buildStoi(map_chemComp, list_eq_reactions);
            // read the logK values
            read_logK(list_eq_reactions); 
            // allocate memory for local residual vector, activity coefficients, and ln concentrations.
            _vec_res           = ogsChem::LocalVector::Zero( _I ); 
            _ln_activity_coeff = ogsChem::LocalVector::Ones( _I ); // by default set activity coeff to 1
            _ln_activity       = ogsChem::LocalVector::Zero( _I ); 
            // allocate memory for local Jacobi matrix
            _mat_Jacobi = ogsChem::LocalMatrix::Zero( _I, _I ); 
		    // allocate memory for vector AI
            if ( this->_I_sec_min > 0 )
                _AI = ogsChem::LocalVector::Zero(_I_sec_min); 
            // flip the initialization flag
            isInitialized = true; 
	    }
    };
    
    /**
      * destructor of the class
      */
    ~chemEqReactSysActivity(void){}; 
    /**
      * whether the reduction scheme has been initialized
      */
	bool IsInitialized(void) {return isInitialized;}; 
	
	/**
      * get the number of components
      */
	size_t get_n_Comp(void) {return _I; }; 

    /**
      * get the number of mobile components
      */
	size_t get_n_Comp_mob(void) {return _I_mob; }; 

	/**
	  * get the number of equilibrium reactions
	  */
	size_t get_n_Eq_React(void) { return _J; };

    /**
      * calculate the total mass using 
      * the concentration vector of basis, and
      * the concentration vector of secondary variables
      */
    void calc_tot_mass(LocalVector & vec_conc_basis, 
                       LocalVector & vec_conc_second, 
                       LocalVector & vec_tot_mass)
    {
        vec_tot_mass = vec_conc_basis - _matStoi.transpose() * vec_conc_second; 
    }; 

    /**
      * calculate the residual of the reaction system using 
      * the concentration vector of all components and
      * the total concentration constrain of the basis species
      */
    void calc_residual(LocalVector & vec_unknowns, 
                       LocalVector & vec_tot_mass_constrain, 
                       LocalVector & vec_residual)
    {
        size_t i; 
        double res_tmp, phi;
        ogsChem::LocalVector c_basis, c_sec_min, c_second; 
        ogsChem::LocalVector ln_activity_basis, ln_activity_sec_mob, ln_activity_sec_sorp; 
        ogsChem::LocalVector vec_ln_a; 
        ogsChem::LocalVector vec_cur_mass_balance;
        ogsChem::LocalVector lnK_min;
        ogsChem::LocalMatrix Stoi_mob, Stoi_sorp, Stoi_min; 

        vec_cur_mass_balance = ogsChem::LocalVector::Zero( _I_basis  ); 
        ln_activity_basis    = ogsChem::LocalVector::Zero( _I_basis );
        c_basis              = ogsChem::LocalVector::Zero( _I_basis  ); 
        c_second             = ogsChem::LocalVector::Zero( _I_second );
        _ln_activity         = ogsChem::LocalVector::Zero( _I        ); 

        c_sec_min  = vec_unknowns.tail(_I_sec_min); 

        // part 0), calculate the concentration of secondary 
        // non-mineral components
        Stoi_mob  = _matStoi.topRows(    _J_mob );
        Stoi_sorp = _matStoi.middleRows( _J_mob, _J_sorp ); 
        Stoi_min  = _matStoi.bottomRows( _J_min );
    
        // TODO, do the activity correction
        this->_activity_model->calc_activity_logC( vec_unknowns, 
                                                   _ln_activity_coeff, 
                                                   _ln_activity ); 
        // now split the unknown vector
        ln_activity_basis = _ln_activity.head(_I_basis);
        // calculate the secondary mobile component concentrations
        ln_activity_sec_mob  = _ln_activity.segment( _I_basis, _I_sec_mob ); 
        // calculate the secondary sorption component concentrations
        ln_activity_sec_sorp = _ln_activity.segment( _I_basis + _I_sec_mob, _I_sec_sorp ); 
        lnK_min = _vec_lnK.tail(_J_min); 
    
        // fill in the basis concentrations
        for (i=0; i < _I_basis ; i++ )
            c_basis(i) = std::exp( ln_activity_basis(i) - _ln_activity_coeff(i) ); 
        // fill in the secondary concentrations
        for (i=0; i < ( _I_sec_mob + _I_sec_sorp ); i++)
            c_second( i ) = std::exp( _ln_activity(_I_basis + i) - _ln_activity_coeff( _I_basis + i ) );
        c_second.tail( _I_sec_min )  = _ln_activity.tail( _I_sec_min );
        
        // start calculating the residuals. 
        // part 1), n_basis mass balance equations
        this->calc_tot_mass( c_basis, c_second, vec_cur_mass_balance ); 
        vec_residual.head( _I_basis ) = vec_tot_mass_constrain - vec_cur_mass_balance; 

        // part 2), the secondary equilibrium reactions, 
        vec_residual.segment( _I_basis, _I_sec_mob ) = ln_activity_sec_mob - _vec_lnK.head( _J_mob ) + Stoi_mob * ln_activity_basis; 
        vec_residual.segment( _I_basis + _I_sec_mob, _I_sec_sorp ) = ln_activity_sec_sorp - _vec_lnK.segment( _J_mob, _J_sorp ) + Stoi_sorp * ln_activity_basis;

        // part 3), n_react_min mineral reactions, 
        // AKA, the "complementary problem".
        for ( i=0; i < _J_min; i++ )
        {
            // attention, this is spectial for mineral reactions
            phi  = -1.0 * lnK_min(i) + Stoi_min.row(i) * ln_activity_basis;
            res_tmp  = std::min( phi, c_sec_min(i) ); 
            vec_residual(_I_basis + _I_sec_mob + _I_sec_sorp + i) = res_tmp; 
        }  // end of for
        

    };  // end of function calc_residual

    /**
      * calculate the Jacobi of the reaction system analytically using 
      * the concentration vector of all components, 
      * the total concentration constrain of the basis species, 
      * and an already calculated residual vector
      */
    void calc_Jacobi(LocalVector & vec_unknowns,
                     LocalVector & vec_tot_mass_constrain,
                     LocalVector & vec_res_base)
    {
        const double epsilon = 1.0e-8; 
        size_t i; 
        LocalVector vec_unknown_tmp; 
        LocalVector vec_res_tmp; 
        // Jacobi matrix is a square matrix, 
        // with n_cols and n_rows equal to the number of concentrations
        size_t n_unknowns(vec_unknowns.rows());
        vec_res_tmp = LocalVector::Zero( n_unknowns ); 
        _mat_Jacobi.setZero(); 

        // using the numerical increment method to construct Jacobi matrix 
        for (i=0; i<n_unknowns; i++)
        {
            // increment the vec_unknown
            vec_unknown_tmp = vec_unknowns; 
            if ( vec_unknown_tmp.norm() < std::numeric_limits<double>::epsilon() )
            {
                vec_unknown_tmp(i) = epsilon;
                this->calc_residual(vec_unknown_tmp, vec_tot_mass_constrain, vec_res_tmp); 
                _mat_Jacobi.col(i) = ( vec_res_tmp - vec_res_base ) / epsilon ;
            }
            else
            {
                vec_unknown_tmp(i) = vec_unknown_tmp(i) + epsilon * vec_unknown_tmp.norm();
                this->calc_residual(vec_unknown_tmp, vec_tot_mass_constrain, vec_res_tmp); 
                _mat_Jacobi.col(i) = ( vec_res_tmp - vec_res_base ) / (epsilon * vec_unknown_tmp.norm());
            }  // end of else
        }  // end of for loop
    #ifdef _DEBUG
	    // debugging--------------------------
	    // std::cout << "Jacobi Matrix: \n"; 
	    // std::cout << _mat_Jacobi << std::endl;
	    // end of debugging-------------------
    #endif
    };

/**
      * solve the equilibrium reaction system 
      * using the Newton iterations, input values are
      * the vector of concentration values
      * the convergence tolerance, 
      * the maximum number of iterations. 
      * and the status of convergence
      */
    void solve_EqSys_Newton(LocalVector & vec_conc, 
                            size_t & result, 
                            size_t & node_idx, 
                            double iter_tol,
                            double rel_tol,
                            double max_iter)
    {
        LocalVector x, x_new;
        LocalVector dx; 
        LocalVector total_mass; 
        LocalVector conc_basis  = LocalVector::Zero( _I_basis ); 
        LocalVector conc_second = LocalVector::Zero( _I_second ); 
        // number of iterations
        size_t i, j, iter;
        double d_norm, d1_norm, d_norm0; 

        x              = LocalVector::Zero( _I ); 
        x_new          = LocalVector::Zero( _I );
        dx             = LocalVector::Zero( _I );
        total_mass     = LocalVector::Zero( _I_basis ); 
        // calculate the bulk composition 
        // in terms of basis species
        conc_basis  = vec_conc.head( _I_basis ); 
        conc_second = vec_conc.tail( _I_second ); 

        this->calc_tot_mass( conc_basis, conc_second, total_mass ); 

        #ifdef _DEBUG
            /*
            std::cout << "Total mass for the basis: \n";
            std::cout << total_mass << "\n"; 
            std::cout << "Conc_basis: \n";
            std::cout << conc_basis << "\n"; 
            std::cout << "Conc_second: \n";
            std::cout << conc_second << std::endl; 
            */
        #endif

        // unknown vector is composed of 
        // aqueous mobile basis species
        for ( i=0; i < ( _I_basis + _I_sec_mob + _I_sec_sorp ); i++ )
            x(i) = std::log( vec_conc(i) ); 
        // and the amount of mineral
		x.segment(_I_basis + _I_sec_mob + _I_sec_sorp, _I_sec_min) = vec_conc.segment(_I_basis + _I_sec_mob + _I_sec_sorp, _I_sec_min);

        // start solving the system
        iter = 0; 
        dx = LocalVector::Ones(this->_I); 
        // now updating the saturation index and minerals
        this->update_minerals( x, total_mass );
        // evaluate the residual
        this->calc_residual(x, total_mass, _vec_res); 
        d_norm0 = _vec_res.norm(); 
        d_norm  = d_norm0; 
        while (true)
        {
            #ifdef _DEBUG
                // display the residual
                // std::cout << "Iteration #" << iter << "||res|| = " << _vec_res.norm() << "||delta_x|| = " << dx.norm() << std::endl; 
            #endif
            // convergence criteria
            if ( d_norm < iter_tol  )
            {
                #ifdef _DEBUG
                    // std::cout << "Newton iteration successfully converged!\n"; 
                #endif
                result = 0;
                // update concentrations
                this->update_concentations( x, vec_conc );
                break;  // break the loop
            }
            else if ( dx.norm() < rel_tol )
            {
                #ifdef _DEBUG
                // std::cout << "Warning, Newton iteration stagnent on Node #" << node_idx << "! Exit the iteration!\n" ; 
                #endif
                result = 0;
                // update concentrations
                this->update_concentations( x, vec_conc );
                break;  // break the loop
            }
            else if ( iter > max_iter )
            {
                #ifdef _DEBUG
                std::cout << "ERROR! Node #" << node_idx  << "Newton iterationan does not converge! Simulation stops!\n"; 
                #endif
                result = 1; 
                return; // stop the program
            }
            // form Jacobian matrix
            this->calc_Jacobi(x, total_mass, _vec_res); 
            // solving for increment
			this->solve_minimization(_mat_Jacobi, _vec_res, dx);

            // increment of unkowns
            this->increment_unknown( x, dx, x_new ); 

            // update the mineral
            this->update_minerals( x_new, total_mass );
            // evaluate residual with x_new
            this->calc_residual(x_new, total_mass, _vec_res); 


            // line search begins
            j = 0; 
            while ( j < max_iter )
            {
                // d1_norm = norm(res,inf);
                d1_norm = _vec_res.norm(); 
            
                if (d1_norm < d_norm)
                    break;
            
                // cut into half
                dx = dx * 0.5;
            
                // increment of unkowns
                this->increment_unknown( x, dx, x_new ); 
                // now updating the saturation index and minerals
                this->update_minerals( x_new, total_mass );
                // evaluate residual with x_new
                this->calc_residual(x_new, total_mass, _vec_res); 
            
                j++;
            }  // end of while        
            d_norm = d1_norm; 

            // d_norm = _vec_res.norm(); 
            x = x_new; 

            // increase the iteration count
            iter++; 
        }
    }; 

	/** HS 20131106: rewrite the minimization solve function
	* This function is doing the job:
	* Minimize norm(dx) on
	* L(d) := { x belong to Real^n , norm(J*dx + b) = min{ norm(J*dx + b, y belong to Real^n }}
	*
	* The numerical algorithm is exactly following page 47 of the dissertation from
	* Joachim Hoffmann (2010) Reactive Transport and Mineral Dissolution /Precipitation in Porous Media:
	* Efficient Solution Algorithms, Benchmark Computations and Existence of Global Solutions.
	* University of Erlangen-Nuernberg.
	*
	* It is an excellant solution of singular J problem!
	*/
	void solve_minimization(ogsChem::LocalMatrix & J, 
                            ogsChem::LocalVector & b,
							ogsChem::LocalVector & dx)
	{
		// step 0: variable definition and initialization
		int n_r, n_rows_M;
		ogsChem::LocalMatrix Q, R, P, B, RB, V, M, tmp;
		ogsChem::LocalVector z, y, y1, y2;
		Eigen::FullPivHouseholderQR<Eigen::MatrixXd> qr_decomp;
		y = ogsChem::LocalVector::Zero(dx.rows());

		// step 1: perform Housholder QR decomposition on J, 
		// so that Q^T * J * P = { R  B }
		//                       { 0  0 }
		qr_decomp.compute(J);
		Q = qr_decomp.matrixQ();
		P = qr_decomp.colsPermutation();
		n_r = qr_decomp.rank();
		n_rows_M = J.cols() - n_r;
		RB = Q.transpose() * J * P;

		if (n_r == J.cols())
		{
			// if n_rank == n_cols, directly solve
			dx = qr_decomp.solve(-b);
		}
		else
		{
			// step 2: split R and B
			R = RB.topLeftCorner(n_r, n_r);
			B = RB.topRightCorner(n_r, RB.cols() - n_r);

			// step 3: if n_rank < n_cols, calculate V, z and y based on R and B. 
			// solve R*V = B
			qr_decomp.compute(R);
			V = qr_decomp.solve(B);
			// Rz = (Q^T *(-b))
			// (I + V^TV)*y2 = V^T * z
			M = ogsChem::LocalMatrix::Identity(n_rows_M, n_rows_M) + V.transpose() * V;
			tmp = (Q.transpose() * (-1.0 * b)).topRows(n_r);
			z = qr_decomp.solve(tmp);
			y2 = M.fullPivHouseholderQr().solve(V.transpose() * z);
			// y1 = z - V*y2
			y1 = z - V * y2;
			// formulate y
			y.head(n_r) = y1;
			y.tail(J.rows() - n_r) = y2;
			// apply permuation
			dx = P * y;
		}
		return;
	};

private:
	/**
      * private flag indicating initialization
      */
	bool isInitialized; 

	/**
      * a list of all kinetic reactions
      */
	std::vector<ogsChem::chemReactionEq*> & _list_eq_reactions; 

    /**
      * stoichiometric matrix S
      */	
	LocalMatrix _matStoi; 

    /**
      * stoichiometric matrix from the input reactions
      */	
	LocalMatrix _matStoi_input;

    /**
      * vector of natural log K values
      */	
    LocalVector _vec_lnK; 
	/**
      * _I is the number of components and _J is the number of reactions
      */
	size_t _I, _J; 

	/**
      * number of mobile, sorption and mineral reactions
      */
	size_t _J_mob, _J_sorp, _J_min;

    /**
      * number of mobile, sorption and mineral components
      */
	size_t _I_mob, _I_sec_mob, _I_sec_sorp, _I_sec_min, _I_kin;

    /**
      * number of basis components
      */
	size_t _I_basis;

    /**
      * number of secondary components
      */
	size_t _I_second;

    /**
      * vector of system residual
      * its size is equal to the number of components
      */
    ogsChem::LocalVector _vec_res;

    /**
      * vector of activity coefficients
      * its size is equal to the number of components
      */
    ogsChem::LocalVector _ln_activity_coeff; 

    /**
      * vector of ln concentrations.
      * notice the c_min are sill in linear scale, not loged. 
      */
    ogsChem::LocalVector _ln_activity; 

    /**
      * vector of staturation index
      * size of this vector is equal to the amount of minerals
      * if AI = 1, meaning the mineral is saturated and present
      * if AI = 0, meaning the mineral is under-saturated and not present
      */
    ogsChem::LocalVector _AI;

    /**
      * matrix of system Jacobi
      * its size is equal to the n_components by n_components
      */
    ogsChem::LocalMatrix _mat_Jacobi; 
    
    /**
      * pointer to the activity model
      */
    ogsChem::chemActivityModelAbstract *_activity_model; 

    /**
      * count the components 
      */
    void countComp(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp)
    {
	    _I_mob = 0; 
        _I_sec_mob = 0; 
	    _I_sec_sorp= 0; 
	    _I_sec_min = 0; 
		_I_kin = 0; 

	    BaseLib::OrderedMap<std::string, ogsChem::ChemComp*>::iterator it; 
	    for( it = map_chemComp.begin(); it != map_chemComp.end(); it++ )
	    {
		    switch ( it->second->getCompType() )
		    {
			case ogsChem::BASIS_COMP:
				_I_mob++; 
				break;
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
				_I_kin++;
				break;
			default:
			    _I_sec_min++;
			    break; 
		    }
	    }

		_I = _I_mob + _I_sec_sorp + _I_sec_min ;
    }; 

    /**
      * count how many mobile, sorption and mineral reactions
      */
    void countReactions(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp, std::vector<ogsChem::chemReactionEq*> & list_eq_reactions)
    {
	    _J_mob = 0; 
	    _J_sorp= 0; 
	    _J_min = 0; 

	    size_t it; 
	    for( it=0; it < list_eq_reactions.size(); it++ )
	    {
            switch ( list_eq_reactions[it]->get_type() )
		    {
            case ogsChem::MOB_EQ_REACT: 
			    _J_mob++; 
			    break;
            case ogsChem::SORP_EQ_REACT: 
			    _J_sorp++;
			    break;
            case ogsChem::MIN_EQ_REACT: 
			    _J_min++;
			    break;
		    default:
			    _J_min++;
			    break; 
		    }  // end of switch
	    }  // end of for

        _J = _J_mob + _J_sorp + _J_min;
    };

    /**
      * construct stoichiometric matrix out of list of components and reactions
      */
    void buildStoi(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp, 
		           std::vector<ogsChem::chemReactionEq*>                & list_eq_reactions)
    {
	    size_t i,j, tmp_idx; 
	    double tmp_stoi; 
	    std::string tmp_str;
	    BaseLib::OrderedMap<std::string, ogsChem::ChemComp*>::iterator tmp_Comp;
	
	    // creat the memory for Stoi matrix
	    _matStoi_input = LocalMatrix::Zero(_I, _J); 
	
	    // based on the reactions, fill in the stoi matrix
	    // loop over all reactions
	    for ( j=0; j < list_eq_reactions.size(); j++ )
	    {	// for each reaction
		    // find each participating components 
		    for ( i=0; i < list_eq_reactions[j]->get_vecCompNames().size(); i++ ){
			    tmp_str  = list_eq_reactions[j]->get_vecCompNames()[i]; 
			    tmp_Comp = map_chemComp.find(tmp_str);
				if ( tmp_Comp == map_chemComp.end() ) // this component not found!
				{
					std::string err_str("The component ");
					err_str.append(tmp_str); 
					err_str.append(" in reaction definition not found in MCP list!"); 
					ERR(err_str.c_str());
					exit(1); 
				}
			    tmp_idx  = tmp_Comp->second->getIndex(); 
			    if ( list_eq_reactions[j]->get_vecStoi().size() > 2 )
			    {
				    // normal reactions
				    tmp_stoi = list_eq_reactions[j]->get_vecStoi()[i];
				    // and put them into Stoi matrix
				    _matStoi_input(tmp_idx,j) = tmp_stoi; 
			    }
			    else  // this is a basis component
			    {
			        _matStoi_input(tmp_idx,j) = 1.0; 
			    }

		    }  // end of for i
	    }  // end of for j

        // number of secondary components
        _I_second = _J ; 
        // number of basis components
        _I_basis = _I - _I_second; 
        // number of secondary mobile components
        _I_sec_mob = _I_mob - _I_basis ; 
        // organize which components are basis
        // and which are secondary
        _matStoi = _matStoi_input.topRows( _I_basis ).transpose(); 

    #ifdef _DEBUG
	    // debugging--------------------------
	    std::cout << "Stoichiometric Matrix S: " << std::endl; 
	    std::cout << _matStoi << std::endl;
	    // end of debugging-------------------
    #endif
    };

    /**
      * reading in the logK values oaf each equilibrium reactions
      */
    void read_logK(std::vector<ogsChem::chemReactionEq*> & list_eq_reactions)
    {
        size_t i; 
        _vec_lnK = ogsChem::LocalVector::Zero(_J); 
        for ( i=0 ; i < list_eq_reactions.size(); i++ )
            _vec_lnK(i) = list_eq_reactions[i]->get_ln_K(); 
    };

    /**
      * increment the unknown with a damping factor
      * to prevent the concentrations falling into negative values
      */    
    void increment_unknown(LocalVector & x_old, 
                           LocalVector & delta_x, 
                           LocalVector & x_new)
    {
        size_t i, n_unknowns;
        double damp_factor; 

        n_unknowns = x_old.rows();
        // increment with a damping factor for the minerals
        for (i=0; i<n_unknowns; i++)
        {
            if ( i < ( _I_basis + _I_sec_mob + _I_sec_sorp ) )
                x_new(i) = x_old(i) + delta_x(i);
            else
            {
                damp_factor = 1.0 / std::max(1.0, -1.33*delta_x(i) / x_old(i) );
                x_new(i) = x_old(i) + damp_factor * delta_x(i);
            }
        }  // end of for

    };  // end of func increment_unknown;

    /**
      * update the concentration of minerals
      */
    void update_minerals(LocalVector & vec_unknowns, LocalVector & mass_constrain)
    {
        size_t i, idx; 
        double cbarmin, phi; 
        ogsChem::LocalMatrix Stoi_min; 
        ogsChem::LocalVector logK_min;
        ogsChem::LocalVector c_sec_min; 
        ogsChem::LocalVector c_basis, c_second; 
        ogsChem::LocalVector ln_activity_basis;
        ogsChem::LocalVector ln_c_basis; 
        ogsChem::LocalVector ln_activity;
        ogsChem::LocalVector ln_activity_coeff;

        Stoi_min = this->_matStoi.bottomRows(_I_sec_min);
        logK_min = this->_vec_lnK.tail(_I_sec_min); 

        c_basis = ogsChem::LocalVector::Zero( _I_basis ); 
		c_second = ogsChem::LocalVector::Zero(_I_second); 
        ln_activity = ogsChem::LocalVector::Zero( _I    ); 
        ln_activity_coeff = ogsChem::LocalVector::Zero( _I ); 
        ln_activity_basis = ogsChem::LocalVector::Zero( _I_basis ); 

        this->_activity_model->calc_activity_logC( vec_unknowns, ln_activity_coeff, ln_activity ); 

        // take the first section which is basis concentration
        ln_activity_basis = ln_activity.head( _I_basis );
        for ( i=0; i < _I_basis; i++)
            c_basis(i) = std::exp( ln_activity(i) - ln_activity_coeff(i) ); 
		for (i = 0; i < (_I_sec_mob + _I_sec_sorp); i++)
			c_second(i) = std::exp(ln_activity(_I_basis + i) - ln_activity_coeff(_I_basis + i)); 
        // notice that we intentially leave the minerals to be zero here. 

        // and the minerals
		c_sec_min = vec_unknowns.segment( _I_basis + _I_sec_mob + _I_sec_sorp, _I_sec_min);

        for ( i=0; i < _I_sec_min; i++ )
        {
            idx = _I_basis + _I_sec_mob + _I_sec_sorp + i; 
            if ( _AI(i) == 1 )
            {
                cbarmin = cal_cbarmin_by_total_mass(i, c_basis, c_second, mass_constrain);
            }  // end of if AI(i)

            phi  = -logK_min(i) + Stoi_min.row(i) * ln_activity_basis;
        
            if ( phi > cbarmin )
            {
                _AI(i) = 0;
                cbarmin = 0.0; 
            }
            else
            {
                if ( _AI(i) == 0 )
                {
                    cbarmin = cal_cbarmin_by_total_mass(i, c_basis, c_second, mass_constrain);
                }
                _AI(i) = 1; 
            }  // end of else

            vec_unknowns(idx) = cbarmin;
        }  // end of for 

    };  // end of function update_minerals

    /**
      * calcuate one particular mineral concentration 
      * by the amount of basis concentration and total mass constrain
      */
    double cal_cbarmin_by_total_mass(size_t idx_min, LocalVector & c_basis, LocalVector & c_second, LocalVector & tot_mass)
    {
        double cbarmin;
        ogsChem::LocalVector conc_min, res_tmp;
        ogsChem::LocalMatrix matStoi_min_trans; 
		conc_min    = ogsChem::LocalVector::Zero(_I_second);
        res_tmp     = ogsChem::LocalVector::Zero( _I_basis );
    
		matStoi_min_trans = _matStoi.bottomRows( _J_min ).transpose();
        res_tmp       = tot_mass - c_basis + _matStoi.transpose() * c_second; 

		conc_min = (-1.0 * matStoi_min_trans).householderQr().solve(res_tmp);
		cbarmin  = conc_min( idx_min );
    return cbarmin;
    };  // end of function cal_cbarmin_by_total_mass

    /**
      * update all concentrations based on p
      */
    void update_concentations(LocalVector & vec_unknowns, LocalVector & vec_concentrations)
    {
        size_t i; 

        for ( i=0; i < ( _I_basis + _I_sec_mob + _I_sec_sorp ) ; i++ )
            vec_concentrations(i) = std::exp( vec_unknowns(i) );
        vec_concentrations.tail( _I_sec_min ) = vec_unknowns.tail( _I_sec_min ); 
    };

};

}  // end of namespace

#endif