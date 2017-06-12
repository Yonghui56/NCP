/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemEqReactSys.cpp
 *
 * Created on 2013-05-26 by Reza Zolfaghari & Haibing Shao
 */

#include "LocalProblem.h"
#include "logog.hpp"
//#define _DEBUG

LocalProblem::LocalProblem(ogsChem::chemReductionGIA* ReductionGIA, MathLib::StepperBulischStoer<Local_ODE_Xi_immob_GIA>* sbs, Local_ODE_Xi_immob_GIA* local_ode_xi_immob_GIA)
: _n_Comp(ReductionGIA->get_n_Comp()), _ReductionGIA(ReductionGIA), _sbs(sbs), _local_ode_xi_immob_GIA(local_ode_xi_immob_GIA), _I_mob(_ReductionGIA->get_n_Comp_mob()), _I_sorp(_ReductionGIA->get_n_Comp_sorb()), _I_min(_ReductionGIA->get_n_Comp_min()), _I_NMin_bar(_ReductionGIA->get_n_Comp_sorb()), _I_kin(_ReductionGIA->get_n_Comp_kin())
    , _J_mob(_ReductionGIA->get_J_mob()), _J_sorp(_ReductionGIA->get_J_sorp()), _J_min(_ReductionGIA->get_J_min())
    , _n_xi_Kin_bar(_ReductionGIA->get_n_xi_Kin_bar())
	, _n_xi_Mob(_ReductionGIA->get_n_xi_Mob()), _n_eta(_ReductionGIA->get_n_eta()), _n_eta_bar(_ReductionGIA->get_n_eta_bar()), _n_xi_Sorp_tilde(_ReductionGIA->get_n_xi_Sorp_tilde()), _n_xi_Min_tilde(_ReductionGIA->get_n_xi_Min_tilde())
    , _n_xi_Sorp(_ReductionGIA->get_n_xi_Sorp()), _n_xi_Min(_ReductionGIA->get_n_xi_Min()), _n_xi_Sorp_bar_li(_ReductionGIA->get_n_xi_Sorp_bar_li()), _n_xi_Sorp_bar_ld(_ReductionGIA->get_n_xi_Sorp_bar_ld()), _n_xi_Kin(_ReductionGIA->get_n_xi_Kin())
	, _mat_c_mob_2_xi_mob(_ReductionGIA->get_matrix_C2Xi()), _mat_c_immob_2_xi_immob(_ReductionGIA->get_matrix_Cbar2XiBar()), _n_xi_Sorp_bar(_ReductionGIA->get_n_xi_Sorp_bar()), _mat_S1min (_ReductionGIA->get_matrix_S1min())
    , _mat_Ald(_ReductionGIA->get_matrix_Ald()),_mat_Ssorp(_ReductionGIA->get_matrix_Ssorp()), _mat_A2kin(_ReductionGIA->get_matrix_A2kin()), _mat_S1mob(_ReductionGIA->get_matrix_S1mob()),_mat_c_mob_2_eta_mob(_ReductionGIA->get_matrix_C2Eta())
    , _n_xi_Min_bar(_ReductionGIA->get_n_xi_Min_bar()), _mat_c_immob_2_eta_immob(_ReductionGIA->get_matrix_C2EtaBar()), _list_kin_reactions(_ReductionGIA->get_list_kin_reactions()), _J_tot_kin(_ReductionGIA->get_n_xi_Kin_total()), _n_xi_global(_ReductionGIA->get_n_xi_global())
	, _n_xi_local(_ReductionGIA->get_n_xi_local())
	,_logk_mob(_ReductionGIA->get_logk_mob()), _logk_sorp(_ReductionGIA->get_logk_sorp()), _logk_min(_ReductionGIA->get_logk_min())
	,_activity_model(ReductionGIA->get_activity_model()) //RZ: 16.12.2013

{

}

void LocalProblem::solve_LocalProblem_Newton_LineSearch(std::size_t & node_idx,
														double dt,
														const double iter_tol,
														const double rel_tol,
														const double max_iter, 
														ogsChem::LocalVector & x,
														ogsChem::LocalVector & vec_eta,
														ogsChem::LocalVector & vec_etabar,
														ogsChem::LocalVector & vec_xi_local,
														ogsChem::LocalVector & vec_xi_global,
														ogsChem::LocalVector & vec_xi_bar_kin_old,
														ogsChem::LocalVector & lnk_mob,
														ogsChem::LocalVector & lnk_sorp,
														ogsChem::LocalVector & lnk_min, 
														ogsChem::LocalVector & vec_AI)
{
    ogsChem::LocalVector x_new, vec_residual;
    ogsChem::LocalVector dx;

    // TODO initialize water content and deltaT

    // number of iterations
    size_t j, iter;
    const double alpha (0.5);
    double d_norm(0.0), d1_norm(0.0);

	// HS: the number of unknowns in the local problem equals
	// to the number of chemical components
    x_new          = ogsChem::LocalVector::Zero( _n_Comp );
	dx             = ogsChem::LocalVector::Ones( _n_Comp );
	_mat_Jacobian  = ogsChem::LocalMatrix::Zero( _n_Comp, _n_Comp);
	vec_residual   = ogsChem::LocalVector::Zero( _n_Comp );

	// vec_tot_mass_constrain contains xi global and eta mobile which acts as a total mass constrain for the local problem.
	_vec_eta            = vec_eta;
	_vec_etabar         = vec_etabar;
	_vec_XiSorpTilde    = vec_xi_global.head(_n_xi_Sorp_tilde);
	_vec_XiMinTilde     = vec_xi_global.segment(_n_xi_Sorp_tilde, _n_xi_Min_tilde);
	_vec_Xikin          = vec_xi_global.tail(_n_xi_Kin);
	_vec_XiBarKin       = vec_xi_local.tail(_n_xi_Kin_bar);
	_vec_XiBarKin_old   = vec_xi_bar_kin_old;

    // start solving the system
    iter = 0; 

    // now updating the saturation index and minerals
    if(_n_xi_Min > 0)
    {
    this->calculate_AI(x,vec_AI);
    this->update_minerals_conc_AI( x, vec_AI );
    }

#ifdef _DEBUG
	// debugging--------------------------
    // std::cout << "x Vector: \n";
    // std::cout << x << std::endl;
	// std::cout << "vec_AI Vector: \n";
    // std::cout << vec_AI << std::endl;
	// end of debugging-------------------
#endif

	// update the value of xikinbar in the vector of unknowns(x_new)
	// HS: deactivate the following block if using backward Euler for the kinetic evaluation. 
	/*
	if (_n_xi_Kin_bar > 0)
		this->ODE_solver(dt, x, _vec_XiBarKin, _vec_XiBarKin_old);
	*/
    // evaluate the residual
	this->calc_residual(dt, x, _vec_XiBarKin_old, vec_residual, _vec_XiBarKin, vec_AI);

    // evaluate norm of residual vector
    d_norm = vec_residual.norm();

#ifdef _DEBUG
	// debugging--------------------------
	// std::cout << "Residual Vector: \n";
	// std::cout << vec_residual << std::endl;
	// end of debugging-------------------
#endif

	// save the previous values
	x_new  =  x;
	while ((iter < max_iter) && (d_norm > iter_tol) && (dx.norm() > rel_tol )){ //RZ: 4.12.2013 // HS 2014Jan14

        // form Jacobian matrix
		this->calc_Jacobian(dt, x, vec_residual, _vec_XiBarKin, vec_AI);

#ifdef _DEBUG
	// debugging--------------------------
	// std::cout << "_mat_Jacobian: \n";
	// std::cout << _mat_Jacobian << std::endl;
	// end of debugging-------------------
#endif

        // solving for increment
        this->solve_minimization(_mat_Jacobian, vec_residual, dx );

#ifdef _DEBUG
	// debugging--------------------------
  	// std::cout << "dx Vector: \n";
    // std::cout << dx << std::endl;
    // end of debugging-------------------
#endif

    // increment of unknowns
    this->increment_unknown( x, dx, x_new ); 

#ifdef _DEBUG
	// debugging--------------------------
	// std::cout << "x_new Vector: \n";
	// std::cout << x_new << std::endl;
	// end of debugging-------------------
#endif

    	// updating the saturation index and minerals
        if(_n_xi_Min > 0)
		this->update_minerals_conc_AI( x_new, vec_AI );

#ifdef _DEBUG
	// debugging--------------------------
//	std::cout << "x_new: \n";
//	std::cout << x_new << std::endl;
	// end of debugging-------------------
#endif
		// HS: deactivate the following block if using backward Euler for the kinetic evaluation. 
		/*
		// update the value of xikinbar in the vector of unknowns(x_new)
		if (_n_xi_Kin_bar > 0)
			this->ODE_solver(dt, x_new, _vec_XiBarKin, _vec_XiBarKin_old);
	    */
        // evaluate residual with x_new
		this->calc_residual(dt, x_new, _vec_XiBarKin_old, vec_residual, _vec_XiBarKin, vec_AI);


#ifdef _DEBUG
	// debugging--------------------------
	// std::cout << "x_new Vector: \n";
	// std::cout << x_new << std::endl;

	// std::cout << "vec_residual: \n";
	// std::cout << vec_residual << std::endl;
	// end of debugging-------------------
#endif

        // line search begins
        j = 0; 
        while ( j < 30 )
        {
            // d1_norm = norm(res,inf);
            d1_norm = vec_residual.norm();
            
            if (d1_norm < d_norm)
                break;
            
            // updating dx
            dx = dx * alpha;

            // increment of unknowns
            this->increment_unknown( x, dx, x_new );
            //RZ08Jan2014: this increment slightly decrease the overall simulation time.
            //x_new = x_new - dx;

            // now updating the saturation index and minerals
            if(_n_xi_Min > 0)
				this->update_minerals_conc_AI( x_new, vec_AI );

            // update the value of xikinbar in the vector of unknowns(x_new)
			if (_n_xi_Kin_bar > 0)
			this->ODE_solver(dt, x_new, _vec_XiBarKin, _vec_XiBarKin_old);

			// evaluate residual with x_new
			this->calc_residual(dt, x_new, _vec_XiBarKin_old, vec_residual, _vec_XiBarKin, vec_AI);

		#ifdef _DEBUG
        	// std::cout << "vec_residual: \n";
        	// std::cout << vec_residual << std::endl;
            // display the residual
            // std::cout << "Line Search Iteration #" << iter << "||res|| = " << d1_norm << "||delta_x|| = " << dx.norm() << std::endl;
		#endif

            j++;
        }  // end of while
        d_norm = d1_norm; 
        x = x_new; 

        // increase the iteration count
        iter++; 
    }
}

void LocalProblem::ODE_solver(double dt,
							  ogsChem::LocalVector & vec_unknowns,
							  ogsChem::LocalVector & vec_Xi_Kin_bar, 
							  ogsChem::LocalVector & vec_Xi_Kin_bar_old)
{
	MathLib::LocalVector conc_Mob        = MathLib::LocalVector::Zero(_I_mob);
	MathLib::LocalVector ln_conc_Mob     = MathLib::LocalVector::Zero(_I_mob);
	MathLib::LocalVector conc_Sorp       = MathLib::LocalVector::Zero(_I_sorp);
	MathLib::LocalVector ln_conc_Sorp    = MathLib::LocalVector::Zero(_I_sorp);
	MathLib::LocalVector vec_xi_kin_rate = MathLib::LocalVector::Zero(_J_tot_kin);
	MathLib::LocalVector conc_Min_bar    = MathLib::LocalVector::Zero(_I_min);
	MathLib::LocalVector conc_Kin_bar    = MathLib::LocalVector::Zero(_I_kin);

    ln_conc_Mob  = vec_unknowns.segment(0, _I_mob );  // ln scale
    // HS: NonMin refers to sorbed concentrations only, 
	// kinetic concentrations are not included.
	ln_conc_Sorp = vec_unknowns.segment(_I_mob, _I_sorp);  // ln scale
	conc_Min_bar = vec_unknowns.segment(_I_mob + _I_sorp, _I_min);  // linear scale
	conc_Kin_bar = vec_unknowns.tail(_I_kin);  // linear scale

    MathLib::LocalVector  loc_eta        	= MathLib::LocalVector::Zero( _n_eta );
    MathLib::LocalVector  loc_eta_bar     	= MathLib::LocalVector::Zero( _n_eta_bar );
    MathLib::LocalVector  loc_xi_global     = MathLib::LocalVector::Zero( _n_xi_global );
    MathLib::LocalVector  loc_xi_local      = MathLib::LocalVector::Zero( _n_xi_local );
    MathLib::LocalVector  local_conc        = MathLib::LocalVector::Zero( _n_Comp );

    // convert the ln mobile conc to mobile conc
	ln_conc_Mob  = vec_unknowns.head(_I_mob);
	ln_conc_Sorp = vec_unknowns.segment(_I_mob, _I_sorp); 
    this->cal_exp_conc_vec(_I_mob,  ln_conc_Mob,  conc_Mob  );  // converted to linear scale
	this->cal_exp_conc_vec(_I_sorp, ln_conc_Sorp, conc_Sorp );  // converted to linear scale
	conc_Min_bar = vec_unknowns.segment(_I_mob+_I_sorp, _I_min);
	conc_Kin_bar = vec_unknowns.tail(_I_kin); 

    local_conc.head(_I_mob)				 	     = conc_Mob;
	local_conc.segment(_I_mob, _I_sorp)          = conc_Sorp;
	local_conc.segment(_I_mob + _I_sorp, _I_min) = conc_Min_bar; 
	local_conc.tail(_I_kin)                      = conc_Kin_bar;

    _ReductionGIA->Conc2EtaXi(local_conc, loc_eta, loc_eta_bar, loc_xi_global, loc_xi_local);
    //using ode solver for vec_XiBarKin
	this->_local_ode_xi_immob_GIA->update_eta_xi(loc_eta, loc_eta_bar, loc_xi_global, loc_xi_local);
	// TODO. strictly saying, here should supply time value, not dt value. although it does not influence the result in monod2d example. 
	if (vec_xi_kin_rate.rows() > 0)
		vec_xi_kin_rate = (*_local_ode_xi_immob_GIA)(dt, loc_xi_local); 

    //_sbs->set_y(Xi_Kin_bar);
	_sbs->set_y(vec_Xi_Kin_bar_old);
    _sbs->set_dydx(vec_xi_kin_rate);

	// solve the local ODE problem for xi kin bar
    _sbs->step(dt, _local_ode_xi_immob_GIA);
    //vec_Xi_Kin_bar_new = _sbs->get_y();
    vec_Xi_Kin_bar = _sbs->get_y();

}  // end of function ode solver


void LocalProblem::calc_residual(double dt,
								 ogsChem::LocalVector & vec_unknowns,
								 ogsChem::LocalVector & vec_xi_Kin_bar_old, 
								 ogsChem::LocalVector & vec_residual,
								 ogsChem::LocalVector & vec_Xi_Kin_bar,
								 ogsChem::LocalVector & vec_AI)
{
	ogsChem::LocalVector ln_conc_Mob, ln_conc_Sorp, conc_Mob, conc_Sorp, conc_Min_bar, conc_Kin_bar, conc_bar;
	conc_Mob     = ogsChem::LocalVector::Zero(_I_mob);
	conc_Sorp    = ogsChem::LocalVector::Zero(_I_sorp);
	conc_Min_bar = ogsChem::LocalVector::Zero(_I_min);
	conc_Kin_bar = ogsChem::LocalVector::Zero(_I_kin);

	ln_conc_Mob = vec_unknowns.head(_I_mob);
	ln_conc_Sorp = vec_unknowns.segment(_I_mob, _I_sorp);
	// convert the ln mobile conc to mobile conc
	this->cal_exp_conc_vec(_I_mob, ln_conc_Mob, conc_Mob);
	// convert the ln sorp conc to sorp conc
	this->cal_exp_conc_vec(_I_sorp, ln_conc_Sorp, conc_Sorp);
	conc_Min_bar = vec_unknowns.segment(_I_mob + _I_sorp, _I_min);
	conc_Kin_bar = vec_unknowns.tail(_I_kin); 

	// Eq. 3.55
	this->residual_conc_Mob (ln_conc_Mob, vec_residual);
    // Eq. 3.56
    this->residual_Eta (conc_Mob, vec_residual);
	// Eq. 3.57 - 58
    if(_n_xi_Sorp > 0)
		this->residual_xi_Sorp_tilde(conc_Mob, conc_Sorp, conc_Min_bar, conc_Kin_bar, vec_residual);
	// Eq. 3.59
	if(_n_xi_Min > 0)
		this->residual_xi_Min_tilde(conc_Mob, conc_Sorp, conc_Min_bar, conc_Kin_bar, vec_residual);
    // Eq. 3.61
	if (_J_sorp > 0)
		this->residual_conc_Sorp(ln_conc_Mob, ln_conc_Sorp, vec_residual);	
	// Eq. 3.62
	this->residual_conc_Min(ln_conc_Mob, conc_Min_bar, vec_residual, vec_AI);
	// Eq. 3.63
	this->residual_Eta_bar(conc_Sorp, conc_Min_bar, conc_Kin_bar, vec_residual);

    if(_n_xi_Kin != 0)
    {
    	// Eq. 3.60
    	this->residual_xi_Kin			(conc_Mob, vec_residual);

		// Eq. 3.64
    	this->residual_xi_KinBar_Eq		(conc_Sorp, conc_Min_bar, conc_Kin_bar, vec_residual); 
    	
		// Eq. 3.65
		// RZ: solve xi kin bar using ode solver after solving local problem.
		this->residual_xi_KinBar_Kin(dt, conc_Mob, conc_Sorp, conc_Min_bar, conc_Kin_bar, vec_Xi_Kin_bar, vec_residual);
    }

}  // end of function calc_residual

void LocalProblem::calc_Jacobian(double dt,
								 ogsChem::LocalVector & vec_x,
							     ogsChem::LocalVector & vec_residual,
							     ogsChem::LocalVector & vec_Xi_Kin_bar,
							     ogsChem::LocalVector & vec_AI)
{
	const double delta_xi = 1.0e-8;  //calcite example
	//const double delta_xi = 1.0e-6;    //monod2d
    int i;
    ogsChem::LocalVector vec_x_incremented, vec_residual_incremented;
    vec_residual_incremented = vec_residual;

	for (i = 0; i < _mat_Jacobian.cols(); i++)
	{
		vec_x_incremented  = vec_x;

		// numerical protection
		if (std::fabs(vec_x(i)) > 1.0e-16)
		{
			vec_x_incremented(i) += delta_xi * std::fabs(vec_x(i));
			this->calc_residual(dt, vec_x_incremented, _vec_XiBarKin_old, vec_residual_incremented, vec_Xi_Kin_bar, vec_AI);
			_mat_Jacobian.col(i) = (vec_residual_incremented - vec_residual) / (delta_xi * std::fabs(vec_x(i)));

		}
		else
		{
			vec_x_incremented(i) += delta_xi;
			this->calc_residual(dt, vec_x_incremented, _vec_XiBarKin_old, vec_residual_incremented, vec_Xi_Kin_bar, vec_AI);
			_mat_Jacobian.col(i) = (vec_residual_incremented - vec_residual ) / delta_xi;
		}
	}



#ifdef _DEBUG
	// debugging--------------------------
	// std::cout << "Jacobi Matrix: \n";
	// std::cout << _mat_Jacobian << std::endl;
	// end of debugging-------------------
#endif
}

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
void LocalProblem::solve_minimization(ogsChem::LocalMatrix & J,
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
	Q        = qr_decomp.matrixQ();
	P        = qr_decomp.colsPermutation();
	n_r      = qr_decomp.rank(); 
	n_rows_M = J.cols() - n_r; 
	RB  = Q.transpose() * J * P; 

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
		V   = qr_decomp.solve(B);
		// Rz = (Q^T *(-b))
		// (I + V^TV)*y2 = V^T * z
		M = ogsChem::LocalMatrix::Identity(n_rows_M, n_rows_M) + V.transpose() * V;
		tmp = (Q.transpose() * (-1.0 * b)).topRows(n_r);
		z   = qr_decomp.solve(tmp);
		y2  = M.fullPivHouseholderQr().solve(V.transpose() * z);
		// y1 = z - V*y2
		y1  = z - V * y2;
		// formulate y
		y.head(n_r) = y1; 
		y.tail(J.rows() - n_r) = y2; 
		// apply permuation
		dx = P * y; 
	}
return; 
}

void LocalProblem::increment_unknown(ogsChem::LocalVector & x_old,
                                       ogsChem::LocalVector & delta_x,
                                       ogsChem::LocalVector & x_new)
{
    size_t i, n_unknowns;
	double tmp_value, damp_factor; 
	n_unknowns = x_new.size(); 
	const double alpha(1.33);
	const double eps(std::numeric_limits<double>::epsilon());

	for (i = 0; i < n_unknowns; i++)
	{
		if (i < (_I_mob + _I_sorp))
			x_new = x_old + delta_x;
		else
		{
			// HS 2013Dec25: adding a safty control here
			if (std::fabs(x_old(i)) < eps)
				tmp_value = -1.0 * alpha * delta_x(i) / eps; 
			else
				tmp_value = -1.0 * alpha * delta_x(i) / x_old(i);

			damp_factor = 1.0 / std::max(1.0, tmp_value );
			x_new(i) = x_old(i) + damp_factor * delta_x(i);
		}
	}

}  // end of func increment_unknown


void LocalProblem::update_minerals_conc_AI(ogsChem::LocalVector & vec_unknowns,
									       ogsChem::LocalVector & vec_AI)
{
    size_t i, idx; 
    double  phi(0.0);
    ogsChem::LocalVector ln_conc_Mob  = ogsChem::LocalVector::Zero(_I_mob);
    ogsChem::LocalVector ln_conc_Sorp = ogsChem::LocalVector::Zero(_I_sorp);
    ogsChem::LocalVector conc_Min_bar = ogsChem::LocalVector::Zero(_I_min);
	ogsChem::LocalVector conc_Kin_bar = ogsChem::LocalVector::Zero(_I_kin); 
    //ogsChem::LocalMatrix mat_S1min_transposed;
    //logConc_Mob = ogsChem::LocalVector::Zero(_I_mob);
    ogsChem::LocalMatrix  mat_S1min_transposed = ogsChem::LocalMatrix::Zero(_mat_S1min.cols(),_mat_S1min.rows());
    mat_S1min_transposed = _mat_S1min.transpose();


    // take the first section which is basis concentration
    ln_conc_Mob    = vec_unknowns.head( _I_mob   );
    // take the sorp component concentrations
	ln_conc_Sorp   = vec_unknowns.segment(_I_mob, _I_sorp);
    // take the mineral part
    conc_Min_bar   = vec_unknowns.segment(_I_mob + _I_sorp, _I_min );
	// take the kinetic part
	conc_Kin_bar   = vec_unknowns.tail(_I_kin); 

    //RZ: 16.12.2013 disable incorporating activity coefficients into reaction constant k and using activities instead of concentrations directly in LMA.
    ogsChem::LocalVector ln_activity       = ogsChem::LocalVector::Zero( _n_Comp    );
    ogsChem::LocalVector ln_activity_coeff = ogsChem::LocalVector::Zero( _n_Comp );
    ogsChem::LocalVector ln_Conc 		   = ogsChem::LocalVector::Zero( _n_Comp );
    ln_Conc.head(_I_mob) 		   		   = ln_conc_Mob;  // the log concentrations of immobile nonmineral and linear conc of minerals are also included but not used.
    this->_activity_model->calc_activity_logC( ln_Conc, ln_activity_coeff, ln_activity );
    //End 16.12.2013

	// _AI = 1 if mineral is present and _AI = 0 if mineral is not present.
    for ( i=0; i < _I_min; i++ )
    {
		idx = _I_mob + _I_sorp + i;

    	//RZ: 16.12.2013
    	phi  = -_logk_min(i) + mat_S1min_transposed.row(i) * ln_activity.head(_I_mob);  //RZ:2.12.2013 include activity model

    	// if mineral is present, calculate mineral concentration
    	if	(vec_AI(i) == 1)
    	{
    		// update mineral concentration
			conc_Min_bar(i) = cal_cbarmin_by_constrain(i, ln_conc_Mob, ln_conc_Sorp, conc_Min_bar, conc_Kin_bar);
    		vec_unknowns(idx) = conc_Min_bar(i);
    	}  // end of if

    	if (phi > conc_Min_bar(i))
    	{
    		vec_AI(i) = 0;
    		conc_Min_bar(i) = 0.0;
    		vec_unknowns(idx) = conc_Min_bar(i);
    	} // end of if
    	else
    	{
    		if(vec_AI(i) == 0)
    		{
				conc_Min_bar(i) = cal_cbarmin_by_constrain(i, ln_conc_Mob, ln_conc_Sorp, conc_Min_bar, conc_Kin_bar);
    			vec_unknowns(idx) = conc_Min_bar(i);
    		} // end of if

    		vec_AI(i) = 1;

    	} // end of else

    }  // end of for loop

}  // end of function update_minerals


double LocalProblem::cal_cbarmin_by_constrain(size_t        idx_min,
	                                          ogsChem::LocalVector & ln_conc_Mob,
                                              ogsChem::LocalVector & ln_conc_Sorp,
                                              ogsChem::LocalVector & conc_Min_bar, 
											  ogsChem::LocalVector & conc_Kin_bar)
{

    double cbarmin (0.0);
	ogsChem::LocalVector temp               = ogsChem::LocalVector::Zero( _n_xi_Sorp_bar + _n_xi_Min_bar + _n_xi_Kin_bar);
	ogsChem::LocalVector conc_bar           = ogsChem::LocalVector::Zero( _n_xi_Sorp_bar + _n_xi_Min_bar + _n_xi_Kin_bar);
    ogsChem::LocalVector Conc_Mob 			= ogsChem::LocalVector::Zero(_I_mob);
    ogsChem::LocalVector conc_Sorp        	= ogsChem::LocalVector::Zero(_I_sorp);
    ogsChem::LocalVector xi_mobile          = ogsChem::LocalVector::Zero(_n_xi_Mob + _n_xi_Sorp + _n_xi_Min + _n_xi_Kin);
    ogsChem::LocalMatrix xi_sorp_bar_Ald    = ogsChem::LocalVector::Zero(_n_xi_Sorp_bar_ld);

    // convert the ln mobile conc to mobile conc
    this->cal_exp_conc_vec(_I_mob, ln_conc_Mob, Conc_Mob);
    xi_mobile = _mat_c_mob_2_xi_mob * Conc_Mob;

    // convert the ln nonmineral immobile conc to nonmineral immobile conc
	this->cal_exp_conc_vec(_I_sorp, ln_conc_Sorp, conc_Sorp);

	conc_bar.head(_I_sorp) = conc_Sorp;
	conc_bar.segment(_I_sorp, _I_min) = conc_Min_bar; 
	conc_bar.tail(conc_Kin_bar.rows()) = conc_Kin_bar;  //TODO: RZ: potential bug

    // xi_sorp_bar
    temp 			= _mat_c_immob_2_xi_immob * conc_bar;
    xi_sorp_bar_Ald = _mat_Ald * temp.segment(_n_xi_Sorp_bar_li,_n_xi_Sorp_bar_ld);
    
    cbarmin 		= - _vec_XiMinTilde(idx_min) + xi_mobile(_n_xi_Mob + _n_xi_Sorp_bar_li + idx_min) - xi_sorp_bar_Ald(idx_min);


return cbarmin;
}  // end of function cal_cbarmin_by_total_mass


// take log of concentration vectors
void LocalProblem::cal_ln_conc_vec(size_t                 idx_size,
								   ogsChem::LocalVector & conc_Mob,
								   ogsChem::LocalVector & ln_conc_Mob)
{
	double tmp_x;
	std::size_t i;


	for (i = 0; i < idx_size; i++)
	{
		tmp_x    = conc_Mob(i);
        if (std::abs(tmp_x) < std::numeric_limits<double>::epsilon() || tmp_x <= 0.0 )
            tmp_x = std::numeric_limits<double>::epsilon();
		ln_conc_Mob(i)  = std::log(tmp_x);
	}

}

// take exponential of natural log concentration vectors
void LocalProblem::cal_exp_conc_vec(size_t    		      idx_size,
									ogsChem::LocalVector & ln_conc,
									ogsChem::LocalVector & Conc)
{
	double tmp_x;
	std::size_t i;

	for (i = 0; i < idx_size; i++)
	{
		tmp_x    	 = ln_conc(i);
		Conc(i)  = std::exp(tmp_x);
	}
}

// Eq. 3.55
void LocalProblem::residual_conc_Mob(ogsChem::LocalVector & ln_conc_Mob,
									 ogsChem::LocalVector & vec_residual)
{

    //RZ: 16.12.2013 disable incorporating activity coefficients into reaction constant k and using activities instead of concentrations directly in LMA.
    ogsChem::LocalVector ln_activity       = ogsChem::LocalVector::Zero( _n_Comp    );
    ogsChem::LocalVector ln_activity_coeff = ogsChem::LocalVector::Zero( _n_Comp );
    ogsChem::LocalVector ln_Conc 		   = ogsChem::LocalVector::Zero( _n_Comp );
    ln_Conc.head(_I_mob) 		   		   = ln_conc_Mob;  // the log concentrations of immobile nonmineral and linear conc of minerals are also included but not used.
    this->_activity_model->calc_activity_logC( ln_Conc, ln_activity_coeff, ln_activity );
    //End 16.12.2013

	//RZ: 16.12.2013
	vec_residual.head(_n_xi_Mob) 	= - _logk_mob + _mat_S1mob.transpose() * ln_activity.head(_I_mob);

	//vec_residual.head(_n_xi_Mob) 	= - _logk_mob + _mat_S1mob.transpose() * ln_conc_Mob;

}

// Eq. 3.56
void LocalProblem::residual_Eta(ogsChem::LocalVector & conc_Mob,
		                        ogsChem::LocalVector & vec_residual)
{
	// eta acts as a constrain
	vec_residual.segment(_n_xi_Mob, _n_eta)  = -1.0 * _vec_eta + _mat_c_mob_2_eta_mob * conc_Mob;

#ifdef _DEBUG
	// debugging--------------------------
	//std::cout << "_vec_eta: \n";
	//std::cout << _vec_eta << std::endl;

	//std::cout << "conc_Mob: \n";
	//std::cout << conc_Mob << std::endl;

	//std::cout << "_mat_c_mob_2_eta_mob: \n";
	//std::cout << _mat_c_mob_2_eta_mob << std::endl;

	//std::cout << "vec_residual: \n";
	//std::cout << vec_residual << std::endl;
	// end of debugging-------------------
#endif

}

// Eq. 3.57 - 58
void LocalProblem::residual_xi_Sorp_tilde(ogsChem::LocalVector & conc_Mob,
										  ogsChem::LocalVector & conc_Sorp,
										  ogsChem::LocalVector & conc_Min_bar,
										  ogsChem::LocalVector & conc_Kin_bar,
										  ogsChem::LocalVector & vec_residual)
{
	ogsChem::LocalVector vec_XiSorp      = ogsChem::LocalVector::Zero(_n_xi_Mob + _n_xi_Sorp + _n_xi_Min + _n_xi_Kin);
	ogsChem::LocalVector vec_XiSorpBarLI = ogsChem::LocalVector::Zero(_n_xi_Sorp_bar + _n_xi_Min_bar + _n_xi_Kin_bar);
	ogsChem::LocalVector conc_bar        = ogsChem::LocalVector::Zero(_I_sorp + _I_min + _I_kin);

	conc_bar.head(_I_sorp) = conc_Sorp;
	conc_bar.segment(_I_sorp, _I_min) = conc_Min_bar; 
	conc_bar.tail(_I_kin) = conc_Kin_bar;

	vec_XiSorp        = _mat_c_mob_2_xi_mob * conc_Mob;
	vec_XiSorpBarLI   = _mat_c_immob_2_xi_immob * conc_bar;
    // HS 09.02.2014: The index of vec_residual was wrong in the next line. 
	// vec_residual.segment(_n_xi_Mob + _n_eta, _n_xi_Sorp_tilde) = - _vec_XiSorpTilde + vec_XiSorp.segment(_n_xi_Mob,_n_xi_Sorp_bar_li ) - vec_XiSorpBarLI.segment(0,_n_xi_Sorp_bar_li);
    vec_residual.segment(_J_mob + _n_eta, _n_xi_Sorp_tilde) = -_vec_XiSorpTilde + vec_XiSorp.segment(_n_xi_Mob, _n_xi_Sorp_bar_li) - vec_XiSorpBarLI.segment(0, _n_xi_Sorp_bar_li);
}

// Eq. 3.59
void LocalProblem::residual_xi_Min_tilde(ogsChem::LocalVector & conc_Mob,
										 ogsChem::LocalVector & conc_Sorp,
										 ogsChem::LocalVector & conc_Min_bar,
										 ogsChem::LocalVector & conc_Kin_bar, 
										 ogsChem::LocalVector & vec_residual)
{
	ogsChem::LocalVector  vec_XiMin     = ogsChem::LocalVector::Zero(_n_xi_Mob + _n_xi_Sorp + _n_xi_Min + _n_xi_Kin);
	ogsChem::LocalVector  vec_XiMinBar  = ogsChem::LocalVector::Zero(_n_xi_Sorp_bar + _n_xi_Min_bar + _n_xi_Kin_bar);
	ogsChem::LocalVector conc_bar       = ogsChem::LocalVector::Zero(_I_sorp + _I_min + _I_kin);
	ogsChem::LocalVector A 				= ogsChem::LocalVector::Zero(_n_xi_Min);
	ogsChem::LocalVector B 				= ogsChem::LocalVector::Zero(_n_xi_Sorp_bar_ld);

	conc_bar.head(_I_sorp) = conc_Sorp;
	conc_bar.segment(_I_sorp, _I_min) = conc_Min_bar;
	conc_bar.tail(_I_kin) = conc_Kin_bar;

	vec_XiMin        = _mat_c_mob_2_xi_mob * conc_Mob;
	vec_XiMinBar     = _mat_c_immob_2_xi_immob * conc_bar;

	// HS: 20131103, TODO: The following lines are very strange. Need to double check! 
	A = vec_XiMin.segment(_n_xi_Mob + _n_xi_Sorp_bar_li, _n_xi_Min);
	B = vec_XiMinBar.tail(_n_xi_Sorp_bar_ld);

	if (B.rows() > 0)
        vec_residual.segment(_J_mob + _n_eta + _n_xi_Sorp_tilde, _n_xi_Min_tilde) = -_vec_XiMinTilde + A - conc_Min_bar - (_mat_Ald * B);
	else
        vec_residual.segment(_J_mob + _n_eta + _n_xi_Sorp_tilde, _n_xi_Min_tilde) = -_vec_XiMinTilde + A - conc_Min_bar;
}

// Eq. 3.60
void LocalProblem::residual_xi_Kin(ogsChem::LocalVector & conc_Mob,
								   ogsChem::LocalVector & vec_residual)
{
	ogsChem::LocalVector conc_tmp  = ogsChem::LocalVector::Zero(_n_xi_Mob + _n_xi_Sorp + _n_xi_Min + _n_xi_Kin);
	conc_tmp        = _mat_c_mob_2_xi_mob * conc_Mob;
    vec_residual.segment(_J_mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Kin) = -_vec_Xikin + conc_tmp.segment(_n_xi_Mob + _n_xi_Sorp + _n_xi_Min, _n_xi_Kin);
}

// Eq. 3.61
void LocalProblem::residual_conc_Sorp(ogsChem::LocalVector & ln_conc_Mob,
									  ogsChem::LocalVector & ln_conc_Sorp,
									  ogsChem::LocalVector & vec_residual)
{

	ogsChem::LocalVector   ln_conc_tmp = ogsChem::LocalVector::Zero(_I_mob + _I_sorp);

	ln_conc_tmp.head(_I_mob)  = ln_conc_Mob ;
	ln_conc_tmp.tail(_I_sorp) = ln_conc_Sorp;

    //RZ: 16.12.2013 disable incorporating activity coefficients into reaction constant k and using activities instead of concentrations directly in LMA.
    ogsChem::LocalVector ln_activity       = ogsChem::LocalVector::Zero( _n_Comp    );
    ogsChem::LocalVector ln_activity_coeff = ogsChem::LocalVector::Zero( _n_Comp );
    ogsChem::LocalVector ln_Conc 		   = ogsChem::LocalVector::Zero( _n_Comp );
    ln_Conc.head(_I_mob + _I_sorp) 		   = ln_conc_tmp;  // the log concentrations of immobile nonmineral and linear conc of minerals are also included but not used.
    this->_activity_model->calc_activity_logC( ln_Conc, ln_activity_coeff, ln_activity );
    //End 16.12.2013
    vec_residual.segment(_J_mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde +_n_xi_Kin, _n_xi_Sorp)  = - _logk_sorp + _mat_Ssorp.transpose() * ln_activity.head(_I_mob + _I_sorp);
}

// Eq. 3.62
void LocalProblem::residual_conc_Min(ogsChem::LocalVector & ln_conc_Mob,
	     	 	 	 	 	 	 	 ogsChem::LocalVector & conc_Min_bar,
	     	 	 	 	 	 	 	 ogsChem::LocalVector & vec_residual,
	     	 	 	 	 	 	 	 ogsChem::LocalVector & vec_AI)
{
	size_t        i, idx;
	double        phi(0.0);
	ogsChem::LocalMatrix   mat_S1minT  = ogsChem::LocalMatrix::Zero(_mat_S1min.cols(), _mat_S1min.rows());

	idx  = _n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin + _n_xi_Sorp;

	mat_S1minT = _mat_S1min.transpose();

    //RZ: 16.12.2013 disable incorporating activity coefficients into reaction constant k and using activities instead of concentrations directly in LMA.
    ogsChem::LocalVector ln_activity       = ogsChem::LocalVector::Zero( _n_Comp    );
    ogsChem::LocalVector ln_activity_coeff = ogsChem::LocalVector::Zero( _n_Comp );
    ogsChem::LocalVector ln_Conc 		   = ogsChem::LocalVector::Zero( _n_Comp );
    ln_Conc.head(_I_mob) 		   		   = ln_conc_Mob;  // the log concentrations of immobile nonmineral and linear conc of minerals are also included but not used.
    this->_activity_model->calc_activity_logC( ln_Conc, ln_activity_coeff, ln_activity );
    //End 16.12.2013

	for (i=0; i < _n_xi_Min; i++)
	{
		phi  = -_logk_min(i) + mat_S1minT.row(i) * ln_activity.head(_I_mob); //RZ: 16.12.2013

		if(vec_AI(i) == 1) //RZ: 3-Nov-13 Must be kept like this!
			vec_residual(idx + i) 	= phi;
		else
			vec_residual(idx + i) 	= conc_Min_bar(i);
	//	vec_residual(idx + i) 	= std::min(phi, conc_Min_bar(i));
	}
}

// Eq. 3.63
void LocalProblem::residual_Eta_bar(ogsChem::LocalVector & conc_Sorp,
									ogsChem::LocalVector & conc_Min_bar,
									ogsChem::LocalVector & conc_Kin_bar, 
									ogsChem::LocalVector & vec_residual)
{
	ogsChem::LocalVector conc_bar  = ogsChem::LocalVector::Zero(_I_sorp + _I_min + _I_kin);

	conc_bar.head(_I_sorp) = conc_Sorp;
	conc_bar.segment(_I_sorp, _I_min) = conc_Min_bar; 
	conc_bar.tail(_I_kin) = conc_Kin_bar;

	// etabar acts as a constrain
    vec_residual.segment(_J_mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin + _n_xi_Sorp + _n_xi_Min, _n_eta_bar) = -_vec_etabar + (_mat_c_immob_2_eta_immob * conc_bar);
}

// Eq. 3.64
void LocalProblem::residual_xi_KinBar_Eq(ogsChem::LocalVector & conc_Sorp,
										 ogsChem::LocalVector & conc_Min_bar,
										 ogsChem::LocalVector & conc_Kin_bar,
										 ogsChem::LocalVector & vec_residual)
{

	ogsChem::LocalVector   conc_tmp  = ogsChem::LocalVector::Zero(_n_xi_Sorp_bar + _n_xi_Min_bar + _n_xi_Kin_bar);
	ogsChem::LocalVector   conc_bar = ogsChem::LocalVector::Zero(_I_sorp + _I_min + _I_kin);

	conc_bar.head(_I_sorp)   = conc_Sorp;
	conc_bar.segment(_I_sorp, _I_min) = conc_Min_bar;
	conc_bar.tail(_I_kin)    = conc_Kin_bar;

	conc_tmp     = _mat_c_immob_2_xi_immob * conc_bar;

    vec_residual.segment(_J_mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin + _n_xi_Sorp + _n_xi_Min + _n_eta_bar, _n_xi_Kin_bar) =
			- _vec_XiBarKin + conc_tmp.tail(_n_xi_Kin_bar);
}

// HS disabled. 

// Eq. 3.65
void LocalProblem::residual_xi_KinBar_Kin(double deltaT,
										  ogsChem::LocalVector & conc_Mob,
										  ogsChem::LocalVector & conc_NonMin_bar,
										  ogsChem::LocalVector & conc_Min_bar,
										  ogsChem::LocalVector & conc_kin_bar, 
										  ogsChem::LocalVector & Xi_Kin_bar,
										  ogsChem::LocalVector & vec_residual)
{
	ogsChem::LocalVector  conc 		   = ogsChem::LocalVector::Zero(_n_Comp);
	ogsChem::LocalVector  vec_rateKin  = ogsChem::LocalVector::Zero(_J_tot_kin);
	// ogsChem::LocalVector  temp_vec1    = ogsChem::LocalVector::Zero(_n_xi_Kin_bar);
	// ogsChem::LocalVector  temp_vec2    = ogsChem::LocalVector::Zero(_n_xi_Kin_bar);
	// ogsChem::LocalVector  temp_vec3    = ogsChem::LocalVector::Zero(_n_xi_Kin_bar);

	conc.head	(_I_mob)  					      = conc_Mob;
	conc.segment(_I_mob, _I_NMin_bar) 		      = conc_NonMin_bar;
	conc.segment(_I_mob + _I_NMin_bar, _n_xi_Min) = conc_Min_bar;
	conc.tail(_I_kin) = conc_kin_bar; 

	this->reaction_rates(conc, vec_rateKin);

	// temp_vec1 = (( Xi_Kin_bar - _vec_XiBarKin_old) / deltaT);
	// temp_vec2 = (_mat_A2kin * vec_rateKin);
	// temp_vec3 =  temp_vec1 - temp_vec2;
	// HS disabled. 2014.03.06
	// vec_residual.segment(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde +  _n_xi_Kin + _n_xi_Sorp + _n_xi_Min + _n_eta_bar + _n_xi_Kin_bar, _n_xi_Kin_bar) = temp_vec3;
	Xi_Kin_bar = deltaT * _mat_A2kin * vec_rateKin + _vec_XiBarKin_old; 

}

//problem specific reaction rates
void LocalProblem::reaction_rates(ogsChem::LocalVector & conc,
								  ogsChem::LocalVector & vec_rateKin)
{
	std::size_t i;
	// then calculate the rates and fill them in the rate vector
	for ( i=0; i < this->_J_tot_kin; i++ )
	{
		// get to the particular kin equation and calculate its rate
		this->_list_kin_reactions[i]->calcReactionRate( conc );
		vec_rateKin(i) = this->_list_kin_reactions[i]->getRate();
	}
}

/*
 * RZ: AI vector contains 1 and 0. 1 if mineral is present and 0 if absent.
 */
//calculate AI vector
void LocalProblem::calculate_AI(ogsChem::LocalVector & vec_unknowns,
								ogsChem::LocalVector & vec_AI)
{
    double phi = 0.0;
    size_t idx, j;
    ogsChem::LocalVector ln_conc_Mob  = ogsChem::LocalVector::Zero(_I_mob);
    ogsChem::LocalVector ln_conc_Sorp = ogsChem::LocalVector::Zero(_I_sorp);
    ogsChem::LocalVector conc_Min_bar = ogsChem::LocalVector::Zero(_I_min);
	ogsChem::LocalVector conc_Kin_bar = ogsChem::LocalVector::Zero(_I_kin);
    //ogsChem::LocalMatrix mat_S1min_transposed;
    //logConc_Mob = ogsChem::LocalVector::Zero(_I_mob);
    ogsChem::LocalMatrix  mat_S1min_transposed = ogsChem::LocalMatrix::Zero(_mat_S1min.cols(),_mat_S1min.rows());
    mat_S1min_transposed = _mat_S1min.transpose();


    // take the first section which is basis concentration
    ln_conc_Mob    = vec_unknowns.head( _I_mob   );
    // take the sorp component concentrations
	ln_conc_Sorp   = vec_unknowns.segment(_I_mob, _I_sorp);
    // take the mineral part
    conc_Min_bar   = vec_unknowns.segment(_I_mob + _I_sorp, _I_min );
	// take the kinetic part
    conc_Kin_bar   = vec_unknowns.tail(_I_kin);

    //RZ: 16.12.2013 disable incorporating activity coefficients into reaction constant k and using activities instead of concentrations directly in LMA.
    ogsChem::LocalVector ln_activity       = ogsChem::LocalVector::Zero( _n_Comp    );
    ogsChem::LocalVector ln_activity_coeff = ogsChem::LocalVector::Zero( _n_Comp );
    ogsChem::LocalVector ln_Conc 		   = ogsChem::LocalVector::Zero( _n_Comp );
    ln_Conc.head(_I_mob) 		   		   = ln_conc_Mob;  // the log concentrations of immobile nonmineral and linear conc of minerals are also included but not used.
    this->_activity_model->calc_activity_logC( ln_Conc, ln_activity_coeff, ln_activity );
    //End 16.12.2013

    for(j =0; j < _I_min; j++)
    {
    	idx = _I_mob + _I_sorp + j;
    	conc_Min_bar(j) = cal_cbarmin_by_constrain(j, ln_conc_Mob, ln_conc_Sorp, conc_Min_bar, conc_Kin_bar);

    	//RZ: 16.12.2013
    	phi  = -_logk_min(j) + mat_S1min_transposed.row(j) * ln_activity.head(_I_mob);  //RZ:2.12.2013 include activity model

    	// if mineral concentration  >= phi ; mineral is present; saturated case; precipitate the mineral.
    	if(conc_Min_bar(j) > phi)  //RZ:9-Dec-13
    	{
    		vec_AI(j) = 1;
    	}// end of if
    	// if mineral concentration < phi : mineral is NOT present; under saturated case; dissolve the mineral
    	else
    	{
    		vec_AI(j) = 0;
    		conc_Min_bar(j) = 0.0;
    		vec_unknowns(idx) = conc_Min_bar(j);
    	}// end of else
    }
}

