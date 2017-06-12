/**
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LocalProblem__EOS->cpp
 *
 * Created on 2014-05-14 by Yonghui HUANG & Haibing Shao
 */


#include "logog.hpp"
#include <algorithm>
#include "math.h"
#include "LocalProblem_EOS_TotalDensity.h"
#include "EOS_TotalDensityForm.h"

/**
  * constructor of the class
  */
LocalProblem_EOS_TotalDensity::LocalProblem_EOS_TotalDensity()
{
	// initialize the corresponding EOS
	// _EOS = new EOS_H2_H2O_ISOT();
	_EOS = new EOS_TotalDensityForm();

	// get the degree of freedom
	N = _EOS->get_n(); 

	Res = ogsChem::LocalVector::Ones(N); 
	Matrix_Jacobian = ogsChem::LocalMatrix::Zero(N, N);
	/*test */
	
}

/**
  * destructor of the class
  */
LocalProblem_EOS_TotalDensity::~LocalProblem_EOS_TotalDensity(void)
{
	BaseLib::releaseObject(_EOS);
}

/**
 * solve the EOS
 * In our case, we have 2 inputs, P and X
 * and 8 outputs, which are the secondary variables.
 * the initial guess will be embedded into the Output vector.
 */
void LocalProblem_EOS_TotalDensity::solve(ogsChem::LocalVector & Input, ogsChem::LocalVector & Output)
{
	ogsChem::LocalVector  vec_rest_var = ogsChem::LocalVector::Zero(4);
	if (Input.size() == 2 && Output.size() == 7)
	{
		
		std::size_t m_flag = 1; 
		//This part is for semi-smooth newton iteration of local problem with complementary condition
		_EOS->set_env_condition(Input);
		U_ini = Output.head(3);
		this->solve_LocalProblem_Newton_LineSearch(m_flag);
		
		if (m_flag == 0)
		{
			Output.head(3) = U_cur;
			_EOS->calc_Rest_SecVar(U_cur, vec_rest_var);
			Output.tail(4) = vec_rest_var;
			//std::cout << U_cur << std::endl;
		}
		else
		{
			WARN("Solving local EOS problem does not converge! \n Using old values as seoncdary varibales. \n"); 
		}
	}
	else
	{
		ERR("When solving local EOS problem, the size of input and out vector is not correct! ");
		exit(1); 
	}

}
//void LocalProblem_EOS::der

/**
  * TODO: describe this function
  */
void LocalProblem_EOS_TotalDensity::solve_minimization(ogsChem::LocalMatrix & J,
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
}

/**
  * Newton iteration with line search
  */
void LocalProblem_EOS_TotalDensity::solve_LocalProblem_Newton_LineSearch(std::size_t & flag)
{
	// flag = 0: converged within max_iter; 
	//      = 1: hit max_iter without converg;
	flag = 1; 

	ogsChem::LocalVector U_pre;// vec_residual;
	ogsChem::LocalVector deltaU;

	// number of iterations
	size_t j(0), Iter(0);//j for line search
	const double neta(0.5);// damping factor for ..
	const double alpha(0.2971);// damping factor for line search
	double d_norm(0.0), d1_norm(0.0);

	/**
	*  Initialize
	*/
	U_cur = ogsChem::LocalVector::Zero(N);
	U_pre = ogsChem::LocalVector::Zero(N);
	deltaU = ogsChem::LocalVector::Ones(N);
	//Matrix_Jacobian = ogsChem::LocalMatrix::Zero(N, N);
	//Res = ogsChem::LocalVector::Zero(N);

	// start solving the system


#ifdef _DEBUG
	// debugging--------------------------
	// std::cout << "Residual Vector: \n";
	// std::cout << vec_residual << std::endl;
	// end of debugging-------------------
#endif

	// save the previous values
	U_cur = U_ini;
	U_pre = U_ini;
	_EOS->eval(U_pre, Res);
	d_norm = Res.norm();
	while (true)
	{

		if (Iter > Iter_tol)
		{
			flag = 1;
			//INFO("residuals : %1.3e", d_norm);
			break; 
		}
		if (d_norm < tot)
		{
			flag = 0; 
			//test for MoMaS benchmark 
			//if (U_cur(0) < 1e-14)
				//U_cur(0) = 0;
			break; 
		}

		// construct Jacobian matrix
		_EOS->calc_Jacobian(U_pre, Matrix_Jacobian);
		
        #ifdef _DEBUG
		// debugging--------------------------
		//std::cout << "Matrix_Jacobian: \n";
		//std::cout << Matrix_Jacobian << std::endl;
		//std::cout << "Residual_Vector: \n";
		//std::cout << Res << std::endl;
		//end of debugging-------------------
        #endif

		//deltaU = Res / Matrix_Jacobian;
		this->solve_minimization(Matrix_Jacobian, Res, deltaU);

        #ifdef _DEBUG
		// debugging--------------------------
		//std::cout << "dx Vector: \n";
		//std::cout << deltaU << std::endl;
		// end of debugging-------------------
        #endif
		
		U_cur = U_pre + neta*deltaU;//  set the damping factor as 0.5

		// evaluate residual with x_new
		_EOS->eval(U_cur, Res);
		// std::cout << "Norm of Residual is "<< std::endl;
		// std::cout << Res.norm() << std::endl;
		j = 0;
		// line search begins
		while (j < Iter_tol)
		{
			// d1_norm = norm(res,inf);
			d1_norm = Res.norm();

			if (d1_norm < d_norm)
				break;

			// updating dx
			deltaU = deltaU * alpha;//

			// increment of unknowns
			//U_cur = U_pre  + deltaU;
			U_cur = U_cur  - deltaU;

			// evaluate residual with x_new
			_EOS->eval(U_cur, Res);

            #ifdef _DEBUG
			// std::cout << "vec_residual: \n";
			// std::cout << Res << std::endl;
			// display the residual
			// std::cout << "Line Search Iteration #" << Iter << "||res|| = " << d1_norm << "||delta_x|| = " << Res.norm() << std::endl;
            #endif

			j++;
		}  // end of while
		d_norm = d1_norm;
		U_pre = U_cur;
		
		// increase the iteration count
		Iter++;
	}
}


