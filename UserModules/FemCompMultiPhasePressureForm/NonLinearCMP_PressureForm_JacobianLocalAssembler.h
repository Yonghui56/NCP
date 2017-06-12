/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NonLinearCMP_2P2C_JacobianLocalAssembler.h
 *
 * 2015-03-21 by Yonghui Huang
 */

/**
  * This file is same as the MassTransportTimeODELocalAssembler.h
  * The difference is, the compound molecular diffusion coefficient is disabled, 
  */

#ifndef NON_LINEAR_CMP_PRESSUREFORM_JACOBIAN_LOCAL_ASSEMBLER_H
#define NON_LINEAR_CMP_PRESSUREFORM_JACOBIAN_LOCAL_ASSEMBLER_H

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Core/Integration/Integration.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Compound.h"
#include "Ogs6FemData.h"
#include "NonLinearCMP_PressureForm_TimeODELocalAssembler.h"

template <class T_NODAL_FUNCTION_SCALAR, class T_FUNCTION_DATA >//,
class NonLinearCMP_PressureForm_JacobianLocalAssembler : public NumLib::IElementWiseTransientJacobianLocalAssembler
{
public:
	typedef MathLib::LocalVector LocalVectorType;
	typedef MathLib::LocalMatrix LocalMatrixType;
	NonLinearCMP_PressureForm_JacobianLocalAssembler(FemLib::LagrangeFeObjectContainer* feObjects, T_FUNCTION_DATA* func_data, const MeshLib::CoordinateSystem &problem_coordinates)
		: _feObjects(*feObjects), _function_data(func_data), _problem_coordinates(problem_coordinates)
    {
		
	};

	virtual ~NonLinearCMP_PressureForm_JacobianLocalAssembler()
    {
		_function_data = NULL; 
		
    };

	
	T_FUNCTION_DATA* get_function_data(void)
	{
		return _function_data;
	}
	
	void setTheta(double v)
	{
		assert(v >= .0 && v <= 1.0);
		_Theta = v;
	}
	/*
	// inexact jacobian assembly
	void assembly(const NumLib::TimeStep &time, const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &, const MathLib::LocalVector & u1, const MathLib::LocalVector & u0, MathLib::LocalMatrix & localJ)
    {
		
		size_t el(0);
		const size_t n_nodes = e.getNumberOfNodes();
		const size_t mat_id = e.getGroupID();
		const std::size_t elem_id = e.getID();
		//clear the local Jacobian Matrix

		MathLib::LocalMatrix TMP_M(2*n_nodes,2* n_nodes);
		FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
		MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];
		double dt = time.getTimeStepSize();
		LocalMatrixType _M = LocalMatrixType::Zero(2 * n_nodes, 2 * n_nodes);
		LocalMatrixType _K = LocalMatrixType::Zero(2 * n_nodes, 2 * n_nodes);
		_M = _function_data->get_elem_M_matrix().at(elem_id);
		_K = _function_data->get_elem_K_matrix().at(elem_id);
		double eps = 1e-7;// can also be defined as numerical way
		TMP_M = (1.0 / dt)*_M + _Theta*_K;
		//LocalMatrixType Local_LHS=(1/dt)*M+
		//LocalVectorType e_vec = LocalVectorType::Zero(2 * n_nodes);
		//_local_assembler=assembleODE(time, e, u1, u0, M, K);
		for (size_t u_idx = 0; u_idx < 2 * n_nodes; u_idx++)
		{
			//clear the epsilon vectoe 
			LocalVectorType e_vec = LocalVectorType::Zero(2 * n_nodes);
			e_vec(u_idx) = eps*std::abs(u1(u_idx));

			localJ.block(0,u_idx,u1.size(),1) = (TMP_M*(u1 - e_vec) - TMP_M*(u1 + e_vec)) / 2 / eps/std::abs(u1(u_idx));
			//debugging--------------------------
			//std::cout << "localJ: \n";
			//std::cout << localJ << std::endl;
			//end of debugging-------------------
			//localJ=
		}	
    }  // end of function assembly
	*/
	/*
	* exact jacobian matrix assembly
	*/
	void assembly(const NumLib::TimeStep &time, const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &, const MathLib::LocalVector & u1, const MathLib::LocalVector & u0, MathLib::LocalMatrix & localJ)
	{
		size_t el(0);
		const size_t n_nodes = e.getNumberOfNodes();
		const size_t mat_id = e.getGroupID();
		const std::size_t elem_id = e.getID();
		const size_t n_dof = u1.size();
		double dt = time.getTimeStepSize();
		TMP_M_1 = MathLib::LocalMatrix::Zero(n_dof, n_dof);
		TMP_M_2 = MathLib::LocalMatrix::Zero(n_dof, n_dof);
		FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
		MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];

		double eps = 1e-7;

		for (size_t u_idx = 0; u_idx < n_dof; u_idx++)
		{
			M1 = MathLib::LocalMatrix::Zero(n_dof, n_dof);
			K1 = MathLib::LocalMatrix::Zero(n_dof, n_dof);
			F1 = MathLib::LocalVector::Zero(n_dof);
			M2 = MathLib::LocalMatrix::Zero(n_dof, n_dof);
			K2 = MathLib::LocalMatrix::Zero(n_dof, n_dof);
			F2 = MathLib::LocalVector::Zero(n_dof);
			//clear local_r
			local_r_1 = MathLib::LocalVector::Zero(n_dof);
			local_r_2 = MathLib::LocalVector::Zero(n_dof);
			//clear the epsilon vectoe 
			MathLib::LocalVector e_vec = MathLib::LocalVector::Zero(n_dof);
			e_vec(u_idx) = eps*(1 + std::abs(u1(u_idx)));
			assemble_jac_ODE(time, e, u1 - e_vec, u0, M1, K1, F1);
			TMP_M_1 = (1.0 / dt)*M1 + _Theta*K1;
			local_r_1 = TMP_M_1*(u1 - e_vec);
			TMP_M_1 = (1.0 / dt) * M1 - (1. - _Theta) * K1;
			local_r_1.noalias() -= TMP_M_1 * u0;
			local_r_1.noalias() -= F1;
			assemble_jac_ODE(time, e, u1 + e_vec, u0, M2, K2, F2);
			TMP_M_2 = (1.0 / dt)*M2 + _Theta*K2;
			local_r_2 = TMP_M_2*(u1 + e_vec);
			TMP_M_2 = (1.0 / dt) * M2 - (1. - _Theta) * K2;
			local_r_2.noalias() -= TMP_M_2 * u0;
			local_r_2.noalias() -= F2;
			localJ.block(0, u_idx, u1.size(), 1) += (local_r_1 - local_r_2) / 2 / eps / (1 + std::abs(u1(u_idx)));
			//debugging--------------------------
			//std::cout << "localJ: \n";
			//std::cout << localJ << std::endl;
			//end of debugging-------------------
			//localJ=
		}
	}
protected:
	virtual void assemble_jac_ODE(const NumLib::TimeStep & /*time*/, const MeshLib::IElement &e, const MathLib::LocalVector & u1, const MathLib::LocalVector & u0, MathLib::LocalMatrix  & localM, MathLib::LocalMatrix & localK, MathLib::LocalMatrix & localF)
	{
		// current element
		FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
		// integration method 
		FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
		// now get the sizes
		const std::size_t n_dim = e.getDimension();
		const std::size_t mat_id = e.getGroupID();
		const std::size_t ele_id = e.getID();
		const std::size_t n_nodes = e.getNumberOfNodes(); // number of connecting nodes
		const std::size_t n_dof = u1.size();
		const std::size_t n_gsp = q->getNumberOfSamplingPoints(); // number of Gauss points
		std::size_t node_id(0);  // index of the node

		// get the instance of porous media class
		MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];
		// get the fluid property for gas phase
		MaterialLib::Fluid* fluid1 = Ogs6FemData::getInstance()->list_fluid[0];
		// get the fluid property for liquid phase
		MaterialLib::Fluid* fluid2 = Ogs6FemData::getInstance()->list_fluid[1];
		// get the first (light) component - hydrogen
		MaterialLib::Compound* component1 = Ogs6FemData::getInstance()->list_compound[0];
		// get the second (heavy) component - water
		MaterialLib::Compound* component2 = Ogs6FemData::getInstance()->list_compound[1];

		const NumLib::TXPosition e_pos(NumLib::TXPosition::Element, e.getID());
		double geo_area = 1.0;
		pm->geo_area->eval(e_pos, geo_area);
		const bool hasGravityEffect = _problem_coordinates.hasZ();//detect the gravity term based on the mesh structure

		if (hasGravityEffect) {
			vec_g = MathLib::LocalVector::Zero(_problem_coordinates.getDimension());
			vec_g[_problem_coordinates.getIndexOfZ()] = 9.81;
		}
		/*
		*SECONDARY VARIABLES
		*/
		// define the secondary variable value at specific one gauss point
		//Scalar value
		S_gp = MathLib::LocalVector::Zero(1, 1);
		dPC_gp = MathLib::LocalVector::Zero(1, 1);
		PC_gp = MathLib::LocalVector::Zero(1, 1);
		PG_gp = MathLib::LocalVector::Zero(1, 1);

		// Loop over each node on this element e
		// calculate the secondary variables' values on each node 
		// calculate the derivative of secondary variables based on P and X on each node
		// store the value in the vector respectively
		M = MathLib::LocalMatrix::Zero(2, 2);//method 2
		D = MathLib::LocalMatrix::Zero(2, 2);//method 2
		H = MathLib::LocalVector::Zero(2);//for gravity term
		tmp = MathLib::LocalMatrix::Zero(1, 1);
		localMass_tmp = MathLib::LocalMatrix::Zero(n_nodes, n_nodes); // for method2
		localDispersion_tmp = MathLib::LocalMatrix::Zero(n_nodes, n_nodes); //for method2
		localDispersion = MathLib::LocalMatrix::Zero(2 * n_nodes, 2 * n_nodes);
		localGravity_tmp = MathLib::LocalVector::Zero(n_nodes);//tmp matrix for gravity 
		
		
		/*
		* Loop for each Gauss Point
		*/
		Kr_L_gp = 0.0;
		Kr_G_gp = 0.0; Lambda_L = 0.0; Lambda_G = 0.0;
		R_s = 0.0;
		for (j = 0; j < n_gsp; j++)
		{
			// calc on each gauss point
			// get the sampling point
			q->getSamplingPoint(j, gp_x);
			// compute basis functions
			fe->computeBasisFunctions(gp_x);
			// get the coordinates of current location 
			fe->getRealCoordinates(real_x);
			// transfer to Gauss point position
			NumLib::TXPosition gp_pos(NumLib::TXPosition::IntegrationPoint, e.getID(), j, real_x);
			double fac = geo_area*fe->getDetJ()*q->getWeight(j);
			// get the porosity
			pm->porosity->eval(gp_pos, poro);
			// get the intrinsic permeability
			pm->permeability->eval(gp_pos, K_perm);
			// evaluate viscosity of each fluid phase
			fluid1->dynamic_viscosity->eval(gp_pos, mu_G);
			fluid2->dynamic_viscosity->eval(gp_pos, mu_L);

			//fluid1->molar_mass->eval(gp_pos, M_G);
			//fluid2->molar_mass->eval(gp_pos, M_L);

			fluid1->density->eval(gp_pos, rho_G_std);
			fluid2->density->eval(gp_pos, rho_L_std);

			// evaluate componential diffusion coefficients
			component1->molecular_diffusion->eval(gp_pos, D_G);
			component2->molecular_diffusion->eval(gp_pos, D_L);

			// evaluation of the shape function
			MathLib::LocalMatrix &Np = *fe->getBasisFunction();
			MathLib::LocalMatrix &dNp = *fe->getGradBasisFunction();



			P_gp = Np*u1.head(n_nodes);// primary variable of LIQUID_PRESSURE on gauss point
			X_gp = Np*u1.tail(n_nodes);// primary variable of MOLAR_FRACTION on gauss point
			R_s = X_gp(0, 0) / rho_g_std;
			PC_gp(0) = X_gp(0, 0) / C_h - P_gp(0, 0);
			PG_gp(0) = std::max(X_gp(0, 0) / C_h, P_gp(0, 0));
			S_gp(0) = pm->getSat_byPC(PC_gp(0));
			dPC_gp(0) = pm->get_diriv_PC(PC_gp(0));

			

			//Calc each entry of the mass matrix
			M(0, 0) = -(C_v / C_h - 1)*poro*X_gp(0, 0)*dPC_gp(0);
			M(0, 1) = poro*(1 + (C_v / C_h - 1)*(S_gp(0, 0) + (X_gp(0, 0) / C_h)*dPC_gp(0)));
			M(1, 0) = poro*rho_l_std*dPC_gp(0);
			M(1, 1) = -poro*(rho_l_std / C_h)*dPC_gp(0);
			//-------------debugging------------------------
			//std::cout << "M=" << std::endl;
			//std::cout << M << std::endl;
			//--------------end debugging-------------------

			//assembly the mass matrix
			for (ii = 0; ii < 2; ii++){
				for (jj = 0; jj < 2; jj++){
					tmp(0, 0) = M(ii, jj);
					localMass_tmp.setZero();
					fe->integrateWxN(j, tmp, localMass_tmp);
					localM.block(n_nodes*ii, n_nodes*jj, n_nodes, n_nodes) += localMass_tmp;
				}
			}
			//-------------debugging------------------------
			//std::cout << "localM=" << std::endl;
			//std::cout << localM << std::endl;
			//-------------End debugging------------------------
			//calculate the relative permeability 
			Kr_L_gp = pm->getKr_L_bySg(S_gp(0)); // change the Relative permeability calculation into porous media file			
			Kr_G_gp = pm->getKr_g_bySg(S_gp(0));
			Lambda_L = K_perm*Kr_L_gp / mu_L;
			Lambda_G = K_perm*Kr_G_gp / mu_G;


			//Calc each entry of the Laplace Matrix
			D(0, 0) = Lambda_L*X_gp(0, 0);
			D(0, 1) = Lambda_G*(C_v / (pow(C_h, 2)))*X_gp(0, 0) + poro*(1 - S_gp(0))*D_G*(F / (R_s + F));
			D(1, 0) = Lambda_L*rho_l_std;
			D(1, 1) = -poro*(1 - S_gp(0))*D_L*(F / (R_s + F) / G)*rho_l_std / rho_g_std;
			//-------------debugging------------------------
			//std::cout << "D=" << std::endl;
			//std::cout << D << std::endl;
			//--------------end debugging-------------------
			//
			for (ii = 0; ii < 2; ii++){
				for (jj = 0; jj < 2; jj++){
					tmp(0, 0) = D(ii, jj);
					localDispersion_tmp.setZero();
					fe->integrateDWxDN(j, tmp, localDispersion_tmp);
					localK.block(n_nodes*ii, n_nodes*jj, n_nodes, n_nodes) += localDispersion_tmp;
				}
			}
			//-------------debugging------------------------
			//std::cout << "localK=" << std::endl;
			//std::cout << localK << std::endl;
			//--------------end debugging-------------------
			//----------------assembly the gravity term--------------------------------
			H(0) = -Lambda_G*pow(C_v, 2)*pow((X_gp(0, 0) / C_h), 2) - Lambda_L*X_gp(0, 0)*(rho_l_std + X_gp(0, 0));
			H(1) = -Lambda_L*rho_l_std*(rho_l_std + X_gp(0, 0));
			if (hasGravityEffect) {
				// F += dNp^T * H* gz
				for (int idx = 0; idx < 2; idx++){
					tmp(0, 0) = H(idx);
					localGravity_tmp.setZero();
					localGravity_tmp = fac * dNp.transpose()*tmp*vec_g;
					localF.block(n_nodes*idx, 0, n_nodes, 1) += localGravity_tmp;

				}
			}

		}
	}
private:
    /**
      * FEM object
      */ 
    FemLib::LagrangeFeObjectContainer _feObjects;

	/**
	  * pointer to the function data class
	  */
	T_FUNCTION_DATA* _function_data; 
	MeshLib::CoordinateSystem _problem_coordinates;
	//LocalMatrixType M = LocalMatrixType::Zero(2 * n_nodes, 2 * n_nodes);
	//LocalMatrixType K = LocalMatrixType::Zero(2 * n_nodes, 2 * n_nodes);
	double _Theta;
	std::size_t i, j, ii, jj;

	double real_x[3], gp_x[3];

	MathLib::LocalVector local_r_1;
	MathLib::LocalVector local_r_2;
	MathLib::LocalMatrix M1, M2, K1, K2, F1, F2, TMP_M_1, TMP_M_2;
	MathLib::LocalVector vec_g; // gravity term 
	MathLib::LocalVector H;//this is for local gravity vector 
	MathLib::LocalVector localGravity_tmp;

	MathLib::LocalVector S_gp; // THE value of saturation on each gauss point
	MathLib::LocalVector dPC_gp;//The value of derivative PC on each gauss point
	MathLib::LocalVector PC_gp;
	MathLib::LocalVector PG_gp;

	MathLib::LocalMatrix P_gp; // the value of primary variable P on each gauss point
	MathLib::LocalMatrix X_gp;// the value of primary variable X on each gauss point



	MathLib::LocalMatrix M;//Local Mass Matrix
	MathLib::LocalMatrix D;//Local Laplace Matrix
	MathLib::LocalMatrix tmp;//method 2


	MathLib::LocalMatrix localMass_tmp; //for method2
	MathLib::LocalMatrix localDispersion_tmp; //for method2
	MathLib::LocalMatrix localDispersion;

	double Kr_L_gp, Kr_G_gp;
	double Lambda_L, Lambda_G;
	double poro, K_perm, mu_L, mu_G, D_G, D_L;
	double rho_G_std;
	double rho_L_std;
	double R_s;

	const double C_h = 1.53e-8;
	const double C_v = 7.939211048841233e-07;
	const double rho_l_std = 1000.0;
	const double rho_g_std = 0.08;
	const double F = 1.388888888888889e+03;
	const double G = 12500.0;
	
};

#endif  // end of ifndef
