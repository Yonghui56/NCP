/**
* Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.com/LICENSE.txt
*
*
* \file NonLinearCMP_2P2C_TimeODELocalAssembler.h
*
* Created on 2015-06-19 by Yonghui HUANG
*/

#ifndef NON_LINEAR_CMP_GLOBALCOMPLEMENTARYFORM_TIME_ODE_LOCAL_ASSEMBLER_H
#define NON_LINEAR_CMP_GLOBALCOMPLEMENTARYFORM_TIME_ODE_LOCAL_ASSEMBLER_H

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Core/Integration/Integration.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "MaterialLib/PorousMedia.h"
#include "Ogs6FemData.h"
#include "EOS_GlobalComplementaryForm.h"




/**
* \brief Local assembly of time ODE components for linear transport in porous media
* This class is same as the MassTransportTimeODELocalAssembler class
* onlye difference is that no compound information is provided and no molecular
* diffusion is included in the assembly.
*/
template <class T, class T_NODAL_FUNCTION_SCALAR, class T_FUNCTION_DATA>
class NonLinearCMP_GlobalComplementaryForm_TimeODELocalAssembler : public T
{
public:
	typedef MathLib::LocalVector LocalVectorType;
	typedef MathLib::LocalMatrix LocalMatrixType;

	NonLinearCMP_GlobalComplementaryForm_TimeODELocalAssembler(FemLib::LagrangeFeObjectContainer* feObjects, T_FUNCTION_DATA* func_data, const MeshLib::CoordinateSystem &problem_coordinates)
		: _feObjects(*feObjects), _function_data(func_data), _problem_coordinates(problem_coordinates)
	{
		_EOS = new EOS_GlobalComplementaryForm();
	};
	
	virtual ~NonLinearCMP_GlobalComplementaryForm_TimeODELocalAssembler() {
		BaseLib::releaseObject(_EOS);
	};
	T_FUNCTION_DATA* get_function_data(void)
	{
		return _function_data;
	}

protected:
	
	virtual void assembleODE(const NumLib::TimeStep & /*time*/, const MeshLib::IElement &e, const LocalVectorType & u1, const LocalVectorType & u0, LocalMatrixType & localM, LocalMatrixType & localK, LocalVectorType & localF)
	{
		// -------------------------------------------------
		// current element
		FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
		// integration method 
		FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
		// now get the sizes
		const std::size_t n_dim = e.getDimension();
		const std::size_t mat_id = e.getGroupID();
		const std::size_t ele_id = e.getID();
		const std::size_t n_nodes = e.getNumberOfNodes(); // number of connecting nodes
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
		const bool hasGravityEffect =  _problem_coordinates.hasZ();//detect the gravity term based on the mesh structure

		if (hasGravityEffect) {
			vec_g = LocalVectorType::Zero(_problem_coordinates.getDimension());
			vec_g[_problem_coordinates.getIndexOfY()] = 9.81;
		}
		

		M = LocalMatrixType::Zero(2, 3);//method 2
		D = LocalMatrixType::Zero(2, 3);//method 2
		H = LocalVectorType::Zero(2);//for gravity term
		tmp = LocalMatrixType::Zero(1, 1);
		Input = LocalVectorType::Zero(3);
		localMass_tmp = LocalMatrixType::Zero(n_nodes, n_nodes); // for method2
		localDispersion_tmp = LocalMatrixType::Zero(n_nodes, n_nodes); //for method2
		localDispersion = LocalMatrixType::Zero(n_nodes, n_nodes);
		LocalJacb = LocalMatrixType::Zero(3*n_nodes, 3 * n_nodes);
		localGravity_tmp = LocalVectorType::Zero(n_nodes);//tmp vect for gravity 
		localRes = LocalVectorType::Zero(n_nodes);//tmp vect for gravity 
		/*
		* Loop for each Gauss Point
		*/
		Kr_L_gp = 0.0;
		Kr_G_gp = 0.0; lambda_L = 0.0; lambda_G = 0.0;
		Lambda_h = 0.0;
		//-------------------Begin assembly on each gauss point---------------------------
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
			LocalMatrixType &Np = *fe->getBasisFunction();
			LocalMatrixType &dNp = *fe->getGradBasisFunction();


			P_gp = Np*u1.head(n_nodes);//gas phase pressure
			X_gp = Np*u1.block(n_nodes, 0, n_nodes, 1);// here we assume to be x_L^h
			S_gp = Np*u1.tail(n_nodes);//gas phase saturation
			Input(0) = P_gp(0, 0);
			Input(1) = X_gp(0, 0);
			Input(2) = S_gp(0, 0);

			double const RT = R*temperature;
			double const PC_gp = pm->getPc_bySat(S_gp(0, 0));
			double const dPC_dSg_gp = pm->Deriv_dPCdS(S_gp(0, 0));
            double const N_G = P_gp(0, 0) / RT;
            double const N_L = rho_l_std / M_L;
            double const x_G_h = 1;// X_gp(0, 0)*N_L / P_gp(0, 0) / H;
            double const d_x_G_h_dX = 0.0;// N_L / P_gp(0, 0) / H;
            double const d_x_G_h_dPG = 0.0;// -X_gp(0, 0)*N_L / P_gp(0, 0) / P_gp(0, 0) / H;

			
            M(0, 0) = poro*(x_G_h*S_gp(0, 0)*(1 / RT) + N_G*S_gp(0, 0)*d_x_G_h_dPG);
            M(0, 1) = poro*(N_G*S_gp(0, 0)*d_x_G_h_dX+ N_L*(1- S_gp(0, 0)));
            M(0, 2) = poro*(N_G* x_G_h - N_L*X_gp(0, 0));
			M(1, 0) = 0.0;
			M(1, 1) = 0.0;
			M(1, 2) = -poro*N_L;//
			
			//-------------debugging------------------------
			//std::cout << "M=" << std::endl;
			// std::cout << M << std::endl;
			//--------------end debugging-------------------

			//assembly the mass matrix
			for (ii = 0; ii < 2; ii++){
				for (jj = 0; jj < 3; jj++){
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
			Kr_L_gp = pm->getKr_L_bySg(S_gp(0, 0)); // change the Relative permeability calculation into porous media file			
			Kr_G_gp = pm->getKr_g_bySg(S_gp(0, 0));
			lambda_L = K_perm*Kr_L_gp / mu_L;
			lambda_G = K_perm*Kr_G_gp / mu_G;
			
			//Calc each entry of the Laplace Matrix
            D(0, 0) = N_L*X_gp(0, 0)*lambda_L + N_G*x_G_h*lambda_G;
			D(0, 1) = poro*(1 - S_gp(0, 0))*D_L*N_L;//
			D(0, 2) = -N_L*X_gp(0, 0)*lambda_L*dPC_dSg_gp;
			D(1, 0) = N_L*lambda_L;
			D(1, 1) = -poro*(1 - S_gp(0, 0))*D_L*N_L;// 
            D(1, 2) = -N_L*lambda_L*dPC_dSg_gp;
			//-------------debugging------------------------
			//std::cout << "D=" << std::endl;
			// std::cout << D << std::endl;
			//--------------end debugging-------------------
			//
			for (ii = 0; ii < 2; ii++){
				for (jj = 0; jj < 3; jj++){
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
            H(0) = 0.0;
            H(1) = 0.0;
			//-------------debugging------------------------
			//std::cout << "H=" << std::endl;
			//std::cout << H << std::endl;
			//--------------end debugging-------------------
			if (hasGravityEffect) {
				// F += dNp^T * H* gz
				
				for (int idx = 0; idx < 2; idx++){
					tmp(0, 0) = H(idx);
					localGravity_tmp.setZero();
					//fe->integrateDWxvec_g(j, tmp, localGravity_tmp, vec_g);
					localGravity_tmp = fac * dNp.transpose()*tmp(0,0)*vec_g;
					localF.block(n_nodes*idx, 0, n_nodes, 1) += localGravity_tmp;

				}
			}
			//----------------end assembly--------------------------------
			//std::cout << "H=" << std::endl;
			//std::cout << localF << std::endl;

		}
		//----------------end assemby gaus---------------------------------
		//----------------calc jacobian for complementary contrains--------------

		//std::cout << "H=" << std::endl;
		//std::cout << LocalJacb<< std::endl;
		_function_data->get_elem_M_matrix()[ele_id] = localM;
		_function_data->get_elem_K_matrix()[ele_id] = localK;
		
		//_function_data->get_elem_J_matrix()[ele_id] = LocalJacb;
	}

private:
	FemLib::LagrangeFeObjectContainer _feObjects;
	MeshLib::CoordinateSystem _problem_coordinates;
	/**
	  * pointer to the function data class
	  */
	T_FUNCTION_DATA* _function_data;
	std::size_t i, j, ii, jj;

	double real_x[3], gp_x[3];

	LocalVectorType vec_g; // gravity term 
	LocalVectorType H;//this is for local gravity vector 
	LocalVectorType localGravity_tmp;
	LocalVectorType localRes;

	LocalVectorType PC;
	LocalVectorType rho_L_h;
	LocalVectorType rho_G_h;
	LocalVectorType dPC_dSg;

	LocalVectorType Input;

	
	LocalMatrixType rho_L_h_gp;
	LocalMatrixType rho_G_h_gp;
	LocalMatrixType PC_gp;
	LocalMatrixType dPC_dSg_gp;
	//LocalMatrixType PG_gp;
	//LocalMatrixType PL_gp;
	LocalMatrixType drho_L_h_dXgp;
	LocalMatrixType drho_L_h_dPgp;
	LocalMatrixType drho_L_h_dSgp;
	//primary variable on gauss point
	LocalMatrixType P_gp;
	LocalMatrixType X_gp;
	LocalMatrixType S_gp;


	LocalMatrixType M;//Local Mass Matrix
	LocalMatrixType D;//Local Laplace Matrix
	LocalMatrixType tmp;//method 2


	LocalMatrixType localMass_tmp; //for method2
	LocalMatrixType localDispersion_tmp; //for method2
	LocalMatrixType localDispersion;
	LocalMatrixType matN;
	LocalMatrixType LocalJacb;

    double const temperature = 295.15;
    double const R = 8.314;
	double Kr_L_gp, Kr_G_gp;
	double lambda_L, lambda_G;
	double poro, K_perm, mu_L, mu_G, D_G, D_L;
	double rho_G_std;
	double rho_L_std;
	double Lambda_h;
	double Var_a;
	const double M_L = 0.018;
	const double M_G = 0.002;

	double RHO_L;// RHO_L=rho_L_std+rho_L_h
	double RHO_G;// RHO_G=rho_G_h+rho_G_w
	const double rho_l_std = 1000.0;
	EOS_GlobalComplementaryForm * _EOS;
	std::size_t isinf;
};



#endif  // end of ifndef