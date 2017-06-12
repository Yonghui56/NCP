/**
* Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.com/LICENSE.txt
*
*
* \file NonLinearCMP_2P2C_TimeODELocalAssembler.h
*
* Created on 2015-03-19 by Yonghui HUANG
*/

#ifndef NON_LINEAR_CMP_PRESSUREFORM_TIME_ODE_LOCAL_ASSEMBLER_H
#define NON_LINEAR_CMP_PRESSUREFORM_TIME_ODE_LOCAL_ASSEMBLER_H

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Core/Integration/Integration.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "MaterialLib/PorousMedia.h"
#include "Ogs6FemData.h"




/**
* \brief Local assembly of time ODE components for linear transport in porous media
* This class is same as the MassTransportTimeODELocalAssembler class
* onlye difference is that no compound information is provided and no molecular
* diffusion is included in the assembly.
*/
template <class T, class T_NODAL_FUNCTION_SCALAR, class T_FUNCTION_DATA>
class NonLinearCMP_PressureForm_TimeODELocalAssembler : public T
{
public:
	typedef MathLib::LocalVector LocalVectorType;
	typedef MathLib::LocalMatrix LocalMatrixType;

	NonLinearCMP_PressureForm_TimeODELocalAssembler(FemLib::LagrangeFeObjectContainer* feObjects, T_FUNCTION_DATA* func_data, const MeshLib::CoordinateSystem &problem_coordinates)
		: _feObjects(*feObjects), _function_data(func_data), _problem_coordinates(problem_coordinates)
	{
		
	};
	
	virtual ~NonLinearCMP_PressureForm_TimeODELocalAssembler() {};

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
			vec_g = LocalVectorType::Zero(_problem_coordinates.getDimension());
			vec_g[_problem_coordinates.getIndexOfZ()] = -9.81;
		}
		/*
		*SECONDARY VARIABLES
		*/
		// define the secondary variable value at specific one gauss point
		//Scalar value
		S_gp = LocalVectorType::Zero(1, 1);
		dPC_gp = LocalVectorType::Zero(1, 1);
		PC_gp = LocalVectorType::Zero(1, 1);
		PG_gp = LocalVectorType::Zero(1, 1);

		// Loop over each node on this element e
		// calculate the secondary variables' values on each node 
		// calculate the derivative of secondary variables based on P and X on each node
		// store the value in the vector respectively
		M = LocalMatrixType::Zero(2, 2);//method 2
		D = LocalMatrixType::Zero(2, 2);//method 2
		H = LocalVectorType::Zero(2);//for gravity term
		tmp = LocalMatrixType::Zero(1, 1);
		localMass_tmp = LocalMatrixType::Zero(n_nodes, n_nodes); // for method2
		localDispersion_tmp = LocalMatrixType::Zero(n_nodes, n_nodes); //for method2
		localDispersion = LocalMatrixType::Zero(2 * n_nodes, 2 * n_nodes);
		localGravity_tmp = LocalVectorType::Zero(n_nodes);//tmp matrix for gravity 
		// Initialize
		_function_data->getS()->setNumberOfIntegationPoints(ele_id, n_gsp);
		_function_data->getPC()->setNumberOfIntegationPoints(ele_id, n_gsp);
		_function_data->getPG()->setNumberOfIntegationPoints(ele_id, n_gsp);
		_function_data->getdPC()->setNumberOfIntegationPoints(ele_id, n_gsp);
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
			LocalMatrixType &Np = *fe->getBasisFunction();
			LocalMatrixType &dNp = *fe->getGradBasisFunction();



			P_gp = Np*u1.head(n_nodes);// primary variable of LIQUID_PRESSURE on gauss point
			X_gp = Np*u1.tail(n_nodes);// primary variable of MOLAR_FRACTION on gauss point
			R_s = X_gp(0, 0) / rho_g_std;
			PC_gp(0) = X_gp(0, 0) / C_h - P_gp(0, 0);
			PG_gp(0) = std::max(X_gp(0, 0) / C_h, P_gp(0, 0));
			S_gp(0) = pm->getSat_byPC(PC_gp(0));
			dPC_gp(0) = pm->get_diriv_PC(PC_gp(0));
			
			_function_data->getS()->setIntegrationPointValue(ele_id, j, S_gp);
			_function_data->getPC()->setIntegrationPointValue(ele_id, j, PC_gp);
			_function_data->getPG()->setIntegrationPointValue(ele_id, j, PG_gp);
			
			//Calc each entry of the mass matrix
			M(0, 0) = -(C_v / C_h - 1)*poro*X_gp(0,0)*dPC_gp(0);
			M(0, 1) = poro*(1 + (C_v / C_h - 1)*(S_gp(0,0) + (X_gp(0,0) / C_h)*dPC_gp(0)));
			M(1, 0) = poro*rho_l_std*dPC_gp(0);
			M(1, 1) = -poro*(rho_l_std / C_h)*dPC_gp(0);
			//-------------debugging------------------------
			std::cout << "M=" << std::endl;
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
		_function_data->get_elem_M_matrix()[ele_id] = localM;
		_function_data->get_elem_K_matrix()[ele_id] = localK;

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

	LocalVectorType S_gp; // THE value of saturation on each gauss point
	LocalVectorType dPC_gp;//The value of derivative PC on each gauss point
	LocalVectorType PC_gp;
	LocalVectorType PG_gp;

	LocalMatrixType P_gp; // the value of primary variable P on each gauss point
	LocalMatrixType X_gp;// the value of primary variable X on each gauss point



	LocalMatrixType M;//Local Mass Matrix
	LocalMatrixType D;//Local Laplace Matrix
	LocalMatrixType tmp;//method 2


	LocalMatrixType localMass_tmp; //for method2
	LocalMatrixType localDispersion_tmp; //for method2
	LocalMatrixType localDispersion;
	
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