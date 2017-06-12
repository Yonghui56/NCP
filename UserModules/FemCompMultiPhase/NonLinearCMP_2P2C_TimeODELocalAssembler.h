/**
* Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.com/LICENSE.txt
*
*
* \file NonLinearCMP_2P2C_TimeODELocalAssembler.h
*
* Created on 2014-05-19 by Yonghui HUANG and Haibing SHAO
*/

#ifndef NON_LINEAR_CMP_2P2C_TIME_ODE_LOCAL_ASSEMBLER_H
#define NON_LINEAR_CMP_2P2C_TIME_ODE_LOCAL_ASSEMBLER_H

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Core/Integration/Integration.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "MaterialLib/PorousMedia.h"
#include "Ogs6FemData.h"

#include "EOS_H2_H2O_ISOT_3x_Complementary.h"
#include "AbstractEOS.h"

/**
* \brief Local assembly of time ODE components for linear transport in porous media
* This class is same as the MassTransportTimeODELocalAssembler class
* onlye difference is that no compound information is provided and no molecular
* diffusion is included in the assembly.
*/
template <class T, class T_NODAL_FUNCTION_SCALAR, class T_FUNCTION_DATA>
class NonLinearCMP_2P2C_TimeODELocalAssembler : public T
{
public:
	typedef MathLib::LocalVector LocalVectorType;
	typedef MathLib::LocalMatrix LocalMatrixType;

	NonLinearCMP_2P2C_TimeODELocalAssembler(FemLib::LagrangeFeObjectContainer* feObjects, T_FUNCTION_DATA* func_data)
		: _feObjects(*feObjects), _function_data(func_data)
	{
		//_EOS = new EOS_H2_H2O_ISOT_3x_Complementary();
	};
	
	virtual ~NonLinearCMP_2P2C_TimeODELocalAssembler() {};

protected:
	
	virtual void assembleODE(const NumLib::TimeStep & /*time*/, const MeshLib::IElement &e, const LocalVectorType & u1, const LocalVectorType & u0, LocalMatrixType & localM, LocalMatrixType & localK, LocalVectorType & localF)
	{

		// -------------------------------------------------
        // current element
		FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
        // integration method 
        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        // now get the sizes
		const std::size_t n_dim   = e.getDimension();
		const std::size_t mat_id  = e.getGroupID();
		const std::size_t ele_id  = e.getID();
        const std::size_t n_nodes = e.getNumberOfNodes(); // number of connecting nodes
        const std::size_t n_gsp   = q->getNumberOfSamplingPoints(); // number of Gauss points
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

        
		PC = LocalVectorType::Zero(n_nodes);
		PL = LocalVectorType::Zero(n_nodes);
		PG = LocalVectorType::Zero(n_nodes);
		X_m = LocalVectorType::Zero(n_nodes);
		X_M = LocalVectorType::Zero(n_nodes);
		N_G = LocalVectorType::Zero(n_nodes);
		N_L = LocalVectorType::Zero(n_nodes);
		S = LocalVectorType::Zero(n_nodes);
		Mat_Sec_Derv = LocalMatrixType::Zero(8, 2);
		dPCdP = LocalVectorType::Zero(n_nodes);
		dPCdX = LocalVectorType::Zero(n_nodes);
		dPGdP = LocalVectorType::Zero(n_nodes);
		dPGdX = LocalVectorType::Zero(n_nodes);
		dPLdP = LocalVectorType::Zero(n_nodes);
		dPLdX = LocalVectorType::Zero(n_nodes);
		dX_mdP = LocalVectorType::Zero(n_nodes);
		dX_mdX = LocalVectorType::Zero(n_nodes);
		dX_MdP = LocalVectorType::Zero(n_nodes);
		dX_MdX = LocalVectorType::Zero(n_nodes);
		dN_GdP = LocalVectorType::Zero(n_nodes);
		dN_GdX = LocalVectorType::Zero(n_nodes);
		dN_LdP = LocalVectorType::Zero(n_nodes);
		dN_LdX = LocalVectorType::Zero(n_nodes);
		dSdP = LocalVectorType::Zero(n_nodes);
		dSdX = LocalVectorType::Zero(n_nodes);

        // define the secondary variable value at specific one point
        PC_gp = LocalMatrixType::Zero(1, 1);
        X_m_gp = LocalMatrixType::Zero(1, 1);
        X_M_gp = LocalMatrixType::Zero(1, 1);
        N_G_gp = LocalMatrixType::Zero(1, 1);
        N_L_gp = LocalMatrixType::Zero(1, 1);
        S_gp = LocalMatrixType::Zero(1, 1);
        P_gp = LocalMatrixType::Zero(1, 1);
        X_gp = LocalMatrixType::Zero(1, 1);


		// Loop over each node on this element e
		// calculate the secondary variables' values on each node 
		// calculate the derivative of secondary variables based on P and X on each node
		// store the value in the vector respectively
		for (i = 0; i < n_nodes; i++)
		{
			node_id = e.getNodeID(i); // e is the current element
			
			PC[i] = _function_data->get_PC()->getValue(node_id);
			PL[i] = _function_data->get_PL()->getValue(node_id);
			PG[i] = _function_data->get_PG()->getValue(node_id);
			X_m[i] = _function_data->get_X_m()->getValue(node_id);
			X_M[i] = _function_data->get_X_M()->getValue(node_id);
			N_G[i] = _function_data->getN_G()->getValue(node_id);
			N_L[i] = _function_data->getN_L()->getValue(node_id);
			S[i] = _function_data->getS()->getValue(node_id);
			Mat_Sec_Derv = _function_data->get_mat_secDer()->getValue(node_id);

			dSdP[i] = Mat_Sec_Derv(0, 0);
			dSdX[i] = Mat_Sec_Derv(0, 1);

			dX_mdP[i] = Mat_Sec_Derv(1, 0);
			dX_mdX[i] = Mat_Sec_Derv(1, 1);

			dX_MdP[i] = Mat_Sec_Derv(2, 0);
			dX_MdX[i] = Mat_Sec_Derv(2, 1);

			dPCdP[i] = Mat_Sec_Derv(3, 0);
			dPCdX[i] = Mat_Sec_Derv(3, 1);

			dPLdP[i] = Mat_Sec_Derv(4, 0);
			dPLdX[i] = Mat_Sec_Derv(4, 1);
			
			dPGdP[i] = Mat_Sec_Derv(5, 0);
			dPGdX[i] = Mat_Sec_Derv(5, 1);

			dN_LdP[i] = Mat_Sec_Derv(6, 0);
			dN_LdX[i] = Mat_Sec_Derv(6, 1);

			dN_GdP[i] = Mat_Sec_Derv(7, 0);
			dN_GdX[i] = Mat_Sec_Derv(7, 1);
		}  // end of for i

		/*
        std::cout << "Mat_Sec_Derv=" << std::endl;
		std::cout << Mat_Sec_Derv << std::endl;
        */

        poro = 0.0; K_perm = 0.0; mu_L = 0.0; mu_G = 0.0; D_G = 0.0; D_L = 0.0; Lambda_L = 0.0; Lambda_G = 0.0;
        Kr_L_gp = 0.0; Kr_G_gp = 0.0;

		C_h = Hen*M_G;
		C_v = M_G / R/T;
		rho_L_std = 0.0;
		rho_G_std = 0.0;
		gamma = 0.0;

        // Initialize the gradient 
        dSgdX_gp = LocalMatrixType::Zero(1, 1);
        dSgdP_gp = LocalMatrixType::Zero(1, 1);
        dX_mdP_Liq_gp = LocalMatrixType::Zero(1, 1); //dX_mdP
        dX_MdP_Gas_gp = LocalMatrixType::Zero(1, 1);//dX_MdP
        dX_MdX_Gas_gp = LocalMatrixType::Zero(1, 1);//dX_MdS
        dX_mdX_Liq_gp = LocalMatrixType::Zero(1, 1);//dX_mdS
        dN_GdP_gp = LocalMatrixType::Zero(1, 1);//dN_GdP
        dN_GdX_gp = LocalMatrixType::Zero(1, 1);//dN_GdS
        dN_LdP_gp = LocalMatrixType::Zero(1, 1);
        dN_LdX_gp = LocalMatrixType::Zero(1, 1);
        dPCdP_gp = LocalMatrixType::Zero(1, 1);
        dPCdX_gp = LocalMatrixType::Zero(1, 1);
		dPGdP_gp = LocalMatrixType::Zero(1, 1);
		dPGdX_gp = LocalMatrixType::Zero(1, 1);
		dPLdP_gp = LocalMatrixType::Zero(1, 1);
		dPLdX_gp = LocalMatrixType::Zero(1, 1);
        M = LocalMatrixType::Zero(2, 2);//method 2
        D = LocalMatrixType::Zero(2, 2);//method 2
        tmp = LocalMatrixType::Zero(1, 1);//method2

        B11_dP = 0.0; B12_dP = 0.0; B13_dP = 0.0; B14_dP = 0.0;
        B11_dX = 0.0; B12_dX = 0.0; B13_dX = 0.0; B14_dX = 0.0;
        B21_dP = 0.0; B22_dP = 0.0; B23_dP = 0.0; B24_dP = 0.0;
        B21_dX = 0.0; B22_dX = 0.0; B23_dX = 0.0; B24_dX = 0.0;

		localMass_tmp = LocalMatrixType::Zero(n_nodes, n_nodes); // for method2
		localDispersion_tmp = LocalMatrixType::Zero(n_nodes, n_nodes); //for method2
		localDispersion = LocalMatrixType::Zero(2*n_nodes, 2*n_nodes);


		// loop for each gauss point
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

			PC_gp = Np*PC;//one scalar value
			X_m_gp = Np*X_m;
			X_M_gp = Np*X_M;
			N_G_gp = Np*N_G;
			N_L_gp = Np*N_L;
			S_gp = Np*S;
			//get the value of primary variable on this gauss point
			P_gp = Np*u1.head(3);
			X_gp = Np*u1.tail(3);
			// get the derivative value on this gauss point
			dSgdP_gp = Np*dSdP;
			dSgdX_gp = Np*dSdX;
			dX_mdP_Liq_gp = Np*dX_mdP;
			dX_mdX_Liq_gp = Np*dX_mdX;
			dX_MdP_Gas_gp = Np*dX_MdP;
			dX_MdX_Gas_gp = Np*dX_MdX;
			dN_LdP_gp = Np*dN_LdP;
			dN_LdX_gp = Np*dN_LdX;
			dN_GdP_gp = Np*dN_GdP;
			dN_GdX_gp = Np*dN_GdX;
			dPCdP_gp = Np*dPCdP;
			dPCdX_gp = Np*dPCdX;
			dPGdP_gp = Np*dPGdP;
			dPGdX_gp = Np*dPGdX;
			dPLdP_gp = Np*dPLdP;
			dPLdX_gp = Np*dPLdX;

			// Define three character function for calculate the gradient
			xi_gp = pm->Characteristic_Function_xi(S_gp(0, 0));
			//alpha_gp = pm->Characteristic_Function_alpha(S_gp(0, 0));
			//beta_gp = pm->Characteristic_Function_Beta(S_gp(0, 0));
            

			M(0, 1) = poro*(S_gp(0, 0)*N_G_gp(0, 0) + (1 - S_gp(0, 0))*N_L_gp(0, 0) + X_gp(0, 0)*(N_G_gp(0, 0)*dSgdX_gp(0, 0)*xi_gp + S_gp(0, 0)*dN_GdX_gp(0, 0)
				+ (1 - S_gp(0, 0))*dN_LdX_gp(0, 0) - xi_gp*dSgdX_gp(0, 0)*N_L_gp(0, 0)));//M12(0, 0))
			M(1,1)= poro*(-( S_gp(0, 0)*N_G_gp(0, 0)+(1 - S_gp(0, 0)) *N_L_gp(0, 0))
				+ (1 - X_gp(0, 0))*(N_G_gp(0, 0)*dSgdX_gp(0, 0)*xi_gp + S_gp(0, 0)*dN_GdX_gp(0, 0) + (1 - S_gp(0, 0))*dN_LdX_gp(0, 0) - xi_gp* dSgdX_gp(0, 0)*N_L_gp(0, 0)));//)M22(0, 0) //modify minus
			M(0, 0) = poro*(X_gp(0, 0)*(xi_gp*dSgdP_gp(0, 0)*N_G_gp(0, 0) + S_gp(0, 0)*dN_GdP_gp(0, 0) + (1 - S_gp(0, 0))*dN_LdP_gp(0, 0) - xi_gp*dSgdP_gp(0, 0)*N_L_gp(0, 0)));//M11(0, 0) 
			M(1, 0) = poro*((1 - X_gp(0, 0))*(xi_gp*dSgdP_gp(0, 0)*N_G_gp(0, 0) + S_gp(0, 0)*dN_GdP_gp(0, 0) + (1 - S_gp(0, 0))*dN_LdP_gp(0, 0) - xi_gp*dSgdP_gp(0, 0)*N_L_gp(0, 0)));//M21(0, 0)
			//-------------debugging------------------------
			//std::cout << "M=" << std::endl;
			//std::cout << M << std::endl;
			//--------------end debugging-------------------

			//LOOP FOR ASSEMBLYing the matrix			
			for (ii = 0; ii < 2; ii++){
				for (jj = 0; jj < 2; jj++){
					tmp(0, 0) = M(ii, jj);
					localMass_tmp.setZero(); 
					fe->integrateWxN(j, tmp, localMass_tmp);
					localM.block(n_nodes*ii, n_nodes*jj, n_nodes, n_nodes) += localMass_tmp;
				}
			}
			
			//calculate the relative permeability 
			Kr_L_gp= pm->getKr_L_bySg(S_gp(0, 0)); // change the Relative permeability calculation into porous media file			
			Kr_G_gp = pm->getKr_g_bySg(S_gp(0, 0));
			Lambda_L = K_perm*Kr_L_gp / mu_L;
			Lambda_G = K_perm*Kr_G_gp / mu_G;
			gamma = pm->getweightedFunc(S_gp(0, 0));

			//diffusion term----Eq.1
			//dp
			B11_dP = N_G_gp(0, 0)*X_M_gp(0, 0)*Lambda_G*dPGdP_gp(0, 0);
				//+ N_G_gp(0, 0)*X_M_gp(0, 0)*Lambda_G*(1 - gamma)*dPCdP_gp(0, 0)
				//- N_G_gp(0, 0)*X_M_gp(0, 0)*Lambda_G*PC_gp(0, 0)*xi_gp*(2 - 2 * S_gp(0, 0))*dSgdP_gp(0, 0);//NG*X_M*K*Krg/mu_g
			B12_dP = N_L_gp(0, 0)*X_m_gp(0, 0)*Lambda_L*dPLdP_gp(0, 0);
				//- N_L_gp(0, 0)*X_m_gp(0, 0)*Lambda_L *gamma*dPCdP_gp(0, 0)
				//- N_L_gp(0, 0)*X_m_gp(0, 0)*Lambda_L*PC_gp(0, 0)*xi_gp*(2 - 2 * S_gp(0, 0))*dSgdP_gp(0, 0);
			B13_dP = N_G_gp(0, 0)*D_G*S_gp(0, 0)*poro*dX_MdP_Gas_gp(0, 0)*xi_gp;
			B14_dP = N_L_gp(0, 0)*D_L*(1 - S_gp(0, 0))*poro*dX_mdP_Liq_gp(0, 0)*xi_gp;
			D(0,0)= B11_dP + B12_dP + B13_dP + B14_dP;//D11(0, 0)

			//dx
			B11_dX = N_G_gp(0, 0)*X_M_gp(0, 0)*Lambda_G*dPGdX_gp(0, 0);
				//(1 - gamma)*dPCdX_gp(0, 0)
				//+ N_G_gp(0, 0)*X_M_gp(0, 0)*Lambda_G*PC_gp(0, 0)*xi_gp*(2 - 2 * S_gp(0, 0))*dSgdX_gp(0, 0);
			B12_dX = N_L_gp(0, 0)*X_m_gp(0, 0)*Lambda_L*dPLdX_gp(0, 0);
				//-N_L_gp(0, 0)*X_m_gp(0, 0)*Lambda_L*gamma*dPCdX_gp(0, 0)
				//- N_L_gp(0, 0)*X_m_gp(0, 0)*Lambda_L *PC_gp(0, 0)*xi_gp*(2 - 2 * S_gp(0, 0))*dSgdX_gp(0, 0);
			B13_dX = N_G_gp(0, 0)*D_G*S_gp(0, 0)*poro*(1 - xi_gp + xi_gp*dX_MdX_Gas_gp(0, 0));
			B14_dX = N_L_gp(0, 0)*D_L*(1 - S_gp(0, 0))*poro*(1 - xi_gp + xi_gp*dX_mdX_Liq_gp(0, 0));
			D(0,1) = B11_dX + B12_dX + B13_dX + B14_dX;//D12(0, 0)

			//Diffusion term ----Eq.2
			//dp
			B21_dP = N_G_gp(0, 0)*(1 - X_M_gp(0, 0))*Lambda_G*dPGdP_gp(0, 0);
				//+ N_G_gp(0, 0)*(1 - X_M_gp(0, 0))*Lambda_G*(1 - gamma)*dPCdP_gp(0, 0)
				//- N_G_gp(0, 0)*(1 - X_M_gp(0, 0))*Lambda_G*PC_gp(0, 0)*xi_gp*(2 - 2 * S_gp(0, 0))*dSgdP_gp(0, 0);//NG*X_M*K*Krg/mu_g
			B22_dP = N_L_gp(0, 0)*(1 - X_m_gp(0, 0))*Lambda_L*dPLdP_gp(0, 0);
				//- N_L_gp(0, 0)*(1-X_m_gp(0, 0))*Lambda_L *gamma*dPCdP_gp(0, 0)
				//- N_L_gp(0, 0)*(1-X_m_gp(0, 0))*Lambda_L *PC_gp(0, 0)*xi_gp*(2 - 2 * S_gp(0, 0))*dSgdP_gp(0, 0);

			B23_dP = -N_G_gp(0, 0)*D_G*S_gp(0, 0)*poro*dX_MdP_Gas_gp(0, 0)*xi_gp;
			B24_dP = -N_L_gp(0, 0)*D_L*(1 - S_gp(0, 0))*poro*dX_mdP_Liq_gp(0, 0)*xi_gp;
			D(1,0) = B21_dP + B22_dP + B23_dP + B24_dP;// D21(0, 0)

			//dx
			B21_dX = N_G_gp(0, 0)*(1 - X_M_gp(0, 0))*Lambda_G*dPGdX_gp(0, 0);
				//(1 - gamma)*dPCdX_gp(0, 0)
				//+ N_G_gp(0, 0)*(1 - X_M_gp(0, 0))*Lambda_G*PC_gp(0, 0)*xi_gp*(2 - 2 * S_gp(0, 0))*dSgdX_gp(0, 0);
			//N_G_gp(0, 0)*(1 - X_M_gp(0, 0))*Lambda_G*dPGdX_gp(0, 0);
			B22_dX = N_L_gp(0, 0)*(1 - X_m_gp(0, 0))*Lambda_L*dPLdX_gp(0, 0);
				//gamma*dPCdX_gp(0, 0)
				//- N_L_gp(0, 0)*X_m_gp(0, 0)*Lambda_L *PC_gp(0, 0)*xi_gp*(2 - 2 * S_gp(0, 0))*dSgdX_gp(0, 0);

			B23_dX = -N_G_gp(0, 0)*D_G*S_gp(0, 0)*poro*(1 - xi_gp + xi_gp*dX_MdX_Gas_gp(0, 0));
			B24_dX = -N_L_gp(0, 0)*D_L*(1 - S_gp(0, 0))*poro*(1 - xi_gp + xi_gp*dX_mdX_Liq_gp(0, 0));

			D(1,1)= B21_dX + B22_dX + B23_dX + B24_dX;//D22(0, 0) 
			
            //-------------debugging------------------------
			//std::cout << "D=" << std::endl;
			//std::cout << D << std::endl;
			//--------------end debugging-------------------
			
			for (ii = 0; ii < 2; ii++){
				for (jj = 0; jj < 2; jj++){
					tmp(0, 0) = D(ii, jj);
					localDispersion_tmp.setZero(); 
					fe->integrateDWxDN(j, tmp, localDispersion_tmp); 
					localK.block(n_nodes*ii, n_nodes*jj, n_nodes, n_nodes) += localDispersion_tmp; 
				}
			}

			//----------debugging---------------------------
			//std::cout << "dshape=" << std::endl;
			//std::cout << *fe->getGradBasisFunction() << std::endl;
			//std::cout << "localK" << std::endl;
			//std::cout << localK << std::endl;
			//----------debugging end-----------------------

		} // end of loop over all gauss points for j
		
		// debugging--------------------------
		//std::cout << "localM: \n";
		//std::cout << localM << std::endl;
		// end of debugging-------------------	

		// debugging--------------------------
		//std::cout << "localK: \n";
		//std::cout << localK << std::endl;
		// end of debugging-------------------		

        // storing the local matrix
		_function_data->get_elem_M_matrix()[ele_id] = localM;
		_function_data->get_elem_K_matrix()[ele_id] = localK;
	}

private:
	FemLib::LagrangeFeObjectContainer _feObjects;

	/**
	  * pointer to the function data class
	  */
	T_FUNCTION_DATA* _function_data;

    /**
      * index numbers
      */
    std::size_t i, j, ii, jj;  

    /**
      * TODO
      */
    double real_x[3], gp_x[3];

    /**
      * TODO
      */
    std::size_t xi_gp;
    std::size_t alpha_gp;
    std::size_t beta_gp;

    /**
      * TODO
      */
    LocalVectorType PC;
    LocalVectorType PL;
    LocalVectorType PG; 
    LocalVectorType X_m; 
    LocalVectorType X_M; 
    LocalVectorType N_G; 
    LocalVectorType N_L; 
    LocalVectorType S; 
    LocalMatrixType Mat_Sec_Derv; 
    LocalVectorType dPCdP; 
    LocalVectorType dPCdX; 
    LocalVectorType dPLdP; 
    LocalVectorType dPLdX; 
	LocalVectorType dPGdP;
	LocalVectorType dPGdX;
    LocalVectorType dX_mdP;
    LocalVectorType dX_mdX;
    LocalVectorType dX_MdP;
    LocalVectorType dX_MdX;
    LocalVectorType dN_GdP;
    LocalVectorType dN_GdX;
    LocalVectorType dN_LdP;
    LocalVectorType dN_LdX;
    LocalVectorType dSdP;
    LocalVectorType dSdX;

    // define the secondary variable value at specific one point
    LocalMatrixType PC_gp; 
    LocalMatrixType X_m_gp; 
    LocalMatrixType X_M_gp; 
    LocalMatrixType N_G_gp; 
    LocalMatrixType N_L_gp; 
    LocalMatrixType S_gp; 
    LocalMatrixType P_gp; 
    LocalMatrixType X_gp; 


    /**
      * TODO
      */
    double poro, K_perm, mu_L, mu_G, D_G, D_L, Lambda_L, Lambda_G;
    double Kr_L_gp, Kr_G_gp;
	double rho_L_std;
	double rho_G_std;
	double C_h;
	double C_v;
	const double R = 8.314;
	const double T = 303;
	const double M_L = 0.01;
	const double M_G = 0.002;
	const double Hen = 7.65e-6;// also assume Henry constant
    double gamma;

    // Initialize the gradient 
    LocalMatrixType dSgdX_gp;
    LocalMatrixType dSgdP_gp;
    LocalMatrixType dX_mdP_Liq_gp; //dX_mdP
    LocalMatrixType  dX_MdP_Gas_gp;//dX_MdP
    LocalMatrixType  dX_MdX_Gas_gp;//dX_MdS
    LocalMatrixType dX_mdX_Liq_gp;//dX_mdS
    LocalMatrixType dN_GdP_gp;//dN_GdP
    LocalMatrixType dN_GdX_gp;//dN_GdS
    LocalMatrixType dN_LdP_gp;
    LocalMatrixType dN_LdX_gp;
    LocalMatrixType dPCdP_gp;
    LocalMatrixType dPCdX_gp;
	LocalMatrixType dPGdP_gp;
	LocalMatrixType dPGdX_gp;
	LocalMatrixType dPLdP_gp;
	LocalMatrixType dPLdX_gp;
    LocalMatrixType M;//method 2
    LocalMatrixType D;//method 2
    LocalMatrixType tmp;//method 2

    LocalMatrixType localMass_tmp; //for method2
    LocalMatrixType localDispersion_tmp; //for method2
    LocalMatrixType localDispersion;

    /**
      * TODO
      */
    double B11_dP, B12_dP, B13_dP, B14_dP;
    double B11_dX, B12_dX, B13_dX, B14_dX;
    double B21_dP, B22_dP, B23_dP, B24_dP;
    double B21_dX, B22_dX, B23_dX, B24_dX;

};



#endif  // end of ifndef