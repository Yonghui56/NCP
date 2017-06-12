/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Ogs5ToOgs6.cpp
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#include "Ogs5ToOgs6.h"

#include "logog.hpp"

#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TimeStepping/TimeStepFunction.h"
#include "NumLib/TimeStepping/TimeStepFunctionNewtonAdaptive.h"
#include "MaterialLib/Fluid.h"
#include "MaterialLib/Fracture.h"
#include "SolutionLib/Fem/FemDirichletBC.h"
#include "SolutionLib/Fem/FemNeumannBC.h"
#include "GeoProcessBuilder.h"
#include "ChemLib/chemReactionKin.h"
#include "ChemLib/chemReductionKin.h"
#include "ChemLib/chemReductionGIA.h"
#include "ChemLib/chemReactionEqMob.h"
#include "ChemLib/chemReactionEqSorp.h"
#include "ChemLib/chemReactionEqMin.h"
#include "ChemLib/chemActivityModelAqDBH.h"
#include "ChemLib/chemActivityModelAqDVS.h"
#include "ChemLib/chemActivityModelUnity.h"
#include "ChemLib/chemEqReactSys.h"
#include "ChemLib/chemKinReactSys.h"

using namespace ogs5;

namespace ogs6
{

namespace Ogs5ToOgs6
{

void convertTimeStepping(const std::vector<CTimeDiscretization*> &td_vector, std::vector<NumLib::ITimeStepFunction*> &tf_vector)
{
    bool done_shared = false;
    for (size_t i=0; i<td_vector.size(); i++)
    {
        const CTimeDiscretization* td = td_vector[i];
        if (td->time_control_name == "NEWTON_ADAPTIVE" )
        {
            NumLib::ITimeStepFunction* tf = new NumLib::TimeStepFunctionNewtonAdaptive(td->time_start, 
                                                                                       td->time_end, 
                                                                                       td->min_time_step, 
                                                                                       td->max_time_step, 
                                                                                       td->time_adapt_tim_vector, 
                                                                                       td->time_adapt_coe_vector);
            tf_vector.push_back(tf);
            // if (!td->time_independence)
            //     done_shared = true;
        }
        // else if (td->time_independence || !done_shared) {
        else {
            NumLib::ITimeStepFunction* tf = new NumLib::TimeStepFunctionVector(td->time_start, td->time_end, td->time_step_vector);
            tf_vector.push_back(tf);
            // if (!td->time_independence)
            //     done_shared = true;
        }
    }
}

void convertFluidProperty(const CFluidProperties &mfp, MaterialLib::Fluid &fluid)
{
    if (mfp.density_model==1) {
        fluid.density = new NumLib::TXFunctionConstant(mfp.rho_0);
		fluid.drho_dp = new NumLib::TXFunctionConstant(mfp.drho_dp);
		fluid.molar_mass = new NumLib::TXFunctionConstant(mfp.molar_mass);
    }

    if (mfp.viscosity_model==1) {
        fluid.dynamic_viscosity = new NumLib::TXFunctionConstant(mfp.my_0);
    }

    if (mfp.heat_capacity_model==1) {
        fluid.specific_heat = new NumLib::TXFunctionConstant(mfp.specific_heat_capacity);
    }

    if (mfp.heat_conductivity_model==1) {
        fluid.thermal_conductivity = new NumLib::TXFunctionConstant(mfp.heat_conductivity);
    }

}

void convertSolidProperty(const CSolidProperties &msp, MaterialLib::Solid &solid)
{
    if (msp.Density_mode==1) {
        solid.density = new NumLib::TXFunctionConstant((*msp.data_Density)(0));
    }

    solid.poisson_ratio = new NumLib::TXFunctionConstant(msp.PoissonRatio);

    if (msp.Youngs_mode==1) {
        solid.Youngs_modulus = new NumLib::TXFunctionConstant((*msp.data_Youngs)(0));
    }

    if (msp.Capacity_mode==1) {
        solid.specific_heat = new NumLib::TXFunctionConstant((*msp.data_Capacity)(0));
    }

    if (msp.Conductivity_mode==1) {
        solid.thermal_conductivity = new NumLib::TXFunctionConstant((*msp.data_Conductivity)(0));
    }
}

void convertPorousMediumProperty(const Ogs5FemData &ogs5fem, const CMediumProperties &mmp, MaterialLib::PorousMedia &pm)
{
    if (mmp.permeability_model==1) {
        pm.hydraulic_conductivity = new NumLib::TXFunctionConstant(mmp.permeability_tensor[0]);
        pm.permeability = new NumLib::TXFunctionConstant(mmp.permeability_tensor[0]);
    }

    if (mmp.porosity_model==1) {
        pm.porosity = new NumLib::TXFunctionConstant(mmp.porosity_model_values[0]);
    }

    if (mmp.storage_model==1) {
        pm.storage = new NumLib::TXFunctionConstant(mmp.storage_model_values[0]);
    }

    if (mmp.mass_dispersion_model==1) {
        pm.dispersivity_long  = new NumLib::TXFunctionConstant(mmp.mass_dispersion_longitudinal); 
        pm.dispersivity_trans = new NumLib::TXFunctionConstant(mmp.mass_dispersion_transverse);
    }

	//Run through all liquid phases
	for (int i=0; i<mmp.num_phases; i++) {

		if (mmp.permeability_saturation_model[i] == 0) {  // this is a curve
            size_t perm_sat_model; 
            perm_sat_model = (size_t)mmp.permeability_saturation_model[i]; 
            pm.perm_saturation_model.push_back(perm_sat_model); 
            // get the curve number
            size_t curve_idx; 
            std::vector<double> vec_x, vec_y;
            curve_idx = (size_t)mmp.perm_saturation_value[i]; 
            // read the curve 
            ogs5fem.kurven_vector[curve_idx-1]->exportCurveValues(vec_x, vec_y);
            // create ITXFunction
            NumLib::FunctionLinear1D* f_perm_sat = NULL;
            MathLib::LinearInterpolation* xy_curve = new                      
            MathLib::LinearInterpolation(vec_x, vec_y);
            f_perm_sat = new NumLib::FunctionLinear1D(xy_curve); 
            // add it into the list
            pm.perm_saturation_curve.push_back( f_perm_sat ); 
        }

	}

	if (mmp.capillary_pressure_model == 0) {  
        // capillary pressure vs saturation curve
        pm.capp_sat_model = (size_t)mmp.capillary_pressure_model; 
        // get the curve number
        size_t curve_idx; 
        std::vector<double> vec_x, vec_y; 
        curve_idx = (size_t)mmp.capillary_pressure_values[0]; 
        // read the curve
        ogs5fem.kurven_vector[curve_idx-1]->exportCurveValues(vec_x, vec_y);
        // create ITXFunction
        MathLib::LinearInterpolation* xy_curve = new 
        MathLib::LinearInterpolation(vec_x, vec_y);
        pm.capp_sat_curve = new NumLib::FunctionLinear1D(xy_curve); 
    }
	if (mmp.capillary_pressure_model == 4) {
		// assume to be van Genuchten Model
		// the index corresponds to the OGS5 
		pm.capp_sat_model = (size_t)mmp.capillary_pressure_model;
		// get the curve number
		pm.Pb_van_Genuchten = mmp.capillary_pressure_values[0];
		pm.S_lr = mmp.capillary_pressure_values[1];
		pm.S_gr = mmp.capillary_pressure_values[2];
		pm.n_van_Genuchten = mmp.capillary_pressure_values[3];
	}

	if (mmp.capillary_pressure_model == 6){
		pm.capp_sat_model = (size_t)mmp.capillary_pressure_model;
		//assume to be Brooks corey Model
		pm.P_entry_Brook = mmp.capillary_pressure_values[0];
		pm.S_lr = mmp.capillary_pressure_values[1];
		pm.S_gr = mmp.capillary_pressure_values[2];
		pm.Lambda_Brook = mmp.capillary_pressure_values[3];
	}

    pm.geo_area = new NumLib::TXFunctionConstant(mmp.geo_area);
}


void convertFractureProperty(const CMediumProperties &mmp, MaterialLib::Fracture &pm)
{
    if (mmp.permeability_model==1) {
        pm.hydraulic_conductivity = new NumLib::TXFunctionConstant(mmp.permeability_tensor[0]);
        pm.permeability = new NumLib::TXFunctionConstant(mmp.permeability_tensor[0]);
    }

    if (mmp.porosity_model==1) {
        pm.porosity = new NumLib::TXFunctionConstant(mmp.porosity_model_values[0]);
    }

    if (mmp.storage_model==1) {
        pm.storage = new NumLib::TXFunctionConstant(mmp.storage_model_values[0]);
    }

    pm.geo_area = new NumLib::TXFunctionConstant(mmp.geo_area);
}
void convertCompoundProperty(const CompProperties &mcp, MaterialLib::Compound &new_cp)
{
    // copy the component name
    new_cp.name = mcp.compname;
    // mobility of component
    new_cp.is_mobile = (mcp.mobil == 1);
    // charge of the component
    new_cp.charge = mcp.charge;
    // checking the type of the component
    new_cp.comp_type = mcp.comp_type; 

    if (mcp.diffusion_model==1) {
        new_cp.molecular_diffusion = new NumLib::TXFunctionConstant(mcp.diffusion_model_values[0]);
    }
}


std::string convertLinearSolverType(int ls_method)
{
    std::string ls_sol_type = "";
    switch (ls_method) {
    case 1:
        ls_sol_type = "GAUSS";
        break;
    case 2:
        ls_sol_type = "BICGSTAB";
        break;
    case 3:
        ls_sol_type = "BICG";
        break;
    case 5:
        ls_sol_type = "CG";
        break;
    case 805:
        ls_sol_type = "PARDISO";
        break;
	case 9:
		ls_sol_type = "GMRES";
		break;
    default:
        break;
    } 
    return ls_sol_type;
}

std::string convertLinearSolverPreconType(int ls_precon)
{
    std::string str = "";
    switch (ls_precon) {
    case 0:
        str = "NONE";
        break;
    case 1:
        str = "JACOBI";
        break;
    case 100:
        str = "ILU";
        break;
    default:
        break;
    } 
    return str;
}

NumLib::TXFunctionType::type convertDistributionType(FiniteElement::DistributionType ogs5_type)
{
    switch (ogs5_type) {
        case FiniteElement::CONSTANT:
            return NumLib::TXFunctionType::CONSTANT;
        case FiniteElement::LINEAR:
            return NumLib::TXFunctionType::GEOSPACE;
        case FiniteElement::FUNCTION:
            return NumLib::TXFunctionType::ANALYTICAL;
        default:
            return NumLib::TXFunctionType::INVALID;
    }

}

bool convert(const Ogs5FemData &ogs5fem, Ogs6FemData &ogs6fem, BaseLib::Options &option)
{

    // -------------------------------------------------------------------------
    // Materials
    // -------------------------------------------------------------------------
    // MFP
    const size_t n_fluid = ogs5fem.mfp_vector.size();
    for (size_t i=0; i<n_fluid; i++)
    {
        CFluidProperties* mfp = ogs5fem.mfp_vector[i];

        MaterialLib::Fluid* fluid = new MaterialLib::Fluid();
        ogs6fem.list_fluid.push_back(fluid);

        convertFluidProperty(*mfp, *fluid);
    }

    // MSP
    for (size_t i=0; i<ogs5fem.msp_vector.size(); i++)
    {
        CSolidProperties* msp = ogs5fem.msp_vector[i];

        MaterialLib::Solid* solid = new MaterialLib::Solid();
        ogs6fem.list_solid.push_back(solid);

        convertSolidProperty(*msp, *solid);
    }

    // MMP
    for (size_t i=0; i<ogs5fem.mmp_vector.size(); i++)
    {
        CMediumProperties* mmp = ogs5fem.mmp_vector[i];
        MaterialLib::IMedium* mmp_ogs6 = NULL;
        if (mmp->is_fracture) {
            MaterialLib::Fracture* frac = new MaterialLib::Fracture();
            convertFractureProperty(*mmp, *frac);
            ogs6fem.list_pm.push_back(NULL);
            mmp_ogs6 = frac;
        } else {
            MaterialLib::PorousMedia* pm = new MaterialLib::PorousMedia();
            convertPorousMediumProperty(ogs5fem, *mmp, *pm);
            ogs6fem.list_pm.push_back(pm);
            mmp_ogs6 = pm;
        }
        ogs6fem.list_medium.push_back(mmp_ogs6);
    }

    // MCP
    for (size_t i=0; i<ogs5fem.cp_vector.size(); i++)
    {
        CompProperties* mcp = ogs5fem.cp_vector[i];
        MaterialLib::Compound* new_cp = new MaterialLib::Compound();
        ogs6fem.list_compound.push_back(new_cp);
        convertCompoundProperty(*mcp, *new_cp);
    }

    // -------------------------------------------------------------------------
    // Geometry
    // -------------------------------------------------------------------------
    ogs6fem.geo_unique_name = ogs5fem.geo_unique_name;
    ogs6fem.geo = ogs5fem.geo_obj;

    // -------------------------------------------------------------------------
    // Mesh
    // -------------------------------------------------------------------------
    if (ogs5fem.list_mesh.size()==0) {
        ERR("***Error: no mesh found in ogs5");
        return false;
    }
    for (size_t i=0; i<ogs5fem.list_mesh.size(); i++) {
        ogs6fem.list_mesh.push_back(ogs5fem.list_mesh[i]);
//        ogs6fem.list_dis_sys.push_back(new DiscreteLib::DiscreteSystem(*ogs5fem.list_mesh[i]));
    }

    // -------------------------------------------------------------------------
    // Time group
    // -------------------------------------------------------------------------
    convertTimeStepping(ogs5fem.time_vector, ogs6fem.list_tim);

    // -------------------------------------------------------------------------
    // User-defined curves
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Individual process and IVBV
    // -------------------------------------------------------------------------
    // PCS
    BaseLib::Options* optPcsData = option.getSubGroup("processList");
    if (optPcsData==NULL)
        optPcsData = option.addSubGroup("processList");
    size_t masstransport_counter = 0;
    if (ogs5fem.pcs_vector.size()==0) {
        ERR("***Error: no PCS found in ogs5");
        return false;
    }
    for (size_t i=0; i<ogs5fem.pcs_vector.size(); i++)
    {
        CRFProcess* rfpcs = ogs5fem.pcs_vector[i];
        std::string pcs_name = FiniteElement::convertProcessTypeToString(rfpcs->getProcessType());
        std::vector<std::string>& var_name = rfpcs->primary_variable_name;

        BaseLib::Options* optPcs = optPcsData->addSubGroup("process");
        optPcs->addOption("type", pcs_name);
        optPcs->addOption("name", pcs_name);

        //Mesh
        optPcs->addOptionAsNum("MeshID", rfpcs->mesh_id);

        //Time
        if (rfpcs->timegroup_id<0)
            rfpcs->timegroup_id = i;
        optPcs->addOptionAsNum("TimeGroupID", rfpcs->timegroup_id);

        // IC
        BaseLib::Options* optIcList = optPcs->addSubGroup("ICList");
        for (size_t i=0; i<ogs5fem.ic_vector.size(); i++)
        {
            CInitialCondition* rfic = ogs5fem.ic_vector[i];
            std::string ic_pcs_name = FiniteElement::convertProcessTypeToString(rfic->getProcessType());
			if ( ic_pcs_name.compare(pcs_name)==0 ) {
				for (size_t j=0; j<var_name.size(); j++) {
					if ( rfic->primaryvariable_name.find(var_name[j])!=std::string::npos) {
						BaseLib::Options* optIc = optIcList->addSubGroup("IC");
						optIc->addOption("Variable", rfic->primaryvariable_name);
						optIc->addOption("GeometryType", rfic->geo_type_name);
						optIc->addOption("GeometryName", rfic->geo_name);
						optIc->addOption("DistributionType", FiniteElement::convertDisTypeToString(rfic->getProcessDistributionType()));
						optIc->addOptionAsNum("DistributionValue", rfic->geo_node_value);
					}  // end of if rfic
				}
			}
        }  // end of for i

        // BC
        BaseLib::Options* optBcList = optPcs->addSubGroup("BCList");
        for (size_t i=0; i<ogs5fem.bc_vector.size(); i++)
        {
            CBoundaryCondition* rfbc = ogs5fem.bc_vector[i];
            std::string bc_pcs_name = FiniteElement::convertProcessTypeToString(rfbc->getProcessType());
            if ( bc_pcs_name.compare(pcs_name)==0 )
				for ( size_t j=0; j<var_name.size(); j++ ) {
					if ( rfbc->primaryvariable_name.find(var_name[j])!=std::string::npos) {
						BaseLib::Options* optBc = optBcList->addSubGroup("BC");
						optBc->addOption("Variable", rfbc->primaryvariable_name);
						optBc->addOption("GeometryType", rfbc->geo_type_name);
						optBc->addOption("GeometryName", rfbc->geo_name);
						NumLib::TXFunctionType::type ogs6dis_type = convertDistributionType(rfbc->getProcessDistributionType());
                        optBc->addOption("DistributionType", NumLib::convertTXFunctionTypeToString(ogs6dis_type));
						switch (rfbc->getProcessDistributionType())
						{
                            case FiniteElement::CONSTANT:
                                optBc->addOptionAsNum("DistributionValue", rfbc->geo_node_value);
                                break;
                            case FiniteElement::LINEAR:
                                {
                                    const size_t n_pt = rfbc->_PointsHaveDistribedBC.size();
                                    optBc->addOptionAsNum("PointSize", n_pt);
                                    for (size_t k=0; k<n_pt; k++) {
                                        BaseLib::Options* optPtList = optBc->addSubGroup("PointValueList");
                                        optPtList->addOptionAsNum("PointID", rfbc->_PointsHaveDistribedBC[k]);
                                        optPtList->addOptionAsNum("Value", rfbc->_DistribedBC[k]);
                                    }
                                }
                                break;
                            case FiniteElement::FUNCTION:
                                {
                                    optBc->addOption("DistributionFunction", rfbc->function_exp);
                                }
                                break;
                            default:
                                //error
                                break;
						}
					}  // end of if rfbc
				}
        }  // end of for i

        //ST
        BaseLib::Options* optStList = optPcs->addSubGroup("STList");
        for (size_t i=0; i<ogs5fem.st_vector.size(); i++)
        {
            CSourceTerm* rfst = ogs5fem.st_vector[i];
            std::string st_pcs_name = FiniteElement::convertProcessTypeToString(rfst->getProcessType());
            if (st_pcs_name.compare(pcs_name)==0 )
				for (size_t j=0; j<var_name.size(); j++ ) {
					if ( rfst->primaryvariable_name.find(var_name[j])!=std::string::npos) {
						BaseLib::Options* optSt = optStList->addSubGroup("ST");
						optSt->addOption("Variable", rfst->primaryvariable_name);
						optSt->addOption("GeometryType", rfst->geo_type_name);
						optSt->addOption("GeometryName", rfst->geo_name);
                        std::string ogs5distype = FiniteElement::convertDisTypeToString(rfst->getProcessDistributionType());
                        std::string ogs6sp_distype = ogs5distype;
                        bool isNeumannBC = (ogs5distype.find("_NEUMANN")!=std::string::npos);
                        if (isNeumannBC) {
                            optSt->addOption("STType", "NEUMANN");
                            ogs6sp_distype = BaseLib::replaceString("_NEUMANN", "", ogs6sp_distype);
                        } else {
                            optSt->addOption("STType", "SOURCESINK");
                        }
                        bool isTransientST = (rfst->CurveIndex>0);
                        std::string ogs6distype;
                        if (isTransientST) {
                            // spatial
                            BaseLib::Options* opStSpace = optSt->addSubGroup("Spatial");
                            opStSpace->addOption("DistributionType", ogs6sp_distype);
                            opStSpace->addOptionAsNum("DistributionValue", rfst->geo_node_value);
                            // temporal
                            BaseLib::Options* opStTime = optSt->addSubGroup("Temporal");
                            std::string ogs6tim_distype = NumLib::convertTXFunctionTypeToString(NumLib::TXFunctionType::T_LINEAR);
                            opStTime->addOption("DistributionType", ogs6tim_distype);
                            BaseLib::Options* opStTimeValue = opStTime->addSubGroup("TimeValue");
                            ogs5::Kurven* curv = ogs5fem.kurven_vector[rfst->CurveIndex-1];
                            for (size_t i_pt=0; i_pt<curv->stuetzstellen.size(); i_pt++) {
                                BaseLib::Options* opStTimeValuePoint = opStTimeValue->addSubGroup("Point");
                                opStTimeValuePoint->addOptionAsNum("Time", curv->stuetzstellen[i_pt]->punkt);
                                opStTimeValuePoint->addOptionAsNum("Value", curv->stuetzstellen[i_pt]->wert);
                            }

                            ogs6distype = NumLib::convertTXFunctionTypeToString(NumLib::TXFunctionType::TX);
                        } else {
                            ogs6distype = ogs6sp_distype;
                            optSt->addOptionAsNum("DistributionValue", rfst->geo_node_value);
                        }
						optSt->addOption("DistributionType", ogs6distype);
					}  // end of if rfst
				}
		}  // end of for i

        // NUM
        BaseLib::Options* optNum = optPcs->addSubGroup("Numerics");
        for (size_t i=0; i<ogs5fem.num_vector.size(); i++)
        {
            CNumerics* rfnum = ogs5fem.num_vector[i];
            if (rfnum->pcs_type_name.compare(pcs_name)==0) {
                // linear solver
                BaseLib::Options* optLS = optNum->addSubGroup("LinearSolver");
                optLS->addOption("solver_type", convertLinearSolverType(rfnum->ls_method));
                optLS->addOption("precon_type", convertLinearSolverPreconType(rfnum->ls_precond));
                //optLS->addOption("error_type", "NONE");
                optLS->addOptionAsNum("error_tolerance", rfnum->ls_error_tolerance);
                optLS->addOptionAsNum("max_iteration_step", rfnum->ls_max_iterations);

                // nonlinear solver
                BaseLib::Options* optNS = optNum->addSubGroup("NonlinearSolver");
                optNS->addOption("solver_type", rfnum->nls_method_name);
                optNS->addOptionAsNum("error_type", rfnum->nls_error_method);
                optNS->addOptionAsNum("error_tolerance", rfnum->nls_error_tolerance[0]);
                optNS->addOptionAsNum("max_iteration_step", rfnum->nls_max_iterations);

                // other stuff
                optNum->addOptionAsNum("TimeTheta", rfnum->ls_theta);
                optNum->addOptionAsNum("GaussPoint", rfnum->ele_gauss_points);
                optNum->addOptionAsNum("FEM_FCT", rfnum->fct_method);
            }
        }

        // some process specific
        if (pcs_name.find("MASS_TRANSPORT")!=std::string::npos) {
            optPcs->addOptionAsNum("CompoundID", masstransport_counter);
            masstransport_counter++;
        }
		
    }  // end of for ogs5fem.pcs_vector.size()

	// HS: add chemical reactions
	size_t n_KinReactions, n_Compound; 
	n_KinReactions = ogs5fem.KinReact_vector.size(); 
	if ( n_KinReactions > 0 ) {  // only create the "optKRC" option when there are kinetic reactions defined. 
		// HS: convert all compounds to ChemComp data structure
		n_Compound = ogs6fem.list_compound.size(); 
		for ( size_t i=0; i < n_Compound ; i++ )
		{
			ogsChem::ChemComp* mChemComp = NULL; 
			// directly use the type of constructor to initialize the ChemComp data structure
			mChemComp = new ogsChem::ChemComp(); 
			// set index
			mChemComp->setIndex(i); 
			// set name
			mChemComp->set_name(ogs6fem.list_compound[i]->name); 
			// set mobility
			if ( ogs6fem.list_compound[i]->is_mobile )
				mChemComp->set_mobility( ogsChem::MOBILE  ); 
			else
			{
				mChemComp->set_mobility(ogsChem::IMMOBILE); 
			}

            // set component type
            if ( ogs6fem.list_compound[i]->comp_type.find("BASIS_COMP") == 0 )
                mChemComp->set_compTyps(ogsChem::BASIS_COMP);
            else if ( ogs6fem.list_compound[i]->comp_type.find("AQ_PHASE_COMP") == 0 )
                mChemComp->set_compTyps(ogsChem::AQ_PHASE_COMP);
            else if ( ogs6fem.list_compound[i]->comp_type.find("GAS_PHASE_COMP") == 0)
                mChemComp->set_compTyps(ogsChem::GAS_PHASE_COMP);
            else if ( ogs6fem.list_compound[i]->comp_type.find("SORPTION_COMP") == 0)
                mChemComp->set_compTyps(ogsChem::SORPTION_COMP); 
            else if ( ogs6fem.list_compound[i]->comp_type.find("MIN_PHASE_COMP") == 0 )
                mChemComp->set_compTyps(ogsChem::MIN_PHASE_COMP);
            else if ( ogs6fem.list_compound[i]->comp_type.find("SS_PHASE_COMP") == 0 )
                mChemComp->set_compTyps(ogsChem::SS_PHASE_COMP);
			else if ( ogs6fem.list_compound[i]->comp_type.find("KIN_COMP") == 0 )
				mChemComp->set_compTyps(ogsChem::KIN_COMP);

            // set charge
            mChemComp->set_charge(ogs6fem.list_compound[i]->charge); 

			ogs6fem.map_ChemComp.insert(ogs6fem.list_compound[i]->name, mChemComp); 
			mChemComp = NULL; 
		}

		// adding the reactions one after another
		for ( size_t i=0; i < n_KinReactions ; i++ )
		{
			// add reactions one after another
			ogs5::CKinReact* rfKinReact = ogs5fem.KinReact_vector[i];
            // check which type of reaction it is. 
            if ( rfKinReact->getType().find( "EQ_REACT" ) != std::string::npos )
            {   // if it is an equilibrium reaction
                // eq_mob?
                if ( rfKinReact->getType().find( "MOB_EQ_REACT" ) == 0 )
                {
                    ogsChem::chemReactionEqMob* mEqReaction; 
                    mEqReaction = new ogsChem::chemReactionEqMob(); 
                    // reading the reaction stoichiometry
                    mEqReaction->readReactionKRC( ogs6fem.map_ChemComp, rfKinReact );
                    // adding the instance into the list
                    ogs6fem.list_eq_reactions.push_back(mEqReaction); 
                }
                // eq_sorp?
                else if ( rfKinReact->getType().find( "SORP_EQ_REACT" ) == 0 )
                {
                    ogsChem::chemReactionEqSorp* mEqReaction; 
                    mEqReaction = new ogsChem::chemReactionEqSorp(); 
                    // reading the reaction stoichiometry
                    mEqReaction->readReactionKRC( ogs6fem.map_ChemComp, rfKinReact );
                    // adding the instance into the list
                    ogs6fem.list_eq_reactions.push_back(mEqReaction); 
                }
                // eq_min?
                else if ( rfKinReact->getType().find( "MIN_EQ_REACT" ) == 0 )
                {
                    ogsChem::chemReactionEqMin* mEqReaction; 
                    mEqReaction = new ogsChem::chemReactionEqMin(); 
                    // reading the reaction stoichiometry
                    mEqReaction->readReactionKRC( ogs6fem.map_ChemComp, rfKinReact );
                    // adding the instance into the list
                    ogs6fem.list_eq_reactions.push_back(mEqReaction); 
                }
            }
            else 
            {   // if it is kinetic reaction
                // creating reaction instance 
                ogsChem::chemReactionKin* mKinReaction; 
                mKinReaction = new ogsChem::chemReactionKin(); 
                // convert the ogs5 KRC data structure into ogs6 Kinetic reactions
                mKinReaction->readReactionKRC( ogs6fem.map_ChemComp, rfKinReact ); 
                // adding the instance of one single kinetic reaction
                ogs6fem.list_kin_reactions.push_back(mKinReaction); 
            }  // end of if else
        }  // end of for


        // using the process name to control which class to initialize
		bool HAVE_KIN_REACT_GIA(false), HAVE_REACT_GIA(false), HAVE_REACT_TRANS_OPS(false);
		for (size_t i = 0; i < ogs5fem.pcs_vector.size(); i++)
		{
			ogs5::FiniteElement::ProcessType pcs_type = ogs5fem.pcs_vector[i]->getProcessType(); 

			if ( pcs_type == ogs5::FiniteElement::ProcessType::KIN_REACT_GIA )
				HAVE_KIN_REACT_GIA = true;
			else if (pcs_type == ogs5::FiniteElement::ProcessType::REACT_GIA)
				HAVE_REACT_GIA = true;
			else if (pcs_type == ogs5::FiniteElement::ProcessType::REACT_TRANS_OPS)
				HAVE_REACT_TRANS_OPS = true; 
		}

        // first case the operator splitting mode
		if (HAVE_REACT_TRANS_OPS)
        {
            // initialize the activity model of the chemical system. 
            ogs5::CKinReactData* ogs5KinReactData = ogs5fem.KinReactData_vector[0]; 
            ogsChem::chemEqReactSysActivity* mEqReactSys; 
            ogsChem::chemKinReactSys* mKinReactSys; 
            if ( ogs5KinReactData->activity_model == 1 ) // Debyl-Huekel
            {
                ogsChem::chemActivityModelAqDBH* mActivityModelAqDBH 
                            = new ogsChem::chemActivityModelAqDBH(ogs6fem.map_ChemComp); 
                mEqReactSys = new ogsChem::chemEqReactSysActivity(ogs6fem.map_ChemComp, 
                                                                  ogs6fem.list_eq_reactions, 
                                                                  mActivityModelAqDBH); 
                mKinReactSys = new ogsChem::chemKinReactSys(ogs6fem.map_ChemComp,
                                                            ogs6fem.list_kin_reactions,
                                                            mActivityModelAqDBH);
                ogs6fem.m_EqReactSys = mEqReactSys;
                ogs6fem.m_KinReactSys = mKinReactSys;
                
            }
            else if ( ogs5KinReactData->activity_model == 2 ) // Davies
            {
                ogsChem::chemActivityModelAqDVS* mActivityModelAqDVS 
                            = new ogsChem::chemActivityModelAqDVS(ogs6fem.map_ChemComp); 
                mEqReactSys = new ogsChem::chemEqReactSysActivity(ogs6fem.map_ChemComp, 
                                                                  ogs6fem.list_eq_reactions, 
                                                                  mActivityModelAqDVS); 
                mKinReactSys = new ogsChem::chemKinReactSys(ogs6fem.map_ChemComp,
                                                            ogs6fem.list_kin_reactions,
                                                            mActivityModelAqDVS);
                ogs6fem.m_EqReactSys = mEqReactSys; 
                ogs6fem.m_KinReactSys = mKinReactSys;
                
            }
            else // by default, use the unity model
            {
                ogsChem::chemActivityModelUnity* mActivityModelUnity 
                            = new ogsChem::chemActivityModelUnity();
                mEqReactSys = new ogsChem::chemEqReactSysActivity(ogs6fem.map_ChemComp, 
                                                                  ogs6fem.list_eq_reactions, 
                                                                  mActivityModelUnity); 
                mKinReactSys = new ogsChem::chemKinReactSys(ogs6fem.map_ChemComp,
                                                            ogs6fem.list_kin_reactions,
                                                            mActivityModelUnity); 
                ogs6fem.m_EqReactSys = mEqReactSys; 
                ogs6fem.m_KinReactSys = mKinReactSys; 
                
            }
            mEqReactSys = NULL; 
            mKinReactSys = NULL; 
        }
        // then is the the KIN_GIA mode
		else if ( HAVE_KIN_REACT_GIA )
        {   // initialize the kin-reduction scheme 
            ogsChem::chemReductionKin* mReductionKinScheme;
            mReductionKinScheme = new ogsChem::chemReductionKin( ogs6fem.map_ChemComp, ogs6fem.list_kin_reactions );
            ogs6fem.m_KinReductScheme = mReductionKinScheme;
            mReductionKinScheme = NULL;
        }
		//  RZ:2.12.2013 then the REDUCT_GIA mode
		else if ( HAVE_REACT_GIA )
        {  // initialize the full reduction scheme
        	ogs5::CKinReactData* ogs5KinReactData = ogs5fem.KinReactData_vector[0];
            ogsChem::chemReductionGIA* mReductionGIAScheme;
            if ( ogs5KinReactData->activity_model == 1 ) // Debyl-Huekel
            {
                ogsChem::chemActivityModelAqDBH* mActivityModelAqDBH = new ogsChem::chemActivityModelAqDBH(ogs6fem.map_ChemComp);

                mReductionGIAScheme = new ogsChem::chemReductionGIA(ogs6fem.map_ChemComp, ogs6fem.list_eq_reactions, ogs6fem.list_kin_reactions, mActivityModelAqDBH);
                ogs6fem.m_GIA_ReductScheme = mReductionGIAScheme;

            }
            else if ( ogs5KinReactData->activity_model == 2 ) // Davies
            {
                ogsChem::chemActivityModelAqDVS* mActivityModelAqDVS = new ogsChem::chemActivityModelAqDVS(ogs6fem.map_ChemComp);
                mReductionGIAScheme = new ogsChem::chemReductionGIA(ogs6fem.map_ChemComp, ogs6fem.list_eq_reactions, ogs6fem.list_kin_reactions, mActivityModelAqDVS);
                ogs6fem.m_GIA_ReductScheme = mReductionGIAScheme;
            }
            else // by default, use the unity model
            {
                ogsChem::chemActivityModelUnity* mActivityModelUnity = new ogsChem::chemActivityModelUnity();
                mReductionGIAScheme = new ogsChem::chemReductionGIA(ogs6fem.map_ChemComp, ogs6fem.list_eq_reactions, ogs6fem.list_kin_reactions, mActivityModelUnity);
                ogs6fem.m_GIA_ReductScheme = mReductionGIAScheme;

            }
            mReductionGIAScheme = NULL;
        }  // end of if else
    }  // end of if ( n_KinReactions > 0 )

    // -------------------------------------------------------------------------
    // Coupling
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Output
    // -------------------------------------------------------------------------

    // OUT
    BaseLib::Options* optOut = option.addSubGroup("outputList");
    for (size_t i=0; i<ogs5fem.out_vector.size(); i++)
    {
        COutput* rfout = ogs5fem.out_vector[i];

        // convert *_X1 to *_X (e.g. DISPLACEMENT_X1 to DISPLACEMENT_X)
        for (size_t j=0; j<rfout->_nod_value_vector.size(); j++) {
            const size_t len_postfix = 3;
            std::string str = rfout->_nod_value_vector[j];
            if (str.length()<len_postfix+1) continue;
            std::string last3 = str.substr(str.length()-len_postfix, len_postfix);
            if (last3.compare("_X1")==0) {
                rfout->_nod_value_vector[j].replace(str.length()-len_postfix, len_postfix, "_X");
            } else if (last3.compare("_Y1")==0) {
                rfout->_nod_value_vector[j].replace(str.length()-len_postfix, len_postfix, "_Y");
            } else if (last3.compare("_Z1")==0) {
                rfout->_nod_value_vector[j].replace(str.length()-len_postfix, len_postfix, "_Z");
            }
        }

        BaseLib::Options* opt = optOut->addSubGroup("output");
        opt->addOption("dataType", rfout->dat_type_name);
        opt->addOption("meshID", "0"); //TODO
        opt->addOption("geoType", rfout->geo_type);
        opt->addOption("geoName", rfout->geo_name);
        opt->addOption("timeType", rfout->tim_type_name);
        opt->addOptionAsNum("timeSteps", rfout->nSteps);
        for (size_t j=0; j<rfout->_nod_value_vector.size(); j++) {
            BaseLib::Options* optVal = opt->addSubGroup("nodeValue");
            optVal->addOption("name", rfout->_nod_value_vector[j]);
        }
        for (size_t j=0; j<rfout->_ele_value_vector.size(); j++) {
            BaseLib::Options* optVal = opt->addSubGroup("elementValue");
            optVal->addOption("name", rfout->_ele_value_vector[j]);
        }
//        opt->addOptionAsArray("MMPValues", rfout->mmp_value_vector);
//        opt->addOptionAsArray("MFPValues", rfout->mfp_value_vector);
    }

    ogs6fem.outController.initialize(option, ogs6fem.output_dir, ogs6fem.project_name, ogs6fem.list_mesh, *ogs6fem.geo, ogs6fem.geo_unique_name);

    return true;
}

};
} //end
