/**
    * \brief Local assembly of time ODE components for RichardsFlow in porous media
    */

template <class T>
class RichardsFlowTimeODELocalAssembler: public T
{
public:
    typedef MathLib::LocalVector LocalVectorType;
    typedef MathLib::LocalMatrix LocalMatrixType;

    RichardsFlowTimeODELocalAssembler(FemLib::LagrangeFeObjectContainer* feObjects, const MeshLib::CoordinateSystem &problem_coordinates)
    : _feObjects(*feObjects), _problem_coordinates(problem_coordinates)
    {};

    virtual ~RichardsFlowTimeODELocalAssembler() {};

protected:
    virtual void assembleODE(const NumLib::TimeStep &/*time*/, const MeshLib::IElement &e, const LocalVectorType &u1,   const LocalVectorType &/*u0*/, LocalMatrixType &localM, LocalMatrixType &localK, LocalVectorType &localF)
    {
		FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
        //const size_t n_dim = e.getDimension();
        size_t mat_id = e.getGroupID();
        // get pointer to corresponding porous media class
        MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];
        // get pointer to corresponding fluid (water) class
        MaterialLib::Fluid* fluid = Ogs6FemData::getInstance()->list_fluid[0];
        // whether include gravity
        // const bool hasGravityEffect = _problem_coordinates.hasZ();
        const bool hasGravityEffect = true;

        MathLib::LocalMatrix mass_mat_coeff    = MathLib::LocalMatrix::Zero(1,1);  // coefficient of mass matrix
        MathLib::LocalMatrix Pw                = MathLib::LocalMatrix::Zero(1,1);  // water pressure 
        MathLib::LocalMatrix local_k_mu;                                           // permeability divided by viscosity
        MathLib::LocalVector vec_g             = MathLib::LocalVector::Zero(_problem_coordinates.getDimension());
        double storage  = 0.0;  // storage term
        double Sw       = 0.0;  // saturation of water
        double rho_w    = 0.0;  // density of water
        double drhow_dp = 0.0;  // derivative of water density over pressure
        double Pc       = 0.0;  // capillary pressure 
        double poro     = 0.0;  // porosity 
        double dSwdPc   = 0.0;  // derivative of water saturation over cap pressure
        double k        = 0.0;  // intrinsic permeability
        double k_rel    = 0.0;  // relative permeability
        double mu       = 0.0;  // dynamic viscosity of water

        // if there is gravity, set up a gravity vector
        if (hasGravityEffect) {
            vec_g[_problem_coordinates.getIndexOfY()] = -9.81;
        }
        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
            q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);
            NumLib::TXPosition gp_pos(NumLib::TXPosition::IntegrationPoint, e.getID(), j, real_x);
            MathLib::LocalMatrix &Np  = *fe->getBasisFunction(); // get basis function
            MathLib::LocalMatrix &dNp = *fe->getGradBasisFunction(); // get gradient basis function
            local_k_mu = MathLib::LocalMatrix::Identity(e.getDimension(), e.getDimension());
            double fac = fe->getDetJ() * q->getWeight(j);

            // geting variables and parameters
            // get porosity
            pm->porosity->eval(gp_pos, poro);
            // get storage
            pm->storage->eval(gp_pos, storage);
            // get primary variable water pressure pw
            // u1 is the current step primary variable. 
            Pw = Np * u1; 
            Pc = -1.0 * Pw(0,0); 
            // get water saturation using pw
            Sw = pm->getSwbyPc(Pc);
            // density of water
            fluid->density->eval(gp_pos, rho_w);
            // get drhow_dp
            fluid->drho_dp->eval(gp_pos, drhow_dp); 
            // viscosity of the fluid
            fluid->dynamic_viscosity->eval(gp_pos, mu);
			// get dSwdPc
            dSwdPc = pm->getdSwdPc( Pc ); 
            // get k_rel
            k_rel = pm->getKrelbySw(Sw,0 );  // 0 stands for aq. phase
            // get intrinsic permeability
            pm->permeability->eval(gp_pos, k); 

            // calculate mass matrix coefficient
            mass_mat_coeff(0,0) = storage * Sw + poro * Sw * drhow_dp - poro * dSwdPc; 
            // multiply shape shape 
            fe->integrateWxN(j, mass_mat_coeff, localM);

            // calculate laplace matrix coefficient
            local_k_mu *= k * k_rel / mu; 
            // multiply dshape dshape
            fe->integrateDWxDN(j, local_k_mu, localK);
            // if includes gravity
            if (hasGravityEffect) { 
                // since no primary vairable involved
                // directly assemble to the Right-Hand-Side
                // F += dNp^T * K * gz
                localF.noalias() += fac * dNp.transpose() * local_k_mu * rho_w * vec_g;
            } // end of if hasGravityEffect
        } // end of for GP

        // testing mass lumping----------------------------
        // std::cout << "Local Mass Matrix before lumping: " << std::endl; 
        // std::cout << localM << std::endl; 
        for (int idx_ml=0; idx_ml < localM.cols(); idx_ml++ )
        {
            double mass_lump_val;
            mass_lump_val = localM.col(idx_ml).sum();
            localM.col(idx_ml).setZero(); 
            localM(idx_ml, idx_ml) = mass_lump_val; 
        }
        // std::cout << "Local Mass Matrix after lumping: " << std::endl; 
        // std::cout << localM << std::endl; 
        // std::cout << "Local K Matrix: " << std::endl; 
        // std::cout << localK << std::endl; 
        //end of testing mass lumping---------------------

    }  // end of assemble ODE


private:
    FemLib::LagrangeFeObjectContainer _feObjects;
	MeshLib::CoordinateSystem _problem_coordinates;
};
