/**
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ReductConc.h
 *
 * Created on    2013-05-17 by Reza Zolfaghari & Haibing Shao
 */

#ifndef REDUCT_CONC_H
#define REDUCT_CONC_H

#include "BaseLib/CodingTools.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerEQSLocalAssemblerWithStorage.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerResidualLocalAssembler.h"
#include "NumLib/Function/DiscreteDataConvergenceCheck.h"
#include "SolutionLib/Fem/FemIVBVProblem.h"
#include "SolutionLib/Fem/TemplateFemEquation.h"
#include "SolutionLib/Fem/SingleStepFEM.h"
#include "ProcessLib/AbstractTransientProcess.h"
#include "MathLib/DataType.h"
#include "NumLib/Nonlinear/DiscreteNRSolverWithStepInitFactory.h"
#include "FemGIARedSolution.h"
#include "SingleStepGIAReduction.h"
#include "NestedGIALocalProbNRIterationStepInitializer.h"
#include "LinearTransportTimeODELocalAssemberML.h"
#include "UserModules/FemKinReactGIA/LinearTransportJacobianLocalAssembler.h"
#include "NonLinearGIATimeODELocalAssembler.h"
#include "NonLinearGIAJacabianLocalAssembler.h"
#include "LocalProblem.h"
#include  "ogs6THMC/ChemLib/chemReductionGIA.h"
#include "TemplateFemEquation_GIAReduct.h"
#include "SingleStepGIAReductionNonlinear.h"
#include "Local_ODE_Xi_immob_GIA.h"


template <class T_DISCRETE_SYSTEM, class T_LINEAR_SOLVER >
class FunctionReductConc
	: public ProcessLib::Process
{
public:
	// input variable is velocity
    enum In { Velocity=0 };
	// no output variable
    enum Out { Concentrations=0 };

	// local matrix and vector
	typedef MathLib::LocalMatrix LocalMatrix;
	typedef MathLib::LocalVector LocalVector;

	typedef FunctionReductConc     MyFunctionData;    // FEM problem class
    typedef T_DISCRETE_SYSTEM      MyDiscreteSystem;  // Discretization
    typedef T_LINEAR_SOLVER        MyLinearSolver;    // linear solver
	
    // memory for discretized concentration vector
    typedef typename FemLib::FemNodalFunctionVector<MyDiscreteSystem>::type MyNodalFunctionVector;
	typedef typename FemLib::FemNodalFunctionScalar<MyDiscreteSystem>::type MyNodalFunctionScalar;

    // local assembler
	// for the linear systems, use the same settings as Mass_Transport
    typedef LinearTransportTimeODELocalAssemblerML<NumLib::ElementWiseTimeEulerEQSLocalAssemblerWithStorage> MyLinearAssemblerType;
    typedef LinearTransportTimeODELocalAssemblerML<NumLib::ElementWiseTimeEulerResidualLocalAssembler> MyLinearResidualAssemblerType;
    typedef LinearTransportJacobianLocalAssembler MyLinearJacobianAssemblerType;                                        

	// to be changed
	typedef NestedGIALocalProbNRIterationStepInitializer<MyFunctionData>                          MyNRIterationStepInitializer;
	typedef NumLib::DiscreteNRSolverWithStepInitFactory<MyNRIterationStepInitializer>             MyDiscreteNonlinearSolverFactory;
	
	/**
      * linear Equation definition
      */
    typedef SolutionLib::TemplateFemEquation<
            MyDiscreteSystem,
            MyLinearSolver,
            MyLinearAssemblerType,
            MyLinearResidualAssemblerType,
            MyLinearJacobianAssemblerType
            > MyLinearEquationType;

    /**
      * linear IVBV problem definition
      */
    typedef SolutionLib::FemIVBVProblem<
            MyDiscreteSystem,
            MyLinearEquationType
            > MyLinearTransportProblemType;

    /** 
      * linear equation solution algorithm definition
      */
    typedef SolutionLib::SingleStepFEM<
            MyLinearTransportProblemType,
            MyLinearSolver
            > MyLinearSolutionType;

	/**
      * nonlinear coupled part solving xi
      * Equation definition
      */
	typedef TemplateFemEquation_GIAReduct<
		    MyDiscreteSystem,
            MyFunctionData,
			MyLinearSolver
	        > MyNonLinearEquationType;
	/**
      * FEM IVBV problem definition
      */
	typedef SolutionLib::FemIVBVProblem<
		    MyDiscreteSystem,
	        MyNonLinearEquationType
	        > MyNonLinearReactiveTransportProblemType;
	
    /**
	  * Solution algorithm definition
	  */
	typedef SingleStepGIAReductionNonlinear<
			MyNonLinearReactiveTransportProblemType,
			MyFunctionData,
			MyLinearSolver, 
			MyDiscreteNonlinearSolverFactory
			> MyNonLinearSolutionType;

    /**
	  * the general reduction problem part
	  */
	typedef SolutionLib::FemGIARedSolution<MyDiscreteSystem> MyGIAReductionProblemType; 

	typedef SolutionLib::SingleStepGIAReduction<MyFunctionData, 
		                                        MyGIAReductionProblemType,
		                                        MyLinearTransportProblemType, 
		                                        MyLinearSolutionType, 
												MyNonLinearReactiveTransportProblemType, 
	                                            MyNonLinearSolutionType> MyGIAReductionSolution;
    typedef typename MyGIAReductionProblemType::MyVariable MyVariableConc;

    FunctionReductConc() 
        : Process("REACT_GIA_REDUCT", 1, 1),
          _feObjects(0),
    	_non_linear_eqs(0),  //RZ
    	_n_xi_global(0), _current_time_step(0),	myNRIterator(0),
    	_solution(0), _n_xi_Kin_bar(0), _I_mob(0), _n_xi_Sorp_bar(0), _n_eta(0), _I_min(0), _non_linear_solution(0), _dofManager(0),
    	_J_tot_kin(0), _n_eta_bar(0), _n_Comp(0), _n_xi_local(0), _vel(0), _n_xi_Sorp_bar_li(0), _n_xi_Sorp(0), _ReductionGIA(0),
    	_I_sorp(0), _dis_sys(0), _nl_sol_dofManager(0), _n_xi_Sorp_tilde(0), myNSolverFactory(0),
    	_n_xi_Mob(0), _non_linear_problem(0), _n_xi_Min_bar(0), _n_xi_Min_tilde(0), _n_xi_Min(0), _n_xi_Kin(0),
    	_theta(0),_tim(0), _problem(0), _n_xi_Sorp_bar_ld(0),_msh_id(0)
    {
        // set default parameter name
		ProcessLib::Process::setInputParameterName(Velocity, "Velocity");
    };
	
    /**
     * destructor, reclaim the memory
     */
    virtual ~FunctionReductConc()
    {
    	BaseLib::releaseObject(_feObjects);
    	BaseLib::releaseObject(_non_linear_eqs);
    	BaseLib::releaseObject(myNRIterator);
    	BaseLib::releaseObject(myNSolverFactory);

    	BaseLib::releaseObject(_problem);
    	BaseLib::releaseObject(_solution);

    	BaseLib::releaseObject(_ReductionGIA);
    	BaseLib::releaseObject(_sbs);
    	BaseLib::releaseObject(_pSolve);

    	BaseLib::releaseObject(_non_linear_eqs);
    	BaseLib::releaseObject(_non_linear_problem);
    	BaseLib::releaseObject(_non_linear_solution);

    	BaseLib::releaseObjectsInStdVector(_linear_problems);
    	BaseLib::releaseObjectsInStdVector(_linear_solutions);

    	BaseLib::releaseObjectsInStdVector(_concentrations);
    	BaseLib::releaseObjectsInStdVector(_eta);
    	BaseLib::releaseObjectsInStdVector(_eta_bar);
    	BaseLib::releaseObjectsInStdVector(_xi_global_cur);
        BaseLib::releaseObjectsInStdVector(_xi_global_pre);
    	BaseLib::releaseObjectsInStdVector(_xi_local_new);
    	BaseLib::releaseObjectsInStdVector(_xi_local_old);

		BaseLib::releaseObjectsInStdVector(_xi_sorp_rates);
		BaseLib::releaseObjectsInStdVector(_xi_min_rates);
		BaseLib::releaseObjectsInStdVector(_xi_kin_rates);
		BaseLib::releaseObjectsInStdVector(_nodal_vec_AI);

    	BaseLib::releaseObjectsInStdVector(_global_vec_Rate);
    	BaseLib::releaseObject(_nl_sol_dofManager);
    	BaseLib::releaseObject(_feObjects);
    	BaseLib::releaseObject(_problem);
    	BaseLib::releaseObject(_solution);
    	BaseLib::releaseObject(_dis_sys);
    	BaseLib::releaseObject(_dofManager);

    };


    /**
	  * initialization of the problem class
	  */
    virtual bool initialize(const BaseLib::Options &option);

    /**
	  * finalize but nothing to do here
	  */
    virtual void finalize() {};

	/**
	  * returns the convergence checker
	  */
    virtual NumLib::IConvergenceCheck* getConvergenceChecker()
    {
        return &_checker;
    }

	/**
	  * function to solve the current time step
	  */
	virtual int solveTimeStep(const NumLib::TimeStep &time)
    {
		INFO("Solving %s...", getProcessName().c_str());
        initializeTimeStep(time);
        getSolution()->solveTimeStep(time);
        updateOutputParameter(time);
        return 0;
    }

    /**
	  * function to suggest the next time step
	  */
    virtual double suggestNext(const NumLib::TimeStep &time_current) 
    {
        return getSolution()->suggestNext(time_current); // is crashing here. why dt is zero?
    }

    /**
	  * called when this problem is awake
	  */
    virtual bool isAwake(const NumLib::TimeStep &time) 
    { 
        return getSolution()->isAwake(time);
    }

    virtual bool accept(const NumLib::TimeStep &time)
    {
        return getSolution()->accept(time);
    }

    /**
	  * called when this time step is accepted
	  */
    virtual void finalizeTimeStep(const NumLib::TimeStep & time)
    {
         output(time);
		 getSolution()->finalizeTimeStep(time);
    };

	void start_node_values_search(MathLib::LocalMatrix &mat_S1_ast,
								  MathLib::LocalVector &optimalXi,
								  MathLib::LocalVector &local_xi_Mob,
								  MathLib::LocalVector &local_xi_Sorp_tilde,
								  MathLib::LocalVector &local_xi_Sorp_bar,
								  MathLib::LocalVector &local_xi_Min_tilde,
								  MathLib::LocalVector &local_xi_Min_bar,
								  MathLib::LocalVector &local_xi_Kin,
								  MathLib::LocalVector &local_xi_Kin_bar,
								  MathLib::LocalVector &local_eta,
								  MathLib::LocalVector &local_etabar,
								  MathLib::LocalVector &local_conc);

	void calculate_node_concentration(MathLib::LocalVector &local_xi_Mob,
									  MathLib::LocalVector &local_xi_Sorp_tilde,
									  MathLib::LocalVector &local_xi_Sorp_bar,
									  MathLib::LocalVector &local_xi_Min_tilde,
									  MathLib::LocalVector &local_xi_Min_bar,
									  MathLib::LocalVector &local_xi_Kin,
									  MathLib::LocalVector &local_xi_Kin_bar,
									  MathLib::LocalVector &loc_eta,
									  MathLib::LocalVector &loc_etabar,
									  MathLib::LocalVector &loc_xi_global,
									  MathLib::LocalVector &loc_xi_local,
									  MathLib::LocalVector &loc_conc,
									  MathLib::LocalMatrix &mat_S1_ast,
									  MathLib::LocalVector &optimalXi,
									  size_t &node_idx);

	void calculate_concentration_using_optimal_xi(MathLib::LocalVector &optimalXi,
			   	   	   	   	   	   	   	   	   	  MathLib::LocalVector &local_xi_Mob,
			   	   	   	   	   	   	   	   	   	  MathLib::LocalVector &local_xi_Sorp_tilde,
			   	   	   	   	   	   	   	   	   	  MathLib::LocalVector &local_xi_Sorp_bar,
			   	   	   	   	   	   	   	   	   	  MathLib::LocalVector &local_xi_Min_tilde,
			   	   	   	   	   	   	   	   	   	  MathLib::LocalVector &local_xi_Min_bar,
			   	   	   	   	   	   	   	   	   	  MathLib::LocalVector &local_xi_Kin,
			   	   	   	   	   	   	   	   	   	  MathLib::LocalVector &local_xi_Kin_bar,
			   	   	   	   	   	   	   	   	   	  MathLib::LocalVector &local_eta,
			   	   	   	   	   	   	   	   	   	  MathLib::LocalVector &local_etabar,
			   	   	   	   	   	   	   	   	   	  MathLib::LocalVector &local_conc);

    /**
	  * set function for eta and xi
	  */
	void set_eta_node_values     	( std::size_t eta_idx,   MyNodalFunctionScalar* new_eta_node_values   );
	void set_eta_bar_node_values    ( std::size_t eta_bar_idx, MyNodalFunctionScalar* new_eta_bar_node_values );
	void set_xi_global_node_values  ( std::size_t xi_global_idx,    MyNodalFunctionScalar* new_xi_global_node_values    );

	/**
      * get function for _nodal_vec_AI
	  */
	ogsChem::LocalVector& get_nodal_vec_AI(std::size_t node_index) { return *(_nodal_vec_AI[node_index]); };
    
    template <class T_X>
    void update_xi_global_cur_nodal_values  ( const T_X & x_new );

	
    void update_xi_local_node_values ( void );

	/**
      * calculate the reaction rates on each node
      */
	void update_node_kin_reaction_rates(void);

	/**
      * update the change of kinetic reaction rate over the change of xi_mob
      */
	void update_node_kin_reaction_drates_dxi(void);

	/**
      * calculate nodal local problem
      */
	void calc_nodal_local_problem(double dt, const double iter_tol, const double rel_tol, const double max_iter);

	ogsChem::chemReductionGIA* getReductionGIA(){return _ReductionGIA;}

	std::vector<MyNodalFunctionScalar*> & get_xi_global_cur() {return _xi_global_cur;}
    std::vector<MyNodalFunctionScalar*> & get_xi_global_pre() {return _xi_global_pre;}
	std::vector<MyNodalFunctionScalar*> & get_xi_local_new() {return _xi_local_new;}
	std::vector<MyNodalFunctionScalar*> & get_eta() {return _eta;}
	std::vector<MyNodalFunctionScalar*> & get_eta_bar() {return _eta_bar;}
	std::vector<MyNodalFunctionScalar*> & get_global_vec_Rate() {return _global_vec_Rate; }
	std::vector<MyNodalFunctionScalar*> & get_concentrations() {return _concentrations; }
	std::vector<MyNodalFunctionScalar*> & get_drates_dxi() {return _drates_dxi; }
	std::vector<MyNodalFunctionScalar*> & get_xi_sorp_rates() { return _xi_sorp_rates; }
	std::vector<MyNodalFunctionScalar*> & get_xi_min_rates() { return _xi_min_rates; }
	std::vector<MyNodalFunctionScalar*> & get_xi_kin_rates() { return _xi_kin_rates; }

	FemLib::LagrangeFeObjectContainer* get_feObjects(){return _feObjects;}
	NumLib::TimeStep const* getTimeStep() {return _current_time_step;}
    /// get the time step function
    /// @return Time step function
    NumLib::ITimeStepFunction* getTimeStepFunction() const {return _tim;};  //RZ: 11.12.2013
    virtual void set_BC_conc_node_values(std::size_t node_idx, std::size_t i_var, double node_value);
    /**
      * convert nodal concentration values to eta and xi
      */
	virtual void convert_conc_to_eta_xi(void);

    MyGIAReductionSolution* getSolution(void) { return _solution; }

    void copy_cur_xi_global_to_pre(void); 

    void copy_cur_xi_local_to_pre(void);
    /**
      * convert nodal eta and xi values to concentrations
      */
	virtual void convert_eta_xi_to_conc(void);

	virtual void update_lnK(void);

	void Perform_Extrapolation();


protected:
    virtual void initializeTimeStep(const NumLib::TimeStep &time);

    /**
      * this function is called to exchange output parameters
      */ 
    virtual void updateOutputParameter(const NumLib::TimeStep &time);

    /**
      * output the result of current solution
      */ 
    virtual void output(const NumLib::TimeStep &time);


private:
    DISALLOW_COPY_AND_ASSIGN(FunctionReductConc);

private:


	DiscreteLib::DofEquationIdTable* _dofManager;

	MyDiscreteSystem *_dis_sys;

	NumLib::ITimeStepFunction* _tim;

	NumLib::TimeStep const* _current_time_step;

	/**
      * linear problems
      */
	std::vector<MyLinearTransportProblemType*> _linear_problems;

    /**
      * linear solutions
      */
    std::vector<MyLinearSolutionType*>         _linear_solutions;

    /**
      * Nonlinear iterator
      */
	MyNRIterationStepInitializer*              myNRIterator;
	
    /**
      * Nonlinear solver factory
      */
    MyDiscreteNonlinearSolverFactory*          myNSolverFactory;
    
    /**
      * nonlinear equation
      */ 
	MyNonLinearEquationType*                  _non_linear_eqs;
    
    /**
      * non-linear problem
      */ 
	MyNonLinearReactiveTransportProblemType * _non_linear_problem;

    /**
      * non-linear solution
      */
	MyNonLinearSolutionType*                  _non_linear_solution;

	/**
      * reduction problem
      */ 
	MyGIAReductionProblemType* _problem;

    /** 
      * reduction solution
      */
	MyGIAReductionSolution*    _solution;

	/**
      * FEM object
      */
    FemLib::LagrangeFeObjectContainer* _feObjects;
    
    double _theta;

    /**
      * velocity function
      */
    NumLib::ITXFunction* _vel;

	/**
      * convergence checker
      */
	NumLib::DiscreteDataConvergenceCheck _checker;

    /**
      * pointer to the reduction scheme. 
      */
	ogsChem::chemReductionGIA* _ReductionGIA;

	MathLib::StepperBulischStoer<Local_ODE_Xi_immob_GIA>* _sbs;

	//pointer to the local problem
    LocalProblem* _pSolve;

    /**
      * concentrations vector
      * including all components in the MCP data structure
      */
	std::vector<MyNodalFunctionScalar*> _concentrations;

    /**
      * nodal eta_mobile values
      */ 
    std::vector<MyNodalFunctionScalar*> _eta;

    /**
      * nodal eta_immobile values
      */ 
	std::vector<MyNodalFunctionScalar*> _eta_bar;
    
    /**
      * nodal xi_global values at the current time step
      */ 
    std::vector<MyNodalFunctionScalar*> _xi_global_cur;

    /**
      * nodal xi_global values at the previous time step
      */ 
    std::vector<MyNodalFunctionScalar*> _xi_global_pre;

    /**
      * nodal xi_local values
      */ 
    std::vector<MyNodalFunctionScalar*> _xi_local_old;

    /**
      * nodal xi_immobile values
      */ 
	std::vector<MyNodalFunctionScalar*> _xi_local_new;
	
	/**
	  * RHS rate values for nonlinear PDE of xi_sorp
	  */
	std::vector<MyNodalFunctionScalar*> _xi_sorp_rates;

	/**
	  * RHS rate values for nonlinear PDE of xi_min
	  */
	std::vector<MyNodalFunctionScalar*> _xi_min_rates;

	/**
	  * RHS rate values for nonlinear PDE of xi_kin
	  */
	std::vector<MyNodalFunctionScalar*> _xi_kin_rates;

	/**
	  * index on each node, whether each mineral is present or not
	  */
	std::vector<ogsChem::LocalVector*> _nodal_vec_AI;

    /**
      * global reaction rate vector
      */
    std::vector<MyNodalFunctionScalar*> _global_vec_Rate;

	/**
      * global _drates_dxi vector
      */
    std::vector<MyNodalFunctionScalar*> _drates_dxi;

    /**
      * degree of freedom equation ID talbe for the nonlinear problem
      */ 
    DiscreteLib::DofEquationIdTable * _nl_sol_dofManager;

    /**
      * the nested ode local problem
      */
    Local_ODE_Xi_immob_GIA*        _local_ode_xi_immob_GIA;

    /**
      * the id of the msh applied in this process
      */ 
    std::size_t _msh_id;

	std::size_t _n_xi_local,  _n_xi_global, _n_eta, _n_eta_bar; //TODO uninitialized

	//TODO no need to make them member variables
	std::size_t _n_xi_Mob, _n_xi_Sorp_tilde,_n_xi_Sorp, _n_xi_Sorp_bar, _n_xi_Sorp_bar_li, _n_xi_Sorp_bar_ld, _n_xi_Min, _n_xi_Min_tilde, _n_xi_Min_bar, _n_xi_Kin, _n_xi_Kin_bar, _I_mob, _I_sorp, _I_min, _I_kin, _J_tot_kin;
    std::size_t _n_Comp;
}; 

#include "ReductConc.hpp"

#endif  // end of ifndef
