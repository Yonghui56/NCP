/**
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ReductConc.h
 *
 * Created on    2013-03-18 by Haibing Shao
 */

#ifndef OPS_CONC_H
#define OPS_CONC_H

#include "BaseLib/CodingTools.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerResidualLocalAssembler.h"
#include "NumLib/Function/DiscreteDataConvergenceCheck.h"
#include "SolutionLib/Fem/FemIVBVProblem.h"
#include "SolutionLib/Fem/TemplateFemEquation.h"
#include "SolutionLib/Fem/SingleStepFEM.h"
#include "ProcessLib/AbstractTransientProcess.h"
#include "MathLib/DataType.h"
#include "ChemLib/chemEqReactSysActivity.h"
#include "NumLib/Nonlinear/DiscreteNRSolverWithStepInitFactory.h"
#include "UserModules/FemKinReactGIA/LinearTransportTimeODELocalAssember.h"
#include "UserModules/FemKinReactGIA/LinearTransportJacobianLocalAssembler.h"
#include "UserModules/FemKinReactGIA/NestedOdeNRIterationStepInitializer.h"
#include "FemReactOPSsolution.h"
#include "SingleStepReactOPS.h"

template <class T_DISCRETE_SYSTEM, class T_LINEAR_SOLVER >
class FunctionOPSConc
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

	typedef FunctionOPSConc        MyFunctionData;    // FEM problem class
    typedef T_DISCRETE_SYSTEM      MyDiscreteSystem;  // Discretization
    typedef T_LINEAR_SOLVER        MyLinearSolver;    // linear solver
	
    // memory for discretized concentration vector
    typedef typename FemLib::FemNodalFunctionVector<MyDiscreteSystem>::type MyNodalFunctionVector;
	typedef typename FemLib::FemNodalFunctionScalar<MyDiscreteSystem>::type MyNodalFunctionScalar;

    // local assembler
	// for the linear systems, use the same settings as Mass_Transport
    typedef LinearTransportTimeODELocalAssembler<NumLib::ElementWiseTimeEulerEQSLocalAssembler> MyLinearAssemblerType;          
    typedef LinearTransportTimeODELocalAssembler<NumLib::ElementWiseTimeEulerResidualLocalAssembler> MyLinearResidualAssemblerType; 
    typedef LinearTransportJacobianLocalAssembler MyLinearJacobianAssemblerType;                                        
	// for the nonlinear part, use different settings
	typedef NonLinearReactiveTransportTimeODELocalAssembler<NumLib::ElementWiseTimeEulerEQSLocalAssembler, MyNodalFunctionScalar>      MyNonLinearAssemblerType; 
	typedef NonLinearReactiveTransportTimeODELocalAssembler<NumLib::ElementWiseTimeEulerResidualLocalAssembler, MyNodalFunctionScalar> MyNonLinearResidualAssemblerType; 
	typedef NonLinearReactiveTransportJacobianLocalAssembler<MyNodalFunctionScalar, MyFunctionData> MyNonLinearJacobianAssemblerType; 
	typedef NestedOdeNRIterationStepInitializer<MyNodalFunctionScalar, MyFunctionData> MyNRIterationStepInitializer;
	typedef NumLib::DiscreteNRSolverWithStepInitFactory<MyNRIterationStepInitializer> MyDiscreteNonlinearSolverFactory; 
	
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
	typedef SolutionLib::TemplateFemEquation<
		    MyDiscreteSystem,
			MyLinearSolver,
			MyNonLinearAssemblerType,
			MyNonLinearResidualAssemblerType,
			MyNonLinearJacobianAssemblerType
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
	typedef SolutionLib::SingleStepFEM<
			MyNonLinearReactiveTransportProblemType,
			MyLinearSolver, 
			MyDiscreteNonlinearSolverFactory
			> MyNonLinearSolutionType;

    /**
	  * the general reduction problem part
	  */
	typedef SolutionLib::FemReactOPSsolution<MyDiscreteSystem> MyReactOPSProblemType; 
	typedef SolutionLib::SingleStepReactOPS<MyFunctionData, 
                                            MyReactOPSProblemType,
                                            MyLinearTransportProblemType, 
                                            MyLinearSolutionType, 
                                            MyNonLinearReactiveTransportProblemType, 
                                            MyNonLinearSolutionType> MyReactOPSSolution; 
    typedef typename MyReactOPSProblemType::MyVariable MyVariableConc;

    FunctionOPSConc() 
        : Process("REACT_TRANS_OPS", 1, 1),
          _feObjects(0)
    {
        // set default parameter name
		ProcessLib::Process::setInputParameterName(Velocity, "Velocity");
    };
	
	/**
	  * destructor, reclaim the memory
	  */
    virtual ~FunctionOPSConc()
    {
        BaseLib::releaseObject(_problem); 
        BaseLib::releaseObjectsInStdVector(_linear_problems);
        BaseLib::releaseObjectsInStdVector(_concentrations);
        BaseLib::releaseObject(_local_eq_react_sys); 
        BaseLib::releaseObject(_local_kin_react_sys); 
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
        return getSolution()->suggestNext(time_current); 
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
    virtual void finalizeTimeStep(const NumLib::TimeStep &time)
    {
         output(time);
		 getSolution()->finalizeTimeStep(time);
    };

	/**
      * update the change of kinetic reaction rate over the change of xi_mob
      */
	void set_concentrations(size_t conc_idx, MyNodalFunctionScalar* new_conc_nodal_values)
    {
        size_t node_idx; 
        double nodal_conc; 

        for (node_idx = _concentrations[0]->getDiscreteData()->getRangeBegin();
             node_idx < _concentrations[0]->getDiscreteData()->getRangeEnd(); 
             node_idx++ )
        {
            nodal_conc = new_conc_nodal_values->getValue(node_idx); 
            if ( nodal_conc < 0.0 )
                nodal_conc = 1.0e-36; 
            
            this->_concentrations[conc_idx]->setValue( node_idx, nodal_conc );        
        }
    
    }; 

    MyNodalFunctionScalar* get_concentrations(size_t idx_conc)
    {
        return _concentrations[idx_conc];
    };

	/**
	  * set boundary node values to the concentration vector
	  */
	void set_BC_conc_node_values(std::size_t node_idx, std::size_t i_var, double node_value); 

	/**
      * calculate nodal chemistry system
      */
	void calc_nodal_react_sys(double dt);

	/**
	  * calculate nodal equilibrium reaction system
	  */
	void calc_nodal_kin_react_sys(double dt);

	/**
	  * return the number of chemical components
	  */
	std::size_t get_n_Comp(void) { return _n_Comp; };

protected:
    virtual void initializeTimeStep(const NumLib::TimeStep &time);

    /**
      * this function is called to exchange output parameters
      */ 
    virtual void updateOutputParameter(const NumLib::TimeStep &time);

    /**
      * get the pointer of solution class for current problem
      */ 
    virtual MyReactOPSSolution* getSolution() {return _solution;};

    /**
      * output the result of current solution
      */ 
    virtual void output(const NumLib::TimeStep &time);

private:
    DISALLOW_COPY_AND_ASSIGN(FunctionOPSConc);

private:
    /**
      * linear problems
      */
	std::vector<MyLinearTransportProblemType*> _linear_problems;

    /**
      * linear solutions
      */
    std::vector<MyLinearSolutionType*>         _linear_solutions;


	/**
      * the local equilibrium reactions system
      */ 
	ogsChem::chemEqReactSysActivity*                   _local_eq_react_sys; 

    /**
      * the local kinetic reactions system
      */
    ogsChem::chemKinReactSys*                          _local_kin_react_sys;

	/**
      * reactive transport operator splitting problem
      */ 
	MyReactOPSProblemType* _problem; 

    /** 
      * reduction solution
      */
	MyReactOPSSolution*    _solution; 

	/**
      * FEM object
      */
    FemLib::LagrangeFeObjectContainer* _feObjects;
    
	/**
      * convergence checker
      */
	NumLib::DiscreteDataConvergenceCheck _checker;

    /**
      * concentrations vector
      * including all components in the MCP data structure
      */
	std::vector<MyNodalFunctionScalar*> _concentrations; 

    /**
      * the id of the msh applied in this process
      */ 
    size_t _msh_id; 

    /**
      * number of chemical components
      */
    size_t _n_Comp; 

    /**
      * number of mobile components
      */
    size_t _n_Comp_mob; 


}; 

#include "OPSConc.hpp"

#endif  // end of ifndef
