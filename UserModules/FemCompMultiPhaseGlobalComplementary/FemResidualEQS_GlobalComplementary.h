/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemResidualEQS.h
 *
 * Created on 2013-08-016 by Yonghui Huang & Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "DiscreteLib/Utils/Tools.h"
#include "NumLib/Function/IFunction.h"
#include "NumLib/TransientAssembler/TransientElementWiseVectorUpdater.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "FemLib/Function/FemFunction.h"

#include "SolutionLib/Fem/FemVariable.h"
#include "SolutionLib/DataType.h"

//#include "UserModules/FemReactGIAReduct/ReductConc.h"
//#include "ReductionGIANodeInfo.h"


/**
 * \brief Template class for FEM residual functions
 *
 * \tparam T_LOCAL_RESIDUAL_ASSEMBLER
 */
template <class T_DIS_SYS, class T_FUNCTION_DATA, class T_LOCAL_RESIDUAL_ASSEMBLER>
class TemplateTransientResidualFEMFunction_Global_Complementary
    : public NumLib::TemplateFunction<SolutionLib::SolutionVector, SolutionLib::SolutionVector>
{
public:
    typedef T_DIS_SYS MyDiscreteSystem;
    typedef SolutionLib::FemVariable MyVariable;
    typedef typename FemLib::FemNodalFunctionScalar<MyDiscreteSystem>::type MyNodalFunctionScalar;
	typedef T_LOCAL_RESIDUAL_ASSEMBLER UserLocalResidualAssembler;
	typedef typename NumLib::TransientElementWiseVectorUpdater<UserLocalResidualAssembler> MyUpdater;
	typedef typename T_DIS_SYS::template MyVectorAssembler<double, MyUpdater>::type MyGlobalAssembler;

    /// constructor
    /// @param problem        Fem problem
    /// @param linear_eqs    Discrete linear equation
	TemplateTransientResidualFEMFunction_Global_Complementary(MyDiscreteSystem* sys,
    		                                        const std::vector<MyVariable*> &list_var,
    		                                        DiscreteLib::DofEquationIdTable* dofManager,
													UserLocalResidualAssembler* asssembler,
    		                                        T_FUNCTION_DATA* function_data)
													: _dis_sys(sys), _dofManager(dofManager), _local_assembler(asssembler),
          _time_step(0), _st(0), _u_n0(0), _list_var(list_var), _function_data(function_data)
          
    {
    };

    ///
	virtual ~TemplateTransientResidualFEMFunction_Global_Complementary()
    {
        _function_data	= NULL;
    };


    ///
    NumLib::TemplateFunction<SolutionLib::SolutionVector,SolutionLib::SolutionVector>* clone() const
    {
		return new TemplateTransientResidualFEMFunction_Global_Complementary
			<T_DIS_SYS, T_FUNCTION_DATA, T_LOCAL_RESIDUAL_ASSEMBLER>(_dis_sys, _list_var, _dofManager, _local_assembler, _function_data);
    }

    /// reset property
    void reset(const NumLib::TimeStep& t, const SolutionLib::SolutionVector* u_n0, SolutionLib::SolutionVector* st)
    {
        _time_step = &t;
        _u_n0 = u_n0;
        _st = st;
    };

    /// evaluate residual
    /// @param u_n1    current results
    /// @param r residual
    void eval(const SolutionLib::SolutionVector &u_n1, SolutionLib::SolutionVector &r);

	const NumLib::TimeStep* getTimeStepObj(void)
	{
		return _time_step;
	}

private:
	void GlobalResidualAssembler(T_FUNCTION_DATA * _function_data,
                                 SolutionLib::SolutionVector & residual_global);

private:
    MyDiscreteSystem *_dis_sys;
    DiscreteLib::DofEquationIdTable* _dofManager;
	UserLocalResidualAssembler *_local_assembler;
    const NumLib::TimeStep* _time_step;
    const SolutionLib::SolutionVector* _st;
    const SolutionLib::SolutionVector *_u_n0;
    std::vector<MyVariable*> _list_var;
    T_FUNCTION_DATA* _function_data;
	MathLib::LocalVector _vec_Res;
    //TODO pass via constructor
    //
    // std::map<size_t, ReductionGIANodeInfo*>* _bc_info;
    NumLib::ITXFunction* _vel;
    FemLib::IFemNumericalIntegration* _q;
    FemLib::IFiniteElement* _fe;
    /**
      * pointer to the activity model
      */
   
};


template <class T_DIS_SYS, class T_FUNCTION_DATA, class T_LOCAL_RESIDUAL_ASSEMBLER>
void TemplateTransientResidualFEMFunction_Global_Complementary<T_DIS_SYS, T_FUNCTION_DATA, T_LOCAL_RESIDUAL_ASSEMBLER>::eval(const SolutionLib::SolutionVector &u_n1, SolutionLib::SolutionVector &r)
{
    // input, output
	MeshLib::IMesh* msh = _dis_sys->getMesh();
    size_t msh_id = _dis_sys->getMesh()->getID();
	const SolutionLib::SolutionVector *u_n = this->_u_n0;
	const NumLib::TimeStep &t_n1 = *this->_time_step;
    // assembly
    r = .0;
    //node based operations
	// calculate the global residual
    //GlobalResidualAssembler(*this->_time_step, u_n1, *_u_n0, r);
	MyUpdater updater(&t_n1, msh, _dofManager, u_n, &u_n1, _local_assembler);
	MyGlobalAssembler assembler(&updater);
	r.construct(*_dofManager, assembler);
	GlobalResidualAssembler(_function_data, r);
    // add source/sink term (i.e. RHS in linear equation)
    if (_st!=0) {
        r += *_st;
    }

    // set residuals to zero for Dirichlet BC
    std::vector<size_t> list_bc1_eqs_id;
    std::vector<double> list_bc1_val;
    for (size_t i=0; i<_list_var.size(); i++) {
        MyVariable* var = _list_var[i];
        std::vector<size_t> var_bc_id;
        std::vector<double> var_bc_val;
        for (size_t j=0; j<var->getNumberOfDirichletBC(); j++) {
        	SolutionLib::FemDirichletBC* bc1 = var->getDirichletBC(j);
            bc1->setup(var->getCurrentOrder());
            DiscreteLib::convertToEqsValues(*_dofManager, i, msh_id, bc1->getListOfBCNodes(), bc1->getListOfBCValues(), list_bc1_eqs_id, list_bc1_val);
        }
    }
    r.setSubvector(list_bc1_eqs_id, .0);
}

template <class T_DIS_SYS, class T_FUNCTION_DATA, class T_LOCAL_RESIDUAL_ASSEMBLER>
void TemplateTransientResidualFEMFunction_Global_Complementary<T_DIS_SYS, T_FUNCTION_DATA, T_LOCAL_RESIDUAL_ASSEMBLER>
::GlobalResidualAssembler(T_FUNCTION_DATA * _function_data,
                               SolutionLib::SolutionVector & residual_global)
{
	
	const size_t nnodes = _dis_sys->getMesh()->getNumberOfNodes();
	_vec_Res = MathLib::LocalVector::Zero(3 * nnodes, 1);
	_vec_Res = _function_data->get_vec_Res();
	for (size_t node_idx = 0; node_idx < nnodes; node_idx++){
		residual_global[3 * node_idx+2] = -_vec_Res[2 * nnodes+node_idx];
		//residual_global[2 * nnodes + node_idx] = -_vec_Res[2 * nnodes + node_idx];
	}
}



