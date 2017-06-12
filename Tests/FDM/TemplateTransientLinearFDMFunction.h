/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TemplateTransientLinearFDMFunction.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "NumLib/Function/IFunction.h"
#include "DiscreteLib/Serial/DiscreteSystem.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "StencilWiseTransientLinearEQSAssembler.h"
#include "FdmFunction.h"
#include "BoundaryConditions.h"


namespace FdmLib
{

typedef DiscreteLib::IDiscreteVector<double> MyFemVector;

/**
 * \brief Template class for transient linear FEM functions
 *
 * \tparam T_FEM_PROBLEM
 * \tparam T_TIME_ODE_ASSEMBLER
 * \tparam T_SPACE_ASSEMBLER
 * \tparam T_LINEAR_SOLVER
 */
template <
    class T_USER_FEM_PROBLEM,
    class T_LOCAL_ASSEMBLER
    >
class TemplateTransientLinearFDMFunction
    : public NumLib::TemplateFunction<MyFemVector, MyFemVector>
{
public:
    typedef T_USER_FEM_PROBLEM UserFemProblem;
    typedef T_LOCAL_ASSEMBLER UserLocalAssembler;

    /// constructor
    /// @param problem        Fem problem
    /// @param linear_eqs    Discrete linear equation
    TemplateTransientLinearFDMFunction(UserFemProblem* problem, UserLocalAssembler* asssembler, DiscreteLib::IDiscreteLinearEquation* linear_eqs)
        : _problem(problem), _local_assembler(asssembler),  _linear_eqs(linear_eqs),
          _t_n1(0), _u_n0(0)
    {
    };

    ///
    virtual ~TemplateTransientLinearFDMFunction() {};

    ///
    NumLib::TemplateFunction<MyFemVector,MyFemVector>* clone() const
    {
        return new TemplateTransientLinearFDMFunction<
                    T_USER_FEM_PROBLEM,
                    T_LOCAL_ASSEMBLER
                    >(_problem, _local_assembler, _linear_eqs);
    }

    /// reset property
    void reset(const NumLib::TimeStep* t, MyFemVector* u_n0)
    {
        this->_t_n1 = const_cast<NumLib::TimeStep*>(t);
        this->_u_n0 = u_n0;
    };

    /// solve linear equations discretized with FEM
    /// @param u0    initial guess
    /// @param u_n1 new results
    void eval(const MyFemVector &u0, MyFemVector &u_n1)
    {
        // input, output
        const NumLib::TimeStep &t_n1 = *this->_t_n1;
        MyFemVector* u_n = this->_u_n0;

        // prepare data
        UserFemProblem* pro = _problem;

        //
        _linear_eqs->initialize(); // should be called before BC

        // setup BC
        for (size_t i=0; i<pro->getNumberOfDirichletBC(); i++) {
            FdmLib::FdmDirichletBC<double> *bc1 = pro->getFdmDirichletBC(i);
            bc1->setup();
            size_t varid = 0; //?
            _linear_eqs->setPrescribedDoF(varid, bc1->getListOfBCNodes(), bc1->getListOfBCValues());
        }
        for (size_t i=0; i<pro->getNumberOfNeumannBC(); i++)
            pro->getFdmNeumannBC(i)->setup();

        //TODO temporally
        std::vector<MyFemVector*> vec_un;
        vec_un.push_back(const_cast<MyFemVector*>(u_n));
        std::vector<MyFemVector*> vec_un1;
        vec_un1.push_back(const_cast<MyFemVector*>(&u_n1));

        // assembly
        StencilWiseTransientLinearEQSAssembler assembler(&t_n1, &vec_un, &vec_un1, _local_assembler);
        _linear_eqs->construct(assembler);

        //apply BC1,2
        for (size_t i=0; i<pro->getNumberOfNeumannBC(); i++) {
            FdmLib::FdmNeumannBC<double, double> *bc2 = pro->getFdmNeumannBC(i);
            size_t varid = 0; //?
            _linear_eqs->addRHS(varid, bc2->getListOfBCNodes(), bc2->getListOfBCValues(), -1.0);
        }

        // solve
        _linear_eqs->setX(u0);
        _linear_eqs->solve();
        _linear_eqs->getX(u_n1);
    }


private:
    UserFemProblem* _problem;
    UserLocalAssembler *_local_assembler;
    DiscreteLib::IDiscreteLinearEquation* _linear_eqs;
    NumLib::TimeStep* _t_n1;
    MyFemVector* _u_n0;
};


} //end
