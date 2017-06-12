/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DiscreteLinearEquation.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"

#include "MathLib/LinAlg/Sparse/Sparsity.h"

#include "MeshLib/Core/IMesh.h"

#include "DiscreteLib/Core/IDiscreteLinearEquation.h"
#include "DiscreteLib/Core/IDiscreteLinearEquationAssembler.h"
#include "DiscreteLib/Utils/DofEquationIdTable.h"
#include "DiscreteLib/Utils/SparsityBuilderDummy.h"

#include "AbstractMeshBasedDiscreteLinearEquation.h"

namespace DiscreteLib
{

/**
 * \brief Implementation of mesh based linear equation classes 
 *        combined with a specific linear solver class
 *
 * \tparam T_LINEAR_SOLVER      Linear solver class
 * \tparam T_SPARSITY_BUILDER   Sparsity builder
 * 
 * - linear equation
 */
template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER=SparsityBuilderDummy>
class DiscreteLinearEquation : public AbstractMeshBasedDiscreteLinearEquation
{
public:
    typedef T_LINEAR_SOLVER MyLinearSolverType;
    typedef T_SPARSITY_BUILDER MySparsityBuilder;

private:
    /// constructor
    /// \param msh
    /// \param linear_solver
    /// \param dofManager
    DiscreteLinearEquation(const MeshLib::IMesh* msh, MyLinearSolverType* linear_solver, DofEquationIdTable* dofManager)
        : AbstractMeshBasedDiscreteLinearEquation(msh, dofManager), _eqs(linear_solver), _do_create_eqs(true)
    {
    };

public:

    static DiscreteLinearEquation<MyLinearSolverType,MySparsityBuilder>* createInstance(IDiscreteSystem &dis_sys, const MeshLib::IMesh* msh, MyLinearSolverType* linear_solver, DofEquationIdTable* dofManager)
    {
        DiscreteLinearEquation<MyLinearSolverType,MySparsityBuilder> *eqs;
        eqs = new DiscreteLinearEquation<MyLinearSolverType,MySparsityBuilder>(msh, linear_solver, dofManager);
        dis_sys.addLinearEquation(eqs);
        return eqs;
    }

    /// initialize
    virtual void initialize()
    {
        DofEquationIdTable* dofManager = getDofMapManger();
        assert(dofManager->getNumberOfVariables()>0);

        if (_do_create_eqs && !_eqs->isCreated()) {
            _do_create_eqs = false;
            MathLib::RowMajorSparsity* sparse = new MathLib::RowMajorSparsity();
            MySparsityBuilder sp_builder(*getMesh(), *dofManager, *sparse);
            AbstractMeshBasedDiscreteLinearEquation::setSparsity(sparse);
            _eqs->create(dofManager->getTotalNumberOfActiveDoFs(), sparse);
        } else {
            _eqs->reset();
        }
        _list_prescribed_dof_id.clear();
        _list_prescribed_values.clear();
    }

    /// construct the linear equation
    virtual void construct(IDiscreteLinearEquationAssembler& assemler)
    {
        assert(getDofMapManger()->getNumberOfVariables()>0);

        assemler.assembly(*(MeshLib::IMesh*)getMesh(), *AbstractMeshBasedDiscreteLinearEquation::getDofMapManger(), *_eqs);
//        //apply 1st bc
//        _eqs->setKnownX(_list_prescribed_dof_id, _list_prescribed_values);
    }

    /// set prescribed dof
    virtual void setPrescribedDoF(size_t varId, std::vector<size_t> &list_discrete_pt_id, std::vector<double> &list_prescribed_values)
    {
        const size_t n = list_discrete_pt_id.size();
        const DofEquationIdTable* dofmap = getDofMapManger(); 
        const size_t msh_id = getMesh()->getID();
        for (size_t i=0; i<n; i++) {
            size_t pt_id = list_discrete_pt_id[i];
            if (dofmap->isActiveDoF(varId, msh_id, pt_id)) {
                _list_prescribed_dof_id.push_back(dofmap->mapEqsID(varId, msh_id, pt_id));
                _list_prescribed_values.push_back(list_prescribed_values[i]);
            }
        }
    }
   
    /// set additional RHS values
    virtual void addRHS(size_t dofId, std::vector<size_t> &list_discrete_pt_id, std::vector<double> &list_rhs_values, double fkt=1.0)
    {
        const size_t n = list_discrete_pt_id.size();
        const DofEquationIdTable* dofmap = getDofMapManger(); 
        const size_t msh_id = getMesh()->getID();
        for (size_t i=0; i<n; i++) {
            size_t pt_id = list_discrete_pt_id[i];
            if (dofmap->isActiveDoF(dofId, msh_id, pt_id)) {
                _eqs->addRHS(dofmap->mapEqsID(dofId, msh_id, pt_id), list_rhs_values[i]*fkt);
            }
        }
    }

    /// add RHS
    virtual void addRHS(const GlobalVectorType &v, double fkt=1.0)
    {
        const size_t n = v.size();
        for (size_t i=0; i<n; i++) {
            _eqs->addRHS(i, v[i]*fkt);
        }
    }

    /// solve
    virtual void solve()
    {
        //apply 1st bc
        _eqs->setKnownX(_list_prescribed_dof_id, _list_prescribed_values);
        _eqs->solve();
    }


    /// get the solution vector
    virtual double* getLocalX()
    {
        return _eqs->getX();
    }

    /// copy global vector of x
    virtual void getGlobalX(std::vector<double> &x)
    {
        x.resize(_eqs->getDimension());
        double *tmp_x = _eqs->getX();
        for (size_t i=0; i<x.size(); i++)
            x[i] = tmp_x[i];
    }

    /// get a global vector of x
    virtual void getX(GlobalVectorType &x)
    {
        double *tmp_x = _eqs->getX();
        for (size_t i=x.getRangeBegin(); i<x.getRangeEnd(); i++)
            x[i] = tmp_x[i];
    };

    ///
    virtual void setX(const GlobalVectorType &x)
    {
        double *tmp_x = _eqs->getX();
        for (size_t i=x.getRangeBegin(); i<x.getRangeEnd(); i++)
            tmp_x[i] = x[i];
    };

    /// return a linear solver object
    MyLinearSolverType* getLinearSolver() {return _eqs;};

private:
    DISALLOW_COPY_AND_ASSIGN(DiscreteLinearEquation);

private:
    MyLinearSolverType* _eqs;
    bool _do_create_eqs;
    std::vector<size_t> _list_prescribed_dof_id;
    std::vector<double> _list_prescribed_values;
	
};

} //end
