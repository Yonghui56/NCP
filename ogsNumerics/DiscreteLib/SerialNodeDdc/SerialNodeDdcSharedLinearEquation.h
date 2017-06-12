/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NodeDDCSerialSharedDiscreteSystem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/BidirectionalMap.h"
#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/Serial/DiscreteSystem.h"
#include "DiscreteLib/Serial/DiscreteLinearEquation.h"
#include "DiscreteLib/DDC/DecomposedDomain.h"
#include "DiscreteLib/DDC/SubDomain.h"
#include "DiscreteLib/Utils/SparsityBuilderDDC.h"
#include "DecomposedVector.h"

namespace DiscreteLib
{

template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
class SerialNodeDdcSharedLinearEquation : public IDiscreteLinearEquation
{
private:
    SerialNodeDdcSharedLinearEquation(IDiscreteSystem &sys, DecomposedDomain* ddc_global, T_LINEAR_SOLVER* shared_linear_solver, DofEquationIdTable* global_dofManager)
    : _ddc_global(ddc_global), _sheared_eqs(shared_linear_solver), _global_dofManager(global_dofManager)
    {
        for (size_t i=0; i<ddc_global->getNumberOfSubDomains(); i++) {
            SubDomain* dom = ddc_global->getSubDomain(i);
            DofEquationIdTable *local_dofManager = new DofEquationIdTable();
            _list_local_eq.push_back(DiscreteLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>::createInstance(sys, dom->getLoalMesh(), shared_linear_solver, local_dofManager));
        }
        _do_create_eqs = true;
    };

public:
    static SerialNodeDdcSharedLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>*
    createInstance(IDiscreteSystem &sys, DecomposedDomain* ddc_global, T_LINEAR_SOLVER* shared_linear_solver, DofEquationIdTable* global_dofManager)
    {
        SerialNodeDdcSharedLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>* eqs;
        eqs = new SerialNodeDdcSharedLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>(sys, ddc_global, shared_linear_solver, global_dofManager);
        sys.addLinearEquation(eqs);
        return eqs;
    }

    virtual ~SerialNodeDdcSharedLinearEquation() 
    {
        //BaseLib::releaseObjectsInStdVector(_list_local_eq);
    };

    /// initialize EQS
    void initialize()
    {
        DofEquationIdTable* global_dofManager = getDofMapManger();
        if (_do_create_eqs) {
            _do_create_eqs = false;
            //create global linear equation
            MathLib::RowMajorSparsity global_sparse;
            SparsityBuilderFromDDC<T_SPARSITY_BUILDER> sp_builder(*_ddc_global, *global_dofManager, global_sparse);
            _sheared_eqs->create(global_dofManager->getTotalNumberOfActiveDoFs(), &global_sparse);
            //local dof
            size_t mesh_id = _ddc_global->getID();
            for (size_t i=0; i<_list_local_eq.size(); i++) {
                SubDomain* dom = _ddc_global->getSubDomain(i);
                MeshLib::IMesh* local_msh = dom->getLoalMesh();
                DofEquationIdTable *local_dofManager = _list_local_eq[i]->getDofMapManger();
                // for each var
                for (size_t j=0; j<global_dofManager->getNumberOfVariables(); j++) {
                    local_dofManager->addVariableDoFs(mesh_id, 0, local_msh->getNumberOfNodes());
                    IEquationIdStorage* local_pt2dof = local_dofManager->getPointEquationIdTable(j, mesh_id);
                    for (size_t k=0; k<local_msh->getNumberOfNodes(); k++) {
                        size_t dof_id = global_dofManager->mapEqsID(mesh_id, j, dom->getGlobalLocalIdMap()->local2global(k));
                        local_pt2dof->set(k, dof_id);
                    }
                }
                // ghost
                if (dom->getGhostList()!=0) {
                    std::set<size_t>* dom_ghosts = dom->getGhostList();
                    std::vector<size_t> vec_ghosts(dom_ghosts->begin(), dom_ghosts->end());
                    local_dofManager->setGhostPoints(mesh_id, vec_ghosts);
                }
            }
        } else {
            _sheared_eqs->reset();
        }
    }

    /// set prescribed dof
    void setPrescribedDoF(size_t dofId, std::vector<size_t> &list_discrete_pt_id, std::vector<double> &list_prescribed_values)
    {
        //assert(_list_prescribed_dof_id.size()==0);
        _list_prescribed_dof_id.clear();
        _list_prescribed_values.clear();

        const size_t n = list_discrete_pt_id.size();
        for (size_t i=0; i<n; i++) {
            size_t pt_id = list_discrete_pt_id[i];
            if (_global_dofManager->isActiveDoF(dofId, 0, pt_id)) {
                _list_prescribed_dof_id.push_back(_global_dofManager->mapEqsID(dofId, 0, pt_id));
                _list_prescribed_values.push_back(list_prescribed_values[i]);
            }
        }
    }

    /// construct 
    void construct(IDiscreteLinearEquationAssembler& assemler) 
    {
        //DofEquationIdTable* dofManager = getDofMapManger();

        // local eqs
        for (size_t i=0; i<_list_local_eq.size(); i++) {
            _list_local_eq[i]->construct(assemler);
            //_global_eqs->printout();
        }

        //apply 1st bc
        _sheared_eqs->setKnownX(_list_prescribed_dof_id, _list_prescribed_values);
        //_global_eqs->printout();
    }

    /// solve
    void solve()
    {
        _sheared_eqs->solve();
    }

    /// get solution
    double* getLocalX()
    {
        return _sheared_eqs->getX();
    }

    void getGlobalX(std::vector<double> &x)
    {
        x.resize(_sheared_eqs->getDimension());
        double *tmp_x = _sheared_eqs->getX();
        for (size_t i=0; i<x.size(); i++)
            x[i] = tmp_x[i];
    }
    void getX(IDiscreteVector<double> &x)
    {
        double *tmp_x = _sheared_eqs->getX();
        for (size_t i=x.getRangeBegin(); i<x.getRangeEnd(); i++)
            x[i] = tmp_x[i];
    };
    ///
    virtual void setX(const GlobalVectorType &x)
    {
        double *tmp_x = _sheared_eqs->getX();
        for (size_t i=x.getRangeBegin(); i<x.getRangeEnd(); i++)
            tmp_x[i] = x[i];
    };
    /// get a Dof map manager
    DofEquationIdTable* getDofMapManger() const
    {
        return _global_dofManager;
    }
    /// set additional RHS values
    void addRHS(size_t dofId, std::vector<size_t> &list_discrete_pt_id, std::vector<double> &list_rhs_values, double fkt)
    {
        const size_t n = list_discrete_pt_id.size();
        for (size_t i=0; i<n; i++) {
            size_t pt_id = list_discrete_pt_id[i];
            if (_global_dofManager->isActiveDoF(dofId, 0, pt_id)) {
                _sheared_eqs->addRHS(_global_dofManager->mapEqsID(dofId, 0, pt_id), list_rhs_values[i]*fkt);
            }
        }
    }

    virtual void addRHS(const GlobalVectorType &v, double fkt=1.0)
    {
        const size_t n = v.size();
        for (size_t i=0; i<n; i++) {
            _sheared_eqs->addRHS(i, v[i]*fkt);
        }
    }

private:
    DecomposedDomain* _ddc_global;
    std::vector<AbstractMeshBasedDiscreteLinearEquation*> _list_local_eq;
    T_LINEAR_SOLVER* _sheared_eqs;
    bool _do_create_eqs;
    DofEquationIdTable* _global_dofManager;
    std::vector<size_t> _list_prescribed_dof_id;
    std::vector<double> _list_prescribed_values;

    DISALLOW_COPY_AND_ASSIGN(SerialNodeDdcSharedLinearEquation);
};

} //end
