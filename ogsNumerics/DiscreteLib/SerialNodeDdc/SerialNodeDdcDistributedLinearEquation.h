/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NodeDDCSerialDistributedDiscreteSystem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
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
class SerialNodeDdcDistributedLinearEquation : public IDiscreteLinearEquation
{
private:
    //TODO is it necessary to pass instance of linear solvers? passing settings seems to be enough.
    SerialNodeDdcDistributedLinearEquation(IDiscreteSystem &sys, DecomposedDomain* ddc_global, T_LINEAR_SOLVER* global_linear_solver, DofEquationIdTable* global_dofManager)
    : _ddc_global(ddc_global), _global_linear_solver(global_linear_solver), _do_create_eqs(true), _global_dofManager(global_dofManager)
    {
        size_t mesh_id = ddc_global->getID();
        for (size_t i=0; i<ddc_global->getNumberOfSubDomains(); i++) {
            SubDomain* dom = ddc_global->getSubDomain(i);
            MeshLib::IMesh* local_msh = dom->getLoalMesh();
            // create local linear solver
            T_LINEAR_SOLVER* local_solver = new T_LINEAR_SOLVER; //TODO option
            // create local DoFmap
            DofEquationIdTable *local_dofManager = new DofEquationIdTable();
            for (size_t j=0; j<global_dofManager->getNumberOfVariables(); j++) {
                local_dofManager->addVariableDoFs(mesh_id, 0, local_msh->getNumberOfNodes());
            }
            if (dom->getGhostList()!=0) {
                std::set<size_t>* dom_ghosts = dom->getGhostList();
                std::vector<size_t> vec_ghosts(dom_ghosts->begin(), dom_ghosts->end());
                local_dofManager->setGhostPoints(mesh_id, vec_ghosts);
            }
            local_dofManager->setNumberingType(DofNumberingType::BY_POINT);
            local_dofManager->construct();
            // create local EQS
            _list_local_eq.push_back(DiscreteLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>::createInstance(sys, dom->getLoalMesh(), local_solver, local_dofManager));
        }
    };

public:
    static SerialNodeDdcDistributedLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>* createInstance(IDiscreteSystem &sys, DecomposedDomain* ddc_global, T_LINEAR_SOLVER* global_linear_solver, DofEquationIdTable* global_dofManager)
    {
        SerialNodeDdcDistributedLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>* eqs;
        eqs = new SerialNodeDdcDistributedLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>(sys, ddc_global, global_linear_solver, global_dofManager);
        sys.addLinearEquation(eqs);
        return eqs;
    }

    virtual ~SerialNodeDdcDistributedLinearEquation() 
    {
        //BaseLib::releaseObjectsInStdVector(_list_local_eq);
    };

    /// initialize EQS
    void initialize()
    {
        //init local 
        for (size_t i=0; i<_list_local_eq.size(); i++) {
            _list_local_eq[i]->initialize();
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
        // local eqs
        for (size_t i=0; i<_list_local_eq.size(); i++) {
            _list_local_eq[i]->construct(assemler);
        }
        //apply 1st bc
    }

    /// solve
    void solve()
    {
        const size_t mesh_id = _ddc_global->getID();
        // prepare global eqs
        DofEquationIdTable* global_dofManager = getDofMapManger();
        if (_do_create_eqs) {
            _do_create_eqs = false;
            //create global linear equation
            MathLib::RowMajorSparsity global_sparse;
            SparsityBuilderFromDDC<T_SPARSITY_BUILDER> sp_builder(*_ddc_global, *global_dofManager, global_sparse);
            _global_linear_solver->create(global_dofManager->getTotalNumberOfActiveDoFs(), &global_sparse);
        } else {
            _global_linear_solver->reset();
        }
        // merge local eqs into a global eqs
        MathLib::ILinearEquation* global_solver = _global_linear_solver;
        for (size_t i=0; i<_list_local_eq.size(); i++) {
            SubDomain* dom = _ddc_global->getSubDomain(i);
            IGlobaLocalMappingTable* pt_id_mapping = dom->getGlobalLocalIdMap();
            DiscreteLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>* local_eqs = _list_local_eq[i];
            DofEquationIdTable* local_dofTable = local_eqs->getDofMapManger();
            MathLib::ILinearEquation* local_solver = local_eqs->getLinearSolver();
            //local_solver->printout();
            MathLib::RowMajorSparsity* local_sp = local_eqs->getSparsity();
            const size_t n_local_pt = dom->getLoalMesh()->getNumberOfNodes();
            const size_t n_var = global_dofManager->getNumberOfVariables(mesh_id);
            // for each point
            for (size_t i_pt=0; i_pt<n_local_pt; i_pt++) {
                if (local_dofTable->isGhostPoint(mesh_id, i_pt)) continue; //skip ghost
                const size_t i_global_pt = pt_id_mapping->local2global(i_pt);
                // for each var
                for (size_t i_var=0; i_var<n_var; i_var++) {
                    // eqs id
                    const size_t local_row_id = local_dofTable->mapEqsID(i_var, mesh_id, i_pt);
                    if (local_row_id==BaseLib::index_npos) continue;
                    const size_t global_row_id = global_dofManager->mapEqsID(i_var, mesh_id, i_global_pt);
                    // add RHS
                    global_solver->addRHS(global_row_id, local_solver->getRHS(local_row_id));
                    // add A
                    const std::set<size_t> &local_sp_row = (*local_sp)[local_row_id];
                    for (std::set<size_t>::const_iterator itr=local_sp_row.begin(); itr!=local_sp_row.end(); ++itr) {
                        const size_t local_column_id = *itr;
                        size_t tmp_var_id, tmp_mesh_id, tmp_pt_id;
                        local_dofTable->mapDoF(local_column_id, tmp_var_id, tmp_mesh_id, tmp_pt_id);
                        size_t global_column_id = global_dofManager->mapEqsID(tmp_var_id, tmp_mesh_id, pt_id_mapping->local2global(tmp_pt_id));
                        global_solver->addA(global_row_id, global_column_id, local_solver->getA(local_row_id, local_column_id)); 
                    }
                }
            }
        }
        //_global_linear_solver->printout();
        //apply 1st bc
        _global_linear_solver->setKnownX(_list_prescribed_dof_id, _list_prescribed_values);
        //_global_eqs->printout();
        // solve
        _global_linear_solver->solve();
    }

    /// get solution
    double* getLocalX()
    {
        return _global_linear_solver->getX();
    }

    void getGlobalX(std::vector<double> &x)
    {
        x.resize(_global_linear_solver->getDimension());
        double *tmp_x = _global_linear_solver->getX();
        for (size_t i=0; i<x.size(); i++)
            x[i] = tmp_x[i];
    }
    void getX(IDiscreteVector<double> &x)
    {
        double *tmp_x = _global_linear_solver->getX();
        for (size_t i=x.getRangeBegin(); i<x.getRangeEnd(); i++)
            x[i] = tmp_x[i];
    };
    ///
    virtual void setX(const GlobalVectorType &x)
    {
        double *tmp_x = _global_linear_solver->getX();
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
                _global_linear_solver->addRHS(_global_dofManager->mapEqsID(dofId, 0, pt_id), list_rhs_values[i]*fkt);
            }
        }
    }

    virtual void addRHS(const GlobalVectorType &v, double fkt=1.0)
    {
        const size_t n = v.size();
        for (size_t i=0; i<n; i++) {
            _global_linear_solver->addRHS(i, v[i]*fkt);
        }
    }

private:
    DecomposedDomain* _ddc_global;
    std::vector<DiscreteLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>*> _list_local_eq;
    T_LINEAR_SOLVER* _global_linear_solver;
    bool _do_create_eqs;
    DofEquationIdTable* _global_dofManager;
    std::vector<size_t> _list_prescribed_dof_id;
    std::vector<double> _list_prescribed_values;

    DISALLOW_COPY_AND_ASSIGN(SerialNodeDdcDistributedLinearEquation);
};

} //end
