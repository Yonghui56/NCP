
#pragma once

#include "MathLib/LinAlg/LinearEquations/LisMPILinearEquation.h"


namespace DiscreteLib
{
/**
 * 
 */
template<class T_SPARSITY_BUILDER>
class LisMPILinearSystem : public AbstractMeshBasedDiscreteLinearEquation
{
public:
    ///
    LisMPILinearSystem(MeshLib::IMesh &msh, MathLib::LisMPILinearEquation &linear_solver, DofEquationIdTable &dofManager)
        : AbstractMeshBasedDiscreteLinearEquation(msh, dofManager), _eqs(&linear_solver), _do_create_eqs(true)
    {

    };

    /// set prescribed dof
    void setPrescribedDoF(size_t dofId, std::vector<size_t> &list_discrete_pt_id, std::vector<double> list_prescribed_values)
    {
        //assert(_list_prescribed_dof_id.size()==0);
        _list_prescribed_dof_id.clear();
        _list_prescribed_values.clear();

//        const DofMap* dofMap = getDofMapManger()->getDofMap(dofId);
//        const size_t n = list_discrete_pt_id.size();
//        for (size_t i=0; i<n; i++) {
//            size_t pt_id = list_discrete_pt_id[i];
//            if (dofMap->isActiveDoF(pt_id)) {
//                _list_prescribed_dof_id.push_back(dofMap->getEqsID(pt_id));
//                _list_prescribed_values.push_back(list_prescribed_values[i]);
//            }
//        }
    }

    /// construct the linear equation
    void construct(IDiscreteLinearEquationAssembler& assemler)
    {
        DofEquationIdTable* dofManager = getDofMapManger();
        assert(dofManager->getNumberOfVariables()>0);

        if (_do_create_eqs) {
            _do_create_eqs = false;
            MathLib::RowMajorSparsity sparse;
            T_SPARSITY_BUILDER sp_builder(*getMesh(), *dofManager, sparse);
            _eqs->create(dofManager->getTotalNumberOfActiveDoFs(), &sparse);
        } else {
            _eqs->reset();
        }
        assemler.assembly(*getMesh(), *dofManager, *_eqs);

        //apply 1st bc
        _eqs->setKnownX(_list_prescribed_dof_id, _list_prescribed_values);
    }

    /// set additional RHS values
    void addRHS(size_t dofId, std::vector<size_t> &list_discrete_pt_id, std::vector<double> list_rhs_values, double fkt=1.0)
    {
//        const DofMap* dofMap = getDofMapManger()->getDofMap(dofId);
//        const size_t n = list_discrete_pt_id.size();
//        for (size_t i=0; i<n; i++) {
//            size_t pt_id = list_discrete_pt_id[i];
//            if (dofMap->isActiveDoF(pt_id)) {
//                _eqs->addRHS(dofMap->getEqsID(pt_id), list_rhs_values[i]*fkt);
//            }
//        }
    }


    /// solve
    void solve()
    {
        _eqs->solve();
    }


    /// get the solution vector
    double* getLocalX()
    {
        return _eqs->getX();
    }

    void getGlobalX(std::vector<double> &x)
    {
        x.resize(_eqs->getDimension());
        double *tmp_x = _eqs->getX();
        for (size_t i=0; i<x.size(); i++)
            x[i] = tmp_x[i];
    }

private:
    MathLib::LisMPILinearEquation* _eqs;
    bool _do_create_eqs;
    std::vector<size_t> _list_prescribed_dof_id;
    std::vector<double> _list_prescribed_values;

    DISALLOW_COPY_AND_ASSIGN(LisMPILinearSystem);
};

} //end
