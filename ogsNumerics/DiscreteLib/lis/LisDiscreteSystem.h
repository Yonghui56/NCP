
#if defined(USE_LIS) && defined(USE_MPI)

#pragma once

#include <vector>
#include <algorithm>

#include "mpi.h"

#include "MathLib/LinAlg/LinearEquations/LisMPILinearEquation.h"

#include "DiscreteLib/Serial/DiscreteSystem.h"
#include "DiscreteLib/LinearEquation/MeshBasedDiscreteLinearEquation.h"
#include "LisMPIDiscreteVector.h"
#include "LisMPILinearSystem.h"


namespace DiscreteLib
{

/**
 * \brief Discrete system based on LIS and MPI
 */
class LisMPIDiscreteSystem : public IDiscreteSystem
{
public:
    LisMPIDiscreteSystem(MeshLib::IMesh &local_msh) : _local_msh(&local_msh)
    {
    };

    virtual ~LisMPIDiscreteSystem()
    {
        Base::releaseObjectsInStdVector(_list_vector);
        Base::releaseObjectsInStdVector(_vec_linear_sys);
    }

    LisMPIDiscreteVector* createVector(MPI_Comm comm, int local_n, int global_n)
    {
        LisMPIDiscreteVector* v = new LisMPIDiscreteVector(comm, local_n, global_n);
        _list_vector.push_back(v);
        return v;
    }

    template<class T_SPARSITY_BUILDER>
    IDiscreteLinearEquation* createLinearEquation(MathLib::LisMPILinearEquation &linear_solver, DofEquationIdTable &dofManager)
    {
        _vec_linear_sys.push_back(new LisMPILinearSystem<T_SPARSITY_BUILDER>(*_local_msh, linear_solver, dofManager));
        //return _vec_linear_sys.size()-1;
        return _vec_linear_sys.back();
    };


private:
    MeshLib::IMesh *_local_msh;
    std::vector<LisMPIDiscreteVector*> _list_vector;
    std::vector<AbstractMeshBasedDiscreteLinearEquation*> _vec_linear_sys;

    DISALLOW_COPY_AND_ASSIGN (LisMPIDiscreteSystem);
};

class LisSolver // : public ILinearEquations
{
public:
    LisSolver() 
    {
        _global_dim = 0;
        _local_dim = 0;
        _dynamic = false;
    }
    virtual ~LisSolver();

    static void initialize(int argc, char* argv[]);
    static void finalize();

    void setOption(const MathLib::LIS_option &option)
    {
        _option = option;
    }
    MathLib::LIS_option &getOption() 
    {
        return _option;
    }

    void createDynamic(size_t local_n, size_t global_n);
    MathLib::SparseTableCRS<int>* createCRS(size_t local_n, size_t global_n);
    void setOption(const Base::Options &option);
    void reset();

    size_t getDimension() const { return _global_dim; }
    size_t getLocalDimension() const { return _local_dim; }


    double getA(size_t rowId, size_t colId);
    void setA(size_t rowId, size_t colId, double v);
    void addA(size_t rowId, size_t colId, double v);
    void addA(std::vector<size_t> &vec_row_pos, std::vector<size_t> &vec_col_pos, MathLib::Matrix<double> &sub_matrix, double fkt=1.0);
    void addA(std::vector<size_t> &vec_pos, MathLib::Matrix<double> &sub_matrix, double fkt=1.0);

    double getRHS(size_t rowId);
    double* getRHS();
    void setRHS(size_t rowId, double v);
    void addRHS(std::vector<size_t> &vec_pos, double *sub_vector, double fkt=1.0);
    void addRHS(size_t rowId, double v);

    double* getX();

    void setKnownX(size_t row_id, double x);
    void setKnownX(const std::vector<size_t> &vec_id, const std::vector<double> &vec_x);
    void solve();

    void getRange(int &is, int &ie) const {is = _is; ie = _ie;};
    MathLib::SparseTableCRS<int>* getCRS() {return &_crs;} ;
    void assembleMatrix();

    size_t createVector()
    {
        LIS_VECTOR u;
        int err = lis_vector_duplicate(_A,&u); CHKERR(err);
        _vec_u.push_back(u);
        return _vec_u.size()-1;
    }
    void destroyVector(size_t i)
    {
        LIS_VECTOR &u = _vec_u[i];
        lis_vector_destroy(u);
        _vec_u.erase(_vec_u.begin()+i);
    }
    void setVectorAll(size_t i, double v)
    {
        LIS_VECTOR &u = _vec_u[i];
        lis_vector_set_all(v, u);
    }
    void matvecToRHS(size_t i)
    {
        LIS_VECTOR &u = _vec_u[i];
        lis_matvec(_A,u,_b);
    }
private:
    MathLib::LIS_option _option;
    LIS_MATRIX _A;
    LIS_VECTOR _b;
    LIS_VECTOR _x;
    int _local_dim;
    int _global_dim;
    std::vector<double> _tmp_b;
    std::vector<double> _tmp_x;
    int _is;
    int _ie;
    MathLib::SparseTableCRS<int> _crs;
    std::vector<LIS_VECTOR> _vec_u;
    bool _dynamic;
};
} //end
#endif
