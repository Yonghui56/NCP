
#pragma once

#include "ILinearEquation.h"
#include "LisInterface.h"

namespace MathLib
{

class LisMPILinearEquation : public MathLib::ILinearEquation
{
public:
    LisMPILinearEquation() 
    {
        _global_dim = 0;
        _local_dim = 0;
        _dynamic = false;
        _created = false;
        _A = 0;
        _b = 0;
        _x = 0;
        _is = 0;
        _ie = 0;
    }
    virtual ~LisMPILinearEquation();

    void setOption(const MathLib::LIS_option &option)
    {
        _option = option;
    }
    MathLib::LIS_option &getOption() 
    {
        return _option;
    }

    void create( size_t local_n, MathLib::RowMajorSparsity* sparse )
    {
        _created = true;
    }

    bool isCreated() const { return _created; };

    void createDynamic(size_t local_n, size_t global_n);

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
    void printout(std::ostream &os=std::cout) const
    {

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
    bool _created;
};
} //end
