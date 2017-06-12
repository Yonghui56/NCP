
#include "LisMPIDiscreteVector.h"

#include <iostream>

namespace DiscreteLib
{

LisMPIDiscreteVector::LisMPIDiscreteVector(MPI_Comm comm, int n_local, int n_global)
{
    _exist_any_ghost_global = false;
    initialize(comm, n_local, n_global);
};

LisMPIDiscreteVector::LisMPIDiscreteVector(MPI_Comm comm, int n_local, int n_global, std::vector<int> &ghost_id) : _ghost_id(ghost_id)
{
    _exist_any_ghost_global = true;
    initialize(comm, n_local, n_global, ghost_id.size());
};

/// destructor
LisMPIDiscreteVector::~LisMPIDiscreteVector() 
{
    lis_vector_destroy(_v);
};

/// check if this vector contain an element with the given global id
bool LisMPIDiscreteVector::has(int global_idx) const
{
    bool found = false;
    if (!isGhost(global_idx)) {
        //std::cout << _local_id << ": " << global_idx << " is not ghost" << std::endl;
        found = (_i_start<=global_idx && global_idx<_i_end);
    } else {
        std::cout << _local_id << ": " << global_idx << " is ghost" << std::endl;
        found = (_ghost_id.end()!=std::find(_ghost_id.begin(), _ghost_id.end(), global_idx));
        //std::cout << _local_id << ": ghost list = ";
        //for (size_t i=0; i<_ghost_id.size(); i++) {
        //    std::cout << _ghost_id[i] << " ";
        //}
        std::cout << std::endl;
    }
    if (!found)
        std::cout << _local_id << ": ***error-> " << global_idx << " is not found" << std::endl;
    return found;
}

void LisMPIDiscreteVector::finishUpdate()
{
    int err = 0;
    // push
    //std::cout << _local_id << ": set values" << std::endl;
    //err = lis_vector_set_values2(LIS_INS_VALUE, _i_start, _local_n, &_data[0], _v);
    err = lis_vector_set_values(LIS_INS_VALUE, _local_n, &_list_index[0], &_data[0], _v);
    CHKERR(err);
    // pull
    if (_exist_any_ghost_global) {
        //std::cout << _local_id << ": get ghost values" << std::endl;
        err = lis_vector_gather(_v, &_temp_all_x[0]);
        CHKERR(err);
        for (size_t i=0; i<_ghost_id.size(); i++) {
            const int ghost_id = _ghost_id[i];
            //std::cout << _local_id << ": ghost_id=" << ghost_id << ", local id=" << access(ghost_id) << std::endl;
            _data[access(ghost_id)] = _temp_all_x[ghost_id];
            //err = lis_vector_get_value(_v, ghost_id, &_data[access(ghost_id)]);
        }
        CHKERR(err);
    }
}

void LisMPIDiscreteVector::getGlobal(double *v)
{
    int err = lis_vector_gather(_v, v);
    CHKERR(err);
}

void LisMPIDiscreteVector::print()
{
    int err = lis_vector_print(_v);
    CHKERR(err);
}

/// initialize this object
void LisMPIDiscreteVector::initialize(MPI_Comm comm, int n_local, int n_global, size_t n_ghost)
{
    MPI_Comm_rank(comm, &_local_id);
    int err = lis_vector_create(comm, &_v);
    err = lis_vector_set_size(_v, n_local, n_global);
    err = lis_vector_get_size(_v, &_local_n, &_global_n);
    err = lis_vector_get_range(_v, &_i_start, &_i_end);
    CHKERR(err);
    DiscreteVector::resize(_local_n+n_ghost);
    _temp_all_x.resize(_global_n); //allocate temporary buffer
    _list_index.resize(_local_n);
    for (int i=0; i<_local_n; i++)
        _list_index[i] = _i_start + i;

    std::cout << _local_id << ": start=" << _i_start << ", end=" << _i_end << std::endl;
}

double LisMPIDiscreteVector::getNorm1()
{
    double val = .0;
    int err = lis_vector_nrm1(_v, &val);
    CHKERR(err);
    return val;
}

double LisMPIDiscreteVector::getNorm2()
{
    double val = .0;
    int err = lis_vector_nrm2(_v, &val);
    CHKERR(err);
    return val;
}

double LisMPIDiscreteVector::getNormInf()
{
    double val = .0;
    int err = lis_vector_nrmi(_v, &val);
    CHKERR(err);
    return val;
}

} // end
