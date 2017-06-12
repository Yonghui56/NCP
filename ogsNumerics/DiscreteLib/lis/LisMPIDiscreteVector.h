
#pragma once

#include <vector>
#include <algorithm>
#include <cassert>

#include "mpi.h"
#include "lis.h"

#include "DiscreteLib/Core/DiscreteVector.h"


namespace DiscreteLib
{

/**
 * \brief Discrete vector based on LIS and MPI
 *
 * This class utilize DiscreteVector and LIS functionality to manipulate distributed vectors.
 * Local data are stored in DiscreteVector and are updated to Lis when finishUpdate() is called.
 * 
 */
class LisMPIDiscreteVector : public DiscreteVector<double>
{
public:
    /// Constructor without ghost elements
    /// @param comm        MPI communicator
    /// @param rank        rank
    /// @param n_local     size of the local vector
    /// @param n_global    size of the global vector 
    LisMPIDiscreteVector(MPI_Comm comm, int n_local, int n_global);

    /// Constructor with ghost elements
    /// @param comm        MPI communicator
    /// @param rank        rank
    /// @param n_local     size of the local vector
    /// @param n_global    size of the global vector 
    /// @param ghost_id    a list of ghost id
    LisMPIDiscreteVector(MPI_Comm comm, int n_local, int n_global, std::vector<int> &ghost_id);

    /// destructor
    virtual ~LisMPIDiscreteVector();

    /// get the global vector size
    size_t getGlobalSize() const {return _global_n;};

    /// get the range begin
    size_t getRangeBegin() const {return _i_start;};
    /// get the range end
    size_t getRangeEnd() const {return _i_end;};

    /// get the number of ghost elements in this vector
    size_t getNumberOfGhostElements() const {return _ghost_id.size();};
    /// get the ghost element id (global) in this vector
    size_t getGhostElementId(size_t i) const {return _ghost_id[i];};
    /// check if the give id can be a ghost element in this vector. This function doesn't check if actually this vector contain the element.
    bool isGhost(int global_idx) const {return !(_i_start<=global_idx && global_idx<_i_end);};

    /// check if this vector contain an element with the given global id
    bool has(int global_idx) const;

    /// access with local index
    double& operator[] (size_t i_local) { return _data[i_local]; }
    /// access with local index
    const double& operator[] (size_t i_local) const { return _data[i_local]; }

    /// access with global index
    double& global (int i) { return _data[access(i)]; }
    /// access with global index
    const double& global (int i) const { return _data[access(i)]; }

    /// this function should be called whenever the vector is updated
    void finishUpdate();

    /// get global vector
    void getGlobal(double *v);

    /// use Lis print function
    void print();

    /// compute norm1
    double getNorm1();
    /// compute norm2
    double getNorm2();
    /// compute infinite norm
    double getNormInf();

private:
    LIS_VECTOR _v;
    std::vector<double> _temp_all_x;
    int _global_n;
    int _local_n;
    int _i_start;
    int _i_end;
    int _local_id;
    std::vector<int> _list_index;
    /// a list of global id of ghost elements
    std::vector<int> _ghost_id;
    bool _exist_any_ghost_global;

    /// initialize this object
    void initialize(MPI_Comm comm, int n_local, int n_global, size_t n_ghost=0);

    /// get local id
    inline int access(int global_idx) const 
    {
        assert (has(global_idx));
        if (!isGhost(global_idx)) {
            return global_idx - _i_start;
        } else {
            return _i_end + (std::find(_ghost_id.begin(), _ghost_id.end(), global_idx) - _ghost_id.begin()) - _i_start;
        }
    };
};

} //end

