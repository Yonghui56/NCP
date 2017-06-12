/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FdmFunction.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "NumLib/Function/TXFunction.h"
#include "MathLib/LinAlg/VectorNorms.h"
#include "DiscreteLib/Serial/DiscreteSystem.h"

namespace FdmLib
{

/**
 * \brief
 *
 */
template<typename Tvalue>
class TemplateFDMFunction : public NumLib::ITXFunction
{
public:
    /// @param msh         Mesh
    TemplateFDMFunction(DiscreteLib::DiscreteSystem &dis, MeshLib::IMesh &msh, Tvalue v0)
    {
        initialize(dis, msh);
        resetNodalValues(v0);
    }

    /// @param msh Mesh
    TemplateFDMFunction(DiscreteLib::DiscreteSystem &dis, MeshLib::IMesh &msh)
    {
        initialize(dis, msh);
    }

    /// @param org source object for copying
    TemplateFDMFunction(const TemplateFDMFunction<Tvalue> &org)
    {
        assign(org);
    }

    ///
    virtual ~TemplateFDMFunction()
    {
    }

    ///
    TemplateFDMFunction &operator=(const TemplateFDMFunction &org)
    {
        assign(org);
        return *this;
    }

    /// make a clone of this object
    /// @return MathLib::IFunction*
    TemplateFDMFunction<Tvalue>* clone() const
    {
        TemplateFDMFunction<Tvalue> *obj = new TemplateFDMFunction<Tvalue>(*this);
        return obj;
    };

    /// get the spatial dimension
    size_t getDimension() const { return _msh->getDimension(); }

    /// get the mesh
    const MeshLib::IMesh* getMesh() const { return _msh; }

    /// get the number of nodes
    size_t getNumberOfNodes() const { return _msh->getNumberOfNodes(); }


    /// evaluate this function at the given point
    void eval(const NumLib::TXPosition &/*pt*/, Tvalue &v)
    {
        throw "eval() is not implemented yet.";
        v = (*_nodal_values)[0];
    };

    /// get nodal value
    Tvalue& getValue(size_t node_id)
    {
        return (*_nodal_values)[node_id];
    }

    /// get nodal values
    DiscreteLib::IDiscreteVector<Tvalue>* getNodalValues()
    {
        return _nodal_values;
    }

    const DiscreteLib::IDiscreteVector<Tvalue>* getNodalValues() const
    {
        return _nodal_values;
    }

    /// set nodal values
    void setNodalValues( Tvalue* x, size_t i_start, size_t n )
    {
        for (size_t i=0; i<n; ++i)
            (*_nodal_values)[i+i_start] = x[i];
        //std::copy(x, x+getNumberOfNodes(), _nodal_values->begin());
    }

    /// set nodal values
    void setNodalValues( const DiscreteLib::IDiscreteVector<Tvalue> &x )
    {
        *_nodal_values = x;
    }

    /// reset nodal values with the given value
    void resetNodalValues (Tvalue &v)
    {
        std::fill(_nodal_values->begin(), _nodal_values->end(), v);
    }

//    double norm_diff(const TemplateFDMFunction<Tvalue> &ref) const
//    {
//        const size_t n = _nodal_values->size();
//        if (n!=ref._nodal_values->size()) {
//            std::cout << "***Warning in TemplateFEMIntegrationPointFunction::norm_diff(): size of two vectors is not same." << std::endl;
//            return .0;
//        }
//
//        const DiscreteLib::DiscreteVector<double>* vec_prev = ref.getNodalValues();
//        const DiscreteLib::DiscreteVector<double>* vec_cur = this->getNodalValues();
//        DiscreteLib::DiscreteVector<double> vec_diff(vec_prev->size());
//        vec_diff = *vec_cur;
//        vec_diff -= *vec_prev;
//        return MathLib::norm_max(vec_diff, vec_diff.size());
//    }

    void printout() const
    {
        std::cout << "nodal_values = ";
        for (size_t i=_nodal_values->getRangeBegin(); i<_nodal_values->getRangeEnd(); ++i)
            std::cout << (*_nodal_values)[i] << " ";
        std::cout << std::endl;
    }
private:
    /// initialize
    void initialize(DiscreteLib::DiscreteSystem &dis, MeshLib::IMesh &msh)
    {
        _discrete_system = &dis;
        _msh = &msh;
        size_t nnodes = msh.getNumberOfNodes();
        _nodal_values = dis.createVector<Tvalue>(nnodes);
    }

    /// Assign this object from the given object
    void assign(const TemplateFDMFunction<Tvalue> &org)
    {
        initialize(*org._discrete_system, *org._msh);
        for (size_t i=org._nodal_values->getRangeBegin(); i<org._nodal_values->getRangeEnd(); ++i)
            (*_nodal_values)[i] = (*org._nodal_values)[i];
        //std::copy(org._nodal_values->begin(), org._nodal_values->end(), _nodal_values->begin());
    }

private:
    DiscreteLib::DiscreteSystem* _discrete_system;
    DiscreteLib::IDiscreteVector<Tvalue>* _nodal_values;
    MeshLib::IMesh* _msh;
};


typedef TemplateFDMFunction<double> FdmFunctionScalar;


class FdmCellVectorFunction : public NumLib::ITXFunction
{
public:
    explicit FdmCellVectorFunction(size_t n)
    {
        _vec.resize(n);
    };
    virtual ~FdmCellVectorFunction() {};

    virtual FdmCellVectorFunction* clone() const
    {
        FdmCellVectorFunction *obj = new FdmCellVectorFunction(_vec.size());
        return obj;
    };

    virtual void eval(const NumLib::TXPosition /*x*/, NumLib::ITXFunction::DataType &val) const
    {
        val = _vec[0];
    }

    void setValue(size_t i, MathLib::LocalVector &v)
    {
        _vec[i] = v;
    }

    void printout() const
    {
        std::cout << "cell_values = ";
        for (size_t i=0; i<_vec.size(); ++i) {
            const MathLib::LocalVector &val1 = _vec[i];
            std::cout << "(";
            for (int j=0; j<val1.size(); ++j) std::cout << val1[j] << " ";
            std::cout << ") ";
        }
        std::cout << std::endl;
    }

private:
    std::vector<MathLib::LocalVector> _vec;
};

} //end
