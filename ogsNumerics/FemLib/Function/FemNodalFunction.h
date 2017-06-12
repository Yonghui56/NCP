/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemNodalFunction.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <algorithm>

#include "BaseLib/CodingTools.h"

#include "MathLib/Vector.h"
#include "MathLib/LinAlg/VectorNorms.h"
#include "MeshLib/Core/IMesh.h"

//#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "DiscreteLib/Core/IDiscreteVector.h"
#include "MathLib/DataType.h"

#include "NumLib/Function/TXFunction.h"
#include "NumLib/Function/ITXDiscreteFunction.h"

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Tools/FemElementObjectContainer.h"


#include "FemLib/Core/PolynomialOrder.h"


namespace FemLib
{

/**
 * \brief Template class for FEM node-based functions (assuming Lagrangian elements)
 *
 * This class represents the following
 * u^h(x) = N(x)*u_i
 *
 * @tparam Tvalue Nodal value type, e.g. double, vector
 */
template<class T_DIS_SYS, typename Tvalue>
class TemplateFEMNodalFunction : public NumLib::ITXDiscreteFunction<Tvalue>
{
public:
    typedef TemplateFEMNodalFunction<T_DIS_SYS, Tvalue> MyClassType;
    typedef T_DIS_SYS MyDiscreteSystem;
    typedef DiscreteLib::IDiscreteVector<Tvalue> MyVector;

    ///
    TemplateFEMNodalFunction()
    : _discrete_system(0), _nodal_values(0), _order(PolynomialOrder::Linear), _feObjects(0)
    {};

//    /// @param dis         Discrete system
//    /// @param order     Polynomial order
//    /// @param v0        initial value
//    TemplateFEMNodalFunction(MyDiscreteSystem &dis, PolynomialOrder::type order, Tvalue v0)
//    {
//        initialize(dis, order);
//        resetNodalValues(v0);
//    }
//
//    /// @param dis         Discrete system
//    /// @param order Polynomial order
//    TemplateFEMNodalFunction(MyDiscreteSystem &dis, PolynomialOrder::type order)
//    {
//        initialize(dis, order);
//    }

    /// @param org source object for copying
    explicit TemplateFEMNodalFunction(const MyClassType &org)
    {
        assign(org);
    }

    ///
    virtual ~TemplateFEMNodalFunction()
    {
        if (_nodal_values != 0) {
            _discrete_system->deleteVector(_nodal_values);
            _nodal_values = 0;
        }
        //BaseLib::releaseObject(_feObjects);
    }

    ///
    TemplateFEMNodalFunction &operator=(const TemplateFEMNodalFunction &org)
    {
        assign(org);
        return *this;
    }

    /// make a clone of this object
    /// @return MathLib::IFunction*
    MyClassType* clone() const
    {
        return new MyClassType(*this);
    };

    /// get this discrete system
    MyDiscreteSystem* getDiscreteSystem() const {return _discrete_system;};

    ///
    size_t getNumberOfNodes() const {return _nodal_values->size();};


    /// evaluate this function at the given point
    virtual void eval(const NumLib::TXPosition x, NumLib::ITXFunction::DataType &v) const
    {
        switch (x.getIdObjectType()) {
        case NumLib::TXPosition::Node:
            {
                v = (*_nodal_values)[x.getId()];
            }
            break;
        default:
            break;
        }
    };

    /// get nodal value
    Tvalue& getValue(size_t node_id)
    {
        return (*_nodal_values)[node_id];
    }

    /// get nodal value
    const Tvalue& getValue(size_t node_id) const
    {
        return (*_nodal_values)[node_id];
    }

    ///
    void setValue(size_t node_id, Tvalue &v)
    {
        (*_nodal_values)[node_id] = v;
    }

    /// get nodal values
    MyVector* getDiscreteData()
    {
        return _nodal_values;
    }

    ///
    const MyVector* getDiscreteData() const
    {
        return _nodal_values;
    }

    /// set nodal values
    void setNodalValues( Tvalue* x, size_t i_start, size_t n )
    {
        for (size_t i=0; i<n; ++i)
            (*_nodal_values)[i+i_start] = x[i];
    }

    /// set nodal values
    void setNodalValues( const MyVector &x )
    {
        *_nodal_values = x;
    }

    /// reset nodal values with the given value
    void resetNodalValues (Tvalue &v)
    {
        *_nodal_values = v;
    }

    /// get Finite element object container
    IFeObjectContainer* getFeObjectContainer() const
    {
        return _feObjects;
    }

    void setFeObjectContainer(IFeObjectContainer *fe)
    {
        _feObjects = fe;
    }

    /// printout internal data for debugging
    void printout() const
    {
        std::cout << "nodal_values = ";
        for (size_t i=_nodal_values->getRangeBegin(); i<_nodal_values->getRangeEnd(); ++i)
            std::cout << (*_nodal_values)[i] << " ";
        std::cout << std::endl;
    }

    /// initialize
    void initialize(MyDiscreteSystem &dis, PolynomialOrder::type order)
    {
        _discrete_system = &dis;
        _order = order;
        _nodal_values = dis.template createVector<Tvalue>(dis.getMesh()->getNumberOfNodes());
        _feObjects = 0;
    }

    template <class T_SYS>
    void initialize(T_SYS &dis, PolynomialOrder::type order, Tvalue v0)
    {
        initialize(dis, order);
        resetNodalValues(v0);
    }

    PolynomialOrder::type getOrder() const {return _order;};

private:
    /// Assign this object from the given object
    void assign(const MyClassType &org)
    {
        initialize(*org._discrete_system, org._order);
        for (size_t i=org._nodal_values->getRangeBegin(); i<org._nodal_values->getRangeEnd(); ++i)
            (*_nodal_values)[i] = (*org._nodal_values)[i];
    }

private:
    MyDiscreteSystem* _discrete_system;
    MyVector* _nodal_values;
    PolynomialOrder::type _order;
    IFeObjectContainer* _feObjects;
};

template <class T_DIS_SYS>
class TemplateFEMNodalFunction<T_DIS_SYS, double> : public NumLib::ITXDiscreteFunction<double>
{
public:
    inline void eval(const NumLib::TXPosition x,  NumLib::ITXFunction::DataType &v) const
    {
        NumLib::ITXFunction::DataType val(1,1);
        val(0,0) = (*_nodal_values)[x.getId()];
        v = val;
    };

public:
    typedef TemplateFEMNodalFunction<T_DIS_SYS, double> MyClassType;
    typedef T_DIS_SYS MyDiscreteSystem;
    typedef DiscreteLib::IDiscreteVector<double> MyVector;

    ///
    TemplateFEMNodalFunction()
    : _discrete_system(0), _nodal_values(0), _order(PolynomialOrder::Linear), _feObjects(0)
    {};

    /// @param org source object for copying
    explicit TemplateFEMNodalFunction(const MyClassType &org)
    {
        assign(org);
    }

    ///
    virtual ~TemplateFEMNodalFunction()
    {
        if (_nodal_values != 0) {
            _discrete_system->deleteVector(_nodal_values);
            _nodal_values = 0;
        }
        //BaseLib::releaseObject(_feObjects);
    }

    ///
    TemplateFEMNodalFunction &operator=(const TemplateFEMNodalFunction &org)
    {
        assign(org);
        return *this;
    }

    /// make a clone of this object
    /// @return MathLib::IFunction*
    MyClassType* clone() const
    {
        return new MyClassType(*this);
    };

    /// get this discrete system
    MyDiscreteSystem* getDiscreteSystem() const {return _discrete_system;};

    ///
    size_t getNumberOfNodes() const {return _nodal_values->size();};



    /// get nodal value
    double getValue(size_t node_id) const
    {
        return (*_nodal_values)[node_id];
    }

    ///
    void setValue(size_t node_id, double v)
    {
        (*_nodal_values)[node_id] = v;
    }

    /// get nodal values
    MyVector* getDiscreteData()
    {
        return _nodal_values;
    }

    ///
    const MyVector* getDiscreteData() const
    {
        return _nodal_values;
    }

    /// set nodal values
    void setNodalValues( double* x, size_t i_start, size_t n )
    {
        for (size_t i=0; i<n; ++i)
            (*_nodal_values)[i+i_start] = x[i];
    }

    /// set nodal values
    void setNodalValues( const MyVector &x )
    {
        *_nodal_values = x;
    }

	/// reset nodal values with the given value
    void resetNodalValues (double v)
    {
        *_nodal_values = v;
    }

    /// get Finite element object container
    IFeObjectContainer* getFeObjectContainer() const
    {
        return _feObjects;
    }

    void setFeObjectContainer(IFeObjectContainer *fe)
    {
        _feObjects = fe;
    }

    /// printout internal data for debugging
    void printout() const
    {
        std::cout << "nodal_values = ";
        for (size_t i=_nodal_values->getRangeBegin(); i<_nodal_values->getRangeEnd(); ++i)
            std::cout << (*_nodal_values)[i] << " ";
        std::cout << std::endl;
    }

    /// initialize
    void initialize(MyDiscreteSystem &dis, PolynomialOrder::type order)
    {
        _discrete_system = &dis;
        _order = order;
        size_t n = dis.getMesh()->getNumberOfNodes();
        _nodal_values = dis.template createVector<double>(n);
        _feObjects = 0;
    }

    template <class T_SYS>
    void initialize(T_SYS &dis, PolynomialOrder::type order, double v0)
    {
        initialize(dis, order);
        resetNodalValues(v0);
    }

    PolynomialOrder::type getOrder() const {return _order;};

private:
    /// Assign this object from the given object
    void assign(const MyClassType &org)
    {
        initialize(*org._discrete_system, org._order);
        for (size_t i=org._nodal_values->getRangeBegin(); i<org._nodal_values->getRangeEnd(); ++i)
            (*_nodal_values)[i] = (*org._nodal_values)[i];
        _feObjects = org._feObjects;
    }

private:
    MyDiscreteSystem* _discrete_system;
    MyVector* _nodal_values;
    PolynomialOrder::type _order;
    IFeObjectContainer* _feObjects;
};

///// evaluate this function at the given point
//template <class T_DIS_SYS>
//inline void TemplateFEMNodalFunction<T_DIS_SYS,double>::eval(const NumLib::TXPosition x,  NumLib::ITXFunction::DataType &v) const
//{
//    NumLib::ITXFunction::DataType val(1,1);
//    val(0,0) = (*_nodal_values)[x.getId()];
//    v = val;
//};

template <class T_DIS_SYS>
struct FemNodalFunctionScalar
{
    typedef TemplateFEMNodalFunction<T_DIS_SYS, double> type;
};

template <class T_DIS_SYS>
struct FemNodalFunctionVector
{
    typedef TemplateFEMNodalFunction<T_DIS_SYS, MathLib::LocalVector> type;
};

template <class T_DIS_SYS>
struct FemNodalFunctionMatrix
{
	typedef TemplateFEMNodalFunction<T_DIS_SYS, MathLib::LocalMatrix> type;
};

//typedef TemplateFEMNodalFunction<double> FemNodalFunctionScalar;
//typedef TemplateFEMNodalFunction<MathLib::LocalVector> FemNodalFunctionVector;

} //end
