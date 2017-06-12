/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemIntegrationPointFunction.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <valarray>
#include <cmath>

#include "MathLib/Vector.h"
#include "MeshLib/Core/IMesh.h"
#include "MathLib/DataType.h"
#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "DiscreteLib/Core/IDiscreteVector.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/Function/ITXDiscreteFunction.h"



namespace FemLib
{

/**
 * \brief Template class for FEM integration point-based functions
 */
template<class T_DIS_SYS, typename Tvalue>
class TemplateFEMIntegrationPointFunction : public NumLib::ITXDiscreteFunction<MathLib::TemplateVectorX<Tvalue> >
{
public:
    typedef MathLib::TemplateVectorX<Tvalue> IntegrationPointVectorType;
    typedef DiscreteLib::IDiscreteVector<IntegrationPointVectorType > MyDiscreteVector;
    typedef T_DIS_SYS MyDiscreteSystem;

    TemplateFEMIntegrationPointFunction()
    : _discrete_system(0), _values(0)
    {};

    virtual ~TemplateFEMIntegrationPointFunction()
    {
        if (_values!=0) {
            _discrete_system->deleteVector(_values);
            _values = 0;
        }
    }

//    ///
//    TemplateFEMIntegrationPointFunction(MyDiscreteSystem* dis, size_t len)
//    {
//        initialize(dis, len);
//    };
//
//    ///
//    TemplateFEMIntegrationPointFunction(MyDiscreteSystem* dis, size_t len, Tvalue &v0)
//    {
//        initialize(dis, len);
//        for (size_t i=0; i<len; i++) {
//            (*_values)[i].resize(1);
//            (*_values)[i][0] = v0;
//        }
//    };

    ///
    explicit TemplateFEMIntegrationPointFunction(const TemplateFEMIntegrationPointFunction &src)
    : _discrete_system(0), _values(0)
    {
        initialize((MyDiscreteSystem*)src._discrete_system, src._values);
        //TODO (*this->_values) = (*src._values);
    };

    /// vector operation: set data
    virtual TemplateFEMIntegrationPointFunction& operator= (const TemplateFEMIntegrationPointFunction &src)
    {
        initialize((MyDiscreteSystem*)src._discrete_system, src._values);
        return *this;
    }

    ///
    TemplateFEMIntegrationPointFunction<T_DIS_SYS, Tvalue>* clone() const
    {
        return new TemplateFEMIntegrationPointFunction<T_DIS_SYS, Tvalue>(*this);
    };

    ///
    const MeshLib::IMesh* getMesh() const
    {
        return this->_discrete_system->getMesh();
    }

    ///
    MyDiscreteSystem* getDiscreteSystem() const {return _discrete_system;};

    ///
    virtual void eval(const NumLib::TXPosition x, NumLib::ITXFunction::DataType &v) const
    {
        switch (x.getIdObjectType()) {
        case NumLib::TXPosition::IntegrationPoint:
            {
                size_t ele_id = x.getId(0);
                size_t gp_id = x.getId(1);
                IntegrationPointVectorType &gp_values = (*_values)[ele_id];
                if (gp_values.size()<gp_id+1) return;
                v = gp_values[gp_id]; 
            }
            break;
        case NumLib::TXPosition::Element:
            {
                // calculate mean value
                size_t ele_id = x.getId();
                IntegrationPointVectorType &gp_values = (*_values)[ele_id];
                if (gp_values.size()==0) return;
                Tvalue val = gp_values[0]; 
                val *= .0;
                for (size_t i=0; i<gp_values.size(); i++)
                    val += gp_values[i];
                val /= gp_values.size();
                v = val;
            }
            break;
        default:
            if (_values->size() == 0) return;
            if ((*_values)[0].size() == 0) return;
            v = (*_values)[0][0];
            break;
        }
    };

    void setIntegrationPointValue( size_t i_e, size_t ip, Tvalue &q )
    {
        assert(ip<(*_values)[i_e].size());
        (*_values)[i_e][ip] = q;
    }

    void setNumberOfIntegationPoints(size_t i_e, size_t n)
    {
        (*_values)[i_e].resize(n);
    }

    const IntegrationPointVectorType& getIntegrationPointValues(size_t i_e) const
    {
        return (*_values)[i_e];
    }

    bool hasIntegrationPointValues(size_t i_e) const
    {
        return (*_values)[i_e].size()>0;
    }

    const MyDiscreteVector* getDiscreteData() const
    {
        return _values;
    }

    MyDiscreteVector* getDiscreteData()
    {
        return _values;
    }

    void printout() const
    {
        std::cout << "integration_pt_values = ";
        for (size_t i=_values->getRangeBegin(); i<_values->getRangeEnd(); ++i) {
            const IntegrationPointVectorType &val1 = (*_values)[i];
            std::cout << "(";
            for (size_t j=0; j<val1.size(); ++j) std::cout << val1[j] << " ";
            std::cout << ") ";
        }
        std::cout << std::endl;
    }

public:
    void initialize(MyDiscreteSystem* dis)
    {
        _discrete_system = dis;
        _values = dis->template createVector<IntegrationPointVectorType>(dis->getMesh()->getNumberOfElements());
    }

    void initialize(MyDiscreteSystem* dis, Tvalue v0)
    {
        initialize(dis);
        for (size_t i=0; i<_values->size(); i++) {
            (*_values)[i].resize(1);
            (*_values)[i][0] = v0;
        }
    }

    void initialize(MyDiscreteSystem* dis, MyDiscreteVector* val)
    {
        initialize(dis);
        (*_values) = (*val);
    }

private:
    MyDiscreteSystem* _discrete_system;
    MyDiscreteVector* _values;
};

template<class T_DIS_SYS>
class TemplateFEMIntegrationPointFunction<T_DIS_SYS, double> : public NumLib::ITXFunction
{
public:
    typedef MathLib::TemplateVectorX<double> IntegrationPointVectorType;
    typedef DiscreteLib::IDiscreteVector<IntegrationPointVectorType > MyDiscreteVector;
    typedef T_DIS_SYS MyDiscreteSystem;

    inline void eval(const NumLib::TXPosition x,  NumLib::ITXFunction::DataType &v) const
    {
        size_t ele_id = x.getId();
        IntegrationPointVectorType &gp_values = (*_values)[ele_id];
        if (gp_values.size()==0) return;
        double val = .0;
        for (size_t i=0; i<gp_values.size(); i++)
            val += gp_values[i];
        val /= gp_values.size();

        v.resize(1,1);
        v(0,0) = val;
        NumLib::ITXFunction::DataType mat(1,1);
    };


    TemplateFEMIntegrationPointFunction()
    : _discrete_system(0), _values(0)
    {};

    virtual ~TemplateFEMIntegrationPointFunction()
    {
        if (_values!=0) {
            _discrete_system->deleteVector(_values);
            _values = 0;
        }
    }

//    ///
//    TemplateFEMIntegrationPointFunction(MyDiscreteSystem* dis, size_t len)
//    {
//        initialize(dis, len);
//    };
//
//    ///
//    TemplateFEMIntegrationPointFunction(MyDiscreteSystem* dis, size_t len, Tvalue &v0)
//    {
//        initialize(dis, len);
//        for (size_t i=0; i<len; i++) {
//            (*_values)[i].resize(1);
//            (*_values)[i][0] = v0;
//        }
//    };

    ///
    explicit TemplateFEMIntegrationPointFunction(const TemplateFEMIntegrationPointFunction &src)
    {
        initialize(src._discrete_system, src._values->size());
        //TODO (*this->_values) = (*src._values);
        _values = 0;
    };

    ///
    TemplateFEMIntegrationPointFunction<T_DIS_SYS, double>* clone() const
    {
        return new TemplateFEMIntegrationPointFunction<T_DIS_SYS, double>(*this);
    };

//    ///
//    const MeshLib::IMesh* getMesh() const
//    {
//        return this->_discrete_system->getMesh();
//    }

    ///
    MyDiscreteSystem* getDiscreteSystem() const {return _discrete_system;};

//    ///
//    virtual void eval(const NumLib::TXPosition x, NumLib::ITXFunction::DataType &v) const
//    {
//        switch (x.getIdObjectType()) {
//        case NumLib::TXPosition::IntegrationPoint:
//            {
//                size_t ele_id = x.getId(0);
//                size_t gp_id = x.getId(1);
//                IntegrationPointVectorType &gp_values = (*_values)[ele_id];
//                if (gp_values.size()==0) return;
//                v = gp_values[gp_id];
//            }
//            break;
//        case NumLib::TXPosition::Element:
//            {
//                // calculate mean value
//                size_t ele_id = x.getId();
//                IntegrationPointVectorType &gp_values = (*_values)[ele_id];
//                if (gp_values.size()==0) return;
//                Tvalue val = gp_values[0];
//                val *= .0;
//                for (size_t i=0; i<gp_values.size(); i++)
//                    val += gp_values[i];
//                val /= gp_values.size();
//                v = val;
//            }
//            break;
//        default:
//            break;
//        }
//    };

    void setIntegrationPointValue( size_t i_e, size_t ip, double &q )
    {
        assert(ip<(*_values)[i_e].size());
        (*_values)[i_e][ip] = q;
    }

    void setNumberOfIntegationPoints(size_t i_e, size_t n)
    {
        (*_values)[i_e].resize(n);
    }

    const IntegrationPointVectorType& getIntegrationPointValues(size_t i_e) const
    {
        return (*_values)[i_e];
    }

    bool hasIntegrationPointValues(size_t i_e) const
    {
        return (*_values)[i_e].size()>0;
    }

    const MyDiscreteVector* getElementValues() const
    {
        return _values;
    }

    void printout() const
    {
        std::cout << "integration_pt_values = ";
        for (size_t i=_values->getRangeBegin(); i<_values->getRangeEnd(); ++i) {
            const IntegrationPointVectorType &val1 = (*_values)[i];
            std::cout << "(";
            for (size_t j=0; j<val1.size(); ++j) std::cout << val1[j] << " ";
            std::cout << ") ";
        }
        std::cout << std::endl;
    }

public:
    void initialize(MyDiscreteSystem* dis)
    {
        _discrete_system = dis;
        _values = dis->template createVector<IntegrationPointVectorType>(dis->getMesh()->getNumberOfElements());
    }

    void initialize(MyDiscreteSystem* dis, double v0)
    {
        initialize(dis);
        for (size_t i=0; i<_values->size(); i++) {
            (*_values)[i].resize(1);
            (*_values)[i][0] = v0;
        }
    }

private:
    MyDiscreteSystem* _discrete_system;
    MyDiscreteVector* _values;
};

template <class T_DIS_SYS>
struct FEMIntegrationPointFunctionScalar
{
    typedef TemplateFEMIntegrationPointFunction<T_DIS_SYS, double> type;
};

template <class T_DIS_SYS>
struct FEMIntegrationPointFunctionVector
{
    typedef TemplateFEMIntegrationPointFunction<T_DIS_SYS, MathLib::LocalVector> type;
};

//typedef TemplateFEMIntegrationPointFunction<LocalVector> FEMIntegrationPointFunctionVector;

} //end
