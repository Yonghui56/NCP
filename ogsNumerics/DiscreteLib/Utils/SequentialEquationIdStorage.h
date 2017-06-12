/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SequentialEquationIdStorage.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <set>
#include <algorithm>
#include <cassert>
#include <iostream>

#include "BaseLib/CodingTools.h"
#include "IEquationIdStorage.h"


namespace DiscreteLib
{

class SequentialEquationIdStorage : public IEquationIdStorage
{
public:
    SequentialEquationIdStorage(size_t pt_id_start, size_t n) : _pt_id_start(pt_id_start), _n(n)
    {
        _dof_begin = 0;
        _delta_per_pt = 1;
    }
    virtual ~SequentialEquationIdStorage() {};

    bool isSequential() const {return true;};

    bool hasKey(size_t pt_id) const
    {
        return (_pt_id_start<=pt_id && pt_id<_pt_id_start+_n);
    }

    bool hasValue(size_t eqs_id) const
    {
        bool in_range = (address(_pt_id_start) <= eqs_id && eqs_id <= address(_pt_id_start+_n-1));
        if (!in_range) return false;
        if ((eqs_id - _dof_begin)%_delta_per_pt!=0) return false;
        return true;
    }

    void key_range(size_t &i_start, size_t &i_end) const
    {
        i_start = _pt_id_start;
        i_end = i_start + _n;
    }

    void activate(size_t pt_id, bool b)
    {
        if (b) {
            _deactive.erase(pt_id);
        } else {
            _deactive.insert(pt_id);
        }
    }
    
    bool isActive(size_t pt_id) const { return _deactive.count(pt_id)==0;};

    void set(size_t /*pt_id*/, long /*eqs_id*/)
    {
        //invalid
        std::cout << "***Error: called invalid functions. SequentiallyMappedAddress::set()." << std::endl;
    }

    size_t setAll(size_t address_start, size_t dn_pt)
    {
        _dof_begin = address_start;
        _delta_per_pt = dn_pt;
        return _dof_begin + (_n-_deactive.size())*_delta_per_pt;
    }

    size_t size() const {return _n;};

    size_t address(size_t pt_id) const
    {
        assert(_pt_id_start<=pt_id && pt_id<_pt_id_start+_n);

        if (_deactive.count(pt_id)>0) return -1;

        if (_pt_id_start<=pt_id && pt_id<_pt_id_start+_n) {
            size_t loc = 0;
            if (_deactive.size()==0) {
                loc = _dof_begin + (pt_id-_pt_id_start)*_delta_per_pt;
            } else {
                size_t inactive_cnt = 0;
                for (size_t i=_pt_id_start; i<pt_id; i++) {
                    if (_deactive.count(i)>0) inactive_cnt++;
                }
                loc = _dof_begin + (pt_id-_pt_id_start-inactive_cnt)*_delta_per_pt;
            }
            return loc;
        } else {
            return BaseLib::index_npos;
        }
    }

    size_t key(size_t /*address_id*/) const
    {
        //TODO
        return BaseLib::index_npos;
    }

private:
    size_t _pt_id_start;
    size_t _n;
    size_t _dof_begin;
    size_t _delta_per_pt;
    std::set<size_t> _deactive;
};

} //end

