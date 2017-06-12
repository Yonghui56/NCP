/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file RandomEquationIdStorage.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <cassert>
#include <iostream>

#include "BaseLib/CodingTools.h"
#include "BaseLib/BidirectionalMap.h"
#include "IEquationIdStorage.h"

namespace DiscreteLib
{

class RandomEquationIdStorage : public IEquationIdStorage
{
public:
    RandomEquationIdStorage(const std::vector<size_t> &sorted_list_pt_id) : _list_pt_id(sorted_list_pt_id)
    {
        for (size_t i=0; i<sorted_list_pt_id.size(); i++)
            set(sorted_list_pt_id[i], 0);
    }
    RandomEquationIdStorage(size_t pt_id_start, size_t n)
    {
        _list_pt_id.resize(n);
        for (size_t i=0; i<n; i++) {
            _list_pt_id[i] = pt_id_start + i;
            set(_list_pt_id[i], 0);
        }
    }
    virtual ~RandomEquationIdStorage() {};

    bool isSequential() const {return false;};

    bool hasKey(size_t pt_id) const
    {
        return _list_pt_id.end() !=  std::find(_list_pt_id.begin(), _list_pt_id.end(), pt_id);
    }

    bool hasValue(size_t eqs_id) const
    {
        return (_map_pt2eqs.countInB(eqs_id)>0);
    }


    void key_range(size_t &i_start, size_t &i_end) const
    {
        i_start = _list_pt_id[0];
        i_end = _list_pt_id.back()+1;
    }

    void activate(size_t pt_id, bool b)
    {
        if (b) {
            _deactive.erase(pt_id);
        } else {
            _deactive.insert(pt_id);
            set(pt_id, -1);
        }
    }

    bool isActive(size_t pt_id) const { return _deactive.count(pt_id)==0;};

    void set(size_t pt_id, long eqs_id)
    {
        _map_pt2eqs.insert(pt_id, eqs_id);
    }

    size_t setAll(size_t dof_start, size_t delta_per_pt)
    {
        size_t last = 0;
        if (_deactive.size()==0) {
            for (size_t i=0; i<_list_pt_id.size(); i++) {
                last = dof_start + i*delta_per_pt;
                set(_list_pt_id[i], last);
            }
        } else {
            size_t counter = 0;
            for (size_t i=0; i<_list_pt_id.size(); i++) {
                const size_t pt_id = _list_pt_id[i];
                if (_deactive.count(pt_id)>0) continue;
                last = dof_start + counter*delta_per_pt;
                set(pt_id, last);
                counter++;
            }
        }
        return last+1;
    }

    size_t size() const {return _map_pt2eqs.size();};

    size_t address(size_t pt_id) const
    {
        if (_map_pt2eqs.countInA(pt_id)==0) return -1;
        return _map_pt2eqs.mapAtoB(pt_id);
    }

    size_t key(size_t address_id) const
    {
        if (_map_pt2eqs.countInB(address_id)==0) return -1;
        return _map_pt2eqs.mapBtoA(address_id);
    }

private:
    BaseLib::BidirectionalMap<size_t, size_t> _map_pt2eqs;
    std::vector<size_t> _list_pt_id;
    std::set<size_t> _deactive;
};

} //end

