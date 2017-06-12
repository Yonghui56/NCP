/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file BidirectionalMap.h
 *
 * Created on 2012-02-24 by Norihiro Watanabe
 */

#pragma once

#include <cstddef>
#include <map>

namespace BaseLib
{

/**
 * \brief A bidrectional map class between unique keys A and B
 *
 * \tparam T1   Value type of the first key
 * \tparam T2   Value type of the second key
 */
template <typename T1, typename T2>
class BidirectionalMap
{
public:
    typedef std::map<T1, T2> MapA;
    typedef std::map<T2, T1> MapB;

    /**
     *
     */
    BidirectionalMap() {};

    /**
     * Copy constructor
     *
     * \param src   source object
     */
    explicit BidirectionalMap(const BidirectionalMap &src)
    {
        _map1 = src._map1;
        _map2 = src._map2;
    }

    /**
     * Copy 
     *
     * \param src   source object
     */
    BidirectionalMap &operator=(const BidirectionalMap &src)
    {
        _map1 = src._map1;
        _map2 = src._map2;
        return *this;
    }

    /**
     * insert a set of keys
     *
     * \param v1    Key1
     * \param v2    Key2
     */
    void insert(const T1 &v1, const T2 &v2)
    {
        _map1[v1] = v2;
        _map2[v2] = v1;
    }

    /**
     * return the number of elements found with key1
     */
    size_t countInA(const T1 &v) const { return _map1.count(v);};

    /**
     * return the number of elements found with key2
     */
    size_t countInB(const T2 &v) const { return _map2.count(v);};

    /**
     * return size of entries
     */
    size_t size() const { return _map1.size(); };


    /**
     * clear entries
     */
    void clear() { _map1.clear(); _map2.clear(); };

    /**
     * return key2 corresponding to the given key1
     */
    const T2& mapAtoB(const T1 &v) const 
    {
        typename MapA::const_iterator itr = _map1.find(v);
        //if (itr!=_map1.end()) 
        return itr->second;
    };

    /**
     * return key1 corresponding to the given key2
     */
    const T1& mapBtoA(const T2 &v) const 
    {
        typename MapB::const_iterator itr = _map2.find(v);
        //if (itr!=_map2.end()) 
        return itr->second;
    };

private:
    MapA _map1;
    MapB _map2;
};

}
