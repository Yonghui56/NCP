/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file MeshElementToFemElementType.h
 *
 * Created on 2012-11-29 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include "IMeshElementToFemElementType.h"

namespace FemLib
{

/**
 * \brief One-to-one mapping from mesh element to FE type
 */
class MeshElementToFemElementType
: public IMeshElementToFemElementType
{
    typedef int MyFeType;
    typedef std::vector<std::vector<MyFeType> > MyTable;
public:
    MeshElementToFemElementType(size_t total_n_ele, size_t max_order)
    {
        _ele2fetype.resize(total_n_ele, std::vector<MyFeType>(max_order));
    }

    virtual ~MeshElementToFemElementType() {};

    /**
     * get FE type
     * @param e         Mesh element
     * @param order     Polynomial order
     * @return FE type
     */
    virtual int getFeType(const MeshLib::IElement& e, const int order) const
    {
        const size_t ele_id = e.getID();
        assert (ele_id < _ele2fetype.size());
        assert ((unsigned)order < _ele2fetype[ele_id].size()+1);

        return _ele2fetype[ele_id][order];
    }

    /**
     * add a mapping between FE type and element shape
     * @param ele_id
     * @param order
     * @param fe_type
     */
    void addFeType(const size_t ele_id, int order, int fe_type)
    {
        assert (ele_id < _ele2fetype.size());
        assert ((unsigned)order < _ele2fetype[ele_id].size()+1);
        _ele2fetype[ele_id][order] = fe_type;
    }

private:
    MyTable _ele2fetype;
};

}
