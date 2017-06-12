/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file MeshElementGroupAndShapeToFemElementType.h
 *
 * Created on 2012-11-27 by Norihiro Watanabe
 */

#pragma once

#include <map>
#include "IMeshElementToFemElementType.h"

namespace FemLib
{

/**
 * \brief Converter from element group ID and shape to FE type
 */
class MeshElementGroupAndShapeToFemElementType
: public IMeshElementToFemElementType
{
    typedef int OrderType;
    typedef size_t GroupIdType;
    typedef std::pair<MeshLib::ElementShape::type, OrderType> MyKeySub;
    typedef std::pair<GroupIdType, MyKeySub> MyKey;
    typedef std::map<MyKey, int> MyTable;
public:
    virtual ~MeshElementGroupAndShapeToFemElementType() {};

    /**
     *get FE type
     * @param e
     * @param order
     * @return
     */
    virtual int getFeType(const MeshLib::IElement& e, const int order) const
    {
        MyKey key(e.getGroupID(), MyKeySub(e.getShapeType(), order));
        MyTable::const_iterator itr = _map_ele2fetype.find(key);
        if (itr != _map_ele2fetype.end())
            return itr->second;
        else
            return -1;
    }

    /**
     * register FE type
     * @param group_id
     * @param shape_type
     * @param order
     * @param fe_type
     */
    void addFeType(size_t group_id, MeshLib::ElementShape::type shape_type, int order, int fe_type)
    {
        MyKey key(group_id, MyKeySub(shape_type, order));
        _map_ele2fetype[key] = fe_type;
    }

private:
    MyTable _map_ele2fetype;
};

}
