/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TXFunctionElementGroupID.h
 *
 * Created on 2012-11-15 by Norihiro Watanabe
 */
 
 #pragma once

#include "NumLib/Function/TXPosition.h"
#include "NumLib/Function/ITXFunction.h"
#include "MeshLib/Core/IMesh.h"

namespace ogs6
{

class TXFunctionElementGroupID : public NumLib::ITXFunction
{
public:
    explicit TXFunctionElementGroupID(const MeshLib::IMesh* msh)
    : _msh(msh)
    {
        ITXFunction::isConst(true);
        ITXFunction::isTemporallyConst(true);
        ITXFunction::isSpatiallyConst(false);
    };

    virtual ~TXFunctionElementGroupID() {};

    virtual void eval(const NumLib::TXPosition x, DataType &v) const
    {
        v.resize(1,1);
        switch (x.getIdObjectType()) {
        case NumLib::TXPosition::Node:
        case NumLib::TXPosition::Element:
            {
                v(0,0) = _msh->getElement(x.getId())->getGroupID();
            }
            break;
        default:
            break;
        }
    }
    
    virtual TXFunctionElementGroupID* clone() const
    {
        return new TXFunctionElementGroupID(_msh);
    }


private:
    const MeshLib::IMesh* _msh;
};

} //end
