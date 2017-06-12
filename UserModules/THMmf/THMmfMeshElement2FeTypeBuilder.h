/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file THMmfMeshElement2FeTypeBuilder.h
 *
 * Created on 2012-11-29 by Norihiro Watanabe
 */

#pragma once

#include <set>
#include "logog.hpp"

#include "MeshLib/Core/IMesh.h"
#include "FemLib/Tools/MeshElementToFemElementType.h"
#include "MaterialLib/IMedium.h"
#include "MaterialLib/MediumType.h"
#include "THMmfFiniteElementType.h"

namespace THMmf
{

class THMmfMeshElement2FeTypeBuilder
{
public:
    static FemLib::MeshElementToFemElementType* construct(const MeshLib::IMesh &msh, const size_t max_order, std::vector<MaterialLib::IMedium*> list_mmp)
    {
        std::set<size_t> set_matno_ie;
        for (size_t i=0; i<list_mmp.size(); i++) {
            if (list_mmp[i]->getMediumType() == MaterialLib::MediumType::Fracture) {
                set_matno_ie.insert(i);
            }
        }

        const size_t n_ele = msh.getNumberOfElements();
        FemLib::MeshElementToFemElementType* ele2fetype = new FemLib::MeshElementToFemElementType(n_ele, max_order);
        for (size_t i=0; i<n_ele; ++i) {
            const MeshLib::IElement* e = msh.getElement(i);
            if (set_matno_ie.count(e->getGroupID()) > 0) {
                switch (e->getShapeType()) {
                    case MeshLib::ElementShape::QUAD:
                        ele2fetype->addFeType(i, 1, THMmfFiniteElementType::IE_QUAD4);
                        break;
                    case MeshLib::ElementShape::TRIANGLE:
                        ele2fetype->addFeType(i, 1, THMmfFiniteElementType::IE_TRI3);
                        break;
                    default:
                        ERR("***Error: Unsupported element type (%d) for interface FE. ", e->getShapeType() );
                        break;
                }
            } else {
                switch (e->getShapeType()) {
                    case MeshLib::ElementShape::LINE:
                        ele2fetype->addFeType(i, 1, THMmfFiniteElementType::LINE2);
                        ele2fetype->addFeType(i, 2, THMmfFiniteElementType::LINE3);
                        break;
                    case MeshLib::ElementShape::QUAD:
                        ele2fetype->addFeType(i, 1, THMmfFiniteElementType::QUAD4);
                        ele2fetype->addFeType(i, 2, THMmfFiniteElementType::QUAD9);
                        break;
                    case MeshLib::ElementShape::TRIANGLE:
                        ele2fetype->addFeType(i, 1, THMmfFiniteElementType::TRI3);
                        ele2fetype->addFeType(i, 2, THMmfFiniteElementType::TRI6);
                        break;
                    case MeshLib::ElementShape::TETRAHEDRON:
                        ele2fetype->addFeType(i, 1, THMmfFiniteElementType::TET4);
                        ele2fetype->addFeType(i, 2, THMmfFiniteElementType::TET10);
                        break;
                    default:
                        ERR("***Error: Unsupported element type (%d) in THMmf. ", e->getShapeType() );
                        break;
                }
            }
        }

        return ele2fetype;
    }
};

}
