/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TemplateFeBase.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MeshLib/Core/IElement.h"
#include "MeshLib/Core/IMesh.h"

#include "IFemElement.h"

namespace FemLib
{

/**
 * \brief Base class for implementation of finite element classes
 *
 * This class owns the following data
 * - Mesh element object
 * - Mesh object
 *
 * \tparam T_FETYPE         Finite element type
 * \tparam N_VARIABLES      The number of variables
 */
template < size_t N_VARIABLES >
class TemplateFeBase : public IFiniteElement
{
public:
    /**
     *
     * @param msh
     */
    explicit TemplateFeBase(MeshLib::IMesh* msh) : _fe_type(-1), _msh(msh), _ele(NULL) {};

    ///
    virtual ~TemplateFeBase() {};

    ///
    void setMesh(MeshLib::IMesh* msh) {_msh = msh;};

    ///
    const MeshLib::IMesh* getMesh() const {return _msh;};

    /// return mesh element
    virtual MeshLib::IElement* getElement() const {return _ele;};

    /// return the number of variables
    virtual size_t getNumberOfVariables() const { return N_VARIABLES; };

    /// return finite element type
    virtual int getFeType() const { return TemplateFeBase<N_VARIABLES>::_fe_type; };

    /// set finite element type ID
    void setFeType(int n) {_fe_type = n;};

protected:
    ///
    void setElement(MeshLib::IElement* e) {_ele = e;};

private:
    int _fe_type;
    MeshLib::IMesh* _msh;
    MeshLib::IElement* _ele;
};

}

