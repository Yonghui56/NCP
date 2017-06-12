/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file MeshBasedProblem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

namespace MeshLib
{
class IMesh;
}


namespace SolutionLib
{

/**
 * \brief Mesh based discrete IVBV problems
 */
class MeshBasedProblem
{
public:
    ///
    explicit MeshBasedProblem(MeshLib::IMesh* msh) : _msh(msh)
    {};

    ///
    virtual ~MeshBasedProblem() {}

    /// set a mesh
    void setMesh(MeshLib::IMesh* msh) {_msh = msh;};

    /// get the mesh
    MeshLib::IMesh* getMesh() const {return _msh;};

private:
    MeshLib::IMesh* _msh;
};

}
