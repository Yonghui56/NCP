/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IStencil.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

namespace FdmLib
{

class IStencil
{
public:
    virtual ~IStencil() {};

    virtual size_t getCentralNodeID() const = 0;

    virtual const std::vector<size_t>& getSurroundingNodes() const = 0;
};


class Stencil5 : public IStencil
{
public:
    Stencil5() : _central_nodeid(0) {};
    virtual ~Stencil5() {};

    virtual void setCentralNodeID(size_t i) { _central_nodeid = i; };
    virtual size_t getCentralNodeID() const { return _central_nodeid; };

    virtual void addSurroundingNode(size_t i) { _surrounding_nodes.push_back(i);};
    virtual const std::vector<size_t>& getSurroundingNodes() const { return _surrounding_nodes;};

private:
    size_t _central_nodeid;
    std::vector<size_t> _surrounding_nodes;
};

}
