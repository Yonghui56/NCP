/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DecomposedDomain.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "BaseLib/CodingTools.h"

#include "DdcEnums.h"

namespace DiscreteLib
{

class SubDomain;

/**
 * \brief Decomposed domain class
 *
 * This class owns
 * - decomposition type
 * - sub-domains
 */
class DecomposedDomain
{
public:
    /**
     * \param ddc_type decomposition type
     */
    explicit DecomposedDomain(DecompositionType::type ddc_type)
    : _id(0), _ddc_type(ddc_type), _n_discrete_pt(0), _n_nodes(0), _n_eles(0)
    { };

    ///
    virtual ~DecomposedDomain()
    {
        BaseLib::releaseObjectsInStdVector(_list_dom);
    }

    /// return this decomposition type
    DecompositionType::type getDecompositionType() const {return _ddc_type;};

    /// return the number of total decomposed objects, e.g. the number of all nodes
    size_t getTotalNumberOfDecomposedObjects() const {return _n_discrete_pt;};

    /// return the number of sub domains
    size_t getNumberOfSubDomains() const {return _list_dom.size();};

    /// add a new sub domain
    size_t addSubDomain(SubDomain* sub);

    /// return a sub domain
    SubDomain* getSubDomain(size_t i) {return _list_dom[i];};

    /// find a sub domain in which the given object exist and return its domain id
    size_t findSubDomainID(size_t global_obj_id) const;

    ///
    size_t getNumberOfNodes() const {return _n_nodes;};

    ///
    size_t getNumberOfElements() const {return _n_eles;};
    
    ///
    size_t getID() const {return _id;};

private:
    DISALLOW_COPY_AND_ASSIGN(DecomposedDomain);

private:
    size_t _id;
    DecompositionType::type _ddc_type;
    std::vector<size_t> _list_dom_start_obj;
    std::vector<SubDomain*> _list_dom;
    size_t _n_discrete_pt;
    size_t _n_nodes;
    size_t _n_eles;
};

} //end

