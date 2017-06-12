/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DiscreteEnums.h
 *
 * Created on 2012-08-17 by Norihiro Watanabe
 */

#pragma once


namespace DiscreteLib
{

/**
 * \brief Discrete System Type
 */
struct DiscreteSystemType
{
    enum type
    {
        Serial,
        SerialShared,
        SerialDistributed,
        OpenMP,
        MPI
    };
};

}
