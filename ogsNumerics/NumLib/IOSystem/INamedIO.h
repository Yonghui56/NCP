/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file INamedIO.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include <vector>

namespace NumLib
{

/**
 * \brief Interface of in/out systems with named keys
 */
class INamedIO
{
public:
    ///
    virtual ~INamedIO() {};

    ///
    virtual void setInputParameterName(size_t i, const std::string& key) = 0;

    ///
    virtual void setOutputParameterName(size_t i, const std::string& key) = 0;

    ///
    virtual bool isValid() const = 0;

    ///
    virtual bool hasInputParameter(const std::string& key) const = 0;

    ///
    virtual bool hasOutputParameter(const std::string& key) const = 0;

    /// get the name of parameter with the given name
    virtual const std::string& getInputParameterName(size_t i) const = 0;

    /// get the name of parameter with the given name
    virtual const std::string& getOutputParameterName(size_t i) const = 0;

    /// get the index of parameter with the given name
    virtual int getInputParameterID(const std::string&  key) const = 0;

    /// get the index of parameter with the given name
    virtual int getOutputParameterID(const std::string&  key) const = 0;

    /// get the number of output parameters
    virtual size_t getNumberOfInputParameterNames() const = 0;

    /// get the number of output parameters
    virtual size_t getNumberOfOutputParameterNames() const = 0;
};

} //end
