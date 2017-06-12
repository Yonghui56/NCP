/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IFunction.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "IClonable.h"

namespace NumLib
{

/**
 *
 */
class IFunction : public IClonable
{
public:
    virtual ~IFunction() {};
};

/**
 *
 */
template<typename Tpos, typename Tval>
class TemplateFunction : public IFunction
{
public:
    virtual ~TemplateFunction() {};
    virtual void eval(const Tpos &x, Tval &val) = 0;
};

/**
 *
 */
template<typename Tpos, typename Tval, class T_FUNCTION>
class AbstractDefaultCloneFunction : public TemplateFunction<Tpos, Tval>
{
public:
    virtual ~AbstractDefaultCloneFunction() {};
    virtual IFunction* clone() const
    {
        return new T_FUNCTION();
    }
};

}
