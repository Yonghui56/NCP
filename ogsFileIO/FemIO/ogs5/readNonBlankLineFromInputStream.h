/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file readNonBlankLineFromInputStream.h
 *
 * Created on 2011-04-19 by Thomas Fisher
 */


#ifndef READNONBLANKLINEFROMINPUTSTREAM_H_
#define READNONBLANKLINEFROMINPUTSTREAM_H_

#include <istream>
#include <string>

/**
 * read a non blank line from given input stream
 * @param in the input stream
 * @return read line into a string
 */
std::string readNonBlankLineFromInputStream(std::istream & in);

#endif /* READNONBLANKLINEFROMINPUTSTREAM_H_ */
