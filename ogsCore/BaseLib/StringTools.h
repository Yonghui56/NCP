/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file StringTools.h
 *
 * Created on 2010-06-16 by Thomas Fischer
 */

#ifndef STRINGTOOLS_H
#define STRINGTOOLS_H

#include <string>
#include <list>
#include <sstream>
#include <fstream>
#include <iostream>
#include <ctype.h>

namespace BaseLib
{

/**
 *   Splits a string into a list of strings.
 *  \param str String to be splitted
 *  \param delim Character indicating that the string should be splitted
 *  \return
 */
std::list<std::string> splitString(const std::string &str, char delim);

/**
 *   Replaces a substring with another in a string
 *  \param searchString Search for this string
 *  \param replaceString Replace with this string
 *  \param stringToReplace Search and replace in this string
 *  \return The modified string
 */
std::string replaceString(const std::string &searchString, const std::string &replaceString, std::string stringToReplace);

/**
 *   Converts a number (double, float, int, ...) into a string
 *  \param d The number to be converted
 *  \return The number as string
 */
template<typename T> std::string number2str(T d)
{
    std::stringstream out;
    out << d;
    return out.str();
}

/**
 *   Converts a string into a number (double, float, int, ...)
 *  Example: size_t number (str2number<size_t> (str));
 *  \param str string to be converted
 *  \return the number
 */
template<typename T> T str2number (const std::string &str)
{
    std::stringstream strs (str, std::stringstream::in | std::stringstream::out);
    T v;
    strs >> v;
    return v;
}

/**
 * Strip whitespace (or other characters) from the beginning and end of a string.
 */
void trim(std::string &str, char ch=' ');

#ifdef MSVC
void correctScientificNotation(std::string filename, size_t precision = 0);
#endif

/**
 * Left padding
 * @param src
 * @param total_len
 * @param paddingChar
 * @return
 */
std::string leftPadding(const std::string &src, size_t total_len, const char paddingChar =' ');

/**
 * Right padding
 * @param src
 * @param total_len
 * @param paddingChar
 * @return
 */
std::string rightPadding(const std::string &src, size_t total_len, const char paddingChar =' ');

/**
 * Left and right padding and make string at center
 * @param src
 * @param total_len
 * @param paddingChar
 * @return
 */
std::string bothPadding(const std::string &src, size_t total_len, const char paddingChar = ' ');

}

#endif //STRINGTOOLS_H
