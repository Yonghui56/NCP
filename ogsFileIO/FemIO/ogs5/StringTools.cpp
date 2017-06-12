/*
 * StringTools.cpp
 *
 *  Created on: Jun 16, 2010
 *      Author: TF
 */

#include "StringTools.h"

namespace ogs5
{

std::list<std::string> splitString(const std::string &str, char delim)
{
	std::list<std::string> strList;
	std::stringstream ss(str);
	std::string item;
	while(getline(ss, item, delim))
		strList.push_back(item);
	return strList;
}

std::string replaceString(const std::string &searchString,
                          const std::string &replaceString,
                          std::string stringToReplace)
{
	std::string::size_type pos = stringToReplace.find(searchString, 0);
	int intLengthSearch = searchString.length();

	while (std::string::npos != pos)
	{
		stringToReplace.replace(pos, intLengthSearch, replaceString);
		pos = stringToReplace.find(searchString, 0);
	}
	return stringToReplace;
}

void trim(std::string &str, char ch)
{
	std::string::size_type pos = str.find_last_not_of(ch);
	if(pos != std::string::npos)
	{
		str.erase(pos + 1);
		pos = str.find_first_not_of(ch);
		if(pos != std::string::npos)
			str.erase(0, pos);
	}
	else
		str.erase(str.begin(), str.end());
}


}
