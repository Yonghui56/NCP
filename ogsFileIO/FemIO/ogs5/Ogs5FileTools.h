/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Ogs5FileTools.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <string>

#include "makros.h"
#include "readNonBlankLineFromInputStream.h"

inline std::string GetUncommentedLine(std::string& line)
{
    std::string zeile = "";
    int i = 0, j = 0;
    //----------------------------------------------------------------------
    i = (int) line.find_first_not_of(" ",0); //Anf���ngliche Leerzeichen ���berlesen, i=Position des ersten Nichtleerzeichens im string
    j = (int) line.find(";",i);           //Nach Kommentarzeichen ; suchen. j = Position des Kommentarzeichens, j=-1 wenn es keines gibt.
    if((i != -1))
        zeile = line.substr(i,j - i);  //Ab erstem nicht-Leerzeichen bis Kommentarzeichen rauskopieren in neuen substring, falls Zeile nicht leer ist
    i = (int) zeile.find_last_not_of(" "); // Suche nach dem letzten Zeichen, dass kein Leerzeichen ist
    if(i >= 0)
    {
        line = zeile.substr(0,i + 1); // Leerzeichen am Ende rausschneiden
        zeile = line;
    }

    return zeile;
}

inline std::ios::pos_type GetNextSubKeyword(std::ifstream* file,std::string* line, bool* keyword)
{
    char buffer[MAX_ZEILE];
    std::ios::pos_type position;
    position = file->tellg();
    *keyword = false;
    std::string line_complete;
    int i,j;
    // Look for next subkeyword
    while(!(line_complete.find("$") != std::string::npos) && (!file->eof()))
    {
        file->getline(buffer,MAX_ZEILE);
        line_complete = buffer;
        if(line_complete.find("#") != std::string::npos)
        {
            *keyword = true;
            return position;
        }
        i = (int) line_complete.find_first_not_of(" ",0);
        j = (int) line_complete.find(";",i); //Nach Kommentarzeichen ; suchen. j = Position des Kommentarzeichens, j=-1 wenn es keines gibt.
        if(j < 0)
            j = (int)line_complete.length();
        //if(j!=i) break;
        if(i != -1)
            *line = line_complete.substr(i,j - i);  //Ab erstem nicht-Leerzeichen bis Kommentarzeichen rauskopieren in neuen substring, falls Zeile nicht leer ist
    }
    return position;
}
