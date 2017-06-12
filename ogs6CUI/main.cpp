/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file main.cpp
 *
 * Created on 2012-07-17 by Norihiro Watanabe
 */

#include <iostream>
#include <exception>
#include "THMCSimulator.h"

int main ( int argc, char *argv[] )
{
    int returncode = 0;
    try {
        ogs6::THMCSimulator sim(argc, argv);
        returncode = sim.execute();
    } catch (const char* e) {
        std::cerr << "EXCEPTION: " << e << std::endl;
    } catch (std::exception& e) {
        std::cerr << "EXCEPTION: " << e.what() << std::endl;
    } catch (...) {
        std::cerr << "EXCEPTION: Unknown exception occurred!" << std::endl;
    }

    return returncode;
}
