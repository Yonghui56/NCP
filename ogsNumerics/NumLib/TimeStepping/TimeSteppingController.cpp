/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TimeSteppingController.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "TimeSteppingController.h"

#include <iostream>

#include "logog.hpp"

namespace NumLib
{

size_t TimeSteppingController::solve(double time_end) 
{
    TimeStep time_current(_time_begin);
    //INFO("\n\n***Start time steps with t0=%d s\n", _time_begin);

    while (time_current.getTime()<time_end) {
        double time_next = _root_subsystems->suggestNext(time_current);
        if (!(time_next > time_current.getTime())) {
            //error
            ERR("error - the suggested next time step is invalid.");
            break;
        }
        if (time_next > time_end) {
            time_next = time_end;
        }
        
        TimeStep t_n1(time_current, time_next);
        INFO("\n#############################################################");
        INFO("Time step %d: t=%f s, dt=%f s ", t_n1.getTimeStepCount(), time_next, t_n1.getTimeStepSize());
        INFO("#############################################################");
        _root_subsystems->solveTimeStep(t_n1);
        if (_root_subsystems->accept(t_n1)) {
            _root_subsystems->finalizeTimeStep(time_next);
            time_current.finalize(time_next);
            doSomethingAfterTimeStepAccepted(time_current);
        }
    }

    return time_current.getTimeStepCount();
}

}
