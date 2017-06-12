/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementVelocity.hpp
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

//#include "ElementVelocity.h"

#include "logog.hpp"
#include "MathLib/Vector.h"
#include "DiscreteLib/Utils/DiscreteSystemContainerPerMesh.h"
#include "MaterialLib/PorousMedia.h"

#include "Ogs6FemData.h"

template <class T>
bool FunctionPressureToHead<T>::initialize(const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    if (femData->list_fluid.size()==0) {
        ERR("***Error: No fluid properties found.");
        return false;
    } else {
        if (femData->list_fluid[0]->density==NULL) {
            ERR("***Error: Fluid density is not configured.");
            return false;
        }
    }

    size_t msh_id = option.getOptionAsNum<size_t>("MeshID");
    MeshLib::IMesh* msh = femData->list_mesh[msh_id];
    MyDiscreteSystem* dis = 0;
    dis = DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(msh);

    _dis = dis;
    _h = new MyNodalFunctionScalar();
    _h->initialize(*dis, FemLib::PolynomialOrder::Linear, .0);
    
    // set initial output
    OutputVariableInfo var(this->getOutputParameterName(Head), _dis->getMesh()->getID(), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _h);
    femData->outController.setOutput(var.name, var); 

    // initial output parameter
    this->setOutput(Head, _h);

    return true;
}

template <class T>
void FunctionPressureToHead<T>::finalizeTimeStep(const NumLib::TimeStep &/*time*/)
{
    //std::cout << "Velocity=" << std::endl;
    //_vel->printout();
    //update data for output
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    OutputVariableInfo var(this->getOutputParameterName(Head), _dis->getMesh()->getID(), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _h);
    femData->outController.setOutput(var.name, var); 

};

template <class T>
int FunctionPressureToHead<T>::solveTimeStep(const NumLib::TimeStep &/*time*/)
{
    INFO("Converting Pressure to Hydraulic head...");

    const MeshLib::IMesh *msh = _dis->getMesh();
    MyNodalFunctionScalar *f_p = (MyNodalFunctionScalar*)getInput(Pressure);
    MyNodalFunctionScalar *f_h = _h;;

    Ogs6FemData* femData = Ogs6FemData::getInstance();
    MaterialLib::Fluid* fluid = femData->list_fluid[0];
    const double g = 9.81;

    for (size_t i=0; i<msh->getNumberOfNodes(); i++) {
        double p = f_p->getValue(i);
        double z = (*msh->getNodeCoordinatesRef(i))[2];
        double rho_f = .0;
        fluid->density->eval(NumLib::TXPosition(NumLib::TXPosition::Node, i), rho_f);
        double h = p / (rho_f * g) - z;
        f_h->setValue(i, h);
    }

    setOutput(Head, f_h);

    return 0;
}

