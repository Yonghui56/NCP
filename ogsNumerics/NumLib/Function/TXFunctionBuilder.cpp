/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TXFunctionBuilder.cpp
 *
 * Created on 2012-11-15 by Norihiro Watanabe
 */

#include "TXFunctionBuilder.h"

#include <string>
#include "logog.hpp"
#include "MathLib/Interpolation/LinearInterpolation.h"
#include "Multiplication.h"

namespace NumLib
{

ITXFunction* TXFunctionBuilder::create(const BaseLib::Options &opDistribution)
{
    const std::string dis_name = opDistribution.getOption("DistributionType");
    NumLib::ITXFunction* f_bc = NULL;
    NumLib::TXFunctionType::type dis_type = NumLib::convertStringToTXFunctionType(dis_name);
    switch (dis_type)
    {
        case NumLib::TXFunctionType::CONSTANT:
            {
                double dis_v = opDistribution.getOptionAsNum<double>("DistributionValue");
                f_bc =  new NumLib::TXFunctionConstant(dis_v);
            }
            break;
        case NumLib::TXFunctionType::GEOSPACE:
            {
//                size_t n_pt = opDistribution.getOptionAsNum<size_t>("PointSize");
//                std::vector<const BaseLib::OptionGroup*> list_opDistLinear = opDistribution.getSubGroupList("PointValue");
//                assert (n_pt == list_opDistLinear.size());
//                for (size_t j=0; j<n_pt; j++) {
//                    size_t pt_id = list_opDistLinear[j]->getOptionAsNum<size_t>("PointID");
//                    double pt_val = list_opDistLinear[j]->getOptionAsNum<double>("Value");
//                }
            }
            break;
        case NumLib::TXFunctionType::ANALYTICAL:
            {
                std::string str_f_exp =opDistribution.getOption("DistributionFunction");
                f_bc = new NumLib::TXFunctionAnalytical(str_f_exp);
            }
            break;
        case NumLib::TXFunctionType::T_LINEAR:
            {
                std::vector<double> vec_t, vec_val;
                {
                    // structure is assumed
                    // + TimeValue
                    //     + Point: Time=0, Value=0
                    //     + Point: Time=10, Value=20
                    const BaseLib::OptionGroup* opDistLinear = opDistribution.getSubGroup("TimeValue");
                    std::vector<const BaseLib::OptionGroup*> list_opDistTime = opDistLinear->getSubGroupList("Point");
                    for (size_t j=0; j<list_opDistTime.size(); j++) {
                        double pt_t = list_opDistTime[j]->getOptionAsNum<double>("Time");
                        double pt_val = list_opDistTime[j]->getOptionAsNum<double>("Value");
                        vec_t.push_back(pt_t);
                        vec_val.push_back(pt_val);
                    }
                }
                MathLib::LinearInterpolation* time_curve = new MathLib::LinearInterpolation(vec_t, vec_val);
                f_bc = new NumLib::TXFunctionTimeCurve<MathLib::LinearInterpolation>(time_curve);
            }
            break;
        case NumLib::TXFunctionType::TX:
            {
                const BaseLib::Options* opTF = opDistribution.getSubGroup("Temporal");
                ITXFunction* fT = create(*opTF);
                const BaseLib::Options* opXF = opDistribution.getSubGroup("Spatial");
                ITXFunction* fX= create(*opXF);
                f_bc = new NumLib::TXCompositFunction<ITXFunction,ITXFunction, Multiplication>(fT, fX);
            }
            break;
        default:
            ERR("***ERROR in FemVariableBuilder. Distribution type %s is not supported.", dis_name.c_str());
            break;
    }
    return f_bc;
}


}
