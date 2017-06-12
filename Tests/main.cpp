/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file main.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include <gtest/gtest.h>
#include "logog.hpp"
#include "logog/include/formatter.hpp"

/**
 * new formatter for logog
 */
class FormatterCustom : public logog::FormatterGCC
{
    virtual TOPIC_FLAGS GetTopicFlags( const logog::Topic &topic )
    {
        return ( Formatter::GetTopicFlags( topic ) &
                 ~( TOPIC_LEVEL_FLAG | TOPIC_FILE_NAME_FLAG | TOPIC_LINE_NUMBER_FLAG ));
    }
};

//int add (int x, int y) {return x+y;};
//
//TEST(Dummy, Test1)
//{
//    ASSERT_EQ(2, add(1, 1));
//}


int main(int argc, char *argv[])
{
    int ret = 0;
    LOGOG_INITIALIZE();
    try {
        logog::Cout out;
        FormatterCustom custom_format;
        out.SetFormatter(custom_format);

#if 0
    argc = 2;
    //argv[1] = "--gtest_filter=Math.Nonlinear*";
    //argv[1] = "--gtest_filter=Num.Discrete*:FEM.*";
    //argv[1] = "--gtest_filter=Discrete.NDDC*";
    //argv[1] = "--gtest_filter=Solution.Coupling*";
    //argv[1] = "--gtest_filter=Math.Matrix*";
    //argv[1] = "--gtest_filter=Coupling.*";
    //argv[1] = "--gtest_filter=Math.SystemOfEqs*";
    //argv[1] = "--gtest_filter=FEM.LIE*";
    //argv[1] = "--gtest_filter=Solution.CouplingF*";
    //argv[1] = "--gtest_filter=*";
#endif
        ::testing::InitGoogleTest(&argc, argv);
        ret = RUN_ALL_TESTS();
    } catch (char* e) {
        std::cerr << e << std::endl;
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
    } catch (...) {
        std::cerr << "Unknown exception occurred!" << std::endl;
    }
    LOGOG_SHUTDOWN();

    return ret;
}

