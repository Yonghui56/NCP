/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TestUtil.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#ifndef _TEST_UTIL_H_
#define _TEST_UTIL_H_

#include <gtest/gtest.h>
#include "GeoLib/Point.h"
#include "MathLib/Vector.h"

inline void ASSERT_DOUBLE_ARRAY_EQ(const double* Expected, const double* Actual, size_t N, double epsilon=1.0e-8) {
    for (size_t i=0; i<N; i++) \
        ASSERT_NEAR(Expected[i], Actual[i], epsilon); 
}

inline void ASSERT_DOUBLE_ARRAY_EQ(double Expected, const double* Actual, size_t N, double epsilon=1.0e-8) {
    for (size_t i=0; i<N; i++) \
        ASSERT_NEAR(Expected, Actual[i], epsilon); 
}

inline void ASSERT_DOUBLE_ARRAY_EQ(const MathLib::Vector2D &Expected, const MathLib::Vector2D* Actual, size_t N, double epsilon=1.0e-8) {
    for (size_t i=0; i<N; i++) 
        for (size_t j=0; j<2; ++j) 
            ASSERT_NEAR(Expected[j], Actual[i][j], epsilon); 
}

inline void ASSERT_DOUBLE_ARRAY_EQ(GeoLib::Point Expected, GeoLib::Point Actual, double epsilon=1.0e-8) {
    for (size_t i=0; i<3; i++) \
        ASSERT_NEAR(Expected[i], Actual[i], epsilon); 
}

template<typename T1, typename T2>
inline void ASSERT_DOUBLE_ARRAY_EQ(T1 &Expected, T2 &Actual, size_t N, double epsilon=1.0e-8) {
    for (size_t i=0; i<N; i++) \
        ASSERT_NEAR(Expected[i], Actual[i], epsilon); 
};

#endif
