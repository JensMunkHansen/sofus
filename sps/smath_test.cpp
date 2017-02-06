#include <sps/sps_export.h>
#include <sps/math.h>
#include <sps/smath.hpp>

#include <float.h>

#include <cstdlib>
#include <iostream>
#include <cassert>

#include <stdint.h>

#include <gtest/gtest.h>

using namespace sps;


TEST(smath_test, test_basis_vectors_double)
{
  euler_t<double> euler; //  = { 0.1, 0.2, 0.3, 0.0 }; Not working on _MSC_VER
  euler.alpha = 0.1;
  euler.beta  = 0.2;
  euler.gamma = 0.3;
  double vec0[4], vec1[4], vec2[4];

  basis_vectors(vec0, vec1, vec2, euler);

  point_t<double> vec = point_t<double>();
  basis_vectors(vec, euler, 0);

  double diff = sqrt(SQUARE(vec[0] - vec0[0]) +   SQUARE(vec[1] - vec0[1]) +   SQUARE(vec[2] - vec0[2]));
  ASSERT_LT(diff, FLT_EPSILON);
  basis_vectors(vec, euler, 1);
  diff = sqrt(SQUARE(vec[0] - vec1[0]) +   SQUARE(vec[1] - vec1[1]) +   SQUARE(vec[2] - vec1[2]));
  ASSERT_LT(diff, FLT_EPSILON);
  basis_vectors(vec, euler, 2);
  diff = sqrt(SQUARE(vec[0] - vec2[0]) +   SQUARE(vec[1] - vec2[1]) +   SQUARE(vec[2] - vec2[2]));
  ASSERT_LT(diff, FLT_EPSILON);
}


int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
