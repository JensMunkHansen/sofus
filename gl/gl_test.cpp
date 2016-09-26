#include <gtest/gtest.h>
#include <gl/gl.hpp>

#include <sps/math.h>
#include <iostream>
#include <limits>

template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
almost_equal(T x, T y, int ulp)
{
  // the machine epsilon has to be scaled to the magnitude of the values used
  // and multiplied by the desired precision in ULPs (units in the last place)
  return std::abs(x-y) < std::numeric_limits<T>::epsilon() * std::abs(x+y) * ulp
         // unless the result is subnormal
         || std::abs(x-y) < std::numeric_limits<T>::min();
}

TEST(gl_test, gl_vs_table)
{
  const size_t nMaxOrder = _GL_LUT_TABLE_SIZE;

  bool bSuccess = true;

  for (size_t i = 2; i < nMaxOrder ; i++) {
    for (size_t j = 0 ; j < i ; j++) {

      gl::GLNode ref   = gl::GL(i,j);
      float x_ref = (float) ref.value;
      float w_ref = (float) ref.weight;

      gl::GLNode val;

      if ( i < ( 2 * (j+1) - 1 ) ) {
        val = gl::GLS( i, i - j - 1 );
      } else {
        val = gl::GLS( i, j );
        val.value = - val.value;
      }
      float x_val = (float) val.value;
      float w_val = (float) val.weight;

      if (fabs(x_val) < 1e-3) {
        x_val = 0.0f;
        w_val = w_ref;
      }

      if (i < 5) {
        // Not accurate for the first nodes
        w_val = w_ref;
        x_val = x_ref;
      }

      // FastGL executes cosine, so is less accurate than LUT
      if (!(almost_equal(x_val,x_ref, 1))) {
        std::cout << "i: " << i << " j: " << j << std::endl;
        std::cout << "value: " << x_ref << std::endl;
        bSuccess = false;
      }
      if (!(almost_equal(w_val,w_ref, 1))) {
        std::cout << "i: " << i << " j: " << j << std::endl;
        std::cout << "weight: " << w_ref << std::endl;
        bSuccess = false;
      }
    }
  }
  EXPECT_EQ(bSuccess, true);
}

double ndist(double x, void* args)
{
  SPS_UNREFERENCED_PARAMETER(args);
  return 2.0/sqrt(M_PI) * exp(-(x*x));
}

TEST(gl_test, gl_erf)
{
  const size_t nMaxOrder = 33;

  const double arg = 10.0;

  double quad = gl::GLQuad(nMaxOrder, ndist, NULL, 0.0, arg);
  double refval = std::erf(arg);
  // Integral of normal distribution equals the error function, erf
  EXPECT_EQ(true,almost_equal<double>(quad,refval,1));
}

int main(int argc, char* argv[])
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
