#include <cstdio>
#include <cstdlib>
#include <sps/trigintrin.h>
#include <cstring>
#include <algorithm> // std::min, std::max
#include <cassert>

#include <gtest/gtest.h>

TEST(trigintrin_test, test_sin_cos_log)
{
  __m128 a = _mm_set1_ps(2.0f);
  __m128 b = _mm_set1_ps(4.0f);
  _mm_arccos_ps(a);
  _mm_arcsin_ps(a);
  _mm_arctan2_ps(a,b);


  ALIGN16_BEGIN float vout[4] ALIGN16_END;

  float max_diff[] = {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f, 0.0f};
  float x,diff;
  __m128 vx,vr;

  for (size_t i = 0 ; i < 1000 ; i++) {
    x  = -float(M_PI)/2.0f + float(i)/1000 * float(M_PI);
    vx = _mm_set1_ps(x);
    vr = _mm_cos_ps(vx);
    memcpy(vout,(void*)&vr,16);
    diff = fabs(vout[0]-cos(x));
    max_diff[3] = std::max<float>(diff, max_diff[3]);
  }

  for (size_t i = 0 ; i < 1000 ; i++) {
    x  = -float(M_PI)/2.0f + float(i)/1000 * float(M_PI);
    vx = _mm_set1_ps(x);
    vr = _mm_sin_ps(vx);
    memcpy(vout,(void*)&vr,16);
    diff = fabs(vout[0]-sin(x));
    max_diff[4] = std::max<float>(diff, max_diff[4]);
  }

  __m128 vr1 = _mm_setzero_ps();
  for (size_t i = 0 ; i < 1000 ; i++) {
      x = -float(M_PI) / 2.0f + float(i) / 1000 * float(M_PI);
    vx = _mm_set1_ps(x);
    _mm_sin_cos_ps(vx,&vr,&vr1);
    memcpy(vout,(void*)&vr,16);
    diff = fabs(vout[0]-sin(x));
    max_diff[5] = std::max<float>(diff, max_diff[5]);
    memcpy(vout,(void*)&vr1,16);
    diff = fabs(vout[0]-cos(x));
    max_diff[6] = std::max<float>(diff, max_diff[6]);

  }

  for (size_t i = 0 ; i < 1000 ; i++) {
    x = exp(float(i+1)/12.0f);
    vx = _mm_set1_ps(x);
    vr = _mm_log_ps(vx);
    //vr = (__m128) xlogf((vfloat)vx);
    memcpy(vout,(void*)&vr,16);
    diff = fabs(vout[0]-logf(x));
    max_diff[7] = std::max<float>(diff, max_diff[7]);
  }

  for (size_t i = 0 ; i < 1000 ; i++) {
    x = float(i);
    vx = _mm_set1_ps(x);
    vr = _mm_cbrtf_ps(vx);
    //vr = xcbrtf(vx);

    memcpy(vout,(void*)&vr,16);
    diff = fabs(vout[0] - cbrtf(x));
    max_diff[8] = std::max<float>(diff, max_diff[8]);
  }


  //  printf("cos error: %f\n",max_diff[3]);
  //  printf("sin error: %f\n",max_diff[4]);
  //  printf("sin_cos error: %f %f\n",max_diff[5],max_diff[6]);
  //  printf("log error: %f\n",max_diff[7]);
  //  printf("cbrtf error: %f\n",max_diff[8]);

  ASSERT_LT( max_diff[0]/2.0f, 6.7e-5);
}

TEST(trigintrin_test, arcsin_arccos_arctan2)
{
  ALIGN16_BEGIN float vout[4] ALIGN16_END;

  float max_diff[] = {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f, 0.0f};
  float x,diff;
  __m128 vx,vr;

  for (size_t i = 0 ; i < 1000 ; i++) {
    x  = -1.0f + float(i)/1000 * 2.0f;
    vx = _mm_set1_ps(x);
    vr = _mm_arccos_ps(vx);
    memcpy(vout,(void*)&vr,16);
    diff = fabs(vout[0]-acos(x));
    if (diff < 3.14)
      max_diff[0] = std::max<float>(diff, max_diff[0]);

    vr = _mm_arcsin_ps(vx);
    memcpy(vout,(void*)&vr,16);
    diff = fabs(vout[0]-asin(x));
    if (diff < 3.14)
      max_diff[1] = std::max<float>(diff, max_diff[1]);
  }

  for (size_t i = 0 ; i < 1000 ; i++) {
    x  = -100.0f + float(i)/1000 * 200.0f;
    vx = _mm_set1_ps(x);
    vr = _mm_arctan2_ps(vx,_m_one_ps);
    memcpy(vout,(void*)&vr,16);
    diff = fabs(vout[0]-atan2(x,1.0f));
    max_diff[2] = std::max<float>(diff, max_diff[2]);
  }

  assert((max_diff[0]/2.0f < 6.7e-5) && "_mm_arcsin_ps too inaccurate");
  assert((max_diff[1]/2.0f < 6.7e-5) && "_mm_arccos_ps too inaccurate");
  assert((max_diff[2]/2.0f < 6.7e-5) && "_mm_arctan2_ps too inaccurate");

  ASSERT_LT( max_diff[0]/2.0f, 6.7e-5);
  ASSERT_LT( max_diff[1]/2.0f, 6.7e-5);
  ASSERT_LT( max_diff[2]/2.0f, 6.7e-5);

  //  printf("arccos error: %f\n",max_diff[0]);
  //  printf("arcsin error: %f\n",max_diff[1]);
  //  printf("arctan error: %f\n",max_diff[2]);
}

TEST(trigintrin_test, test_exp)
{
  ALIGN16_BEGIN float vout[4] ALIGN16_END;

  float max_diff[] = {0.0f};
  float x,diff;
  __m128 vx,vr;

  for (size_t i = 0 ; i < 1000 ; i++) {
    x = float(i+1)/100.0f;
    vx = _mm_set1_ps(x);
    vr = _mm_exp_ps(vx);
    memcpy(vout,(void*)&vr,16);
    diff = fabs(vout[0]-exp(x));
    max_diff[0] = std::max<float>(diff, max_diff[0]);
  }
  assert((max_diff[0] < 2.1e-3) && "_mm_exp_ps too inaccurate");

  ASSERT_LT( max_diff[0], 2.1e-3);

  //  printf("exp error: %f\n",max_diff[0]);
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
