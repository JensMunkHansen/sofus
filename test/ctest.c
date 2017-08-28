#include <fnm/config.h>
#include <fnm/if_fnm.h>

#include <sps/multi_malloc.h>

#include <stdint.h>
#include <stdio.h>

#include <complex.h>

#include <setjmp.h>
#include <cmocka.h>

int fnm_ctest_pressure_linear_array()
{

//! [LinearArrayC example]
  const float f0 = 1.0e6f;
  const float c  = 1500.0f;
  const size_t nElements = 64; // 128
  const float kerf   = 5.0e-4f;
  const float width  = 3e-3f;
  const float height = 50e-3f;

  size_t nDiv = 4;

  const size_t nx = 170;
  const size_t nz = 250;

  float* pos = (float*) malloc(nx*nz*3*sizeof(float));

  const float d = (width+kerf)*nx;

  const float dx = (1.5f * d) / nx;
  const float dz = (2.0f * d) / nz;

  float focus[3] = {0,0,d};

  float* xs = (float*) malloc(nx*sizeof(float));
  float* zs = (float*) malloc(nz*sizeof(float));

  float wx = (nx-1.0f) / 2;

  size_t i,j;
  for (i = 0 ; i < nx ; i++) {
    xs[i] = ( (float)i - wx) * dx;
  }
  for (i = 0 ; i < nz ; i++) {
    zs[i] = (float)i * dz;
  }

  Aperture* a = NULL;
  int err = ApertureLinearCreate(&a, nElements, width, kerf, height);

#ifdef FNM_CLOSURE_FUNCTIONS
  ApertureRwFloatParamSet(a,F0,&f0,0);
  ApertureRwFloatParamSet(a,C,&c,0);
  ApertureRwFloatParamSet(a,Focus, &focus[0], 1);
#endif
  // TODO: Narrow down interface
  ApertureNDivWSet(a,nDiv);
  ApertureNDivHSet(a,nDiv);
  ApertureNThreadsSet(a,1);

  ApertureFocusingTypeSet(a,Rayleigh);

  for (i = 0 ; i < nx ; i++) {
    for (j = 0 ; j < nz ; j++) {
      pos[3*(i*nz+j) + 0] = xs[i];
      pos[3*(i*nz+j) + 1] = 0.0f;
      pos[3*(i*nz+j) + 2] = zs[j];
    }
  }

  size_t nresults = 0;
  float complex* results = NULL;

#ifndef FNM_DOUBLE_SUPPORT
  err = ApertureCalcCwFieldRef(a, pos, nx*nz, 3, &results, &nresults);
#else
  err = ApertureCalcCwFast(a, pos, nx*nz, 3, &results, &nresults);
#endif

  // Clean-up
  err = ApertureDestroy(a);

  if (results) {
    free(results);
  }

  free(pos);
  free(xs);
  free(zs);

  //! [LinearArrayC example]
  SPS_UNREFERENCED_PARAMETER(err);

  return 0;
}

void dummy()
{

}
int main(int argc, char* argv[])
{

  const struct CMUnitTest tests[] = {
#ifdef FNM_CLOSURE_FUNCTIONS
    unit_test(fnm_ctest_pressure_linear_array),
#endif
    unit_test(dummy),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
