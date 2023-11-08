#include <gtest/gtest.h>

#include <sps/memory>
#ifdef _WIN32
# include <sps/win32/memory>
#endif
#include <sps/unix/memory>
#include <fnm/config.h>
#include <fnm/fnm_data.hpp>
#include <fnm/fnm.hpp>
#include <fnm/circular.hpp>

#include <iostream>

// No leaks
TEST(fnm_test, test_nothing) {
  EXPECT_EQ(1, 1);
}

class A {
 public:
  A() {
    m_data = sps::unique_aligned_array_create<float>(20);
  }
 private:
  sps::unique_aligned_array<float> m_data;
};

class B {
 public:
  B() : m_data(2,2) {
    m_data = sps::unique_aligned_multi_array<sps::element_rect_t<float>, 2>(4,4);
  }
 private:
  sps::unique_aligned_multi_array<sps::element_rect_t<float>, 2> m_data;
};

// No leaks
TEST(fnm_test, unique_members) {
  A a;
  B b;
}

class dummy {
 public:
  dummy() : m_data() {}
  dummy(sps::nix::unique_aligned_multi_array<float,2>&& data) /* noexcept */ : m_data(std::move(data)) {}
  sps::nix::unique_aligned_multi_array<float,2> m_data;
};

TEST(fnm_test, multi_arrays) {
  auto arr1 = sps::unique_aligned_multi_array<float, 2>(3,4);

  for (size_t i = 0 ; i < 3 ; i++)
    for (size_t j = 0 ; j < 4 ; j++)
      arr1[i][j] = float(i*4 + j);

  // Example using non-contiguous memory
  std::unique_ptr<int*, std::function<void(int**)>> x(
        new int*[10](),
  [](int** x) {
    std::for_each(x, x + 10, std::default_delete<int[]>());
    delete[] x;
  });

  for (size_t row = 0; row < 10; ++row) {
    (x.get())[row] = new int[10];
  }

  size_t k = 10;
  for (size_t i = 0 ; i < k ; i++)
    for (size_t j = 0 ; j < 10 ; j++)
      x.get()[i][j] = int(i*10 + j);

  // Example using contiguous memory
  std::unique_ptr<int*, std::function<void(int**)>> y(new int*[k](),
  [](int** x) {
    _mm_free(&(x[0][0]));
    delete[] x;
  });
  y.get()[0] = (int*) _mm_malloc(k*10*sizeof(int),16);
  for (size_t row = 1; row < k; ++row) {
    (y.get())[row] = &(y.get()[0][0]);
  }

  for (size_t i = 0 ; i < k ; i++)
    for (size_t j = 0 ; j < 10 ; j++)
      y.get()[i][j] = int(i*10 + j);

  auto arr3 = sps::nix::unique_aligned_multi_array_create<float, 2>(10,10);
  for (size_t i = 0 ; i < k ; i++) {
    for (size_t j = 0 ; j < 10 ; j++) {
      arr3[i][j] = float(i*10 + j);
      arr3(i,j) =  float(i*10 + j);
    }
  }
  EXPECT_EQ(1,1);
}

TEST(fnm_test, allocate) {
  fnm::ApertureData<float> a;
  fnm::Aperture<float> b;

  fnm::Aperture<float> c(10,1.0f,0.0f,2.0f);

  fnm::Aperture<float> d;

  EXPECT_EQ(1,1);
}

TEST(fnm_test, elements_get) {
  auto aperture = fnm::Aperture<float>(64, 1.0f, 0.0f, 1.0f);
  const sps::element_rect_t<float>** elements = NULL;
  size_t nElements, nSubElements;
  aperture.ElementsRefGet(&nElements, &nSubElements, elements);
  EXPECT_EQ(elements[0][0].center[0], -31.5f);
}

// Leaks 48 bytes
TEST(fnm_test, scalar_float_set) {
#ifdef FNM_CLOSURE_FUNCTIONS

  const size_t nElements = 12;
  auto aperture = fnm::Aperture<float>(nElements, 1.0f, 0.0f, 1.0f);

  float alpha = 4.0f;

  aperture.RwFloatParamSet(fnm::RwParamType::Alpha, &alpha, 0);

  float beta = 7.0f;
  aperture.RwFloatParamSet(fnm::RwParamType::Beta, &beta, 0);

  float* fDelays = (float*)malloc(12 * sizeof(float));

  fDelays[2] = 2.0f;
  aperture.RwFloatParamSet(fnm::RwParamType::ElementDelays,
                           fDelays, 1, nElements);
  if (fDelays) {
    free(fDelays);
  }

  float* pfDelays = NULL;
  size_t nOutElements;
  aperture.RwFloatParamGet(fnm::RwParamType::ElementDelays,
                           &pfDelays, 1, &nOutElements);

  if (pfDelays) {
    free(pfDelays);
    pfDelays = NULL;
  }

  float* pfAlpha = NULL;

  aperture.RwFloatParamGet(fnm::RwParamType::Alpha, &pfAlpha, 0, &nOutElements);

  EXPECT_EQ(alpha, *pfAlpha);

  if (pfAlpha) {
    free(pfAlpha);
  }

  float* focus = (float*) malloc(3*sizeof(float));
  focus[0] = 1.0f;
  focus[1] = 2.0f;
  focus[2] = 3.0f;
  aperture.RwFloatParamSet(fnm::RwParamType::Focus, focus, 1, 3);

  free(focus);

  aperture.RwFloatParamGet(fnm::RwParamType::Focus, &pfDelays, 1, &nOutElements);
  if (pfDelays) {
    free(pfDelays);
  }

  aperture.RwBooleanParamSet0D(fnm::RwParamType::AttenuationEnabled, true);

  bool enabled = aperture.AttenuationEnabledGet();
  SPS_UNREFERENCED_PARAMETER(enabled);
#endif
}

TEST(fnm_test, pressure_linear_array) {
//! [LinearArray example]
  const float f0 = 1.0e6f;
  const float c  = 1500;
  const size_t nElements = 128;
  const float kerf   = 5.0e-4f;
  const float width  = 3e-3f;
  const float height = 50e-3f;

  size_t nDiv = 4;

  const size_t nx = 170;
  const size_t nz = 250;

  auto pos = sps::unique_aligned_array_create<float>(nx*nz*3);

  const float d = (width+kerf)*nx;

  const float dx = (1.5f * d) / nx;
  const float dz = (2.0f * d) / nz;

  float focus[3] = {0, 0, d};

  std::vector<float> xs(nx);
  std::vector<float> zs(nz);

  float wx = static_cast<float>(nx-1) / 2;
  for (size_t i = 0 ; i < nx ; i++) {
    xs[i] = (static_cast<float>(i) - wx) * dx;
  }
  for (size_t i = 0 ; i < nz ; i++) {
    zs[i] = static_cast<float>(i) * dz;
  }

  fnm::Aperture<float> a(nElements, width, kerf, height);

  a.F0Set(f0);
  a.CSet(c);
  a.NDivWSet(nDiv);
  a.NDivHSet(nDiv);
  a.NThreadsSet(8);
  a.FocusSet(focus);
  a.FocusingTypeSet(fnm::FocusingType::Rayleigh);

  for (size_t i = 0 ; i < nx ; i++) {
    for (size_t j = 0 ; j < nz ; j++) {
      pos[3*(i*nz+j) + 0] = xs[i];
      pos[3*(i*nz+j) + 1] = 0.0f;
      pos[3*(i*nz+j) + 2] = zs[j];
    }
  }

  // Output variables
  size_t nresults = 0;
  std::complex<float>* results = NULL;

  // Compute pressure
  int err;
  err = a.CalcCwFast(pos.get(), nx*nz, 3, &results, &nresults);
  SPS_UNREFERENCED_PARAMETER(err);

  // Clean-up
  free(results);
  //! [LinearArray example]
}

TEST(fnm_test, circular_array) {
  fnm::CircularAperture<float> aperture(1.2f);
}

int main(int argc, char* argv[]) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
