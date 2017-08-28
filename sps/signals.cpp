/**
 * @file   signals.cpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Fri May 29 02:43:18 2015
 *
 * @brief
 *
 *
 */

/*
 *  This file is part of SOFUS.
 *
 *  SOFUS is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  SOFUS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with SOFUS.  If not, see <http://www.gnu.org/licenses/>.
 */

// TODO: Implement pruned FFTs
//       Use FFTW's h2h
//         fftw_plan_r2r_1d(N, out, in, FFTW_HC2R / FFTW_R2HC, FFTW_PATIENT);
//         fftw_execute_r2r(plan_r2hc, in, out);

// For a halfcomplex array, hc[n], the kth component thus has its real
// part in hc[k] and its imaginary part in hc[n-k], with the exception
// of k == 0 or n/2 (the latter only if n is even)—in these two cases,
// the imaginary part is zero due to symmetries of the real-input DFT

#include <sps/signals.hpp>
#include <sps/cenv.h>
#include <sps/mm_malloc.h>
#include <sps/debug.h>
#include <cstring>   // memset
#include <algorithm> // min/max

#include <sps/math.h>
#include <sps/extintrin.h>

#include <stdint.h>
#include <emmintrin.h>

#include <cstring>

#include <fftw3.h>

#include <algorithm>
#include <cassert>

#include <complex>

#include <iostream>

#include <mutex>
std::mutex g_plan_mutex;

// TODO: Solve order of destruction sequence

#include <cstdio>
#define DEBUG_PLAN 0

// TODO: Avoid static const members with in-class initializers
#ifdef _MSC_VER
# pragma warning(push)
# pragma warning(disable:4127)
#endif

namespace sps {

// Consider using FFTW_MEASURE or FFTW_PATIENT (overwrites input)

// Erases the input arrays
// FFTW_MEASURE
// Even more optimal
//#define SPS_FFTW_PLAN FFTW_PATIENT
//#define SPS_FFTW_PLAN FFTW_EXHAUSTIVE

#define SPS_FFTW_FAST FFTW_ESTIMATE
#define SPS_FFTW_ACCURATE FFTW_PATIENT

  template <typename T>
  class Signal1DPlan {
  public:
    static Signal1DPlan& Instance();
    //    static void Wisdom(size_t** nFFTSamples, size_t* nFFTs);
  private:
    Signal1DPlan();
    ~Signal1DPlan();
    Signal1DPlan(Signal1DPlan const&) = default;
    Signal1DPlan& operator=(Signal1DPlan const&) = default;
  };

  template<>
  class Signal1DPlan<float> {
  public:

    // TODO: Consider storing plans for non-aligned memory
    static const size_t nFFTLengths = 31;
    static __THREAD fftwf_plan forward[Signal1DPlan<float>::nFFTLengths];
    static __THREAD fftwf_plan backward[Signal1DPlan<float>::nFFTLengths];

    static const bool measure = false;

    static const bool reuse = true;

    static void DestroyPlans()
    {
      for (size_t i = 0 ; i < nFFTLengths ; i++) {
        if (forward[i]) {
          fftwf_destroy_plan(forward[i]);
          forward[i] = NULL;
        }
        if (backward[i]) {
          fftwf_destroy_plan(backward[i]);
          backward[i] = NULL;
        }
      }
    }

    static Signal1DPlan& Instance()
    {
      // A singleton cannot be made thread_local
      static Signal1DPlan singleton;
      return singleton;
    }

    fftwf_plan& Forward(size_t size, float* in, std::complex<float>* out)
    {
      // Check power of two
      assert(size && ((size & (size-1))) ==0);
      size_t index = (size_t)log2(size);

      if (Signal1DPlan<float>::reuse == false) {
        fftwf_destroy_plan(forward[index]);
        forward[index] = NULL;
      }

      if (forward[index]) {
        return forward[index];
      } else {
        std::lock_guard<std::mutex> guard(g_plan_mutex);

        if (Signal1DPlan<float>::measure) {
          float* dummy = (float*) _mm_malloc(size*sizeof(float),16);
          memcpy(dummy, in, size*sizeof(float));
          forward[index] = fftwf_plan_dft_r2c_1d((int)size, dummy, reinterpret_cast<fftwf_complex*>(out), SPS_FFTW_ACCURATE);
          _mm_free(dummy);
        } else {
          forward[index] = fftwf_plan_dft_r2c_1d((int)size, in, reinterpret_cast<fftwf_complex*>(out), SPS_FFTW_FAST);
        }
        return forward[index];
      }
    }

    fftwf_plan& Backward(size_t size, std::complex<float>* in, float* out)
    {
      // Check power of two
      assert(size && ((size & (size-1))) ==0);
      size_t index = (size_t)log2(size);

      if (Signal1DPlan<float>::reuse == false) {
        fftwf_destroy_plan(backward[index]);
        backward[index] = NULL;
      }


      if (backward[index]) {
        return backward[index];
      } else {
        std::lock_guard<std::mutex> guard(g_plan_mutex);

        if (Signal1DPlan<float>::measure) {

          size_t nComplex = size/2 + 1;
          size_t nbytes = 16 * ((nComplex+1)*sizeof(std::complex<float>) + 15) / 16;
          std::complex<float>* dummy = (std::complex<float>*) _mm_malloc(nbytes,16);
          memcpy(dummy, in, nComplex*sizeof(std::complex<float>));
          backward[index] = fftwf_plan_dft_c2r_1d((int)size, reinterpret_cast<fftwf_complex*>(dummy), out, SPS_FFTW_ACCURATE);
          _mm_free(dummy);
        } else {
          backward[index] = fftwf_plan_dft_c2r_1d((int)size, reinterpret_cast<fftwf_complex*>(in), out, SPS_FFTW_FAST);
        }
        return backward[index];
      }
    }

  private:
    Signal1DPlan()
    {
#if USE_FFTW_THREADS
      int state = fftw_init_threads();
      fftw_plan_with_nthreads(4);
#endif
    }
    ~Signal1DPlan()
    {
#if USE_FFTW_THREADS
      fftw_cleanup_threads();
#else
      fftw_cleanup();
#endif

#if DEBUG_PLAN
      FILE* fp = fopen("out.txt","w");
#endif
      for (size_t i = 0 ; i < nFFTLengths ; i++) {
        if (forward[i]) {
#if DEBUG_PLAN
          fftwf_fprint_plan(forward[i], fp);
#endif
          fftwf_destroy_plan(forward[i]);
          forward[i] = NULL;
        }
        if (backward[i]) {
          fftwf_destroy_plan(backward[i]);
          backward[i] = NULL;
        }
      }
#if DEBUG_PLAN
      fclose(fp);
#endif
    }
    Signal1DPlan(Signal1DPlan const&) = default;
    Signal1DPlan& operator=(Signal1DPlan const&) = default;
  };

  __THREAD fftwf_plan Signal1DPlan<float>::forward[Signal1DPlan<float>::nFFTLengths]  = {NULL};
  __THREAD fftwf_plan Signal1DPlan<float>::backward[Signal1DPlan<float>::nFFTLengths] = {NULL};

  template<>
  class Signal1DPlan<double> {
  public:
    static const size_t nFFTLengths = 31;
    static __THREAD fftw_plan forward[Signal1DPlan<double>::nFFTLengths];
    static __THREAD fftw_plan backward[Signal1DPlan<double>::nFFTLengths];

    static Signal1DPlan& Instance()
    {
      // Apparently, the singleton cannot be made thread_local
      static Signal1DPlan singleton;
      return singleton;
    }
    /*
    void Wisdom(size_t** nFFTSamples, size_t* nFFTs) {
      *nFFTs = nFFTLengths;
      size_t* _nFFTSamples = (size_t*) malloc(Signal1DPlan<double>::nFFTLengths * sizeof(size_t));

      for (size_t i = 0 ; i < nFFTLengths ; i++) {
        if (forward[i] != NULL) {
          (*nFFTSamples)[i] = pow(2.0f, (double)i);
        }
      }
      *nFFTSamples = _nFFTSamples;
    }
    */

    fftw_plan& Forward(size_t size, double* in, std::complex<double>* out)
    {
      // Check power of two
      assert(size && ((size & (size-1))) ==0);
      size_t index = (size_t)log2(size);
      if (forward[index]) {
        return forward[index];
      } else {
        std::lock_guard<std::mutex> guard(g_plan_mutex);
        forward[index] = fftw_plan_dft_r2c_1d((int)size, in, reinterpret_cast<fftw_complex*>(out), SPS_FFTW_FAST);
        return forward[index];
      }
    }

    fftw_plan& Backward(size_t size, std::complex<double>* in, double* out)
    {
      // Check power of two
      assert(size && ((size & (size-1))) ==0);
      size_t index = (size_t)log2(size);
      if (backward[index]) {
        return backward[index];
      } else {
        std::lock_guard<std::mutex> guard(g_plan_mutex);
        backward[index] = fftw_plan_dft_c2r_1d((int)size, reinterpret_cast<fftw_complex*>(in), out, SPS_FFTW_FAST);
        return backward[index];
      }
    }

  private:
    Signal1DPlan() {}
    ~Signal1DPlan()
    {
      for (size_t i = 0 ; i < nFFTLengths ; i++) {
        if (forward[i]) {
          fftw_destroy_plan(forward[i]);
          forward[i] = NULL;
        }
        if (backward[i]) {
          fftw_destroy_plan(backward[i]);
          backward[i] = NULL;
        }
      }
    }
    Signal1DPlan(Signal1DPlan const&) = default;
    Signal1DPlan& operator=(Signal1DPlan const&) = default;
  };

  __THREAD fftw_plan Signal1DPlan<double>::forward[Signal1DPlan<double>::nFFTLengths] = {NULL};
  __THREAD fftw_plan Signal1DPlan<double>::backward[Signal1DPlan<double>::nFFTLengths] = {NULL};






  template class Signal1DPlan<float>;
  template class Signal1DPlan<double>;

  template <typename T>
  STATIC_INLINE_BEGIN T* _mm_padarray(const T* input, const size_t len, const size_t newlen)
  {
    T* output = (T*) _mm_malloc(newlen*sizeof(T),16);
#if 0
    memset(output, 0, newlen*sizeof(T));
    memcpy(output, (void*)input, len*sizeof(T));
#else
    memcpy(output, (void*)input, len*sizeof(T));
    memset(&output[len], 0, (newlen-len)*sizeof(T));
#endif
    return output;
  }

  template <typename T>
  void DivideArray(T *Data, size_t NumEl, T Divisor)
  {
    size_t n;
    for(n = 0; n < NumEl; n++)
      Data[n] /= Divisor;
  }

  template <>
  void DivideArray<float>(float *Data, size_t NumEl, float Divisor)
  {

    // TODO: Consider using _mm_loadu_ps for the first (if any)

    const size_t align = 16;
    const uintptr_t mask = ~(uintptr_t)(align - 1);

    float *aligned_ptr = (float *)(((uintptr_t)Data+align-1) & mask); // mask is 0x0F

    // Offset before first aligned address
    size_t offset = (size_t) ((uintptr_t)aligned_ptr - (uintptr_t) Data) / sizeof(float);

    const float inv_Divisor = 1.0f / Divisor;

    // Number of unaligned values before alignment (if any)
    int nCount = (int) std::min<size_t>(NumEl,offset);

    int n = 0;
    for (n=0 ; n < nCount ; n++) {
      Data[n] *= inv_Divisor;
    }

    const __m128 divisor = _mm_set1_ps(inv_Divisor);

    // Number of aligned values after first aligned address (TODO: Verify)
    nCount = int( NumEl - (align/sizeof(float) - 1) - offset);

    // Shouldn't it be
    nCount = int((align/sizeof(float)) * ((NumEl - offset) / (align / sizeof(float))));

    // Only works for alignment equal to 16 (SSE2,SSE4)
    for (; n < nCount ; n += align/sizeof(float)) {
      __m128 va = _mm_load_ps(&aligned_ptr[n]);
      va = _mm_mul_ps(va, divisor);
      _mm_store_ps(&aligned_ptr[n],va);
    }

    // Remaining unaligned values (if any)
    for(; n < (int)NumEl; n++)
      Data[n] *= inv_Divisor;
  }

  template <>
  void DivideArray<double>(double *Data, size_t NumEl, double Divisor)
  {

    const size_t align = 16;
    const uintptr_t mask = ~(uintptr_t)(align - 1);

    double *aligned_ptr = (double *)(((uintptr_t)Data+align-1) & mask);

    int n = 0;

    size_t offset = (size_t) ((uintptr_t)aligned_ptr - (uintptr_t) Data) / sizeof(double);

    const double inv_Divisor = 1.0 / Divisor;

    // Number of unaligned values before alignment (if any)
    int nCount = int(std::min<size_t>(NumEl,offset));

    // Shouldn't it be
    nCount = int((align/sizeof(double)) * ((NumEl - offset) / (align / sizeof(double))));

    for (n=0 ; n < nCount ; n++) {
      Data[n] *= inv_Divisor;
    }

    const __m128d divisor = _mm_set1_pd(inv_Divisor);

    // Number of aligned values after first aligned address
    nCount = (int) NumEl - (align/sizeof(double) - 1) - (int) offset;

    // Divide aligned values (continue from above)
    for (; n < nCount ; n += align/sizeof(double)) {
      __m128d va = _mm_load_pd(&aligned_ptr[n]);
      va = _mm_mul_pd(va, divisor);
      _mm_store_pd(&aligned_ptr[n],va);
    }

    // Remaining unaligned values (if any)
    for(; n < (int)NumEl; n++)
      Data[n] *= inv_Divisor;
  }

  template <class T>
  signal1D<T>::signal1D() : data(NULL), offset(0), ndata(0), nbytes(0) {}

  template <class T>
  signal1D<T>::signal1D(size_t ndata) : data(NULL), offset(0), ndata(ndata)
  {
    nbytes = 16 * (ndata * sizeof(T) + 15) / 16;
    data = (T*) _mm_malloc(nbytes,16);
#ifndef NDEBUG
    memset(data,0,nbytes);
#endif
  }

  template <class T>
  signal1D<T>::signal1D(size_t ndata, size_t _nbytes) : data(NULL), offset(0), ndata(ndata), nbytes(0)
  {
    nbytes = 16 * (ndata * sizeof(T) + 15) / 16;
    nbytes = std::max<size_t>(_nbytes,nbytes);
    data = (T*) _mm_malloc(nbytes,16);
#ifndef NDEBUG
    memset(data,0,nbytes);
#endif
  }

  template <class T>
  signal1D<T>::~signal1D()
  {
    if (data) {
      _mm_free(data);
      data = NULL;
    }
  }


#if defined(__GNUG__) || (defined(_MSC_VER) && (_MSC_VER >= 1800))

  template <typename T>
  bool pack_r2c(const sps::signal1D<T> &a, sps::signal1D<std::complex<T> >& c)
  {

    size_t nbytes = 16 * (a.ndata * sizeof(T) + 15) / 16;
    c.ndata = a.ndata / 2;

    if (c.nbytes < nbytes) {
      c.nbytes = nbytes;
      if (c.data) {
        _mm_free(c.data);
        c.data = NULL;
      }
      c.data = (std::complex<T>*) _mm_malloc(c.nbytes, 16);
    }
    memcpy(c.data,a.data,a.ndata * sizeof(T));
    return true;
  }

  template <typename T>
  bool unpack_c2r(const sps::signal1D<std::complex<T> >& a, sps::signal1D<T> &c)
  {

    size_t nbytes = 16 * (a.ndata * sizeof(std::complex<T>) + 15) / 16;
    c.ndata = 2*a.ndata;

    if (c.nbytes < nbytes) {
      c.nbytes = nbytes;
      if (c.data) {
        _mm_free(c.data);
        c.data = NULL;
      }
      c.data = (T*) _mm_malloc(c.nbytes, 16);
    }
    memcpy(c.data, a.data, a.ndata * sizeof(std::complex<T>));
    return true;
  }

  template <> bool
  fft<float>(const sps::signal1D<float> &a, const size_t &n, sps::signal1D<std::complex<float> >& c)
  {

    c.ndata = n/2 + 1;
    c.offset = 0;

    // TODO: Make 2 a divisor (pair-wise multiplication)
    //            4 a divisor (four-wise multiplication - filters)

    // One extra sample for pair-wise multiplication of complex data
    size_t nbytes = 16 * ((c.ndata+1)*sizeof(std::complex<float>) + 15) / 16;

    bool retval = true, pad_a = false;

    float* _a = NULL;

    while(retval) {
      // Allocate if necessary
      if (c.data) {
        if (c.nbytes < nbytes) {
          _mm_free(c.data);
          c.nbytes = nbytes;
          c.data = (std::complex<float>*) _mm_malloc(c.nbytes,16);
        }
      } else {
        c.nbytes = nbytes;
        c.data = (std::complex<float>*) _mm_malloc(c.nbytes,16);
      }

      // Pad array if necessary
      pad_a = a.nbytes < n * sizeof(float);
      _a = pad_a ? _mm_padarray<float>(a.data, a.ndata, n) : a.data;
      assert(((uintptr_t)_a & 0x0F) == 0 && "Data must be aligned");

      // No need for zeroing output (TEST - NOW WE ZERO OUTPUT)
      //memset(c.data,0,c.ndata * sizeof(std::complex<float>));

      Signal1DPlan<float>& p = Signal1DPlan<float>::Instance();
      fftwf_execute_dft_r2c(p.Forward(n, _a, c.data), _a, reinterpret_cast<fftwf_complex*>(c.data));
      break;
    }
    return true;
  }

  template <> bool
  ifft<float>(const sps::signal1D<std::complex<float> >& a, const size_t &n, sps::signal1D<float> &c)
  {

    c.ndata = n;
    c.offset = 0;

    size_t nbytes = 16 * (n*sizeof(float) + 15) / 16;

    bool retval = true, pad_a = false;

    std::complex<float>* _a = NULL;

    while(retval) {
      // Allocate if necessary
      if (c.data) {
        if (c.nbytes != nbytes) {
          _mm_free(c.data);
          c.nbytes = nbytes;
          c.data = (float*) _mm_malloc(c.nbytes,16);
        }
      } else {
        c.nbytes = nbytes;
        c.data = (float*) _mm_malloc(c.nbytes,16);
      }

      // Pad array if necessary
      pad_a = a.nbytes < (n/2 + 1) * sizeof(std::complex<float>);
      _a = pad_a ? _mm_padarray<std::complex<float> >(a.data, a.ndata, n/2 + 1) : a.data;
      assert(((uintptr_t)_a & 0x0F) == 0 && "Data must be aligned");

      // Not necessary
      // memset(c.data,0,sizeof(float) * n);

      Signal1DPlan<float>& p = Signal1DPlan<float>::Instance();
      fftwf_execute_dft_c2r(p.Backward(n,_a,c.data), reinterpret_cast<fftwf_complex*>(_a), c.data);

      DivideArray<float>(c.data, n, float(n));

      break;
    }
    return true;
  }

  template <> bool
  fft<double>(const sps::signal1D<double> &a, const size_t &n, sps::signal1D<std::complex<double> >& c)
  {

    c.ndata = n/2 + 1;
    c.offset = 0;

    // One extra sample for pair-wise multiplication of complex data
    size_t nbytes = 16 * ((c.ndata+1)*sizeof(std::complex<double>) + 15) / 16;

    bool retval = true, pad_a = false;

    double* _a = NULL;

    while(retval) {
      if (c.data) {
        if (c.nbytes != nbytes) {
          _mm_free(c.data);
          c.nbytes = nbytes;
          c.data = (std::complex<double>*) _mm_malloc(c.nbytes,16);
        }
      } else {
        c.nbytes = nbytes;
        c.data = (std::complex<double>*) _mm_malloc(c.nbytes,16);
      }

      pad_a = a.nbytes < n * sizeof(double);
      _a = pad_a ? _mm_padarray<double>(a.data, a.ndata, n) : a.data;
      assert(((uintptr_t)_a & 0x0F) == 0 && "Data must be aligned");

      // Not necessary
      //      memset(c.data,0,c.ndata * sizeof(std::complex<double>));

      Signal1DPlan<double>& p = Signal1DPlan<double>::Instance();
      fftw_execute_dft_r2c(p.Forward(n,_a,c.data), _a, reinterpret_cast<fftw_complex*>(c.data));
      break;
    }
    return true;
  }

  template <> bool
  ifft<double>(const sps::signal1D<std::complex<double> >& a, const size_t &n, sps::signal1D<double> &c)
  {

    c.ndata = n;
    c.offset = 0;

    size_t nbytes = 16 * (n*sizeof(double) + 15) / 16;

    bool retval = true, pad_a = false;

    std::complex<double>* _a = NULL;

    while(retval) {
      if (c.data) {
        if (c.nbytes != nbytes) {
          _mm_free(c.data);
          c.nbytes = nbytes;
          c.data = (double*) _mm_malloc(c.nbytes,16);
        }
      } else {
        c.nbytes = nbytes;
        c.data = (double*) _mm_malloc(c.nbytes,16);
      }

      pad_a = a.nbytes < (n/2 + 1) * sizeof(std::complex<double>);
      _a = pad_a ? _mm_padarray<std::complex<double> >(a.data, a.ndata, n/2 + 1) : a.data;
      assert(((uintptr_t)_a & 0x0F) == 0 && "Data must be aligned");

      // Not necessary
      //      memset(c.data,0,sizeof(double) * n);

      Signal1DPlan<double>& p = Signal1DPlan<double>::Instance();
      fftw_execute_dft_c2r(p.Backward(n,_a,c.data), reinterpret_cast<fftw_complex*>(_a), c.data);

      DivideArray<double>(c.data, n, double(n));

      break;
    }
    return true;
  }

#endif

  template <> bool
  conv_fft<float>(const sps::signal1D<float> &a,
                  const sps::signal1D<float> &b,
                  sps::signal1D<float>& c)
  {

    float *_a = NULL, *_b = NULL;
    fftwf_complex *fft_a = NULL, *fft_b = NULL;

    bool pad_a = false, pad_b = false;

    size_t n_a = a.ndata;
    size_t n_b = b.ndata;
    size_t n = next_power_two<size_t>(n_a+n_b-1);

    debug_print("na: %zu, nb: %zu, n: %zu\n",n_a,n_b,n);

    bool retval = true;

    if ((!a.data || !b.data ) || (a.ndata==0 || b.ndata==0)) {
      debug_print("Uninitialized data\n");
      retval = false;
    }
    while (retval) {

      // If the arrays are long enough, we don't need to copy data
      pad_a = a.nbytes < n * sizeof(float);
      // Ups
      _a = pad_a ? _mm_padarray<float>(a.data, n_a, n) : a.data;
      assert(((uintptr_t)_a & 0x0F) == 0 && "Data must be aligned");

      pad_b = b.nbytes < n * sizeof(float);
      _b = pad_b ? _mm_padarray<float>(b.data, n_b, n) : b.data;
      assert(((uintptr_t)_b & 0x0F) == 0 && "Data must be aligned");

      // One extra complex point added to use fast complex multiply (we multiply two at a time)
      fft_a = (fftwf_complex*) _mm_malloc(sizeof(fftwf_complex) * ((n/2+1) + 1),16);
      fft_b = (fftwf_complex*) _mm_malloc(sizeof(fftwf_complex) * ((n/2+1) + 1),16);

      // TEST TO ZERO ARRAYS
      //memset(fft_a, 0, sizeof(fftwf_complex) * ((n/2+1) + 1));
      //memset(fft_b, 0, sizeof(fftwf_complex) * ((n/2+1) + 1));

      // No need to zero output array
      c.ndata = n_a+n_b-1;
      c.offset = a.offset + b.offset;

      size_t nbytes = 16 * (n*sizeof(float) + 15) / 16; // Added 16 for complex multiply (removed)

      // IS IT AN ISSUE THAT WE REALLOC OUTPUT
      if (c.data) {
        assert(((uintptr_t)c.data & 0x0F) == 0 && "Data must be aligned");
        if (c.nbytes < nbytes) { // Was != nbytes
          _mm_free(c.data);
          c.nbytes = nbytes;
          c.data = (float*) _mm_malloc(c.nbytes,16);
        }
      } else {
        c.nbytes = nbytes;
        c.data = (float*) _mm_malloc(c.nbytes,16);
      }

      memset(c.data,0,c.nbytes);
      Signal1DPlan<float>& p = Signal1DPlan<float>::Instance();

      fftwf_plan forward = p.Forward(n,_a,(std::complex<float>*)fft_a);
      fftwf_plan backward= p.Backward(n,(std::complex<float>*)fft_a,c.data);
      fftwf_execute_dft_r2c(forward, _a, fft_a);
      fftwf_execute_dft_r2c(forward, _b, fft_b); // Error here???

      //#define CONV_FFT_USE_COMPLEX 1

#ifdef CONV_FFT_USE_COMPLEX
      // It is incredible slow
      std::complex<float> * cfft_a = reinterpret_cast<std::complex<float>* >(fft_a);
      std::complex<float> * cfft_b = reinterpret_cast<std::complex<float>* >(fft_b);

      for (size_t i = 0 ; i < (n/2+1) ; i++) {
        cfft_a[i] = cfft_a[i] * cfft_b[i] / float(n);
      }
      fftwf_execute_dft_c2r(backward, fft_a, c.data);

#elif defined(CONV_FFT_USE_SCALARS)
      for (size_t i = 0 ; i < (n/2+1) ; i++) {
        // For some reason this is way faster than using <complex>
        float real = (fft_a[i][0]*fft_b[i][0] - fft_a[i][1]*fft_b[i][1]);
        float imag = (fft_a[i][0]*fft_b[i][1] + fft_a[i][1]*fft_b[i][0]);
        fft_a[i][0] = real;
        fft_a[i][1] = imag;
      }
      fftwf_execute_dft_c2r(backward, fft_a, c.data);

      // Why n-1
      DivideArray<float>(c.data, n-1, float(n));

#else
      __m128 divisor = _mm_set1_ps(1.0f / float(n));

      for (size_t i = 0 ; i < (n/2+1)+1 ; i+=2) {
        __m128 vec_a = _mm_load_ps((float*)&fft_a[i]);
        __m128 vec_b = _mm_load_ps((float*)&fft_b[i]);
        vec_b = _mm_mul_ps(
              _mm_mulcmplx_ps(vec_a,vec_b),
              divisor);
        _mm_store_ps((float*)&fft_a[i][0],vec_b);
      }
#endif
      fftwf_execute_dft_c2r(backward, fft_a, c.data);
      break;
    }

    if (pad_a)
      _mm_free(_a);
    _a = NULL;

    if (pad_b)
      _mm_free(_b);
    _b = NULL;

    if (fft_a)
      _mm_free(fft_a);
    fft_a = NULL;

    if (fft_b)
      _mm_free(fft_b);
    fft_b = NULL;

    return retval;
  }

  template <> bool
  conv_fft<double>(const sps::signal1D<double> &a,
                   const sps::signal1D<double> &b,
                   sps::signal1D<double>& c)
  {

    double *_a = NULL, *_b = NULL;

    fftw_complex *fft_a = NULL, *fft_b = NULL;

    bool pad_a = false, pad_b = false;

    size_t n_a = a.ndata;
    size_t n_b = b.ndata;
    size_t n = next_power_two<size_t>(n_a+n_b-1);

    bool retval = true;

    if (!a.data || !b.data )
      retval = false;

    while (retval) {
      pad_a = a.nbytes < n * sizeof(double);
      _a = pad_a ? _mm_padarray<double>(a.data, n_a, n) : a.data;

      pad_b = b.nbytes < n * sizeof(double);
      _b = pad_b ? _mm_padarray<double>(b.data, n_b, n) : b.data;

      // Here FFTW is used (no need to memset)
      fft_a = (fftw_complex*) _mm_malloc(sizeof(fftw_complex) * (n/2+1),16);
      fft_b = (fftw_complex*) _mm_malloc(sizeof(fftw_complex) * (n/2+1),16);

      size_t nbytes = 16 * (n*sizeof(double) + 15) / 16;

      c.ndata = n_a+n_b-1;
      c.offset = a.offset + b.offset;

      // Allocate if needed
      if (c.data) {
        if (c.nbytes != nbytes) {
          _mm_free(c.data);
          c.nbytes = nbytes;
          c.data = (double*) _mm_malloc(c.nbytes,16);
        }
      } else {
        c.nbytes = nbytes;
        c.data = (double*) _mm_malloc(c.nbytes,16);
      }

      Signal1DPlan<double>& p = Signal1DPlan<double>::Instance();
      fftw_plan forward = p.Forward(n,_a,(std::complex<double>*)fft_a);
      fftw_plan backward= p.Backward(n,(std::complex<double>*)fft_a, c.data);

      fftw_execute_dft_r2c(forward, _a, fft_a);
      fftw_execute_dft_r2c(forward, _b, fft_b);

      // Complex multiplication
      for (size_t i = 0 ; i < (n/2+1) ; i++) {
        double real = fft_a[i][0]*fft_b[i][0] - fft_a[i][1]*fft_b[i][1];
        double imag = fft_a[i][0]*fft_b[i][1] + fft_a[i][1]*fft_b[i][0];
        fft_a[i][0] = real;
        fft_a[i][1] = imag;
      }

      fftw_execute_dft_c2r(backward, fft_a, c.data);

      // Why n-1
      DivideArray<double>(c.data, n-1, double(n));
      break;
    }
    if (pad_a)
      _mm_free(_a);
    _a = NULL;
    if (pad_b)
      _mm_free(_b);
    _b = NULL;

    if (fft_a)
      _mm_free(fft_a);
    fft_a = NULL;
    if (fft_b)
      _mm_free(fft_b);
    fft_b = NULL;

    return retval;
  }

  template <typename T>
  bool conv_fft_fs(const T& fs,
                   const sps::signal1D<T> &a,
                   const sps::signal1D<T> &b,
                   sps::signal1D<T>& c)
  {

    const T invfs = T(1.0) / fs;
    bool retval = conv_fft(a,b,c);

    for (size_t i = 0; i < c.ndata ; i++) {
      c.data[i] = c.data[i] * invfs;
    }

    return retval;
  }

  template <class T>
  msignal1D<T>::msignal1D() : offset(0), ndata(0), nbytes(0), m_data(NULL) {}

  template <class T>
  msignal1D<T>::msignal1D(size_t ndata) : msignal1D(ndata, 16 * (ndata * sizeof(T) + 15) / 16) {}

  template <class T>
  msignal1D<T>& msignal1D<T>::operator=(const msignal1D& a)
  {
    this->m_data = a.m_data;
    this->offset = a.offset;
    this->ndata  = a.ndata;
    this->nbytes = a.nbytes;
    return *this;
  }

  template <class T>
  msignal1D<T>::msignal1D(const sps::msignal1D<T>& a)
  {
    this->m_data = a.m_data;
    this->offset = a.offset;
    this->ndata  = a.ndata;
    this->nbytes = a.nbytes;
  }

  template <class T>
  msignal1D<T>::msignal1D(size_t _ndata, size_t _nbytes) : offset(0),
    ndata(_ndata),
    nbytes(16 * ( std::max<size_t>(_ndata * sizeof(T),_nbytes) + 15 ) / 16), m_data(NULL)
  {
    m_data.reset( (T*) _mm_malloc(nbytes,16), [=](T *p) {
      _mm_free(p);
      p=NULL;
    });
#ifndef NDEBUG
    memset(m_data.get(),0,nbytes);
#endif
  }
  /*
  template <class T>
  const T& msignal1D<T>::operator[](const size_t index) const
  {
    return this->m_data.get()[index];
  }

  template <class T>
  inline T& msignal1D<T>::operator[](const size_t index)
  {
    return this->m_data.get()[index];
  }
  */

  template <class T>
  void msignal1D<T>::reset(size_t _ndata, size_t _nbytes)
  {
    // Maximum of current and number of elements
    _ndata  = std::max<size_t>(_ndata,ndata);

    // Maximum of current and new size
    _nbytes = std::max<size_t>(_nbytes, _ndata*sizeof(T));

    // Re-alloc if needed
    if (_nbytes > nbytes) {
      nbytes = _nbytes;
      m_data.reset( (T*) _mm_malloc(nbytes,16), [=](T *p) {
        _mm_free(p);
        p=NULL;
      });
    }
    // Reset data
    ndata = _ndata;
    memset(m_data.get(),0,ndata*sizeof(T));
  }
  template <class T>
  void msignal1D<T>::reverse()
  {
    std::shared_ptr<T> data1 = m_data;
    m_data.reset( (T*) _mm_malloc(nbytes,16), [=](T *p) {
      _mm_free(p);
      p=NULL;
    });
    for (size_t i = 0 ; i < ndata ; i++) {
      m_data.get()[i] = data1.get()[ndata-1-i];
    }
  }

  template <class T>
  void msignal1D<T>::pad(size_t _nbytes)
  {
    if ((nbytes > 0) && (_nbytes > nbytes)) {
      // Shallow copy
      std::shared_ptr<T> data1 = m_data;
      m_data.reset( (T*) _mm_malloc(_nbytes,16), [=](T *p) {
        _mm_free(p);
        p=NULL;
      });
      memset(m_data.get(),0,_nbytes);
      memcpy(m_data.get(), data1.get(), ndata*sizeof(T));
      nbytes = _nbytes;
    }
  }

  template <class T>
  msignal1D<T>::~msignal1D()
  {
    if (ndata) {
      if (m_data) {
        m_data.reset((T*)NULL);
      }
    }
  }

  template <class T>
  void msignal1D<T>::scale(const T& s)
  {
    for (size_t i = 0 ; i < ndata ; i++) {
      m_data.get()[i] *= s;
    }
  }

  template <>
  void msignal1D<float>::scale(const float& s)
  {
    __m128 multiplier = _mm_set1_ps(s);
    for (size_t i = 0 ; i < ndata ; i += 4) {
      __m128 va = _mm_load_ps(&(m_data.get()[i]));
      va = _mm_mul_ps(va, multiplier);
      _mm_store_ps(&(m_data.get()[i]),va);
    }
  }


// TODO: Overwrites output if length is 1 (optimization)
  template<typename T>
  bool mconv_fft_fs(const T& fs,
                    const msignal1D<T>& a, const msignal1D<T>& b, msignal1D<T>& c)
  {

    const T invfs = T(1.0) / fs;
    bool retval = mconv_fft(a,b,c);
    retval = true;
    T val;

    // Always scale using 1/fs
    if (retval) {
      for (size_t i = 0; i < c.ndata ; i++) {
        val = c.m_data.get()[i] * invfs;
        c.m_data.get()[i] = val;
      }
    }

    return retval;
  }

  template <typename T>
  bool conv_fft(const sps::signal1D<T> &a,
                sps::signal1D<T> &b)
  {
    // In-place of second argument
    return true;
  }

  template <>
  bool conv_fft_out_in(const signal1D<float> &a,
                       signal1D<float> &b)
  {

    float *_a = NULL, *_b = NULL;
    fftwf_complex *fft_a = NULL;

    bool pad_a = false;

    size_t n_a = a.ndata;
    size_t n_b = b.ndata;
    size_t n = next_power_two<size_t>(n_a+n_b-1);

    bool retval = true;

    if (!a.data || !b.data )
      retval = false;

    // Internal samples needed = n/2+1 complex samples + one extra to make an even number, we use SIMD
    size_t nInternal = 2 * (n/2 + 2);

    while (retval) {
      pad_a = a.nbytes < n * sizeof(float);
      _a = pad_a ? _mm_padarray(a.data, n_a, n) : a.data;

      // In-place
      if (b.nbytes < nInternal * sizeof(float)) {
        // We need to reallocate
        float* pTmp = b.data;
        b.nbytes = 16 * (nInternal * sizeof(float) + 15) / 16;
        b.data = (float*) _mm_malloc(b.nbytes,16);
        memset(b.data,0,nInternal*sizeof(float));
        memcpy(b.data,(void*)pTmp,n_b*sizeof(float));
        _mm_free(pTmp);
      }
      _b =  b.data;

      // One extra point added to use fast complex multiply (we multiply two at a time)
      fft_a = (fftwf_complex*) _mm_malloc(sizeof(fftwf_complex) * ((n/2+1) + 1),16);

      b.ndata = n_a+n_b-1;
      b.offset = a.offset + b.offset;

      Signal1DPlan<float>& p = Signal1DPlan<float>::Instance();

      fftwf_plan forward = p.Forward(n,_a, (std::complex<float>*)fft_a);
      fftwf_plan backward = p.Backward(n,(std::complex<float>*)_b,_b);

      fftwf_execute_dft_r2c(forward, _a, fft_a);
      fftwf_execute_dft_r2c(forward, _b, (fftwf_complex*)_b);

      __m128 divisor = _mm_set1_ps(1.0f / float(n));

      for (size_t i = 0 ; i < (n/2+1)+1 ; i+=2) {
        __m128 vec_a = _mm_load_ps((float*)&fft_a[i]);
        __m128 vec_b = _mm_load_ps((float*)&_b[2*i]);
		vec_b = _mm_mul_ps(_mm_mulcmplx_ps(vec_a, vec_b),divisor);
        _mm_store_ps((float*)&_b[2*i], vec_b);
      }

      fftwf_execute_dft_c2r(backward, (fftwf_complex*)_b, _b);
      break;
    }
    if (pad_a)
      _mm_free(_a);
    _a = NULL;

    if (fft_a) {
      _mm_free(fft_a);
      fft_a = NULL;
    }

    return retval;

  }

  template <>
  bool conv_fft_out_in(const signal1D<double> &a,
                       signal1D<double> &b)
  {
    SPS_UNREFERENCED_PARAMETERS(a,b);
    assert(false && "Not implemented");
    return false;
  }

  template <class T>
  void conv_ansi(const T* Signal, size_t SignalLen,
                 const T* Kernel, size_t KernelLen,
                 T* Result)
  {
    size_t n;

    for (n = 0; n < SignalLen + KernelLen - 1; n++) {
      size_t kmin, kmax, k;

      Result[n] = 0;

      kmin = (n >= KernelLen - 1) ? n - (KernelLen - 1) : 0;
      kmax = (n < SignalLen - 1) ? n : SignalLen - 1;

      for (k = kmin; k <= kmax; k++) {
        Result[n] += Signal[k] * Kernel[n - k];
      }
    }
  }

  template <typename T>
  bool conv(const signal1D<T> &a,
            const signal1D<T> &b,
            signal1D<T>& c)
  {

    // check validity of params
    if ((a.data == NULL) || (b.data==NULL)) {
      return false;
    }

    int na = (int)a.ndata;
    int nb = (int)b.ndata;

    T* _a = a.data;
    T* _b = b.data;

    size_t nbytes = 16 * ((na+nb-1)*sizeof(T) + 15) / 16;

    // Allocate if needed
    if (c.data) {
      if (c.nbytes != nbytes) {
        _mm_free(c.data);
        c.nbytes = nbytes;
        c.data = (T*) _mm_malloc(c.nbytes,16);
      }
    } else {
      c.nbytes = nbytes;
      c.data = (T*) _mm_malloc(c.nbytes,16);
    }

    memset(c.data,0,(na+nb-1)*sizeof(T));
    c.ndata = na+nb-1;
    T* out = c.data;

#if 0
    // Works if nb > na
    int i, j, k;

    // convolution from out[0] to out[nb-2]
    for(i = 0; i < nb - 1; ++i) {
      out[i] = T(0);
      for(j = 0, k = i ; j <= i; ++j, --k)
        out[i] += _a[j] * _b[k];
    }

    // convolution from out[nb-1] to out[na-1]
    for(i = nb-1; i < na; ++i) {
      out[i] = T(0);
      for(j = i, k = 0; k < nb; --j, ++k) // AV (trouble)
        out[i] += _a[j] * _b[k];
    }

    // convolution from out[na] to out[na + nb-2] (last)
    for(i = 0; i < nb - 1; ++i) {
      out[na+i] = T(0);
      for(j = 0, k = i + 1; j < (nb - 1)-i; j++,k++) {
        // i=0,j=0,k=1 (issue)
        if (i+j > 0) {
          out[na+i] += _a[na-nb+k] * _b[nb-1-j];
        }
      }
    }
#else
    conv_ansi(_a, na, _b, nb, out);
#endif
    c.offset = a.offset + b.offset;

    return true;
  }

  template<typename T>
  bool mconv_fft(const msignal1D<T>& a, const msignal1D<T>& b, signal1D<T>& c)
  {
    bool retval = true;

    size_t na = a.ndata;
    size_t nb = b.ndata;

    bool bMulA = false;
    bool bMulB = false;

    T scale = T(1.0);

    // < 3, we have an issue
    if (na+nb == 0) {
      return false;
    }

    debug_print("na: %zu, nb: %zu\n",na,nb);
    if (na < 2) {
      if (nb < 2) {
        // A is multiplied by b[0]
        bMulA = true;
      } else {
        bMulB = true;
      }
    } else if (nb < 2) {
      bMulA = true;
    }

    if (bMulA && bMulB) {
      // What about bMulA and bMulB !!!
      exit(-1);
    }

    if (bMulA) {
      if (c.data) {
        if (c.nbytes < a.ndata*sizeof(T)) {
          _mm_free(c.data);
          c.data = (T*) _mm_malloc(a.ndata*sizeof(T),16);
        }
      } else {
        c.data = (T*) _mm_malloc(a.ndata*sizeof(T),16);
      }
      if (nb == 1) {
        scale = b.m_data.get()[0];
      }

      for (size_t i = 0 ; i < a.ndata ; i++) {
        c.data[i] = a.m_data.get()[i] * scale;
      }
      return retval;
    }

    if (bMulB) {
      if (c.data) {
        if (c.nbytes < b.ndata*sizeof(T)) {
          _mm_free(c.data);
          c.data = (T*) _mm_malloc(b.ndata*sizeof(T),16);
        }
      } else {
        c.data = (T*) _mm_malloc(b.ndata*sizeof(T),16);
      }
      if (na == 1) {
        scale = a.m_data.get()[0];
      }

      debug_print("multiplier: %f\n", scale);

      for (size_t i = 0 ; i < b.ndata ; i++) {
        c.data[i] = b.m_data.get()[i] * scale;
      }
      return retval;
    }

    // Unmanaged signals (temporary)
    signal1D<T> _a, _b;
    _a.data   = a.m_data.get();
    _a.ndata  = a.ndata;
    _a.offset = a.offset;
    _a.nbytes = a.nbytes;

    _b.data   = b.m_data.get();
    _b.ndata  = b.ndata;
    _b.offset = b.offset;
    _b.nbytes = b.nbytes;

    // Convolve the signals
    retval = conv_fft<T>(_a,_b,c);

    // Prevent destruction of data
    _a.data = NULL;
    _b.data = NULL;

    return retval;
  }

// WRONG

// Make scaling more elegant
  template<typename T>
  bool mconv_fft(const msignal1D<T>& a, const msignal1D<T>& b, msignal1D<T>& c)
  {

    bool retval = true;

    size_t na = a.ndata;
    size_t nb = b.ndata;

    if (na+nb == 0) {
      // Consider exiting if any is zero
      return false;
    }

#if 1
    // [1] * [] = []

    size_t n = next_power_two<size_t>(na+nb-1);

    size_t nbytes = 16 * (n*sizeof(T) + 15) / 16;

    if (c.nbytes < nbytes) {
      c.nbytes = nbytes;
      c.m_data = std::shared_ptr<T>( (T*) _mm_malloc(c.nbytes,16), [=](T *p) {
        _mm_free(p);
        p = NULL;
      });
    }
    c.ndata  = na + nb - 1;
    c.offset = a.offset + b.offset;

#else

    // FIX THIS (UGLY) ask if na and nb > 1 else if na==1 or nb==1 or both
    if (na < 2) {
      if (nb < 2) {
        if (b.ndata == 1) {
          if (b.data.get()[0] == T(1.0)) {
            // TODO: Copy values
            c = a;
            return retval;
          } else {
            if (a.ndata == 1) {
              // Scale
              c.ndata  = na + nb - 1;
              c.nbytes = 16;
              c.offset = a.offset + b.offset;
              c.data = std::shared_ptr<T>( (T*) _mm_malloc(c.nbytes,16), [=](T *p) {
                _mm_free(p);
                p = NULL;
              });
              c.data.get()[0] = a.data.get()[0] * b.data.get()[0];
              return retval;
            } else {
              c = 1;
              return retval;
            }
          }
        } else {
          c = a;
          return retval;
        }
      } else {
        // Scale
        if (a.ndata == 1) {
          c.ndata  = na + nb - 1;
          c.nbytes = 16 * (c.ndata*sizeof(T) + 15) / 16;
          c.offset = a.offset + b.offset;
          c.data = std::shared_ptr<T>( (T*) _mm_malloc(c.nbytes,16), [=](T *p) {
            _mm_free(p);
            p = NULL;
          });
          for (size_t i = 0 ; i < b.ndata ; i++) {
            c.data.get()[i] = b.data.get()[i] * a.data.get()[0];
          }
          return retval;
        } else {
          c = b;
          return retval;
        }
      }
    } else if (nb < 2) {
      // TODO: Scale
      c = a;
      return retval;
    }

    // FFT points
    size_t n = next_power_two<size_t>(na+nb-1);

    size_t nbytes = 16 * (n*sizeof(T) + 15) / 16; // Added +16 for extra space complex multiply

    // TODO: TODO: TODO: Make sure that no memory is allocated by conv_fft for output signal
    c.ndata  = na + nb - 1;
    c.nbytes = nbytes;
    c.offset = a.offset + b.offset;
    c.data = std::shared_ptr<T>( (T*) _mm_malloc(c.nbytes,16), [=](T *p) {
      _mm_free(p);
      p = NULL;
    });

#endif

    // Unmanaged signals (temporary)
    signal1D<T> _a, _b, _c;
    _a.data   = a.m_data.get();
    _a.ndata  = a.ndata;
    _a.offset = a.offset;
    _a.nbytes = a.nbytes;

    _b.data   = b.m_data.get();
    _b.ndata  = b.ndata;
    _b.offset = b.offset;
    _b.nbytes = b.nbytes;

    _c.data   = c.m_data.get();
    _c.ndata  = c.ndata;
    _c.offset = c.offset;
    _c.nbytes = c.nbytes;

    // Convolve the signals
    retval = conv_fft<T>(_a,_b,_c); // Error - when not enough memory

    // Prevent destruction of data
    _a.data = NULL;
    _b.data = NULL;
    _c.data = NULL;

    return retval;
  }

  template<typename T>
  bool mfft(const msignal1D<T>& a, const size_t& n, msignal1D<std::complex<T> >& c)
  {

    bool retval = true;

    c.ndata  = n/2 + 1;
    c.offset = 0;
    // We often multiple two complex numbers at a time (consider c.ndata+1
    size_t nbytes = 16 * ( (c.ndata+1)*sizeof(std::complex<T>) + 15) / 16;

    if (c.nbytes < nbytes) {
      c.nbytes = nbytes;
      c.m_data = std::shared_ptr<std::complex<T>>( (std::complex<T>*) _mm_malloc(c.nbytes,16), [=](std::complex<T> *p) {
        _mm_free(p);
        p = NULL;
      });
    }

    // Unmanaged signals (temporary)
    signal1D<T> _a;
    _a.data   = a.m_data.get();
    _a.ndata  = a.ndata;
    _a.offset = a.offset;
    _a.nbytes = a.nbytes;

    signal1D<std::complex<T> > _c;
    _c.data   = c.m_data.get();
    _c.ndata  = c.ndata;
    _c.offset = c.offset;
    _c.nbytes = c.nbytes;

    retval = fft<T>(_a, n, _c);

    // Prevent destruction of data
    _a.data = NULL;
    _c.data = NULL;

    return retval;
  }

  template<typename T>
  bool mifft(const msignal1D<std::complex<T> >& a, const size_t &n, msignal1D<T>& c)
  {

    bool retval = true;

    c.ndata  = n;
    c.offset = 0;
    size_t nbytes = 16 * (c.ndata*sizeof(T) + 15) / 16;

    if (c.nbytes < nbytes) {
      c.nbytes = nbytes;
      c.m_data = std::shared_ptr<T>( (T*) _mm_malloc(c.nbytes,16), [=](T *p) {
        _mm_free(p);
        p = NULL;
      });
    }

    // Unmanaged signals (temporary)
    signal1D<std::complex<T> > _a;
    _a.data   = a.m_data.get();
    _a.ndata  = a.ndata;
    _a.offset = a.offset;
    _a.nbytes = a.nbytes;

    signal1D<T> _c;
    _c.data   = c.m_data.get();
    _c.ndata  = c.ndata;
    _c.offset = c.offset;
    _c.nbytes = c.nbytes;

    // Convolve the signals
    retval = ifft<T>(_a,n,_c);

    // Prevent destruction of data
    _a.data = NULL;
    _c.data = NULL;

    return retval;
  }

  template void SPS_EXPORT DivideArray<float>(float *Data, size_t NumEl, float Divisor);
  template void SPS_EXPORT DivideArray<double>(double *Data, size_t NumEl, double Divisor);

  template struct signal1D<float>;
  template struct signal1D<double>;

  template struct signal1D<std::complex<float> >;
  template struct signal1D<std::complex<double> >;

//template struct Signal1DPlan<float>;
//template struct Signal1DPlan<double>;

// template void Signal1DPlan<float>::Wisdom(size_t** nFFTSamples, size_t* nFFTs);
//template void Signal1DPlan::Wisdom<double>(size_t** nFFTSamples, size_t* nFFTs);

  template bool SPS_EXPORT conv<float>(const signal1D<float> &a,
                                       const signal1D<float> &b,
                                       signal1D<float>& c);

  template bool SPS_EXPORT conv<double>(const signal1D<double> &a,
                                        const signal1D<double> &b,
                                        signal1D<double>& c);

  template bool SPS_EXPORT conv_fft<float>(const signal1D<float> &a,
      const signal1D<float> &b,
      signal1D<float>& c);

  template bool SPS_EXPORT conv_fft_fs<float>(const float &fs,
      const signal1D<float> &a,
      const signal1D<float> &b,
      sps::signal1D<float>& c);

  template bool SPS_EXPORT conv_fft<double>(const sps::signal1D<double> &a,
      const sps::signal1D<double> &b,
      sps::signal1D<double>& c);

  template bool SPS_EXPORT conv_fft_fs<double>(const double &fs,
      const signal1D<double> &a,
      const signal1D<double> &b,
      signal1D<double>& c);

  template bool SPS_EXPORT conv_fft_out_in<float>(const sps::signal1D<float> &a,
      sps::signal1D<float> &b);

  template bool SPS_EXPORT conv_fft_out_in<double>(const sps::signal1D<double> &a,
      sps::signal1D<double> &b);

// TEST adding SPS_EXPORT
//template class std::shared_ptr<float>;
//template class std::shared_ptr<double>;


  /* TEST TO EXPORT */
  template class msignal1D<float>;
  template class msignal1D<double>;

  template class msignal1D<std::complex<float> >;
  template class msignal1D<std::complex<double> >;

  template class std::complex<float>;
  template class std::complex<double>;

  template bool SPS_EXPORT mconv_fft<float>(const msignal1D<float> &a,
      const msignal1D<float> &b,
      msignal1D<float>& c);

  template bool SPS_EXPORT mconv_fft<double>(const msignal1D<double> &a,
      const msignal1D<double> &b,
      msignal1D<double>& c);

  template bool SPS_EXPORT mconv_fft<float>(const msignal1D<float> &a,
      const msignal1D<float> &b,
      signal1D<float>& c);

  template bool SPS_EXPORT mconv_fft<double>(const msignal1D<double> &a,
      const msignal1D<double> &b,
      signal1D<double>& c);

  template bool SPS_EXPORT mconv_fft_fs<float>(const float& fs,
      const msignal1D<float> &a,
      const msignal1D<float> &b,
      msignal1D<float>& c);

  template bool SPS_EXPORT mconv_fft_fs<double>(const double& fs,
      const msignal1D<double> &a,
      const msignal1D<double> &b,
      msignal1D<double>& c);

  template bool SPS_EXPORT fft<float>(const signal1D<float>& a, const size_t &n, signal1D<std::complex<float> >& c);
  template bool SPS_EXPORT fft<double>(const signal1D<double>& a, const size_t &n, signal1D<std::complex<double> >& c);
  template bool SPS_EXPORT ifft<float>(const signal1D<std::complex<float> >& a, const size_t &n, signal1D<float>& c);
  template bool SPS_EXPORT ifft<double>(const signal1D<std::complex<double> >& a, const size_t &n, signal1D<double>& c);

  template bool SPS_EXPORT pack_r2c<double>(const signal1D<double>& a, signal1D<std::complex<double> >& c);
  template bool SPS_EXPORT unpack_c2r<double>(const signal1D<std::complex<double> >& a, signal1D<double>& c);
  template bool SPS_EXPORT pack_r2c<float>(const signal1D<float>& a, signal1D<std::complex<float> >& c);
  template bool SPS_EXPORT unpack_c2r<float>(const signal1D<std::complex<float> >& a, signal1D<float>& c);

  template bool SPS_EXPORT mfft<float>(const msignal1D<float>& a, const size_t &n, msignal1D<std::complex<float> >& c);
  template bool SPS_EXPORT mfft<double>(const msignal1D<double>& a, const size_t &n, msignal1D<std::complex<double> >& c);
  template bool SPS_EXPORT mifft<float>(const msignal1D<std::complex<float> >& a, const size_t &n, msignal1D<float>& c);
  template bool SPS_EXPORT mifft<double>(const msignal1D<std::complex<double> >& a, const size_t &n, msignal1D<double>& c);

}

#ifdef _MSC_VER
# pragma warning(pop)
#endif


/*

// Maybe this is better

bool conv_1d2(const float_type* in, size_t n_in, float_type* out,
              const float_type* kernel, size_t kernelSize) {

  size_t i, j, k;

  // check validity of params
  if (!in || !out || !kernel)
    return false;

  // convolution from out[0] to out[kernelSize-2]
  for(i = 0; i < kernelSize - 1; ++i) {
    out[i] = 0;
    for(j = 0, k = i ; j <= i; ++j, --k)
      out[i] += in[j] * kernel[k];
  }

  // convolution from out[kernelSize-1] to out[n_in-1]
  for(i = kernelSize-1; i < n_in; ++i) {
    out[i] = 0;
    for(j = i, k = 0; k < kernelSize; --j, ++k)
      out[i] += in[j] * kernel[k];
  }

  // convolution from out[n_in] to out[n_in + kernelSize-2] (last)
  for(i = 0; i < kernelSize - 1; ++i) {
    out[n_in+i] = 0;
    for(j = 0, k= i + 1; j < (kernelSize - 1)-i; j++,k++)
      out[n_in+i] += in[n_in-kernelSize+k] * kernel[kernelSize-1-j];
  }

  return true;
}


*/














/*

fftw_plan fftw_plan_many_dft(int rank, const int *n, int howmany,
                                  fftw_complex *in, const int *inembed,
                                  int istride, int idist,
                                  fftw_complex *out, const int *onembed,
                                  int ostride, int odist,
                                  int sign, unsigned flags);
This routine plans multiple multidimensional complex DFTs, and it extends the fftw_plan_dft routine (see Complex DFTs) to compute howmany transforms, each having rank rank and size n. In addition, the transform data need not be contiguous, but it may be laid out in memory with an arbitrary stride. To account for these possibilities, fftw_plan_many_dft adds the new parameters howmany, {i,o}nembed, {i,o}stride, and {i,o}dist. The FFTW basic interface (see Complex DFTs) provides routines specialized for ranks 1, 2, and 3, but the advanced interface handles only the general-rank case.

howmany is the number of transforms to compute. The resulting plan computes howmany transforms, where the input of the k-th transform is at location in+k*idist (in C pointer arithmetic), and its output is at location out+k*odist. Plans obtained in this way can often be faster than calling FFTW multiple times for the individual transforms. The basic fftw_plan_dft interface corresponds to howmany=1 (in which case the dist parameters are ignored).

Each of the howmany transforms has rank rank and size n, as in the basic interface. In addition, the advanced interface allows the input and output arrays of each transform to be row-major subarrays of larger rank-rank arrays, described by inembed and onembed parameters, respectively. {i,o}nembed must be arrays of length rank, and n should be elementwise less than or equal to {i,o}nembed. Passing NULL for an nembed parameter is equivalent to passing n (i.e. same physical and logical dimensions, as in the basic interface.)

The stride parameters indicate that the j-th element of the input or output arrays is located at j*istride or j*ostride, respectively. (For a multi-dimensional array, j is the ordinary row-major index.) When combined with the k-th transform in a howmany loop, from above, this means that the (j,k)-th element is at j*stride+k*dist. (The basic fftw_plan_dft interface corresponds to a stride of 1.)

For in-place transforms, the input and output stride and dist parameters should be the same; otherwise, the planner may return NULL.

Arrays n, inembed, and onembed are not used after this function returns. You can safely free or reuse them.

Examples: One transform of one 5 by 6 array contiguous in memory:

        int rank = 2;
        int n[] = {5, 6};
        int howmany = 1;
        int idist = odist = 0; // unused because howmany = 1
        int istride = ostride = 1; // array is contiguous in memory
        int *inembed = n, *onembed = n;
Transform of three 5 by 6 arrays, each contiguous in memory, stored in memory one after another:

        int rank = 2;
        int n[] = {5, 6};
        int howmany = 3;
        int idist = odist = n[0]*n[1]; // = 30, the distance in memory
                                       //   between the first element
                                       //   of the first array and the
                                       //   first element of the second array
        int istride = ostride = 1; // array is contiguous in memory
        int *inembed = n, *onembed = n;
Transform each column of a 2d array with 10 rows and 3 columns:

        int rank = 1; // not 2: we are computing 1d transforms
        int n[] = {10}; // 1d transforms of length 10
        int howmany = 3;
        int idist = odist = 1;
        int istride = ostride = 3; // distance between two elements in
                                   //   the same column
        int *inembed = n, *onembed = n;


fftw_plan fftw_plan_many_dft_r2c(int rank, const int *n, int howmany,
                                      double *in, const int *inembed,
                                      int istride, int idist,
                                      fftw_complex *out, const int *onembed,
                                      int ostride, int odist,
                                      unsigned flags);
     fftw_plan fftw_plan_many_dft_c2r(int rank, const int *n, int howmany,
                                      fftw_complex *in, const int *inembed,
                                      int istride, int idist,
                                      double *out, const int *onembed,
                                      int ostride, int odist,
                                      unsigned flags);

Like fftw_plan_many_dft, these two functions add howmany, nembed,
stride, and dist parameters to the fftw_plan_dft_r2c and
fftw_plan_dft_c2r functions, but otherwise behave the same as the
basic interface.


The interpretation of howmany, stride, and dist are the same as for
fftw_plan_many_dft, above. Note that the stride and dist for the real
array are in units of double, and for the complex array are in units
of fftw_complex.

If an nembed parameter is NULL, it is interpreted as what it would be
in the basic interface, as described in Real-data DFT Array
Format. That is, for the complex array the size is assumed to be the
same as n, but with the last dimension cut roughly in half. For the
real array, the size is assumed to be n if the transform is
out-of-place, or n with the last dimension “padded” if the transform
is in-place.

If an nembed parameter is non-NULL, it is interpreted as the physical
size of the corresponding array, in row-major order, just as for
fftw_plan_many_dft. In this case, each dimension of nembed should be
>= what it would be in the basic interface (e.g. the halved or padded
n).

Arrays n, inembed, and onembed are not used after this function
returns. You can safely free or reuse them.


fftw_plan fftw_plan_r2r_1d(int n, double *in, double *out,
                                fftw_r2r_kind kind, unsigned flags);
     fftw_plan fftw_plan_r2r_2d(int n0, int n1, double *in, double *out,
                                fftw_r2r_kind kind0, fftw_r2r_kind kind1,
                                unsigned flags);
     fftw_plan fftw_plan_r2r_3d(int n0, int n1, int n2,
                                double *in, double *out,
                                fftw_r2r_kind kind0,
                                fftw_r2r_kind kind1,
                                fftw_r2r_kind kind2,
                                unsigned flags);
     fftw_plan fftw_plan_r2r(int rank, const int *n, double *in, double *out,
                             const fftw_r2r_kind *kind, unsigned flags);

Plan a real input/output (r2r) transform of various kinds in zero or
more dimensions, returning an fftw_plan (see Using Plans).

Once you have created a plan for a certain transform type and
parameters, then creating another plan of the same type and
parameters, but for different arrays, is fast and shares constant data
with the first plan (if it still exists).

The planner returns NULL if the plan cannot be created. A non-NULL
plan is always returned by the basic interface unless you are using a
customized FFTW configuration supporting a restricted set of
transforms, or for size-1 FFTW_REDFT00 kinds (which are not defined).

Arguments

rank is the dimensionality of the transform (it should be the size of
the arrays *n and *kind), and can be any non-negative integer. The
‘_1d’, ‘_2d’, and ‘_3d’ planners correspond to a rank of 1, 2, and 3,
respectively. A rank of zero is equivalent to a copy of one number
from input to output.  n, or n0/n1/n2, or n[rank], respectively, gives
the (physical) size of the transform dimensions. They can be any
positive integer.  Multi-dimensional arrays are stored in row-major
order with dimensions: n0 x n1; or n0 x n1 x n2; or n[0] x n[1] x
... x n[rank-1]. See Multi-dimensional Array Format.  FFTW is
generally best at handling sizes of the form 2a 3b 5c 7d 11e 13f,where
e+f is either 0 or 1, and the other exponents are arbitrary. Other
sizes are computed by means of a slow, general-purpose algorithm
(which nevertheless retains O(n log n) performance even for prime
sizes). (It is possible to customize FFTW for different array sizes;
see Installation and Customization.) Transforms whose sizes are powers
of 2 are especially fast.  For a REDFT00 or RODFT00 transform kind in
a dimension of size n, it is n-1 or n+1, respectively, that should be
factorizable in the above form.  in and out point to the input and
output arrays of the transform, which may be the same (yielding an
in-place transform). These arrays are overwritten during planning,
unless FFTW_ESTIMATE is used in the flags. (The arrays need not be
initialized, but they must be allocated.)  kind, or kind0/kind1/kind2,
or kind[rank], is the kind of r2r transform used for the corresponding
dimension. The valid kind constants are described in Real-to-Real
Transform Kinds. In a multi-dimensional transform, what is computed is
the separable product formed by taking each transform kind along the
corresponding dimension, one dimension after another.  flags is a
bitwise OR (‘|’) of zero or more planner flags, as defined in Planner
Flags.

*/

/* Local variables: */
/* indent-tab-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* indent-tabs-mode: nil */
/* End: */
