/**
 * @file   signals.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Fri May 29 02:44:47 2015
 * 
 * @brief  
 * 
 * 
 */

// TODO: Make an interface

#pragma once

#ifndef SPS_EXPORT
# include <sps/sps_export.h>
#endif

#include <cstddef>
#include <sps/mm_malloc.h>
#include <cstring>   // memset
#include <algorithm> // min/max
#include <memory>    // shared_ptr

// TODO: Check for availability of <complex>
#if defined(__GNUG__) || (defined(_MSC_VER) && (_MSC_VER >= 1800))
# include <complex>
#else
# include <complex>
#endif

namespace sps {

  template <typename T>
  void DivideArray(T *Data, size_t NumEl, T Divisor);
  
  template<typename T>
  struct SPS_EXPORT msignal1D {
    msignal1D() : offset(0), ndata(0), nbytes(0), data(NULL) {}

    msignal1D(size_t ndata) : msignal1D(ndata, 16 * (ndata * sizeof(T) + 15) / 16) {}
      
    msignal1D(size_t _ndata, size_t _nbytes) : offset(0),
                                               ndata(_ndata),
      nbytes(16 * ( std::max<size_t>(_ndata * sizeof(T),_nbytes) + 15 ) / 16), data(NULL) {
      data.reset( (T*) _mm_malloc(nbytes,16), [=](T *p) {_mm_free(p); p=NULL;});
#ifndef NDEBUG
      memset(data.get(),0,nbytes);
#endif
    }

    /** 
     * Reset signal, set it to zero
     * 
     * @param _ndata 
     * @param _nbytes 
     */
    void reset(size_t _ndata = 0, size_t _nbytes = 0) {
      // Maximum of current and number of elements
      _ndata  = std::max<size_t>(_ndata,ndata);

      // Maximum of current and new size
      _nbytes = std::max<size_t>(_nbytes, _ndata*sizeof(T));

      // Re-alloc if needed
      if (_nbytes > nbytes) {
        nbytes = _nbytes;
        data.reset( (T*) _mm_malloc(nbytes,16), [=](T *p) {_mm_free(p); p=NULL;});
      }
      // Reset data
      ndata = _ndata;
      memset(data.get(),0,ndata*sizeof(T));
    }

    // TODO: Consider padding to number of elements
    void pad(size_t _nbytes) {
      if ((nbytes > 0) && (_nbytes > nbytes)) {
        // Shallow copy
        std::shared_ptr<T> data1 = data;
        data.reset( (T*) _mm_malloc(_nbytes,16), [=](T *p) {_mm_free(p); p=NULL;});
        memset(data.get(),0,_nbytes);
        memcpy(data.get(), data1.get(), ndata*sizeof(T));
        nbytes = _nbytes;
      }
    }
    ~msignal1D() {
      if (ndata) {
        if (data) {
          data.reset((T*)NULL);
        }
        //        printf("Freed signal");
      }
    }
    
    int offset;              // Offset can be negative
    size_t ndata;            // Number of samples
    size_t nbytes;           // Number of bytes - may exceed sizeof(T) * ndata
    std::shared_ptr<T> data; // Data
  };


  /** 
   * Convolution of managed signals with data held by shared_ptr's
   * 
   * @param a Input 
   * @param b Input 
   * @param c Output : length is len(a) + len(b) - 1
   * 
   * @return 
   */  
  template<typename T>
  bool mconv_fft(const msignal1D<T>& a, const msignal1D<T>& b, msignal1D<T>& c);

  /** 
   * Convolution of managed signals with data held by
   * shared_ptr's. The result is multiplied by 1/fs to get units
   * right.
   * 
   * @param fs Sampling frequency 
   * @param a Input
   * @param b Input
   * @param c Ouput : length is len(a) + len(b) - 1
   * 
   * @return 
   */  
  template<typename T>
  bool mconv_fft_fs(const T& fs, const msignal1D<T>& a, const msignal1D<T>& b, msignal1D<T>& c);

  template <typename T>
  class SPS_EXPORT Signal1DPlan {
  public:
    static Signal1DPlan& Instance();
    //    static void Wisdom(size_t** nFFTSamples, size_t* nFFTs);
  private:
    Signal1DPlan();
    ~Signal1DPlan();
    Signal1DPlan(Signal1DPlan const&) = default;
    Signal1DPlan& operator=(Signal1DPlan const&) = default;
  };

  template <typename T>
  struct signal1D {
    signal1D() : data(NULL), offset(0), ndata(0), nbytes(0) {}
    signal1D(size_t ndata) : data(NULL), offset(0), ndata(ndata) {
      nbytes = 16 * (ndata * sizeof(T) + 15) / 16;
      data = (T*) _mm_malloc(nbytes,16);
#ifndef NDEBUG
      memset(data,0,nbytes);
#endif
    }
    signal1D(const signal1D& other) = default;
 
    signal1D(size_t ndata, size_t _nbytes) : data(NULL), offset(0), ndata(ndata), nbytes(0) {
      nbytes = 16 * (ndata * sizeof(T) + 15) / 16;
      nbytes = std::max<size_t>(_nbytes,nbytes);
      data = (T*) _mm_malloc(nbytes,16);
#ifndef NDEBUG
      memset(data,0,nbytes);
#endif
    }
    ~signal1D() {
      if (data) {
        _mm_free(data);
        data = NULL;
      }
    }

    T* data;
    int offset;    // Offset can be negative
    size_t ndata;  // Number of samples
    size_t nbytes; // Number of bytes - may exceed sizeof(T) * ndata
  };

  /** 
   * 
   * 
   * @param a 
   * @param b 
   * @param c unmanaged output - length is a.ndata + b.ndata - 1
   * 
   * @return 
   */
  template<typename T>
  bool mconv_fft(const msignal1D<T>& a, const msignal1D<T>& b, signal1D<T>& c);
  
  /** 
   * Convolution of 1D signals (using FFT)
   * 
   * @param a Input A - length is a.ndata
   * @param b Input B - length is b.ndata
   * @param c Output  - length is a.ndata + b.ndata - 1
   * 
   * @return 
   */
  template <typename T>
  bool SPS_EXPORT conv_fft(const sps::signal1D<T> &a,
                           const sps::signal1D<T> &b,
                           sps::signal1D<T>& c);

  /** 
   * Convolution
   * 
   * @param fs 
   * @param a 
   * @param b 
   * @param c 
   * 
   * @return 
   */
  template <typename T>
  bool conv_fft_fs(const T& fs,
                   const sps::signal1D<T> &a,
                   const sps::signal1D<T> &b,
                   sps::signal1D<T>& c);

  /** 
   * In-place convolution (2nd argument)
   * 
   * @param a 
   * @param b 
   * 
   * @return 
   */  
  template <typename T>
  bool conv_fft_out_in(const signal1D<T> &a,
                       signal1D<T> &b);

  template <typename T>
  bool SPS_EXPORT pack_r2c(const sps::signal1D<T> &a, sps::signal1D<std::complex<T> >& c);

  template <typename T>
  bool SPS_EXPORT unpack_c2r(const sps::signal1D<std::complex<T> >& a, sps::signal1D<T> &c);
  
  /** 
   * FFT of signal
   * 
   * @param a input, length is a.ndata
   * @param c output, length is the 
   * 
   * @return 
   */
  template <typename T>
  bool SPS_EXPORT fft(const sps::signal1D<T> &a, const size_t &n, sps::signal1D<std::complex<T> >& c);

  /** 
   * FFT of signal
   * 
   * @param a input, length is a.ndata
   * @param c output, length is the 
   * 
   * @return 
   */
  template <typename T>
  bool SPS_EXPORT ifft(const sps::signal1D<std::complex<T> >& a, const size_t &n, sps::signal1D<T> &c);
  
#if defined(__GNUG__) || (defined(_MSC_VER) && (_MSC_VER >= 1800))

  template <typename T>
  bool mfft(const sps::msignal1D<T> &a, const size_t &n, sps::msignal1D<std::complex<T> >& c);

  template <typename T>
  bool mifft(const sps::msignal1D<std::complex<T> >& a, const size_t &n, sps::msignal1D<T> &c);
  
#endif

  /** 
   * Convolution of 1D signals (there is a bug)
   * 
   * @param a Input A - length is na
   * @param b Input B - length is nb (must be less than na)
   * @param c Output  - length is a.ndata + b.ndata - 1
   * 
   * @return 
   */
  template <typename T>
  bool SPS_EXPORT conv(const sps::signal1D<T> &a,
                       const sps::signal1D<T> &b,
                       sps::signal1D<T>& c);
}

/* Local variables: */
/* indent-tab-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* indent-tabs-mode: nil */
/* End: */
