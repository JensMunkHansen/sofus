/**
 * @file   signals.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Fri May 29 02:44:47 2015
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

// TODO: Create interface without shared_ptr member or use PIMPL

#pragma once

#include <sps/sps_export.h>

#include <cstddef>
#include <memory>    // shared_ptr
#include <complex>

namespace std {
  // Not nice
  template class SPS_EXPORT std::shared_ptr<float>;
  template class SPS_EXPORT std::shared_ptr<double>;
  template class SPS_EXPORT std::shared_ptr<std::complex<float> >;
  template class SPS_EXPORT std::shared_ptr<std::complex<double> >;

#ifdef _MSC_VER
  template class SPS_EXPORT shared_ptr<float>;
  template class SPS_EXPORT shared_ptr<double>;
#endif
}

namespace sps {

  template <typename T>
  void SPS_EXPORT DivideArray(T *Data, size_t NumEl, T Divisor);

  template<typename T>
  class SPS_EXPORT msignal1D {
  public:
    msignal1D();

    msignal1D(size_t ndata);

    msignal1D(size_t _ndata, size_t _nbytes);

#ifndef SWIG_VERSION
    msignal1D& operator=(const sps::msignal1D<T>& a);

    msignal1D(const sps::msignal1D<T>& a);

    inline T& operator[](size_t index)
    {
      return this->m_data.get()[index];
    }

    inline const T& operator[](size_t index) const
    {
      return this->m_data.get()[index];
    }

    inline T* get()
    {
      return this->m_data.get();
    }
#endif
    void scale(const T& s);
    /**
     * Reset signal, set it to zero
     *
     * @param _ndata
     * @param _nbytes
     */
    void reset(size_t _ndata = 0, size_t _nbytes = 0);

    // TODO: Consider padding to number of elements
    void pad(size_t _nbytes);

    /**
     * Reverse the data
     *
     */
    void reverse();

    ~msignal1D();

    // TODO: Consider using nsamples and (ndata >= sizeof(T)*nsamples)
    int offset;              // Offset can be negative
    size_t ndata;            // Number of samples
    size_t nbytes;           // Number of bytes - may exceed sizeof(T) * ndata

    std::shared_ptr<T>& data()
    {
      return this->m_data;
    }
  private:
  public:

    // Should not be here
    std::shared_ptr<T> m_data; // Data
  };

  template <typename T>
  struct SPS_EXPORT signal1D {
    signal1D();
    signal1D(size_t ndata);
    signal1D(const signal1D& other) = default;
    signal1D(size_t ndata, size_t _nbytes);
    ~signal1D();

    T* data;
    int offset;    // Offset can be negative
    size_t ndata;  // Number of samples
    size_t nbytes; // Number of bytes - may exceed sizeof(T) * ndata
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
  bool SPS_EXPORT mconv_fft(const msignal1D<T>& a, const msignal1D<T>& b, msignal1D<T>& c);

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
  bool SPS_EXPORT mconv_fft_fs(const T& fs, const msignal1D<T>& a, const msignal1D<T>& b, msignal1D<T>& c);

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
  bool SPS_EXPORT mconv_fft(const msignal1D<T>& a, const msignal1D<T>& b, signal1D<T>& c);

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
  bool SPS_EXPORT conv_fft_fs(const T& fs,
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
   * @param a real input, length is a.ndata
   * @param c complex output, length is a.ndata/2 + 1
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
  bool SPS_EXPORT mfft(const sps::msignal1D<T> &a, const size_t &n, sps::msignal1D<std::complex<T> >& c);

  template <typename T>
  bool SPS_EXPORT mifft(const sps::msignal1D<std::complex<T> >& a, const size_t &n, sps::msignal1D<T> &c);

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
