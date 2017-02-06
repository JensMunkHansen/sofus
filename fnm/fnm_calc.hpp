/**
 * @file   fnm_calc.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Mon Jun 13 08:33:33 2016
 *
 * @brief Function used for Fast-Nearfield-Method
 *
 *
 */

#pragma once

#include <fnm/fnm.hpp>
#include <fnm_data.hpp>
#include <sps/stdlib.h>
#include <complex>

namespace fnm {

#if HAVE_PTHREAD_H
  extern pthread_t threads[N_MAX_THREADS];
  extern pthread_attr_t attr;
#endif

  /**
   *
   *  Field above a corner of an element with dimensions \f$s\f$ and \f$l\f$ computed using Gauss-Legendre integration
   *
   *  \f{eqnarray*}{H_{s,l}(z;k) =
   *   - i/(2\pi k) ( &s \int_0^l (exp(-ik\sqrt{\sigma^2+z^2+s^2}) - exp(-ikz)) / (\sigma^2+s^2) d\sigma + \\
   *                  &l \int_0^s (exp(-ik\sqrt{\sigma^2+z^2+l^2}) - exp(-ikz)) / (\sigma^2+l^2) d\sigma) \f}
   *
   *  Reference implementation, see <a href="http://www.ncbi.nlm.nih.gov/pmc/articles/pmid/15139602/">Rapid calculations of time-harmonic nearfield pressures produced by rectangular pistons, J. McGough, 2004</a>
   *
   * @param s width (x-dimension)
   * @param l height (y-dimension)
   * @param z distance to element
   * @param k wave-number
   * @param uxs abcissa coordinates for x-coordinate
   * @param uweights weights for x-coordinate
   * @param nUs number of x-coordinates
   * @param vxs abcissa coordinates for y-coordinate
   * @param vweights  weights for y-coordinate
   * @param nVs number of y-coordinates
   *
   * @return Complex field value
   */
  template <class T>
  std::complex<T> CalcHz(const T& s,
                         const T& l,
                         const T& z,
                         const T& k,
                         const T* uxs,
                         const T* uweights,
                         const size_t nUs,
                         const T* vxs,
                         const T* vweights,
                         const size_t nVs);


  /**
   *
   * Spatial impulse response computed at the positions pos for a
   * frequency f0. The frequency f0 as well as the positional
   * information for the aperture is contained in the data structure.
   \verbatim
               ^
               |
         +-----+-----+-----------+ (x,y)
         |     |     |           |
         +-----+-----+-----------+
         |     |     |           |
         |     |     |           |
         |     |     |           |
         +     o-----+-----------+---->
         |           |           |
     hh  |           |           |
         |           |           |
         +-----+-----+-----------+
                 hw

   \endverbatim
   *
   * (Scalar reference implementation)
   *
   * @param data          Structure holding element positions and f0
   * @param pos           Positions
   * @param nPositions    # of positions
   * @param odata         complex output
   * @return error code
   */
  template <class T>
  int CalcCwFieldRef(const ApertureData<T>& data,
                     const T* pos, const size_t nPositions,
                     std::complex<T>** odata);

  /**
   *
   *
   * @param data
   * @param pos
   * @param nPositions
   * @param odata
   *
   * @return
   */
  template <class T>
  int CalcCwThreaded(const ApertureData<T>* data,
                     const T* pos, const size_t nPositions,
                     std::complex<T>** odata);

  /**
   * Naive integral, no re-use of parts and many abcissa are needed when projection is far away from rectangle
   *
   * @param data
   * @param pos
   * @param nPositions
   * @param odata
   */
  template <class T>
  void CalcCwField(const ApertureData<T>& data,
                   const T* pos, const size_t nPositions,
                   std::complex<T>** odata);

  /**
   * The positions must all equal the focus point and the number
   * should match the number of elements. By doing so, a set of
   * phases is computed for focusing.
   *
   * Used by \ref FocusUpdate
   *
   * @param aperture data
   * @param pos
   * @param nPositions
   * @param odata
   *
   * @return
   */
  template <class T>
  int CalcCwFocus(const ApertureData<T>& data,
                  const T* pos, const size_t nPositions,
                  std::complex<T>** odata);
}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
