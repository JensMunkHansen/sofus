/**
 * @file   fnm_calc.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Mon Jun 13 08:33:33 2016
 *
 * @brief Function used for Fast-Nearfield-Method
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

#pragma once

#include <fnm_data.hpp>
#include <sps/stdlib.h>
#include <sps/progress.hpp>
#include <complex>

namespace fnm {
  template <class T>
  struct GLQuad2D;

#if HAVE_PTHREAD_H
  extern pthread_t threads[N_MAX_THREADS];
  extern pthread_attr_t attr;
#endif

  /** @defgroup fnm_calc_functions FNM Calculation functions
   *  @ingroup fnm
   * @{
   */

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
   * @param gl Abcissae and weights for integration
   *
   * @return Complex field value
   */
  template <class T>
  std::complex<T> CalcHz(const T& s,
                         const T& l,
                         const T& z,
                         const T& k,
                         const GLQuad2D<T>* gl);


  /**
   *
   * Spatial impulse response computed at the positions pos for a
   * frequency f0. The principle of superposition is used as shown
   * below and the ranges of integration are reduced to give the most
   * accurate result with the least amount of abcissas. This is a
   * reference implementation. Eight integrals are evaluated, which
   * can be reduced to four. See @ref CalcCwFieldFourRef
   *
   * The frequency f0 as well as the positional information for the
   * aperture is contained in the @ref ApertureData structure
   * data.
   *
   * Superposition principle:
   \verbatim
               ^
               |
         +-----+-----+-----------+ (x,y)
         |     |     |           |
         +-----+-----+-----------+
         |     |     |           |
         |     |     |           |
         |     |     |           |        =
         +     o-----+-----------+---->
         |           |           |
     hh  |           |           |
         |           |           |
         +-----+-----+-----------+
                 hw

               ^                             ^                                     ^                                      ^
               |                             |                                     |                                      |
         +-----------------------+ (x,y)     |     +-----------+ (x,y)       +-----+-----+-----------+ (x,y)              |     +-----------+ (x,y)
         |     |                 |           |     |           |             |     |                 |                    |     |           |
         +     |                 |           |     |           |             +-----+-----------------+              +     |     +-----------+
         |     |                 |           |     |           |             |     |                                |     |
         |     |                 |           |     |           |         hh  |     |                            hh  |     |
         |     |                 |       -   |     |           |       -     |     |                        +       |     |
         |     o-----------------|---->      o-----+-----------+---->        +     o---------------------->         +     o---------------------->
         |                       |                 |           |                                                    |
     hh  |                       |                 |           |                                                    |
         |                       |                 |           |                                                    |
         +-----------+-----------+                 +-----------+                                                    +-----+-----+
                 hw                            hw                                                                           hw

   \endverbatim
   *
   * \f{eqnarray}{
   * H(x,y,z,k) &= H_{hw+|x|,hh+|y|}(z;k) - H_{|x|-hw,hh+|y|}(z;k) - H_{hw+|x|,|y|-hh}(z;k) + H_{|x|-hw,|y|-hh}(z;k)\nonumber\\
   * &= H_{hw+|x|,hh+|y|}(z;k) + H_{hw-|x|,hh+|y|}(z;k) + H_{hw+|x|,hh-|y|}(z;k) + H_{hw-|x|,hh-|y|}(z;k)\nonumber\\
   * &= -\frac{i}{2\pi k} \Big( (hh+|y|)\int_{0}^{|x|+hw} M(hh+|y|,\sigma, z, k)d\sigma +(hw+|x|)\int_{0}^{|y|+hh} M(hw+|x|,\sigma, z, k)d\sigma \nonumber\\
   * &\phantom{= -\frac{i}{2\pi k} \Big(} -(hh+|y|)\int_{0}^{|x|-hw} M(hh+|y|,\sigma, z, k)d\sigma - (hw-|x|)\int_{0}^{|y|-hh} M(hw-|x|,\sigma, z, k)d\sigma\nonumber\\
   * &\phantom{= -\frac{i}{2\pi k} \Big(} - (hh-|y|)\int_{0}^{|x|-hw} M(hh-|y|,\sigma, z, k)d\sigma -(hw+|x|)\int_{0}^{|y|-hh} M(hw+|x|,\sigma, z, k)d\sigma \nonumber\\
   * &\phantom{= -\frac{i}{2\pi k} \Big(} + (hh-|y|)\int_{0}^{|x|+hw} M(hh-|y|,\sigma, z, k)d\sigma + (hw-|x|)\int_{0}^{|y|+hh} M(hw-|x|,\sigma, z, k)d\sigma\nonumber\\
   * &= -\frac{i}{2\pi k} \Big( (hh+|y|)\int_{|x|-hw}^{|x|+hw} M(hh+|y|,\sigma, z, k)d\sigma + (hh-|y|)\int_{|x|-hw}^{|x|+hw} M(hh-|y|,\sigma, z, k)d\sigma \nonumber\\
   * &\phantom{= -\frac{i}{2\pi k} \Big(} +  (hw+|x|)\int_{|y|-hh}^{|y|+hh} M(hw+|x|,\sigma, z, k)d\sigma + (hw-|x|)\int_{|y|-hh}^{|y|+hh} M(hw-|x|,\sigma, z, k)d\sigma\Big),
   * \f}
   *
   * where \f$M(s,\sigma,z,k)\f$ is the kernel function
   *\f[
   * M(s,\sigma,z,k) = \frac{exp(-ik\sqrt{\sigma^2+z^2+s^2}) - exp(-ikz)}{\sigma^2+s^2}
   \f]
   *
   * @param sysparm
   * @param data          Structure holding element positions and f0
   * @param pos           Positions
   * @param nPositions    # of positions
   * @param odata         complex output
   * @return error code
   */
#ifdef USE_PROGRESS_BAR
  template <class T>
  int CalcCwFieldRef(const sysparm_t<T>* sysparm,
                     const ApertureData<T>* data,
                     const T* pos, const size_t nPositions,
                     std::complex<T>** odata,
                     sps::ProgressBarInterface* pBar);
#else
  template <class T>
  int CalcCwFieldRef(const sysparm_t<T>* sysparm,
                     const ApertureData<T>* data,
                     const T* pos, const size_t nPositions,
                     std::complex<T>** odata,
                     void* pBar);
#endif

  /**
   * CalcCwFieldFourRef
   *
   * Same as @ref CalcCwFieldRef except the number of integrals is reduced to four
   *
   * @param sysparm
   * @param data
   * @param pos
   * @param nPositions
   * @param odata
   *
   * @return
   */
  template <class T>
  int CalcCwFieldFourRef(const sysparm_t<T>* sysparm,
                         const ApertureData<T>* data,
                         const T* pos, const size_t nPositions,
                         std::complex<T>** odata);

  /*** @} */

  /**
   * Compute CW response at multiple positions. The range of
   * integration is naive so it is accurate when projections lie
   * inside an element, but requires a huge amount of abcissas to
   * get a usuable result, when projections lie outside an element.
   *
   * @param sysparm
   * @param data
   * @param pos
   * @param nPositions
   * @param odata
   *
   * @return
   */
  template <class T>
  int CalcCwFocusNaiveFast(const sysparm_t<T>* sysparm,
                           const ApertureData<T>& data,
                           const T* pos, const size_t nPositions,
                           std::complex<T>** odata);

  /**
   * CalcSingleFast. TODO: Implement this for double
   *
   * @param s1            lower limit
   * @param s2            upper limit
   * @param l             adjacent edge
   * @param z             z-coordinate
   * @param k             wave numer
   * @param uxs           abcissae
   * @param uweights      weights
   * @param nUs           number of abcissae
   *
   * @return Complex field value
   */
  template <class T>
  STATIC_INLINE_BEGIN
  std::complex<T> CalcSingleFast(const T& s1,
                                 const T& s2,
                                 const T& l,
                                 const T& z,
                                 const T& k,
                                 const T* uxs,
                                 const T* uweights,
                                 const size_t nUs);

  /**
   * Multi-threaded version of @ref CalcCwFieldRef. This is the function to use.
   *
   * @param sysparm
   * @param data
   * @param pos
   * @param nPositions
   * @param odata
   *
   * @return
   */
  template <class T>
  int CalcCwThreaded(const fnm::sysparm_t<T>* sysparm,
                     const ApertureData<T>* data,
                     const T* pos, const size_t nPositions,
                     std::complex<T>** odata, sps::ProgressBarInterface* pbar);

  /**
   * Naive integral, no re-use of parts and many abcissa are needed when projection is far away from rectangle.
   *
   * Consider renaming to CalcCwNaive
   *
   * @param sysparm
   * @param data
   * @param pos
   * @param nPositions
   * @param odata
   */
  template <class T>
  void CalcCwField(const sysparm_t<T>* sysparm,
                   const ApertureData<T>& data,
                   const T* pos, const size_t nPositions,
                   std::complex<T>** odata);

  /**
   * The positions must all equal the focus point and the number
   * should match the number of elements. By doing so, a set of phases
   * is computed for focusing. The function uses @ref CalcHzFast for
   * phase calculations.
   *
   * @param sysparm
   * @param data
   * @param pos
   * @param nPositions
   * @param odata
   *
   * @return
   */
  template <class T>
  int CalcCwFocus(const sysparm_t<T>* sysparm,
                  const ApertureData<T>& data,
                  const T* pos, const size_t nPositions,
                  std::complex<T>** odata);

  /**
   * Reference implementation of the function above. Uses CalcCwFieldRef for phase calculations
   *
   * @param sysparm
   * @param data
   * @param pos
   * @param nPositions
   * @param odata
   *
   * @return
   */
  template <class T>
  int CalcCwFocusRef(const sysparm_t<T>* sysparm,
                     const ApertureData<T>& data,
                     const T* pos, const size_t nPositions,
                     std::complex<T>** odata);
}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
