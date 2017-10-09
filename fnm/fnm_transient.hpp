/**
 * @file   fnm_transient.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Thu May  4 20:25:38 2017
 *
 * @brief
 *
 *
 */
#pragma once

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

#include <fnm/fnm_types.hpp>
#include <sps/math.h>
#include <sps/smath_types.hpp>

#include <sofus/sofus_types.hpp>

namespace fnm {
#ifdef _MSC_VER
  template<class T>
  class ApertureData;
#else
  template<class T>
  struct ApertureData;
#endif
}

namespace fnm {
  template <class T, template <typename> class A>
  T TransientSingleRect(const sysparm_t<T>* sysparm,
                        const ApertureData<T>* data,
                        const T* pos, const size_t nPositions,
                        T** odata, size_t* nSamples, int mask = 0x1F);

  /**
   * Note: Uses scaled Gauss-Legendre weights and abcissae
   *
   * @param element
   * @param limit   Projection and limits for sphere intersecting element
   * @param uv      Gauss-Legendre weights and abcissae scaled by hw (for u) and hh (for v)
   *
   * @return
   */
  template <class T>
  T FourDirect(const sps::element_rect_t<T>* element,
               const sofus::proj_limit_dist_t<T>* limit,
               const GLQuad2D<T>* uv);

  /**
   * Compute direct response for rectangle
   *
   * @param sysparm
   * @param element
   * @param scale
   * @param f0
   * @param pld
   * @param uv
   * @param delay
   * @param iSampleSignalStart
   * @param nSamples
   * @param odata
   */
  template <class T, template <typename> class A>
  void DirectWaveSingle(const sysparm_t<T>* sysparm,
                        const sps::element_rect_t<T>* element,
                        const T& scale,
                        const T& f0,
                        const sofus::proj_limit_dist_t<T>* pld,
                        const GLQuad2D<T>* uv,
                        const T delay,
                        const int iSampleSignalStart,
                        const size_t nSamples,
                        T* odata);

  /**
   * Compute edge response for rectangle
   *
   * Note: Uses scaled Gauss-Legendre weights and abcissae
   *
   * @param sysparm
   * @param element
   * @param scale
   * @param iEdge
   * @param pld
   * @param pGL
   * @param f0
   * @param delay
   * @param iSampleSignalStart
   * @param nSamples
   * @param odata
   */
  template <typename T, template <typename> class A>
  void EdgeResponse(const sysparm_t<T>* sysparm,
                    const sps::element_rect_t<T>* element,
                    const T& scale,
                    const size_t iEdge,
                    const sofus::proj_limit_dist_t<T>* pld,
                    const GLQuad1D<T>* pGL,
                    const T& f0,
                    const T& delay,
                    const int& iSampleSignalStart,
                    const size_t& nSamples,
                    T* odata);

}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
