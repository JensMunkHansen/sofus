/**
 * @file   circular_calc.hpp
 * @author Jens Munk Hansen <jmh@debian9laptop.parknet.dk>
 * @date   Mon Aug 21 12:58:21 2017
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

#include <fnm/fnm_export.h>
#include <fnm/config.h>
#include <sps/config.h>
#include <sps/smath.hpp>

#include <fnm/circular.hpp>

namespace fnm {
  /**
   * Transient field computed. The result is in units of [Pa]. The tone burst used is +/- 1.0 V.
   *
   * @param sysparm
   * @param data Data describing the transducer
   * @param pos
   * @param nPositions
   * @param odata
   * @param nSamples
   *
   * @return
   */
  template <class T, template <typename> class A>
  T CalcCircularTransientFieldRefZeroDelay(const fnm::sysparm_t<T>* sysparm,
      const fnm::CircularApertureData<T>& data,
      const T* pos, const size_t nPositions,
      T** odata, size_t* nSamples);

  /**
   * Compute CW field from circular transducer
   *
   * @param sysparm
   * @param data
   * @param pos
   * @param nPositions
   * @param odata
   * @return success
  */
  template <class T>
  int CalcCircularCwFieldRef(const fnm::sysparm_t<T>* sysparm,
                             const fnm::CircularApertureData<T>& data,
                             const T* pos,
                             const size_t nPositions,
                             std::complex<T>** odata);
}
