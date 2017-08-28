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
  template <class T>
  T CalcCircularTransientFieldRef(const fnm::sysparm_t<T>* sysparm,
                                  const fnm::CircularApertureData<T>& data,
                                  const T* pos, const size_t nPositions,
                                  T** odata, size_t* nSamples);
  template <class T>
  T CalcCircularTransientFieldRefZeroDelay(const fnm::sysparm_t<T>* sysparm,
      const fnm::CircularApertureData<T>& data,
      const T* pos, const size_t nPositions,
      T** odata, size_t* nSamples);

  template <class T>
  int CalcCircularCwFieldRef(const fnm::sysparm_t<T>* sysparm,
                             const fnm::CircularApertureData<T>& data,
                             const T* pos,
                             const size_t nPositions,
                             std::complex<T>** odata);
}
