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
#ifdef FNM_PULSED_WAVE
  template <class T, template <typename> class A>
  T TransientSingleRect(const sysparm_t<T>* sysparm,
                        const ApertureData<T>* data,
                        const T* pos, const size_t nPositions,
                        T** odata, size_t* nSamples, int mask = 0x1F);

  template <class T>
  T CalcFdTransientRef(const sysparm_t<T>* sysparm,
                       const ApertureData<T>* data,
                       const T* pos, const size_t nPositions, const size_t nDim,
                       T** odata, size_t* nSignals, size_t* nSamples);
#endif
}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
