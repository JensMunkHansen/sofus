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

#include <fnm/circular.hpp>
#include <fnm/circular_data.hpp>
#include <fnm/circular_calc.hpp>
#include <fnm/fnm_transient.hpp>
#include <fnm/fnm.hpp>
#include <fnm/fnm_basis.hpp>

#include <sps/smath.hpp>
#include <sps/msignals.hpp>

#include <gl/gl.hpp>

#include <cassert>

#include <vector>

// TODO: Invent structur with (r,z, dist, fStart, fStop)

namespace fnm {

  template <class T>
  const size_t& CircularAperture<T>::NMaxDivAGet() const
  {
    return this->m_sysparm->nDivA;
  }

  template <class T>
  int CircularAperture<T>::NMaxDivASet(const size_t& nDiv)
  {
    int retval = -1;
    if (nDiv > 1) {
      this->m_sysparm->nDivA = nDiv;
      retval = 0;
    }
    return retval;
  }

  template <class T>
  const T& CircularAperture<T>::F0Get() const
  {
    return m_data->m_f0;
  }

  template <class T>
  const T& CircularAperture<T>::FsGet() const
  {
    return m_sysparm->fs;
  }

  template <class T>
  void CircularAperture<T>::FsSet(const T& fs)
  {
    m_sysparm->fs = fs;
  }

  template <class T>
  const T& CircularAperture<T>::CGet() const
  {
    return m_sysparm->c;
  }

  template <class T>
  void CircularAperture<T>::CSet(const T& c)
  {
    m_sysparm->c = c;
  }

  template <class T>
  const T& CircularAperture<T>::DensityGet() const
  {
    return m_sysparm->rho;
  }

  template <class T>
  void CircularAperture<T>::DensitySet(const T& rho)
  {
    m_sysparm->rho = rho;
  }

  template <class T>
  const T& CircularAperture<T>::RadiusGet() const
  {
    return m_data->m_element.circle.radius;
  }

  template <class T>
  void CircularAperture<T>::RadiusSet(const T& radius)
  {
    m_data->m_element.circle.radius = radius;
  }

  template <class T>
  const T& CircularAperture<T>::DelayGet() const
  {
    return m_data->m_delay;
  }

  template <class T>
  void CircularAperture<T>::DelaySet(const T& delay)
  {
    m_data->m_delay = delay;
  }

  template <class T>
  void CircularAperture<T>::F0Set(const T& f0)
  {
    // Make static variable
    const T eps = std::numeric_limits<T>::epsilon();
    assert(f0 > eps);
    if (f0 > eps) {
      // Invalidate focus
      m_data->m_f0 = f0;
    }
  }

  template <class T>
  const T& CircularAperture<T>::GridSectorScaleGet() const
  {
    return m_data->m_gridSectorScale;
  }

  template <class T>
  void CircularAperture<T>::GridSectorScaleSet(const T& gss)
  {
    if (gss >= T(0.0) && gss <= T(1.0)) {
      m_data->m_gridSectorScale = gss;
    }
  }

  template <class T>
  CircularAperture<T>::CircularAperture()
  {
    // Aligned data segment
    m_data = (CircularApertureData<T>*) _mm_malloc(sizeof(CircularApertureData<T>), 4*sizeof(T));

    // Consider moving initializer
    new (this->m_data) CircularApertureData<T>();

    // Assign pointer to global sysparm
    this->m_sysparm = Aperture<T>::DefaultSysParmGet();
  }

  template <class T>
  CircularAperture<T>::CircularAperture(T radius) : CircularAperture<T>()
  {
    m_data->m_element.circle.radius = radius;
    m_data->Initialize();
  }

  template <class T>
  CircularAperture<T>::~CircularAperture()
  {
    if (m_data) {
      m_data->~CircularApertureData();
      _mm_free(m_data);
      m_data = nullptr;
    }

    // Free local system parameters if they are set
    if (m_sysparm != Aperture<T>::DefaultSysParmGet()) {
      delete m_sysparm;
      m_sysparm = nullptr;
    }
  }

  template <class T>
  const T& CircularAperture<T>::WGet() const
  {
    return m_sysparm->w;
  }

  template <class T>
  void CircularAperture<T>::WSet(const T& w)
  {
    m_sysparm->w = w;
  }

  template <class T>
  T CircularAperture<T>::CalcTransientFieldRef(const T* pos, const size_t nPositions, const size_t nDim,
      T** odata, size_t* nSignals, size_t* nSamples)
  {
    T retval = T(0.0);
    assert(nDim == 3);

    if (nDim != 3) {
      *odata = NULL;
      *nSignals = 0;
      *nSamples = 0;
      retval = T(-1.0);
      return retval;
    }

    // TODO: Call HanningWeightedPulse is excitation
    T tStart = CalcCircularTransientFieldRefZeroDelay<T, ToneBurst>(this->m_sysparm, *this->m_data,
               pos, nPositions,
               odata, nSamples);
    *nSignals = nPositions;
    return tStart;
  }


  template <class T>
  int CircularAperture<T>::CalcCwFieldRef(const T* pos, const size_t nPositions, const size_t nDim,
                                          std::complex<T>** odata, size_t* nOutPositions)
  {
    int retval = 0;
    assert(nDim == 3);
    if (nDim != 3) {
      *odata = NULL;
      *nOutPositions = 0;
      retval = -1;
      return retval;
    }

    retval = CalcCircularCwFieldRef<T>(this->m_sysparm, *this->m_data, pos, nPositions, odata);
    if (retval == 0) {
      *nOutPositions = nPositions;
    } else {
      *nOutPositions = 0;
    }

    return retval;
  }

#ifdef FNM_DOUBLE_SUPPORT
  template class CircularAperture<double>;
#endif
  template class CircularAperture<float>;
}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
