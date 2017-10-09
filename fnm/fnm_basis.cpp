/**
 * @file   fnm_basis.cpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Thu May  4 20:25:38 2017
 *
 * @brief
 *
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

#include <fnm/config.h>
#include <fnm/fnm_basis.hpp>

namespace fnm {
  template <class T>
  T (*const Identity<T>::SpatialBasisFunction[Identity<T>::nTerms])(T,T,T) = {
    Identity<T>::sbf0,
  };

  template <class T>
  T (*const Identity<T>::TemporalBasisFunction[Identity<T>::nTerms])(T,T,T) = {
    Identity<T>::tbf0,
  };

  template <class T>
  Identity<T>::Identity()
  {
    std::fill(m_fTerms, m_fTerms + Identity<T>::nTerms, T(0.0));
  }

  template <class T>
  inline void Identity<T>::ResetSpatial()
  {
    std::fill(m_fTerms, m_fTerms + Identity<T>::nTerms, T(0.0));
  }

  template <class T>
  inline void Identity<T>::UpdateSpatial(T factor, T tau, T W, T f0)
  {
    for (size_t iTerm = 0 ; iTerm < Identity<T>::nTerms ; iTerm++) {
      this->m_fTerms[iTerm] += factor * Identity<T>::SpatialBasisFunction[iTerm](tau, W, f0);
    }
  }
  /*
  template <class T>
  inline T Identity<T>::EvaluateTSD(T t, T W, T f0) const
  {
    T result = T(0.0);
    for (size_t iTerm = 0 ; iTerm < Identity<T>::nTerms ; iTerm++) {
      result += this->m_fTerms[iTerm] * Identity<T>::TemporalBasisFunction[iTerm](t,W,f0);
    }
    return result;
  }
  */

  template <class T>
  T (*const ToneBurst<T>::SpatialBasisFunction[ToneBurst<T>::nTerms])(T,T,T) = {
    ToneBurst<T>::sbf0,
    ToneBurst<T>::sbf1,
  };

  template <class T>
  T (*const ToneBurst<T>::TemporalBasisFunction[ToneBurst<T>::nTerms])(T,T,T) = {
    ToneBurst<T>::tbf0,
    ToneBurst<T>::tbf1,
  };

  template <class T>
  ToneBurst<T>::ToneBurst()
  {
    std::fill(m_fTerms, m_fTerms + ToneBurst<T>::nTerms, T(0.0));
  }

  template <class T>
  inline void ToneBurst<T>::ResetSpatial()
  {
    std::fill(m_fTerms, m_fTerms + ToneBurst<T>::nTerms, T(0.0));
  }

  template <class T>
  inline void ToneBurst<T>::UpdateSpatial(T factor, T tau, T W, T f0)
  {
    for (size_t iTerm = 0 ; iTerm < ToneBurst<T>::nTerms ; iTerm++) {
      this->m_fTerms[iTerm] += factor * ToneBurst<T>::SpatialBasisFunction[iTerm](tau, W, f0);
    }
  }
  /*
  template <class T>
  inline T ToneBurst<T>::EvaluateTSD(T t, T W, T f0) const
  {
    T result = T(0.0);
    for (size_t iTerm = 0 ; iTerm < ToneBurst<T>::nTerms ; iTerm++) {
      result += this->m_fTerms[iTerm] * ToneBurst<T>::TemporalBasisFunction[iTerm](t,W,f0);
    }
    return result;
  }
  */

  template <class T>
  T (*const HanningWeightedPulse<T>::SpatialBasisFunction[HanningWeightedPulse<T>::nTerms])(T,T,T) = {
    HanningWeightedPulse<T>::sbf0,
    HanningWeightedPulse<T>::sbf1,
    HanningWeightedPulse<T>::sbf2,
    HanningWeightedPulse<T>::sbf3,
    HanningWeightedPulse<T>::sbf4,
    HanningWeightedPulse<T>::sbf5
  };

  template <class T>
  T (*const HanningWeightedPulse<T>::TemporalBasisFunction[HanningWeightedPulse<T>::nTerms])(T,T,T) = {
    HanningWeightedPulse<T>::tbf0,
    HanningWeightedPulse<T>::tbf1,
    HanningWeightedPulse<T>::tbf2,
    HanningWeightedPulse<T>::tbf3,
    HanningWeightedPulse<T>::tbf4,
    HanningWeightedPulse<T>::tbf5
  };

  template <class T>
  HanningWeightedPulse<T>::HanningWeightedPulse()
  {
    std::fill(m_fTerms, m_fTerms + HanningWeightedPulse<T>::nTerms, T(0.0));
  }

  template <class T>
  inline void HanningWeightedPulse<T>::ResetSpatial()
  {
    std::fill(m_fTerms, m_fTerms + HanningWeightedPulse<T>::nTerms, T(0.0));
  }

  template <class T>
  inline void HanningWeightedPulse<T>::UpdateSpatial(T factor, T tau, T W, T f0)
  {
    for (size_t iTerm = 0 ; iTerm < HanningWeightedPulse<T>::nTerms ; iTerm++) {
      this->m_fTerms[iTerm] += factor * HanningWeightedPulse<T>::SpatialBasisFunction[iTerm](tau, W, f0);
    }
  }
  /*
  template <class T>
  T HanningWeightedPulse<T>::EvaluateTSD(T t, T W, T f0) const
  {
    T result = T(0.0);
    for (size_t iTerm = 0 ; iTerm < HanningWeightedPulse<T>::nTerms ; iTerm++) {
      result += this->m_fTerms[iTerm] * HanningWeightedPulse<T>::TemporalBasisFunction[iTerm](t,W,f0);
    }
    return result;
  }
  */

#ifdef FNM_DOUBLE_SUPPORT
  template class Identity<double>;
  template class ToneBurst<double>;
  template class HanningWeightedPulse<double>;
#endif
  template class Identity<float>;
  template class ToneBurst<float>;
  template class HanningWeightedPulse<float>;
}
