/**
 * @file   fnm_basis.hpp
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
#pragma once

#include <sps/cenv.h>
#include <sps/math.h>
#include <sps/cmath>
#include <cstddef>
#include <algorithm> // for std::fill

// TODO: Vectorize all basis functions and evaluate

namespace fnm {
  /*! \brief ToneBurst
   *
   * @tparam T floating point type
   *
   * Tone burst separated using time-space decomposition for solving PDEs
   */
  template<class T>
  class ToneBurst {
  public:
    static const size_t nTerms = 2;

    STATIC_INLINE_BEGIN T Evaluate(T t, T W, T f0) STATIC_INLINE_END {
      return (fabs(t / W - T(0.5)) < T(0.5)) ? sps::sin<T>(T(M_2PI)*f0*t) : T(0.0);
    }
    STATIC_INLINE_BEGIN T sbf0(T tau, T W, T f0) STATIC_INLINE_END {
      SPS_UNREFERENCED_PARAMETER(W);
      return sps::cos<T>(T(M_2PI)*f0*tau);
    }
    STATIC_INLINE_BEGIN T sbf1(T tau, T W, T f0) STATIC_INLINE_END {
      SPS_UNREFERENCED_PARAMETER(W);
      return sps::sin<T>(T(M_2PI)*f0*tau);
    }
    STATIC_INLINE_BEGIN T tbf0(T t, T W, T f0) STATIC_INLINE_END {
      SPS_UNREFERENCED_PARAMETER(W);
      return sps::sin<T>(T(M_2PI)*f0*t);
    }
    STATIC_INLINE_BEGIN T tbf1(T t, T W, T f0) STATIC_INLINE_END {
      SPS_UNREFERENCED_PARAMETER(W);
      return -sps::cos<T>(T(M_2PI)*f0*t);
    }

    static T(*const SpatialBasisFunction[nTerms])(T, T, T);
    static T(*const TemporalBasisFunction[nTerms])(T, T, T);

    /** @name Introducing state for evaluation
     *
     */
    ///@{

    ToneBurst();

    inline void ResetSpatial();

    // TODO: Check if it is faster to use a static function add argument T (coeff*)[nTerms]
    inline void UpdateSpatial(T factor, T tau, T W, T f0);

    inline T EvaluateTSD(T t, T W, T f0) const
    {
      T result = T(0.0);
      for (size_t iTerm = 0 ; iTerm < ToneBurst<T>::nTerms ; iTerm++) {
        result += this->m_fTerms[iTerm] * ToneBurst<T>::TemporalBasisFunction[iTerm](t,W,f0);
      }
      return result;
    }
    ///@}

  private:
    T m_fTerms[nTerms];
  };

  template<class T>
  class Identity {
  public:
    static const size_t nTerms = 1;

    STATIC_INLINE_BEGIN T Evaluate(T t, T W, T f0) STATIC_INLINE_END {
      SPS_UNREFERENCED_PARAMETER(f0);
      return (fabs(t / W - T(0.5)) < T(0.5)) ? T(1.0) : T(0.0);
    }
    STATIC_INLINE_BEGIN T sbf0(T tau, T W, T f0) STATIC_INLINE_END {
      SPS_UNREFERENCED_PARAMETERS(tau, W, f0);
      return T(1.0);
    }
    STATIC_INLINE_BEGIN T tbf0(T t, T W, T f0) STATIC_INLINE_END {
      SPS_UNREFERENCED_PARAMETERS(t, W, f0);
      return T(1.0);
    }
    static T(*const SpatialBasisFunction[nTerms])(T, T, T);
    static T(*const TemporalBasisFunction[nTerms])(T, T, T);

    /** @name Introducing state for evaluation
     *
     */
    ///@{

    Identity();

    inline void ResetSpatial();

    inline void UpdateSpatial(T factor, T tau, T W, T f0);

    inline T EvaluateTSD(T t, T W, T f0) const
    {
      T result = T(0.0);
      for (size_t iTerm = 0 ; iTerm < Identity<T>::nTerms ; iTerm++) {
        result += this->m_fTerms[iTerm] * Identity<T>::TemporalBasisFunction[iTerm](t,W,f0);
      }
      return result;
    }
    ///@}
  private:
    T m_fTerms[nTerms];
  };


  // Not the derivative. We could compute a strip otherwise, we need to differentiate at the end

  /*! \brief HanningWeightedPulse
   *
   * @tparam T floating point type
   *
   * Hanning-weighted pulse separated using time-space decomposition for solving PDEs
   */
  template<class T>
  class HanningWeightedPulse {
  public:
    static const size_t nTerms = 6;

    STATIC_INLINE_BEGIN T Evaluate(T t, T W, T f0) STATIC_INLINE_END {
      return (fabs(t / W - T(0.5)) < T(0.5)) ?
      T(0.5)*(T(1.0) - sps::cos<T>(T(M_2PI)*t / W)) * sps::sin<T>(T(M_2PI)*f0*t) : T(0.0);
    }

    STATIC_INLINE_BEGIN T sbf0(T tau, T W, T f0) STATIC_INLINE_END {
      SPS_UNREFERENCED_PARAMETER(W);
      return sps::cos<T>(T(M_2PI)*f0*tau);
    }
    STATIC_INLINE_BEGIN T sbf1(T tau, T W, T f0) STATIC_INLINE_END {
      SPS_UNREFERENCED_PARAMETER(W);
      return sps::sin<T>(T(M_2PI)*f0*tau);
    }
    STATIC_INLINE_BEGIN T sbf2(T tau, T W, T f0) STATIC_INLINE_END {
      return sps::cos<T>(T(M_2PI)*tau / W)*sps::cos<T>(T(M_2PI)*f0*tau);
    }
    STATIC_INLINE_BEGIN T sbf3(T tau, T W, T f0) STATIC_INLINE_END {
      return sps::cos<T>(T(M_2PI)*tau / W)*sps::sin<T>(T(M_2PI)*f0*tau);
    }
    STATIC_INLINE_BEGIN T sbf4(T tau, T W, T f0) STATIC_INLINE_END {
      return sps::sin<T>(T(M_2PI)*tau / W)*sps::cos<T>(T(M_2PI)*f0*tau);
    }
    STATIC_INLINE_BEGIN T sbf5(T tau, T W, T f0) STATIC_INLINE_END {
      return sps::sin<T>(T(M_2PI)*tau / W)*sps::sin<T>(T(M_2PI)*f0*tau);
    }
    STATIC_INLINE_BEGIN T tbf0(T t, T W, T f0) STATIC_INLINE_END {
      SPS_UNREFERENCED_PARAMETER(W);
      return T(0.5)*sps::sin<T>(T(M_2PI)*f0*t);
    }
    STATIC_INLINE_BEGIN T tbf1(T t, T W, T f0) STATIC_INLINE_END {
      SPS_UNREFERENCED_PARAMETER(W);
      return -T(0.5)*sps::cos<T>(T(M_2PI)*f0*t);
    }
    STATIC_INLINE_BEGIN T tbf2(T t, T W, T f0) STATIC_INLINE_END {
      return -T(0.5)*sps::cos<T>(T(M_2PI)*t / W)*sps::sin<T>(T(M_2PI)*f0*t);
    }
    STATIC_INLINE_BEGIN T tbf3(T t, T W, T f0) STATIC_INLINE_END {
      return T(0.5)*sps::cos<T>(T(M_2PI)*t / W)*sps::cos<T>(T(M_2PI)*f0*t);
    }
    STATIC_INLINE_BEGIN T tbf4(T t, T W, T f0) STATIC_INLINE_END {
      return -T(0.5)*sps::sin<T>(T(M_2PI)*t / W)*sps::sin<T>(T(M_2PI)*f0*t);
    }
    STATIC_INLINE_BEGIN T tbf5(T t, T W, T f0) STATIC_INLINE_END {
      return T(0.5)*sps::sin<T>(T(M_2PI)*t / W)*sps::cos<T>(T(M_2PI)*f0*t);
    }
    static T(*const SpatialBasisFunction[nTerms])(T, T, T);
    static T(*const TemporalBasisFunction[nTerms])(T, T, T);

    /** @name Introducing state for evaluation
     *
     */
    ///@{

    HanningWeightedPulse();

    inline void ResetSpatial();

    inline void UpdateSpatial(T factor, T tau, T W, T f0);

    inline T EvaluateTSD(T t, T W, T f0) const
    {
      T result = T(0.0);
      for (size_t iTerm = 0 ; iTerm < HanningWeightedPulse<T>::nTerms ; iTerm++) {
        result += this->m_fTerms[iTerm] * HanningWeightedPulse<T>::TemporalBasisFunction[iTerm](t,W,f0);
      }
      return result;
    }
    ///@}
  private:
    T m_fTerms[nTerms];
  };
}
/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
