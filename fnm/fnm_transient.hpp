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
  /*! \brief ToneBurst
   *
   * @tparam T floating point type
   *
   * Tone burst separated using time-space decomposition for solving PDEs
   */
  template<class T>
  struct ToneBurst {
    static const size_t nTerms = 2;

    STATIC_INLINE_BEGIN T eval(T t, T W, T f0) STATIC_INLINE_END {
      return (fabs(t / W - T(0.5)) < T(0.5)) ? sin(T(M_2PI)*f0*t) : T(0.0);
    }
    STATIC_INLINE_BEGIN T sbf0(T tau, T W, T f0) STATIC_INLINE_END {
      SPS_UNREFERENCED_PARAMETER(W);
      return cos(T(M_2PI)*f0*tau);
    }
    STATIC_INLINE_BEGIN T sbf1(T tau, T W, T f0) STATIC_INLINE_END {
      SPS_UNREFERENCED_PARAMETER(W);
      return sin(T(M_2PI)*f0*tau);
    }
    STATIC_INLINE_BEGIN T tbf0(T t, T W, T f0) STATIC_INLINE_END {
      SPS_UNREFERENCED_PARAMETER(W);
      return sin(T(M_2PI)*f0*t);
    }
    STATIC_INLINE_BEGIN T tbf1(T t, T W, T f0) STATIC_INLINE_END {
      SPS_UNREFERENCED_PARAMETER(W);
      return -cos(T(M_2PI)*f0*t);
    }
    static T(*const SpatialBasisFunction[nTerms])(T, T, T);
    static T(*const TemporalBasisFunction[nTerms])(T, T, T);
  };

  // Not the derivative. We could compute a strip otherwise, we need to differentiate at the end

  /*! \brief HanningWeightedPulse
   *
   * @tparam T floating point type
   *
   * Hanning-weighted pulse separated using time-space decomposition for solving PDEs
   */
  template<class T>
  struct HanningWeightedPulse {
    static const size_t nTerms = 6;

    STATIC_INLINE_BEGIN T eval(T t, T W, T f0) STATIC_INLINE_END {
      return (fabs(t / W - T(0.5)) < T(0.5)) ?
      T(0.5)*(T(1.0) - cos(T(M_2PI)*t / W)) * sin(T(M_2PI)*f0*t) : T(0.0);
    }

    STATIC_INLINE_BEGIN T sbf0(T tau, T W, T f0) STATIC_INLINE_END {
      SPS_UNREFERENCED_PARAMETER(W);
      return cos(T(M_2PI)*f0*tau);
    }
    STATIC_INLINE_BEGIN T sbf1(T tau, T W, T f0) STATIC_INLINE_END {
      SPS_UNREFERENCED_PARAMETER(W);
      return sin(T(M_2PI)*f0*tau);
    }
    STATIC_INLINE_BEGIN T sbf2(T tau, T W, T f0) STATIC_INLINE_END {
      return cos(T(M_2PI)*tau / W)*cos(T(M_2PI)*f0*tau);
    }
    STATIC_INLINE_BEGIN T sbf3(T tau, T W, T f0) STATIC_INLINE_END {
      return cos(T(M_2PI)*tau / W)*sin(T(M_2PI)*f0*tau);
    }
    STATIC_INLINE_BEGIN T sbf4(T tau, T W, T f0) STATIC_INLINE_END {
      return sin(T(M_2PI)*tau / W)*cos(T(M_2PI)*f0*tau);
    }
    STATIC_INLINE_BEGIN T sbf5(T tau, T W, T f0) STATIC_INLINE_END {
      return sin(T(M_2PI)*tau / W)*sin(T(M_2PI)*f0*tau);
    }
    STATIC_INLINE_BEGIN T tbf0(T t, T W, T f0) STATIC_INLINE_END {
      SPS_UNREFERENCED_PARAMETER(W);
      return T(0.5)*sin(T(M_2PI)*f0*t);
    }
    STATIC_INLINE_BEGIN T tbf1(T t, T W, T f0) STATIC_INLINE_END {
      SPS_UNREFERENCED_PARAMETER(W);
      return -T(0.5)*cos(T(M_2PI)*f0*t);
    }
    STATIC_INLINE_BEGIN T tbf2(T t, T W, T f0) STATIC_INLINE_END {
      return -T(0.5)*cos(T(M_2PI)*t / W)*sin(T(M_2PI)*f0*t);
    }
    STATIC_INLINE_BEGIN T tbf3(T t, T W, T f0) STATIC_INLINE_END {
      return T(0.5)*cos(T(M_2PI)*t / W)*cos(T(M_2PI)*f0*t);
    }
    STATIC_INLINE_BEGIN T tbf4(T t, T W, T f0) STATIC_INLINE_END {
      return -T(0.5)*sin(T(M_2PI)*t / W)*sin(T(M_2PI)*f0*t);
    }
    STATIC_INLINE_BEGIN T tbf5(T t, T W, T f0) STATIC_INLINE_END {
      return T(0.5)*sin(T(M_2PI)*t / W)*cos(T(M_2PI)*f0*t);
    }
    static T(*const SpatialBasisFunction[nTerms])(T, T, T);
    static T(*const TemporalBasisFunction[nTerms])(T, T, T);
  };

  template <class T>
  T TransientSingleRect(const sysparm_t<T>* sysparm,
                        const ApertureData<T>* data,
                        const T* pos, const size_t nPositions,
                        T** odata, size_t* nSamples);

  template <class T>
  T CalcFdTransientRef(const sysparm_t<T>* sysparm,
                       const ApertureData<T>* data,
                       const T* pos, const size_t nPositions, const size_t nDim,
                       T** odata, size_t* nSignals, size_t* nSamples);
}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
