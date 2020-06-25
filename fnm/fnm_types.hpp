/**
 * @file   fnm_types.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Tue Oct 18 20:17:24 2016
 *
 * @brief Structure types used for field simulations
 *
 * Copyright 2017 Jens Munk Hansen
 */

#pragma once

#include <fnm/config.h>
#include <fnm/fnm_export.h>
#include <fnm/fnm_types.h>

#include <sps/cenv.h>

#if FNM_PULSED_WAVE
# include <sofus/sofus_types.hpp>
#endif

#include <sps/globals.hpp>

#include <cstddef>

namespace fnm {

// TODO(JMH): Use this construct for handling system parameters, i.e. sysparm_t
template <class T>
struct FNM_EXPORT bla : public sps::Default<bla<T> > {
  T k;
  bla();
 private:
  //bla& operator=(const bla& rhs) = delete;
  friend class sps::Default<bla<T> >;
};

/*! \brief Sysparm structure
 *
 *
 * A structure containing global simulation parameters.
 *
 * TODO: Separate into frequency domain and time domain parameters
 */
template <class T>
struct FNM_EXPORT sysparm_t {
  sysparm_t();

  /// Speed of sound
  T c;

  /// Sample frequency
  T fs;

  /// Center frequency (for temporal response)
  T f0;

  /// Frequency dependent attenuation
  T att;

  /// Frequency independent attenuation
  T beta;

  /// Use attenuation
  bool use_att;

  /// Density [kg/m3]
  T rho;

  /// Number of threads
  size_t nThreads;

  /** @name Frequency domain parameters
   *
   */
  ///@{

  /// Number of width abcissas
  size_t nDivW;

  /// Number of height abcissas
  size_t nDivH;

  /// Number of angular abcissas
  size_t nDivA;

  /// Maximum number of angular sectors
  size_t nMaxSectors;

  /// Width of Hanning-weighted sinusoid pulse
  T w;

  ///@}

#if FNM_PULSED_WAVE
  /** @name Time domain parameters
   *
   */
  ///@{

  /// Time-domain calculation type
  sofus::TimeDomainCalcType timeDomainCalcType;

  /// Time-domain integration order
  sofus::TimeDomainIntOrder timeDomainIntOrder;

  /// Enforce positive delays, when derived from a focus point
  bool enforcePositiveDelays;

  /// Soft-baffle
  bool soft_baffle;
  ///@}
#endif
};

typedef FNM_ExcitationTypeNS::ExcitationType_Value ExcitationType;

typedef FNM_FocusingTypeNS::FocusingType_Value FocusingType;

typedef FNM_ApodizationTypeNS::ApodizationType_Value ApodizationType;

typedef RwParamTypeNS::RwParamType_Value RwParamType;

typedef FNM_TypeNS::Type_Value Type;

/*! \brief GL (Gauss-Legendre) structure
 *
 * Weight and abcissas for 1D integration
 */
template <class T>
struct GLQuad1D {
  T* xs;       ///< abcissas
  T* ws;       ///< weights
  size_t nx;   ///< number of abcissas
};

/*! \brief GL (Gauss-Legendre) structure
 *
 * Weight and abcissas for 2D integration
 */
template <class T>
struct GLQuad2D {
  GLQuad1D<T> u;
  GLQuad1D<T> v;
};
}  // namespace fnm
/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
