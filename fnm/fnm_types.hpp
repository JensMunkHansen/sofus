/**
 * @file   fnm_types.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Tue Oct 18 20:17:24 2016
 *
 * @brief
 *
 *
 */

#pragma once

#include <fnm/config.h>
#include <fnm/fnm_export.h>

#include <sps/cenv.h>
#include <sps/memory>
#include <sps/smath.hpp> // sps::point_t and sps::euler_t, now element_t

#if FNM_PULSED_WAVE
# include <sofus/sofus_types.hpp>
#endif
namespace fnm {

  /*! \brief Sysparm structure
   *
   *
   * A structure containing global simulation parameters.
   *
   * TODO: Separate into frequency domain and time domain parameters
   */
  template <class T>
  struct FNM_EXPORT sysparm_t {
    /// Speed of sound
    T c;

    /// Number of width abcissas
    size_t nDivW;
    /// Number of height abcissas
    size_t nDivH;

    T att;
    T beta;
    bool use_att;

#if FNM_PULSED_WAVE
    /// Time-domain calculation type
    sofus::TimeDomainCalcType timeDomainCalcType;

    /// Bum
    sofus::PulsedWaveIntOrder pulseWaveIntOrder;
#endif
  };

  struct FNM_EXPORT FocusingTypeNS {
    enum Value {
      Rayleigh          = 0x00, ///< Rayleigh integral is solved to fix phase of complex signal
      Pythagorean       = 0x01, ///< Distance from center of element is used to fix phase
      Delays            = 0x02,
      FocusingTypeCount = 0x03, ///< Used for invalid focusing type
    };
  };

  // TODO: Renaming to time-domain / freq. domain, when Goertzel is implemeted and we can do
  // time-space decomposition.
  typedef FocusingTypeNS::Value FocusingType;

}
