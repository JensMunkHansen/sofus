/**
 * @file   fnm_types.hpp
 * @author Jens Munk Hansen <jmh@jmhlaptop>
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

namespace fnm {

#if 0
  /*! \brief Element representation
   *
   *  The element representation is used no matter if the
   *  PP(Plane-parallel)-wave approximation is used or the more
   *  accurate analytical responses.
   */
#if defined(_WIN32)
  template <typename T>
  struct FNM_EXPORT element_t : sps::win32::base_align16
#else
  template <typename T>
  ALIGN16_BEGIN struct FNM_EXPORT ALIGN16_END element_t
#endif
  {
    ALIGN16_BEGIN sps::point_t<T> center ALIGN16_END;
    ALIGN16_BEGIN sps::euler_t<T> euler ALIGN16_END;
    ALIGN16_BEGIN T hw ALIGN16_END;  ///< Half width
    ALIGN16_BEGIN T hh ALIGN16_END;  ///< Half height
    ALIGN16_BEGIN T dummy[2] ALIGN16_END;

    //@{ Cached variables.

    /*
      There is no need for specializing using either __m128 or
      __m256. The data need to be loaded into a register in all cases.
    */
    ALIGN16_BEGIN T normal[4] ALIGN16_END;
    ALIGN16_BEGIN T uvector[4] ALIGN16_END;
    ALIGN16_BEGIN T vvector[4] ALIGN16_END;
    ALIGN16_BEGIN T vertices[3][4] ALIGN16_END;
    //@}
  };
#endif

  /*! \brief Sysparm structure
   *
   *
   * A structure containing global simulation parameters
   */
  template <class T>
  struct FNM_EXPORT sysparm_t {
    /// Speed of sound
    T c;
    /// Number of width abcissas
    size_t nDivW;
    /// Number of height abcissas
    size_t nDivH;
  };

  struct FNM_EXPORT FocusingTypeNS {
    enum Value {
      Rayleigh          = 0x00, ///< Rayleigh integral is solved to fix phase of complex signal
      Pythagorean       = 0x01, ///< Distance from center of element is used to fix phase
      FocusingTypeCount = 0x02,
    };
  };
  typedef FocusingTypeNS::Value FocusingType;
}
