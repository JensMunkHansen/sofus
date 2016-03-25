#pragma once

#include <sps/smath.hpp>
#include <fnm/fnm_export.h>

namespace fnm {
  /*! \brief Element representation
   *
   *  The element representation is used no matter if the
   *  PP(Plane-parallel)-wave approximation is used or the more
   *  accurate analytical responses.
   */
  #ifdef _MSC_VER
  // Not aligned
  template <typename T>
  struct FNM_EXPORT element_t
  #else
  // There is no explicit requirement for alignment of the object
  template <typename T>
  ALIGN16_BEGIN struct ALIGN16_END element_t
  #endif
  {
  // There is no explicit requirement for alignment of the members
    ALIGN16_BEGIN sps::point_t<T> center ALIGN16_END;
    ALIGN16_BEGIN sps::euler_t<T> euler ALIGN16_END;
    ALIGN16_BEGIN T hw ALIGN16_END;  /* Half width  */
    ALIGN16_BEGIN T hh ALIGN16_END;  /* Half height */
    ALIGN16_BEGIN T dummy[2] ALIGN16_END;
  };

  template <typename T>
  struct FNM_EXPORT sysparm_t {
    T c;
    size_t nDivW;
    size_t nDivH;
  };
}
