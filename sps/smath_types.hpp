#pragma once


#ifdef _MSC_VER
// TODO: Try out: alignas (16) std::array<char, 10>
# include <sps/aligned_array.hpp>
# include <sps/memory>
#else
# include <array>
#endif

#ifdef __GNUG__
// Check C++11
# include <cstdalign>
#endif

#include <limits>

#include <cstring>   // memcpy

#ifdef _WIN32
namespace std {
  template class SPS_EXPORT aligned_array<float,4>;
  template class SPS_EXPORT aligned_array<double,4>;
}
#endif

/** @addtogroup SPS */
/*@{*/

namespace sps {

  /*! \brief Point type (aligned)
   *
   *
   *  Point as an aligned array
   */
#ifdef _MSC_VER
  template <typename T>
  struct SPS_EXPORT point_t : public std::aligned_array<T,4> {};
#else

  //template<class T>
  //using point_t alignas(16) = std::array<T,4>;
  template <typename T>
  ALIGN16_BEGIN struct ALIGN16_END SPS_EXPORT point_t : public std::array<T,4> {};
#endif

  /*! \brief Rectangle
   *
   *
   *  Rectangle structure as an array of \ref sps::point_t
   */
#ifdef _MSC_VER
  template <typename T>
  struct SPS_EXPORT rect_t : public std::aligned_array<point_t<T> ,4> {};
#else
  template <typename T>
  struct rect_t : public std::array<point_t<T> ,4> {};
#endif


  /*! \brief Bounding box
   *
   *
   *  The bounding box is used for fast parallel computation and to
   *  optimize the cache usage.
   */
  template <typename T>
  struct SPS_EXPORT bbox_t {
    /// Boundary minimum
    point_t<T> min;
    /// Boundary maximum
    point_t<T> max;
  };

#ifdef _MSC_VER

#pragma warning(push)
#pragma warning(disable : 4324)
  ALIGN16_BEGIN struct ALIGN16_END sps_base_align16 {
    /*
      Empty classes are "special" because the standard requires that no complete object have size 0,
      so a class with no data members or bases always requires some padding if it doesn't have a vtable.
     */
    /*
     Apparently, MSVC creates a struct of size 16, if a struct of size 16 is derived from a struct of
     size 0, so we should remove this and ignore warnings instead.

    char dummy[16];

    */
  };
  ALIGN32_BEGIN struct ALIGN32_END sps_base_align32 {
    /*
     Apparently, MSVC creates a struct of size 16, if a struct of size 16 is derived from a struct of
     size 0, so we should remove this and ignore warnings instead.

    char dummy[32];

    */
  };
#pragma warning(pop)

  /*! \brief Euler angles
   *
   *
   *  Euler angles, the z-x-z' convention is used.
   */
  template <typename T>
  struct SPS_EXPORT euler_t : sps_base_align16 {
    /// alpha
    T alpha;
    /// beta
    T beta;
    /// gamma
    T gamma;
    /// Dummy used by alignment
    T dummy;
  };
#else
  template <typename T>
  ALIGN16_BEGIN struct ALIGN16_END euler_t {
    /// alpha
    T alpha;
    /// beta
    T beta;
    /// gamma
    T gamma;
    /// Dummy used by alignment
    T dummy;
  };
#endif

  /*! \brief Element representation
   *
   *  The element representation is used no matter if the
   *  PP(Plane-parallel)-wave approximation is used or the more
   *  accurate analytical responses.
   */
#if defined(_WIN32)
  template <typename T>
  struct SPS_EXPORT element_t : sps::win32::base_align16
#else
  template <typename T>
  ALIGN16_BEGIN struct SPS_EXPORT ALIGN16_END element_t
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
}
/*@}*/
