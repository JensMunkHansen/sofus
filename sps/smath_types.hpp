/**
 * @file   smath_types.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Sat Mar 18 13:55:23 2017
 *
 * @brief  Mathematical structures
 *
 *
 */
#pragma once

#include <sps/sps_export.h>

#ifdef _MSC_VER
# include <sps/aligned_array.hpp>
#else
# include <array>
#endif

#include <sps/memory> // For sps::aligned<>

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
  template <typename T>
  // Must be 32-byte aligned for double. TODO: Make alignment depend on T
#if 1 /* DOUBLE SUPPORT */
struct __attribute__((aligned(4*sizeof(T)))) SPS_EXPORT point_t : public std::array<T,4> {
#else
struct __attribute__((aligned(16))) SPS_EXPORT point_t : public std::array<T,4> {
#endif
    point_t()
  {
    (*this)[0] = (*this)[1] = (*this)[2] = (*this)[3] = T(0.0);
  }
  };
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

  /*! \brief Euler angles
   *
   *
   *  Euler angles, the z-x-z' convention is used.
   */
  template <typename T>
  struct SPS_EXPORT euler_t : sps::aligned<4*sizeof(T)> {
    /// alpha
    T alpha;
    /// beta
    T beta;
    /// gamma
    T gamma;
    /// Dummy used by alignment
    T dummy;
  };

  /*! \brief Element representation
   *
   *  The element representation is used no matter if the
   *  PP(Plane-parallel)-wave approximation is used or the more
   *  accurate analytical responses.
   */

  // TODO: Consider inheriting from sps::aligned<4*sizeof(T)>
#if defined(_WIN32)
# if 0 /* DOUBLE SUPPORT */
  template <typename T>
  struct SPS_EXPORT element_t : sps::win32::base_align32
# else
  template <typename T>
  struct SPS_EXPORT element_t : sps::win32::base_align16
# endif
#else
  template <typename T>
  struct SPS_EXPORT __attribute__((aligned(4*sizeof(T)))) element_t
#endif
  {

#ifndef _MSC_VER /* DOUBLE SUPPORT */
    ALIGN32_BEGIN sps::point_t<T> center ALIGN32_END; ///< Center position
    ALIGN32_BEGIN sps::euler_t<T> euler ALIGN32_END;  ///< Euler angles
#else
    ALIGN16_BEGIN sps::point_t<T> center ALIGN16_END; ///< Center position
    ALIGN16_BEGIN sps::euler_t<T> euler ALIGN16_END;  ///< Euler angles
#endif
    T hw;       ///< Half width
    T hh;       ///< Half height
    T dummy[2]; ///< Dummy for alignment

    /** @name Cached variables
     *
     */
    ///@{

    /// Normal vector
    ALIGN16_BEGIN T normal[4] ALIGN16_END;
    /// First basis vector
    ALIGN16_BEGIN T uvector[4] ALIGN16_END;
    /// Second basis vector
    ALIGN16_BEGIN T vvector[4] ALIGN16_END;
    /// Four vertices, xyz major-order
    ALIGN16_BEGIN T vertices[3][4] ALIGN16_END;

    ///@}
  };
#ifdef _MSC_VER
  template struct SPS_EXPORT bbox_t<float>;
  template struct SPS_EXPORT bbox_t<double>;
#endif
}
/*@}*/

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
