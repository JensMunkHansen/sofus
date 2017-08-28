/**
 * @file   smath_types.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Sat Mar 18 13:55:23 2017
 *
 * @brief  Mathematical structures
 *
 *
 */
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
  struct SPS_EXPORT point_t : public std::aligned_array<T,4>
#else
  // Must be 32-byte aligned for double.
  template <typename T>
struct __attribute__((aligned(4*sizeof(T)))) SPS_EXPORT point_t : public std::array<T,4>
#endif
  {
    static const point_t xaxis;
    static const point_t yaxis;
    static const point_t zaxis;
    point_t() = default;
    // Neede for initialization of constant vectors
    point_t(const T& a, const T& b, const T&c)
    {
      (*this)[0] = a;
      (*this)[1] = b;
      (*this)[2] = c;
    }
  };

  template <typename T>
  const point_t<T> point_t<T>::xaxis = point_t<T>(T(1.0),T(0.0),T(0.0));

  template <typename T>
  const point_t<T> point_t<T>::yaxis = point_t<T>(T(0.0),T(1.0),T(0.0));

  template <typename T>
  const point_t<T> point_t<T>::zaxis = point_t<T>(T(0.0),T(0.0),T(1.0));

  /*! \brief Vector aliased as point
   *
   */
  template <typename T>
  using vector_t = point_t<T>;

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
struct __attribute__((aligned(4*sizeof(T)))) SPS_EXPORT rect_t : public std::array<point_t<T> ,4> {};
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

  template <typename T>
  struct SPS_EXPORT circle_t : sps::aligned<4*sizeof(T)> {
    sps::point_t<T> center; ///< Center position
    sps::euler_t<T> euler;  ///< Euler angles
    T radius;               ///< Radius
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
  struct SPS_EXPORT element_rect_t : sps::win32::base_align32
# else
  template <typename T>
  struct SPS_EXPORT element_rect_t : sps::win32::base_align16
# endif
#else
  template <typename T>
  struct SPS_EXPORT __attribute__((aligned(4*sizeof(T)))) element_rect_t
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

  template <typename T>
  struct SPS_EXPORT element_circular_t : sps::aligned<4*sizeof(T)> {
    sps::circle_t<T> circle; // Circular element

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
