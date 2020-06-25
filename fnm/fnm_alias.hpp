/**
 * @file   fnm_alias.hpp
 * @author Jens Munk Hansen <jmh@debian9laptop.parknet.dk>
 * @date   Thu May 24 19:44:02 2018
 *
 * @brief
 *
 *
 */
#pragma once

#include <sps/smath.hpp>

namespace fnm {
/** Type-alias for a unique two-dimensional array of elements */
template <class T>
using element_array =
  sps::unique_aligned_multi_array<sps::element_rect_t<T>, 2>;

template <class T>
using point_array =
  sps::unique_aligned_array<sps::point_t<T> >;

// Function alias used in this compilation unit
template <class T>
constexpr element_array<T>(*element_array_create)(size_t m, size_t n) =
  &sps::unique_aligned_multi_array_create<sps::element_rect_t<T>, 2>;

template <class T>
constexpr element_array<T>(*point_array_create)(size_t m) =
  &sps::unique_aligned_array_create<sps::point_t<T> >;

}  // namespace fnm

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
