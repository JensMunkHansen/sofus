/**
 * @file   fnm_delays.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Wed Apr  4 22:12:35 2018
 *
 * @brief
 *
 *
 */
#pragma once

#include <sps/cenv.h>
#include <sps/smath_types.hpp>

namespace fnm {
template <class T>
int PythagoreanFocusDelays(
  const sps::point_t<T>* pElementPositions, const size_t& nElements,
  const sps::point_t<T>& focus, const sps::point_t<T>& centerFocus,
  const size_t& nDelays, T* pDelays);


template <class T>
int PytDelays0(const sps::unique_aligned_array<sps::point_t<T> >& positions,
               const size_t& nElements,
               const sps::point_t<T>& centerFocus,
               const sps::point_t<T>& focus0,
               const sps::point_t<T>& focus1,
               const T& c,
               sps::unique_aligned_array<T>& delays, bool allPositive = false);

template <class T>
int PytDelays1(const sps::unique_aligned_array<sps::point_t<T> >& positions,
               const size_t& nElements,
               const sps::point_t<T>& centerFocus,
               const sps::point_t<T>& focus0,
               const sps::point_t<T>& focus1,
               const T& c,
               sps::unique_aligned_array<T>& delays, bool allPositive = false);

template <class T>
int PytDelays2(const sps::unique_aligned_array<sps::point_t<T> >& positions,
               const size_t& nElements,
               const sps::point_t<T>& centerFocus,
               const sps::point_t<T>& focus0,
               const sps::point_t<T>& focus1,
               const T& c,
               sps::unique_aligned_array<T>& delays, bool allPositive = false);

// The function PytDelays2 (Pythagorean Delays 2) is aliased as PytDelays
SPS_ALIAS_TEMPLATE_FUNCTION(PytDelays, PytDelays2)

}  // namespace fnm
