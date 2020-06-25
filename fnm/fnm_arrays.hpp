/**
 * @file   fnm_arrays.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Sat Mar 18 13:17:09 2017
 *
 * @brief  Geometries for transducer arrays
 *
 * Copyright 2017 Jens Munk Hansen
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

#include <sps/memory>
#include <cstdlib>
#include <sps/smath.hpp>

#include <fnm/fnm_types.hpp>
#include <fnm/fnm_alias.hpp>

namespace fnm {

/**
 * Create Matrix array
 *
 * @param[in] nRows
 * @param[in] nCols
 * @param[in] rowWidth
 * @param[in] rowKerf
 * @param[in] colWidth
 * @param[in] colKerf
 * @param[out] elements
 *
 * @return
 */
template <class T>
void MatrixArray(
  const size_t nRows,
  const size_t nCols,
  const T rowWidth,
  const T rowKerf,
  const T colWidth,
  const T colKerf,
  element_array<T> &&elements);


/**
 * Compute focused linear array
 *
 * @param[in] nElements
 * @param[in] nSubH
 * @param[in] nSubW
 * @param[in] width
 * @param[in] kerf
 * @param[in] height
 * @param[in] eFocus
 * @param[in] arcPlacement 0 = outside, 1 = inside
 * @param[out] elements
 */
template <class T>
void FocusedLinearArray(
  const size_t nElements,
  const size_t nSubH,  // Elevation
  const size_t nSubW,
  const T width,
  const T kerf,
  const T height,
  const T eFocus,
  const int arcPlacement,
  element_array<T> &&elements);

/**
 * Compute focused convex array
 *
 * @param[in] nElements
 * @param[in] nSubH        sub-elements: elevation
 * @param[in] nSubW        sub-elements: azimuth
 * @param[in] width
 * @param[in] kerf
 * @param[in] height
 * @param[in] radius
 * @param[in] eFocus       elevation focus
 * @param[in] arcPlacement 0 = outside, 1 = inside
 * @param[out] elements
 *
 * @return
 */
template <class T>
void FocusedConvexArray(
  const size_t nElements,
  const size_t nSubH,  // Elevation
  const size_t nSubW,
  const T width,
  const T kerf,
  const T height,
  const T radius,
  const T eFocus,
  const int arcPlacement,
  element_array<T> &&elements);

}  // namespace fnm
/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
