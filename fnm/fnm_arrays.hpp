/**
 * @file   fnm_arrays.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Sat Mar 18 13:17:09 2017
 *
 * @brief  Geometries for transducer arrays
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

#include <cstdlib>

#include <sps/memory>
#include <sps/smath.hpp>

namespace fnm {

  /**
   * Create Matrix array
   *
   * @param nRows
   * @param nCols
   * @param rowWidth
   * @param rowKerf
   * @param colWidth
   * @param colKerf
   * @param elements     output
   *
   * @return
   */
  template <class T>
  void MatrixArray(const size_t nRows,
                   const size_t nCols,
                   const T rowWidth,
                   const T rowKerf,
                   const T colWidth,
                   const T colKerf,
                   sps::deleted_aligned_multi_array<sps::element_t<T>, 2> &&elements);


  /**
   * Compute focused linear array
   *
   * @param nElements
   * @param nSubH
   * @param nSubW
   * @param width
   * @param kerf
   * @param height
   * @param eFocus
   * @param arcPlacement 0 = outside, 1 = inside
   * @param elements     output
   */
  template <class T>
  void FocusedLinearArray(const size_t nElements,
                          const size_t nSubH, // Elevation
                          const size_t nSubW,
                          const T width,
                          const T kerf,
                          const T height,
                          const T eFocus,
                          const int arcPlacement,
                          sps::deleted_aligned_multi_array<sps::element_t<T>, 2> &&elements);

  /**
   * Compute focused convex array
   *
   * @param nElements
   * @param nSubH        sub-elements: elevation
   * @param nSubW        sub-elements: azimuth
   * @param width
   * @param kerf
   * @param height
   * @param radius
   * @param eFocus       elevation focus
   * @param arcPlacement 0 = outside, 1 = inside
   * @param elements     output
   *
   * @return
   */
  template <class T>
  void FocusedConvexArray(const size_t nElements,
                          const size_t nSubH, // Elevation
                          const size_t nSubW,
                          const T width,
                          const T kerf,
                          const T height,
                          const T radius,
                          const T eFocus,
                          const int arcPlacement,
                          sps::deleted_aligned_multi_array<sps::element_t<T>, 2> &&elements);

}
/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
