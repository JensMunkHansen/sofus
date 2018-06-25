#include <fnm/config.h>
#include <fnm/fnm_arrays.hpp>
#include <sps/smath.hpp>
#include <cassert>

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

namespace fnm {

  template <class T>
  void MatrixArray(const size_t nRows,
                   const size_t nCols,
                   const T rowWidth,
                   const T rowKerf,
                   const T colWidth,
                   const T colKerf,
                   sps::deleted_aligned_multi_array<sps::element_rect_t<T>, 2> &&elements)
  {
    elements      = sps::deleted_aligned_multi_array_create<sps::element_rect_t<T>,2>(nRows*nCols,1);
    T rowPitch = rowWidth + rowKerf;
    T colPitch = colWidth + colKerf;
    T hh = T(0.5) * colWidth;
    T hw = T(0.5) * rowWidth;

    for (size_t iRow=0 ; iRow < nRows ; iRow++) {
      sps::point_t<T> center = sps::point_t<T>(); // Only to avoid warnings
      center[1] = (T(iRow) - T(0.5)*(T(nRows)-T(1.0)))*rowPitch;
      center[2] = T(0.0);
      for (size_t iCol=0 ; iCol < nCols ; iCol++) {
        auto& element = elements[iRow*nCols + iCol][0];

        sps::point_t<T> center1 = center;
        center1[0] =  (T(iCol) - T(0.5)*(T(nCols)-T(1.0)))*colPitch;
        element.hh     = hh;
        element.hw     = hw;
        element.center = center1;
        element.euler.alpha = T(0.0);
        element.euler.beta  = T(0.0);
        element.euler.gamma = T(0.0);
      }
    }
  }

  template <class T>
  void FocusedLinearArray(const size_t nElements,
                          const size_t nSubH, // Elevation (minimum is 1)
                          const size_t nSubW, // Sub-elements in azimuth (=1)
                          const T width,
                          const T kerf,
                          const T height,
                          const T eFocus,
                          const int arcPlacement,
                          sps::deleted_aligned_multi_array<sps::element_rect_t<T>, 2> &&elements)
  {
    assert(nSubH>0);

    elements      = sps::deleted_aligned_multi_array_create<sps::element_rect_t<T>,2>(nElements,nSubH);

    T focus = eFocus;
    T halfHeight = T(0.5)*height;

    T R = sqrt(SQUARE(focus)+SQUARE(halfHeight));

    T elSector = T(2.0)*atan2(halfHeight,focus); /* equals pi if focus is 0.0 */

    T dEl = elSector / T(nSubH);
    T pitch = width + kerf;
    T hw  = T(0.5)*(pitch - kerf);

    T chordLength = T(2.0) * R * sin(T(0.5)*dEl);
    T tanLength;

    if (nSubH == 1) {
      /* We do not allow focusing using a single element - design decision */
      tanLength   = height;
      chordLength = height;
    } else {
      if (eFocus < std::numeric_limits<T>::epsilon()) {
        /* We can focus at 0.0, but decide not to. */
        /* This can be skipped if we allow focusing at 0.0 */
        dEl = T(0.0);
        tanLength   = height / T(nSubH);
        chordLength = height / T(nSubH);
      } else {
        /* Infinity, when focus is 0.0 and nSubH == 1 */
        tanLength = T(2.0) * R * tan(T(0.5)*dEl);
      }
    }

    T hh = T(0.0);

    if (arcPlacement == 0) {
      // Outer
      hh = T(0.5)*tanLength;
    } else {
      // Inner
      hh = T(0.5)*chordLength;
    }

    T elAngle     = T(0.0);
    T lastElAngle = T(0.0);

    // Get z-coordinate right, we cannot focus using a single element
    if (nSubH == 1) {
      R = focus;
    }

    // We decide not to allow focusing at 0.0
    if (eFocus < std::numeric_limits<T>::epsilon()) {
      /* This can be skipped if we allow focusing at 0.0 */
      R = focus;
    }

    if (arcPlacement == 0) {
      // Placement outside arc (default)
      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {

        for (size_t iSubH = 0 ; iSubH < nSubH ; iSubH++) {
          elAngle = (T(iSubH) - T(0.5)*(T(nSubH)-T(1.0))) * dEl;
          sps::point_t<T> center = sps::point_t<T>(); // Only to avoid warnings
          center[0] = (T(iElement) - T(0.5)*(T(nElements)-T(1.0)))*pitch;
          if (eFocus < std::numeric_limits<T>::epsilon()) {
            /* This can be skipped if we allow focusing at 0.0 */
            center[1] = (T(iSubH) - T(0.5)*(T(nSubH)-T(1.0))) * T(2.0) * hh;
          } else {
            center[1] = sin(elAngle)*R;
          }
          // For z-coordinate R must equal focus if nSub == 1
          center[2] = - (cos(elAngle)*R - focus)  +  (R-focus);/* Added + (R-focus) */
          for (size_t iSubW = 0 ; iSubW < nSubW ; iSubW++) {

            sps::point_t<T> center1 = center;
            center1[0] = center1[0] + T(2.0) * hw/T(nSubW) * (T(iSubW) - T(0.5)*(T(nSubW)-T(1.0)) );

            elements[iElement][iSubH].hh     = hh;
            elements[iElement][iSubH].hw     = hw / T(nSubW);
            elements[iElement][iSubH].center = center1;
            elements[iElement][iSubH].euler.alpha = T(0.0);
            elements[iElement][iSubH].euler.beta  = elAngle;
            elements[iElement][iSubH].euler.gamma = T(0.0);
          }
        }
      }
    } else {
      // Placement inside arc
      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
        for (size_t iSubH = 0 ; iSubH < nSubH ; iSubH++) {
          if (iSubH == 0) {
            lastElAngle = (T(iSubH) - T(0.5)*T(nSubH))*dEl;
          }
          elAngle = (T(iSubH+1) - T(0.5)*T(nSubH))*dEl;

          sps::point_t<T> center = sps::point_t<T>();
          center[0] = (T(iElement) - T(0.5)*(T(nElements)-T(1.0)))*pitch;
          if (eFocus < std::numeric_limits<T>::epsilon()) {
            /* This can be skipped if we allow focusing at 0.0 */
            center[1] = (T(iSubH) - T(0.5)*(T(nSubH)-T(1.0))) * T(2.0) * hh;
          } else {
            center[1] = T(0.5)*(sin(elAngle)+sin(lastElAngle))*R;
          }
          center[2] = -(T(0.5)*(cos(elAngle)+cos(lastElAngle))*R-focus);

          for (size_t iSubW = 0 ; iSubW < nSubW ; iSubW++) {
            sps::point_t<T> center1 = center;
            center1[0] = center1[0] + T(2.0)*hw/T(nSubW)*(T(iSubW)-T(0.5)*(T(nSubW)-T(1.0)));

            elements[iElement][iSubH].hh     = hh;
            elements[iElement][iSubH].hw     = hw/T(nSubW);
            elements[iElement][iSubH].center = center1;
            elements[iElement][iSubH].euler.alpha = T(0.0);
            elements[iElement][iSubH].euler.beta  = T(0.5)*(elAngle+lastElAngle);
            elements[iElement][iSubH].euler.gamma = T(0.0);
          }
          lastElAngle = elAngle;
        }
      }
    }
  }

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
                          sps::deleted_aligned_multi_array<sps::element_rect_t<T>, 2> &&elements)
  {
    SPS_UNREFERENCED_PARAMETER(nSubW);
    assert(nSubW==1);
    assert(arcPlacement==0);

    elements      = sps::deleted_aligned_multi_array_create<sps::element_rect_t<T>,2>(nElements,nSubH);

    T halfHeight = T(0.5)*height;
    T focus = eFocus;
    T pitch = width + kerf;
    if (focus < T(0.0))
      focus = T(0.0);

    T azR = radius;

    // Arc-length from center to center
    T azArcLength = pitch * (T(nElements)-T(1.0));

    T azSegment = azArcLength / azR;

    T dAz = azSegment / std::max<T>(T(nElements) - T(1.0),T(1.0));

    T azTanLength   = T(2.0) * T(azR) * tan(T(0.5)*T(dAz));

    T elR = sqrt(SQUARE(focus)+SQUARE(T(0.5)*height));

    // Sector from outer edge to outer edge
    T elSector = T(2.0) * atan2(halfHeight,focus);

    T dEl =  elSector / T(nSubH);

    T elChordLength = T(2.0) * elR * sin(T(0.5)*dEl);

    T elTanLength = T(0.0);
    if (focus != 0.0) {
      elTanLength   = T(2.0) * elR * tan(T(0.5)*dEl);
    } else {
      elTanLength   = height;
    }

    T hh;
    if (arcPlacement == 0) {
      // Outer
      hh = T(0.5)*elTanLength;
    } else {
      // Inner
      hh = T(0.5)*elChordLength;
    }

    SPS_UNREFERENCED_PARAMETER(azTanLength);
    T hw = T(0.5)*(pitch - kerf);

    if (arcPlacement == 0) {
      // Outer radius
      for (size_t iEl = 0 ; iEl < nSubH ; iEl++) {
        // Rotation in elevation
        T elAngle = (T(iEl) - T(0.5)*(T(nSubH)-T(1.0)))*dEl;

        sps::point_t<T> center;
        center[0] = T(0.0);
        center[1] = sin(elAngle)*elR;
        center[2] = -(cos(elAngle)*elR - elR)+azR;
        for (size_t iAz = 0 ; iAz < nElements ; iAz++) {

          T azAngle =  (T(iAz) - T(0.5)*(T(nElements)-T(1.0)))*dAz;
          // Rotation origin in azimuth
          sps::euler_t<T> euler;
          euler.alpha = azAngle;
          euler.beta  = T(0.0);
          euler.gamma = T(0.0);

          sps::point_t<T> center1;
          // Rotate about origin
          sps::basis_rotate<T, sps::EulerIntrinsicYXY>(center,euler,&center1);
          center1[2] = center1[2] - azR;

          elements[iAz][iEl].hh     = hh;
          elements[iAz][iEl].hw     = hw;
          elements[iAz][iEl].center = center1;
          elements[iAz][iEl].euler.alpha = azAngle;
          elements[iAz][iEl].euler.beta  = elAngle;
          elements[iAz][iEl].euler.gamma = T(0.0);
        }
      }
    }
  }

  template void
  FocusedLinearArray(const size_t nElements,
                     const size_t nSubH, // Elevation
                     const size_t nSubW,
                     const float pitch,
                     const float kerf,
                     const float height,
                     const float eFocus,
                     const int arcPlacement,
                     sps::deleted_aligned_multi_array<sps::element_rect_t<float>, 2> &&elements);

  template void
  FocusedConvexArray(const size_t nElements,
                     const size_t nSubH, // Elevation
                     const size_t nSubW,
                     const float pitch,
                     const float kerf,
                     const float height,
                     const float radius,
                     const float eFocus,
                     const int arcPlacement,
                     sps::deleted_aligned_multi_array<sps::element_rect_t<float>, 2> &&elements);

  template void
  MatrixArray(const size_t nRows,
              const size_t nCols,
              const float rowWidth,
              const float rowKerf,
              const float colWidth,
              const float colKerf,
              sps::deleted_aligned_multi_array<sps::element_rect_t<float>, 2> &&elements);

#ifdef FNM_DOUBLE_SUPPORT
  template void
  FocusedLinearArray(const size_t nElements,
                     const size_t nSubH, // Elevation
                     const size_t nSubW,
                     const double pitch,
                     const double kerf,
                     const double height,
                     const double eFocus,
                     const int arcPlacement,
                     sps::deleted_aligned_multi_array<sps::element_rect_t<double>, 2> &&elements);

  template void
  FocusedConvexArray(const size_t nElements,
                     const size_t nSubH, // Elevation
                     const size_t nSubW,
                     const double pitch,
                     const double kerf,
                     const double height,
                     const double radius,
                     const double eFocus,
                     const int arcPlacement,
                     sps::deleted_aligned_multi_array<sps::element_rect_t<double>, 2> &&elements);

  template void
  MatrixArray(const size_t nRows,
              const size_t nCols,
              const double rowWidth,
              const double rowKerf,
              const double colWidth,
              const double colKerf,
              sps::deleted_aligned_multi_array<sps::element_rect_t<double>, 2> &&elements);

#endif

}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
