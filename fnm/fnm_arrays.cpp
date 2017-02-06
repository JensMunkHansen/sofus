#include <fnm/fnm_arrays.hpp>
#include <sps/smath.hpp>

#include <cassert>
namespace fnm {
  template <class T>
  void FocusedLinearArray(const size_t nElements,
                          const size_t nSubH, // Elevation
                          const size_t nSubW,
                          const T pitch,
                          const T kerf,
                          const T height,
                          const T efocus,
                          const int arcPlacement,
                          sps::deleted_aligned_multi_array<sps::element_t<T>, 2> &&elements)
  {

    // TODO: Support nSubW != 1
    assert(nSubW==1);
    elements      = sps::deleted_aligned_multi_array_create<sps::element_t<T>,2>(nElements,nSubH);

    T focus = efocus;
    T halfHeight = T(0.5)*height;

    T R = sqrt(SQUARE(focus)+SQUARE(halfHeight));

    T elSector = T(2.0)*atan2(halfHeight,focus);

    T dEl = elSector / T(nSubH);
    T hw  = T(0.5)*(pitch - kerf);

    T chordLength = T(2.0) * R * sin(T(0.5)*dEl);
    T tanLength   = T(2.0) * R * tan(T(0.5)*dEl);

    T hh = T(0.0);

    if (arcPlacement == 0) {
      hh = T(0.5)*tanLength;
    } else {
      hh = T(0.5)*chordLength;
    }

    T elAngle     = T(0.0);
    T lastElAngle = T(0.0);

    if (arcPlacement == 0) {
      // Placement outside arc
      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {

        for (size_t iSubH = 0 ; iSubH < nSubH ; iSubH++) {
          elAngle = (T(iSubH) - T(0.5)*(T(nSubH)-T(1.0)))*dEl;
          sps::point_t<T> center;
          center[0] = (T(iElement) - T(0.5)*(T(nElements)-T(1.0)))*pitch;
          center[1] = sin(elAngle)*R;
          center[2] = -(cos(elAngle)*R-focus);

          for (size_t iSubW = 0 ; iSubW < nSubW ; iSubW++) {
            // TODO: Consider allowing nSubW > 1
            sps::point_t<T> center1 = center;
            center1[0] = center1[0] + T(2.0)*hw/T(nSubW)*(T(iSubW)-T(0.5)*(T(nSubW)-T(1.0)));

            elements[iElement][iSubH].hh     = hh;
            elements[iElement][iSubH].hw     = hw/T(nSubW);
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

          sps::point_t<T> center;
          center[0] = (T(iElement) - T(0.5)*(T(nElements)-T(1.0)))*pitch;
          center[1] = T(0.5)*(sin(elAngle)+sin(lastElAngle))*R;
          center[2] = -(T(0.5)*(cos(elAngle)+cos(lastElAngle))*R-focus);

          for (size_t iSubW = 0 ; iSubW < nSubW ; iSubW++) {
            // TODO: Consider allowing nSubW > 1
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
                          const T pitch,
                          const T kerf,
                          const T height,
                          const T radius,
                          const T efocus,
                          const int arcPlacement,
                          sps::deleted_aligned_multi_array<sps::element_t<T>, 2> &&elements)
  {

    assert(nSubW==1);
    assert(arcPlacement==0);

    elements      = sps::deleted_aligned_multi_array_create<sps::element_t<T>,2>(nElements,nSubH);

    T halfHeight = T(0.5)*height;
    T focus = efocus;
    if (focus < T(0.0))
      focus = T(0.0);

    T azR = radius;

    // Arc-length from center to center
    T azArcLength = pitch * (T(nElements)-T(1.0));

    T azSegment = azArcLength / azR;

    T dAz = azSegment / std::max<T>(T(nElements) - T(1.0),T(1.0));

    T azTanLength   = T(2.0) * T(azR) * tan(T(0.5)*T(dAz));

    //azAngles = (np.r_[0:opt.nElements] - (opt.nElements - 1.0)/2) * dAz

    T elR = sqrt(SQUARE(focus)+SQUARE(T(0.5)*height));

    // Sector from outer edge to outer edge
    T elSector = T(2.0) * atan2(halfHeight,focus);

    T dEl =  elSector / T(nSubH);

    T elChordLength = T(2.0) * elR * sin(T(0.5)*dEl);

    T elTanLength = T(0.0);
    if (focus != 0.0) {
      elTanLength   = T(2.0) * elR * tan(T(0.5)*dEl);
    } else {
      elTanLength   = T(0.5)*height;
    }

    T hh;
    if (arcPlacement == 0) {
      // Outer
      hh = T(0.5)*elTanLength;
    } else {
      // Inner
      //elAngles = (np.r_[0:(opt.nSubH+1)] - opt.nSubH/2.0) * dEl
      hh = T(0.5)*elChordLength;
    }

    T hw = T(0.5)*azTanLength;

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
          sps::basis_rotate<T>(center,euler,center1);
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

  template void FNM_EXPORT FocusedLinearArray(const size_t nElements,
      const size_t nSubH, // Elevation
      const size_t nSubW,
      const float pitch,
      const float kerf,
      const float height,
      const float efocus,
      const int arcPlacement,
      sps::deleted_aligned_multi_array<sps::element_t<float>, 2> &&elements);

  template void FNM_EXPORT FocusedConvexArray(const size_t nElements,
      const size_t nSubH, // Elevation
      const size_t nSubW,
      const float pitch,
      const float kerf,
      const float height,
      const float radius,
      const float efocus,
      const int arcPlacement,
      sps::deleted_aligned_multi_array<sps::element_t<float>, 2> &&elements);


}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
