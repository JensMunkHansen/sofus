/**
 * @file   fnm_delays.cpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Wed Apr  4 22:18:44 2018
 *
 * @brief
 *
 * Copyright 2018 Jens Munk Hansen
 */

#include <cmath>
#include <fnm/fnm_delays.hpp>
#include <sps/smath_types.hpp>
#include <sps/smath.hpp>

namespace fnm {

template <class T>
void InteriorEllipsis(const sps::point_t<T>* pPositions,
                      const size_t& nElements,
                      sps::ellipsis_t<T>* pEllipsis) {

  sps::point_t<T> xdcCenter{T(0.0), T(0.0), T(0.0)};
  for (size_t iElement=0 ; iElement < nElements ; ++iElement) {
    xdcCenter += pPositions[iElement];
  }

  T hh = T(0.0);
  T hw = T(0.0);

  for (size_t iElement=0 ; iElement < nElements ; ++iElement) {
    hh = std::max<T>(hh, abs(pPositions[iElement][1] - xdcCenter[1]));
    hw = std::max<T>(hw, abs(pPositions[iElement][0] - xdcCenter[0]));
  }

  pEllipsis->hh = hh;
  pEllipsis->hw = hw;
  pEllipsis->center = xdcCenter;
}

// Never pass smart pointer by value or reference!!!! Pass the underlying array instead by ref

// TODO: Make explicit specialization using std::copysignf<T> for float
template <class T>
int PytDelays0(const sps::unique_aligned_array<sps::point_t<T> >& positions,
               const size_t& nElements,
               const sps::point_t<T>& centerFocus,
               const sps::point_t<T>& focus0,
               const sps::point_t<T>& focus1,
               const T& c,
               sps::unique_aligned_array<T>& delaysOut, bool allPositive) {
  SPS_UNREFERENCED_PARAMETER(centerFocus);
  const T signFocus = std::copysign<T>(T(1.0), focus0[2]);
  const auto& pos = positions;

  // Pythagorean
  auto delays = sps::unique_aligned_array_create<T>(nElements);
  T maxDelay = - std::numeric_limits<T>::max();
  T minDelay = std::numeric_limits<T>::max();
  if (sps::norm(focus1 - sps::point_t<T>({T(0.0), T(0.0), T(0.0)})) > 1e-6) {
    // Axi-quad focusing
    T right_x = focus1[0] - pos[0][0];
    T right_z = focus1[2] - pos[0][2];
    T l_right = sqrt(SQUARE(right_x) + SQUARE(right_z));

    T left_x = focus1[0] - pos[(nElements-1)][0];
    T left_z = focus1[2] - pos[(nElements-1)][2];
    T l_left = sqrt(SQUARE(left_x) + SQUARE(left_z));
    for (size_t iElement=0 ; iElement < nElements ; iElement++) {

      // Axi-quad - not working when focus is behind aperture
      T ele_x = pos[iElement][0];
      T diff_x = ele_x - focus0[0];
      T ele_z = pos[iElement][2];
      T diff_z = ele_z - focus0[2];
      T delay;
      if (diff_x * right_z - diff_z * right_x < 0) {
        delay =
          (diff_x * right_x + diff_z * right_z) / l_right / c;
      } else if (diff_x * left_z - diff_z * left_x > 0) {
        delay =
          (diff_x * left_x + diff_z * left_z) / l_left / c;
      } else {
        delay =
          - sqrt(SQUARE(diff_x) + SQUARE(diff_z)) / c;
      }
      minDelay = std::min<T>(delay, minDelay);
      maxDelay = std::max<T>(maxDelay, delay);
      delays[iElement] = delay;
    }

    for (size_t iElement=0 ; iElement < nElements ; iElement++) {
      // These delays are by design positive
      delaysOut[iElement] = delays[iElement] - minDelay;
      // Shift to negative
      delaysOut[iElement] +=
        allPositive && signFocus > T(0.0) ? T(0.0) : -(maxDelay-minDelay);
    }
  } else {
    // Quadratic focusing
    for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
      delays[iElement] =
        norm(pos[iElement] - focus0) / c;
      maxDelay = std::max<T>(maxDelay, delays[iElement]);
      minDelay = std::min<T>(minDelay, delays[iElement]);
    }
    for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
      // These delays are negative
      delaysOut[iElement] =
        signFocus*(minDelay - delays[iElement]);
      delaysOut[iElement] +=
        allPositive && signFocus > T(0.0) ? (maxDelay - minDelay) : T(0.0);
    }
  }
  return 0;
}

template <class T>
int PytDelays1(const sps::unique_aligned_array<sps::point_t<T> >& positions,
               const size_t& nElements,
               const sps::point_t<T>& centerFocus,
               const sps::point_t<T>& focus0,
               const sps::point_t<T>& focus1,
               const T& c,
               sps::unique_aligned_array<T>& delaysOut, bool allPositive) {
  const T signFocus = std::copysign<T>(T(1.0), focus0[2]);
  const auto& pos = positions;

  // Pythagorean
  auto delays = sps::unique_aligned_array_create<T>(nElements);
  T maxDelay = - std::numeric_limits<T>::max();
  T minDelay = std::numeric_limits<T>::max();

  if (sps::norm(focus1 - sps::point_t<T>({T(0.0), T(0.0), T(0.0)})) > 1e-6) {
    // Axi-quad focusing - the array is anticipated rectangular and
    // axis-aligned
    sps::point_t<T> xdcCenter = {T(0.0), T(0.0), T(0.0)};
    T hh = T(0.0);
    T hw = T(0.0);

    for (size_t iElement=0 ; iElement < nElements ; ++iElement) {
      hh = std::max<T>(hh, abs(pos[iElement][1] - xdcCenter[1]));
      hw = std::max<T>(hw, abs(pos[iElement][0] - xdcCenter[0]));
    }

    sps::ellipsis_t<T> ellipsis;
    ellipsis.hh = hh;
    ellipsis.hw = hw;
    ellipsis.center = xdcCenter;

    for (size_t iElement=0 ; iElement < nElements ; iElement++) {
      // TODO(JMH): Figure this out in 2D - diagonals are wrong
      sps::point_t<T> rightArcPoint;
      sps::point_t<T> leftArcPoint;

      // TODO(JMH): Consider finding nearest boundary instead of left/right
      tan_point_ellipsis(ellipsis, fabs(pos[iElement][1]), -fabs(pos[iElement][0]), &rightArcPoint);
      tan_point_ellipsis(ellipsis, -fabs(pos[iElement][1]), fabs(pos[iElement][0]), &leftArcPoint);

      sps::point_t<T> r2f = focus1 - rightArcPoint;
      sps::point_t<T> l2f = focus1 - leftArcPoint;

      T df2r = norm(r2f);
      T df2l = norm(l2f);

      auto& ele = pos[iElement];

      sps::point_t<T> e2f = ele - focus0;

      T delay;
      // TODO(JMH): Sign of projection
      //
      sps::point_t<T> c2f = focus0 - ellipsis.center;

      if (dot(sps::cross<T>(e2f, l2f), sps::cross<T>(l2f, c2f)) < T(0.0)) {
        delay = dot(e2f, l2f) / df2l / c;
      } else {
        if (dot(sps::cross<T>(e2f, r2f), sps::cross<T>(c2f, r2f)) > T(0.0)) {
          delay = dot(e2f, r2f) / df2r / c;
        } else {
          delay =
            -norm(e2f) / c;
        }
      }
      minDelay = std::min<T>(delay, minDelay);
      maxDelay = std::max<T>(maxDelay, delay);
      delays[iElement] = delay;
    }
    for (size_t iElement = 0; iElement < nElements; iElement++) {
      // These delays are by design positive
      delaysOut[iElement] = delays[iElement] - minDelay;
      // Shift to negative
      delaysOut[iElement] +=
        allPositive && signFocus > T(0.0) ? T(0.0) : -(maxDelay-minDelay);
    }
  } else {
    // Quadratic focusing - respecting center focus
    T centerDelay = norm(centerFocus - focus0) / c;

    for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
      T delay =
        signFocus * (centerDelay - norm(pos[iElement] - focus0) / c);
      minDelay = std::min<T>(delay, minDelay);
      maxDelay = std::max<T>(maxDelay, delay);
      delays[iElement] = delay;
    }
    for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
      delaysOut[iElement] = delays[iElement];
      delaysOut[iElement] +=
        allPositive ? (maxDelay - minDelay) : T(0.0);
    }
    for (size_t iElement=0 ; iElement < nElements ; iElement++) {
      // These delays are by design positive
      delaysOut[iElement] = delays[iElement] - minDelay;
      // Shift to negative
      delaysOut[iElement] +=
        allPositive && signFocus > T(0.0) ? T(0.0) : -(maxDelay-minDelay);
    }
  }
  return 0;
}

// Pythagorean focused and defocused (works)
// Axi-Quad focused (works)
// Axi-Quad defocused (not working)
template <class T>
int PytDelays2(const sps::unique_aligned_array<sps::point_t<T> >& positions,
               const size_t& nElements,
               const sps::point_t<T>& centerFocus,
               const sps::point_t<T>& focus0,
               const sps::point_t<T>& focus1,
               const T& c,
               sps::unique_aligned_array<T>& delaysOut, bool allPositive) {
  const T signFocus = std::copysign<T>(T(1.0), focus0[2]);
  const auto& pos = positions;

  // Pythagorean
  auto delays = sps::unique_aligned_array_create<T>(nElements);
  T maxDelay = - std::numeric_limits<T>::max();
  T minDelay = std::numeric_limits<T>::max();

  T centerDelay = norm(centerFocus - focus0) / c;

  if (sps::norm(focus1 - sps::point_t<T>({T(0.0), T(0.0), T(0.0)})) > 1e-6) {
    // Find interior ellipsis
    // TODO(JMH): Support subapertures
    sps::ellipsis_t<T> ellipsis;
    InteriorEllipsis(positions.get(), nElements, &ellipsis);

    // Anticipate center of ellipsis equals center focus
    sps::point_t<T> c2f = focus0 - ellipsis.center;

    for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
      sps::point_t<T> arcPoints[2];

      auto& ele = pos[iElement];

#if 0
      tan_point_ellipsis(ellipsis, ele[1], -ele[0], &arcPoints[0]);
      tan_point_ellipsis(ellipsis, -ele[1], ele[0], &arcPoints[1]);
#else
      sps::element_rect_t<T> rect;
      rect.hh = ellipsis.hh;
      rect.hw = ellipsis.hw;
      rect.center = ellipsis.center;

      intcp_line_rect(rect, ele[1], -ele[0], &arcPoints[0]);
      intcp_line_rect(rect, -ele[1], ele[0], &arcPoints[1]);
#endif
      sps::point_t<T> arc2focus[2];

      arc2focus[0] = focus1 - arcPoints[0];
      arc2focus[1] = focus1 - arcPoints[1];

      T distFocus2arc[2];

      distFocus2arc[0] = norm(arc2focus[0]);
      distFocus2arc[1] = norm(arc2focus[1]);

      sps::point_t<T> e2f = ele - focus0;

      T delay;

      auto arcDistances = {norm(ele-arcPoints[0]), norm(ele-arcPoints[1])};
      auto arcDistMin = std::min_element(arcDistances.begin(), arcDistances.end());
      size_t iNear = std::distance(arcDistances.begin(), arcDistMin);

      if (dot(sps::cross<T>(e2f, arc2focus[iNear]), sps::cross<T>(arc2focus[iNear], c2f)) < T(0.0)) {
        delay = signFocus*(centerDelay + dot(e2f, arc2focus[iNear]) / distFocus2arc[iNear] / c);
      } else {
        delay = signFocus*(centerDelay - norm(e2f) / c);
      }
      minDelay = std::min<T>(delay, minDelay);
      maxDelay = std::max<T>(maxDelay, delay);
      delays[iElement] = delay;
    }
    // How, if focus is behind
    // delaysOut[iElement] +=
    //    allPositive && signFocus > T(0.0) ? T(0.0) : -(maxDelay-minDelay);
  } else {
    // Quadratic focusing - respecting center focus

    for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
      T delay =
        signFocus * (centerDelay - norm(pos[iElement] - focus0) / c);
      minDelay = std::min<T>(delay, minDelay);
      maxDelay = std::max<T>(maxDelay, delay);
      delays[iElement] = delay;
    }
  }
  for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
    delaysOut[iElement] = delays[iElement];
    delaysOut[iElement] -=
        allPositive ? minDelay : T(0.0);
  }
  return 0;
}

// Not used anywhere!!!
template <class T>
int PythagoreanFocusDelays(
  const sps::point_t<T>* pElementPositions, const size_t& nElements,
  const sps::point_t<T>& focus, const sps::point_t<T>& centerFocus,
  const size_t& nDelays, T* pDelays) {
  T c = 1540.0;

  auto delays = sps::unique_aligned_array_create<T>(nElements);

  T maxDelay = - std::numeric_limits<T>::max();
  T minDelay =  std::numeric_limits<T>::max();

  for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
    delays[iElement] = norm(pElementPositions[iElement] - focus) / c;
    maxDelay = std::max<T>(maxDelay, delays[iElement]);
    minDelay = std::min<T>(minDelay, delays[iElement]);
  }

  // Using center_focus for computing minDelay (TEST - should we do this)
  minDelay = norm(focus - centerFocus) / c;

  for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
    pDelays[iElement] = minDelay - delays[iElement];
#if 0
    // TODO(JEM): Verify if it is okay, if delays and phases are opposite
    for (size_t jElement = 0 ; jElement < nSubElements ; jElement++) {
      m_data->m_phases[iElement*nSubElements+jElement] =
        T(M_2PI) * (delays[iElement] - maxDelay) * m_data->m_f0;
    }
#endif
  }
  return 0;
}

template void
InteriorEllipsis(const sps::point_t<float>* pPositions,
                 const size_t& nElements,
                 sps::ellipsis_t<float>* pEllipsis);

template int
PytDelays0(const sps::unique_aligned_array<sps::point_t<float> >& positions,
           const size_t& nElements,
           const sps::point_t<float>& centerFocus,
           const sps::point_t<float>& focus0,
           const sps::point_t<float>& focus1,
           const float& c,
           sps::unique_aligned_array<float>& delays, bool allPositive);

template int
PytDelays1(const sps::unique_aligned_array<sps::point_t<float> >& positions,
           const size_t& nElements,
           const sps::point_t<float>& centerFocus,
           const sps::point_t<float>& focus0,
           const sps::point_t<float>& focus1,
           const float& c,
           sps::unique_aligned_array<float>& delays, bool allPositive);
template int
PytDelays2(const sps::unique_aligned_array<sps::point_t<float> >& positions,
           const size_t& nElements,
           const sps::point_t<float>& centerFocus,
           const sps::point_t<float>& focus0,
           const sps::point_t<float>& focus1,
           const float& c,
           sps::unique_aligned_array<float>& delays, bool allPositive);
}  // namespace fnm
