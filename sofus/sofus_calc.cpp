/**
 * @file   sofus_calc.cpp
 * @author Jens Munk Hansen <jmh@debian9laptop.parknet.dk>
 * @date   Tue Sep 26 19:43:42 2017
 *
 * @brief
 *
 * Copyright 2017 Jens Munk Hansen
 */

// TODO(JEM): Split projection into (compute limits)  + (projection)

#include <sofus/config.h>
#include <sps/debug.h>
#include <sps/profiler.h>

#include <sofus/sofus_calc.hpp>
#include <sofus/sofus_pulses.hpp>

#include <fnm/fnm_data.hpp>

namespace sofus {

template<class T>
void ComputeBoxTimes(const sysparm_t<T>& sysparm,
                     const sps::bbox_t<T>& box,
                     const T* points, const size_t nPoints,
                     const T* delays, const T* apodizations,
                     const size_t nElements,
                     T* tStart,
                     int* iStartSample, size_t* nSamples) {

  const T c  = sysparm.c;
  const T fs = sysparm.fs;

  T delay_min = T(0.0);
  T delay_max = T(0.0);

#if 0
  // Temporarily we ignore apodization for timing to ensure equal tstart
  sps::unique_aligned_array<T> apodizations0 = sps::unique_aligned_array_create<T>(nElements);
  std::fill(apodizations0.get(), apodizations0.get() + nElements, T(1.0));

  // Find minimum and maximum delay (one-way, e.g. transmit)
  std::pair<T, T> delay_minmax =
    sps::minmax_delay(delays, apodizations0.get(), nElements);
#else
  std::pair<T, T> delay_minmax =
    sps::minmax_delay(delays, apodizations, nElements);
#endif
  delay_min = delay_minmax.first;
  delay_max = delay_minmax.second;

  debug_print("delay_min: %g, delay_max: %g\n", delay_min, delay_max);
  T dist_close = T(0.0);
  T dist_far   = T(0.0);

  if (nPoints > 1) {
    // Simple box surrounding the scatters
    sps::bbox_t<T> scatter_box = sps::bbox_t<T>();
    sps::compute_bounding_box3(points, nPoints, &scatter_box);

    // Uses the 8 corner points
    sps::dists_most_distant_and_closest(scatter_box,
                                        box,
                                        &dist_close,
                                        &dist_far);
  } else {

    sps::point_t<T> point;
    point[0] = points[0];
    point[1] = points[1];
    point[2] = points[2];

    sps::point_t<T> point_close = sps::nearest_point_on_bbox(point, box);
    sps::point_t<T> point_far   = sps::farthest_point_on_bbox(point, box);

    // Distance to nearest and farthest point
    dist_close  = sps::dist_point_to_point(point, point_close);
    dist_far    = sps::dist_point_to_point(point, point_far);
  }

  T _tStart = dist_close / c + delay_min;
  T tEnd = dist_far / c + delay_max;

  // Length of spatial response (1 rounding + 1 range)
  size_t _nSamples = (size_t) ceil(T(2.0) + (tEnd - _tStart)*fs);

  *nSamples = _nSamples;
  // Redundant parameter, tStart is enough
  *iStartSample = static_cast<int>(floor(_tStart*fs));
  *tStart = _tStart;
}

template void
SOFUS_EXPORT ComputeBoxTimes(const sysparm_t<float>& sysparm,
                             const sps::bbox_t<float>& box,
                             const float* points, const size_t nPoints,
                             const float* delays, const float* apodizations,
                             const size_t nElements,
                             float* tStart,
                             int* iSampleStart, size_t* nSamples);
}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
