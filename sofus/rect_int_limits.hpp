/**
 * @file   rect_int_limits.hpp
 * @author  <jens.munk.hansen@gmail.com>
 * @date   Mon Aug 10 16:05:41 2015
 *
 * @brief
 *
 *
 */

#pragma once

// Cache line is 64 bytes, 16 floats

// TODO: Use T* and restrict in prototypes and __builtin_assume_aligned inside function bodies
// TODO: Replace _mm256_dp_pd with _mm256_hsum_pd

#include <sofus/config.h>
#include <sofus/sofus_types.hpp>

#include <sps/cenv.h>
#include <sps/smath.hpp>
#include <sps/extintrin.h>
#include <sps/debug.h>

#include <cstring>
#include <cassert>

namespace sofus {

  /**
   *
   *
   * @param sysparm
   * @param element
   * @param point
   * @param delay
   * @param limits
   *
   * @return
   */
  template <class T>
  STATIC_INLINE_BEGIN int
  calcProjectionAndLimits(const sysparm_t<T>& __restrict sysparm,
                          const sps::element_rect_t<T>& __restrict element,
                          const sps::point_t<T>& __restrict point,
                          const T& __restrict delay,
                          proj_limit_dist_t<T>* __restrict limits) STATIC_INLINE_END;

  /**
   * Compute distance to plane, projection, distance to vertices and sample interval
   *
   * @param sysparm
   * @param element
   * @param point
   * @param delay
   * @param u            Projection, u-coordinate [m]
   * @param v            Projection, v-coordinate [m]
   * @param dist2plane   Distance to plane        [m]
   * @param arrivalTimes Arrivaltimes             [s]
   * @param vdists       Distance to vertices     [m]
   *
   * @return
   */
  template <class T>
  STATIC_INLINE_BEGIN int
  calcProjectionAndBoundaries(const sysparm_t<T>& sysparm,
                              const sps::element_rect_t<T>& element,
                              const sps::point_t<T>& point,
                              const T& delay,
                              T* u,
                              T* v,
                              T* dist2plane,
                              T (*arrivalTimes)[3],
                              T (*vdists)[4]) STATIC_INLINE_END;

  /**
   * Compute distance to plane, projection, distance to vertices and sample interval
   *
   * @param sysparm
   * @param element
   * @param point
   * @param delay
   * @param v_vuvu       Projection           [m]
   * @param v_dist2plane Distance to plane    [m]
   * @param v_hwhw       Half width/height    [m]
   * @param v_vdists     Distance to vertices [m]
   * @param fSampleStart Fractional sample start
   * @param fSampleStop  Fractional sample stop
   *
   * @return
   */
  template <typename T>
  STATIC_INLINE_BEGIN bool
  calcProjectionAndIntegrationLimitsSIMD(const sysparm_t<T>& sysparm,
                                         const sps::element_rect_t<T>& element,
                                         const sps::point_t<T>& point,
                                         const T delay,
                                         __m128 *v_vuvu,
                                         __m128 *v_dist2plane,
                                         __m128 *v_hwhw,
                                         __m128 *v_vdists,
                                         T* fSampleStart,
                                         T* fSampleStop) STATIC_INLINE_END;

  template <typename T>
  STATIC_INLINE_BEGIN bool
  calcProjectionAndIntegrationLimitsSIMD2(const sysparm_t<T>& sysparm,
                                          const sps::element_rect_t<T>& element,
                                          const sps::point_t<T>& point,
                                          const T delay,
                                          __m256d *v_vuvu,
                                          __m256d *v_dist2plane,
                                          __m256d *v_hwhw,
                                          __m256d *v_vdists,
                                          T* fSampleStart,
                                          T* fSampleStop) STATIC_INLINE_END;


  template <>
  inline bool
  calcProjectionAndIntegrationLimitsSIMD(const sysparm_t<double>& sysparm,
                                         const sps::element_rect_t<double>& element,
                                         const sps::point_t<double>& point,
                                         const double delay,
                                         __m128 *v_vuvu,
                                         __m128 *v_dist2plane,
                                         __m128 *v_hwhw,
                                         __m128 *v_vdists,
                                         double* fSampleStart,
                                         double* fSampleStop)
  {
    SPS_UNREFERENCED_PARAMETERS(sysparm,
                                element,
                                point,
                                delay,
                                v_vuvu,
                                v_dist2plane,
                                v_hwhw,
                                v_vdists,
                                fSampleStart,
                                fSampleStop);

    // Note: TODO: Consider using aligned T* arguments and use __m256 to hold 4 vertex distances
    assert(false && "Makes no sense for double precision: Avoid __m128 types on interface");
    return true;
  }

  template <>
  inline bool
  calcProjectionAndIntegrationLimitsSIMD(const sysparm_t<float>& sysparm,
                                         const sps::element_rect_t<float>& element,
                                         const sps::point_t<float>& point,
                                         const float delay,
                                         __m128 *v_vuvu,                       /* Response    */
                                         __m128 *v_dist2plane,                 /* Response    */
                                         __m128 *v_hwhw,                       /* Trivial     */
                                         __m128 *v_vdists,                     /* Response    */
                                         float* fSampleStart,                  /* Integration */
                                         float* fSampleStop)                   /* Integration */
  {

    const float _invc       = 1.0f / sysparm.c;             // [s/m]
    const float hh          = element.hh;                   // [m]
    const float hw          = element.hw;                   // [m]

    // Distance to plane, and projections u,v
    ALIGN16_BEGIN float distuv[4] ALIGN16_END;

    __m128 v_r2p = _mm_sub_ps(
                     _mm_load_ps((float*)&point[0]),
                     _mm_load_ps((float*)&element.center[0]));

    // Distance to plane (TODO: Consider alias to __m128*
    __m128 v_dXvu   = _mm_fabs_ps(_mm_dp_ps(_mm_load_ps((float*)&element.normal[0]),v_r2p,0x71));

    // Projection (u,v)
    v_dXvu   = _mm_add_ps(v_dXvu,_mm_dp_ps(_mm_load_ps((float*)&element.uvector[0]),v_r2p,0x78));
    v_dXvu   = _mm_add_ps(v_dXvu,_mm_dp_ps(_mm_load_ps((float*)&element.vvector[0]),v_r2p,0x74));

    // Store results
    _mm_store_ps(&distuv[0],v_dXvu);

    // TODO: Consider introducing w = distuv[0]
    float u = distuv[3];
    float v = distuv[2];

    // Time to plane
    float planetime = (distuv[0] * _invc) + delay;

    // Only works for float, but it is faster
    const __m128 v_fs          = _mm_broadcast_ss((float*)&sysparm.fs);
    const __m128 v_invc        = _mm_broadcast_ss((float*)&_invc);
    const __m128 v_delay       = _mm_broadcast_ss((float*)&delay);

    const __m128 v_pmmp        = _mm_set_ps(1.0f,-1.0f,-1.0f,1.0f);

    // Vertices of the elements
    __m128 v_vertices_x;
    __m128 v_vertices_y;
    __m128 v_vertices_z;

    // Read cached vertices
    v_vertices_x = _mm_load_ps(element.vertices[0]);
    v_vertices_y = _mm_load_ps(element.vertices[1]);
    v_vertices_z = _mm_load_ps(element.vertices[2]);

    // Point
#if HAVE_IMMINTRIN_H
    __m128 v_px = _mm_broadcast_ss((float*)&point[0]);
    __m128 v_py = _mm_broadcast_ss((float*)&point[1]);
    __m128 v_pz = _mm_broadcast_ss((float*)&point[2]);
#else
    __m128 v_px = _mm_set1_ps((float)point[0]);
    __m128 v_py = _mm_set1_ps((float)point[1]);
    __m128 v_pz = _mm_set1_ps((float)point[2]);
#endif

    // Point to vertices
    __m128 v_pv_x = _mm_sub_ps(v_vertices_x,v_px);
    __m128 v_pv_y = _mm_sub_ps(v_vertices_y,v_py);
    __m128 v_pv_z = _mm_sub_ps(v_vertices_z,v_pz);

    // Vertex dists
    *v_vdists =
      _mm_sqrt_ps(
        _mm_add_ps(
          _mm_add_ps(
            _mm_square_ps(v_pv_x),
            _mm_square_ps(v_pv_y)),
          _mm_square_ps(v_pv_z)));

    // Vertex times
    __m128 v_vtimes     = _mm_add_ps(_mm_mul_ps(*v_vdists,v_invc),v_delay);

#if 0
    debug_print("v_dXvu[0]: %f, v_dXvu[1]: %f, v_dXvu[2]: %f, v_dXvu[3]: %f\n",
                v4f(v_dXvu).f32[0],v4f(v_dXvu).f32[1],v4f(v_dXvu).f32[2],v4f(v_dXvu).f32[3]);
#endif
    *v_vuvu       = _mm_movehl_ps(v_dXvu,v_dXvu);

#if 0
    // f2_test_lines 0.0, -0.001250, 0.0, -0.001250
    debug_print("v_vuvu[0]: %f, v_vuvu[1]: %f, v_vuvu[2]: %f, v_vuvu[3]: %f\n",
                v4f(*v_vuvu).f32[0],v4f(*v_vuvu).f32[1],v4f(*v_vuvu).f32[2],v4f(*v_vuvu).f32[3]);
#endif
#if HAVE_ZMMINTRIN_H
    *v_dist2plane = _mm_broadcastss_ps(v_dXvu);
#elif HAVE_IMMINTRIN_H
    *v_dist2plane = _mm_broadcast_ss((float*)&distuv[0]);
#else
    *v_dist2plane = _mm_set1_ps(distuv[0]);
#endif

    *v_hwhw       = _mm_set_ps(hw,hh,hw,hh);

    // Edge times (from scatter point to the edges)
#if USE_PROJECTIONS
    __m128 v_etimes =
      _mm_add_ps(
        _mm_mul_ps(
          _mm_sqrt_ps(
            _mm_add_ps(
              _mm_square_ps(*v_dist2plane),
              _mm_square_ps(
                _mm_sub_ps(
                  _mm_mul_ps(
                    v_pmmp,
                    *v_hwhw),
                  *v_vuvu)))),
          v_invc),
        v_delay);
#else
    // Compute expensive cross-products
    __m128 v_etimes = _mm_setzero_ps();
#pragma error "Fix me"
#endif

    __m128 v_f_etimes = _mm_mul_ps(v_etimes,v_fs);
    __m128 v_f_vtimes = _mm_mul_ps(v_vtimes,v_fs);

    // Consider arctan + min(edge_time,v_time_1,v_time_2)


    ALIGN16_BEGIN float f_etimes[4] ALIGN16_END;
    ALIGN16_BEGIN float f_vtimes[4] ALIGN16_END;

    _mm_store_ps((float*)&f_etimes[0],v_f_etimes);
    _mm_store_ps((float*)&f_vtimes[0],v_f_vtimes);

    // Time to reach nearest edge
    float s_etimes_min = _mm_min1_ps(v_f_etimes);

    // Time to reach nearest vertex
    float s_vtimes_min = _mm_min1_ps(v_f_vtimes);

    // Time to reach farthest vertex (last sample - any case)
    float s_vtimes_max = _mm_max1_ps(v_f_vtimes);

    // Sampled time
    float s_planetime = planetime * sysparm.fs;

    // Start sample when reaching plane before any edge (inside)
    float sSamplePlane = 0.0f;

    // Start sample when arc length are computed
    float sSampleLines = 0.0f;

    // Stop sample (excluded)
    float sStopSample = s_vtimes_max;

    // First loop is excluded unless we are inside
    sSamplePlane = sStopSample;

    // Find boundary samples
    if ((fabs(u) > hw) && (fabs(v) > hh)) {
      // Outside and in the vertices (a vertex is hit first)
      sSampleLines = s_vtimes_min;
    } else if (u > hw) {
      // Outside
      sSampleLines = f_etimes[3];
    } else if (u < -hw) {
      // Outside
      sSampleLines = f_etimes[1];
    } else {
      if (v > hh) {
        // Outside
        sSampleLines = f_etimes[0];
      } else if (v < -hh) {
        // Outside
        sSampleLines = f_etimes[2];
      } else {
        // Inside
        sSampleLines = s_etimes_min;
        sSamplePlane = s_planetime;
      }
    }

    *fSampleStop = sStopSample;
    *fSampleStart = std::min<float>(sSampleLines,sSamplePlane);

    return true;
  }

  template <>
  inline bool
  calcProjectionAndIntegrationLimitsSIMD2(const sysparm_t<double>& sysparm,
                                          const sps::element_rect_t<double>& element,
                                          const sps::point_t<double>& point,
                                          const double delay,
                                          __m256d *v_vuvu,                       /* Response    */
                                          __m256d *v_dist2plane,                 /* Response    */
                                          __m256d *v_hwhw,                       /* Trivial     */
                                          __m256d *v_vdists,                     /* Response    */
                                          double* fSampleStart,                  /* Integration */
                                          double* fSampleStop)                   /* Integration */
  {

    const double _invc       = 1.0 / sysparm.c;              // [s/m]
    const double hh          = element.hh;                   // [m]
    const double hw          = element.hw;                   // [m]

    // Distance to plane, and projections u,v
    ALIGN32_BEGIN double distuv[4] ALIGN32_END;

    assert(((uintptr_t)&point[0] & 0x1F) == 0 && "Data must be aligned");
    assert(((uintptr_t)&element.center[0] & 0x1F) == 0 && "Data must be aligned");

    __m256d v_r2p = _mm256_sub_pd(
                      _mm256_load_pd((double*)&point[0]),
                      _mm256_load_pd((double*)&element.center[0]));

    assert(((uintptr_t)&element.normal[0] & 0x1F) == 0 && "Data must be aligned");

    __m256d v_dXvu   = _mm256_fabs_pd(_mm256_dp_pd(_mm256_load_pd((double*)&element.normal[0]),v_r2p,0x71));

    // Projection (u,v)
    assert(((uintptr_t)&element.uvector[0] & 0x1F) == 0 && "Data must be aligned");
    assert(((uintptr_t)&element.vvector[0] & 0x1F) == 0 && "Data must be aligned");
    v_dXvu   = _mm256_add_pd(v_dXvu,_mm256_dp_pd(_mm256_load_pd((double*)&element.uvector[0]),v_r2p,0x78));
    v_dXvu   = _mm256_add_pd(v_dXvu,_mm256_dp_pd(_mm256_load_pd((double*)&element.vvector[0]),v_r2p,0x74));

    // Store results
    _mm256_store_pd(&distuv[0],v_dXvu);

    // TODO: Consider introducing w = distuv[0]
    double u = distuv[3];
    double v = distuv[2];

    // Time to plane
    double planetime = (distuv[0] * _invc) + delay;

    // Only works for double, but it is faster
    const __m256d v_fs          = _mm256_broadcast_sd((double*)&sysparm.fs);
    const __m256d v_invc        = _mm256_broadcast_sd((double*)&_invc);
    const __m256d v_delay       = _mm256_broadcast_sd((double*)&delay);

    const __m256d v_pmmp        = _mm256_set_pd(1.0,-1.0,-1.0,1.0);

    // Vertices of the elements
    __m256d v_vertices_x;
    __m256d v_vertices_y;
    __m256d v_vertices_z;

    // Read cached vertices
    v_vertices_x = _mm256_load_pd(element.vertices[0]);
    v_vertices_y = _mm256_load_pd(element.vertices[1]);
    v_vertices_z = _mm256_load_pd(element.vertices[2]);

    // Point
#if HAVE_IMMINTRIN_H
    __m256d v_px = _mm256_broadcast_sd((double*)&point[0]);
    __m256d v_py = _mm256_broadcast_sd((double*)&point[1]);
    __m256d v_pz = _mm256_broadcast_sd((double*)&point[2]);
#else
    __m256d v_px = _mm256_set1_pd((double)point[0]);
    __m256d v_py = _mm256_set1_pd((double)point[1]);
    __m256d v_pz = _mm256_set1_pd((double)point[2]);
#endif

    // Point to vertices
    __m256d v_pv_x = _mm256_sub_pd(v_vertices_x,v_px);
    __m256d v_pv_y = _mm256_sub_pd(v_vertices_y,v_py);
    __m256d v_pv_z = _mm256_sub_pd(v_vertices_z,v_pz);

    // Vertex dists
    *v_vdists =
      _mm256_sqrt_pd(
        _mm256_add_pd(
          _mm256_add_pd(
            _mm256_square_pd(v_pv_x),
            _mm256_square_pd(v_pv_y)),
          _mm256_square_pd(v_pv_z)));

    // Vertex times
    __m256d v_vtimes     = _mm256_add_pd(_mm256_mul_pd(*v_vdists,v_invc),v_delay);

    debug_print("v_dXvu[0]: %f, v_dXvu[1]: %f, v_dXvu[2]: %f, v_dXvu[3]: %f\n",
                v4d(v_dXvu).f64[0],v4d(v_dXvu).f64[1],v4d(v_dXvu).f64[2],v4d(v_dXvu).f64[3]);

    *v_vuvu       = _mm256_movehl_pd(v_dXvu,v_dXvu);

    // f2_double 0.00 -0.001250, 0.015 0.00
    debug_print("v_vuvu[0]: %f, v_vuvu[1]: %f, v_vuvu[2]: %f, v_vuvu[3]: %f\n",
                v4d(*v_vuvu).f64[0],v4d(*v_vuvu).f64[1],v4d(*v_vuvu).f64[2],v4d(*v_vuvu).f64[3]);

#if 0 //HAVE_ZMMINTRIN_H
    *v_dist2plane = _mm256_broadcastsd_pd(v_dXvu);
#elif 1 // HAVE_IMMINTRIN_H
    *v_dist2plane = _mm256_broadcast_sd((double*)&distuv[0]);
#else
    *v_dist2plane = _mm256_set1_pd(distuv[0]);
#endif

    *v_hwhw       = _mm256_set_pd(hw,hh,hw,hh);

    // Edge times (from scatter point to the edges)
#if USE_PROJECTIONS
    __m256d v_etimes =
      _mm256_add_pd(
        _mm256_mul_pd(
          _mm256_sqrt_pd(
            _mm256_add_pd(
              _mm256_square_pd(*v_dist2plane),
              _mm256_square_pd(
                _mm256_sub_pd(
                  _mm256_mul_pd(
                    v_pmmp,
                    *v_hwhw),
                  *v_vuvu)))),
          v_invc),
        v_delay);
#else
    // Compute expensive cross-products
    __m256d v_etimes = _mm256_setzero_pd();
#pragma error "Fix me"
#endif

    __m256d v_f_etimes = _mm256_mul_pd(v_etimes,v_fs);
    __m256d v_f_vtimes = _mm256_mul_pd(v_vtimes,v_fs);

    // Consider arctan + min(edge_time,v_time_1,v_time_2)


    ALIGN32_BEGIN double f_etimes[4] ALIGN32_END;
    ALIGN32_BEGIN double f_vtimes[4] ALIGN32_END;

    _mm256_store_pd((double*)&f_etimes[0],v_f_etimes);
    _mm256_store_pd((double*)&f_vtimes[0],v_f_vtimes);

    // Time to reach nearest edge
    double s_etimes_min = _mm256_min1_pd(v_f_etimes);

    // Time to reach nearest vertex
    double s_vtimes_min = _mm256_min1_pd(v_f_vtimes);

    // Time to reach farthest vertex (last sample - any case)
    double s_vtimes_max = _mm256_max1_pd(v_f_vtimes);

    // Sampled time
    double s_planetime = planetime * sysparm.fs;

    // Start sample when reaching plane before any edge (inside)
    double sSamplePlane = 0.0;

    // Start sample when arc length are computed
    double sSampleLines = 0.0;

    // Stop sample (excluded)
    double sStopSample = s_vtimes_max;

    // First loop is excluded unless we are inside
    sSamplePlane = sStopSample;

    // Find boundary samples
    if ((fabs(u) > hw) && (fabs(v) > hh)) {
      // Outside and in the vertices (a vertex is hit first)
      sSampleLines = s_vtimes_min;
    } else if (u > hw) {
      // Outside
      sSampleLines = f_etimes[3];
    } else if (u < -hw) {
      // Outside
      sSampleLines = f_etimes[1];
    } else {
      if (v > hh) {
        // Outside
        sSampleLines = f_etimes[0];
      } else if (v < -hh) {
        // Outside
        sSampleLines = f_etimes[2];
      } else {
        // Inside
        sSampleLines = s_etimes_min;
        sSamplePlane = s_planetime;
      }
    }

    *fSampleStop = sStopSample;
    *fSampleStart = std::min<double>(sSampleLines,sSamplePlane);

    return true;
  }

  template <typename T>
  STATIC_INLINE_BEGIN
  int calcArrivalTimesSIMD(const sysparm_t<T>& sysparm,
                           const sps::element_rect_t<T>& element,
                           const sps::point_t<T>& point,
                           const T& delay,
                           __m128 *v_hwhw,
                           __m128 *v_vuvu,
                           __m128 *v_dist2plane,
                           __m128 *v_vdists,
                           T (*arrivalTimes)[9]);

  template <>
  inline int
  calcArrivalTimesSIMD(const sysparm_t<double>& sysparm,
                       const sps::element_rect_t<double>& element,
                       const sps::point_t<double>& point,
                       const double& delay,
                       __m128 *v_hwhw,
                       __m128 *v_vuvu,
                       __m128 *v_dist2plane,
                       __m128 *v_vdists,
                       double (*arrivalTimes)[9])
  {
    SPS_UNREFERENCED_PARAMETERS(sysparm, element, point, delay, v_hwhw, v_vuvu, v_dist2plane, v_vdists, arrivalTimes);
    assert(false && "Makes no sense for double precision: Avoid __m128 types on interface");
    return 0;
  }


  template <>
  inline int
  calcArrivalTimesSIMD(const sysparm_t<float>& sysparm,
                       const sps::element_rect_t<float>& element,
                       const sps::point_t<float>& point,
                       const float& delay,
                       __m128 *v_hwhw,
                       __m128 *v_vuvu,
                       __m128 *v_dist2plane,
                       __m128 *v_vdists,
                       float (*arrivalTimes)[9])
  {

    const float _invc       = 1.0f / sysparm.c;//Aperture<float>::_c;   // [s/m]
    const float hh          = element.hh;                   // [m]
    const float hw          = element.hw;                   // [m]

    const __m128 v_invc     = _mm_broadcast_ss((float*)&_invc);
    const __m128 v_delay    = _mm_broadcast_ss((float*)&delay);

    const __m128 v_pmmp     = _mm_set_ps(1.0f,-1.0f,-1.0f,1.0f);

    int nDiscontinuities = 0;

    __m128 v_r2p = _mm_sub_ps(
                     _mm_load_ps((float*)&point[0]),
                     _mm_load_ps((float*)&element.center[0]));

    // Distance to plane
    v4f v_dXvu;
    v_dXvu.v   = _mm_fabs_ps(_mm_dp_ps(_mm_load_ps(element.normal),v_r2p,0x71));

    // Projection (u,v)
    v_dXvu.v   = _mm_add_ps(v_dXvu.v,_mm_dp_ps(_mm_load_ps(element.uvector),v_r2p,0x78));
    v_dXvu.v   = _mm_add_ps(v_dXvu.v,_mm_dp_ps(_mm_load_ps(element.vvector),v_r2p,0x74));

    // Vertices of the elements
    __m128 v_vertices_x;
    __m128 v_vertices_y;
    __m128 v_vertices_z;

    // Read cached vertices
    v_vertices_x = _mm_load_ps(element.vertices[0]);
    v_vertices_y = _mm_load_ps(element.vertices[1]);
    v_vertices_z = _mm_load_ps(element.vertices[2]);

    // Point
    __m128 v_px = _mm_broadcast_ss((float*)&point[0]);
    __m128 v_py = _mm_broadcast_ss((float*)&point[1]);
    __m128 v_pz = _mm_broadcast_ss((float*)&point[2]);

    // Point to vertices
    __m128 v_pv_x = _mm_sub_ps(v_vertices_x,v_px);
    __m128 v_pv_y = _mm_sub_ps(v_vertices_y,v_py);
    __m128 v_pv_z = _mm_sub_ps(v_vertices_z,v_pz);

    // Vertex dists
    *v_vdists =
      _mm_sqrt_ps(
        _mm_add_ps(
          _mm_add_ps(
            _mm_square_ps(v_pv_x),
            _mm_square_ps(v_pv_y)),
          _mm_square_ps(v_pv_z)));

    // Vertex times
    __m128 v_vtimes     = _mm_add_ps(_mm_mul_ps(*v_vdists,v_invc),v_delay);

    *v_vuvu             = _mm_movehl_ps(v_dXvu.v,v_dXvu.v);

    debug_print("v_vuvu[0]: %f, v_vuvu[1]: %f, v_vuvu[2]: %f, v_vuvu[3]: %f\n",
                v4f(*v_vuvu).f32[0],v4f(*v_vuvu).f32[1],v4f(*v_vuvu).f32[2],v4f(*v_vuvu).f32[3]);

#if HAVE_ZMMINTRIN_H
    *v_dist2plane = _mm_broadcastss_ps(v_dXvu.v);
#elif HAVE_IMMINTRIN_H
    *v_dist2plane = _mm_broadcast_ss((float*)&v_dXvu.f32[0]);
#else
    *v_dist2plane = _mm_set1_ps(v_dXvu.f32[0]);
#endif

    float f_planetime = v_dXvu.f32[0] * _invc + delay;

    float u = v_dXvu.f32[3];
    float v = v_dXvu.f32[2];

    *v_hwhw       = _mm_set_ps(hw,hh,hw,hh);

    // Edge times (from scatter point to the edges)
#if USE_PROJECTIONS
    __m128 v_etimes =
      _mm_add_ps(
        _mm_mul_ps(
          _mm_sqrt_ps(
            _mm_add_ps(
              _mm_square_ps(*v_dist2plane),
              _mm_square_ps(
                _mm_sub_ps(
                  _mm_mul_ps(
                    v_pmmp,
                    *v_hwhw),
                  *v_vuvu)))),
          v_invc),
        v_delay);
#else
    // Compute expensive cross-products
    __m128 v_etimes = _mm_setzero_ps();
#pragma error "Fix me"
#endif

    v4f v_f_etimes, v_f_vtimes;
    v_f_etimes.v = v_etimes;//_mm_mul_ps(v_etimes,v_fs);
    v_f_vtimes.v = v_vtimes;//_mm_mul_ps(v_vtimes,v_fs);

    std::fill_n(*arrivalTimes, 9, float(INT_MAX));

    // We always need the vertices
    memcpy(*arrivalTimes,&v_f_vtimes.f32[0],4*sizeof(float));

    nDiscontinuities = 4;

    // Experiment with this (also accumulate multiple contributions at discontinuities)

    // Find boundary samples
    if (fabs(u) > element.hw) {
      if (fabs(v) < element.hh) {
        (*arrivalTimes)[4] = v_f_etimes.f32[1];
        (*arrivalTimes)[5] = v_f_etimes.f32[3];
        nDiscontinuities = 6;
      }
    } else if (fabs(v) > element.hh) {
      (*arrivalTimes)[4] = v_f_etimes.f32[0];
      (*arrivalTimes)[5] = v_f_etimes.f32[2];
      nDiscontinuities = 6;
    } else {
      // Inside
      memcpy(*arrivalTimes + 4,&v_f_etimes.f32[0],4*sizeof(float));
      (*arrivalTimes)[8] = f_planetime;
      nDiscontinuities = 9;
    }

    std::sort(*arrivalTimes, *arrivalTimes + 9, std::less<float>());

    return nDiscontinuities;
  }

  template <class T>
  inline int
  calcProjectionAndBoundaries(const sysparm_t<T>& sysparm,
                              const sps::element_rect_t<T>& element,
                              const sps::point_t<T>& point,
                              const T& delay,
                              T* u,
                              T* v,
                              T* dist2plane,
                              T (*arrivalTimes)[3],
                              T (*vdists)[4])
  {

    const T _invc      = T(1.0) / sysparm.c;         // [s/m]

    // Avoid aliasing
    const T hh = element.hh;
    const T hw = element.hw;

    sps::point_t<T> vertices[4];

    // Compute vertices (if needed)
    for (size_t i_xyz = 0 ; i_xyz < 3 ; i_xyz++) {
      vertices[0][i_xyz] = element.vertices[i_xyz][0];
      vertices[1][i_xyz] = element.vertices[i_xyz][1];
      vertices[2][i_xyz] = element.vertices[i_xyz][2];
      vertices[3][i_xyz] = element.vertices[i_xyz][3];
    }

    // Vertex distances
    (*vdists)[0] = norm(point - vertices[0]);
    (*vdists)[1] = norm(point - vertices[1]);
    (*vdists)[2] = norm(point - vertices[2]);
    (*vdists)[3] = norm(point - vertices[3]);

    // Vertex times
    T vtimes[4];
    vtimes[0] = (*vdists)[0] * _invc + delay;
    vtimes[1] = (*vdists)[1] * _invc + delay;
    vtimes[2] = (*vdists)[2] * _invc + delay;
    vtimes[3] = (*vdists)[3] * _invc + delay;

    sps::point_t<T> r2p = point - element.center;

    sps::point_t<T> hh_dir, hw_dir, normal;
    sps::basis_vectors<T, sps::EulerIntrinsicYXY>(hw_dir, element.euler, 0);
    sps::basis_vectors<T, sps::EulerIntrinsicYXY>(hh_dir, element.euler, 1);
    sps::basis_vectors<T, sps::EulerIntrinsicYXY>(normal, element.euler, 2);

    // Distance to plane
    *dist2plane = fabs(dot(normal,r2p));

    // Time to plane
    T timePlane = (*dist2plane * _invc) + delay;

    // Projection onto plane
    *u = dot(hw_dir,r2p);
    *v = dot(hh_dir,r2p);

    // Edge times (from scatter point to the edges)
    T etimes[4];

    // Distances computed using Pythagorean formula
    etimes[0] = sqrt(SQUARE(*dist2plane) + SQUARE( hh - *v)) * _invc + delay;
    etimes[1] = sqrt(SQUARE(*dist2plane) + SQUARE(-hw - *u)) * _invc + delay;
    etimes[2] = sqrt(SQUARE(*dist2plane) + SQUARE(-hh - *v)) * _invc + delay;
    etimes[3] = sqrt(SQUARE(*dist2plane) + SQUARE( hw - *u)) * _invc + delay;

    std::pair<T*,T*> dist_range;

    dist_range = std::minmax_element(etimes,
                                     etimes + 4,
                                     std::less<T>());
    // Time to reach nearest edge
    T etimes_min = *dist_range.first;

    dist_range = std::minmax_element(vtimes,
                                     vtimes + 4,
                                     std::less<T>());

    // Time to reach nearest vertex
    T vtimes_min = *dist_range.first;

    // Time to reach farthest vertex (last sample - any case)
    T vtimes_max = *dist_range.second;

    // Last time
    (*arrivalTimes)[2] = vtimes_max;

    bool inside = false;

    T timeLines;

    // Find boundary samples
    if ((fabs(*u) > hw) && (fabs(*v) > hh)) {
      // Outside and in the vertices
      timeLines = vtimes_min;
    } else if (*u > hw) {
      // Outside
      timeLines = etimes[3];
    } else if (*u < -hw) {
      // Outside
      timeLines = etimes[1];
    } else {
      if (*v > hh) {
        // Outside
        timeLines = etimes[0];
      } else if (*v < -hh) {
        // Outside
        timeLines = etimes[2];
      } else {
        // Inside
        inside = true;
        timeLines = etimes_min;
      }
    }

    if (inside) {
      // Arrival time to transducer plane (if inside)
      (*arrivalTimes)[0] = timePlane;
    } else {
      (*arrivalTimes)[0] = timeLines;
    }
    (*arrivalTimes)[1] = timeLines;

    return 0;
  }

  template <>
  inline int
  calcProjectionAndLimits(const sysparm_t<double>& __restrict sysparm,
                          const sps::element_rect_t<double>& __restrict element,
                          const sps::point_t<double>& __restrict point,
                          const double& __restrict delay,
                          proj_limit_dist_t<double>* __restrict limits)
  {
    double arrivalTimes[3];
    double vdists[4];
    calcProjectionAndBoundaries(sysparm,
                                element,
                                point,
                                delay,
                                &limits->u,
                                &limits->v,
                                &limits->dist2plane,
                                &arrivalTimes,
                                &vdists);

    limits->fSampleStart = arrivalTimes[0] * sysparm.fs;
    limits->fSampleStop  = arrivalTimes[2] * sysparm.fs;
    return 0;
  }

  template <>
  inline int
  calcProjectionAndLimits(const sysparm_t<float>& __restrict sysparm,
                          const sps::element_rect_t<float>& __restrict element,
                          const sps::point_t<float>& __restrict point,
                          const float& __restrict delay,
                          proj_limit_dist_t<float>* __restrict limits)
  {
    v4f v_vuvu;
    v4f v_dist2plane;
    v4f v_hwhw;
    v4f v_vdists;

    float fSampleStart = 0.0f;
    float fSampleStop = 0.0f;

    calcProjectionAndIntegrationLimitsSIMD(sysparm, element, point, delay,
                                           &v_vuvu.v,
                                           &v_dist2plane.v,
                                           &v_hwhw.v,
                                           &v_vdists.v,
                                           &fSampleStart,
                                           &fSampleStop);
#if 0
    debug_print("u: %f, v: %f, d2p: %f, h: %f, w: %f\n",
                v_vuvu.f32[1],v_vuvu.f32[0], v_dist2plane.f32[0], v_hwhw.f32[1], v_hwhw.f32[0]);
    debug_print("d[0]: %f,d[1]: %f,d[2]: %f,d[3]: %f\n",
                v_vdists.f32[0], v_vdists.f32[1], v_vdists.f32[2], v_vdists.f32[3]);
    debug_print("fs: %f, c: %f, delay: %f\n",
                sysparm.fs, sysparm.c, delay);
#endif
    // Stupid
    limits->u = v_vuvu.f32[1];
    limits->v = v_vuvu.f32[0];
    limits->dist2plane = _mm_cvtss_f32(v_dist2plane.v);
#if 0
    debug_print("fSampleStop: %f, fSampleStop: %f\n", fSampleStart, fSampleStop);
#endif
    limits->fSampleStart = fSampleStart;
    limits->fSampleStop  = fSampleStop;
    _mm_store_ps(&(limits->vdists[0]),v_vdists.v);

    return 0;
  }

}

/*
  The streaming load intrinsic (mm_stream_load_si128) performs the load
  "using a non-temporal memory hint" (according to the Intel Intrinsics
  Guide). This means that the value loaded will not cause anything to be
  evicted from the cache.

  This is useful if you are assembling a lot of data together that you
  are going to operate on immediately and not look at again for a "long"
  time. Most commonly this happens during streaming operations. I have
  used it when I know I am performing a simple operation on a large data
  set, where I know the data will quickly get evicted from the cache
  anyway. Operations such as memcpy also fall under this category.

  The non-streaming load (mm_load_si128) will retrieve the value and it
  will be subject to normal caching rules. It may evict old cache
  entries if needed, and will be able to be retrieved from the cache
  until it is evicted.

  If you expect to use the data again before a normal cache eviction
  would occur, then the non-streaming load is preferred. If you are
  operating on a large data set where a given piece of data is not
  expected to be accessed again before it would have been kicked out of
  the cache, the streaming load is preferred.

*/

/* Local variables: */
/* indent-tab-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* indent-tabs-mode: nil */
/* End: */
