/**
 * @file   sofus_types.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Sat Mar 18 13:21:37 2017
 *
 * @brief  Structure types used for time-domain field simulations
 *
 *
 */
#pragma once

#include <sofus/config.h>
#include <sofus/sofus_export.h>
#include <sofus/sofus_types.h>

#include <sps/memory>
#include <cstddef>

#include <sps/smath_types.hpp>

//! SOFUS interfaces and implementations
namespace sofus {
typedef SOFUS_TimeDomainCalcTypeNS::Value TimeDomainCalcType;

typedef SOFUS_TimeDomainIntOrderNS::Value TimeDomainIntOrder;

typedef SOFUS_ImpulseTypeNS::ImpulseType_Value ImpulseType;

template <typename T>
struct SOFUS_EXPORT sysparm_t {
  sysparm_t() : c(T(1500.0)), fs(T(70e6)), att(T(0.0)), beta(T(0.0)),
    f0(T(3.5e6)), use_att(false), soft_baffle(false), normalize(true),
    nThreads(1), timeDomainCalcType(TimeDomainCalcType::Sphere),
    timeDomainIntOrder(TimeDomainIntOrder::Fourth) {}
  T c;               ///< Propagation speed
  T fs;              ///< Sample frequency
  T att;             ///< Freq. dependent att
  T beta;            ///< Freq. independent att
  T f0;              ///< Center frequency
  bool use_att;      ///< Use attenuation
  bool soft_baffle;  ///< Soft baffle boundary
  bool normalize;    ///< Normalize
  size_t nThreads;   ///< Number of threads

  TimeDomainCalcType timeDomainCalcType;  ///< Time-domain propagator type, ray, plane or sphere
  TimeDomainIntOrder timeDomainIntOrder;  ///< Order used for strip-integration
};

// Projection, Limits and Distances (plane and vertices).
// Spherical: u,  v, dist2plane,                dist_vertices, fStart, fStop
// Plane:    d1, d2, dist2plane(used for soft), not needed,    fCenter

template <typename T>
struct SOFUS_EXPORT proj_limit_dist_t : public sps::aligned<4*sizeof(T)> {
  ALIGN16_BEGIN T vdists[4] ALIGN16_END; ///< Distances to vertices [m]
  T u;           ///< Projection of point onto element (u-coordinate) or position vector onto basis vector u [m]
  T v;           ///< Projection of point onto element (v-coordinate) or position vector onto basis vector v [m]
  T dist2plane;  ///< Used for computing radii or soft-baffle for planar response [m]
  T fSampleStart; ///< Sample start
  T fSampleStop;  ///< Sample stop
};

// TODO: Check how list are created in rcu and redefine this properly
template <typename T>
class SOFUS_EXPORT focus_line_t : public sps::aligned<4*sizeof(T)> {
 public:
  ALIGN16_BEGIN sps::point_t<T> focus ALIGN16_END; ///< Focal point if not dynamic
  ALIGN16_BEGIN sps::euler_t<T> euler ALIGN16_END; ///< Direction
  ALIGN16_BEGIN sps::point_t<T> ALIGN16_END direction;
  T time;           ///< Time after which this is valid
  int dynamic[2];   ///< Dynamic focusing
  T* delays;        ///< Vector of delays
  // Cached variables
  class focus_line_t<T>* next_focus_line;
};
}  // namespace sofus

#include <ostream>
template <typename T>
std::ostream& operator<<(std::ostream& os, const sofus::sysparm_t<T>& sp) {
  os << " c: " << sp.c
     << " fs: " << sp.fs
     << " att: " << sp.att
     << " beta: " << sp.beta
     << " f0: " << sp.f0
     << " use_att: " << sp.use_att
     << " soft_baffle: " << sp.soft_baffle
     << " normalize: " << sp.normalize
     << " nThreads: " << sp.nThreads
     << " tdCalcType" << (int) sp.timeDomainCalcType
     << " tdIntOrder: " << (int) sp.timeDomainIntOrder << std::endl;
  return os;
}






/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
