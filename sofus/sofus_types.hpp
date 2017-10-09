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

#include <cstddef>

//! SOFUS interfaces and implementations
namespace sofus {
  typedef SOFUS_TimeDomainCalcTypeNS::Value TimeDomainCalcType;

  typedef SOFUS_TimeDomainIntOrderNS::Value TimeDomainIntOrder;

  typedef SOFUS_ImpulseTypeNS::ImpulseType_Value ImpulseType;

  typedef SOFUS_ExcitationTypeNS::ExcitationType_Value ExcitationType;

  template <typename T>
  struct SOFUS_EXPORT sysparm_t {
    T c;               ///< Propagation speed
    T fs;              ///< Sample frequency
    T att;             ///< Freq. dependent att
    T beta;            ///< Freq. independent att
    T f0;              ///< Center frequency
    bool use_att;      ///< Use attenuation
    bool soft_baffle;  ///< Soft baffle boundary
    bool normalize;    ///< Normalize
    size_t nThreads;   ///< Number of threads

    TimeDomainCalcType timeDomainCalcType;
    TimeDomainIntOrder timeDomainIntOrder;
  };

  // Projection, Limits and Distances (plane and vertices).
  // Spherical: u,  v, dist2plane,                dist_vertices, fStart, fStop
  // Plane:    d1, d2, dist2plane(used for soft), not needed,    fCenter

  // TODO: Make functions filling out this structure for plane and spherical responses
#ifdef __GNUG__
  template <typename T>
  ALIGN16_BEGIN struct ALIGN16_END SOFUS_EXPORT proj_limit_dist_t
#else
  template <typename T>
  struct SOFUS_EXPORT proj_limit_dist_t
#endif
  {
    ALIGN16_BEGIN T vdists[4] ALIGN16_END; ///< Distances to vertices [m]
    T u;           ///< Projection of point onto element (u-coordinate) or position vector onto basis vector u [m]
    T v;           ///< Projection of point onto element (v-coordinate) or position vector onto basis vector v [m]
    T dist2plane;  ///< Used for computing radii or soft-baffle for planar response [m]
    T fSampleStart; ///< Sample start
    T fSampleStop;  ///< Sample stop
  };


}
/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
