/**
 * @file   sofus_types.h
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Sun Apr 23 15:27:06 2017
 *
 * @brief  ANSI-C types for interface
 *
 *
 */

#pragma once

#include <sofus/config.h>
#include <sofus/sofus_export.h>

#ifdef __cplusplus
extern "C" {
#endif

struct SOFUS_EXPORT SOFUS_TimeDomainCalcTypeNS {
  enum Value {
    Ray    = 0x00,
    Plane  = 0x01,
    Sphere = 0x02,
    TimeDomainCalcTypeCount = 0x03,
  } _Value;
};

struct SOFUS_EXPORT SOFUS_TimeDomainIntOrderNS {
  enum Value {
    Midpoint                = 0x00, // Not implemented
    Trapez                  = 0x01,
    Simpson                 = 0x02,
    BoolesThird             = 0x03,
    Fourth                  = 0x04,
    Fifth                   = 0x05,
    TimeDomainIntOrderCount = 0x06,
  } _Value;
};

struct SOFUS_EXPORT SOFUS_ImpulseTypeNS {
  enum ImpulseType_Value {
    ImpulseTypeNonParametric = 0x00, ///< Specify time-domain temporal impulse response
    ImpulseTypeGaussian      = 0x01, ///< Gaussian temporal response computed from f0 and bandwidth
    ImpulseTypeCount         = 0x02, ///< Unused
  } _ImpulseType_Value;
};

struct SOFUS_EXPORT SOFUS_ExcitationTypeNS {
  enum ExcitationType_Value {
    ExcitationTypeNonParametric        = 0x00, ///< Specify time-domain excitation
    ExcitationTypeToneBurst            = 0x01, ///< Excitation computed from f0 and nCycles
    ExcitationTypeHanningWeightedPulse = 0x02, ///< Hanning weigted excitation computed from f0 and nCycles
    ExcitationTypeCount                = 0x03, ///< Unused
  } _ExcitationType_Value;
};

#ifdef __cplusplus
}
#endif

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
