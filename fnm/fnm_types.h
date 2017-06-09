/**
 * @file   fnm_types.h
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Sun Apr 23 14:56:32 2017
 *
 * @brief  ANSI-C types for interface
 *
 *
 */
#pragma once

#include <fnm/fnm_export.h>

#ifdef __cplusplus
extern "C" {
#endif

struct FNM_EXPORT FNM_FocusingTypeNS {
  enum FocusingType_Value {
    Rayleigh          = 0x00, ///< Rayleigh integral is solved to fix phase of complex signal
    Pythagorean       = 0x01, ///< Distance from center of element is used to fix phase
    Delays            = 0x02, ///< Explicit delays
    FocusingTypeCount = 0x03, ///< Used for invalid focusing type
  } _FocusingType_Value;
};

struct FNM_EXPORT RwParamTypeNS {
  enum RwParamType_Value {
    ElementDelays      = 0x00,
    AttenuationEnabled = 0x01,
    Alpha              = 0x02,
    Beta               = 0x03,
    Positions          = 0x04,
    Focus              = 0x05,
    CenterFocus        = 0x06,
    FocusingType       = 0x07,
    NThreads           = 0x08,
    F0                 = 0x09,
    W                  = 0x0A,
    SysParm            = 0x0B,
    C                  = 0x0C,
    Elements           = 0x0D,
    SubElements        = 0x0E,
    Apodization        = 0x0F,
    Fs                 = 0x10,
    Normalize          = 0x11,
    Excitation         = 0x12,
    Impulse            = 0x13,
    NDivW              = 0x14,
    NDivH              = 0x15,
    RwParamTypeCount   = 0x16,
  } _RwParamType_Value;
};

struct FNM_EXPORT FNM_TypeNS {
  enum Type_Value {
    Float          = 0x00,
    Bool           = 0x01,
    Int32          = 0x02, /* enum */
    SizeT          = 0x03,
    Struct         = 0x04,
    TypeValueCount = 0x05,
  } _Type_Value;
};

#ifdef __cplusplus
}
#endif

