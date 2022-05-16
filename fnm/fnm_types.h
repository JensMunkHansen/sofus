/**
 * @file   fnm_types.h
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Sun Apr 23 14:56:32 2017
 *
 * @brief  ANSI-C types for interface
 *
 * Copyright 2017 Jens Munk Hansen
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

#include <fnm/fnm_export.h>
#include <fnm/config.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
  Enums are only visible for MATLAB if included without recursion
  and unfortunately not if inside a struct
*/
enum MyMatlabEnum2 {
  First = 0,
  Last = 1
};

#pragma pack(4)
typedef struct ApertureProperty {
  int id;
  char name[32];
} ApertureProperty;
#pragma pack()

struct FNM_EXPORT FNM_ExcitationTypeNS {
  enum ExcitationType_Value {
    ExcitationTypeNonParametric        = 0x00, ///< Specify time-domain excitation
    ExcitationTypeDelta                = 0x01, ///< Delta pulse
    ExcitationTypeToneBurst            = 0x02, ///< Excitation computed from f0 and nCycles
    ExcitationTypeHanningWeightedPulse = 0x03, ///< Hanning weigted excitation computed from f0 and nCycles
    ExcitationTypeCount                = 0x04, ///< Unused
  } _ExcitationType_Value;
};

struct FNM_EXPORT FNM_FocusingTypeNS {
  enum FocusingType_Value {
    Rayleigh          = 0x00, ///< Rayleigh integral is solved to fix phase of complex signal
    Pythagorean       = 0x01, ///< Distance from center of element is used to fix phase
    Delays            = 0x02, ///< Explicit delays
    FocusingTypeCount = 0x03, ///< Used for invalid focusing type
  } _FocusingType_Value;
};

// TODO(JMH): Consider adding an enum outside the struct also
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
    NDivW              = 0x10,  ///< Number of width abcissas
    NDivH              = 0x11,  ///< Number of height abcissas
    Fs                 = 0x12,
#if FNM_PULSED_WAVE
    Normalize          = 0x13,
    Excitation         = 0x14,
    Impulse            = 0x15,
    Focus2             = 0x16,
    RwParamTypeCount   = 0x17,
#else
    RwParamTypeCount   = 0x12,
#endif
  } _RwParamType_Value;
};

struct FNM_EXPORT FNM_TypeNS {
  enum Type_Value {
    Float          = 0x00,  ///< Float type
    Bool           = 0x01,
    Int32          = 0x02,  ///< Integer type (32-bit)
    SizeT          = 0x03,
    Struct         = 0x04,
    TypeValueCount = 0x05,
  } _Type_Value;
};

struct FNM_EXPORT FNM_ApodizationTypeNS {
  enum ApodizationType_Value {
    ApodizationTypeNonParametric = 0x00, ///< Specify element-wise apodization
    ApodizationTypeRectangular   = 0x01, ///< Rectangular window using an f-number
    ApodizationTypeHamming       = 0x02, ///< Hamming window
    ApodizationTypeCount         = 0x03, ///< Unused
  } _ApodizationType_Value;
};

#ifdef __cplusplus
}
#endif

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
