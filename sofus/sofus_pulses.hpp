/**
 * @file   sofus_pulses.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Tue Oct 18 04:42:54 2016
 *
 * @brief  Used for time-domain simulations
 *
 * Copyright 2017, Jens Munk Hansen
 */

// TODO(JEM): Figure out how to expose msignal across DLL interface
//       Introduce setters for excitation and impulse (store as unique pointers)

#pragma once

#include <sofus/config.h>
#include <sofus/sofus_export.h>
#include <sofus/sofus_types.h>
#include <fnm/fnm_types.h>
#include <sps/msignals.hpp>      // Consider forward declaration

namespace sofus {

struct SOFUS_EXPORT TransmissionTypeNS {
  enum Value {
    Transmit              = 0x00,  ///< Transmitting
    Receive               = 0x01,  ///< Receiving
    TransmissionTypeCount = 0x02,
  };
};
typedef TransmissionTypeNS::Value TransmissionType;

typedef SOFUS_ImpulseTypeNS::ImpulseType_Value ImpulseType;

typedef FNM_ExcitationTypeNS::ExcitationType_Value ExcitationType;

template <class T>
class SOFUS_EXPORT AperturePulses {
 public:

  /** Reference level used for temporal responses */
  static const T temporalReferenceLevel;

  /** Bandwidth reference level, -6 dB */
  static const T bandWidthReferenceLevel;

  // TODO(JEM): Avoid statics
  static bool bInPlaceScale;

  // TODO(JEM): Avoid statics
  static void InPlaceScale(bool enable = false);

  /**
   * Compute temporal pulse for time-domain simulations. The
   * returned temporal response and scaling is the convolution and
   * product of \ref TemporalPulseGet applied to each of the pulse
   * objects.
   *
   * @param pXmtPulse  Pulse object for transmitting aperture
   * @param pRcvPulse  Pulse object for receiving aperture
   * @param temporal   Temporal response sampled at fs
   * @param factor     Scale factor to avoid convolution of scalars
   */
  static void TemporalPulseCalc(const AperturePulses<T>* pXmtPulse,
                                const AperturePulses<T>* pRcvPulse,
                                sps::msignal1D<T>* temporal,
                                T* factor);

  void DeltaPulseSet(const T& RT);


  /**
   * Set the temporal transducer response to a Gaussian pulse
   *
   * @param fs
   * @param f0
   * @param bw
   */
  void GaussPulseSet(const T& fs, const T& f0, const T& bw);


  // TODO: Consider introducing object, which can be parametric or non-parametric
  void ToneBurstSet(const T& fs, const T& f0, const size_t& nCycles);

  /**
   * Default ctor
   *
   *
   * @return
   */
  AperturePulses();

  /**
   * Default dtor
   *
   *
   * @return
   */
  ~AperturePulses() = default;

  /**
   * Compute temporal signal and possible scaling factor.
   *
   * @param transmissionType Transmit or Receive aperture
   * @param temporal         Temporal response, convolution of impulse and excitation (if transmit).
   * @param factor           Scaling factor (avoid convolutions with scalars)
   *
   * @return Error code
   */
  int TemporalPulseGet(const TransmissionType& transmissionType,
                       sps::msignal1D<T>* temporal, T* factor) const;

  void BandWidthSet(const T& value);

  void FsSet(const T& fs);
 private:
 public:
  ImpulseType m_impulseType;        ///< Parametric or non-parametric impulse
  ExcitationType m_excitationType;  ///< Parametric or non-parametric excitation
  size_t m_nCycles;                 ///< Number of excitation cycles at f0
  T m_fs;                           ///< Sampling frequency
  T m_f0;                           ///< Center frequency
  T m_impulseBandwidth;             ///< Bandwidth of temporal response
  sps::msignal1D<T> excitation;     ///< Excitation
  sps::msignal1D<T> impulse;        ///< Temporal impulse response of transducer
};

}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
