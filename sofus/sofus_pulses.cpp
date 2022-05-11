#include <sofus/config.h>

#include <sofus/sofus_pulses.hpp>
#include <sps/math.h>
#include <sps/memory>
#include <sps/debug.h>
#include <string.h>

namespace sofus {

// Keep this false. We get issues with FFT if not
template <class T>
bool AperturePulses<T>::bInPlaceScale = false;

template <class T>
const T AperturePulses<T>::temporalReferenceLevel = T(-120.0);

template <class T>
const T AperturePulses<T>::bandWidthReferenceLevel = T(-6.0);

template <class T>
void AperturePulses<T>::FsSet(const T& value) {
  this->m_fs = value;

  // Update parametric impulse
  if (this->m_impulseType == sofus::ImpulseType::ImpulseTypeGaussian) {
    this->GaussPulseSet(value, this->m_f0, this->m_impulseBandwidth);
  }
}

template <class T>
void AperturePulses<T>::BandWidthSet(const T& value) {
  // Update parametric impulse
  if (this->m_impulseType == sofus::ImpulseType::ImpulseTypeGaussian) {
    this->GaussPulseSet(this->m_fs, this->m_f0, value);
  }
}

template <class T>
void AperturePulses<T>::InPlaceScale(bool enable) {
  AperturePulses<T>::bInPlaceScale = enable;
}

template <class T>
void AperturePulses<T>::TemporalPulseCalc(
  const AperturePulses<T>* pXmtPulse,
  const AperturePulses<T>* pRcvPulse,
  sps::msignal1D<T>* temporal,
  T* factor) {
  sps::msignal1D<T> xmt_temporal = sps::msignal1D<T>();
  sps::msignal1D<T> rcv_temporal = sps::msignal1D<T>();
  T xmt_scale = T(1.0);
  T rcv_scale = T(1.0);

  pXmtPulse->TemporalPulseGet(TransmissionType::Transmit, &xmt_temporal, &xmt_scale);
  pRcvPulse->TemporalPulseGet(TransmissionType::Receive, &rcv_temporal, &rcv_scale);

  *factor = xmt_scale * rcv_scale;

  if (xmt_temporal.m_data) {
    // Transmit response
    if (rcv_temporal.m_data) {
      // Transmit and receive response
      if (AperturePulses<T>::bInPlaceScale) {
        sps::mconv_fft_fs<T>(pXmtPulse->m_fs, xmt_temporal,
                             rcv_temporal, *temporal);
      } else {
        // We postponse a factor 1/fs
        *factor *= T(1.0) / pXmtPulse->m_fs;
        sps::mconv_fft<T>(xmt_temporal, rcv_temporal, *temporal);
      }
    } else {
      // Only transmit response
      *temporal = xmt_temporal;
    }
  } else {
    if (rcv_temporal.m_data) {
      // Only receive response
      *temporal = rcv_temporal;
    }
  }
}


template <class T>
AperturePulses<T>::AperturePulses() : 
  m_impulseType(ImpulseType::ImpulseTypeNonParametric),
  m_excitationType(ExcitationType::ExcitationTypeNonParametric) {
  this->m_f0 = 1e6;
  this->m_fs = 1e6;
}

// HanningWeightedPulseSet

template <class T>
void AperturePulses<T>::ToneBurstSet(const T& fs, const T& f0,
                                     const size_t& nCycles) {
  this->m_fs = fs;
  this->m_f0 = f0;
  this->m_nCycles = nCycles;
  this->m_excitationType = ExcitationType::ExcitationTypeToneBurst;

  size_t nData = (size_t) round(nCycles * fs / f0);

  excitation.ndata  = nData;
  excitation.offset = 0;
  excitation.m_data   = NULL;

  if (nData > 0) {

    auto data = sps::unique_aligned_array_create<T>(nData);

    T td = (T(2.0) * T(M_PI) * T(nCycles)) / T(nData);

    for (size_t i = 0 ; i < nData ; i++) {
      data[i] = sin(T(i)*td);
    }

    excitation.m_data  =
    std::shared_ptr<T>( (T*) _mm_malloc(sizeof(T)*nData,16), [=](T *p) {
      _mm_free(p);
    });
    memcpy((void*)excitation.m_data.get(), data.get(), nData*sizeof(T));
  }
}

// TODO(JMH): 10-90% rise time, bw = 0.35 / RT
template <class T>
void AperturePulses<T>::DeltaPulseSet(const T& RT) {
  // Rise-time
  SPS_UNREFERENCED_PARAMETER(RT);
  this->impulse.offset = 0;
  this->impulse.m_data = NULL;
  this->impulse.ndata = 1;

  this->impulse.m_data  =
  std::shared_ptr<T>( (T*) _mm_malloc(sizeof(T)*1,16), [=](T *p) {
    _mm_free(p);
  });
  this->impulse.ndata = 1;
  this->impulse.m_data.get()[0] = T(1.0);
}

template <class T>
void AperturePulses<T>::GaussPulseSet(const T& fs, const T& f0, const T& bw) {

  // Update members
  this->m_impulseBandwidth = bw;
  this->m_fs               = fs;
  this->m_f0               = f0;
  this->m_impulseType      = ImpulseType::ImpulseTypeGaussian;

  const T fc = f0;
  const T invfs = T(1.0) / fs;

  debug_print("fs: %f, f0: %f, bw: %f\n", fs, f0, bw);

  T r = pow(T(10.0), (sofus::AperturePulses<T>::bandWidthReferenceLevel / T(20.0))); // Ref level (fraction of max peak)
  T fv = - SQUARE(bw)*SQUARE(fc) / (T(8.0) * log(r));                                // Variance is fv, mean is fc

  size_t nData = 1;
  sps::unique_aligned_array<T> impulse_response;

  // Bandwidth greater than 1%
  if (bw > 0.005) {
    // Determine corresponding time - domain parameters :
    T tv = T(1.0) / (T(4.0) * SQUARE(T(M_PI)) * fv);                // variance is tv, mean is 0

    // Pulse cutoff time
    T tpr = sofus::AperturePulses<T>::temporalReferenceLevel;
    T delta = pow(T(10.0), (tpr / T(20.0)));                          // Ref level(fraction of max peak)
    T tc = sqrt(- T(2.0) * tv * log(delta));                          // Pulse cutoff time

    debug_print("tv: %f, tpr: %f, tc: %f\n", tv,tpr,tc);

    nData = (size_t) (T(2.0) * tc * fs) + 1U;

    debug_print("nData: %zu\n", nData);
    impulse_response = sps::unique_aligned_array_create<T>(nData);

    for (size_t i = 0; i < nData; i++) {
      T t = -tc + T(i) * invfs;
      impulse_response[i] = exp(- SQUARE(t) / (T(2.0)*tv)) * cos(T(2.0) * T(M_PI) * fc * t); // In-phase signal
    }
  } else {
    impulse_response = sps::unique_aligned_array_create<T>(nData);
    impulse_response[0] = T(1.0);
    nData = 0;
  }

  this->impulse.offset = 0;
  this->impulse.m_data = NULL;
  this->impulse.ndata = 0;

  if (nData > 0) {
    this->impulse.m_data  =
    std::shared_ptr<T>( (T*) _mm_malloc(sizeof(T)*nData,16), [=](T *p) {
      _mm_free(p);
    });
    this->impulse.ndata = nData;
    memcpy((void*)this->impulse.m_data.get(), impulse_response.get(), nData *sizeof(T));
  }
}

template <class T>
int AperturePulses<T>::TemporalPulseGet(const sofus::TransmissionType& transmissionType,
                                        sps::msignal1D<T>* temporal, T* factor) const {
  int err = 0;
  *factor = T(1.0);

  const size_t nImpulseSamples    = this->impulse.ndata;

  // Applies to both transmit and receive
  if ( nImpulseSamples == 1 ) {
    *factor *= (this->impulse[0] / this->m_fs);
  }

  if (transmissionType == TransmissionType::Transmit) {
    const size_t nExcitationSamples = this->excitation.ndata;

    if ( nExcitationSamples == 1 ) {
      *factor *= (this->excitation[0] / this->m_fs);
    }

    if (( nExcitationSamples > 1 ) && ( nImpulseSamples > 1 )) {
      if (AperturePulses<T>::bInPlaceScale) {
        sps::mconv_fft_fs<T>(this->m_fs, this->excitation, this->impulse, *temporal);
        *factor = T(1.0);
      } else {
        sps::mconv_fft<T>(this->excitation, this->impulse, *temporal);
        *factor = T(1.0) / this->m_fs;
      }
    } else if ( nImpulseSamples > 1 ) {
      *temporal = this->impulse;
    } else if (nExcitationSamples > 1) {
      *temporal = this->excitation;
    }
  } else if (transmissionType == TransmissionType::Receive) {
    if ( nImpulseSamples > 1 ) {
      *temporal = this->impulse;
    }
  }
  return err;
}

template class AperturePulses<float>;
template class AperturePulses<double>;
}
