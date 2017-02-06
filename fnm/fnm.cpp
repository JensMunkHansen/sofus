/**
 * @file   fnm.cpp
 * @author Jens Munk Hansen <jmh@jmhlaptop>
 * @date   Sun Jun 12 18:57:55 2016
 *
 * @brief  Source file for FNM utilities
 *
 *
 */

// TODO: Reuse results (skipped),
//       Limits independent of z-coordinate (could be re-used)
//       Prevent one thread to eat all messages (or send index of data)
//       Support arbitrary number of threads (done)
//       Use static or unnamed namespace to use internal linkage (done)

#include <fnm/fnm.hpp>
#include <fnm/config.h>
#include <fnm/FnmSIMD.hpp> // Used by CalcCwFocus
#include <fnm/fnm_data.hpp>
#include <fnm/fnm_types.hpp>

#include <fnm/fnm_common.hpp>
#include <fnm/fnm_calc.hpp>
#include <fnm/fnm_arrays.hpp>

#if FNM_PULSED_WAVE
# include <sofus/sofus_types.hpp>
# include <sofus/sofus_calc.hpp>
#endif

#include <sps/memory>            // sps::deleted_aligned_array_create
#include <sps/smath.hpp>
#include <sps/sps_threads.hpp>
#include <sps/sps_mqueue.hpp>

#include <sps/mm_malloc.h>
#include <sps/cerr.h>
#include <sps/profiler.h>
#include <sps/extintrin.h>
#include <sps/trigintrin.h>
#include <sps/debug.h>
#include <sps/multi_malloc.hpp>

#include <sps/progress.hpp>

#include <gl/gl.hpp>

#include <string.h>
#include <new>
#include <stdexcept>

#include <assert.h>

#include <fftw3.h>

#ifdef _MSC_VER
#include <BaseTsd.h>
#endif

#ifndef MSGMAX
# define MSGMAX 8192
#endif

namespace fnm {

  /////////////////////////////////////////////////
  // Static content declared in interface
  /////////////////////////////////////////////////

#if defined(C99) && !defined(__STRICT_ANSI__)
  template <class T>
  sysparm_t<T> Aperture<T>::_sysparm =
    (sysparm_t<T>)
  {
    .c = 1500.0, .nDivW = 16, .DivH = 16,
     .att = 0.0,
      .beta = 0.0,
       .use_att = false
#if FNM_PULSED_WAVE
                  ,
                  .timeDomainCalcType = sofus::TimeDomainCalcType::Sphere,
                   .pulseWaveIntOrder  = sofus::PulsedWaveIntOrder::Fourth
#endif
  };
#else
  template <class T>
  sysparm_t<T> Aperture<T>::_sysparm = {T(1500.0), size_t(16), size_t(16),
                                        T(0.0),
                                        T(0.0),
                                        bool(false)
#if FNM_PULSED_WAVE
                                        ,
                                        sofus::TimeDomainCalcType(sofus::TimeDomainCalcType::Sphere),
                                        sofus::PulsedWaveIntOrder(sofus::PulsedWaveIntOrder::Fourth)
#endif
                                       };
#endif

  template <class T>
  const T Aperture<T>::Neper_dB = T(0.11512925464970231);

  template <class T>
  const T Aperture<T>::dB_Neper = T(8.685889638065035);

#if FNM_PULSED_WAVE
  template <class T>
  T Aperture<T>::fs = 100e6;
#endif
  
  template <class T>
  size_t Aperture<T>::nthreads = 4;

  ////////////////////////////////
  // Implementation of interface
  ////////////////////////////////
  template <class T>
  Aperture<T>::Aperture()
  {
    // Aligned data segment
    m_data = (ApertureData<T>*) _mm_malloc(sizeof(ApertureData<T>),16);
    new (this->m_data) ApertureData<T>();

#ifdef USE_PROGRESS_BAR
    // Progress bar
    m_pbar = NULL;
#endif

    // Data needed for time-domain pulsed waves
#ifdef FNM_PULSED_WAVE
    m_pulses = new sofus::AperturePulses<T>();
#endif
  }

  template <class T>
  Aperture<T>::Aperture(const size_t nElements, const T width, const T kerf, const T height) : Aperture<T>()
  {
    // If I set m_data->m_elements, it works
    const size_t nSubElements = 1;

    if (nElements * nSubElements > 0) {
      auto elements = sps::deleted_aligned_multi_array<sps::element_t<T>, 2>(nElements, nSubElements);

      T wx = (T(nElements)-T(1))/2;

      T pitch = width + kerf;

      for (size_t i = 0 ; i < nElements ; i++) {
        elements[i][0].center[0] = (T(i) - wx)*pitch;
        elements[i][0].center[1] = T(0.0);
        elements[i][0].center[2] = T(0.0);
        memset((void*)&elements[i][0].euler,0,sizeof(sps::euler_t<T>));
        elements[i][0].hh           = height/2;
        elements[i][0].hw           = width/2;
        elements[i][0].euler.alpha  = T(0.0);
        elements[i][0].euler.beta   = T(0.0);
        elements[i][0].euler.gamma  = T(0.0);
      }

      m_data->ElementsSet(std::forward<sps::deleted_aligned_multi_array<sps::element_t<T>,2> >(elements),
                          nElements,
                          nSubElements);
    }
  }

  template <class T>
  Aperture<T>::Aperture(const size_t nElements, const T width, const T kerf, const T height,
                        const size_t nSubH, const T focus) : Aperture<T>()
  {

    const size_t nSubW = 1;
    const size_t nSubElements = nSubH * nSubW;

    if (nElements * nSubElements > 0) {
      auto elements = sps::deleted_aligned_multi_array<sps::element_t<T>, 2>();

      T pitch = width + kerf;

      FocusedLinearArray(nElements,nSubH,nSubW,
                         pitch, kerf, height,
                         focus,
                         0,
                         std::move(elements));

      m_data->ElementsSet(std::forward<sps::deleted_aligned_multi_array<sps::element_t<T>,2> >(elements),
                          nElements,
                          nSubElements);
    }
  }

  template <class T>
  Aperture<T>::Aperture(const size_t nElements, const T width, const T kerf, const T height,
                        const T radius,
                        const size_t nSubH, const T focus) : Aperture<T>()
  {

    const size_t nSubW = 1;
    const size_t nSubElements = nSubH * nSubW;

    if (nElements * nSubElements > 0) {
      auto elements = sps::deleted_aligned_multi_array<sps::element_t<T>, 2>();

      T pitch = width + kerf;

      FocusedConvexArray(nElements,nSubH,nSubW,
                         pitch, kerf, height, radius,
                         focus,
                         0,
                         std::move(elements));

      m_data->ElementsSet(std::forward<sps::deleted_aligned_multi_array<sps::element_t<T>,2> >(elements),
                          nElements,
                          nSubElements);
    }
  }

  template <class T>
  Aperture<T>::~Aperture()
  {
    if (m_data) {
      m_data->~ApertureData();
      _mm_free(m_data);
      m_data = NULL;
    }
#ifdef FNM_PULSED_WAVE
    if (m_pulses) {
      delete m_pulses;
      m_pulses = nullptr;
    }
#endif
  }

  template <class T>
  void Aperture<T>::AlphaSet(const T& value)
  {
    Aperture<T>::_sysparm.att = value;
  }

  template <class T>
  const T& Aperture<T>::AlphaGet() const
  {
    return Aperture<T>::_sysparm.att;
  }

  template <class T>
  void Aperture<T>::BetaSet(const T& value)
  {
    Aperture<T>::_sysparm.beta = value;
  }

  template <class T>
  const T& Aperture<T>::BetaGet() const
  {
    return Aperture<T>::_sysparm.beta;
  }

  template <class T>
  void Aperture<T>::AttenuationEnabledSet(const bool& iEnabled)
  {
    Aperture<T>::_sysparm.use_att = iEnabled;
  }

  template <class T>
  const bool& Aperture<T>::AttenuationEnabledGet() const
  {
    return Aperture<T>::_sysparm.use_att;
  }

  template <class T>
  const size_t& Aperture<T>::NThreadsGet() const
  {
    return Aperture<T>::nthreads;
  }

  template <class T>
  void Aperture<T>::NThreadsSet(const size_t &nThreads)
  {
    assert(nThreads <= N_MAX_THREADS);
    if (nThreads <= N_MAX_THREADS)
      Aperture<T>::nthreads = nThreads;
  }

  template <class T>
  void Aperture<T>::PositionsGet(T **out, size_t* nElements, size_t* nParams) const
  {

    const size_t _nElements          = m_data->m_nelements;
    const size_t arrSize             = _nElements*3*sizeof(T);
    const auto& positions            = m_data->m_pos;

    // The function allocates
    T* arr = (T*) malloc(arrSize);
    if (arr) {
      memset(arr,0,arrSize);
      for (size_t i = 0 ; i < _nElements ; i++) {
        arr[i*3 + 0] = positions[i][0];
        arr[i*3 + 1] = positions[i][1];
        arr[i*3 + 2] = positions[i][2];
      }
    }

    // Assign outputs
    *nElements = m_data->m_nelements;
    *nParams   = 3;
    *out = arr;
  }

  template <class T>
  void Aperture<T>::ExtentGet(T** coordinates, size_t* nDim, size_t* nLimits) const
  {

    const size_t arrSize = 6;
    T* arr = (T*)malloc(arrSize*sizeof(T));
    memset(arr,0,arrSize);

    sps::bbox_t<T> bbox;
    this->m_data->ExtentGet(bbox);

    typedef T extent_t[3][2];
    extent_t* extent = (extent_t*) arr;

    for(size_t iDim = 0 ; iDim < 3 ; iDim++) {
      (*extent)[iDim][0] = bbox.min[iDim];
      (*extent)[iDim][1] = bbox.max[iDim];
    }

    *coordinates = arr;
    *nDim        = 3;
    *nLimits     = 2;
  }

  template <class T>
  T Aperture<T>::AreaGet() const
  {
    return m_data->AreaGet();
  }

  template <class T>
  void Aperture<T>::PositionsSet(const T* pos, const size_t nPositions, const size_t nDim)
  {
    const size_t _nElements          = m_data->m_nelements;
    sps::point_t<T>* positions       = m_data->m_pos.get();
    if (_nElements*3 == nPositions * nDim) {
      // Consider creating and moving an unique_ptr here
      for (size_t i = 0 ; i < _nElements ; i++) {
        positions[i][0] = pos[i*3 + 0];
        positions[i][1] = pos[i*3 + 1];
        positions[i][2] = pos[i*3 + 2];
      }
    }
  }

  template <class T>
  void Aperture<T>::FocusSet(const T iFocus[3])
  {
    if ((iFocus[0] != m_data->m_focus[0]) || (iFocus[0] != m_data->m_focus[0]) || (iFocus[0] != m_data->m_focus[0])) {
      m_data->m_focus_valid = FocusingType::FocusingTypeCount;
    }
    memcpy(&m_data->m_focus[0],&iFocus[0],3*sizeof(T));
  }

  template <class T>
  void Aperture<T>::FocusUpdate()
  {

    if ((m_data->m_focus_valid == m_data->m_focus_type)) {
      debug_print("We cannot skip this if elements,c, or f0 are changed");
      return;
    }

    const size_t nElements    = m_data->m_nelements;
    const size_t nSubElements = m_data->m_nsubelements;

    if (m_data->m_focus_type == FocusingType::Pythagorean) {

      // Geometric focusing - allocate here, since we consider removing m_data->m_delays
      auto delays = sps::deleted_aligned_array_create<T>(nElements);
      T maxDelay = T(0.0);
      T minDelay =  std::numeric_limits<T>::max();
      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
        delays[iElement] = norm(m_data->m_pos[iElement] - m_data->m_focus) / Aperture<T>::_sysparm.c;
        maxDelay = std::max<T>(maxDelay,delays[iElement]);
        minDelay = std::min<T>(minDelay,delays[iElement]);
      }
      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
        // Verify if we can positive delays
        m_data->m_delays[iElement] = minDelay - delays[iElement];

        for (size_t jElement = 0 ; jElement < nSubElements ; jElement++) {
          m_data->m_phases[iElement*nSubElements+jElement] = T(M_2PI) * (delays[iElement] - maxDelay) * m_data->m_f0;
        }
      }
    } else if (m_data->m_focus_type == FocusingType::Rayleigh) {

      // Rayleigh focusing (CW)
      size_t nPositions = nElements;
      auto positions = sps::deleted_aligned_array_create<T>(3*nPositions);
      for (size_t iPosition = 0 ; iPosition < nPositions ; iPosition++) {
        positions[3*iPosition]   = m_data->m_focus[0];
        positions[3*iPosition+1] = m_data->m_focus[1];
        positions[3*iPosition+2] = m_data->m_focus[2];
      }

      std::complex<T>* pFieldValues = NULL;
      size_t nFieldValues;

      // Note, this function allocates!!!
      this->CalcCwFocus(positions.get(),nPositions,3,
                        &pFieldValues, &nFieldValues);

      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
        for (size_t jElement = 0 ; jElement < nSubElements ; jElement++) {
          m_data->m_phases[iElement*nSubElements+jElement] =
            - std::arg<T>(pFieldValues[iElement]);
        }
      }
      free(pFieldValues);
    } else if (m_data->m_focus_type == FocusingType::Delays) {
      // Update phases using delays
      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
        for (size_t jElement = 0 ; jElement < nSubElements ; jElement++) {
          m_data->m_phases[iElement*nSubElements+jElement] =
            std::fmod<T>(m_data->m_f0*m_data->m_delays[iElement], T(1.0))*T(M_2PI);
        }
      }
    }

    // Set validity
    m_data->m_focus_valid = m_data->m_focus_type;
  }

  template <class T>
  void Aperture<T>::ApodizationSet(const T* data, size_t nData)
  {
    const size_t nElements = m_data->m_nelements;

    if (nData == m_data->m_nelements) {
      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
        m_data->m_apodizations[iElement] = data[iElement];
      }
    }
  }

  template <class T>
  void Aperture<T>::ApodizationGet(T** data, size_t* nData) const
  {
    size_t nElements = m_data->m_nelements;
    *data = (T*) malloc(nElements*sizeof(T));
    if (nElements > 0) {
      memcpy(*data,m_data->m_apodizations.get(),nElements*sizeof(T));
      *nData = nElements;
    } else {
      *nData = 0;
    }
  }

  template <class T>
  void Aperture<T>::PhasesGet(T** odata, size_t* nData)
  {

    size_t nElements = m_data->m_nelements;
    size_t nSubElements = m_data->m_nsubelements;
    *odata = (T*) malloc(nElements*sizeof(T));
    if (nElements > 0) {
      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
        (*odata)[iElement] = m_data->m_phases[iElement*nSubElements];
      }
      *nData = nElements;
    } else {
      *nData = 0;
    }
  }

  template <class T>
  void Aperture<T>::SubPhasesGet(T** odata, size_t* nData)
  {

    size_t nElements = m_data->m_nelements;
    size_t nSubElements = m_data->m_nsubelements;
    *odata = (T*) malloc(nElements*nSubElements*sizeof(T));
    if (nElements > 0) {
      memcpy(*odata, m_data->m_phases.get(), nElements*nSubElements*sizeof(T));
      *nData = nElements*nSubElements;
    } else {
      *nData = 0;
    }
  }

  template <class T>
  void Aperture<T>::DelaysGet(T** odata, size_t* nData)
  {

    size_t nElements = m_data->m_nelements;
    *odata = (T*) malloc(nElements*sizeof(T));
    if (nElements > 0) {
      memcpy(*odata,m_data->m_delays.get(),nElements*sizeof(T));
      *nData = nElements;
    } else {
      *nData = 0;
    }
  }

  template <class T>
  void Aperture<T>::DelaysSet(const T* data, const size_t nData)
  {
    const size_t nElements = m_data->m_nelements;

    if (nData == m_data->m_nelements) {
      memcpy(m_data->m_delays.get(), data, sizeof(T)*nElements);
      // Set focusing type to using delays
      m_data->m_focus_type = FocusingType::Delays;
    }
  }

  template <class T>
  void Aperture<T>::FocusGet(T oFocus[3]) const
  {
    memcpy(oFocus,&m_data->m_focus[0],3*sizeof(T));
  }

  template <class T>
  int Aperture<T>::FocusingTypeGet() const
  {
    return m_data->m_focus_type;
  }

  template <class T>
  void Aperture<T>::FocusingTypeSet(const int iFocusingType)
  {
    m_data->m_focus_type = iFocusingType;
  }

  template <class T>
  const T& Aperture<T>::F0Get() const
  {
    return m_data->m_f0;
  }

  template <class T>
  void Aperture<T>::F0Set(const T f0)
  {
    // Make static variable
    const T eps = std::numeric_limits<T>::epsilon();
    assert(f0 > eps);
    if (f0 > eps) {
      // Invalidate focus
      m_data->m_focus_valid = FocusingType::FocusingTypeCount;
      m_data->m_f0 = f0;
    }
  }

  template <class T>
  const sysparm_t<T> Aperture<T>::SysParmGet() const
  {
    return Aperture<T>::_sysparm;
  }

  template <class T>
  void Aperture<T>::SysParmSet(const sysparm_t<T> *arg)
  {
    Aperture<T>::_sysparm.c     = arg->c;
    Aperture<T>::_sysparm.nDivH = arg->nDivH;
    Aperture<T>::_sysparm.nDivW = arg->nDivW;

    // Time-domain parameters

#if FNM_PULSED_WAVE
    // This is propagator
    Aperture<T>::_sysparm.timeDomainCalcType = arg->timeDomainCalcType;
    Aperture<T>::_sysparm.pulseWaveIntOrder = arg->pulseWaveIntOrder;
#endif
  }

#ifdef FNM_PULSED_WAVE
  template <class T>
  const T& Aperture<T>::FsGet() const
  {
    return Aperture<T>::fs;
  }

  template <class T>
  void Aperture<T>::FsSet(const T fs)
  {
    // Make static variable
    assert(fs > T(1.0));
    if (fs > T(1.0)) {
      Aperture<T>::fs = fs;
      m_pulses->fs = fs;
    }
  }

  template <class T>
  void Aperture<T>::ImpulseGet(T** data, size_t* nData) const
  {
    *nData = m_pulses->impulse.ndata;
    *data  = m_pulses->impulse.m_data.get();
  }

  template <class T>
  void Aperture<T>::ImpulseSet(const T* data,
                               const size_t nData)
  {
    m_pulses->impulse.ndata = nData;
    m_pulses->impulse.offset = 0;
    m_pulses->impulse.m_data = NULL;

    if (nData > 0) {
      m_pulses->impulse.m_data  = std::shared_ptr<T>( (T*) _mm_malloc(sizeof(T)*nData,16), [=](T *p) {
        _mm_free(p);
      });
      memcpy((void*)m_pulses->impulse.m_data.get(), data, nData*sizeof(T));
    }
  }

  template <class T>
  void Aperture<T>::ExcitationGet(T** data, size_t* nData) const
  {
    *nData = m_pulses->excitation.ndata;
    *data  = m_pulses->excitation.m_data.get();
  }

  template <class T>
  void Aperture<T>::ExcitationSet(const T* data,
                                  const size_t nData)
  {
    m_pulses->excitation.ndata  = nData;
    m_pulses->excitation.offset = 0;
    m_pulses->excitation.m_data   = NULL;

    if (nData > 0) {
      m_pulses->excitation.m_data  = std::shared_ptr<T>( (T*) _mm_malloc(sizeof(T)*nData,16), [=](T *p) {
        _mm_free(p);
      });
      memcpy((void*)m_pulses->excitation.m_data.get(), data, nData*sizeof(T));
    }
  }
#endif

  template <class T>
  const T& Aperture<T>::CGet()  const
  {
    return Aperture<T>::_sysparm.c;
  }

  template <class T>
  void Aperture<T>::CSet (const T c)
  {
    Aperture<T>::_sysparm.c = c;
    // Invalidate focus
    m_data->m_focus_valid = FocusingType::FocusingTypeCount;
  }

  template <class T>
  const size_t& Aperture<T>::NElementsGet() const
  {
    return m_data->m_nelements;
  }

  template <class T>
  const size_t& Aperture<T>::NSubElementsGet() const
  {
    return m_data->m_nsubelements;
  }

  template <class T>
  const size_t& Aperture<T>::NDivWGet() const
  {
    return Aperture<T>::_sysparm.nDivW;
  }

  template <class T>
  void Aperture<T>::NDivWSet(const size_t nDivW)
  {
    Aperture<T>::_sysparm.nDivW = nDivW;
  }

  template <class T>
  const size_t& Aperture<T>::NDivHGet() const
  {
    return Aperture<T>::_sysparm.nDivH;
  }

  template <class T>
  void Aperture<T>::NDivHSet(const size_t nDivH)
  {
    Aperture<T>::_sysparm.nDivH = nDivH;
  }

  template <class T>
  void Aperture<T>::RectanglesGet(T** out, size_t* nElements,
                                  size_t* nSubElements, size_t* nParams) const
  {

    const size_t nCornerCoordinates = 12;

    const size_t _nElements       = m_data->m_nelements;
    const size_t _nSubElements    = m_data->m_nsubelements; // Length of vectors
    const size_t arrSize          = _nElements * _nSubElements * nCornerCoordinates * sizeof(T);

    // The function allocates
    T* arr = (T*) malloc(arrSize);
    if (arr) {
      memset(arr,0,arrSize);

      const auto& rectangles = m_data->m_rectangles;

      for (size_t iElement = 0 ; iElement < _nElements ; iElement++) {
        for (size_t iSubElement = 0 ; iSubElement < _nSubElements ; iSubElement++) {
          for (size_t iCorner = 0 ; iCorner < Aperture<T>::nVerticesPerElement ; iCorner++) {
            for (size_t iXYZ = 0 ; iXYZ < 3 ; iXYZ++) {
              arr[iElement * _nSubElements * Aperture<T>::nVerticesPerElement * 3 +
                  iSubElement * Aperture<T>::nVerticesPerElement * 3 + iCorner*3 + iXYZ] =
                    rectangles[iElement*_nSubElements +
                               iSubElement][iCorner][iXYZ];
            }
          }
        }
      }

      *nElements    = m_data->m_nelements;
      *nSubElements = m_data->m_nsubelements;
      *nParams      = nCornerCoordinates;

      *out          = arr;
    }
  }

  template <class T>
  void Aperture<T>::ElementsGet(T** out, size_t* nElements,
                                size_t* nParams) const
  {
    // TODO: Use common function for ElementsGet and SubElementsGet

    const size_t nElePosParams = Aperture<T>::nElementPosParameters;

    const size_t _nElements             = m_data->m_nelements;
    const size_t nSubElementsPerElement = 1;
    const size_t arrSize                = _nElements * nSubElementsPerElement * nElePosParams * sizeof(T);

    // The function allocates
    T* arr = (T*) malloc(arrSize);
    if (arr) {
      memset(arr,0,arrSize);

      const auto& elements = m_data->m_elements;

      for (size_t iElement = 0 ; iElement < _nElements ; iElement++) {
        arr[iElement * nSubElementsPerElement * nElePosParams] = elements[iElement][0].hw;

        arr[iElement * nSubElementsPerElement * nElePosParams + 1] = elements[iElement][0].hh;

        memcpy(&arr[iElement * nSubElementsPerElement * nElePosParams + 2],
               &(elements[iElement][0].center[0]),
               sizeof(sps::point_t<T>));
        arr[iElement * nSubElementsPerElement * nElePosParams+5] = elements[iElement][0].euler.alpha;
        arr[iElement * nSubElementsPerElement * nElePosParams+6] = elements[iElement][0].euler.beta;
        arr[iElement * nSubElementsPerElement * nElePosParams+7] = elements[iElement][0].euler.gamma;
      }

      *nElements    = m_data->m_nelements;
      *nParams      = Aperture<T>::nElementPosParameters; // No need to introduce static variable
      *out          = arr;
    }
  }

  template <class T>
  void Aperture<T>::SubElementsGet(T** out, size_t* nElements,
                                   size_t* nSubElements, size_t* nParams) const
  {

    const size_t nElePosParams = Aperture<T>::nElementPosParameters;

    const size_t _nElements             = m_data->m_nelements;
    const size_t nSubElementsPerElement = m_data->m_nsubelements;
    const size_t arrSize                = _nElements * nSubElementsPerElement * nElePosParams * sizeof(T);

    // The function allocates
    T* arr = (T*) malloc(arrSize);

    if (arr) {
      memset(arr,0,arrSize);

      const auto& elements = m_data->m_elements;

      for (size_t iElement = 0 ; iElement < _nElements ; iElement++) {
        for (size_t jElement = 0 ; jElement < nSubElementsPerElement ; jElement++) {

          arr[iElement * nSubElementsPerElement * nElePosParams +
              jElement * nElePosParams + 0] = elements[iElement][jElement].hw;

          arr[iElement * nSubElementsPerElement * nElePosParams +
              jElement * nElePosParams + 1] = elements[iElement][jElement].hh;

          memcpy(&arr[iElement * nSubElementsPerElement * nElePosParams + jElement*nElePosParams + 2],
                 &elements[iElement][jElement].center[0],
                 sizeof(sps::point_t<T>));
          arr[iElement * nSubElementsPerElement * nElePosParams + jElement*nElePosParams+5] =
            elements[iElement][jElement].euler.alpha;
          arr[iElement * nSubElementsPerElement * nElePosParams + jElement*nElePosParams+6] =
            elements[iElement][jElement].euler.beta;
          arr[iElement * nSubElementsPerElement * nElePosParams + jElement*nElePosParams+7] =
            elements[iElement][jElement].euler.gamma;
        }
      }

      // Assign output
      *nElements    = m_data->m_nelements;
      *nSubElements = m_data->m_nsubelements;
      *nParams      = Aperture<T>::nElementPosParameters;
      *out          = arr;
    }
  }

  template <class T>
  bool Aperture<T>::ElementsSet(const T* pos, const size_t nPositions, const size_t nDim) throw (std::runtime_error)
  {
    bool retval = this->SubElementsSet(pos,nPositions,1,nDim);
    if (!retval) {
      throw std::runtime_error("Elements improperly formatted");
    }
    return retval;
  }

  // Screws up or???
  template <class T>
  bool Aperture<T>::SubElementsSet(const T* pos, const size_t nElements,
                                   const size_t nSubElementsPerElement, const size_t nDim)
  {
    bool retval = true;

    // We are not allowed to set 0 elements
    if ((nDim != Aperture<T>::nElementPosParameters) || (nSubElementsPerElement < 1)) {
      // TODO: Consider throwing
      retval = false;
      return retval;
    }

    // Invalidate focus
    m_data->m_focus_valid = FocusingType::FocusingTypeCount;

    sps::deleted_aligned_multi_array<sps::element_t<T>, 2> elements =
      sps::deleted_aligned_multi_array_create<sps::element_t<T>, 2>(nElements,nSubElementsPerElement);

    const size_t nElePosParams = Aperture<T>::nElementPosParameters;

    for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
      for (size_t jElement = 0 ; jElement < nSubElementsPerElement ; jElement++) {

        elements[iElement][jElement].hw = pos[iElement * nSubElementsPerElement * nElePosParams +
                                              jElement * nElePosParams + 0];
        // Error was here
        elements[iElement][jElement].hh = pos[iElement * nSubElementsPerElement * nElePosParams +
                                              jElement * nElePosParams + 1];

        memcpy(&elements[iElement][jElement].center[0],
               &pos[iElement * nSubElementsPerElement * nElePosParams + jElement*nElePosParams + 2],
               3*sizeof(T));
        elements[iElement][jElement].euler.alpha =
          pos[iElement * nSubElementsPerElement * nElePosParams + jElement*nElePosParams+5];
        elements[iElement][jElement].euler.beta =
          pos[iElement * nSubElementsPerElement * nElePosParams + jElement*nElePosParams+6];
        elements[iElement][jElement].euler.gamma =
          pos[iElement * nSubElementsPerElement * nElePosParams + jElement*nElePosParams+7];
      }
    }

    // Set and initialize elements
    this->m_data->ElementsSet(std::forward<sps::deleted_aligned_multi_array<sps::element_t<T>,2> >(elements),
                              nElements, nSubElementsPerElement);

    return retval;
  }

#ifdef USE_PROGRESS_BAR
  template <class T>
  void Aperture<T>::ProgressBarSet(sps::ProgressBarInterface* pbar)
  {
    this->m_pbar = pbar;
  }
#endif

#ifdef FNM_PULSED_WAVE
  // Calculations
  template <class T>
  T Aperture<T>::CalcPwField(const T* pos, const size_t nPositions, const size_t nDim,
                             T** odata, size_t* nSignals, size_t* nSamples)
  {

    assert(nDim == 3);
    if (nDim != 3) {
      *odata = NULL;
      *nSignals = 0;
      *nSamples = 0;
      return T(0.0);
    }

    if (m_data->m_focus_type == FocusingType::Rayleigh) {
      fprintf(stderr, "Rayleigh focusing type is not supported for pulsed wave calculations\n");
    }
    this->FocusUpdate();

    sofus::sysparm_t<T> sysparm;
    sysparm.c   = Aperture<T>::_sysparm.c;
    sysparm.fs  = Aperture<T>::fs;

    sysparm.att = Aperture<T>::_sysparm.att;
    sysparm.beta = Aperture<T>::_sysparm.beta;
    sysparm.use_att = Aperture<T>::_sysparm.use_att;

    T retval = sofus::CalcPwField(sysparm,
                                  m_data,
                                  m_pulses,
                                  Aperture<T>::_sysparm.timeDomainCalcType,
                                  Aperture<T>::_sysparm.pulseWaveIntOrder,
                                  pos, nPositions, nDim,
                                  odata, nSignals, nSamples,
                                  m_pbar);

    return retval;
  }

  template <class T>
  T Aperture<T>::CalcPwFieldThreaded(const T* pos, const size_t nPositions, const size_t nDim,
                                     T** odata, size_t* nSignals, size_t* nSamples)
  {

    assert(nDim == 3);
    if (nDim != 3) {
      *odata = NULL;
      *nSignals = 0;
      *nSamples = 0;
      return T(0.0);
    }

    if (m_data->m_focus_type != FocusingType::Pythagorean) {
      fprintf(stderr, "Focusing type must be Pythagorean for pulsed wave calculations\n");
    }
    this->FocusUpdate();

    sofus::sysparm_t<T> sysparm;
    sysparm.c        = Aperture<T>::_sysparm.c;
    sysparm.fs       = Aperture<T>::fs;
    sysparm.nthreads = Aperture<T>::nthreads;

    sysparm.att      = Aperture<T>::_sysparm.att;
    sysparm.beta     = Aperture<T>::_sysparm.beta;
    sysparm.use_att  = Aperture<T>::_sysparm.use_att;

    // Allocates
    T retval = sofus::CalcPwFieldThreaded(sysparm,
                                          m_data,
                                          m_pulses,
                                          Aperture<T>::_sysparm.timeDomainCalcType,
                                          Aperture<T>::_sysparm.pulseWaveIntOrder,
                                          pos, nPositions, nDim,
                                          odata, nSignals, nSamples,
                                          m_pbar);

    return retval;
  }
#endif

  template <class T>
  int Aperture<T>::CalcCwField(const T* pos, const size_t nPositions, const size_t nDim,
                               std::complex<T>** odata, size_t* nOutPositions)
  {

    assert(nDim == 3);
    if (nDim != 3) {
      *odata = NULL;
      *nOutPositions = 0;
      return -1;
    }

    size_t nScatterers = nPositions;
    *odata = (std::complex<T>*) malloc(nScatterers*sizeof(std::complex<T>));
    *nOutPositions = nScatterers;

    this->FocusUpdate();

    fnm::CalcCwField<T>(*this->m_data,
                        pos, nPositions,
                        odata);
    return 0;
  }

  // Optimal sampling for integral (reduced integration path)
  template <class T>
  int Aperture<T>::CalcCwFieldRef(const T* pos, const size_t nPositions, const size_t nDim,
                                  std::complex<T>** odata, size_t* nOutPositions)
  {

    assert(nDim == 3);
    if (nDim != 3) {
      *odata = NULL;
      *nOutPositions = 0;
      return -1;
    }

    // Allocate data
    size_t nScatterers = nPositions;
    *odata = (std::complex<T>*) malloc(nScatterers*sizeof(std::complex<T>));
    *nOutPositions = nScatterers;

    this->FocusUpdate();

    int retval = fnm::CalcCwFieldRef<T>(*this->m_data,
                                        pos, nPositions,
                                        odata);
    return retval;
  }

#ifndef FNM_DOUBLE_SUPPORT

  template <class T>
  int Aperture<T>::CalcCwFast(const T* pos,
                              const size_t nPositions,
                              const size_t nDim,
                              std::complex<T>** odata,
                              size_t* nOutPositions)
  {

    int retval;
    assert(nDim == 3);
    if (nDim != 3) {
      *odata = NULL;
      *nOutPositions = 0;
      retval = -1;
      return retval;
    }

    *odata = (std::complex<T>*) malloc(nPositions*sizeof(std::complex<T>));
    memset(*odata, 0, nPositions*sizeof(std::complex<T>));
    *nOutPositions = nPositions;

    this->FocusUpdate();

    retval = fnm::CalcCwThreaded<T>(this->m_data, pos, nPositions, odata);

    return retval;
  }

  template <class T>
  int Aperture<T>::CalcCwFocus(const T* pos, const size_t nPositions, const size_t nDim,
                               std::complex<T>** odata, size_t* nOutPositions)
  {
    int retval = 0;
    assert(nDim == 3);
    if (nDim != 3) {
      *odata = NULL;
      *nOutPositions = 0;
      retval = -1;
      return retval;
    }

    *odata = (std::complex<T>*) malloc(nPositions*sizeof(std::complex<T>));
    memset(*odata, 0, nPositions*sizeof(std::complex<T>));
    *nOutPositions = nPositions;

    retval = fnm::CalcCwFocus<T>(*this->m_data,pos,nPositions,odata);

    return retval;
  }

  // Sub-optimal integration range (Old function). TODO: Remove or move to fnm_calc
  template <class T>
  int Aperture<T>::CalcCwField2(const T* pos, const size_t nPositions, const size_t nDim,
                                std::complex<T>** odata, size_t* nOutPositions)
  {

    assert(nDim == 3);
    if (nDim != 3) {
      *odata = NULL;
      *nOutPositions = 0;
      return -1;
    }
    size_t nScatterers = nPositions;
    *odata = (std::complex<T>*) malloc(nScatterers*sizeof(std::complex<T>));
    memset(*odata, 0, nScatterers*sizeof(std::complex<T>));
    *nOutPositions = nScatterers;

    this->FocusUpdate();

    const T lambda = Aperture<T>::_sysparm.c / m_data->m_f0;

    const T k = T(M_2PI)/lambda;

    const size_t nElements = m_data->m_nelements;

    const size_t nSubElements = m_data->m_nsubelements;

    const T* apodizations = m_data->m_apodizations.get();

    T apodization;


    // Need weights and abcissa values
    size_t nDivW = Aperture<T>::_sysparm.nDivW;
    size_t nDivH = Aperture<T>::_sysparm.nDivH;

    auto uweights = sps::deleted_aligned_array<T>();
    auto vweights = sps::deleted_aligned_array<T>();
    auto uxs      = sps::deleted_aligned_array<T>();
    auto vxs      = sps::deleted_aligned_array<T>();

    CalcWeightsAndAbcissae(std::move(uxs),std::move(uweights),
                           std::move(vxs),std::move(vweights));

    sps::deleted_aligned_multi_array<sps::element_t<T>,2>& elements = m_data->m_elements;

    sps::point_t<T> projection;

    for (size_t iPoint = 0 ; iPoint < nPositions ; iPoint++) {

      __m128 vec_point = _mm_set_ps(0.0f, float(pos[iPoint*3 + 2]),float(pos[iPoint*3 + 1]),float(pos[iPoint*3]));

      std::complex<T> final = std::complex<T>(T(0.0),T(0.0));

      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {

        apodization = apodizations[iElement];

        if (apodization != 0.0) {

          for (size_t jElement = 0 ; jElement < nSubElements ; jElement++) {

            const sps::element_t<T>& element = elements[iElement][jElement];

            std::complex<T> result;

            // TODO: Move element_t to smath.hpp and make projection(point_t<T>, element_t<T>)->point_t<T>
            assert(((uintptr_t)&element.center[0] & 0x0F) == 0 && "Data must be aligned");
            __m128 vec_r2p = _mm_sub_ps(vec_point, _mm_load_ps((float*)&element.center[0]));
            _mm_store_ss((float*)&projection[0],_mm_dp_ps(_mm_load_ps(&element.uvector[0]), vec_r2p,0x71));
            _mm_store_ss((float*)&projection[1],_mm_dp_ps(_mm_load_ps(&element.vvector[0]), vec_r2p,0x71));
            _mm_store_ss((float*)&projection[2],_mm_fabs_ps(_mm_dp_ps(_mm_load_ps(&element.normal[0]), vec_r2p,0x71)));

            result = CalcHzAll<T>(element, projection, k,
                                  uxs.get(), uweights.get(), nDivW,
                                  vxs.get(), vweights.get(), nDivH);
            final = final + apodization *
                    result * exp(std::complex<T>(0,m_data->m_phases[iElement*nSubElements+jElement]));
          }
        }
      }
      (*odata)[iPoint] = final;
    }
    return 0;
  }

  template <class T>
  void Aperture<T>::ManagedAllocation(std::complex<T>** outTest, size_t* nOutTest)
  {
    const size_t _nData = 100000000;
    *nOutTest = _nData;
    *outTest = (std::complex<T>*) malloc(_nData*sizeof(std::complex<T>));
  }

  template <class T>
  int KSpace(const T dx, const size_t Nx, const T dy, const size_t Ny, const T lambda, T** data, bool shift=true)
  {
    int retval = 0;

    T k2max = SQUARE(lambda / (T(2.0)*dx)) + SQUARE(lambda / (T(2.0)*dy));

    retval = (k2max <= T(1.0)) ? 0 : -1;

    if (retval == 0) {
      T k = T(M_2PI)/lambda;
      T kyDelta = lambda/T(Ny*dy);
      T kxDelta = lambda/T(Nx*dx);

      if (shift) {
        T kyStart = floor(-T(Ny)/T(2));
        T kxStart = floor(-T(Nx)/T(2.0));
        size_t iix, iiy;
        for (size_t iy = 0 ; iy < (Ny+1)/2 ; iy++) {
          T ky2 = SQUARE((kyStart + iy) * kyDelta);
          for (size_t ix = 0; ix < (Nx+1)/2 ; ix++) {
            T kx = (kxStart + T(ix)) * kxDelta;
            data[iy][ix] = k * sqrt(T(1.0) - (SQUARE(kx)+ky2));
          }
          iix = 0;
          for (size_t ix = (Nx+1)/2 ; ix < Nx ; ix++) {
            T kx = T(iix) * kxDelta;
            data[iy][ix] = k * sqrt(T(1.0) - (SQUARE(kx)+ky2));
            iix++;
          }
        }

        iiy = 0;
        for (size_t iy = (Ny+1)/2 ; iy < Ny ; iy++) {
          T ky2 = SQUARE(T(iiy) * kyDelta);
          for (size_t ix = 0; ix < (Nx+1)/2 ; ix++) {
            T kx = T(kxStart + ix) * kxDelta;
            data[iy][ix] = k * sqrt(T(1.0) - (SQUARE(kx)+ky2));
          }
          iix = 0;
          for (size_t ix = (Nx+1)/2 ; ix < Nx ; ix++) {
            T kx = T(iix) * kxDelta;
            data[iy][ix] = k * sqrt(T(1.0) - (SQUARE(kx)+ky2));
            iix++;
          }
          iiy++;
        }
      }
    }
    return retval;
  }

  template <class T>
  int Aperture<T>::CalcAsa(const T* y0, const size_t nx, const size_t ny,
                           // Vector of z coordinates
                           const T dx, const T dy,
                           const size_t Nx, const size_t Ny,
                           // Output
                           std::complex<T>** p1, size_t *onx, size_t *ony, size_t *onz)
  {
    int retval = 0;

    const T c    = Aperture<T>::_sysparm.c;
    const T f0   = this->m_data->m_f0;

    /* Propagator 'P' */
    bool bRestrict = true;

    /* Attenuation */
    const T beta = T(0.5);
    const T dBperNeper = T(20.0) * log10(exp(T(1.0)));
    const T eps = std::numeric_limits<T>::epsilon();
    const T alpha = beta / dBperNeper * T(100.0) * f0 / T(1.0e6);

    const T lambda = c/f0;

    // Write zero output
    *onx = Nx;
    *ony = Ny;
    *onz = 1;
    *p1 = (std::complex<T>*) malloc(Nx*Ny*1*sizeof(std::complex<T>));

    SPS_UNREFERENCED_PARAMETERS(c,f0,beta,bRestrict,dBperNeper,eps,lambda,alpha);

    while (retval == 0) {
      if (!((dx > eps) && (dy > eps))) {
        retval = -1;
        printf("dx or dy too small\n");
        break;
      }
      if (!(Nx && ((Nx & (Nx-1))) ==0) || !(Ny && ((Ny & (Ny-1)))==0)) {
        // Must be power of two
        retval = -2;
        printf("Nx and Ny must be power of 2: %zu, %zu\n",Nx,Ny);
        break;
      }
      if (!((Nx > nx) && (Ny > ny))) {
        printf("Nx and Ny must be larger than nx and ny\n");
        retval = -3;
        break;
      }

      T k2max = SQUARE(lambda / (T(2.0)*dx)) + SQUARE(lambda / (T(2.0)*dy));
      if (k2max > 1.0) {
        retval = -4;
        printf("dx and dy too small for the given lambda\n");
        break;
      }
      //////////////////////////////////////////////
      // Fourier transform of pressure or velocity
      //////////////////////////////////////////////

      // TODO: Use r2c and use half the values
      fftwf_complex* input = (fftwf_complex*) _mm_malloc(sizeof(fftwf_complex) * Nx * Ny, 16);
      memset(input,0,sizeof(fftwf_complex) * Nx * Ny);

      for (size_t j = 0 ; j < ny ; j++) {
        for (size_t i = 0 ; i < nx ; i++) {
          input[j*Nx + i][0] = float(y0[j*nx + i]);
        }
      }
      fftwf_complex* Y0 = (fftwf_complex*) _mm_malloc(sizeof(fftwf_complex) * Nx * Ny, 16);

      fftwf_plan forward = fftwf_plan_dft_2d((int)Ny, (int)Nx, input, Y0, -1, FFTW_ESTIMATE);
      fftwf_execute_dft(forward,(fftwf_complex *) input, (fftwf_complex *) Y0);
      fftwf_destroy_plan(forward);

      /////////////////////////////////////////////////////////////////////
      // K-space coordinates: Asymmetric for both N even and odd (shiftet)
      /////////////////////////////////////////////////////////////////////

      T** kzspace = (T**) _mm_multi_malloc<T,16>(2,Nx,Ny);
      KSpace<T>(dx,Nx,dy,Ny,lambda,kzspace);

      /////////////////////////////////////////////////////////////////////
      // Propagator (update for each z value)
      /////////////////////////////////////////////////////////////////////

      fftwf_complex* Hprop = (fftwf_complex*) _mm_malloc(sizeof(fftwf_complex) * Nx * Ny, 16);
      memset(Hprop,0,sizeof(fftwf_complex) * Nx * Ny);

      for (size_t iz = 0 ; iz < 1 ; iz++) {

        v4f carg;
        v4f sarg;
        __m128 arg;
        __m128 v_z = _mm_set1_ps(float(0.001));

        // Forward: H = np.conj(np.exp(1j * z * kzspace))
        // Backward: H = np.exp(-1j*z * kzspace) * kxsq_ysq <= 1;
        for (size_t j = 0 ; j < Nx ; j++) {
          for (size_t i = 0 ; i < Ny ; i+=4) {
            arg = _mm_load_ps((float*)&(kzspace[j][i]));
            arg = _mm_neg_ps(_mm_mul_ps(arg,v_z));
            _mm_sin_cos_ps(arg,&sarg.v,&carg.v);

            // Consider shuffling
            Hprop[j*Nx+i][0]   = carg.f32[0];
            Hprop[j*Nx+i+1][0] = carg.f32[1];
            Hprop[j*Nx+i+2][0] = carg.f32[2];
            Hprop[j*Nx+i+3][0] = carg.f32[3];
            Hprop[j*Nx+i][1]   = sarg.f32[0];
            Hprop[j*Nx+i+1][1] = sarg.f32[1];
            Hprop[j*Nx+i+2][1] = sarg.f32[2];
            Hprop[j*Nx+i+3][1] = sarg.f32[3];
            // We could multiply by Y0 in here

            // H = np.conj(np.exp(1j * z * kzspace));
          }
        }
      }


      // TODO: Use _mm_mulcmplx_ps(const __m128 &a, const __m128 &b)

      // iarg = Y0 * H
      // newpress = ifft2(iarg, (Ny,Nx))
      _mm_free(Hprop);

      // Writing to misaligned memory (Works)
      //fftwf_plan plan = fftwf_plan_dft_2d(Ny, Nx, input, (fftwf_complex *)(*p1), -1, FFTW_ESTIMATE);
      //fftwf_execute_dft(plan,(fftwf_complex *) input, (fftwf_complex *) (*p1));

      // fftw(nx*realsize(ny,in,out),-1,nx*ny) // out-of-place
      // Output is upper half of complex transform


#if 1
      for (size_t j = 0 ; j < Ny ; j++) {
        for (size_t i = 0 ; i < Nx ; i++) {
          (*p1)[j*Nx + i].real(T(Y0[j*Nx + i][0]));
          (*p1)[j*Nx + i].imag(T(Y0[j*Nx + i][1]));
          //(*p1)[j*Nx + i].real(T(kzspace[j][i]));
          //(*p1)[j*Nx + i].imag(T(0.0));
        }
      }
#endif

      _mm_multi_free<T,16>(kzspace, 2);
      _mm_free(input);
      _mm_free(Y0);
      break;
    }

    return retval;
  }
#endif

#ifdef FNM_DOUBLE_SUPPORT
  template class Aperture<double>;
#endif

  template class Aperture<float>;
}

// Template instantiation


/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */

