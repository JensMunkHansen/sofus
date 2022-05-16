/**
* @file   fnm.cpp
* @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
* @date   Thu Feb 15 18:44:28 2018
*
* @brief
*
*  This file is part of SOFUS.
*
*  SOFUS is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  SOFUS is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with SOFUS. If not, see <http://www.gnu.org/licenses/>.
*
* Copyright 2017 Jens Munk Hansen
*/

#include <assert.h>
#include <fftw3.h>
#include <sps/mm_malloc.h>
#include <sps/memory>           // sps::unique_aligned_array_create
#include <sps/cerr.h>
#include <sps/profiler.h>

#include <fnm/config.h>

// Temporarily included for CalcAsa
#include <sps/trigintrin.h>
#include <sps/debug.h>
#include <sps/smath_types.hpp>

#include <fnm/fnm.hpp>
#include <fnm/fnm_data.hpp>
#include <fnm/fnm_types.hpp>

#include <fnm/fnm_common.hpp>     // fnm::CalcWeightsAndAbcissae
#include <fnm/fnm_arrays.hpp>     // fnm::FocusedLinearArray, fnm::FocusedConvexArray

#include <fnm/fnm_delays.hpp>     // fnm::PytDelays

#include <fnm/fnm_calc.hpp>       // fnm::Calc routines

# include <fnm/fnm_transient.hpp>  // fnm::CalcFd routines
# include <fnm/fnm_transient_threaded.hpp>  // fnm::CalcFd routines

#include <fnm/fnm_basis.hpp>

#if FNM_PULSED_WAVE
# include <fnm/fnm_echo.hpp>
# include <sofus/sofus_types.hpp>
# include <sofus/sofus_calc.hpp>
# include <sofus/sofus_pulses.hpp>
#endif

#include <new>
#include <stdexcept>
#include <algorithm>

// Closure test
#include <array>
#include <tuple>

#include <sps/sps_threads.hpp>  // getncpus()
#include <sps/smath.hpp>

#include <sps/multi_malloc.hpp>
#include <sps/globals.hpp>

#include <sps/progress.hpp>

#include <gl/gl.hpp>

#ifdef _MSC_VER
#include <BaseTsd.h>
#endif

// __declspec(dllimport) is recognized as a synonym for __attribute__ ((dllimport))  for CYGWIN

/*
  A way to provide static templated structures without introducing
  singletons to the proving interface. The pointer is an atomic and
  is initialized using acquire-release semantics See <sps/globals.hpp> for
  details.
 */
template <typename T, template <typename> class V>
std::atomic<V<T>*> sps::globalstruct<T, V >::pVar{nullptr};


#if defined(__GNUC__)
# if !defined(__CYGWIN__)
#  include <sps/strace.hpp>
# endif

void SofusInit()     __attribute__((constructor(101)));
void SofusDestroy()  __attribute__((destructor(101)));

#elif defined(_WIN32)

void SofusInit();
void SofusDestroy();

#endif

void SofusInit() {
  printf("\nInitializing SOFUS\n\n");
  // TODO(JEM): Avoid deprecated __malloc_hook__
#if defined(__GNUC__) && !defined(__CYGWIN__)
  // #  if !defined(NDEBUG)
# if defined(SPS_STRACE)
  sps::STrace::Instance().Enable();
# endif
#endif
}

void SofusDestroy() {
# ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning(disable : 4099)
# endif

  printf("\nFreeing SOFUS\n\n");
  if (sps::globalstruct<float, fnm::sysparm_t>::pVar) {
    delete sps::globalstruct<float, fnm::sysparm_t>::pVar;
  }
# ifdef FNM_DOUBLE_SUPPORT
  if (sps::globalstruct<double, fnm::sysparm_t>::pVar) {
    delete sps::globalstruct<double, fnm::sysparm_t>::pVar;
  }
# endif

# ifdef _MSC_VER
#  pragma warning(pop)
# endif
}

#if defined(_WIN32)

BOOL APIENTRY DllMain(HANDLE hModule,
                      DWORD  ulReasonForCall,
                      LPVOID lpReserved) {
  SPS_UNREFERENCED_PARAMETERS(hModule, ulReasonForCall, lpReserved);

  switch (ulReasonForCall) {
  case DLL_PROCESS_ATTACH:
    SofusInit();
    break;
  case DLL_THREAD_ATTACH:
    break;
  case DLL_THREAD_DETACH:
    break;
  case DLL_PROCESS_DETACH:
    SofusDestroy();
    break;
  }
  return TRUE;
}
#endif

namespace fnm {

std::mutex myMutex;  ///< Mutex used in this compilation unit

/////////////////////////////////////////////////
// Constant content declared in interface
/////////////////////////////////////////////////

template <class T>
const T Aperture<T>::Neper_dB = T(0.11512925464970231);

template <class T>
const T Aperture<T>::dB_Neper = T(8.685889638065035);

/////////////////////////////////////////////////
// Static content declared in interface
/////////////////////////////////////////////////

#if FNM_PULSED_WAVE
// TODO(JEM): Remove static variable
template <class T>
bool Aperture<T>::normalize = true;
#endif

#ifdef FNM_CLOSURE_FUNCTIONS
template <class T>
ApertureData<T>* Aperture<T>::DataGet() {
  return this->m_pData;
}
#endif

template <class T>
sysparm_t<T>* Aperture<T>::DefaultSysParmGet() {
  ////////////////////////////////////////////////
  // Thread-safe using Acquire-Release semantics
  ////////////////////////////////////////////////
  sysparm_t<T>* pGlobalVar =
    sps::globalstruct<T, fnm::sysparm_t>::pVar.load(std::memory_order_acquire);
  if (!pGlobalVar) {
    std::lock_guard<std::mutex> myLock(myMutex);
    pGlobalVar =
      sps::globalstruct<T, fnm::sysparm_t>::
      pVar.load(std::memory_order_relaxed);
    if (!pGlobalVar) {
      pGlobalVar = new sysparm_t<T>();
      pGlobalVar->nThreads = getncpus();
      sps::globalstruct<T, fnm::sysparm_t>::
      pVar.store(pGlobalVar, std::memory_order_release);
    }
  }
  return pGlobalVar;
}

////////////////////////////////
// Implementation of interface
////////////////////////////////
template <class T>
Aperture<T>::Aperture() : m_sysparm(NULL) {
  // Assign pointer to global sysparm
  this->m_sysparm = Aperture<T>::DefaultSysParmGet();

  // Assign point to singleton
  this->m_pTest = fnm::bla<T>::InstanceGet();

  // Aligned data segment
  m_pData =
    static_cast<ApertureData<T>*>(
      SPS_MM_MALLOC(sizeof(ApertureData<T>), 4*sizeof(T)));
  new (this->m_pData) ApertureData<T>();

  // Progress bar
  m_pbar = nullptr;

  // Data needed for time-domain pulsed waves
  this->m_pData->m_pulses->m_fs = m_sysparm->fs;
}

template <class T>
Aperture<T>::Aperture(
  const size_t nElements, const T width, const T kerf,
  const T height) : Aperture<T>() {
  const size_t nSubH = 1;
  const size_t nSubW = 1;
  const size_t nSubElements = nSubH * nSubW;
  const T focus = T(0.0);

  if (nElements * nSubElements > 0) {
    auto elements = element_array<T>();

    FocusedLinearArray(nElements, nSubH, nSubW,
                       width, kerf, height,
                       focus, 0, std::move(elements));

    m_pData->ElementsSet(std::forward<element_array<T> >(elements),
                         nElements,
                         nSubElements);
  }
}

template <class T>
Aperture<T>::~Aperture() {
  if (m_pData) {
    m_pData->~ApertureData();
    _mm_free(m_pData);
    m_pData = nullptr;
  }

  // Free local system parameters if they are set
  if (m_sysparm != Aperture<T>::DefaultSysParmGet()) {
    delete m_sysparm;
  }
  if (m_pTest != bla<T>::InstanceGet()) {
    delete m_pTest;
  }
}

template <class T>
void Aperture<T>::Rotate(const T iEuler[3], const T iFocus[3]) {
  sps::point_t<T> reference;
  sps::euler_t<T> euler;
  euler.alpha = iEuler[0];
  euler.beta  = iEuler[1];
  euler.gamma = iEuler[2];

  memcpy(&reference[0], iFocus, 3*sizeof(T));

  this->m_pData->Rotate(reference, euler);
}

template <class T>
void Aperture<T>::SubToElements() {
  size_t nElements = 0, nSubElements = 0, nPositions = 0, nParams = 0;
  T* pos = NULL;
  T* delays = NULL;

  // Allocates
  this->DelaysGet(&delays, &nPositions);

  // Allocates
  this->MultiElementsGet(1, &pos, &nElements, &nSubElements, &nParams);

  assert(nParams == 8);

  debug_print("nElements: %zu, nSubElements: %zu\n", nElements, nSubElements);

  size_t nNewElements = nElements * nSubElements;

  this->SubElementsSet(pos, nNewElements, 1, 8);

  T* newDelays;

  if (NULL ==
      (newDelays = static_cast<T*>(SPS_MALLOC(nNewElements*sizeof(T))))) {
    goto SubToElementsError;
  }

  for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
    for (size_t jElement = 0 ; jElement < nSubElements ; jElement++) {
      newDelays[iElement*nSubElements + jElement] = delays[iElement];
    }
  }

  this->DelaysSet(newDelays, nNewElements);

SubToElementsError:
  // Freeing NULL is defined as a no-op, so no need to check pointers
  if (pos) {
    free(pos);
  }
  if (delays) {
    free(delays);
  }
  if (newDelays) {
    free(newDelays);
  }
}

template <class T>
int
Aperture<T>::ArrayCreate(fnm::Aperture<T> **obj) {
  *obj = new Aperture<T>();
  return 0;
}

template <class T>
int
Aperture<T>::MatrixArrayCreate(
  Aperture<T> **obj,
  const size_t nRows,
  const size_t nCols,
  const T rowWidth,
  const T rowKerf,
  const T colWidth,
  const T colKerf) {
  const size_t nSubElements = 1;
  const size_t nElements = nRows*nCols;

  *obj = NULL;

  int retval = -1;

  if (nElements * nSubElements > 0) {
    *obj = new Aperture<T>();

    auto elements = element_array<T>();

    MatrixArray(nRows, nCols, rowWidth, rowKerf,
                colWidth, colKerf, std::move(elements));

    (*obj)->m_pData->ElementsSet(std::forward<element_array<T> >(elements),
                                 nElements,
                                 nSubElements);
    retval = 0;
  }

  return retval;
}

template <class T>
int
Aperture<T>::FocusedLinearArrayCreate(
  fnm::Aperture<T> **obj, const size_t nElements,
  const T width, const T kerf, const T height,
  const size_t nSubH, const T focus) {
  const size_t nSubW = 1;
  const size_t nSubElements = nSubH * nSubW;

  // Null object in case of zero elements
  *obj = NULL;

  int retval = -1;

  if (nElements * nSubElements > 0) {
    *obj = new Aperture<T>();

    auto elements = element_array<T>();

    FocusedLinearArray(nElements, nSubH, nSubW,
                       width, kerf, height,
                       focus,
                       0, /* 0: Outer placement, 1: Inner placement */
                       std::move(elements));

    (*obj)->m_pData->ElementsSet(std::forward<element_array<T> >(elements),
                                 nElements,
                                 nSubElements);
    retval = 0;
  }

  return retval;
}

template <class T>
int Aperture<T>::
FocusedConvexArrayCreate(
  fnm::Aperture<T> **obj,
  const size_t nElements, const T width, const T kerf, const T height,
  const T radius, const size_t nSubH, const T focus) {
  // No sub-elements in azimuth (width)
  const size_t nSubW = 1;
  const size_t nSubElements = nSubH * nSubW;

  // Null object in case of zero elements
  *obj = NULL;

  int retval = -1;

  if (nElements * nSubElements > 0) {
    retval = 0;
    *obj = new Aperture<T>();

    auto elements = element_array<T>();

    FocusedConvexArray(nElements, nSubH, nSubW,
                       width, kerf, height, radius,
                       focus,
                       0, /* 0: Outer placement, 1: Inner placement */
                       std::move(elements));  // elements is now an r-value

    // Why not use std::move here. If elements were an argument to this
    // function, forward would be correct
    (*obj)->m_pData->ElementsSet(std::forward<element_array<T> >(elements),
                                 nElements,
                                 nSubElements);
  }
  return retval;
}

template <class T>
void Aperture<T>::InPlaceOP(T* ioData, size_t nIOdata) {
  SPS_UNREFERENCED_PARAMETER(nIOdata);
  printf("ioData[0]: %f\n", ioData[0]);
  ioData[0] = T(1.0);
}


// One could read back sysparm, modify and write back
template <class T>
void Aperture<T>::AlphaSet(const T& value) {
  this->m_sysparm->att = value;
}

// Returning by argument also
template <class T>
const T& Aperture<T>::AlphaGet() const {
  return this->m_sysparm->att;
}

template <class T>
void Aperture<T>::BetaSet(const T& value) {
  this->m_sysparm->beta = value;
}

template <class T>
const T& Aperture<T>::BetaGet() const {
  return this->m_sysparm->beta;
}

template <class T>
void Aperture<T>::AttenuationEnabledSet(const bool& iEnabled) {
  this->m_sysparm->use_att = iEnabled;
}

template <class T>
const bool& Aperture<T>::AttenuationEnabledGet() const {
  return this->m_sysparm->use_att;
}

template <class T>
const size_t& Aperture<T>::NThreadsGet() const {
  return this->m_sysparm->nThreads;
}

template <class T>
int Aperture<T>::NThreadsSet(const size_t &nThreads) {
  assert(nThreads <= N_MAX_THREADS);
  size_t _nthreads = std::max<size_t>(1U, nThreads);
  int retval = -1;
  if (_nthreads <= N_MAX_THREADS) {
    this->m_sysparm->nThreads = _nthreads;
    retval = 0;
  }
  return retval;
}

template <class T>
void Aperture<T>::PositionsGet(
  T **out, size_t* nElements, size_t* nParams) const {
  const size_t _nElements = m_pData->m_nelements;
  const size_t arrSize    = _nElements*3*sizeof(T);
  const auto& positions   = m_pData->m_pos;

  // The function allocates
  T* arr = static_cast<T*>(SPS_MALLOC(arrSize));
  if (arr) {
    memset(arr, 0, arrSize);
    for (size_t i = 0 ; i < _nElements ; i++) {
      arr[i*3 + 0] = positions[i][0];
      arr[i*3 + 1] = positions[i][1];
      arr[i*3 + 2] = positions[i][2];
    }
  }

  // Assign outputs
  *nElements = m_pData->m_nelements;
  *nParams   = 3;
  *out = arr;
}

template <class T>
void Aperture<T>::
ExtentGet(T** coordinates, size_t* nDim, size_t* nLimits) const {
  const size_t arrSize = 6;

  *coordinates = nullptr;
  *nDim = 0;
  *nLimits = 0;

  T* arr = static_cast<T*>(SPS_MALLOC(arrSize*sizeof(T)));

  if (arr) {
    memset(arr, 0, arrSize);

    sps::bbox_t<T> bbox;
    this->m_pData->ExtentGet(&bbox);

    typedef T extent_t[3][2];
    extent_t* extent = reinterpret_cast<extent_t*>(arr);

    for (size_t iDim = 0 ; iDim < 3 ; iDim++) {
      (*extent)[iDim][0] = bbox.min[iDim];
      (*extent)[iDim][1] = bbox.max[iDim];
    }

    *coordinates = arr;
    *nDim        = 3;
    *nLimits     = 2;
  }
}

template <class T>
T Aperture<T>::AreaGet() const {
  return m_pData->AreaGet();
}

template <class T>
int Aperture<T>::
PositionsSet(const T* pos, const size_t nPositions, const size_t nDim) {
  const size_t _nElements          = m_pData->m_nelements;
  int retval                       = -1;
  sps::point_t<T>* positions       = m_pData->m_pos.get();
  if ((_nElements == nPositions) && (3 == nDim)) {
    // Consider creating and moving a unique_ptr here
    for (size_t i = 0 ; i < _nElements ; i++) {
      positions[i][0] = pos[i*3 + 0];
      positions[i][1] = pos[i*3 + 1];
      positions[i][2] = pos[i*3 + 2];
    }
    retval = 0;
  }
  return retval;
}

template <class T>
void Aperture<T>::FocusSet(const T iFocus[3]) {
  bool changed = false;
  if ((iFocus[0] != m_pData->m_focus[0]) ||
      (iFocus[1] != m_pData->m_focus[1]) ||
      (iFocus[2] != m_pData->m_focus[2])) {
    changed = true;
    m_pData->m_focus_valid = FocusingType::FocusingTypeCount;
  }
  memcpy(&m_pData->m_focus[0], &iFocus[0], 3*sizeof(T));

  if (changed) {
    if (m_pData->m_apodization_type > fnm::ApodizationType::ApodizationTypeNonParametric) {
      // Update apodization using f-number
    }
  }
}

template <class T>
void Aperture<T>::Focus2Get(T oFocus[3]) const {
  memcpy(oFocus, &m_pData->m_focus2[0], 3*sizeof(T));
}

template <class T>
void Aperture<T>::Focus2Set(const T iFocus[3]) {
  bool changed = false;
  if ((iFocus[0] != m_pData->m_focus2[0]) ||
      (iFocus[1] != m_pData->m_focus2[1]) ||
      (iFocus[2] != m_pData->m_focus2[2])) {
    changed = true;
    m_pData->m_focus_valid = FocusingType::FocusingTypeCount;
  }
  memcpy(&m_pData->m_focus2[0], &iFocus[0], 3*sizeof(T));

  if (changed) {
    if (m_pData->m_apodization_type > fnm::ApodizationType::ApodizationTypeNonParametric) {
      // Update apodization using f-number
    }
  }
}

template <class T>
const T& Aperture<T>::XmtFNumberGet() const {
  return m_pData->m_fnumber;
}

// TODO: Don't do this!!!
template <class T>
void Aperture<T>::XmtFNumberSet(const T& fnumber) {
  m_pData->m_fnumber = fnumber;
  if (m_pData->m_apodization_type > ApodizationType::ApodizationTypeNonParametric) {
    // Parametric - not good approach
    sps::point_t<T> direction;
    sps::point_t<T> r2p = this->m_pData->m_focus - this->m_pData->m_center_focus;
    T depth = sps::norm(r2p);

    if (depth > std::numeric_limits<T>::epsilon()) {
      // Avoid division by zero
      direction[0] = r2p[0] / depth;
      direction[1] = r2p[1] / depth;
      direction[2] = r2p[2] / depth;

      // Simple strategy (distance to element - not good)
      this->m_pData->ApodizationSet(
        direction, depth,
        static_cast<ApodizationType>(this->m_pData->m_apodization_type));
    }
  }
}

template <class T>
void Aperture<T>::FocusGet(T oFocus[3]) const {
  memcpy(oFocus, &m_pData->m_focus[0], 3*sizeof(T));
}

template <class T>
void Aperture<T>::CenterFocusSet(const T iFocus[3]) {
  if ((iFocus[0] != m_pData->m_center_focus[0]) ||
      (iFocus[1] != m_pData->m_center_focus[1]) ||
      (iFocus[2] != m_pData->m_center_focus[2])) {
    m_pData->m_focus_valid = FocusingType::FocusingTypeCount;
  }
  memcpy(&m_pData->m_center_focus[0], &iFocus[0], 3*sizeof(T));
}

template <class T>
void Aperture<T>::CenterFocusGet(T oFocus[3]) const {
  memcpy(oFocus, &m_pData->m_center_focus[0], 3*sizeof(T));
}

template <class T>
int Aperture<T>::ApodizationSet(const T* data, size_t nData) {
  const size_t nElements = m_pData->m_nelements;
  int retval = -1;

#if 0
  // Issue here with forwarding (perhaps)
  if (nData == nElements) {
    sps::unique_aligned_array<T> apodization = sps::unique_aligned_array_create<T>(nData);

    for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
      apodization[iElement] = data[iElement];
    }
    // TODO(JMH): Make this work
    retval = m_pData->ApodizationSet(std::forward<sps::unique_aligned_array<T> >(apodization), nData);
  }
#else
  if (nData == nElements) {
    for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
      m_pData->m_sensitivities.get()[iElement] = data[iElement];
    }
    retval = 0;
  }
#endif
  return retval;
}

template <class T>
void Aperture<T>::ApodizationGet(T** data, size_t* nData) const {
  const size_t nElements = m_pData->m_nelements;

  *data = static_cast<T*>(SPS_MALLOC(nElements*sizeof(T)));
  // Perfectly valid to copy zero elements, so no need to do this
  if (nElements > 0 && *data) {
    memcpy(*data, m_pData->m_sensitivities.get(), nElements*sizeof(T));
    *nData = nElements;
  } else {
    *nData = 0;
  }
}

template <class T>
void Aperture<T>::PhasesGet(T** odata, size_t* nData) const {
  const size_t nElements    = m_pData->m_nelements;
  const size_t nSubElements = m_pData->m_nsubelements;

  *odata = static_cast<T*>(SPS_MALLOC(nElements*sizeof(T)));
  if (nElements > 0) {
    for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
      (*odata)[iElement] = m_pData->m_phases[iElement*nSubElements];
    }
    *nData = nElements;
  } else {
    *nData = 0;
  }
}

template <class T>
void Aperture<T>::SubPhasesGet(T** odata, size_t* nData) const {
  const size_t nElements = m_pData->m_nelements;
  const size_t nSubElements = m_pData->m_nsubelements;

  *odata = static_cast<T*>(SPS_MALLOC(nElements*nSubElements*sizeof(T)));
  if (nElements > 0) {
    memcpy(*odata, m_pData->m_phases.get(), nElements*nSubElements*sizeof(T));
    *nData = nElements*nSubElements;
  } else {
    *nData = 0;
  }
}

template <class T>
void Aperture<T>::DelaysGet(T** odata, size_t* nData) const {
  const T* pDelays = nullptr;
  size_t nElements = 0;
  m_pData->DelaysRefGet(&nElements, pDelays);
  *odata = static_cast<T*>(SPS_MALLOC(nElements*sizeof(T)));

  if (nElements > 0) {
    memcpy(*odata, pDelays, nElements*sizeof(T));
    *nData = nElements;
  } else {
    *nData = 0;
  }
}

template <class T>
int Aperture<T>::DelaysSet(const T* data, const size_t nData) {
  const size_t nElements = m_pData->m_nelements;
  int retval = -1;
  if (nData == m_pData->m_nelements) {
    memcpy(m_pData->m_delays.get(), data, sizeof(T)*nElements);
    // Set focusing type to using delays
    m_pData->m_focus_type = FocusingType::Delays;
    retval = 0;
  }
  return retval;
}

template <class T>
int Aperture<T>::FocusingTypeGet() const {
  return m_pData->m_focus_type;
}

template <class T>
void Aperture<T>::FocusingTypeSet(const int iFocusingType) {
  m_pData->m_focus_type = iFocusingType;
}

template <class T>
int Aperture<T>::ApodizationTypeGet() const {
  return m_pData->m_apodization_type;
}

template <class T>
void Aperture<T>::ApodizationTypeSet(const int iApodizationType) {
  m_pData->m_apodization_type = iApodizationType;
}

template <class T>
const T& Aperture<T>::WGet() const {
  return m_sysparm->w;
}

template <class T>
void Aperture<T>::WSet(const T& w) {
  m_sysparm->w = w;
}

template <class T>
const T& Aperture<T>::F0Get() const {
  return m_pData->m_f0;
}

// TODO: Distinguish between f0 and xdcf0 (or fc) (introduce excitation object)
template <class T>
void Aperture<T>::F0Set(const T& f0) {
  // Make static variable
  const T eps = std::numeric_limits<T>::epsilon();
  assert(f0 > eps);
  if (f0 > eps) {
    // Invalidate focus
    m_pData->m_focus_valid = ApodizationType::ApodizationTypeCount;
    m_pData->m_f0 = f0;
    m_sysparm->f0 = f0;
  }
}

template <class T>
const sysparm_t<T> Aperture<T>::SysParmGet() const {
  // A copy is returned, so a new structure is created, when setting this again
  return (*this->m_sysparm);
}

template <class T>
void Aperture<T>::SysParmSet(const sysparm_t<T> *arg) {
  if (this->m_sysparm != Aperture<T>::DefaultSysParmGet()) {
    delete this->m_sysparm;
  }
  this->m_sysparm = new sysparm_t<T>();
  (*this->m_sysparm) = *arg;

  if (m_pTest != bla<T>::InstanceGet()) {
    delete this->m_pTest;
  }
  this->m_pTest = new bla<T>();
  (*this->m_pTest) = bla<T>();
}

template <class T>
const T& Aperture<T>::DensityGet() const {
  return m_sysparm->rho;
}

template <class T>
void Aperture<T>::DensitySet(const T& rho) {
  m_sysparm->rho = rho;
}

template <class T>
int Aperture<T>::ExcitationTypeGet() const {
  return this->m_pData->m_pulses->m_excitationType;
}

template <class T>
void Aperture<T>::ExcitationTypeSet(const int iExcitationType) {
  this->m_pData->m_pulses->m_excitationType = (sofus::ExcitationType) iExcitationType;
}

template <class T>
void fnm::Aperture<T>::BandWidthSet(const T& value) {
  this->m_pData->m_pulses->BandWidthSet(value);
}

template <class T>
const T& fnm::Aperture<T>::BandWidthGet() const {
  return this->m_pData->m_pulses->m_impulseBandwidth;
}

template <class T>
const T& Aperture<T>::FCGet() const {
  return this->m_pData->m_pulses->m_f0;
}

template <class T>
void Aperture<T>::FCSet(const T& f0) {
  const T eps = std::numeric_limits<T>::epsilon();

  // Dependency graph: (f0, fs, bw, parametric) -> new pulse

  assert(f0 > eps);
  if (f0 > eps) {
    this->m_pData->m_pulses->m_f0 = f0;
    // Update parametric impulse
    if (this->m_pData->m_pulses->m_impulseType == sofus::ImpulseType::ImpulseTypeGaussian) {
      // FC is f0 of pulse
      this->m_pData->m_pulses->GaussPulseSet(this->FsGet(), this->FCGet(), this->BandWidthGet());
    } else {
      this->m_pData->m_pulses->DeltaPulseSet(T(0.0));
    }
  }
}

template <class T>
const T& Aperture<T>::FsGet() const {
  return m_sysparm->fs;
}

template <class T>
int Aperture<T>::FsSet(const T& value) {
  int retval = -1;
  if (value > T(1.0)) {
    m_sysparm->fs = value;
    // Sysparm and pulse must be in sync
    this->m_pData->m_pulses->FsSet(value);
    retval = 0;
  }
  return retval;
}

#if FNM_PULSED_WAVE
template <class T>
const bool& Aperture<T>::NormalizeGet() const {
  return Aperture<T>::normalize;
}

template <class T>
void Aperture<T>::NormalizeSet(const bool& value) {
  Aperture<T>::normalize = value;
}

template <class T>
void Aperture<T>::
FocusLinesCreate(size_t nLines, sofus::FocusLineList<T>** obj) const {
  const size_t nElements = m_pData->m_nelements;
  *obj = new sofus::FocusLineList<T>(nLines, nElements);
}


template <class T>
int Aperture<T>::ImpulseTypeGet() const {
  return this->m_pData->m_pulses->m_impulseType;
}

template <class T>
void Aperture<T>::ImpulseTypeSet(const int iImpulseType) {
  sofus::AperturePulses<T>* pPulse = this->PulsesGet();
  pPulse->m_impulseType = (sofus::ImpulseType) iImpulseType;
  if (this->m_pData->m_pulses->m_impulseType == sofus::ImpulseType::ImpulseTypeGaussian) {
    this->m_pData->m_pulses->GaussPulseSet(this->FsGet(), this->FCGet(), this->BandWidthGet());
  }
}


template <class T>
T Aperture<T>::CalcMatchedFilter(
  const Aperture<T>* other,
  const T iFocus[3],
  T** odata, size_t* nSignals, size_t* nSamples,
  int** sizeAndOffsets, size_t* nFilters, size_t* nTwo,
  T** data, size_t* nData) {

  if (m_pData->m_focus_type == FocusingType::Rayleigh) {
    fprintf(stderr, "Rayleigh focusing type is not supported for pulsed wave calculations\n");
    fprintf(stderr, "Focusing type is reset to Pythagorean\n");
    m_pData->m_focus_type = FocusingType::Pythagorean;
  }

  this->FocusUpdate();

  sofus::sysparm_t<T> sysparm;
  sysparm.c         = m_sysparm->c;
  sysparm.fs = m_sysparm->fs;
  sysparm.normalize = Aperture<T>::normalize;
  sysparm.att       = m_sysparm->att;
  sysparm.beta      = m_sysparm->beta;
  sysparm.use_att   = m_sysparm->use_att;
  sysparm.f0        = m_pData->m_f0;
  sysparm.timeDomainCalcType = m_sysparm->timeDomainCalcType;
  sysparm.timeDomainIntOrder = m_sysparm->timeDomainIntOrder;
  sysparm.soft_baffle = false;

  sps::point_t<T> point;
  point[0] = iFocus[0];
  point[1] = iFocus[1];
  point[2] = iFocus[2];

  T tStart = sofus::CalcMatchedFilter(sysparm,
                                      m_pData, other->m_pData,
                                      point,
                                      odata, nSignals, nSamples,
                                      sizeAndOffsets, nFilters, nTwo,
                                      data, nData);
  return tStart;
}

template <class T>
T Aperture<T>::CalcSmfApply(
  const Aperture<T>* other,
  const T iFocus[3],
  const T tStart,
  const T* iData, const size_t nChannels, const size_t nSamples,
  T** data, size_t* nData) {
  if (m_pData->m_focus_type == FocusingType::Rayleigh) {
    fprintf(stderr, "Rayleigh focusing type is not supported for pulsed wave calculations\n");
    fprintf(stderr, "Focusing type is reset to Pythagorean\n");
    m_pData->m_focus_type = FocusingType::Pythagorean;
  }

  this->FocusUpdate();

  sofus::sysparm_t<T> sysparm;
  sysparm.c         = m_sysparm->c;
  sysparm.fs = m_sysparm->fs;
  sysparm.normalize = Aperture<T>::normalize;
  sysparm.att       = m_sysparm->att;
  sysparm.beta      = m_sysparm->beta;
  sysparm.use_att   = m_sysparm->use_att;
  sysparm.f0        = m_pData->m_f0;
  sysparm.timeDomainCalcType = m_sysparm->timeDomainCalcType;
  sysparm.timeDomainIntOrder = m_sysparm->timeDomainIntOrder;
  sysparm.soft_baffle = false;

  sps::point_t<T> point;
  point[0] = iFocus[0];
  point[1] = iFocus[1];
  point[2] = iFocus[2];

  return sofus::CalcSmfApply(sysparm,
                             m_pData,
                             other->m_pData,
                             point,
                             tStart,
                             iData, nChannels, nSamples, data, nData);
}

template <class T>
T Aperture<T>::CalcScat(
  const Aperture<T>* other,
  const T* pos, const size_t nPositions, const size_t nDim,
  const T* data, const size_t nData,
  T** odata, size_t* nSignals, size_t* nSamples) {
  assert(nDim == 3);

  if ((nDim != 3) || (!pos) || nPositions == 0) {
    *odata = NULL;
    *nSignals = 0;
    *nSamples = 0;
    return T(0.0);
  }

  // TODO: Find better way - agree on content of sysparm or split in two
  sofus::sysparm_t<T> sysparm;
  sysparm.c         = m_sysparm->c;
  sysparm.fs        = m_sysparm->fs;
  sysparm.normalize = Aperture<T>::normalize;

  sysparm.att     = m_sysparm->att;
  sysparm.beta    = m_sysparm->beta;
  sysparm.use_att = m_sysparm->use_att;
  sysparm.f0      = m_pData->m_f0;
  sysparm.timeDomainCalcType = m_sysparm->timeDomainCalcType;
  sysparm.timeDomainIntOrder = m_sysparm->timeDomainIntOrder;
  sysparm.soft_baffle = false;

  return sofus::CalcScat(&sysparm,
                         this->m_pData,
                         other->m_pData,
                         pos, nPositions, nDim,
                         data, nData,
                         odata, nSignals, nSamples);
}

#if 1 // __GNUG__ || defined(_MSC_VER) && (_MSC_VER >= 1900)
template <class T>
T Aperture<T>::CalcScatThreaded(
  const Aperture<T>* other,
  const sofus::FocusLineList<T>* pXmtLines,
  const T* pos, const size_t nPositions, const size_t nDim,
  const T* data, const size_t nData,
  T** odata, size_t* nSignals, size_t* nSamples) {
  assert(nDim == 3);

  if ( (nDim != 3) || (!pos) || nPositions == 0 ) {
    *odata = NULL;
    *nSignals = 0;
    *nSamples = 0;
    return T(0.0);
  }

  // TODO: Find better way - agree on content of sysparm or split in two
  sofus::sysparm_t<T> sysparm;
  sysparm.c         = m_sysparm->c;
  sysparm.fs = m_sysparm->fs;
  sysparm.normalize = Aperture<T>::normalize;

  sysparm.att     = m_sysparm->att;
  sysparm.beta    = m_sysparm->beta;
  sysparm.use_att = m_sysparm->use_att;
  sysparm.f0      = m_pData->m_f0;
  sysparm.timeDomainCalcType = m_sysparm->timeDomainCalcType;
  sysparm.timeDomainIntOrder = m_sysparm->timeDomainIntOrder;
  sysparm.soft_baffle = false;
  sysparm.nThreads = m_sysparm->nThreads;

  return sofus::CalcScatThreaded(&sysparm,
                                 pXmtLines,
                                 this->m_pData,
                                 other->m_pData,
                                 this->m_pData->m_pulses,
                                 other->m_pData->m_pulses,
                                 pos, nPositions, nDim,
                                 data, nData,
                                 odata, nSignals, nSamples);
}
#endif

template <class T>
T Aperture<T>::CalcPwEcho(
  const Aperture<T>* other,
  const T* pos, const size_t nPositions, const size_t nDim,
  const T* data, const size_t nData,
  T** odata, size_t* nSignals, size_t* nSamples) {
  assert(nDim == 3);

  if ((nDim != 3) || (!pos) || nPositions == 0) {
    *odata = NULL;
    *nSignals = 0;
    *nSamples = 0;
    return T(0.0);
  }

  if (m_pData->m_focus_type == FocusingType::Rayleigh) {
    fprintf(stderr, "Rayleigh focusing type is not supported for pulsed wave calculations\n");
    fprintf(stderr, "Focusing type is reset to Pythagorean\n");
    m_pData->m_focus_type = FocusingType::Pythagorean;
  }
  this->FocusUpdate();

  // TODO: Find better way - agree on content of sysparm or split in two
  sofus::sysparm_t<T> sysparm;
  sysparm.c         = m_sysparm->c;
  sysparm.fs = m_sysparm->fs;
  sysparm.normalize = Aperture<T>::normalize;

  sysparm.att     = m_sysparm->att;
  sysparm.beta    = m_sysparm->beta;
  sysparm.use_att = m_sysparm->use_att;
  sysparm.f0      = m_pData->m_f0;
  sysparm.timeDomainCalcType = m_sysparm->timeDomainCalcType;
  sysparm.timeDomainIntOrder = m_sysparm->timeDomainIntOrder;
  sysparm.soft_baffle = false;
  T retval = sofus::CalcPwEcho(sysparm,
                               m_pData,
                               other->m_pData,
                               this->m_pData->m_pulses,
                               other->m_pData->m_pulses,
                               pos, nPositions, nDim,
                               data, nData,
                               odata, nSignals, nSamples, this->m_pbar);
  return retval;
}



template <class T>
void Aperture<T>::ImpulseGet(T** data, size_t* nData) const {
  *nData = this->m_pData->m_pulses->impulse.ndata;
  *data  = this->m_pData->m_pulses->impulse.m_data.get();
}

// TODO: Consider adding a method, which takes data and type (use default)
template <class T>
int Aperture<T>::ImpulseSet(const T* data,
                            const size_t nData) {
  this->m_pData->m_pulses->impulse.ndata = nData;
  this->m_pData->m_pulses->impulse.offset = 0;
  this->m_pData->m_pulses->impulse.m_data = NULL;
  int retval = -1;

  if (nData > 0) {
    this->m_pData->m_pulses->impulse.m_data  =
    std::shared_ptr<T>(static_cast<T*>(SPS_MM_MALLOC(sizeof(T)*nData,16)), [=](T *p) {
      _mm_free(p);
    });
    memcpy((void*)this->m_pData->m_pulses->impulse.m_data.get(), data, nData*sizeof(T));
    this->m_pData->m_pulses->m_impulseType = sofus::ImpulseType::ImpulseTypeNonParametric;
    retval = 0;
  }
  return retval;
}

template <class T>
void Aperture<T>::ExcitationGet(T** data, size_t* nData) const {
  // A view is returned. We cannot delete this in the wrapper
  *nData = this->m_pData->m_pulses->excitation.ndata;
  *data  = this->m_pData->m_pulses->excitation.m_data.get();
}

template <class T>
int Aperture<T>::ExcitationSet(const T* data,
                               const size_t nData) {
  int retval = -1;
  this->m_pData->m_pulses->excitation.ndata  = nData;
  this->m_pData->m_pulses->excitation.offset = 0;
  this->m_pData->m_pulses->excitation.m_data = NULL;

  if (nData > 0) {
    this->m_pData->m_pulses->excitation.m_data  =
      std::shared_ptr<T>( static_cast<T*>(SPS_MM_MALLOC(sizeof(T)*nData, 16)),
    [=](T *p) {
      _mm_free(p);
    });
    memcpy(this->m_pData->m_pulses->excitation.m_data.get(), data, nData*sizeof(T));
    retval = 0;
  } else {
    this->m_pData->m_pulses->excitation.m_data = nullptr;
  }
  return retval;
}
#endif

template <class T>
const T& Aperture<T>::CGet()  const {
  return this->m_sysparm->c;
}

template <class T>
void Aperture<T>::CSet(const T& c) {
  this->m_sysparm->c = c;
  // Invalidate focus
  m_pData->m_focus_valid = FocusingType::FocusingTypeCount;
}

template <class T>
const size_t& Aperture<T>::NElementsGet() const {
  return m_pData->m_nelements;
}

template <class T>
const size_t& Aperture<T>::NSubElementsGet() const {
  return m_pData->m_nsubelements;
}

template <class T>
const size_t& Aperture<T>::NDivWGet() const {
  return this->m_sysparm->nDivW;
}

template <class T>
int Aperture<T>::NDivWSet(const size_t& nDivW) {
  int retval = -1;
  size_t _nDivW = nDivW;
  if (nDivW > 1) {
    if (nDivW % 2 == 1) {
      printf("\nNDivW rounded up to be even\n");
      _nDivW = _nDivW + 1;
    }
    this->m_sysparm->nDivW = _nDivW;
    retval = 0;
  }
  return retval;
}

template <class T>
const size_t& Aperture<T>::NDivHGet() const {
  return this->m_sysparm->nDivH;
}

template <class T>
int Aperture<T>::NDivHSet(const size_t& nDivH) {
  int retval = -1;
  size_t _nDivH = nDivH;

  if (nDivH > 1) {
    if (nDivH % 2 == 1) {
      printf("\nNDivH rounded up to be even\n");
      _nDivH = _nDivH + 1;
    }
    this->m_sysparm->nDivH = _nDivH;
    retval = 0;
  }
  return retval;
}

template <class T>
void Aperture<T>::RectanglesGet(T** out, size_t* nElements,
                                size_t* nSubElements, size_t* nParams) const {
  const size_t nCornerCoordinates = 12;

  const size_t _nElements    = m_pData->m_nelements;
  const size_t _nSubElements = m_pData->m_nsubelements;
  const size_t arrSize       = _nElements * _nSubElements * nCornerCoordinates * sizeof(T);

  // The function allocates
  T* arr = static_cast<T*>(SPS_MALLOC(arrSize));
  if (arr) {
    memset(arr, 0, arrSize);

    const auto& rectangles = m_pData->m_rectangles;

    for (size_t iElement = 0 ; iElement < _nElements ; iElement++) {
      for (size_t iSubElement = 0 ; iSubElement < _nSubElements ; iSubElement++) {
        for (size_t iCorner = 0 ; iCorner < Aperture<T>::nVerticesPerElement ; iCorner++) {
          // We need the order [0,1,3,2]
          size_t iiCorner = iCorner;
          if (iCorner == 2) {
            iiCorner = 3;
          } else if (iCorner == 3) {
            iiCorner = 2;
          }
          for (size_t iXYZ = 0 ; iXYZ < 3 ; iXYZ++) {
            arr[iElement * _nSubElements * Aperture<T>::nVerticesPerElement * 3 +
                         iSubElement * Aperture<T>::nVerticesPerElement * 3 + iCorner*3 + iXYZ] =
                  rectangles[iElement*_nSubElements +
                             iSubElement][iiCorner][iXYZ];
          }
        }
      }
    }

    *nElements    = m_pData->m_nelements;
    *nSubElements = m_pData->m_nsubelements;
    *nParams      = nCornerCoordinates;
    *out          = arr;
  }
}

template <class T>
int Aperture<T>::ElementsRefGet(size_t* nElements, size_t* nSubElements,
                                const sps::element_rect_t<T>**& elements) const {
  return m_pData->ElementsRefGet(nElements, nSubElements, elements);
}

template <class T>
int Aperture<T>::ApodizationsRefGet(size_t* nElements, const T*& apodizations) const {
  return m_pData->ApodizationsRefGet(nElements, apodizations);
}

template <class T>
void Aperture<T>::ElementsGet(T** out, size_t* nElements,
                              size_t* nParams) const {
  size_t nSubElements = 0;
  // The first sub-element per element is returned (it may not be central)
  int err = this->MultiElementsGet(0, out, nElements, &nSubElements, nParams);
  SPS_UNREFERENCED_PARAMETERS(err);
}

template <class T>
void Aperture<T>::SubElementsGet(T** out, size_t* nElements,
                                 size_t* nSubElements, size_t* nParams) const {
  int err = this->MultiElementsGet(1, out, nElements, nSubElements, nParams);
  SPS_UNREFERENCED_PARAMETERS(err);
}

template <class T>
int Aperture<T>::MultiElementsGet(int sub,
                                  T** out, size_t* nElements,
                                  size_t* nSubElements, size_t* nParams) const {
  const size_t nElePosParams = Aperture<T>::nElementPosParameters;

  const size_t _nElements             = m_pData->m_nelements;
  const size_t _nSubElements          = m_pData->m_nsubelements;

  // Size of output
  size_t nOutputElements              = _nElements;
  size_t nOutputSubElementsPerElement = _nSubElements;

  // Reduce size if no sub-elements are wanted
  if (sub == 0) {
    nOutputSubElementsPerElement = 1;
  }

  // Size of output array
  const size_t arrSize =
    nOutputElements * nOutputSubElementsPerElement * nElePosParams * sizeof(T);

  // Error code
  int retval = -1;

  T* arr = static_cast<T*>(SPS_MALLOC(arrSize));

  if (arr) {
    memset(arr, 0, arrSize);

    size_t __nElements, __nSubElements;
    const sps::element_rect_t<T>** elements = NULL;

    retval = m_pData->ElementsRefGet(&__nElements, &__nSubElements, elements);

    for (size_t iElement = 0 ; iElement < nOutputElements ; iElement++) {
      for (size_t jElement = 0 ; jElement < nOutputSubElementsPerElement ; jElement++) {
        arr[iElement * nOutputSubElementsPerElement * nElePosParams +
            jElement * nElePosParams + 0] = elements[iElement][jElement].hw;

        arr[iElement * nOutputSubElementsPerElement * nElePosParams +
            jElement * nElePosParams + 1] = elements[iElement][jElement].hh;

        memcpy(&arr[iElement * nOutputSubElementsPerElement * nElePosParams + jElement*nElePosParams + 2],
               &elements[iElement][jElement].center[0],
               sizeof(sps::point_t<T>));
        arr[iElement * nOutputSubElementsPerElement * nElePosParams + jElement*nElePosParams+5] =
          elements[iElement][jElement].euler.alpha;
        arr[iElement * nOutputSubElementsPerElement * nElePosParams + jElement*nElePosParams+6] =
          elements[iElement][jElement].euler.beta;
        arr[iElement * nOutputSubElementsPerElement * nElePosParams + jElement*nElePosParams+7] =
          elements[iElement][jElement].euler.gamma;
      }
    }

    // Assign output
    *nElements    = nOutputElements;
    *nSubElements = nOutputSubElementsPerElement;
    *nParams      = Aperture<T>::nElementPosParameters;
    *out          = arr;
  }
  return retval;
}

template <class T>
int Aperture<T>::ElementsSet(
  const T* pos, const size_t nPositions, const size_t nDim) {
  int retval = this->SubElementsSet(pos, nPositions, 1, nDim);
  if (retval) {
    // Causes SWIG to leak memory
    throw std::runtime_error("Elements improperly formatted");
  }
  return retval;
}

template <class T>
int Aperture<T>::SubElementsSet(
  const T* pos, const size_t nElements,
  const size_t nSubElementsPerElement, const size_t nDim) {
  int retval = 0;

  // We are not allowed to set 0 elements
  if ((nDim != Aperture<T>::nElementPosParameters) ||
      (nSubElementsPerElement < 1)) {
    retval = -1;
    return retval;
  }

  // Invalidate focus
  m_pData->m_focus_valid = FocusingType::FocusingTypeCount;

  element_array<T> elements =
    sps::unique_aligned_multi_array_create<sps::element_rect_t<T>, 2>(nElements, nSubElementsPerElement);

  const size_t nElePosParams = Aperture<T>::nElementPosParameters;

  for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
    for (size_t jElement = 0 ; jElement < nSubElementsPerElement ; jElement++) {
      elements[iElement][jElement].hw =
        pos[iElement * nSubElementsPerElement * nElePosParams +
                     jElement * nElePosParams + 0];
      // Error was here
      elements[iElement][jElement].hh =
        pos[iElement * nSubElementsPerElement * nElePosParams +
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
  this->m_pData->ElementsSet(
    std::forward<element_array<T> >(elements),
    nElements, nSubElementsPerElement);

  return retval;
}

#ifdef USE_PROGRESS_BAR
template <class T>
void Aperture<T>::ProgressBarSet(sps::ProgressBarInterface* pbar) {
  this->m_pbar = pbar;
}
#endif

#ifdef FNM_CLOSURE_FUNCTIONS

template <class T>
class FunctionInfo {
 public:
  RwParamType rw;
  const char name[32];
  Type type;
  size_t nDims;
  bool allocates;
  std::function<void(Aperture<T>*)> get_ptr;
  std::function<void(Aperture<T>*)> set_ptr;
};

template <class T>
class FunctionSelectTable {
 public:
  static FunctionInfo<T> Info[];

  /* Thread local variables: TODO: Move to object or data object */
  static __THREAD size_t idims[3];
  static __THREAD size_t* odims[3];
};

template <class U>
class FunctionArguments {
 public:
  static __THREAD U value;
  static __THREAD U* ptr;        ///< Non-const pointer for getting
  static const __THREAD U* iPtr; ///< Const pointer for setting
};

template <class U>
__THREAD U FunctionArguments<U>::value = U(0);

template <class U>
__THREAD U* FunctionArguments<U>::ptr = nullptr;

template <class U>
__THREAD const U* FunctionArguments<U>::iPtr = nullptr;

template class FunctionArguments<bool>;
template class FunctionArguments<int>;
template class FunctionArguments<size_t>;
template class FunctionArguments<float>;
template class FunctionArguments<double>;

template <class T>
int Aperture<T>::
ParameterInfoGet(const RwParamType& param, Type* type, size_t* nDims) {
  int retval = -1;
  if (param < RwParamType::RwParamTypeCount) {
    *type  = FunctionSelectTable<T>::Info[param].type;
    *nDims = FunctionSelectTable<T>::Info[param].nDims;
    return 0;
  }
  return retval;
}

template <class T>
int Aperture<T>::
ParameterInfoListGet(ApertureProperty** ppPropertyInfo) {
  *ppPropertyInfo =
    (ApertureProperty*) malloc(RwParamType::RwParamTypeCount*sizeof(ApertureProperty));
  for (int i = 0 ; i < RwParamType::RwParamTypeCount ; ++i) {
    (*ppPropertyInfo)[i].id = (RwParamType) i;
    strncpy(
      (*ppPropertyInfo)[i].name,
      FunctionSelectTable<T>::Info[i].name, 31);
    (*ppPropertyInfo)[i].name[31] = '\0';
  }
  return RwParamType::RwParamTypeCount;
}

template <class T>
int Aperture<T>::
ParameterInfoListDestroy(ApertureProperty* ppPropertyInfo) {
  free(ppPropertyInfo);
  return 0;
}

template <class T>
__THREAD size_t FunctionSelectTable<T>::idims[3] = {0, 0, 0};

template <class T>
__THREAD size_t* FunctionSelectTable<T>::odims[3] = {nullptr, nullptr, nullptr};

template <class T>
FunctionInfo<T> FunctionSelectTable<T>::Info[] = {
  {
    RwParamType::ElementDelays, "ElementDelays", Type::Float, 1, true,
    //      TODO: Make this work using RwFloatParamGetUsingThis (use std::ref() or std::cref())
    //      [](Aperture<T>* a)->void{a->DelaysGet(&(a->m_pData->c_ptr), a->m_pData->c_odims[0]);},
    //      [](Aperture<T>* a)->void{a->DelaysSet(a->m_pData->c_iPtr, a->m_pData->c_idims[0]); }
    [](Aperture<T>* a)->void {a->DelaysGet(&(fnm::FunctionArguments<T>::ptr), fnm::FunctionSelectTable<T>::odims[0]);},
    [](Aperture<T>* a)->void {a->DelaysSet(fnm::FunctionArguments<T>::iPtr, fnm::FunctionSelectTable<T>::idims[0]); }
  },
  {
    RwParamType::AttenuationEnabled, "AttenuationEnabled", Type::Bool, 0, false,
    [](Aperture<T>* a)->void {*FunctionArguments<bool>::ptr = a->AttenuationEnabledGet();},
    [](Aperture<T>* a)->void {a->AttenuationEnabledSet(FunctionArguments<bool>::value);},
  },
  {
    RwParamType::Alpha, "Alpha", Type::Float, 0, false,
    [](Aperture<T>* a)->void{*FunctionArguments<T>::ptr = a->AlphaGet();},
    [](Aperture<T>* a)->void{a->AlphaSet(FunctionArguments<T>::value);}
  },
  {
    RwParamType::Beta, "Beta", Type::Float, 0, false,
    [](Aperture<T>* a)->void{*FunctionArguments<T>::ptr = a->BetaGet();},
    [](Aperture<T>* a)->void{a->BetaSet(FunctionArguments<T>::value);}
  },
  {
    RwParamType::Positions, "Positions", Type::Float, 2, true,
    [](Aperture<T>* a)->void{
      a->PositionsGet(&(FunctionArguments<T>::ptr),
                      FunctionSelectTable<T>::odims[0],
                      FunctionSelectTable<T>::odims[1]);
    },
    [](Aperture<T>* a)->void{
      a->PositionsSet(FunctionArguments<T>::iPtr,
                      FunctionSelectTable<T>::idims[0],
                      FunctionSelectTable<T>::idims[1]);
    },
  },
  {
    RwParamType::Focus, "Focus", Type::Float, 1, true,
    [](Aperture<T>* a)->void{
      typedef T tre[3];
      tre* output = (tre*) malloc(3*sizeof(T));
      a->FocusGet(*output);
      FunctionArguments<T>::ptr = (T*) output;
      // Need to set output dimensions
      *FunctionSelectTable<T>::odims[0] = 3;
    },
    [](Aperture<T>* a)->void{
      typedef T tre[3];
      const tre* output = (tre*) FunctionArguments<T>::iPtr;
      assert(FunctionSelectTable<T>::idims[0] == 3);
      a->FocusSet(*output);
    },
  },
  {
    RwParamType::CenterFocus, "CenterFocus", Type::Float, 1, true,
    [](Aperture<T>* a)->void{
      typedef T tre[3];
      tre* output = (tre*) malloc(3*sizeof(T));
      a->CenterFocusGet(*output);
      FunctionArguments<T>::ptr = (T*) output;
      // Need to set output dimensions
      *FunctionSelectTable<T>::odims[0] = 3;
    },
    [](Aperture<T>* a)->void{
      typedef T tre[3];
      const tre* output = (tre*) FunctionArguments<T>::iPtr;
      assert(FunctionSelectTable<T>::idims[0] == 3);
      a->CenterFocusSet(*output);
    },
  },
  {
    RwParamType::FocusingType, "FocusingType", Type::Int32, 0, false,
    [](Aperture<T>* a)->void{*FunctionArguments<int>::ptr = a->FocusingTypeGet();},
    [](Aperture<T>* a)->void{a->FocusingTypeSet(FunctionArguments<int>::value);}
  },
  {
    RwParamType::NThreads, "NThreads", Type::SizeT, 0, false,
    [](Aperture<T>* a)->void{*FunctionArguments<size_t>::ptr = a->NThreadsGet();},
    [](Aperture<T>* a)->void{a->NThreadsSet(FunctionArguments<size_t>::value);}
  },
  {
    RwParamType::F0, "F0", Type::Float, 0, false,
    [](Aperture<T>* a)->void{*FunctionArguments<T>::ptr = a->F0Get();},
    [](Aperture<T>* a)->void{a->F0Set(FunctionArguments<T>::value);}
  },
  {
    RwParamType::W, "W", Type::Float, 0, false,
    [](Aperture<T>* a)->void{*FunctionArguments<T>::ptr = a->WGet();},
    [](Aperture<T>* a)->void{a->WSet(FunctionArguments<T>::value);}
  },
  {
    RwParamType::SysParm, "SysParm", Type::Struct,  0, false, nullptr, nullptr
  },
  {
    RwParamType::C, "C", Type::Float, 0, false,
    [](Aperture<T>* a)->void{*FunctionArguments<T>::ptr = a->CGet();},
    [](Aperture<T>* a)->void{a->CSet(FunctionArguments<T>::value);}
  },
  {
    RwParamType::Elements, "Elements", Type::Float, 2, true,
    [](Aperture<T>* a)->void{
      a->ElementsGet(&(FunctionArguments<T>::ptr),
                     FunctionSelectTable<T>::odims[0],
                     FunctionSelectTable<T>::odims[1]);
    },
    [](Aperture<T>* a)->void{
      a->ElementsSet(FunctionArguments<T>::iPtr,
                     FunctionSelectTable<T>::idims[0],
                     FunctionSelectTable<T>::idims[1]);
    },
  },
  {
    RwParamType::SubElements, "SubElements", Type::Float, 3, true,
    [](Aperture<T>* a)->void{
      a->SubElementsGet(&(FunctionArguments<T>::ptr),
                        FunctionSelectTable<T>::odims[0],
                        FunctionSelectTable<T>::odims[1],
                        FunctionSelectTable<T>::odims[2]);
    },
    [](Aperture<T>* a)->void{
      a->SubElementsSet(FunctionArguments<T>::iPtr,
                        FunctionSelectTable<T>::idims[0],
                        FunctionSelectTable<T>::idims[1],
                        FunctionSelectTable<T>::idims[2]);
    },
  },
  {
    RwParamType::Apodization, "Apodization", Type::Float, 1, false,
    [](Aperture<T>* a)->void{a->ApodizationGet(&(FunctionArguments<T>::ptr), FunctionSelectTable<T>::odims[0]);},
    [](Aperture<T>* a)->void{a->ApodizationSet(FunctionArguments<T>::iPtr, FunctionSelectTable<T>::idims[0]);}
  },
  {
    RwParamType::NDivW, "NDivW", Type::SizeT, 0, false,
    [](Aperture<T>* a)->void{*FunctionArguments<size_t>::ptr = a->NDivWGet();},
    [](Aperture<T>* a)->void{a->NDivWSet(FunctionArguments<size_t>::value);}
  },
  {
    RwParamType::NDivH, "NDivH", Type::SizeT, 0, false,
    [](Aperture<T>* a)->void{*FunctionArguments<size_t>::ptr = a->NDivHGet();},
    [](Aperture<T>* a)->void{a->NDivHSet(FunctionArguments<size_t>::value);}
  },
#if FNM_PULSED_WAVE
  {
    RwParamType::Fs, "Fs", Type::Float, 0, false,
    [](Aperture<T>* a)->void{*FunctionArguments<T>::ptr = a->FsGet();},
    [](Aperture<T>* a)->void{a->FsSet(FunctionArguments<T>::value);}
  },
  {
    RwParamType::Normalize, "Normalize", Type::Bool, 0, false,
    [](Aperture<T>* a)->void {*FunctionArguments<bool>::ptr = a->NormalizeGet();},
    [](Aperture<T>* a)->void {a->NormalizeSet(FunctionArguments<bool>::value);},
  },
  {
    RwParamType::Excitation, "Excitation", Type::Float, 1, false,
    [](Aperture<T>* a)->void{a->ExcitationGet(&(FunctionArguments<T>::ptr), FunctionSelectTable<T>::odims[0]);},
    [](Aperture<T>* a)->void{a->ExcitationSet(FunctionArguments<T>::iPtr, FunctionSelectTable<T>::idims[0]);}
  },
  {
    RwParamType::Impulse, "Impulse", Type::Float, 1, false,
    [](Aperture<T>* a)->void{a->ImpulseGet(&(FunctionArguments<T>::ptr), FunctionSelectTable<T>::odims[0]);},
    [](Aperture<T>* a)->void{a->ImpulseSet(FunctionArguments<T>::iPtr, FunctionSelectTable<T>::idims[0]);}
  },
  {
    RwParamType::Focus2, "Focus2", Type::Float, 1, true,
    [](Aperture<T>* a)->void{
      typedef T tre[3];
      tre* output = reinterpret_cast<tre*>(SPS_MALLOC(3*sizeof(T)));
      a->Focus2Get(*output);
      FunctionArguments<T>::ptr = (T*) output;
      // Need to set output dimensions
      *FunctionSelectTable<T>::odims[0] = 3;
    },
    [](Aperture<T>* a)->void{
      typedef T tre[3];
      const tre* output = reinterpret_cast<const tre*>(FunctionArguments<T>::iPtr);
      assert(FunctionSelectTable<T>::idims[0] == 3);
      a->Focus2Set(*output);
    },
  },
#endif
};

#ifdef FNM_CLOSURE_FUNCTIONS

template <class T>
int Aperture<T>::RwFloatParamGet(int fsel, T** oMultiData, size_t nDim, ...) {
  va_list args;
  va_start(args, nDim);
  int retval = this->RwFloatParamGet(fsel, oMultiData, nDim, args);
  va_end(args);
  return retval;
}

template <class T>
int Aperture<T>::RwFloatParamSet(int fsel, const T* iMultiData, size_t nDim, ...) {
  va_list args;
  va_start(args, nDim);
  int retval = this->RwFloatParamSet(fsel, iMultiData, nDim, args);
  va_end(args);
  return retval;
}

template <class T>
int Aperture<T>::RwFloatParamSet(
  int pSel, const T* f, size_t nDim, va_list args) {
  int retval = -1;

  for (size_t i = 0 ; i < 3 ; i++) {
    FunctionSelectTable<T>::idims[i] = 0;
  }

  auto pred = [&](const FunctionInfo<T>& item) {
    return item.rw == pSel;
  };

#ifdef __GNUG__
  auto it = std::find_if(std::begin(FunctionSelectTable<T>::Info),
                         std::end(FunctionSelectTable<T>::Info), pred);
  if (it != std::end(FunctionSelectTable<T>::Info))
#else
  auto it = std::find_if(FunctionSelectTable<T>::Info,
                         FunctionSelectTable<T>::Info + 22, pred);
  if (it != (FunctionSelectTable<T>::Info + 22))
#endif
  {
#ifdef __GNUG__
    // Index for parameters
    int fsel =
      static_cast<int>(std::distance(std::begin(FunctionSelectTable<T>::Info), it));
#else
    int fsel =
      static_cast<int>(std::distance(FunctionSelectTable<T>::Info, it));
#endif
    // Verify correct number of dimensions
    if ((FunctionSelectTable<T>::Info[fsel].nDims != nDim) || (nDim > 3)) {
      // How to set empty array
      printf("Wrong dimensions: %zu\n", nDim);
      retval = -1;
      return retval;
    }

    // Read input dimensions
    for (size_t i = 0; i < nDim ; i++) {
      FunctionSelectTable<T>::idims[i] = va_arg(args, size_t);
    }

    // No need to allocate for setting a scalar
    FunctionArguments<T>::value = *f;

    // Pointer used for inputs of 1, 2 dimensions
    FunctionArguments<T>::iPtr = f;

    // Check that the parameters are floating point (i.e. T)
    if (FunctionSelectTable<T>::Info[fsel].type == Type::Float) {
      FunctionSelectTable<T>::Info[fsel].set_ptr(this);
      retval = 0;
    } else {
      retval = -1;
    }
  }

  return retval;
}

template <class T>
int Aperture<T>::RwFloatParamGet(int pSel, T** f, size_t nDim, va_list args) {
  debug_print("nDim: %zu\n", nDim);
  int retval = 0;

  // Reset dimensions and pointers
  for (size_t i = 0 ; i < 3 ; i++) {
    if (FunctionSelectTable<T>::odims[i]) {
      *(FunctionSelectTable<T>::odims[i]) = 0;
    }
    FunctionSelectTable<T>::odims[i] = nullptr;
  }
  FunctionArguments<T>::ptr = nullptr;

#ifdef __GNUC__
  auto pred = [&](const FunctionInfo<T>& item) {
    return item.rw == pSel;
  };

  auto it = std::find_if(std::begin(FunctionSelectTable<T>::Info),
                         std::end(FunctionSelectTable<T>::Info), pred);

  if (it != std::end(FunctionSelectTable<T>::Info)) {
    // Index for parameters
    int fsel =
      static_cast<int>(std::distance(std::begin(FunctionSelectTable<T>::Info), it));
#else
  {
    int fsel = pSel;
#endif

    // Verify correct number of dimensions
    if ((FunctionSelectTable<T>::Info[fsel].nDims != nDim) || (nDim > 3)) {
      *f = NULL;
      retval = -1;
      return retval;
    }

    // Read pointers for output dimensions
    debug_print("nDim: %zu\n", nDim);
    for (size_t i = 0; i < nDim ; i++) {
      FunctionSelectTable<T>::odims[i] = va_arg(args, size_t*);
    }

    if (FunctionSelectTable<T>::Info[fsel].type == Type::Float) {
      // Dimensions 1 or 2
      if (FunctionSelectTable<T>::Info[fsel].nDims > 0) {
        // Call function
        FunctionSelectTable<T>::Info[fsel].get_ptr(this);

        // Assign output
        *f = FunctionArguments<T>::ptr;

        auto _odims = FunctionSelectTable<T>::odims;
        debug_print("odims[0]: %zu, odims[1]: %zu, odims[2]: %zu\n",
                    _odims[0] != nullptr ? *(_odims[0]) : 0,
                    _odims[1] != nullptr ? *(_odims[1]) : 0,
                    _odims[2] != nullptr ? *(_odims[2]) : 0);

        // If the function doesn't allocate and returns a view, allocate and copy
        if (!FunctionSelectTable<T>::Info[fsel].allocates) {
          // Compute size
          size_t nData = 1;
          for (size_t i = 0 ; i < nDim ; i++) {
            nData *= *(FunctionSelectTable<T>::odims[i]);
          }
          // Allocate
          *f = static_cast<T*>(SPS_MALLOC(sizeof(T)*nData));
          // Copy
          debug_print("FunctionArguments<T>::ptr: %f\n", *FunctionArguments<T>::ptr);
          memcpy(*f, FunctionArguments<T>::ptr, nData*sizeof(T));
        }
      } else {
        // Scalar, we always allocate
        *f = static_cast<T*>(SPS_MALLOC(sizeof(T)));
        FunctionArguments<T>::ptr = *f;

        // Call function
        FunctionSelectTable<T>::Info[fsel].get_ptr(this);
      }
    }
  }

  return retval;
}

template <class T>
int Aperture<T>::RwFloatParamGetUsingThis(
  int pSel, T** f, size_t nDim, va_list args) {
  debug_print("nDim: %zu\n", nDim);
  int retval = 0;

  // Reset dimensions and pointers
  for (size_t i = 0 ; i < 3 ; i++) {
    if (this->m_pData->c_odims[i]) {
      *(this->m_pData->c_odims[i]) = 0;
    }
    this->m_pData->c_odims[i] = nullptr;
  }
  this->m_pData->c_ptr = nullptr;

#ifdef __GNUC__
  auto pred = [&](const FunctionInfo<T>& item) {
    return item.rw == pSel;
  };

  auto it = std::find_if(std::begin(FunctionSelectTable<T>::Info),
                         std::end(FunctionSelectTable<T>::Info), pred);

  if (it != std::end(FunctionSelectTable<T>::Info)) {
    // Index for parameters
    int fsel = (int)std::distance(std::begin(FunctionSelectTable<T>::Info), it);
#else
  {
    int fsel = pSel;
#endif
    // Verify correct number of dimensions
    if ((FunctionSelectTable<T>::Info[fsel].nDims != nDim) || (nDim > 3)) {
      *f = NULL;
      retval = -1;
      return retval;
    }

    // Read pointers for output dimensions
    debug_print("nDim: %zu\n", nDim);
    for (size_t i = 0; i < nDim ; i++) {
      this->m_pData->c_odims[i] = va_arg(args, size_t*);
    }

    if (FunctionSelectTable<T>::Info[fsel].type == Type::Float) {
      // Dimensions 1 or 2
      if (FunctionSelectTable<T>::Info[fsel].nDims > 0) {
        // Call function
        FunctionSelectTable<T>::Info[fsel].get_ptr(this);

        // Assign output
        *f = this->m_pData->c_ptr;

        debug_print("c_odims[0]: %zu, c_odims[1]: %zu, c_odims[2]: %zu\n",
                    this->m_pData->c_odims[0] != nullptr ? *(this->m_pData->c_odims[0]) : 0,
                    this->m_pData->c_odims[1] != nullptr ? *(this->m_pData->c_odims[1]) : 0,
                    this->m_pData->c_odims[2] != nullptr ? *(this->m_pData->c_odims[2]) : 0);

        //  If the function doesn't allocate and returns a view, allocate and copy
        if (!FunctionSelectTable<T>::Info[fsel].allocates) {
          // Compute size
          size_t nData = 1;
          for (size_t i = 0 ; i < nDim ; i++) {
            nData *=  *(this->m_pData->c_odims[i]);
          }
          // Allocate
          *f = static_cast<T*>(SPS_MALLOC(sizeof(T)*nData));
          // Copy
          debug_print("this->m_pData->c_ptr: %f\n", *this->m_pData->c_ptr);
          memcpy(*f, this->m_pData->c_ptr, nData*sizeof(T));
        }
      } else {
        // Update static variables
        *f = static_cast<T*>(SPS_MALLOC(sizeof(T)));
        this->m_pData->c_ptr = *f;

        // Call function
        FunctionSelectTable<T>::Info[fsel].get_ptr(this);
      }
    }
  }

  return retval;
}
#endif

template <class T>
void Aperture<T>::RwBooleanParamSet0D(int fsel, const bool& value) {
  if (FunctionSelectTable<T>::Info[fsel].type == Type::Bool) {
    FunctionArguments<bool>::value = value;
    FunctionSelectTable<T>::Info[fsel].set_ptr(this);
  }
}

template <class T>
void Aperture<T>::RwSizeTParamSet0D(int fsel, const size_t& value) {
  if (FunctionSelectTable<T>::Info[fsel].type == Type::SizeT) {
    FunctionArguments<size_t>::value = value;
    FunctionSelectTable<T>::Info[fsel].set_ptr(this);
  }
}

template <class T>
void Aperture<T>::RwIntegerParamSet0D(int fsel, const int& value) {
  if (FunctionSelectTable<T>::Info[fsel].type == Type::Int32) {
    FunctionArguments<int>::value = value;
    FunctionSelectTable<T>::Info[fsel].set_ptr(this);
  }
}
#endif

template <class T>
void Aperture<T>::FocusUpdateRef() {
  const size_t nElements    = m_pData->m_nelements;
  const size_t nSubElements = m_pData->m_nsubelements;
  const T c = m_sysparm->c;

  if (m_pData->m_focus_type == FocusingType::Pythagorean) {
    // Geometric focusing - allocate here,
    // since we consider removing m_pData->m_delays
    auto delays = sps::unique_aligned_array_create<T>(nElements);

    T maxDelay = - std::numeric_limits<T>::max();
    T minDelay =  std::numeric_limits<T>::max();

    for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
      delays[iElement] =
        sps::norm(m_pData->m_pos[iElement] - m_pData->m_focus) / c;
      maxDelay = std::max<T>(maxDelay, delays[iElement]);
      minDelay = std::min<T>(minDelay, delays[iElement]);
    }
    // Using center_focus for computing minDelay (TEST - should we do this)
    minDelay = norm(m_pData->m_focus - m_pData->m_center_focus) / c;

    for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
      m_pData->m_delays[iElement] = minDelay - delays[iElement];

      // TODO(JEM): Verify if it is okay, if delays and phases are opposite
      for (size_t jElement = 0 ; jElement < nSubElements ; jElement++) {
        m_pData->m_phases[iElement*nSubElements+jElement] =
          T(M_2PI) * (delays[iElement] - maxDelay) * m_pData->m_f0;
      }
    }
  } else if (m_pData->m_focus_type == FocusingType::Rayleigh) {
    // Rayleigh focusing (CW)
    size_t nPositions = nElements;
    auto positions = sps::unique_aligned_array_create<T>(3*nPositions);
    for (size_t iPosition = 0 ; iPosition < nPositions ; iPosition++) {
      positions[3*iPosition]   = m_pData->m_focus[0];
      positions[3*iPosition+1] = m_pData->m_focus[1];
      positions[3*iPosition+2] = m_pData->m_focus[2];
    }

    std::complex<T>* pFieldValues = NULL;
    size_t nFieldValues;

    // Note, this function allocates!!!
    this->CalcCwFocusRef(positions.get(), nPositions, 3,
                         &pFieldValues, &nFieldValues);

    for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
      for (size_t jElement = 0 ; jElement < nSubElements ; jElement++) {
        T phase = -std::arg<T>(pFieldValues[iElement]);
        m_pData->m_phases[iElement*nSubElements+jElement] = phase;
        debug_print("phase: %f\n", phase);
      }
    }
    if (pFieldValues) {
      free(pFieldValues);
    }
  } else if (m_pData->m_focus_type == FocusingType::Delays) {
    // Update phases using delays
    printf("Update phases\n");
    for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
      for (size_t jElement = 0 ; jElement < nSubElements ; jElement++) {
        m_pData->m_phases[iElement*nSubElements+jElement] =
          - std::fmod<T>(m_pData->m_f0*m_pData->m_delays[iElement],
                         T(1.0)) * T(M_2PI);
      }
    }
  }

  // Set validity
  m_pData->m_focus_valid = m_pData->m_focus_type;
}

template <class T>
void Aperture<T>::FocusUpdate() {

  // Ups, we need to introduce requested and actual
  if ((m_pData->m_focus_valid == m_pData->m_focus_type) && false) {
    debug_print("We cannot skip this if elements, c, or f0 are changed");
    return;
  }

  const size_t nElements    = m_pData->m_nelements;
  const size_t nSubElements = m_pData->m_nsubelements;
  const T c = m_sysparm->c;

  const sps::point_t<T>& focus = m_pData->m_focus;
  const sps::point_t<T>& focus2 = m_pData->m_focus2;
  const sps::point_t<T>& centerFocus = m_pData->m_center_focus;
  const T signFocus = std::copysign<T>(T(1.0), focus[2]);

  const auto& pos = m_pData->m_pos;

  if (m_pData->m_focus_type == FocusingType::Pythagorean) {
    // TODO(JMH): Rename this
    PytDelays(pos, nElements,
              centerFocus,
              focus,
              focus2,
              c,
              m_pData->m_delays,
             false);

    for (size_t iElement=0 ; iElement < nElements ; iElement++) {
      for (size_t jElement = 0 ; jElement < nSubElements ; jElement++) {
        // Phases are set for each sub-element
        m_pData->m_phases[iElement*nSubElements+jElement] =
          std::fmod<T>(-m_pData->m_f0*m_pData->m_delays[iElement], T(1.0))*T(M_2PI);
      }
    }
  } else if (m_pData->m_focus_type == FocusingType::Rayleigh) {
    // Rayleigh focusing (CW)
    size_t nPositions;
#if FNM_PHASED_FOCUS
    nPositions = nElements*nSubElements;
#else
    nPositions = nElements;
#endif
    auto positions = sps::unique_aligned_array_create<T>(3*nPositions);
    for (size_t iPosition = 0 ; iPosition < nPositions ; iPosition++) {
      positions[3*iPosition]   = focus[0];
      positions[3*iPosition+1] = focus[1];
      positions[3*iPosition+2] = focus[2];
    }

    std::complex<T>* pFieldValues = NULL;
    size_t nFieldValues;

    this->CalcCwFocus(positions.get(), nPositions, 3,
                      &pFieldValues, &nFieldValues);

    for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
      for (size_t jElement = 0 ; jElement < nSubElements ; jElement++) {
        T phase;
#if FNM_PHASED_FOCUS
        // Experimenting with setting phase per sub-element
        phase = - std::arg<T>(pFieldValues[iElement*nSubElements + jElement]);
#else
        phase = - std::arg<T>(pFieldValues[iElement]);
#endif
        m_pData->m_phases[iElement*nSubElements+jElement] = signFocus * phase;
        debug_print("phase: %f\n", phase);
      }
    }
    free(pFieldValues);
  } else if (m_pData->m_focus_type == FocusingType::Delays) {
    // Update phases using delays
    for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
      for (size_t jElement = 0 ; jElement < nSubElements ; jElement++) {
        // [-2pi;2pi]
        m_pData->m_phases[iElement*nSubElements+jElement] =
          - std::fmod<T>(m_pData->m_f0*m_pData->m_delays[iElement],
                         T(1.0))*T(M_2PI);
      }
    }
  }

  // Set validity
  m_pData->m_focus_valid = m_pData->m_focus_type;
}

template <class T>
int Aperture<T>::CalcCwFieldNaive(
  const T* pos, const size_t nPositions, const size_t nDim,
  std::complex<T>** odata, size_t* nOutPositions) {
  assert(nDim == 3);
  if (nDim != 3) {
    *odata = NULL;
    *nOutPositions = 0;
    return -1;
  }

  size_t nScatterers = nPositions;
  *odata = (std::complex<T>*) malloc(nScatterers*sizeof(std::complex<T>));
  *nOutPositions = nScatterers;

  // Works if replaced by FocusUpdateRef
  this->FocusUpdate();

  fnm::CalcCwField<T>(this->m_sysparm, *this->m_pData,
                      pos, nPositions,
                      odata);
  return 0;
}

// Even better integration range
template <class T>
int Aperture<T>::CalcCwFieldFourRef(
  const T* pos, const size_t nPositions, const size_t nDim,
  std::complex<T>** odata, size_t* nOutPositions) {
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

  // This causes deviations
  this->FocusUpdateRef();

  int tStart = fnm::CalcCwFieldFourRef<T>(this->m_sysparm, this->m_pData,
                                          pos, nPositions,
                                          odata);
  return tStart;
}

template <class T>
T Aperture<T>::CalcPwFnmThreaded(
  const T* pos, const size_t nPositions, const size_t nDim,
  T** odata, size_t* nSignals, size_t* nSamples, int mask) {
  assert(nDim == 3);
  if (nDim != 3) {
    *odata = NULL;
    *nSignals = 0;
    *nSamples = 0;
    return T(0.0);
  }

  if ((m_pData->m_focus_type != FocusingType::Pythagorean) &&
      (m_pData->m_focus_type != FocusingType::Delays)) {
    fprintf(stderr, "Focusing type must be Pythagorean or set using "
            "delays for pulsed wave calculations\n");
  }
  this->FocusUpdate();

  T tstart;
  if (ExcitationTypeGet() == ExcitationType::ExcitationTypeHanningWeightedPulse) {
    tstart = fnm::CalcPwFnmThreaded<T, HanningWeightedPulse>(
               this->m_sysparm, this->m_pData,
               pos, nPositions, odata, nSamples, mask, m_pbar);
  } else if (ExcitationTypeGet() == ExcitationType::ExcitationTypeToneBurst) {
    tstart = fnm::CalcPwFnmThreaded<T, ToneBurst>(
               this->m_sysparm, this->m_pData,
               pos, nPositions, odata, nSamples, mask, m_pbar);
  } else {
    tstart = fnm::CalcPwFnmThreaded<T, Identity>(
               this->m_sysparm, this->m_pData,
               pos, nPositions, odata, nSamples, mask, m_pbar);
  }
  *nSignals = nPositions;
  return tstart;
}

template <class T>
T Aperture<T>::CalcTransientSingleElementNoDelay(
  const T* pos, const size_t nPositions, const size_t nDim,
  T** odata, size_t* nSignals, size_t* nSamples, int mask) {
  assert(nDim == 3);
  if (nDim != 3) {
    *odata = NULL;
    *nSignals = 0;
    *nSamples = 0;
    return T(0.0);
  }
  T tstart = T(0.0);
  if (this->ExcitationTypeGet() == ExcitationType::ExcitationTypeHanningWeightedPulse) {
    tstart = fnm::TransientSingleRect<T, HanningWeightedPulse>(
                                                               this->m_sysparm, this->m_pData, pos, nPositions, odata, nSamples, mask);
  } else {
    tstart = fnm::TransientSingleRect<T, ToneBurst>(
                                                    this->m_sysparm, this->m_pData, pos, nPositions, odata, nSamples, mask);
  }
  *nSignals = nPositions;
  return tstart;
}

// Optimal sampling for integral (reduced integration path)
template <class T>
int Aperture<T>::CalcCwFieldRef(
  const T* pos, const size_t nPositions, const size_t nDim,
  std::complex<T>** odata, size_t* nOutPositions) {

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

  // This causes deviations
  this->FocusUpdateRef();

  int tStart = fnm::CalcCwFieldRef<T>(this->m_sysparm, this->m_pData,
                                      pos, nPositions,
                                      odata, m_pbar);
  return tStart;
}

template <class T>
int Aperture<T>::CalcCwEcho(const Aperture<T>* pOther,
                            const T* pos, const size_t nPositions, const size_t nDim,
                            std::complex<T>** odata, size_t* nOutPositions) {
  int retval = 0;

  assert(nDim == 3);
  if (nDim != 3) {
    *odata = NULL;
    *nOutPositions = 0;
    retval = -1;
    return retval;
  }

  // std::complex<T>* pFieldData  = (std::complex<T>*) malloc(nPositions*sizeof(std::complex<T>));
  // memset(pFieldData, 0, nPositions*sizeof(std::complex<T>));
  std::complex<T>* pFieldData  = (std::complex<T>*) calloc(nPositions, sizeof(std::complex<T>));

  this->FocusUpdate();

  retval = fnm::CalcCwThreaded<T>(this->m_sysparm, this->m_pData, pos, nPositions, &pFieldData, this->m_pbar);

  const size_t nElements = pOther->m_pData->m_nelements;

  //*odata = (std::complex<T>*) malloc(nElements*sizeof(std::complex<T>));
  //memset(*odata, 0, nElements*sizeof(std::complex<T>));
  *odata = (std::complex<T>*) calloc(nElements, sizeof(std::complex<T>));
  *nOutPositions = nElements;

  retval |= fnm::CalcCwBackThreaded<T>(
              this->m_sysparm, pOther->m_pData, pos, nPositions,
              pFieldData, nPositions,
              odata, this->m_pbar);

  return retval;
}

template <class T>
int Aperture<T>::CalcCwBack(const T* pos,
                            const size_t nPositions,
                            const size_t nDim,
                            const std::complex<T>* pFieldValues,
                            const size_t nComplexValues,
                            std::complex<T>** odata,
                            size_t* nOutPositions) {
  int retval = 0;

  assert(nDim == 3);
  if (nDim != 3  || nPositions != nComplexValues) {
    *odata = NULL;
    *nOutPositions = 0;
    retval = -1;
    return retval;
  }

  const size_t nElements = this->m_pData->m_nelements;

  //*odata = (std::complex<T>*) malloc(nElements*sizeof(std::complex<T>));
  //memset(*odata, 0, nElements*sizeof(std::complex<T>));
  *odata = (std::complex<T>*) calloc(nElements, sizeof(std::complex<T>));

  *nOutPositions = nElements;

  retval = fnm::CalcCwBackThreaded<T>(
             this->m_sysparm, this->m_pData, pos, nPositions,
             pFieldValues, nComplexValues,
             odata, this->m_pbar);

  return retval;
}

template <class T>
int Aperture<T>::CalcCwTimeReversal(
  const Aperture<T>* pOther,
  const T* pos, const size_t nPositions, const size_t nDim,
  const std::complex<T>* pFieldValues, const size_t nComplexValues,
  std::complex<T>** odata, size_t* nOutPositions, size_t* nSamples) {
  int retval = 0;
  assert(nDim == 3);
  const size_t nElements0 = this->m_pData->m_nelements;
  const size_t nElements1 = pOther->m_pData->m_nelements;

  if (nDim != 3  || nElements0 != nComplexValues || nElements0 != nElements1) {
    *odata = NULL;
    *nOutPositions = 0;
    retval = -1;
    return retval;
  }

  /*
  *odata =
    static_cast<std::complex<T>*>(malloc(nPositions*nPositions*sizeof(std::complex<T>)));
  memset(*odata, 0, nPositions*nPositions*sizeof(std::complex<T>));
  */

  *odata =
    static_cast<std::complex<T>*>(calloc(nPositions*nPositions, sizeof(std::complex<T>)));

  *nOutPositions = nPositions;
  *nSamples = nPositions;

  retval = fnm::CalcCwTimeReversal<T>(
             this->m_sysparm, this->m_pData, pOther->m_pData,
             pos, nPositions,
             pFieldValues, nComplexValues,
             odata, this->m_pbar);

  return retval;
}

template <class T>
int Aperture<T>::CalcCwFast(const T* pos,
                            const size_t nPositions,
                            const size_t nDim,
                            std::complex<T>** odata,
                            size_t* nOutPositions) {
  int retval = 0;
  assert(nDim == 3);
  if (nDim != 3) {
    *odata = NULL;
    *nOutPositions = 0;
    retval = -1;
    return retval;
  }

  //*odata = (std::complex<T>*) malloc(nPositions*sizeof(std::complex<T>));
  //memset(*odata, 0, nPositions*sizeof(std::complex<T>));
  *odata = (std::complex<T>*) calloc(nPositions, sizeof(std::complex<T>));
  *nOutPositions = nPositions;

  this->FocusUpdate();

  retval = fnm::CalcCwThreaded<T>(this->m_sysparm, this->m_pData, pos, nPositions, odata, this->m_pbar);

  // Call back function with complex scatter weight


  return retval;
}

template <class T>
int Aperture<T>::CalcCwFocus(
  const T* pos, const size_t nPositions, const size_t nDim,
  std::complex<T>** odata, size_t* nOutPositions) const {
  int retval = 0;
  assert(nDim == 3);
  if (nDim != 3) {
    *odata = NULL;
    *nOutPositions = 0;
    retval = -1;
    return retval;
  }

  // *odata = (std::complex<T>*) malloc(nPositions*sizeof(std::complex<T>));
  // memset(*odata, 0, nPositions*sizeof(std::complex<T>));
  *odata = (std::complex<T>*) calloc(nPositions, sizeof(std::complex<T>));

  *nOutPositions = nPositions;

  retval = fnm::CalcCwFocus<T>(this->m_sysparm, *this->m_pData, pos, nPositions, odata);

  return retval;
}

template <class T>
int Aperture<T>::CalcCwFocusRef(
  const T* pos, const size_t nPositions, const size_t nDim,
  std::complex<T>** odata, size_t* nOutPositions) {
  int retval = 0;
  assert(nDim == 3);
  if (nDim != 3) {
    *odata = NULL;
    *nOutPositions = 0;
    retval = -1;
    return retval;
  }

  // *odata = (std::complex<T>*) malloc(nPositions*sizeof(std::complex<T>));
  // memset(*odata, 0, nPositions*sizeof(std::complex<T>));
  *odata = (std::complex<T>*) calloc(nPositions, sizeof(std::complex<T>));
  *nOutPositions = nPositions;

  // Phases for each element (or sub-element)
  retval = fnm::CalcCwFocusRef<T>(this->m_sysparm, *this->m_pData, pos, nPositions, odata);

  return retval;
}

#ifndef FNM_DOUBLE_SUPPORT

// Sub-optimal integration range (Old function). TODO: Remove or move to fnm_calc
template <class T>
int Aperture<T>::CalcCwFieldNaiveFast(
  const T* pos, const size_t nPositions, const size_t nDim,
  std::complex<T>** odata, size_t* nOutPositions) {
  int retval = 0;
  assert(nDim == 3);
  if (nDim != 3) {
    *odata = NULL;
    *nOutPositions = 0;
    retval = -1;
    return retval;
  }

  size_t nScatterers = nPositions;
  //*odata = (std::complex<T>*) malloc(nScatterers*sizeof(std::complex<T>));
  //memset(*odata, 0, nScatterers*sizeof(std::complex<T>));
  *odata = (std::complex<T>*) calloc(nScatterers, sizeof(std::complex<T>));
  *nOutPositions = nScatterers;

  this->FocusUpdate();

  retval =
    fnm::CalcCwFocusNaiveFast<T>(this->m_sysparm, *this->m_pData, pos, nPositions, odata);
  return retval;
}

template <class T>
int KSpace(
  const T dx, const size_t Nx, const T dy, const size_t Ny,
  const T lambda, T** data, bool shift = true) {
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
int Aperture<T>::CalcAsa(
  const T* y0, const size_t nx, const size_t ny,
  // Vector of z coordinates
  const T dx, const T dy,
  const size_t Nx, const size_t Ny,
  // Output
  std::complex<T>** p1, size_t *onx, size_t *ony, size_t *onz) {

  int retval = 0;

  const T c    = m_sysparm->c;
  const T f0   = this->m_pData->m_f0;

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
    fftwf_complex* input = (fftwf_complex*) SPS_MM_MALLOC(sizeof(fftwf_complex) * Nx * Ny, 16);
    memset(input, 0, sizeof(fftwf_complex) * Nx * Ny);

    for (size_t j = 0 ; j < ny ; j++) {
      for (size_t i = 0 ; i < nx ; i++) {
        input[j*Nx + i][0] = float(y0[j*nx + i]);
      }
    }
    fftwf_complex* Y0 = (fftwf_complex*) SPS_MM_MALLOC(sizeof(fftwf_complex) * Nx * Ny, 16);

    fftwf_plan forward = fftwf_plan_dft_2d((int)Ny, (int)Nx, input, Y0, -1, FFTW_ESTIMATE);
    fftwf_execute_dft(forward,(fftwf_complex *) input, (fftwf_complex *) Y0);
    fftwf_destroy_plan(forward);

    /////////////////////////////////////////////////////////////////////
    // K-space coordinates: Asymmetric for both N even and odd (shifted)
    /////////////////////////////////////////////////////////////////////

    T** kzspace = (T**) _mm_multi_malloc<T,16>(2,Nx,Ny);
    KSpace<T>(dx,Nx,dy,Ny,lambda,kzspace);

    /////////////////////////////////////////////////////////////////////
    // Propagator (update for each z value)
    /////////////////////////////////////////////////////////////////////

    fftwf_complex* Hprop = (fftwf_complex*) SPS_MM_MALLOC(sizeof(fftwf_complex) * Nx * Ny, 16);
    memset(Hprop, 0, sizeof(fftwf_complex) * Nx * Ny);

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
template class FNM_EXPORT Aperture<double>;
#endif

template class FNM_EXPORT Aperture<float>;

#if defined(FNM_CLOSURE_FUNCTIONS) && defined(__GNUG__)
static_assert(RwParamType::RwParamTypeCount == sizeof(FunctionSelectTable<float>::Info) / sizeof(FunctionInfo<float>),
              "Update table FunctionSelectTable, when adding new parameters in fnm_types.h");
#endif

}  // namespace fnm

// Template instantiation

/*
#ifdef __cplusplus
#define IS_C(x)   extern "C" x ;
#define IS_CPP(x) x ;
#else
#define IS_C(x)   x ;
#define IS_CPP(x)
#endif
With this type of header:

IS_C   (void ArrayList_insert(ArrayList *arrlst, void *data, int i))
IS_CPP (void ArrayList_insert(ArrayList *arrlst, char *data, int i))
IS_CPP (void ArrayList_insert(ArrayList *arrlst, Buffer *data, int i))
*/



/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
