/**
 * @file   fnm_transient_threaded.cpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Mon Sep 25 21:15:36 2017
 *
 * @brief
 *
 * Copyright 2017 Jens Munk Hansen
 *
 */
#include <fnm/config.h>
#include <fnm/fnm_types.hpp>
#include <fnm/fnm_common.hpp>
#include <fnm/fnm_transient_threaded.hpp>

#include <fnm/fnm_basis.hpp>
#include <fnm/fnm_data.hpp>
#include <fnm/fnm_response.hpp>
#include <fnm/fnm_profiling.hpp>

#include <sps/queue.hpp>
#include <sps/smath_types.hpp>

#ifdef FNM_PULSED_WAVE
# include <sofus/rect_int_limits.hpp> // calcProjectionAndLimits
# include <sofus/sofus_calc.hpp>      // ComputeBoxTimes
#endif

#include <sps/profiler.h>
#include <sps/progress.hpp>
#include <sps/sps_threads.hpp>

#include <sps/msignals.hpp>

#ifdef __GNUC__
# include <sched.h>
#endif

namespace fnm {

#ifdef HAVE_PTHREAD_H
// External - defined in fnm_calc.cpp
extern pthread_t threads[N_MAX_THREADS];
extern pthread_attr_t attr;
#endif

template <class T>
struct PwFnmThreadArg {
  /// Pointer to system parameters
  const fnm::sysparm_t<T>* pSysparm;
  /// Pointer to aperture data
  const ApertureData<T>* pApertureData;
  /// Mask for enabling direct and edge responses
  int mask;
  /// Start index of output
  size_t iOutputBegin;
  /// Stop index of output
  size_t iOutputEnd;
  /// Number of positions (all threads)
  size_t nPositions;
  /// Position data
  const T* pPos;
  /// Output
  T* pField;
  /// Size of output data (without temporal - there isn't any)
  size_t nData;
  /// Integration start
  int iSampleStart;
  /// Gauss-Legendre coefficienter
  const GLQuad2D<T>* pGL;
  /// Pointer to queue
  sps::queue<float>* pQueue;
  /// Thread ID
  size_t threadId;
  /// CPU ID
  int cpuId;
};

template <class T,  template <typename> class A>
T CalcPwFnmThreaded(const sysparm_t<T>* pSysparm,
                    const ApertureData<T>* pData,
                    const T* pPos, const size_t nPositions,
                    T** pOutData, size_t* pNSamples,
                    int mask,
                    sps::ProgressBarInterface* pBar) {
  SPS_UNREFERENCED_PARAMETER(pBar);

  ProfilerStart();

  // Initialize output pointer and dimensions
  *pOutData = nullptr;
  *pNSamples = 0;

  // Validate input
  if ((!pPos) || nPositions == 0) {
    return T(0.0);
  }

  // Box encapsulating the aperture
  sps::bbox_t<T> box;
  pData->ExtentGet(&box);

  size_t _nElements = 0;

  const T* pApodizations = nullptr;
  const T* pDelays       = nullptr;

  pData->ApodizationsRefGet(&_nElements, pApodizations);
  pData->DelaysRefGet(&_nElements, pDelays);

#if SPS_DEBUG
  debug_print("nElements: %zu\n", _nElements);
  std::cout << "[";
  for (size_t i = 0 ; i < _nElements ; i++) {
    std::cout << pDelays[i] << " ";
  }
  std::cout << "]" << std::endl;
#endif

  size_t nSamples;  // Number of samples required for most extreme point
  T tStart;        // Time for first sample
  int sStart;      // First absolute sample index

  // TODO(JEM): Get rid of one of the sysparm_t<T> structures
  sofus::sysparm_t<T> sysparm;
  sysparm.fs = pSysparm->fs;
  sysparm.c  = pSysparm->c;
  sysparm.nThreads = pSysparm->nThreads;

  auto uws = sps::unique_aligned_array<T>();
  auto vws = sps::unique_aligned_array<T>();
  auto uxs = sps::unique_aligned_array<T>();
  auto vxs = sps::unique_aligned_array<T>();

  /********************************************
   * Compute abcissas and weights             *
   ********************************************/

  size_t _nSubElements;
  const sps::element_rect_t<T>** pElements = nullptr;
  pData->ElementsRefGet(&_nElements, &_nSubElements, pElements);

  // TODO(JEM): Make this without scale (scale assumes all share hh, hw)
  CalcWeightsAndAbcissaeScaled(pSysparm, pElements[0][0],
                               std::move(uxs), std::move(uws),
                               std::move(vxs), std::move(vws));

  GLQuad2D<T> uv;
  uv.u.xs     = uxs.get();
  uv.u.ws     = uws.get();
  uv.u.nx     = pSysparm->nDivW;
  uv.v.xs     = vxs.get();
  uv.v.ws     = vws.get();
  uv.v.nx     = pSysparm->nDivH;

  // Was this changed???
  sofus::ComputeBoxTimes(
    sysparm,
    box,
    pPos, nPositions,
    pDelays, pApodizations, _nElements,
    &tStart,
    &sStart,
    &nSamples);

  // TODO: Take care if sStart is negative
  debug_print("\033[31;1msStart: %d\033[0m\n", sStart);

  // Add length of pulse to nSamples
  nSamples = nSamples + size_t(pSysparm->w * pSysparm->fs) + 1;

  // Allocate output
  T* _data = static_cast<T*>(malloc(nPositions*nSamples*sizeof(T)));
  memset(_data, 0, nPositions*nSamples*sizeof(T));

#ifndef HAVE_THREAD
# error Multi-threading required
#else
  int nproc = 0;
# if defined(_WIN32)
  unsigned int threadID = 0U;
  HANDLE threads[N_MAX_THREADS];
# endif

  nproc = getncpus();

  setvbuf(stdout, NULL, _IONBF, 0);
  setvbuf(stderr, NULL, _IONBF, 0);

# ifdef HAVE_PTHREAD_H
  CallErr(pthread_attr_init, (&fnm::attr));
  CallErr(pthread_attr_setdetachstate, (&fnm::attr, PTHREAD_CREATE_JOINABLE));
# endif

  sps::queue<float>* progressQueue = new sps::queue<float>();

  PwFnmThreadArg<T> threadarg[N_MAX_THREADS];

  // Populate structs for threads
  for (size_t i=0 ; i < sysparm.nThreads ; i++) {
    threadarg[i].pSysparm      = pSysparm;
    threadarg[i].pApertureData = pData;
    threadarg[i].mask          = mask;
    threadarg[i].iOutputBegin  = 0+i*(nPositions/sysparm.nThreads);
    threadarg[i].iOutputEnd    =
      (nPositions/sysparm.nThreads)+i*(nPositions/sysparm.nThreads);
    threadarg[i].nPositions    =
      nPositions / sysparm.nThreads;
    // Output and indexing
    threadarg[i].pField        = _data;
    threadarg[i].iSampleStart  = sStart;
    threadarg[i].nData         = nSamples;

    // Input
    threadarg[i].pPos          = pPos;

    threadarg[i].threadId       = i;
    threadarg[i].cpuId          = i % nproc;

    threadarg[i].pGL            = &uv;
    if (i == 0) {
      threadarg[i].pQueue = progressQueue;
    }
    if (i == (sysparm.nThreads-1)) {
      threadarg[i].iOutputEnd   = nPositions;
    }
  }

  for (size_t i=0 ; i< sysparm.nThreads ; i++) {
# if defined(HAVE_PTHREAD_H)
    CallErr(pthread_create,
            (&fnm::threads[i],
             &fnm::attr,
             &CalcPwFnmThreadFunc<T, A>,
             &threadarg[i]));
# elif defined(_WIN32)
    threads[i] =
      (HANDLE)_beginthreadex(NULL, 0,
                             &CalcPwFnmThreadFunc<T, A>,
                             &threadarg[i], 0, &threadID);
# endif
  }

  // Monitor progress
  float fProgress = 0.0f;

# ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning(disable:4127)
# endif
  while (true) {
    fProgress = progressQueue->pop();
    if (fProgress == 100.0f) {
      break;
    }
# ifdef USE_PROGRESS_BAR
    if (pBar) {
      pBar->show(fProgress);
    }
    if (g_bExit) {
      break;
    }
# endif
# ifdef __GNUC__
    sched_yield();
# elif defined(_WIN32)
    SwitchToThread();
# endif
  }
# ifdef _MSC_VER
#  pragma warning(pop)
# endif

  for (size_t i = 0; i < sysparm.nThreads; i++) {
# if defined(HAVE_PTHREAD_H)
    CallErr(pthread_join, (threads[i], NULL));
# elif defined(_WIN32)
    if (threads[i]) {
      WaitForSingleObject(threads[i], INFINITE);
    }
# endif
  }

  g_bExit = false;

  delete progressQueue;

  ProfilerStop();

# if defined(HAVE_PTHREAD_H)
  CallErr(pthread_attr_destroy, (&fnm::attr));
# endif
#endif

  // Shared stuff

  // Assign output pointers
  *pOutData = static_cast<T*>(_data);
  *pNSamples = nSamples;

  return tStart;
}


template <class T,  template <typename> class A>
#if defined(HAVE_PTHREAD_H)
void* CalcPwFnmThreadFunc(void* ptarg)
#else
unsigned int __stdcall CalcPwFnmThreadFunc(void *ptarg)
#endif
{
  fnm::PwFnmThreadArg<T>* pThreadArg = reinterpret_cast<fnm::PwFnmThreadArg<T>*>(ptarg);

  const fnm::sysparm_t<T>* pSysparm = pThreadArg->pSysparm;

  const ApertureData<T>* pApertureData = pThreadArg->pApertureData;
  const int mask                       = pThreadArg->mask;
  const GLQuad2D<T>* pGL               = pThreadArg->pGL;
  const T* pos                         = pThreadArg->pPos;
  T* odata                             = pThreadArg->pField;
  size_t nData                         = pThreadArg->nData;
  int iSampleSignalStart               = pThreadArg->iSampleStart;

#ifdef HAVE_THREAD
  setcpuid(pThreadArg->cpuId);
#endif
  const size_t nPositions          = pThreadArg->nPositions;
  sps::queue<float>* pQueue        = pThreadArg->pQueue;

  const T amplitude = T(1.0);

  // Initial time
  double tProfStart = 0.0;

  if (pThreadArg->threadId == 0) {
    tProfStart = sps::profiler::time();
  }

  size_t nCallbackPeriod = 1;
  double duration = 0.0;

  for (size_t iPoint = pThreadArg->iOutputBegin ;
       iPoint < pThreadArg->iOutputEnd ; iPoint++) {
    // Why are we calling this function in sofus_calc_threaded
    // ComputeBoundingTimesCompact(pSysparm, box, &pos[iPoint*3], 1, delays, apodizations, nElements, sStart, &signal);


    ALIGN32_BEGIN sps::point_t<T> point ALIGN32_END;

    point[0] = pos[iPoint*3 + 0];
    point[1] = pos[iPoint*3 + 1];
    point[2] = pos[iPoint*3 + 2];

    FnmResponse<T, A>(pSysparm, pApertureData,
                      pGL,
                      // TODO(JMH): Add here SIMD abcissas
                      amplitude, point,
                      iSampleSignalStart,
                      nData,
                      &odata[iPoint*nData],
                      mask);

    if (pThreadArg->threadId==0) {
      // Check performance after 9 points
      if (iPoint == 9) {
        duration = sps::profiler::time() - tProfStart;
        nCallbackPeriod = (size_t) (10.0 / duration);
        nCallbackPeriod = std::max<size_t>(nCallbackPeriod, 1);
      }

      if ( (iPoint > 9) && (iPoint % nCallbackPeriod == 0) ) {
        float val = 100.0f * float(iPoint) / nPositions;
        pQueue->push(val);
# ifdef __GNUC__
        sched_yield();
# elif defined(_WIN32)
        SwitchToThread();
# endif
      }
    }
    if (g_bExit) {
      break;
    }
  }
  debug_print("Thread done\n");

  // Will always happen, even for nPoints < 9
  if (pThreadArg->threadId == 0) {
    pQueue->push(100.0f);
  }
#if HAVE_PTHREAD_H
  pthread_exit(NULL);
#else
  return 0;
#endif
}

template float
CalcPwFnmThreaded<float, ToneBurst>(
  const sysparm_t<float>* pSysparm,
  const ApertureData<float>* pData,
  const float* pPos, const size_t nPositions,
  float** pOutData, size_t* pNSamples, int mask,
  sps::ProgressBarInterface* pBar);

template float
CalcPwFnmThreaded<float, HanningWeightedPulse>(
  const sysparm_t<float>* pSysparm,
  const ApertureData<float>* pData,
  const float* pPos, const size_t nPositions,
  float** pOutData, size_t* pNSamples, int mask,
  sps::ProgressBarInterface* pBar);

template float
CalcPwFnmThreaded<float, Identity>(
  const sysparm_t<float>* pSysparm,
  const ApertureData<float>* pData,
  const float* pPos, const size_t nPositions,
  float** pOutData, size_t* pNSamples, int mask,
  sps::ProgressBarInterface* pBar);

#if FNM_DOUBLE_SUPPORT
template double
CalcPwFnmThreaded<double, ToneBurst>(
  const sysparm_t<double>* pSysparm,
  const ApertureData<double>* pData,
  const double* pPos, const size_t nPositions,
  double** pOutData, size_t* pNSamples,
  int mask,
  sps::ProgressBarInterface* pBar);

template double
CalcPwFnmThreaded<double, HanningWeightedPulse>(
  const sysparm_t<double>* pSysparm,
  const ApertureData<double>* pData,
  const double* pPos, const size_t nPositions,
  double** pOutData, size_t* pNSamples,
  int mask,
  sps::ProgressBarInterface* pBar);

template double
CalcPwFnmThreaded<double, Identity>(
  const sysparm_t<double>* pSysparm,
  const ApertureData<double>* pData,
  const double* pPos, const size_t nPositions,
  double** pOutData, size_t* pNSamples,
  int mask,
  sps::ProgressBarInterface* pBar);
#endif

}



/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
