/**
 * @file   sofus_calc.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Wed Oct 19 21:41:07 2016
 *
 * @brief
 *
 *
 */
#include <sofus/sofus_export.h>
#include <sofus/config.h>
#include <sps/config.h>
#include <sps/smath.hpp>
#include <sps/queue.hpp>

#include <sofus/sofus_types.hpp>
#include <sofus/sofus_pulses.hpp>

/** @defgroup sofus SIMD Optimized Fast Ultrasound Simulator
 *  @brief This module is used for computing responses using time-domain methods
 *
 *  The methods here include some of the methods available using the well-known
 * <a href="http://field-ii.dk/"> Field II</a> by Jørgen Arendt Jensen
 *
 *
 */

// If used by FNM
#include <fnm/fnm_data.hpp>

namespace sofus {
// TODO: Consider using a forward declaration
template <class T>
using ApertureData = fnm::ApertureData<T>;
}

namespace sps {
#ifdef USE_PROGRESS_BAR
class ProgressBarInterface;
#endif
template <class T>
class msignal1D;
}

/*!
 * @addtogroup sofus
 * @{
 */


namespace sofus {

#ifdef HAVE_PTHREAD_H
extern pthread_t threads[N_MAX_THREADS];
extern pthread_attr_t attr;

class dummy {
  static int myStatic;
};

#endif

/**
 * Compute boundary for signals at points from an aperture confined to a box. This can be
 * used for allocating memory for a one-way response.
 *
 * @param sysparm
 * @param box A box encapsulating the aperture
 * @param points
 * @param nPoints
 * @param delays
 * @param apodizations
 * @param nElements
 * @param tStart
 * @param iStartSample
 * @param nSamples
 */
template<class T>
void ComputeBoxTimes(const sysparm_t<T>& sysparm,
                     const sps::bbox_t<T>& box,
                     const T* points, const size_t nPoints,
                     const T* delays, const T* apodizations,
                     const size_t nElements,
                     T* tStart,
                     int* iStartSample, size_t* nSamples);

#if 0
/**
 * Used for CalcPwEcho
 *
 * @param sysparm
 * @param scatter_box
 * @param xmt
 * @param rcv
 * @param startSample
 * @param nSamples
 */
template<class T>
void ComputeBoxBoxTimes(const sysparm_t<T>& sysparm,
                        const sps::bbox_t<T>& scatter_box,
                        const ApertureData<T>* xmt,
                        const ApertureData<T>* rcv,
                        int* startSample, size_t* nSamples);

/**
 * Used for CalcPwEcho
 *
 * @param sysparm
 * @param scatter_box
 * @param xmt
 * @param rcv
 * @param startSample
 * @param nSamples
 */
template<class T>
void ComputeBoxBoxWithoutTimes(const sysparm_t<T>& sysparm,
                               const sps::bbox_t<T>& scatter_box,
                               const ApertureData<T>* xmt,
                               const ApertureData<T>* rcv,
                               int* startSample, size_t* nSamples);

/**
 *
 *
 * @param data
 * @param sysparm
 * @param box
 * @param points
 * @param nPoints
 * @param delays
 * @param apodizations
 * @param nElements
 * @param tStart
 * @param nSamples
 */
template<class T>
void ComputeBoxTimesFarField(const ApertureData<T>* data,
                             const sysparm_t<T>& sysparm,
                             const sps::bbox_t<T>& box,
                             const T* points, const size_t nPoints,
                             const T* delays, const T* apodizations,
                             const size_t nElements,
                             T* tStart, size_t* nSamples);

/**
 *
 *
 * @param sysparm
 * @param pXmt
 * @param pRcv
 * @param box
 * @param point
 * @param tStart
 * @param nSamples
 */
template <class T>
void ComputeBoundingTimes(const sysparm_t<T>& sysparm,
                          const ApertureData<T>* pXmt,
                          const ApertureData<T>* pRcv,
                          const sps::bbox_t<T>& box,
                          const sps::point_t<T>* point,
                          T* tStart, size_t* nSamples);

/**
 * Compact version - works for a single scatter
 *
 * @param sysparm
 * @param box
 * @param points
 * @param nPoints
 * @param delays
 * @param apodizations
 * @param nElements
 * @param reference
 * @param signal
 */
template<class T>
void ComputeBoundingTimesCompact(const sysparm_t<T>* sysparm,
                                 const sps::bbox_t<T>& box,
                                 const T* points, const size_t nPoints,
                                 const T* delays,
                                 const T* apodizations,
                                 const size_t nElements,
                                 const int reference,
                                 sps::msignal1D<T>* signal);

/** @defgroup sofus_calc_functions SOFUS Calculation functions
 *  @ingroup sofus
 * @{
 */

template<class T>
T CalcPwBackField(const sysparm_t<T>& sysparm,
                  const fnm::ApertureData<T>* pApertureData,
                  const AperturePulses<T>* pPulses,
                  const T* pos, const size_t nPositions, const size_t nDim,
                  T** odata, size_t* nChannels, size_t* nSamples,
                  sps::ProgressBarInterface* pBar);


/**
 * CalcMatchedFilter
 *
 * Compute spatial matched filters
 *
 * @param sysparm    System parameters
 * @param data0      Data for transmit aperture (element positions, delays, apodizations)
 * @param data1      Data for receive aperture (element positions)
 * @param point      Location
 * @param odata      Output responses
 * @param nChannels  Number of channels
 * @param nSamples   Number of samples (length of longest response)
 * @param sizeAndOffsets Response length and offsets relative to earliest response
 * @param nFilters   Number of filters = number of receive channels
 * @param nTwo       Length and offsets (2 quantities per filter)
 * @param data       Time offsets (tStart for each filter [s])
 * @param nData      Number of offsets = number of receive channels
 *
 * @return
 */
template <class T>
T CalcMatchedFilter(const sofus::sysparm_t<T>& sysparm,
                    const ApertureData<T>* data0,
                    const ApertureData<T>* data1,
                    const sps::point_t<T>& point,
                    T** odata, size_t* nChannels, size_t* nSamples,
                    int** sizeAndOffsets, size_t* nFilters, size_t* nTwo,
                    T** data, size_t* nData);

/**
 * CalcSmfApply
 *
 * @param sysparm
 * @param data0
 * @param data1
 * @param point
 * @param tStart
 * @param iData
 * @param nChannels
 * @param nSamples
 * @param data
 * @param nData
 *
 * @return
 */
template <class T>
T CalcSmfApply(const sofus::sysparm_t<T>& sysparm,
               const ApertureData<T>* data0,
               const ApertureData<T>* data1,
               const sps::point_t<T>& point,
               const T tStart,
               const T* iData, const size_t nChannels, const size_t nSamples,
               T** data, size_t* nData);

/**
 *
 *
 * @param sysparm
 * @param data0
 * @param data1
 * @param pulses0
 * @param pulses1
 * @param pos
 * @param nPositions
 * @param nDim
 * @param amplitudes
 * @param nScatterers
 * @param odata
 * @param nChannels
 * @param nSamples
 *
 * @return
 */
template <class T>
T CalcPwEcho(const sofus::sysparm_t<T>& sysparm,
             const ApertureData<T>* data0,
             ApertureData<T>* data1,
             const AperturePulses<T>* pulses0,
             const AperturePulses<T>* pulses1,
             const T* pos, const size_t nPositions, const size_t nDim,
             const T* amplitudes, const size_t nScatterers,
             T** odata, size_t* nChannels, size_t* nSamples,
             sps::ProgressBarInterface* pBar);

template <class T>
T CalcScat(const sofus::sysparm_t<T>* pSysparm,
           const ApertureData<T>* pXmtData,
           const ApertureData<T>* pRcvData,
           const T* pos, const size_t nPositions, const size_t nDim,
           const T* amplitudes, const size_t nScatterers,
           T** odata, size_t* nLines, size_t* nSamples);

template <class T>
T CalcScatThreaded(const sofus::sysparm_t<T>* pSysparm,
                   const sofus::FocusLineList<T>* pXmtLines,
                   const ApertureData<T>* pXmtData,
                   const ApertureData<T>* pRcvData,
                   const AperturePulses<T>* pPulses0,
                   const AperturePulses<T>* pPulses1,
                   const T* pos, const size_t nPositions, const size_t nDim,
                   const T* amplitudes, const size_t nScatterers,
                   T** odata, size_t* nLines, size_t* nSamples);


/**
 *
 *
 * @param sysparm
 * @param data
 * @param pos
 * @param nPositions
 * @param nDim
 * @param odata
 * @param nSignals
 * @param nSamples
 * @param pBar
 *
 * @return
 */
#ifdef USE_PROGRESS_BAR
template <class T>
T CalcPwField(const sysparm_t<T>& sysparm,
              const ApertureData<T>* data,
              const T* pos, const size_t nPositions, const size_t nDim,
              T** odata, size_t* nSignals, size_t* nSamples,
              sps::ProgressBarInterface* pBar);
#else
template <class T>
T CalcPwField(const sysparm_t<T>& sysparm,
              const ApertureData<T>* data,
              const T* pos, const size_t nPositions, const size_t nDim,
              T** odata, size_t* nSignals, size_t* nSamples,
              void* pBar);
#endif

/**
 *
 *
 * @param sysparm
 * @param data
 * @param pulses
 * @param pos
 * @param nPositions
 * @param nDim
 * @param odata
 * @param nSignals
 * @param nSamples
 *
 * @return
 */
template <class T>
T CalcPwFieldRef(const sysparm_t<T>& sysparm,
                 const ApertureData<T>* data,
                 const AperturePulses<T>* pulses,
                 const T* pos, const size_t nPositions, const size_t nDim,
                 T** odata, size_t* nSignals, size_t* nSamples);

/**
 *
 *
 * @param sysparm
 * @param data
 * @param pulses
 * @param pos
 * @param nPositions
 * @param nDim
 * @param odata
 * @param nSignals
 * @param nSamples
 * @param pBar
 *
 * @return
 */
#ifdef USE_PROGRESS_BAR
template <class T>
T CalcPwFieldThreaded(const sysparm_t<T>& sysparm,
                      const ApertureData<T>* data,
                      const AperturePulses<T>* pulses,
                      const T* pos, const size_t nPositions, const size_t nDim,
                      T** odata, size_t* nSignals, size_t* nSamples,
                      sps::ProgressBarInterface* pBar);
#else
template <class T>
T CalcPwFieldThreaded(const sysparm_t<T>& sysparm,
                      const ApertureData<T>* data,
                      const AperturePulses<T>* pulses,
                      const T* pos, const size_t nPositions, const size_t nDim,
                      T** odata, size_t* nSignals, size_t* nSamples,
                      void* pBar);
#endif

/*** @} */ // end of group sofus_calc_functions

/**
 *
 *
 * @param ptarg
 *
 * @return
 */
template <class T>
void* CalcPwFieldThreadFunc(void* ptarg);

#endif

}  // namespace sofus

/*! @} End of sofus Group */


/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
