/**
 * @file   fnm.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Tue Jun  7 23:52:37 2016
 *
 * @brief  Contains Aperture class with methods using the
 *         fast nearfield method (FNM)
 *
 */
#pragma once

// TODO: Introduce header with interface structs

/** @defgroup FNM Fast Near-field Method
 * @brief This module is used for computing responses using the fast nearfield method (FNM)
 *
 */

#include <fnm/fnm_export.h>
#include <fnm/config.h>
#include <fnm/fnm_types.hpp>

#ifdef FNM_PULSED_WAVE
#include <sofus/sofus_pulses.hpp>
#endif

#include <complex>

//! Fast Nearfield Method interfaces and implementations
namespace fnm {
  /// Forward-declare ApertureData
  template<class T>
  class ApertureData;
}

namespace std {
  class runtime_error;
}

namespace sps {
  class ProgressBarInterface;
}

/*!
 * @addtogroup FNM
 * @{
 */

namespace fnm {

  /*! \brief Aperture class
   *
   *
   * A class representing an aperture
   */
  template <class T>
  class FNM_EXPORT Aperture {

  public:
    /**
     * Constructor
     *
     */
    Aperture();

    /**
     * Constructor for linear array
     *
     * @param nElements
     * @param width
     * @param kerf
     * @param height
     *
     * @return
     */
    Aperture(const size_t nElements, const T width, const T kerf, const T height);

    /**
     * Destructor
     *
     *
     * @return
     */
    ~Aperture();

    //a{
    /** Read-only attributes */

    /**
     * Get element count
     *
     * @return
     */
    const size_t& NElementsGet() const;

    /**
     * Get sub-element count
     *
     * @return
     */
    const size_t& NSubElementsGet() const;

    void ExtentGet(T** coordinates, size_t* nDim, size_t* nLimits) const;

    T AreaGet() const;

    /**
     * Get phases of elements
     *
     * @param data
     * @param nData
     */
    void PhasesGet(T** data, size_t* nData);

    /**
     * Get delays of elements
     *
     * @param data
     * @param nData
     */
    void DelaysGet(T** data, size_t* nData);

    /**
     * Get rectangles for display
     *
     * @param out
     * @param nElements
     * @param nSubElements
     * @param nParams
     */
    void RectanglesGet(T** out, size_t* nElements,
                       size_t* nSubElements, size_t* nParams) const;

    //a}

    //a{
    /** Read-write attributes */

    /**
     * Get element position data
     *
     * @param out
     * @param nElements
     * @param nParams
     */
    void PositionsGet(T** out, size_t* nElements, size_t* nParams) const;

    /**
     * Set element position data (may differ from center-most sub-element)
     *
     * @param pos
     * @param nPositions
     * @param nDim
     */
    void PositionsSet(const T* pos, const size_t nPositions, const size_t nDim);

    /**
     * Get focus point or virtual source
     *
     * @param oFocus
     */
    void FocusGet(T oFocus[3]) const;

    /**
     * Set focus point or virtual source
     *
     * @param iFocus
     */
    void FocusSet(const T iFocus[3]);

    /**
     * Get focus type used
     *
     * @param oFocus
     */
    int FocusingTypeGet() const;

    /**
     * Set focus type used
     *
     * @param iFocus
     */
    void FocusingTypeSet(const int iFocusingType);

    /**
     * Get number of threads
     *
     *
     * @return # of threads
     */
    const size_t& NThreadsGet() const;

    /**
     * Set number of threads
     *
     * @param nThreads
     */
    void NThreadsSet(const size_t &nThreads);

    /**
     * Get center frequency
     *
     * @return
     */
    const T& F0Get() const;

    /**
     * Set center frequency
     *
     * @param f0
     */
    void F0Set(const T f0);

    void SysParmSet(const sysparm_t<T> *arg);

#ifdef FNM_PULSED_WAVE
    /**
     * Get sampling frequency
     *
     * @return
     */
    const T& FsGet() const;

    /**
     * Set sampling frequency
     *
     * @param f0
     */
    void FsSet(const T f0);

    /**
     * Get excitation (reference to or view of)
     *
     * @param data
     * @param nData
     */
    void ExcitationGet(T** data, size_t* nData) const;

    /**
     * Set excitation
     *
     * @param data
     * @param nData
     */
    void ExcitationSet(const T* data,
                       const size_t nData);

    /**
     * Get impulse (reference to or view of)
     *
     * @param data
     * @param nData
     */
    void ImpulseGet(T** data, size_t* nData) const;

    /**
     * Set impulse
     *
     * @param data
     * @param nData
     */
    void ImpulseSet(const T* data,
                    const size_t nData);

#endif

    /**
     * Get speed of sound
     *
     *
     * @return
     */
    const T& CGet()  const;

    /**
     * Set speed of sound
     *
     * @param c
     */
    void CSet(const T c);

    /**
     * Get element definitions
     *
     * @param out
     * @param nElements
     * @param nParams
     */
    void ElementsGet(T** out, size_t* nElements,
                     size_t* nParams) const;

    bool ElementsSet(const T* pos, const size_t nPositions, const size_t nDim) throw (std::runtime_error);

    /*! \fn bool Aperture::ElementsSet(const T* pos, const size_t nPositions, const size_t nDim)
     *  \brief Set element positions
     *  \param pos  Input data T[nElements][8]
     *  \param nPositions
     *  \param nDim Must equal 8
     *  \exception std::runtime_error nDim != 8
     *  \return true on success
     */

    /**
     * Get sub-element positions
     *
     * @param out
     * @param nElements
     * @param nSubElements
     * @param nParams
     */
    void SubElementsGet(T** out, size_t* nElements,
                        size_t* nSubElements, size_t* nParams) const;

    /**
     * Set sub-element positions
     *
     * @param pos
     * @param nElements
     * @param nSubElementsPerElement
     * @param nDim
     *
     * @return true on success
     */
    bool SubElementsSet(const T* pos, const size_t nElements, const size_t nSubElementsPerElement, const size_t nDim);

    /**
     * Get number of width abcissas
     *
     *
     * @return # of abcissa in the width dimension
     */
    const size_t& NDivWGet() const;

    /**
     * Set number of width abcissas
     *
     * @param nDivW
     */
    void NDivWSet(const size_t nDivW);

    /**
     * Get number of height abcissas
     *
     *
     * @return
     */
    const size_t& NDivHGet() const;

    /**
     * Set number of height abicissas
     *
     * @param nDivH
     */
    void NDivHSet(const size_t nDivH);

    /**
     * Get apodization
     *
     * @param data
     * @param nData
     */
    void ApodizationGet(T** data, size_t* nData) const;

    /**
     * Set apodization
     *
     * @param data
     * @param nData
     */
    void ApodizationSet(const T* data, const size_t nData);

    void ProgressBarSet(sps::ProgressBarInterface* pbar);
    //a}

    /**
     * Update phases after setting focus. This function is used by all
     * methods for adjusting phases for focusing.
     *
     */
    void FocusUpdate();

    //@{ Calculation functions

#ifdef FNM_PULSED_WAVE
    /**
     *
     *
     * @param pos
     * @param nPositions
     * @param nDim
     * @param odata
     * @param nSignals
     * @param nSamples
     *
     * @return
     */
    T CalcPwField(const T* pos, const size_t nPositions, const size_t nDim,
                  T** odata, size_t* nSignals, size_t* nSamples);
#endif

    /**
     * Compute CW response at multiple positions. Reference
     * implementation. The ranges of integration are reduced to give
     * the most accurate result with less abcissas.
     *
     * @param pos
     * @param nPositions
     * @param nDim
     * @param odata
     * @param nOutPositions
     */
    int CalcCwFieldRef(const T* pos, const size_t nPositions, const size_t nDim,
                       std::complex<T>** odata, size_t* nOutPositions);

#ifndef FNM_DOUBLE_SUPPORT
    /**
     * Compute CW response at multiple positions (uses SIMD). The
     * ranges of integration are reduced to give the most accurate
     * result with less abcissas.
     *
     * @param pos
     * @param nPositions
     * @param nDim
     * @param odata
     * @param nOutPositions
     */
    int CalcCwFast(const T* pos, const size_t nPositions, const size_t nDim,
                   std::complex<T>** odata, size_t* nOutPositions);
#endif

    /**
     * Compute CW response at multiple positions. The range of
     * integration is naive so it is accurate when projections lie
     * inside an element, but requires a huge amount of abcissas to
     * get a usuable result, when projections lie outside an element.
     *
     * TODO: Remove when CalcCwFieldRef is fixed
     *
     * @param pos
     * @param nPositions
     * @param nDim
     * @param odata
     * @param nOutPositions
     */
    int CalcCwField(const T* pos, const size_t nPositions, const size_t nDim,
                    std::complex<T>** odata, size_t* nOutPositions);

#ifndef FNM_DOUBLE_SUPPORT
    /**
     * Same as CalcCwField, but uses SIMD.
     *
     * TODO: Remove when CalcCwFieldRef is fixed
     *
     * @param pos
     * @param nPositions
     * @param nDim
     * @param odata
     * @param nOutPositions
     *
     * @return
     */
    int CalcCwField2(const T* pos, const size_t nPositions, const size_t nDim,
                     std::complex<T>** odata, size_t* nOutPositions);
#endif

    //@}

    //a{
    /** Experimental functions */

#ifndef FNM_DOUBLE_SUPPORT
    int CalcAsa(const T* y0,const size_t nx,const size_t ny,
                const T dx,const T dy,
                const size_t Nx,const size_t Ny,
                std::complex<T>** p1, size_t *onx, size_t *ony, size_t *onz);
#endif

    //a}

    //a{
    /** Static variables */
    static const size_t nVerticesPerElement = 4; /**< Number of corners for an element */

    static const size_t nElementPosParameters = 8; /**< Number of parameters for an element */

    /// System parameters (We could make a long list of functions friends of Aperture)
    static sysparm_t<T> _sysparm;

#ifdef FNM_PULSED_WAVE
    static T fs;
#endif
    //a}
  private:
    void ManagedAllocation(std::complex<T>** outTest, size_t* nOutTest);


    /**
     * The positions must all equal the focus point and the number
     * should match the number of elements. By doing so, a set of
     * phases is computed for focusing.
     *
     * Used by \ref FocusUpdate
     *
     * @param pos
     * @param nPositions
     * @param nDim
     * @param odata
     * @param nOutPositions
     *
     * @return
     */
    int CalcCwFocus(const T* pos, const size_t nPositions, const size_t nDim,
                    std::complex<T>** odata, size_t* nOutPositions);

    /// Number of threads
    static size_t nthreads;

    /**
     * Initialize elements. This function must be called whenever
     * elements are changed.
     *
     */
    void initElements();

    // TODO: Move these to .cpp translation unit, must take pointer to Aperture<T>

#ifndef FNM_DOUBLE_SUPPORT
    /**
     * Thread function for computing harmonic response
     *
     * @param ptarg
     *
     * @return
     */
#if defined(HAVE_PTHREAD_H)
# ifdef HAVE_MQUEUE_H
    void* CalcThreadFunc(void* ptarg);
# endif
    // TODO: Invalid read of size 4
    void* CalcThreaded(void* ptarg);
#elif defined(_WIN32)
    unsigned int __stdcall CalcThreaded(void *ptarg);
#endif
#endif

    ///< Data
    ApertureData<T>* m_data;

    sps::ProgressBarInterface* m_pbar;
#ifdef FNM_PULSED_WAVE
    sofus::AperturePulses<T>* m_pulses;
#endif
  };
}
/*! @} End of FNM Group */


/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
