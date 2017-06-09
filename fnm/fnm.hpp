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

/** @defgroup fnm Fast Near-field Method
 *  @brief This module is used for computing pressures and transients using
 *         the Fast Nearfield Method (FNM)
 *
 *  Continuous Wave (CW) pressure fields and gradient are computed in
 *  frequency domain using an expression, which is equivalent to the
 *  Fourier transform of the Rayleigh integral. The procedure is in
 *  the literature referred to as the Fast Near Field Method \cite
 *  mcgough2004rapid .
 *
 */

#include <sps/config.h>
#include <fnm/config.h>
#include <fnm/fnm_export.h>
#include <fnm/fnm_types.hpp>

#include <sps/smath_types.hpp> // Consider forward declaring element_t

#include <complex>

#include <cstdarg>

#ifdef _MSC_VER
# pragma warning(push)
# pragma warning(disable:4290)
#endif

//! Fast Nearfield Method interfaces and implementations
namespace fnm {
  template<class T>
  class ApertureData;
}

namespace std {
  class runtime_error;
}

namespace sps {
#ifdef USE_PROGRESS_BAR
  class ProgressBarInterface;
#endif
}

#ifdef FNM_PULSED_WAVE
namespace sofus {
  template <class T>
  class AperturePulses;
}
#endif

/*!
 * @addtogroup fnm
 * @{
 */

namespace fnm {

#ifdef USE_PROGRESS_BAR
  // TODO: Define functions using ProgIF
  typedef sps::ProgressBarInterface ProgIF;
#else
  typedef void ProgIF;
#endif

  /*! \brief Aperture class
   *
   * @tparam T floating point type
   *
   * A class representing an aperture. Consider using (backslash)nosubgrouping
   */
  template <class T>
  class FNM_EXPORT Aperture {

  public:

    /**
     * Static constructor for empty array
     *
     * @param obj
     *
     * @return
     */
    static int ArrayCreate(Aperture<T> **obj);

    /**
     * Static constructor for elevation focused linear array
     *
     * @param obj
     * @param nElements
     * @param width
     * @param kerf
     * @param height
     * @param nSubH      Number of sub-elements used for elevation
     * @param focus      Elevation focus depth
     *
     * @return
     */
    static int FocusedLinearArrayCreate(Aperture<T> **obj, const size_t nElements,
                                        const T width, const T kerf, const T height,
                                        const size_t nSubH, const T focus);

    /**
     * Create Matrix array - all elements are elements (no sub-elements)
     *
     * @param obj
     * @param nRows
     * @param nCols
     * @param rowWidth
     * @param rowKerf
     * @param colWidth
     * @param colKerf
     *
     * @return
     */
    static int MatrixArrayCreate(Aperture<T> **obj,
                                 const size_t nRows,
                                 const size_t nCols,
                                 const T rowWidth,
                                 const T rowKerf,
                                 const T colWidth,
                                 const T colKerf);

    /**
     * Static constructor for elevation focused convex array
     *
     * @param obj
     * @param nElements
     * @param width
     * @param kerf
     * @param height
     * @param radius
     * @param nSubH
     * @param focus
     *
     * @return
     */
    static int FocusedConvexArrayCreate(fnm::Aperture<T> **obj, const size_t nElements,
                                        const T width, const T kerf, const T height,
                                        const T radius, const size_t nSubH, const T focus);

    /**
     * @brief Constructor.
     *
     *  The constructor can be used to create arbitrary transducers
     *  using @ref ElementsSet or @ref SubElementsSet after
     *  construction.
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

    /** @name Accessors to read-only attributes
     *
     */
    ///@{
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

    /**
     * Get extent of aperture
     *
     * @param coordinates data
     * @param nDim    3 dimensions
     * @param nLimits Min and maximum
     */
    void ExtentGet(T** coordinates, size_t* nDim, size_t* nLimits) const;

    /**
     * Get area of aperture
     *
     * @return
     */
    T AreaGet() const;

    /**
     * Get phases of elements
     *
     * @param data
     * @param nData
     */
    void PhasesGet(T** data, size_t* nData) const;

    /**
     * Get phases of elements
     *
     * @param data
     * @param nData
     */
    void SubPhasesGet(T** data, size_t* nData) const;

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

    ///@}

    /** @name Read-only references for external use
     *
     */
    ///@{
    /**
     * Get reference to (sub)-elements
     *
     * @param nElements
     * @param nSubElements
     * @param elements
     *
     * @return
     */
    int ElementsRefGet(size_t* nElements, size_t* nSubElements,
                       const sps::element_t<T>**& elements) const;

    /**
     * Get reference to apodizations
     *
     * @param nElements
     * @param apodizations
     *
     * @return
     */
    int ApodizationsRefGet(size_t* nElements, const T*& apodizations) const;

    ///@}

    /** @name Accessor (getters) and Mutators (setters)
     *
     */
    ///@{
    /**
     * Get delays of elements
     *
     * @param data
     * @param nData
     */
    void DelaysGet(T** data, size_t* nData) const;

    /**
     * Set delays of elements
     *
     * @param data
     * @param nData
     */
    void DelaysSet(const T* data, const size_t nData);

    /**
     * Enable attenuation
     *
     * @param iEnabled
     */
    void AttenuationEnabledSet(const bool& iEnabled);

    /**
     * Is attenuation enabled
     *
     *
     * @return
     */
    const bool& AttenuationEnabledGet() const;


    /**
     * Set attenuation alpha parameter
     *
     * @param value
     */
    void AlphaSet(const T& value);

    /**
     * Get attenuation alpha parameter
     *
     * @return
     */
    const T& AlphaGet() const;

    /**
     * Set beta attenuation value.
     *
     * Signals are attenuated according to:
     *
     * \f$
     * h(x,y,z) \mapsto \exp(-(f-f_0)*\beta\, d\,) h(x,y,z),
     * \f$
     *
     * where \f$ d\f$ is the distance. Note the unit is Neper and
     * not dB. Multiply an attenuation in dB with \ref Aperture<T>::Neper_dB
     *
     * @param value [Neper/(Hz m)]
     */
    void BetaSet(const T& value);

    /**
     * Get beta attenuation value
     *
     * @return
     */
    const T& BetaGet() const;

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
     * Get center focus point
     *
     * @param oFocus
     */
    void CenterFocusGet(T oFocus[3]) const;

    /**
     * Set center focus point
     *
     * @param iFocus
     */
    void CenterFocusSet(const T iFocus[3]);

    /**
     * Get focus type used
     *
     */
    int FocusingTypeGet() const;

    /**
     * Set focus type used
     *
     * @param iFocusingType
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
    void F0Set(const T& f0);

    /**
     * Get width of pulse [s]
     *
     * @return
     */
    const T& WGet() const;

    /**
     * Set width of pulse [s]
     *
     * @param w
     */
    void WSet(const T& w);

    /**
     * Get system parameters (TODO: Avoid grouping)
     *
     *
     * @return
     */
    const sysparm_t<T> SysParmGet() const;

    /**
     * Set the system parameters
     *
     * @param arg
     */
    void SysParmSet(const sysparm_t<T> *arg);

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
    void CSet(const T& c);

    /**
     * Get element definitions
     *
     * @param out
     * @param nElements
     * @param nParams
     */
    void ElementsGet(T** out, size_t* nElements,
                     size_t* nParams) const;

    void ElementsSet(const T* pos, const size_t nPositions, const size_t nDim);// throw (std::runtime_error);

    /*! \fn bool Aperture::ElementsSet(const T* pos, const size_t nPositions, const size_t nDim)
     *  \brief Set element positions
     *  \param pos  Input data T[nElements][8]
     *  \param nPositions
     *  \param nDim Must equal 8
     *  \exception std::runtime_error nDim != 8
     *  \return true on success
     */

    /**
     * Get sub-element definitions
     *
     * @param out
     * @param nElements
     * @param nSubElements
     * @param nParams
     */
    void SubElementsGet(T** out, size_t* nElements,
                        size_t* nSubElements, size_t* nParams) const;

    /**
     * Set sub-element definitions
     *
     * @param pos
     * @param nElements
     * @param nSubElementsPerElement
     * @param nDim
     *
     * @return true on success
     */
    bool SubElementsSet(const T* pos, const size_t nElements,
                        const size_t nSubElementsPerElement, const size_t nDim);

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

    ///@}

#if FNM_PULSED_WAVE
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
    void FsSet(const T& f0);

    /**
     * Get normalization state
     *
     *
     * @return
     */
    const bool& NormalizeGet() const;

    /**
     * Enable normalization of convolutions
     *
     * @param normalize
     */
    void NormalizeSet(const bool& normalize);

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
    void NDivWSet(const size_t& nDivW);

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
    void NDivHSet(const size_t& nDivH);

#ifdef USE_PROGRESS_BAR
    /**
     * Set the progress bar used for notification
     *
     * @param pbar
     */
    void ProgressBarSet(sps::ProgressBarInterface* pbar);
#endif

#if !(defined(_MSC_VER) && _MSC_VER < 1900)
# ifndef _SWIG_WIN32

#  ifndef SWIG_VERSION
    int RwFloatParamSet(int fsel, const T* f, size_t nDim, va_list args);
    int RwFloatParamGet(int fsel, T** f, size_t nDim, va_list args);
#  endif
    int RwFloatParamSet(int fsel, const T* f, size_t nDim, ...);
    int RwFloatParamGet(int fsel, T** f, size_t nDim, ...);

    void RwBooleanParamSet0D(int fsel, const bool& value);
    void RwSizeTParamSet0D(int fsel, const size_t& value);
    void RwIntegerParamSet0D(int fsel, const int& value);

    int ParameterInfoGet(const RwParamType& param, ScalarType* type, size_t* nDims);
# endif
#endif

    /** @name Static variables
     *
     */
    ///@{

    static const size_t nVerticesPerElement = 4; /**< Number of corners for an element */

    static const size_t nElementPosParameters = 8; /**< Number of parameters for an element */

#if FNM_PULSED_WAVE
    static T fs;           ///< Sample frequency
    static bool normalize; ///< Normalize responses
#endif
    static size_t nthreads; ///< Number of threads

#if (defined(SWIG_VERSION) && (SWIG_VERSION > 0x030000)) && defined(__GNUG__)
    // Supported by SWIG3.0
    static constexpr T Neper_dB = 0.11512925464970231;
    static constexpr T dB_Neper = 8.685889638065035;
#else
    static const T Neper_dB; /**< Nepers per dB, 1 / ( 20 log(e)) */
    static const T dB_Neper; /**< dB per Neper, 20 log(e) */
#endif

    ///@}


    /**
     * Update phases after setting focus. This function is used by all
     * methods for adjusting phases for focusing.
     *
     */
    void FocusUpdate();

    /**
     * Reference implementation of fnm::FocusUpdate
     *
     */
    void FocusUpdateRef();

#ifdef FNM_PULSED_WAVE
    /**
     * Experimental attempt for hanning-weighted pulsed wave
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
    T CalcFdTransientRef(const T* pos, const size_t nPositions, const size_t nDim,
                         T** odata, size_t* nSignals, size_t* nSamples);
#endif

    /** @name FNM Calculation functions
     *
     */
    ///@{

    /**
     * Compute CW response at multiple positions. Reference
     * implementation. See @ref fnm::CalcCwFieldRef.
     *
     * @param pos
     * @param nPositions
     * @param nDim
     * @param odata
     * @param nOutPositions
     *
     * @return
     */
    int CalcCwFieldRef(const T* pos, const size_t nPositions, const size_t nDim,
                       std::complex<T>** odata, size_t* nOutPositions);


    /**
     * Compute CW response at multiple positions. Number of integrals is
     * reduced to four
     *
     * @param pos
     * @param nPositions
     * @param nDim
     * @param odata
     * @param nOutPositions
     *
     * @return
     */
    int CalcCwFieldFourRef(const T* pos, const size_t nPositions, const size_t nDim,
                           std::complex<T>** odata, size_t* nOutPositions);

    //

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
     *
     * @return
     */
    int CalcCwFast(const T* pos, const size_t nPositions, const size_t nDim,
                   std::complex<T>** odata, size_t* nOutPositions);

#ifndef FNM_DOUBLE_SUPPORT
    /**
     * Same as CalcCwField, but uses SIMD.
     *
     * @param pos
     * @param nPositions
     * @param nDim
     * @param odata
     * @param nOutPositions
     *
     * @return
     */
    int CalcCwFieldNaiveFast(const T* pos, const size_t nPositions, const size_t nDim,
                             std::complex<T>** odata, size_t* nOutPositions);
#endif

    /**
     * Compute CW response at multiple positions. The range of
     * integration is naive so it is accurate when projections lie
     * inside an element, but requires a huge amount of abcissas to
     * get a usuable result, when projections lie outside an element.
     *
     * @param pos
     * @param nPositions
     * @param nDim
     * @param odata
     * @param nOutPositions
     */
    int CalcCwFieldNaive(const T* pos, const size_t nPositions, const size_t nDim,
                         std::complex<T>** odata, size_t* nOutPositions);

    /*** @} */

#ifdef FNM_PULSED_WAVE


    /** @name SOFUS accellerated calculation functions
     *
     */
    ///@{

    /**
     *
     *
     * @param other
     * @param iFocus
     * @param odata
     * @param nSignals
     * @param nSamples
     * @param sizeAndOffsets
     * @param nFilters
     * @param nTwo
     * @param data
     * @param nData
     *
     * @return
     */
    T CalcMatchedFilter(const Aperture<T>* other,
                        const T iFocus[3],
                        T** odata, size_t* nSignals, size_t* nSamples,
                        int** sizeAndOffsets, size_t* nFilters, size_t* nTwo,
                        T** data, size_t* nData);

    /**
     *
     *
     * @param other
     * @param iFocus
     * @param tStart
     * @param iData
     * @param nChannels
     * @param nSamples
     * @param data
     * @param nData
     *
     * @return
     */
    T CalcSmfApply(const Aperture<T>* other,
                   const T iFocus[3],
                   const T tStart,
                   const T* iData, const size_t nChannels, const size_t nSamples,
                   T** data, size_t* nData);

    /**
     *
     *
     * @param other
     * @param pos
     * @param nPositions
     * @param nDim
     * @param data
     * @param nData
     * @param odata
     * @param nSignals
     * @param nSamples
     *
     * @return
     */
    T CalcPwEcho(const Aperture<T>* other,
                 const T* pos, const size_t nPositions, const size_t nDim,
                 const T* data, const size_t nData,
                 T** odata, size_t* nSignals, size_t* nSamples);

    /**
     * Compute pulsed-wave field
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

    /**
     * Reference implementation of Aperture::CalcPwField
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
    T CalcPwFieldRef(const T* pos, const size_t nPositions, const size_t nDim,
                     T** odata, size_t* nSignals, size_t* nSamples);


# ifndef FNM_DOUBLE_SUPPORT

    /**
     * Threaded version of Aperture::CalcPwField. \bug: Sometimes to little memory is allocated
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
    T CalcPwFieldThreaded(const T* pos, const size_t nPositions, const size_t nDim,
                          T** odata, size_t* nSignals, size_t* nSamples);
# endif
#endif
    /*** @} */

    /** @name Developers corner (do not touch)
     *
     */
    ///@{
    void InPlaceOP(T* ioData, size_t nIOdata);
    ///@}


#ifndef FNM_DOUBLE_SUPPORT
    /** @name Angular spectrum approach (work in progress)
     *
     */
    ///@{

    /**
     *
     *
     * @param y0
     * @param nx
     * @param ny
     * @param dx
     * @param dy
     * @param Nx
     * @param Ny
     * @param p1
     * @param onx
     * @param ony
     * @param onz
     *
     * @return
     */
    int CalcAsa(const T* y0,const size_t nx,const size_t ny,
                const T dx,const T dy,
                const size_t Nx,const size_t Ny,
                std::complex<T>** p1, size_t *onx, size_t *ony, size_t *onz);
    ///@}

#endif

  private:
    /**
     * Get default system parameters. The system parameters are shared among
     * all apertures, unless system paramters are changed using
     * Aperture::SysParmGet and Aperture::SysParmSet for a specific array.
     *
     *
     * @return
     */
    static sysparm_t<T>* DefaultSysParmGet();

    /**
     * Internal function for getting elements or sub-elements. The
     * function is used by Aperture::ElementsGet and
     * Aperture::SubElementsGet.
     *
     * @param sub 0 for elements, 1 for subelements
     * @param out
     * @param nElements
     * @param nSubElements
     * @param nParams
     *
     * @return
     */
    int MultiElementsGet(int sub,
                         T** out, size_t* nElements,
                         size_t* nSubElements, size_t* nParams) const;

    //    Okay with inline mallocs in Info array
    //    void FocusGetAllocate(T** data, size_t* nData) const;

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
                    std::complex<T>** odata, size_t* nOutPositions) const;

    /**
     * Reference implementation of Aperture::CalcCwFocus
     *
     * @param pos
     * @param nPositions
     * @param nDim
     * @param odata
     * @param nOutPositions
     *
     * @return
     */
    int CalcCwFocusRef(const T* pos, const size_t nPositions, const size_t nDim,
                       std::complex<T>** odata, size_t* nOutPositions);

  private:
    /// PIMPL design (data)

    ///< Sysparm
    sysparm_t<T>* m_sysparm;

    ///< Data
    ApertureData<T>* m_data;

    // TODO: Use typedef ProgIF
#ifdef USE_PROGRESS_BAR
    sps::ProgressBarInterface* m_pbar;
#else
    void* m_pbar;
#endif

#ifdef FNM_PULSED_WAVE
    sofus::AperturePulses<T>* m_pulses;
#endif
  };
}
/*! @} End of FNM Group */

#ifdef _MSC_VER
#pragma warning(pop)
#endif

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
