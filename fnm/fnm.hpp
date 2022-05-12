/**
 * @file   fnm.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Tue Jun  7 23:52:37 2016
 *
 * @brief  Contains Aperture class with methods using the
 *         fast nearfield method (FNM)
 *
 * Copyright 2017, Jens Munk Hansen
 *
 */
#pragma once

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

/** @defgroup fnm Fast Near-field Method
 *  @brief This module is used for computing pressures and transients using
 *         the Fast Nearfield Method (FNM)
 *
 *  Continuous Wave (CW) pressure fields and gradient are computed in
 *  frequency domain using an expression, which is equivalent to the
 *  Fourier transform of the Rayleigh integral. The procedure is in
 *  the literature referred to as the Fast Near Field Method \cite mcgough2004rapid .
 *
 */

#include <sps/config.h>
#include <fnm/config.h>
#include <fnm/fnm_export.h>

#include <complex>
#include <cstdarg>
#include <atomic>
#include <fnm/fnm_types.hpp>

#ifdef FNM_PULSED_WAVE
namespace sofus {
template <class T>
class AperturePulses;
}
# include <sofus/sofus_types.hpp>
# include <sofus/sofus_focus_line.hpp>
#endif

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

template <typename T>
struct element_rect_t;
}


/*!
 * @addtogroup fnm
 * @{
 */

namespace fnm {

#ifdef USE_PROGRESS_BAR
typedef sps::ProgressBarInterface ProgIF;
#else
typedef void ProgIF;
#endif

#ifdef WIN32
# define TEMPLATE_EXTERN
#else
# define TEMPLATE_EXTERN extern
#endif

/*! \brief Aperture class
 *
 * @tparam T floating point type
 *
 * A class representing an aperture. Consider using (backslash)nosubgrouping
 */
template <class T>
class Aperture {
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
   * @param[out] obj
   * @param[in] nElements
   * @param[in] width
   * @param[in] kerf
   * @param[in] height
   * @param[in] nSubH      Number of sub-elements used for elevation (default=1)
   * @param[in] focus      Elevation focus depth (0.0 means no elevation focus)
   *
   * @return
   */
  static int FocusedLinearArrayCreate(Aperture<T> **obj, const size_t nElements,
                                      const T width, const T kerf,
                                      const T height, const size_t nSubH,
                                      const T focus);

  /**
   * Create Matrix array - all elements are elements (no sub-elements)
   *
   * @param[out] obj
   * @param[in] nRows
   * @param[in] nCols
   * @param[in] rowWidth
   * @param[in] rowKerf
   * @param[in] colWidth
   * @param[in] colKerf
   *
   * @return
   */
  static int MatrixArrayCreate(Aperture<T> **obj,
                               const size_t nRows, const size_t nCols,
                               const T rowWidth, const T rowKerf,
                               const T colWidth, const T colKerf);

  /**
   * Static constructor for elevation focused convex array
   *
   * @param[out] obj
   * @param[in] nElements
   * @param[in] width
   * @param[in] kerf
   * @param[in] height
   * @param[in] radius
   * @param[in] nSubH
   * @param[in] focus
   *
   * @return
   */
  static int FocusedConvexArrayCreate(
    fnm::Aperture<T> **obj, const size_t nElements,
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
   * @param[in] nElements
   * @param[in] width
   * @param[in] kerf
   * @param[in] height
   *
   * @deprecated The user is encouraged to use Aperture::FocusedLinearArrayCreate instead
   *
   * @see Aperture::FocusedLinearArrayCreate
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

  void SubToElements();

  /**
   * Rotate the array around a point
   *
   * @param iEuler Euler angles using convention 'yxy'
   * @param iFocus Point to rotate about
   */
  void Rotate(const T iEuler[3], const T iFocus[3]);

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
   * @param[out] coordinates data
   * @param[out] nDim    3 dimensions
   * @param[out] nLimits Min and maximum
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
   * @param[out] data
   * @param[out] nData
   */
  void PhasesGet(T** data, size_t* nData) const;

  /**
   * Get phases of elements
   *
   * @param[out] data
   * @param[out] nData
   */
  void SubPhasesGet(T** data, size_t* nData) const;

  /**
   * Get rectangles for display
   *
   * @param[out] out
   * @param[out] nElements
   * @param[out] nSubElements
   * @param[out] nParams
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
   * @param[out] nElements
   * @param[out] nSubElements
   * @param[out] elements
   *
   * @return
   */
  int ElementsRefGet(size_t* nElements, size_t* nSubElements,
                     const sps::element_rect_t<T>**& elements) const;

  // Google requires all input parameters to be values or const references and
  // output parameters to be pointers. Google never allows non-const reference
  // parameters - except when required by convention, e.g. std::swap

  /**
   * Get reference to apodizations
   *
   * @param[out] nElements
   * @param[out] apodizations
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
   * @param[out] data
   * @param[out] nData
   */
  void DelaysGet(T** data, size_t* nData) const;

  /**
   * Set delays of elements
   *
   * @param[in] data
   * @param[in] nData
   */
  int DelaysSet(const T* data, const size_t nData);

  /**
   * Enable attenuation
   *
   * @param[in] iEnabled
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
   * @param[in] value
   */
  void AlphaSet(const T& value);

  /**
   * Get attenuation alpha parameter
   *
   * @return alpha
   */
  const T& AlphaGet() const;

  /**
   * Set beta attenuation value.
   *
   * Signals are attenuated according to:
   *
   * \f{align}{
   * h(x,y,z) &\mapsto \exp{(-f \beta \, d\,)} h(x,y,z)\nonumber\\
   * &\mapsto \exp{(-(f-f_0)\beta \, d\,)}\exp{(-\beta f0 \, d\,)} h(x,y,z)\nonumber\\
   * &\mapsto \exp{(-(f-f_0)\beta \, d\,)}\exp{(-\alpha\, d\,)} h(x,y,z),\nonumber
   * \f}
   *
   * where \f$ d\f$ is the distance. Note the unit is Neper and
   * not dB. Multiply an attenuation in dB with \ref Aperture<T>::Neper_dB
   *
   * Example:
   *
   * 0.54 dB / (cm MHz), f0 = 2 MHz
   *
   * beta = 0.54 * 100 / 1e6 * 0.11 Neper / dB = 6.22e-6 Neper/(Hz m)
   *
   * alpha = 2e6 * 0.54 * 100 / 1e6 * 0.11 Neper/ dB = 12.43 Neper/m
   *
   * @param[in] value [Neper/(Hz m)]
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
   * @param[out] out
   * @param[out] nElements
   * @param[out] nParams
   */
  void PositionsGet(T** out, size_t* nElements, size_t* nParams) const;

  /**
   * Set element position data (may differ from center-most sub-element)
   *
   * @param[in] pos
   * @param[in] nPositions
   * @param[in] nDim
   */
  int PositionsSet(const T* pos, const size_t nPositions, const size_t nDim);

  /**
   * Get focus point or virtual source
   *
   * @param[out] oFocus
   */
  void FocusGet(T oFocus[3]) const;

  /**
   * Set focus point or virtual source
   *
   * @param[in] iFocus
   */
  void FocusSet(const T iFocus[3]);

  /**
   * Get second focus point
   *
   * @param[out] oFocus
   */
  void Focus2Get(T oFocus[3]) const;

  /**
   * Set second focus point
   *
   * @param[in] iFocus
   */
  void Focus2Set(const T iFocus[3]);

  /**
   * Get center focus point
   *
   * @param[out] oFocus
   */
  void CenterFocusGet(T oFocus[3]) const;

  /**
   * Set center focus point
   *
   * @param[in] iFocus
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
   * @param[in] iFocusingType
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
   * @param[in] nThreads
   */
  int NThreadsSet(const size_t &nThreads);

  /**
   * Get excitation center frequency
   *
   * @return
   */
  const T& F0Get() const;

  /**
   * Set excitation center frequency
   *
   * @param[in] f0
   */
  void F0Set(const T& f0);

  /**
   * Get width of pulse [s]
   *
   * @return
   */
  const T& WGet() const;

  /**
   * Set width of pulse [s]. TODO: Figure this out
   *
   * @param[in] w
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
   * @param[in] arg
   */
  void SysParmSet(const sysparm_t<T> *arg);

#ifdef FNM_PULSED_WAVE
  sofus::AperturePulses<T>* PulsesGet();
#endif

  const T& DensityGet() const;
  void DensitySet(const T& rho);

  /**
   * Get sampling frequency
   *
   * @return
   */
  const T& FsGet() const;

  /**
   * Set sampling frequency
   *
   * @param[in] fs
   */
  int FsSet(const T& fs);


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
   * @param[in] c
   */
  void CSet(const T& c);

  /**
   * Get transmit f-number
   *
   *
   * @return
   */
  const T& XmtFNumberGet() const;

  /**
   * Set transmit f-number
   *
   * @param[in] fnumber
   */
  void XmtFNumberSet(const T& fnumber);

  /**
   * Get element definitions
   *
   * @param[out] out
   * @param[out] nElements
   * @param[out] nParams
   */
  void ElementsGet(T** out, size_t* nElements,
                   size_t* nParams) const;

  int ElementsSet(const T* pos, const size_t nPositions, const size_t nDim);

  /*! \fn bool Aperture::ElementsSet(const T* pos, const size_t nPositions, const size_t nDim)
   *  \brief Set element positions
   *  \param[in] pos  Input data T[nElements][8]
   *  \param[in] nPositions
   *  \param[in] nDim Must equal 8
   *  \exception std::runtime_error nDim != 8
   *  \return true on success
   */

  /**
   * Get sub-element definitions
   *
   * @param[out] out
   * @param[out] nElements
   * @param[out] nSubElements
   * @param[out] nParams
   */
  void SubElementsGet(T** out, size_t* nElements,
                      size_t* nSubElements, size_t* nParams) const;

  /**
   * Set sub-element definitions
   *
   * @param[in] pos
   * @param[in] nElements
   * @param[in] nSubElementsPerElement
   * @param[in] nDim
   *
   * @return true on success
   */
  int SubElementsSet(const T* pos, const size_t nElements,
                     const size_t nSubElementsPerElement, const size_t nDim);

  /**
   * Get apodization
   *
   * @param[out] data
   * @param[out] nData
   */
  void ApodizationGet(T** data, size_t* nData) const;

  /**
   * Set apodization
   *
   * @param[in] data
   * @param[in] nData
   */
  int ApodizationSet(const T* data, const size_t nData);

  /**
   * Get apodization type
   *
   */
  int ApodizationTypeGet() const;

  /**
   * Set apodization
   *
   * @param[in] iApodizationType
   */
  void ApodizationTypeSet(const int iApodizationType);

  ///@}

  /**
   * Get excitation type
   *
   *
   * @return
   */
  int ExcitationTypeGet() const;

  /**
   * Set Excitation type (0: ToneBurst, 1: Hamming-weighted pulse)
   *
   * @param iExcitationType
   */
  void ExcitationTypeSet(const int iExcitationType);

  /**
   * Get transducer center frequency used for computing impulse response
   *
   * @return
   */
  const T& FCGet() const;

  /**
   * Set pulse center frequency
   *
   * @param[out] fc
   */
  void FCSet(const T& fc);

  /**
   * Get bandwidth
   *
   * @return
   */
  const T& BandWidthGet() const;

  /**
   * Set bandwidth
   *
   * @param[in] value
   */
  void BandWidthSet(const T& value);


#if FNM_PULSED_WAVE
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
   * @param[in] value
   */
  void NormalizeSet(const bool& value);

  /**
   * Create a number of lines, which can be used for focusing
   *
   * @param[out] nLines
   * @param[out] obj
   */
  void FocusLinesCreate(size_t nLines, sofus::FocusLineList<T>** obj) const;


  /**
   * Get impulse type
   *
   * @return
   */
  int ImpulseTypeGet() const;

  /**
   * Set impulse type
   *
   * @param iImpulseType
   */
  void ImpulseTypeSet(const int iImpulseType);



  /**
   * Get excitation (reference to or view of)
   *
   * @param[out] data
   * @param[out] nData
   */
  void ExcitationGet(T** data, size_t* nData) const;

  /**
   * Set excitation
   *
   * @param[in] data
   * @param[in] nData
   */
  int ExcitationSet(const T* data,
                    const size_t nData);

  /**
   * Get impulse (reference to or view of)
   *
   * @param[out] data
   * @param[out] nData
   */
  void ImpulseGet(T** data, size_t* nData) const;

  /**
   * Set impulse
   *
   * @param[in] data
   * @param[in] nData
   */
  int ImpulseSet(const T* data,
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
   * @param[in] nDivW
   */
  int NDivWSet(const size_t& nDivW);

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
   * @param[in] nDivH
   */
  int NDivHSet(const size_t& nDivH);

#ifdef USE_PROGRESS_BAR
  /**
   * Set the progress bar used for notification
   *
   * @param[in] pbar
   */
  void ProgressBarSet(sps::ProgressBarInterface* pbar);
#endif

#ifdef FNM_CLOSURE_FUNCTIONS

  ApertureData<T>* DataGet();

#  ifndef SWIG_VERSION
  int RwFloatParamSet(int pSel, const T* f, size_t nDim, va_list args);
  int RwFloatParamGet(int pSel, T** f, size_t nDim, va_list args);
  // Using this pointer (not working)
  int RwFloatParamGetUsingThis(int pSel, T** f, size_t nDim, va_list args);
#  endif

  int RwFloatParamSet(int fsel, const T* iMultiData, size_t nDim, ...);
  int RwFloatParamGet(int fsel, T** oMultiData, size_t nDim, ...);

  void RwBooleanParamSet0D(int fsel, const bool& value);
  void RwSizeTParamSet0D(int fsel, const size_t& value);

  void RwIntegerParamSet0D(int fsel, const int& value);

  int ParameterInfoGet(const RwParamType& param, Type* type,
                       size_t* nDims);

  static int ParameterInfoListGet(ApertureProperty** ppPropertyInfo);

  static int ParameterInfoListDestroy(ApertureProperty* ppPropertyInfo);

  // TODO(JEM): Expose list of parameters
#endif

  /** @name Static variables
   *
   */
  ///@{

  static const size_t nVerticesPerElement = 4; /**< Number of corners for an element */

  static const size_t nElementPosParameters = 8; /**< Number of parameters for an element */

#if FNM_PULSED_WAVE
  // TODO: Not needed
  static bool normalize;  ///< Normalize responses
#endif

#if (defined(SWIG_VERSION) && (SWIG_VERSION > 0x030000)) && defined(__GNUG__)
  // Supported by SWIG3.0
  static constexpr T Neper_dB = 0.11512925464970231;
  static constexpr T dB_Neper = 8.685889638065035;
#else
  static const T Neper_dB;  /**< Nepers per dB, 1 / ( 20 log(e)) */
  static const T dB_Neper;  /**< dB per Neper, 20 log(e) */
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
   * Calculate pulsed wave field using fast near-field method
   *
   * Experimental SSE4 version
   *
   * @param[in] pos
   * @param[in] nPositions
   * @param[in] nDim
   * @param[out] odata
   * @param[out] nSignals
   * @param[out] nSamples
   * @param[in] mask
   *
   * @return
   */
  T CalcPwFnmVectorized(const T* pos, const size_t nPositions,
                        const size_t nDim,
                        T** odata, size_t* nSignals, size_t* nSamples,
                        int mask = 0x1F);


  /**
   *
   *
   * @param[in] pos
   * @param[in] nPositions
   * @param[in] nDim
   * @param[out] odata
   * @param[out] nSignals
   * @param[out] nSamples
   * @param[in] mask
   *
   * @return
   */
  T CalcPwFnmVecThreaded(const T* pos,const size_t nPositions,
                         const size_t nDim,
                         T** odata,size_t* nSignals,size_t* nSamples,
                         int mask = 0x1F);
#endif

  /** @name FNM Calculation functions
   *
   */
  ///@{

  /**
   * Compute CW response at multiple positions. Reference
   * implementation. See @ref fnm::CalcCwFieldRef.
   *
   *
   * @param[in] pos
   * @param[in] nPositions
   * @param[in] nDim
   * @param[out] odata
   * @param[out] nOutPositions
   *
   * @return
   */
  int CalcCwFieldRef(const T* pos, const size_t nPositions, const size_t nDim,
                     std::complex<T>** odata, size_t* nOutPositions);


  /**
   * Compute CW response at multiple positions. Number of integrals is
   * reduced to four
   *
   * @param[in] pos
   * @param[in] nPositions
   * @param[in] nDim
   * @param[out] odata
   * @param[out] nOutPositions
   *
   * @return
   */
  int CalcCwFieldFourRef(const T* pos, const size_t nPositions,
                         const size_t nDim,
                         std::complex<T>** odata, size_t* nOutPositions);

  //

  int CalcCwEcho(const Aperture<T>* pOther,
                 const T* pos, const size_t nPositions, const size_t nDim,
                 std::complex<T>** odata, size_t* nOutPositions);


  int CalcCwBack(const T* pos,
                 const size_t nPositions,
                 const size_t nDim,
                 const std::complex<T>* pFieldValues,
                 const size_t nComplexValues,
                 std::complex<T>** odata,
                 size_t* nOutPositions);

  int CalcCwTimeReversal(const Aperture<T>* pOther,
                         const T* pos, const size_t nPositions, const size_t nDim,
                         const std::complex<T>* pFieldValues, // At elements
                         const size_t nComplexValues,
                         std::complex<T>** odata, size_t* nOutPositions, size_t* nSamples);
  /**
   * Compute CW response at multiple positions (uses SIMD). The
   * ranges of integration are reduced to give the most accurate
   * result with less abcissas.
   *
   * @param[in] pos
   * @param[in] nPositions
   * @param[in] nDim
   * @param[out] odata
   * @param[out] nOutPositions
   *
   * @return
   */
  int CalcCwFast(const T* pos, const size_t nPositions, const size_t nDim,
                 std::complex<T>** odata, size_t* nOutPositions);

#ifndef FNM_DOUBLE_SUPPORT
  /**
   * Same as CalcCwField, but uses SIMD.
   *
   * @param[in] pos
   * @param[in] nPositions
   * @param[in] nDim
   * @param[out] odata
   * @param[out] nOutPositions
   *
   * @return
   */
  int CalcCwFieldNaiveFast(const T* pos, const size_t nPositions,
                           const size_t nDim,
                           std::complex<T>** odata, size_t* nOutPositions);
#endif

  /**
   * Compute CW response at multiple positions. The range of
   * integration is naive so it is accurate when projections lie
   * inside an element, but requires a huge amount of abcissas to
   * get a usuable result, when projections lie outside an element.
   *
   * @param[in] pos
   * @param[in] nPositions
   * @param[in] nDim
   * @param[out] odata
   * @param[out] nOutPositions
   */
  int CalcCwFieldNaive(const T* pos, const size_t nPositions, const size_t nDim,
                       std::complex<T>** odata, size_t* nOutPositions);

  /*** @} */


  /**
   * Calculate pulsed wave field using fast near-field method.
   *
   * @param[in] pos
   * @param[in] nPositions
   * @param[in] nDim
   * @param[out] odata
   * @param[out] nSignals
   * @param[out] nSamples
   * @param[in] mask
   *
   * @return
   */
  T CalcPwFnmThreaded(const T* pos, const size_t nPositions,
                      const size_t nDim,
                      T** odata, size_t* nSignals, size_t* nSamples,
                      int mask = 0x1F);

  /**
   * Compute transient for first element.
   *
   * Mask bits: 0x01: Direct,
   *            0x02: Edge (q0->q3)
   *            0x04: Edge (q0->q1)
   *            0x08: Edge (q1->q2)
   *            0x10: Edge (q2->q3)
   *            0x1F: All contributions (default)
   *
   * @param[in] pos
   * @param[in] nPositions
   * @param[in] nDim
   * @param[out] odata
   * @param[out] nSignals
   * @param[out] nSamples
   * @param[in] mask
   *
   * @return
   */
  T CalcTransientSingleElementNoDelay(const T* pos, const size_t nPositions,
                                      const size_t nDim,
                                      T** odata, size_t* nSignals,
                                      size_t* nSamples, int mask = 0x1F);

#ifdef FNM_PULSED_WAVE

  T CalcPwBackScat(const T* pos, const size_t nPositions,
                   const size_t nDim,
                   const T* data, const size_t nData,  // Amplitudes
                   T** odata, size_t* nSignals, size_t* nSamples);

  /** @name SOFUS accellerated calculation functions
   *
   */
  ///@{

  /**
   *
   *
   * @param[in] other
   * @param[in] iFocus
   * @param[in] odata
   * @param[in] nSignals
   * @param[out] nSamples
   * @param[out] sizeAndOffsets
   * @param[out] nFilters
   * @param[out] nTwo
   * @param[out] data
   * @param[out] nData
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
   * @param[in] other
   * @param[in] iFocus
   * @param[in] tStart
   * @param[in] iData
   * @param[in] nChannels
   * @param[in] nSamples
   * @param[out] data
   * @param[out] nData
   *
   * @return
   */
  T CalcSmfApply(const Aperture<T>* other,
                 const T iFocus[3],
                 const T tStart,
                 const T* iData, const size_t nChannels,
                 const size_t nSamples,
                 T** data, size_t* nData);

  /**
   * Rename to sim image
   *
   * @param[in] other
   * @param[in] pos
   * @param[in] nPositions
   * @param[in] nDim
   * @param[in] data
   * @param[in] nData
   * @param[out] odata
   * @param[out] nSignals
   * @param[out] nSamples
   *
   * @return
   */
  T CalcScat(const Aperture<T>* other, const T* pos, const size_t nPositions,
             const size_t nDim,
             const T* data, const size_t nData,  // Amplitudes
             T** odata, size_t* nSignals, size_t* nSamples);

  T CalcScatThreaded(const Aperture<T>* other,
                     const sofus::FocusLineList<T>* xmtLines,
                     const T* pos, const size_t nPositions, const size_t nDim,
                     const T* data, const size_t nData,  // Amplitudes
                     T** odata, size_t* nSignals, size_t* nSamples);
  /**
   *
   *
   * @param[in] other
   * @param[in] pos
   * @param[in] nPositions
   * @param[in] nDim
   * @param[in] data
   * @param[in] nData
   * @param[out] odata
   * @param[out] nSignals
   * @param[out] nSamples
   *
   * @return
   */
  T CalcPwEcho(const Aperture<T>* other,
               const T* pos, const size_t nPositions, const size_t nDim,
               const T* data, const size_t nData,
               T** odata, size_t* nSignals, size_t* nSamples);

  T CalcPwFnmEcho(
    const Aperture<T>* other,
    const T* pos, const size_t nPositions, const size_t nDim,
    const T* data, const size_t nData,
    T** odata, size_t* nSignals, size_t* nSamples);

  /**
   * Compute pulsed-wave field
   *
   * @param[in] pos
   * @param[in] nPositions
   * @param[in] nDim
   * @param[out] odata
   * @param[out] nSignals
   * @param[out] nSamples
   *
   * @return
   */
  T CalcPwField(const T* pos, const size_t nPositions, const size_t nDim,
                T** odata, size_t* nSignals, size_t* nSamples);

  /**
   * Reference implementation of Aperture::CalcPwField
   *
   * @param[in] pos
   * @param[in] nPositions
   * @param[in] nDim
   * @param[out] odata
   * @param[out] nSignals
   * @param[out] nSamples
   *
   * @return
   */
  T CalcPwFieldRef(const T* pos, const size_t nPositions, const size_t nDim,
                   T** odata, size_t* nSignals, size_t* nSamples);


  /**
   * Threaded version of Aperture::CalcPwField. \bug: Sometimes to little memory is allocated
   *
   * @param[in] pos
   * @param[in] nPositions
   * @param[in] nDim
   * @param[out] odata
   * @param[out] nSignals
   * @param[out] nSamples
   *
   * @return
   */
  T CalcPwFieldThreaded(const T* pos, const size_t nPositions,
                        const size_t nDim,
                        T** odata, size_t* nSignals, size_t* nSamples);
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
   * @param[in] y0
   * @param[in] nx
   * @param[in] ny
   * @param[in] dx
   * @param[in] dy
   * @param[in] Nx
   * @param[in] Ny
   * @param[out] p1
   * @param[out] onx
   * @param[out] ony
   * @param[out] onz
   *
   * @return
   */
  int CalcAsa(const T* y0, const size_t nx, const size_t ny,
              const T dx, const T dy,
              const size_t Nx, const size_t Ny,
              std::complex<T>** p1, size_t *onx, size_t *ony, size_t *onz);
  ///@}

#endif

  /**
  * Get default system parameters. The system parameters are shared among
  * all apertures, unless system paramters are changed using
  * Aperture::SysParmGet and Aperture::SysParmSet for a specific array.
  *
  *
  * @return
  */
  static sysparm_t<T>* DefaultSysParmGet();

 private:
  /**
   * Internal function for getting elements or sub-elements. The
   * function is used by Aperture::ElementsGet and
   * Aperture::SubElementsGet.
   *
   * @param[in] sub 0 for elements, 1 for subelements
   * @param[out] out
   * @param[out] nElements
   * @param[out] nSubElements
   * @param[out] nParams
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
   * @param[in] pos
   * @param[in] nPositions
   * @param[in] nDim
   * @param[out] odata
   * @param[out] nOutPositions
   *
   * @return
   */
  int CalcCwFocus(const T* pos, const size_t nPositions, const size_t nDim,
                  std::complex<T>** odata, size_t* nOutPositions) const;

  /**
   * Reference implementation of Aperture::CalcCwFocus
   *
   * @param[in] pos
   * @param[in] nPositions
   * @param[in] nDim
   * @param[out] odata
   * @param[out] nOutPositions
   *
   * @deprecated
   *
   * @return
   */
  int CalcCwFocusRef(const T* pos, const size_t nPositions, const size_t nDim,
                     std::complex<T>** odata, size_t* nOutPositions);

 private:
  /// PIMPL data design (all pointers are opaque)

  ///< Test structure (replacing m_sysparm)
  bla<T>* m_pTest;

  ///< Sysparm
  sysparm_t<T>* m_sysparm;

  ///< Data
  ApertureData<T>* m_pData;

  ProgIF* m_pbar;

};
}  // namespace fnm
/*! @} End of FNM Group */

#ifdef _MSC_VER
#pragma warning(pop)
#endif

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
