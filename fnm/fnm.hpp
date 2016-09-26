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

/** @defgroup FNM Fast Near-field Method
 * @brief This module is used for computing responses using the fast nearfield method (FNM)
 *
 */

#include <fnm/fnm_export.h>
#include <fnm/config.h>

#include <complex>

//! Fast Nearfield Method interfaces and implementations
namespace fnm {
  /// Forward-declare ApertureData
  template<class T>
  struct ApertureData;
}

namespace std {
  class runtime_error;
}

/*!
 * @addtogroup FNM
 * @{
 */

namespace fnm {

  /*! \brief Sysparm structure
   *
   *
   * A structure containing global simulation parameters
   */
  template <typename T>
  struct FNM_EXPORT sysparm_t {
    /// Speed of sound
    T c;
    /// Number of width abcissas
    size_t nDivW;
    /// Number of height abcissas
    size_t nDivH;
  };

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
     * Get phases of elements
     *
     * @param phases
     * @param nPhases
     */
    void PhasesGet(T** phases, size_t* nPhases);


    void RectanglesGet(T** out, size_t* nElements,
                       size_t* nSubElements, size_t* nParams) const;

    // TODO: Group properly
    //a}


    //a{
    /** Read-write attributes */

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

    //a}

    /**
     * Compute CW response at multiple positions (uses SIMD)
     *
     * @param pos
     * @param nPositions
     * @param nDim
     * @param odata
     * @param nOutPositions
     */
    int CalcCwFast(const T* pos, const size_t nPositions, const size_t nDim,
                   std::complex<T>** odata, size_t* nOutPositions);

    /**
     * Update phases after setting focus. This is called automatically.
     *
     */
    void FocusUpdate();

    // Experimental undocumented functions
    void CalcCwField(const T* pos, const size_t nPositions, const size_t nDim,
                     std::complex<T>** odata, size_t* nOutPositions);

    /**
     * Compute CW response at multiple positions. Reference implementation.
     *
     * @param pos
     * @param nPositions
     * @param nDim
     * @param odata
     * @param nOutPositions
     */
    void CalcCwFieldRef(const T* pos, const size_t nPositions, const size_t nDim,
                        std::complex<T>** odata, size_t* nOutPositions);

    void CalcCwField2(const T* pos, const size_t nPositions, const size_t nDim,
                      std::complex<T>** odata, size_t* nOutPositions);

    void CalcCwFocus(const T* pos, const size_t nPositions, const size_t nDim,
                     std::complex<T>** odata, size_t* nOutPositions);

    //a{
    /** Static variables */
    static const size_t nVerticesPerElement = 4; /**< Number of corners for an element */

    static const size_t nElementPosParameters = 8; /**< Number of parameters for an element */

    /// System parameters (We could make a long list of functions friends of Aperture)
    static sysparm_t<T> _sysparm;
    //a}
  private:

    /// Number of threads
    static size_t nthreads;

    /**
     * Initialize elements. This function must be called whenever
     * elements are changed.
     *
     */
    void initElements();

    // TODO: Move these to .cpp translation unit, must take pointer to Aperture<T>

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
    void* CalcThreaded(void* ptarg);
#elif defined(_WIN32)
    unsigned int __stdcall CalcThreaded(void *ptarg);
#endif

    ///< Data
    ApertureData<T>* m_data;
  };
}
/*! @} End of FNM Group */


/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
