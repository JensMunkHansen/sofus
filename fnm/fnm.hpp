#pragma once

#include <fnm/fnm_export.h>
#include <fnm/config.h>
#include <sps/smath.hpp>

#include <complex>
#include <vector>

#define N_MAX_THREADS 8

namespace fnm {
  /// Forward-declare ApertureData
  template<class T>
  struct FNM_EXPORT ApertureData;

  template <typename T>
  struct FNM_EXPORT element_t;

  template <typename T>
  struct FNM_EXPORT sysparm_t;
}

namespace fnm {

  template <class T>
  class FNM_EXPORT Aperture {
  private:
    
    struct thread_arg {
      size_t iPointBegin;
      size_t iPointEnd;
      const T* pos;
      std::complex<T>* field;
      T k;
      T* uxs;
      T* uweights;
      size_t nDivU;
      T* vxs;
      T* vweights;
      size_t nDivV;
      T* normals;
      T* uvectors;
      T* vvectors;
      size_t thread_id;
      int cpu_id;
    };
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
     * @param pitch
     */
    Aperture(const size_t nElements, const T width, const T kerf, const T height);

    ~Aperture();

    void PositionsGet(T** pos, size_t* nElements, size_t* nParams) const;

    void PhasesGet(T** phases, size_t* nPhases);

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

    const size_t& NThreadsGet() const;

    void NThreadsSet(const size_t &nThreads);

    /**
     * Get element count
     *
     * @return
     */
    const size_t& NElementsGet() const;

    const size_t& NSubElementsGet() const;
    
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
    void CSet (const T c);

    /**
     * Set element positions
     *
     * @param pos  Input data T[nElements][3]
     * @param nElements
     * @param nDim Must equal 3
     *
     * @return
     */
    bool ElementsSet(const T* pos, const size_t nElements, const size_t nDim) throw (std::runtime_error);

    /** 
     * Set sub-element positions
     * 
     * @param pos 
     * @param nElements 
     * @param nSubElementsPerElement
     * @param nDim 
     * 
     * @return 
     */
    bool SubElementsSet(const T* pos, const size_t nElements, const size_t nSubElementsPerElement, const size_t nDim);
    
    const size_t& NDivWGet() const;

    void NDivWSet(const size_t nDivW);

    const size_t& NDivHGet() const;

    void NDivHSet(const size_t nDivH);

    void RectanglesGet(T** outRectangles, size_t* nElements,
                       size_t* nSubElements, size_t* nCornerCoordinates) const;
    
    void CalcCwField(const T* pos, const size_t nPositions, const size_t nDim,
                     std::complex<T>** odata, size_t* nOutPositions);
    void CalcCwFieldRef(const T* pos, const size_t nPositions, const size_t nDim,
                        std::complex<T>** odata, size_t* nOutPositions);

    void CalcCwField2(const T* pos, const size_t nPositions, const size_t nDim,
                      std::complex<T>** odata, size_t* nOutPositions);

    void CalcCwFocus(const T* pos, const size_t nPositions, const size_t nDim,
                     std::complex<T>** odata, size_t* nOutPositions);

    void CalcCwFast(const T* pos, const size_t nPositions, const size_t nDim,
                    std::complex<T>** odata, size_t* nOutPositions);

    //@{
    /** Static variables */
    static const size_t nVerticesPerElement = 4; /**< Number of corners for an element */

    static const size_t nElementPosParameters = 8; /**< Number of parameters for an element */
    //@}
    
    /// System parameters
    static sysparm_t<T> _sysparm;

    static size_t nthreads;

#ifdef HAVE_PTHREAD_H
    static pthread_t threads[N_MAX_THREADS];
    static pthread_attr_t attr;
# ifdef HAVE_MQUEUE_H
    static bool threads_initialized;
    static mqd_t mqd_master;
    static mqd_t mqd_client;
# endif
#endif

    
  private:
    static thread_arg threadarg[N_MAX_THREADS];
    /**
     * Initialize elements. This function must be called whenever
     * elements are changed.
     *
     */
    void initElements();

#if defined(HAVE_PTHREAD_H)
    void* CalcThreaded(void* ptarg);
# ifdef HAVE_MQUEUE_H
    void CalcThreaded(void* ptarg);
# endif
#elif defined(_WIN32)
    unsigned int __stdcall CalcThreaded(void *ptarg);
#endif
    ApertureData<T>* m_data;
  };
}
    
template <class T>
std::complex<T> CalcHz(const T& s,
                       const T& l,
                       const T& z,
                       const T& k,
                       const T* uxs,
                       const T* uweights,
                       const size_t nUs,
                       const T* vxs,
                       const T* vweights,
                       const size_t nVs);

template <class T>
STATIC_INLINE_BEGIN std::complex<T> CalcSingle(const T& s1,
                                               const T& s2,
                                               const T& l,
                                               const T& z,
                                               const T& k,
                                               const T* uxs,
                                               const T* uweights,
                                               const size_t nUs) STATIC_INLINE_END;

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
