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
#include <fnm/fnm_calc.hpp>

#ifdef FNM_PULSED_WAVE
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

#ifdef HAVE_MQUEUE_H
# include <mqueue.h>
#endif

#ifdef _MSC_VER
#include <BaseTsd.h>
#endif

#ifndef MSGMAX
# define MSGMAX 8192
#endif

namespace fnm {

  /////////////////////////////////////////////////
  // Types visible for this compilation unit only
  /////////////////////////////////////////////////
  template <class T>
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
    size_t thread_id;
    int cpu_id;
  };

#if HAVE_THREAD
  // Specializations contain static arrays
  template<class T>
  struct ApertureThreadArgs {
    static thread_arg<T> args[N_MAX_THREADS];
  };
#endif

#ifdef HAVE_PTHREAD_H
# ifdef HAVE_MQUEUE_H
  template <class T>
  struct ApertureQueue {
    static bool threads_initialized;
    static mqd_t mqd_master;
    static mqd_t mqd_client;
  };
# endif
  pthread_t threads[N_MAX_THREADS];
  static pthread_attr_t attr = {{0}};
#endif

  /////////////////////////////////////////////////
  // Static content declared in interface
  /////////////////////////////////////////////////

#if defined(C99) && !defined(__STRICT_ANSI__)
  template <class T>
  sysparm_t<T> Aperture<T>::_sysparm =
    (sysparm_t<T>)
  {
    .c = 1500.0, .nDivW = 16, .DivH = 16
  };
#else
  template <class T>
  sysparm_t<T> Aperture<T>::_sysparm = {T(1500.0), size_t(16), size_t(16)};
#endif

#ifdef FNM_PULSED_WAVE
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
    m_data = (ApertureData<T>*) _mm_malloc(sizeof(ApertureData<T>),16);
    new (this->m_data) ApertureData<T>();
    m_pbar = NULL;
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

    // Experimental stuff with queues
#ifdef HAVE_MQUEUE_H
    size_t i;
    unsigned int prio = 0;

    if (ApertureQueue<T>::threads_initialized) {
      // Signal threads to exit
      CallErr(ApertureQueue<T>::mqd_client = mq_open,
              ("/jmh-client", O_WRONLY | O_CREAT,
               S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
               NULL));

      for (i=0; i<Aperture<T>::nthreads; i++) {
        prio = i;
        if (mq_send(ApertureQueue<T>::mqd_client,"EXIT",4,prio)==-1)
          perror ("mq_send()");
      }

      // Wait for threads to exit. TODO: Do we need this???
      for (i = 0; i < Aperture<T>::nthreads; i++) {
        CallErr(pthread_join,(threads[i],NULL));
      }

      mq_close(ApertureQueue<T>::mqd_client);

      ApertureQueue<T>::threads_initialized = false;
      debug_print("Threads have exited\n");
      CallErr(pthread_attr_destroy,(&attr));
    }
#endif
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
    if (m_data->m_focus_valid == m_data->m_focus_type) {
      return;
    }

    const size_t nElements = m_data->m_nelements;

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
        m_data->m_phases[iElement] = T(M_2PI) * (delays[iElement] - maxDelay) * m_data->m_f0;
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

      // This function allocates!!! TODO: Move CalcCwFocus to fnm_calc
      this->CalcCwFocus(positions.get(),nPositions,3,&pFieldValues,&nFieldValues);

      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
        m_data->m_phases[iElement] = - std::arg<T>(pFieldValues[iElement]);
      }
      free(pFieldValues);
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
    *odata = (T*) malloc(nElements*sizeof(T));
    if (nElements > 0) {
      memcpy(*odata,m_data->m_phases.get(),nElements*sizeof(T));
      *nData = nElements;
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
      m_data->m_f0 = f0;
    }
  }

  template <class T>
  void Aperture<T>::SysParmSet(const sysparm_t<T> *arg)
  {
    Aperture<T>::_sysparm.c = arg->c;
    Aperture<T>::_sysparm.nDivH = arg->nDivH;
    Aperture<T>::_sysparm.nDivW = arg->nDivW;
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

    // Static is needed to avoid temporaries (if a view was returned). We allocate, so this is not a temporary
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

      auto& elements = m_data->m_elements;

      for (size_t iElement = 0 ; iElement < _nElements ; iElement++) {
        for (size_t jElement = 0 ; jElement < nSubElementsPerElement ; jElement++) {

          arr[iElement * nSubElementsPerElement * nElePosParams +
              jElement * nElePosParams + 0] = elements[iElement][jElement].hw;

          arr[iElement * nSubElementsPerElement * nElePosParams +
              jElement * nElePosParams + 1] = elements[iElement][jElement].hh;

          memcpy(&arr[iElement * nSubElementsPerElement * nElePosParams + jElement*nElePosParams + 2],
                 &elements[iElement][jElement].center[0],
                 sizeof(sps::point_t<T>));
          arr[iElement * nSubElementsPerElement * nElePosParams + jElement*nElePosParams+5] = elements[iElement][jElement].euler.alpha;
          arr[iElement * nSubElementsPerElement * nElePosParams + jElement*nElePosParams+6] = elements[iElement][jElement].euler.beta;
          arr[iElement * nSubElementsPerElement * nElePosParams + jElement*nElePosParams+7] = elements[iElement][jElement].euler.gamma;
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

  template <class T>
  bool Aperture<T>::SubElementsSet(const T* pos, const size_t nElements,
                                   const size_t nSubElementsPerElement, const size_t nDim)
  {
    bool retval = true;
    // We are allowed to set 0 elements
    if ((nDim != Aperture<T>::nElementPosParameters) || (nSubElementsPerElement < 1)) {
      retval = false;
      return retval;
    }

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

  template <class T>
  void Aperture<T>::ProgressBarSet(sps::ProgressBarInterface* pbar)
  {
    this->m_pbar = pbar;
  }

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

    if (m_data->m_focus_type != FocusingType::Pythagorean) {
      fprintf(stderr, "Focusing type must be Pythagorean for pulsed wave calculations\n");
    }
    this->FocusUpdate();

    sofus::sysparm_t<T> sysparm;
    sysparm.c  = Aperture<T>::_sysparm.c;
    sysparm.fs = Aperture<T>::fs;

    // Allocates
    T retval = sofus::CalcPwField(sysparm,
                                  m_data,
                                  m_pulses,
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

    const T* apodizations = m_data->m_apodizations.get();

    T apodization;


    // Need weights and abcissa values

    size_t nDivW = Aperture<T>::_sysparm.nDivW;
    size_t nDivH = Aperture<T>::_sysparm.nDivH;

    // TODO: For FDTSD store these
    auto uweights = sps::deleted_aligned_array_create<T>(nDivW);
    auto vweights = sps::deleted_aligned_array_create<T>(nDivH);
    auto uxs      = sps::deleted_aligned_array_create<T>(nDivW);
    auto vxs      = sps::deleted_aligned_array_create<T>(nDivH);

    memset(uxs.get(),0,sizeof(T)*nDivW);
    memset(vxs.get(),0,sizeof(T)*nDivH);

    memset(uweights.get(),0,sizeof(T)*nDivW);
    memset(vweights.get(),0,sizeof(T)*nDivH);

    // Common weights and abcissa
    for (size_t i = 0 ; i < nDivW ; i++) {
      gl::GLNode qp = gl::GL(nDivW,i);
      uxs[i]        = T(qp.value);
      uweights[i]   = T(qp.weight);
    }

    for (size_t i = 0 ; i < nDivH ; i++) {
      gl::GLNode qp = gl::GL(nDivH, i);
      vxs[i]      = T(qp.value);
      vweights[i] = T(qp.weight);
    }

    sps::deleted_aligned_multi_array<sps::element_t<T>,2>& elements = m_data->m_elements;

    sps::point_t<T> projection;

    for (size_t iPoint = 0 ; iPoint < nPositions ; iPoint++) {

      __m128 vec_point = _mm_set_ps(0.0f, float(pos[iPoint*3 + 2]),float(pos[iPoint*3 + 1]),float(pos[iPoint*3]));

      std::complex<T> final = std::complex<T>(T(0.0),T(0.0));

      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {

        apodization = apodizations[iElement];

        if (apodization != 0.0) {

          for (size_t jElement = 0 ; jElement < this->m_data->m_nsubelements ; jElement++) {

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
            final = final + apodization * result * exp(std::complex<T>(0,m_data->m_phases[iElement]));
          }
        }
      }
      (*odata)[iPoint].real(final.real());
      (*odata)[iPoint].imag(final.imag());
    }
    return 0;
  }
#endif

#ifndef FNM_DOUBLE_SUPPORT
  // Reduced integral and fast (the one to use). TODO: Move to fnm_calc and use ApertureData<T>
  template <class T>
  int Aperture<T>::CalcCwFast(const T* pos,
                              const size_t nPositions,
                              const size_t nDim,
                              std::complex<T>** odata,
                              size_t* nOutPositions)
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

    // This leaks
    this->FocusUpdate();

    const T lambda = Aperture<T>::_sysparm.c / m_data->m_f0;
    const T k = T(M_2PI)/lambda;

    size_t nDivW = Aperture<T>::_sysparm.nDivW;
    size_t nDivH = Aperture<T>::_sysparm.nDivH;

    // TODO: For FDTSD store these
    auto uweights = sps::deleted_aligned_array_create<T>(nDivW);
    auto vweights = sps::deleted_aligned_array_create<T>(nDivH);
    auto uxs      = sps::deleted_aligned_array_create<T>(nDivW);
    auto vxs      = sps::deleted_aligned_array_create<T>(nDivH);

    memset(uxs.get(),0,sizeof(T)*nDivW);
    memset(vxs.get(),0,sizeof(T)*nDivH);

    memset(uweights.get(),0,sizeof(T)*nDivW);
    memset(vweights.get(),0,sizeof(T)*nDivH);

    // Common weights and abcissa
    for (size_t i = 0 ; i < nDivW ; i++) {
      gl::GLNode qp = gl::GL(nDivW,i);
      uxs[i]        = T(qp.value);
      uweights[i]   = T(qp.weight);
    }

    for (size_t i = 0 ; i < nDivH ; i++) {
      gl::GLNode qp = gl::GL(nDivH, i);
      vxs[i]      = T(qp.value);
      vweights[i] = T(qp.weight);
    }

    // Threaded or not

#ifndef HAVE_THREAD
    fnm::thread_arg<T> threadarg;
    threadarg.iPointBegin = 0;
    threadarg.iPointEnd   = nPositions;
    threadarg.pos         = pos;
    threadarg.field       = *odata;
    threadarg.k           = k;
    threadarg.uxs         = uxs.get();
    threadarg.uweights    = uweights.get();
    threadarg.nDivU       = nDivW; // Works for threaded
    threadarg.vxs         = vxs.get();
    threadarg.vweights    = vweights.get();
    threadarg.nDivV       = nDivH;
    threadarg.thread_id   = 0;
    threadarg.cpu_id      = 0;

# ifdef _WIN32
    retval                = (unsigned int)CalcThreaded((void*)&threadarg);
# else
    void* thread_retval   = CalcThreaded((void*)&threadarg);
    SPS_UNREFERENCED_PARAMETER(thread_retval);
# endif
#else

    int nproc;

# if defined(HAVE_PTHREAD_H)
# elif defined(_WIN32)
    unsigned int threadID;
    uintptr_t threads[N_MAX_THREADS];
# endif

    nproc = getncpus();

# ifdef HAVE_MQUEUE_H
    sps::pthread_launcher<Aperture<T>,
        &Aperture<T>::CalcThreadFunc> launcher[N_MAX_THREADS];
# else
    sps::pthread_launcher<Aperture<T>,
        &Aperture<T>::CalcThreaded> launcher[N_MAX_THREADS];
# endif

    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

# ifdef HAVE_PTHREAD_H
    // TODO: Shoudl we do this over and over if using message queues???
#  ifdef HAVE_MQUEUE_H
    if ((!ApertureQueue<T>::threads_initialized) && !(nthreads>N_MAX_THREADS)) {
      CallErr(pthread_attr_init,(&attr));
      //CallErr(pthread_attr_setdetachstate, (&attr, PTHREAD_CREATE_JOINABLE));
      CallErr(pthread_attr_setdetachstate, (&attr, PTHREAD_CREATE_DETACHED));
    }
#  else
    CallErr(pthread_attr_init,(&attr));
    CallErr(pthread_attr_setdetachstate, (&attr, PTHREAD_CREATE_JOINABLE));
#  endif
# endif

    debug_print("nthreads: %zu\n", nthreads);

    //HERE
    // Populate structs for threads
    for (size_t i=0 ; i < nthreads ; i++) {
      ApertureThreadArgs<T>::args[i].iPointBegin = 0+i*(nPositions/nthreads);
      ApertureThreadArgs<T>::args[i].iPointEnd   = (nPositions/nthreads)+i*(nPositions/nthreads);
      ApertureThreadArgs<T>::args[i].pos         = pos;
      ApertureThreadArgs<T>::args[i].field       = (*odata);
      ApertureThreadArgs<T>::args[i].k           = k;
      ApertureThreadArgs<T>::args[i].uxs         = uxs.get();
      ApertureThreadArgs<T>::args[i].uweights    = uweights.get();
      ApertureThreadArgs<T>::args[i].nDivU       = nDivW;
      ApertureThreadArgs<T>::args[i].vxs         = vxs.get();
      ApertureThreadArgs<T>::args[i].vweights    = vweights.get();
      ApertureThreadArgs<T>::args[i].nDivV       = nDivH;
      ApertureThreadArgs<T>::args[i].thread_id   = i;
      ApertureThreadArgs<T>::args[i].cpu_id      = ((int) i) % nproc;
      if (i==(nthreads-1))
        ApertureThreadArgs<T>::args[i].iPointEnd = nPositions;
    }

# ifdef HAVE_MQUEUE_H
    // Message queues
    struct mq_attr qattr;

    char buf[MSGMAX];
    unsigned int prio;

    // Create two message queues (could use socket_pair)
    CallErrReturn(ApertureQueue<T>::mqd_master = mq_open,
                  ("/jmh-master", O_RDWR | O_CREAT,
                   S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
                   NULL), EXIT_FAILURE);
    CallErrReturn(mq_close,(ApertureQueue<T>::mqd_master), EXIT_FAILURE);

    CallErrReturn(ApertureQueue<T>::mqd_client = mq_open,
                  ("/jmh-master", O_RDWR | O_CREAT,
                   S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
                   NULL), EXIT_FAILURE);
    CallErrReturn(mq_close,(ApertureQueue<T>::mqd_client), EXIT_FAILURE);

    /* Initiate launcher, note that the threadarg are global and updated
       in between calls to CalcThreaded */
    for (size_t i=0 ; i < nthreads ; i++) {
      launcher[i] =
        sps::pthread_launcher<Aperture<T>,&Aperture<T>::CalcThreadFunc>(this,
            &ApertureThreadArgs<T>::args[i]);
    }

    // Should only happen once
    if ((!ApertureQueue<T>::threads_initialized) && !(nthreads>N_MAX_THREADS)) {

      /* Clear message queues (and create if they don't exist), TODO: Do
         this in two stages */
      CallErrReturn(mq_clear,("/jmh-master"),EXIT_FAILURE);
      CallErrReturn(mq_clear,("/jmh-client"),EXIT_FAILURE);
      debug_print("Queues cleared:\n");

      /* Initiate threads */
      for (size_t i=0; i<nthreads; i++) {
        CallErr(pthread_create,
                (&threads[i], &attr,
                 sps::launch_member_function<sps::pthread_launcher<Aperture<T>,
                 &Aperture<T>::CalcThreadFunc> >,
                 &launcher[i]));
      }
      ApertureQueue<T>::threads_initialized = true;
    }

    // ERROR HERE
    /* Signal idling threads to run */
    CallErrReturn(ApertureQueue<T>::mqd_client = mq_open,
                  ("/jmh-client", O_WRONLY,
                   S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
                   NULL), false);

    for (size_t i=0 ; i < nthreads ; i++) {
      // Placed on queue with decreasing order of priority
      prio = i;
      if (mq_send(ApertureQueue<T>::mqd_client,"RUN",3,prio)==-1)
        perror ("mq_send()");
      debug_print("Send: RUN:\n");
    }
    mq_close(ApertureQueue<T>::mqd_client);

    /* Wait for threads to finish */
    CallErrReturn(ApertureQueue<T>::mqd_master = mq_open,
                  ("/jmh-master", O_RDONLY,
                   S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
                   NULL), false);

    mq_getattr(ApertureQueue<T>::mqd_master, &qattr);

    size_t n_to_go = nthreads;
    size_t n_idle  = 0;

    /* Wait for threads to finish */
    while (true) {
      debug_print("Master waiting\n");
      ssize_t len = mq_receive(ApertureQueue<T>::mqd_master, &buf[0], qattr.mq_msgsize, &prio);
      debug_print("Master received\n");
      buf[len] = '\0';
      debug_print("%s\n",buf);
      if (!strncmp(buf,"DONE",4)) {
        n_to_go--;
      } else if (!strncmp(buf,"READY",5)) {
        n_idle++;
      } else {
        debug_print("Unexpected message: %s\n",buf);
        retval = false;
        break;
      }
      if (n_to_go==0)
        break;
    }

    debug_print("We are done\n");

    // Are we really done???
    mq_close(ApertureQueue<T>::mqd_master);
# else
    /* Without message queues (slower) */
    for (size_t i=0; i<nthreads; i++) {
      launcher[i] =
        sps::pthread_launcher<Aperture<T>,&Aperture<T>::CalcThreaded>(this,
            &ApertureThreadArgs<T>::args[i]);
    }
    for (size_t i=0; i<nthreads; i++) {
#  if defined(HAVE_PTHREAD_H)
      CallErr(pthread_create,
              (&threads[i], &attr,
               sps::launch_member_function<sps::pthread_launcher<Aperture<T>,
               &Aperture<T>::CalcThreaded> >,
               &launcher[i]));
#  elif defined(_WIN32)
      threads[i] =
        (uintptr_t) _beginthreadex(NULL, 0,
                                   sps::launch_member_function<sps::pthread_launcher<Aperture<T>,
                                   &Aperture<T>::CalcThreaded> >,
                                   &launcher[i], 0, &threadID );
#  endif
    }

    for (size_t i = 0; i < nthreads; i++) {
#  if defined(HAVE_PTHREAD_H)
      CallErr(pthread_join,(threads[i],NULL));
#  elif defined(_WIN32)
      WaitForSingleObject((HANDLE) threads[i], INFINITE );
      //GetExitCodeThread(threads[i],&exitCode);
      /*
        The exit value specified in the ExitThread or TerminateThread function.
        The return value from the thread function.
        The exit value of the thread's process.
       */
#  endif
    }
    // Without message queues we destroy attributes
#  if defined(HAVE_PTHREAD_H)
    CallErr(pthread_attr_destroy,(&attr));
#  endif
# endif
#endif
    return retval;
  }
#endif

  // TODO: Move to fnm_calc.cpp
  template <class T>
  int Aperture<T>::CalcCwFocus(const T* pos, const size_t nPositions, const size_t nDim,
                               std::complex<T>** odata, size_t* nOutPositions)
  {

    assert(nDim == 3);
    if (nDim != 3) {
      *odata = NULL;
      *nOutPositions = 0;
      return -1;
    }
    *odata = (std::complex<T>*) malloc(nPositions*sizeof(std::complex<T>));
    memset(*odata, 0, nPositions*sizeof(std::complex<T>));
    *nOutPositions = nPositions;

    const T lambda = Aperture<T>::_sysparm.c / m_data->m_f0;

    const T k = T(M_2PI)/lambda;

    const size_t nElements = m_data->m_nelements;

    const T* apodizations = m_data->m_apodizations.get();

    T apodization;

    // Need weights and abcissa values

    size_t nDivW = Aperture<T>::_sysparm.nDivW;
    size_t nDivH = Aperture<T>::_sysparm.nDivH;

    // TODO: For FDTSD store these
    auto uweights = sps::deleted_aligned_array_create<T>(nDivW);
    auto vweights = sps::deleted_aligned_array_create<T>(nDivH);
    auto uxs      = sps::deleted_aligned_array_create<T>(nDivW);
    auto vxs      = sps::deleted_aligned_array_create<T>(nDivH);

    memset(uxs.get(),0,sizeof(T)*nDivW);
    memset(vxs.get(),0,sizeof(T)*nDivH);

    memset(uweights.get(),0,sizeof(T)*nDivW);
    memset(vweights.get(),0,sizeof(T)*nDivH);

    // Common weights and abcissa
    for (size_t i = 0 ; i < nDivW ; i++) {
      gl::GLNode qp = gl::GL(nDivW,i);
      uxs[i]        = T(qp.value);
      uweights[i]   = T(qp.weight);
    }

    for (size_t i = 0 ; i < nDivH ; i++) {
      gl::GLNode qp = gl::GL(nDivH, i);
      vxs[i]      = T(qp.value);
      vweights[i] = T(qp.weight);
    }

    auto& elements = m_data->m_elements;

    sps::point_t<T> projection;

    // TODO: Update to work for double

    // vectors, pos, output
    for (size_t iPoint = 0 ; iPoint < nPositions ; iPoint++) {

#ifdef FNM_DOUBLE_SUPPORT
      sps::point_t<T> point;
      memcpy(&point[0],&pos[iPoint*3],3*sizeof(T));
#else
      __m128 vec_point = _mm_set_ps(0.0f, float(pos[iPoint*3 + 2]),float(pos[iPoint*3 + 1]),float(pos[iPoint*3]));
#endif

      size_t iElement = iPoint % nElements;

      apodization = apodizations[iElement];

      std::complex<T> final = std::complex<T>(T(0.0),T(0.0));

      if (apodization != 0.0) {

        for (size_t jElement = 0 ; jElement < this->m_data->m_nsubelements ; jElement++) {

          const sps::element_t<T>& element = elements[iElement][jElement];

          std::complex<T> result;

#ifdef FNM_DOUBLE_SUPPORT
          sps::point_t<T> r2p = point - element.center;
          sps::point_t<T> hh_dir,hw_dir,normal;

          // Projection onto plane
          projection[0] = dot(hw_dir,r2p);
          projection[1] = dot(hh_dir,r2p);
          projection[2] = fabs(dot(normal,r2p));
#else
          assert(((uintptr_t)&element.center[0] & 0x0F) == 0 && "Data must be aligned");
          assert(((uintptr_t)&element.uvector[0] & 0x0F) == 0 && "Data must be aligned");
          assert(((uintptr_t)&element.vvector[0] & 0x0F) == 0 && "Data must be aligned");
          __m128 vec_r2p = _mm_sub_ps(vec_point, _mm_load_ps((float*)&element.center[0]));

          // Use vectors stored with elements (not working, if not set using ElementsSet)
          _mm_store_ss((float*)&projection[0],_mm_dp_ps(_mm_load_ps((float*)&element.uvector[0]), vec_r2p,0x71));
          _mm_store_ss((float*)&projection[1],_mm_dp_ps(_mm_load_ps((float*)&element.vvector[0]), vec_r2p,0x71));
          _mm_store_ss((float*)&projection[2],_mm_fabs_ps(_mm_dp_ps(_mm_load_ps((float*)&element.normal[0]), vec_r2p,0x71)));
#endif
          result = CalcHzFast<T>(element, projection, k,
                                 uxs.get(), uweights.get(), nDivW,
                                 vxs.get(), vweights.get(), nDivH);
          T real = result.real();
          T imag = result.imag();

          T carg = cos(m_data->m_phases[iElement]);
          T sarg = sin(m_data->m_phases[iElement]);
          final.real(final.real() + real*carg - imag*sarg);
          final.imag(final.imag() + real*sarg + imag*carg);
        }
      }
      (*odata)[iPoint].real(final.real());
      (*odata)[iPoint].imag(final.imag());
    }
    return 0;
  }

// How do we pass arguments to running thread
#ifdef HAVE_MQUEUE_H
  template <class T>
  void* Aperture<T>::CalcThreadFunc(void *ptarg)
  {

    char buf[MSGMAX];
    thread_arg* threadarg = NULL;
    //cpu_set_t set;
    struct mq_attr qattr;
    unsigned int prio;

    threadarg = (thread_arg*) ptarg;

    int thread_id = threadarg->thread_id;

    // Priority set to thread ID
    prio = thread_id;

    //setcpuid(threadarg->cpu_id);

    debug_print("Queue initialized: iPointBegin: %zu\n", threadarg->iPointBegin);
    mqd_t lmqd_master;
    mqd_t lmqd_client;
    CallErrReturn(lmqd_master = mq_open,
                  ("/jmh-master",O_WRONLY,
                   S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
                   NULL), NULL);

    // Post Message
    if (mq_send(lmqd_master,"READY",5,prio)==-1)
      perror ("mq_send(): READY");
    debug_print("Send: READY\n");

    int state = 0; // Ready
    void* thread_retval = NULL;

    CallErrReturn(lmqd_client = mq_open,
                  ("/jmh-client",O_RDONLY,
                   S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
                   NULL), NULL);

    mq_getattr (lmqd_client, &qattr);

    qattr.mq_flags &= ~O_NONBLOCK;
    mq_setattr(lmqd_client,&qattr,0);
    int rett = 0;

    //ssize_t mq_timedreceive

    // TODO: Receive only with a given priority
    while(mq_receive(lmqd_client, &buf[0], qattr.mq_msgsize, &prio) != -1) {
      debug_print("Message received\n");
      if (errno == EAGAIN) {
        debug_print("Apparently non-blocking\n");
        pthread_exit(NULL);
      }

      if (!strcmp(buf,"RESET"))
        state = 0;
      else if (!strncmp(buf,"RUN",3)) {
        state = 1;
      } else if (!strncmp(buf,"EXIT",4)) {
        state = 3;
        // Break while loop
        break;
      } else {
        printf("Unknown message: %s", buf);
      }

      switch (state) {
      case 0:
      case 1:
        debug_print("Call: CalcThreaded\n");

        // TESTME (use thread ID gettid())
        //threadarg = Aperture<float>::threadarg[

        thread_retval = this->CalcThreaded((void*)threadarg);
        state = 2;
      case 2:
        debug_print("Send: DONE\n");
        rett = mq_send(lmqd_master,"DONE", 4, prio);
        if (rett==-1) {
          printf("Error sending DONE\n");
          perror ("mq_send()");
        }
        state = 0;
        break; // break-out of switch only
      case 3:
        // Should exit
        break;
      default:
        break;
      }
    };

    // TODO: Should I close the master
    CallErrReturn(mq_close,(lmqd_client),NULL);
    CallErrReturn(mq_close,(lmqd_master),NULL);
    // We do not need to join, right?
    pthread_exit(thread_retval);
  }
#endif


#ifndef FNM_DOUBLE_SUPPORT
#if defined(HAVE_PTHREAD_H)
# ifdef HAVE_MQUEUE_H
  template <class T>
  void* Aperture<T>::CalcThreaded(void* ptarg) // Verify
# else
  template <class T>
  void* Aperture<T>::CalcThreaded(void* ptarg)
# endif
#else
  template <class T>
  unsigned int __stdcall Aperture<T>::CalcThreaded(void *ptarg)
#endif
  {
    fnm::thread_arg<T>* pThreadArg = reinterpret_cast<fnm::thread_arg<T>*>(ptarg);

    const T* uxs                 = pThreadArg->uxs;
    const T* vxs                 = pThreadArg->vxs;
    const T* uweights            = pThreadArg->uweights;
    const T* vweights            = pThreadArg->vweights;
    const T* pos                 = pThreadArg->pos;
    const size_t nDivW           = pThreadArg->nDivU;
    const size_t nDivH           = pThreadArg->nDivV;
    std::complex<T>* odata       = pThreadArg->field;

#ifndef HAVE_MQUEUE_H
# ifdef HAVE_THREAD
    setcpuid(pThreadArg->cpu_id);
# endif
#endif

    const size_t nElements = m_data->m_nelements;
    auto& elements = m_data->m_elements;
    T apodization;
    T k = pThreadArg->k;

    const T* apodizations = m_data->m_apodizations.get();

    ALIGN16_BEGIN sps::point_t<T> projection ALIGN16_END;

#ifdef _MSC_VER
    debug_print("iPointBegin: %Iu, iPointEnd: %Iu\n", pThreadArg->iPointBegin, pThreadArg->iPointEnd);
#else
    debug_print("iPointBegin: %zu, iPointEnd: %zu\n", pThreadArg->iPointBegin, pThreadArg->iPointEnd);
#endif

#if FNM_ENABLE_ATTENUATION
    T alpha = T(0.5) * T(100) / T(1e6) * m_data->m_f0 * T(0.1151);
#endif

    // vectors, pos, output
    for (size_t iPoint = pThreadArg->iPointBegin ; iPoint < pThreadArg->iPointEnd ; iPoint++) {

      __m128 vec_point =
        _mm_set_ps(0.0f,
                   float(pos[iPoint*3 + 2]),
                   float(pos[iPoint*3 + 1]),
                   float(pos[iPoint*3 + 0]));

      std::complex<T> final = std::complex<T>(T(0.0),T(0.0));

      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {

        apodization = apodizations[iElement];

        if (apodization != 0.0) {

          for (size_t jElement = 0 ; jElement <  m_data->m_nsubelements ; jElement++) {

            const auto& element = elements[iElement][jElement];
            std::complex<T> result;
            assert( ((uintptr_t)&element.center[0] & 0xF) == 0);
            assert( ((uintptr_t)&element.uvector[0] & 0xF) == 0);
            assert( ((uintptr_t)&element.vvector[0] & 0xF) == 0);
            __m128 vec_r2p = _mm_sub_ps(
                               vec_point,
                               _mm_load_ps((float*)&element.center[0]));
            _mm_store_ss((float*)&projection[0],
                         _mm_dp_ps(_mm_load_ps(&element.uvector[0]),
                                   vec_r2p,0x71));
            _mm_store_ss((float*)&projection[1],
                         _mm_dp_ps(_mm_load_ps(&element.vvector[0]),
                                   vec_r2p,0x71));
            _mm_store_ss((float*)&projection[2],
                         _mm_fabs_ps(_mm_dp_ps(_mm_load_ps(&element.normal[0]),
                                               vec_r2p,0x71)));

            // TODO: Un-roll by a factor of 4
            result = CalcHzFast<T>(element, projection, k,
                                   uxs, uweights, nDivW,
                                   vxs, vweights, nDivH);

            T real = apodization * result.real();
            T imag = apodization * result.imag();

            T carg = cos(m_data->m_phases[iElement]);
            T sarg = sin(m_data->m_phases[iElement]);

            // TODO: Fix attenuation (HERE). It is not working
#if FNM_ENABLE_ATTENUATION
            T dist = (T) _mm_cvtss_f32(_mm_sqrt_ps(_mm_dp_ps(vec_r2p,vec_r2p,0x71)));
            T factor = exp(-dist*alpha);
            final.real(final.real() + factor*(real*carg - imag*sarg));
            final.imag(final.imag() + factor*(real*sarg + imag*carg));
#else
            final.real(final.real() + real*carg - imag*sarg);
            final.imag(final.imag() + real*sarg + imag*carg);
#endif
          }
        }
      }
      odata[iPoint].real(final.real());
      odata[iPoint].imag(final.imag());
    }
    debug_print("Thread done\n");
#if HAVE_PTHREAD_H
# ifdef HAVE_MQUEUE_H
    return NULL;
# else
    pthread_exit(NULL);
# endif
#else
    return 0;
#endif
  }
#endif

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

#ifndef FNM_DOUBLE_SUPPORT
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

  // Specialize static variable
  template <class T>
  thread_arg<T> ApertureThreadArgs<T>::args[N_MAX_THREADS];

#ifdef HAVE_MQUEUE_H
  // Static variables
  template <class T>
  bool  ApertureQueue<T>::threads_initialized = false;

  template <class T>
  mqd_t ApertureQueue<T>::mqd_master = 0;

  template <class T>
  mqd_t ApertureQueue<T>::mqd_client = 0;
#endif

  //  template <class T>
  //  sysparm_t<T> Aperture<float>::_sysparm;

#ifdef HAVE_MQUEUE_H
  template class ApertureQueue<float>;
#endif

#ifdef FNM_DOUBLE_SUPPORT
# ifdef HAVE_MQUEUE_H
  template class ApertureQueue<double>;
# endif
  template class Aperture<double>;
#endif

  template class Aperture<float>;
}

// Template instantiation
#ifdef HAVE_MQUEUE_H
template class sps::pthread_launcher<fnm::Aperture<float>,&fnm::Aperture<float>::CalcThreadFunc>;
#endif

#ifndef FNM_DOUBLE_SUPPORT
template class sps::pthread_launcher<fnm::Aperture<float>,&fnm::Aperture<float>::CalcThreaded>;
#endif

#ifdef FNM_DOUBLE_SUPPORT
# ifdef HAVE_MQUEUE_H
template class sps::pthread_launcher<fnm::Aperture<double>,&fnm::Aperture<double>::CalcThreadFunc>;
# endif
#endif

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */

