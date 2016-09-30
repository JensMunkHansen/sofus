/**
 * @file   fnm.cpp
 * @author Jens Munk Hansen <jmh@jmhlaptop>
 * @date   Sun Jun 12 18:57:55 2016
 *
 * @brief  Source file for FNM utilities
 *
 *
 */

// TODO: Reuse results,
//       Limits independent of z-coordinate
//       Repeated integral
//       Consider complex f0 or using k
//       Prevent one thread to eat all messages (or send index of data)
//       Support arbitrary number of threads
//       prio must be pointer for mq_receive
//       Use static or unnamed namespace to use internal linkage

#include <fnm/fnm.hpp>
#include <fnm/config.h>
#include <fnm/FnmSIMD.hpp>
#include <fnm/fnm_data.hpp>
#include <fnm/fnm_calc.hpp>

#include <sps/smath.hpp>
#include <sps/sps_threads.hpp>
#include <sps/sps_mqueue.hpp>

#include <sps/mm_malloc.h>
#include <sps/cerr.h>
#include <sps/profiler.h>
#include <sps/extintrin.h>
#include <sps/trigintrin.h>
#include <sps/smath.hpp>
#include <sps/debug.h>

#include <gl/gl.hpp>

#include <string.h>
#include <memory>
#include <new>
#include <stdexcept>

#include <assert.h>

#ifdef HAVE_MQUEUE_H
# include <mqueue.h>
#endif

#ifdef _MSC_VER
#include <BaseTsd.h>
#endif

#ifdef GL_REF
# include "gauss_legendre.h"
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
    T* normals;
    T* uvectors;
    T* vvectors;
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
#endif

#if defined(HAVE_THREAD) && HAVE_PTHREAD_H
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

  template <class T>
  size_t Aperture<T>::nthreads = 4;

  ////////////////////////////////
  // Implementation of interface
  ////////////////////////////////
  template <class T>
  Aperture<T>::Aperture() : m_data(NULL)
  {
    m_data = (ApertureData<T>*) _mm_malloc(sizeof(ApertureData<T>),16);
    new (this->m_data) ApertureData<T>();
    initElements();
  }

  template <class T>
  Aperture<T>::Aperture(const size_t nElements, const T width, const T kerf, const T height) : Aperture<T>()
  {
    // TODO: Move to ApertureData
    m_data->m_nelements = nElements;
    m_data->m_nsubelements = 1;
    if (m_data->m_elements) {
      delete m_data->m_elements;
      m_data->m_elements = NULL;
    }
    m_data->m_elements =
      new std::vector< std::vector<element_t<T> > >(m_data->m_nelements,
          std::vector<element_t<T> >(m_data->m_nsubelements));

    if (m_data->m_apodizations) {
      _mm_free(m_data->m_apodizations);
      m_data->m_apodizations = NULL;
    }
    m_data->m_apodizations = (T*) _mm_malloc(m_data->m_nelements*sizeof(T),16);

    if (m_data->m_phases) {
      _mm_free(m_data->m_phases);
      m_data->m_phases = NULL;
    }
    m_data->m_phases = (T*) _mm_malloc(m_data->m_nelements*sizeof(T),16);
    if (m_data->m_phases) {
      memset(m_data->m_phases,0,m_data->m_nelements*sizeof(T));
    }

    if (m_data->m_pos) {
      _mm_free(m_data->m_pos);
      m_data->m_pos = NULL;
    }
    m_data->m_pos = (sps::point_t<T>*) _mm_malloc(nElements*sizeof(sps::point_t<T>),16);

    T wx = T(nElements-1)/2;

    T pitch = width + kerf;

    std::vector<std::vector<element_t<T> > >& elements = *m_data->m_elements;

    for (size_t i = 0 ; i < m_data->m_nelements ; i++) {
      m_data->m_apodizations[i] = T(1.0);
      elements[i][0].center[0] = (T(i) - wx)*pitch;
      elements[i][0].center[1] = T(0.0);
      elements[i][0].center[2] = T(0.0);

      m_data->m_pos[i][0] = elements[i][0].center[0];
      m_data->m_pos[i][1] = elements[i][0].center[1];
      m_data->m_pos[i][2] = elements[i][0].center[2];

      memset((void*)&elements[i][0].euler,0,sizeof(sps::euler_t<T>));
      elements[i][0].hh           = height/2;
      printf("width: %f\n",width);
      elements[i][0].hw           = width/2;
      elements[i][0].euler.alpha  = T(0.0);
      elements[i][0].euler.beta   = T(0.0);
      elements[i][0].euler.gamma  = T(0.0);
    }
    initElements();
  }

  template <class T>
  Aperture<T>::~Aperture()
  {
    debug_print("~Aperture()\n");
    if (m_data) {
      _mm_free(m_data);
      m_data = NULL;
    }
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
# if defined(HAVE_PTHREAD_H)
      CallErr(pthread_attr_destroy,(&attr));
# endif
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
    const sps::point_t<T>* positions = m_data->m_pos;

    // The function allocates
    T* arr = (T*) malloc(arrSize);
    if (arr && positions) {
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
  void Aperture<T>::PositionsSet(const T* pos, const size_t nPositions, const size_t nDim)
  {
    const size_t _nElements          = m_data->m_nelements;
    sps::point_t<T>* positions = m_data->m_pos;
    if (_nElements*3 == nPositions * nDim) {
      for (size_t i = 0 ; i < _nElements ; i++) {
        positions[i][0] = pos[i*3 + 0];
        positions[i][0] = pos[i*3 + 1];
        positions[i][0] = pos[i*3 + 2];
      }
    }
  }

  template <class T>
  void Aperture<T>::FocusSet(const T iFocus[3])
  {
    if ((iFocus[0] != m_data->m_focus[0]) || (iFocus[0] != m_data->m_focus[0]) || (iFocus[0] != m_data->m_focus[0])) {
      m_data->_focus_validity = FocusingType::FocusingTypeCount;
    }
    memcpy(&m_data->m_focus[0],&iFocus[0],3*sizeof(T));
  }

  template <class T>
  void Aperture<T>::FocusUpdate()
  {

    size_t nElements = m_data->m_nelements;

    if (m_data->m_focus_type == FocusingType::Pythagorean) {
      std::unique_ptr<T[]> delays = std::unique_ptr<T[]>(new T[nElements]);
      T maxDelay = T(0.0);
      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
        delays[iElement] = norm(m_data->m_pos[iElement] - m_data->m_focus) / Aperture<T>::_sysparm.c;
        maxDelay = std::max<T>(maxDelay,delays[iElement]);
      }
      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
        m_data->m_phases[iElement] = T(M_2PI) * (delays[iElement] - maxDelay) * m_data->m_f0;
      }
    } else if (m_data->m_focus_type == FocusingType::Rayleigh) {
      // Right for CW only
      size_t nPositions = nElements;
      std::unique_ptr<T[]> positions = std::unique_ptr<T[]>(new T[3*nPositions]);
      for (size_t iPosition = 0 ; iPosition < nPositions ; iPosition++) {
        positions[3*iPosition]   = m_data->m_focus[0];
        positions[3*iPosition+1] = m_data->m_focus[1];
        positions[3*iPosition+2] = m_data->m_focus[2];
      }

      std::complex<T>* pFieldValues = NULL;
      size_t nFieldValues;
      this->CalcCwFocus(positions.get(),nPositions,3,&pFieldValues,&nFieldValues);

      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
        m_data->m_phases[iElement] = - std::arg<T>(pFieldValues[iElement]);
      }
      free(pFieldValues);
    }
  }

  template <class T>
  void Aperture<T>::PhasesGet(T** phases, size_t* nPhases)
  {

    size_t nElements = m_data->m_nelements;
    *phases = (T*) malloc(nElements*sizeof(T));
    if (nElements > 0) {
      memcpy(*phases,m_data->m_phases,nElements*sizeof(T));
      *nPhases = nElements;
    } else {
      *nPhases = 0;
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
    m_data->m_f0 = f0;
  }

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
  void Aperture<T>::initElements()
  {

    // Note: The height dh is used in the y-direction (if euler = (0,0,0))

    size_t iElement, jElement;
    size_t i_xyz;

    const size_t nElements       = m_data->m_nelements;
    const size_t nSubElements    = m_data->m_nsubelements;

    // Not used right now
    const size_t nSubH           = 1;
    const size_t nSubW           = 1;
    const size_t nSubSubElements = 1;
    size_t i_subh, i_subw;

    // Elements
    std::vector<std::vector<element_t<T> > >& elements       = *m_data->m_elements;

    sps::point_t<T> h_dir, w_dir, normal;

    // Delta height and width
    T dh, dw;

    // Lower left position reference
    T pos_ll;

    if (m_data->m_rectangles) {
      delete [] m_data->m_rectangles;
      m_data->m_rectangles = NULL;
    }

    // Should be max(nSubElements) * nSubH * nSubW * nElements
    m_data->m_rectangles =
      new sps::rect_t<T>[nElements*nSubElements*nSubSubElements];

    for (iElement=0 ; iElement < nElements ; iElement++) {
      // Set apodizations
      m_data->m_apodizations[iElement] = T(1.0);

      // TODO: Keep this (introduce possiblity to set positions, without this they can be uninitialized)

      // Center element (reference for delays)
      if (nSubElements % 2 == 1) {
        // Odd number of elements
        m_data->m_pos[iElement] = elements[iElement][nSubElements/2].center;
      } else {
        // Even number of elements (TODO: On an edge)
        m_data->m_pos[iElement] = T(0.5) * (elements[iElement][nSubElements-1].center + elements[iElement][0].center);
      }

      for (jElement=0 ; jElement < nSubElements ; jElement++) {
        const element_t<T>& element = elements[iElement][jElement];

        printf("half width: %f\n",element.hw);
        
        // Hack (not working for double anyway)
        SPS_UNREFERENCED_PARAMETER(normal);
        sps::euler_t<float> euler;
        memcpy((void*)&euler,(void*)&element.euler,sizeof(sps::euler_t<float>));

        __m128 uvector, vvector, normalvector;

        basis_vectors_ps((float*)&uvector, (float*)&vvector, (float*)&normalvector, euler);
        _mm_store_ps((float*)&w_dir[0],uvector);
        _mm_store_ps((float*)&h_dir[0],vvector);

        dh = 2*element.hh/nSubH;
        dw = 2*element.hw/nSubW;

        // Maximum width or height used for memory allocation
        // Loop over x,y,z
        for (i_xyz=0 ; i_xyz < 3 ; i_xyz++) {

          // Corner position of rectangle (lower left)
          pos_ll =
            element.center[i_xyz]
            - element.hh*h_dir[i_xyz]
            - element.hw*w_dir[i_xyz];

          // Compute coordinates for corners of mathematical elements
          // TODO: Consider using lists of lists of array[nSubH][nSubW]
          for (i_subh=0 ; i_subh < nSubH ; i_subh++) {
            for (i_subw=0 ; i_subw < nSubW ; i_subw++) {
              m_data->m_rectangles[iElement*nSubSubElements*nSubElements +
                                   jElement*nSubSubElements +
                                   nSubH*i_subw + i_subh][0][i_xyz] =
                                     pos_ll + i_subh*dh*h_dir[i_xyz]+i_subw*dw*w_dir[i_xyz];

              m_data->m_rectangles[iElement*nSubSubElements*nSubElements +
                                   jElement*nSubSubElements +
                                   nSubH*i_subw + i_subh][1][i_xyz] =
                                     pos_ll + i_subh*dh*h_dir[i_xyz]+(i_subw+1)*dw*w_dir[i_xyz];

              m_data->m_rectangles[iElement*nSubSubElements*nSubElements +
                                   jElement*nSubSubElements +
                                   nSubH*i_subw + i_subh][2][i_xyz] =
                                     pos_ll + (i_subh+1)*dh*h_dir[i_xyz]+i_subw*dw*w_dir[i_xyz];

              m_data->m_rectangles[iElement*nSubSubElements*nSubElements +
                                   jElement*nSubSubElements +
                                   nSubH*i_subw + i_subh][3][i_xyz] =
                                     pos_ll + (i_subh+1)*dh*h_dir[i_xyz]+(i_subw+1)*dw*w_dir[i_xyz];
            } // for i_subw
          } // for i_subh
        } // for i_xyz
      } // jElement
    } // iElement
  }

  template <class T>
  void Aperture<T>::RectanglesGet(T** out, size_t* nElements,
                                  size_t* nSubElements, size_t* nParams) const
  {

    // Needed to avoid temporaries (if a view is returned). We
    // allocate, so this is not a temporary
    static size_t _nCornerCoordinates = 12;

    const size_t _nElements       = m_data->m_nelements;
    const size_t _nSubElements    = m_data->m_nsubelements; // Length of vectors
    const size_t arrSize          = _nElements * _nSubElements * _nCornerCoordinates * sizeof(T);

    // The function allocates
    T* arr = (T*) malloc(arrSize);
    if (arr) {
      memset(arr,0,arrSize);

      const sps::rect_t<T>* rectangles = m_data->m_rectangles;

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
      *nParams      = _nCornerCoordinates;

      *out          = arr;
    }
  }

  template <class T>
  void Aperture<T>::ElementsGet(T** out, size_t* nElements,
                                size_t* nParams) const
  {
    printf("half width: %f\n", (*m_data->m_elements)[0][0].hw);
    printf("half height: %f\n", (*m_data->m_elements)[0][0].hh);

    // TODO: Use common function for ElementsGet and SubElementsGet

    const size_t nElePosParams = Aperture<T>::nElementPosParameters;

    const size_t _nElements             = m_data->m_nelements;
    const size_t nSubElementsPerElement = 1;
    const size_t arrSize                = _nElements * nSubElementsPerElement * nElePosParams * sizeof(T);

    // The function allocates
    T* arr = (T*) malloc(arrSize);
    if (arr) {
      memset(arr,0,arrSize);

      std::vector<std::vector<element_t<T> > >& elements = *m_data->m_elements;

      for (size_t iElement = 0 ; iElement < _nElements ; iElement++) {
        arr[iElement * nSubElementsPerElement * nElePosParams] = elements[iElement][0].hw;

        arr[iElement * nSubElementsPerElement * nElePosParams + 1] = elements[iElement][0].hh;

        printf("half width: %f\n", elements[iElement][0].hw);
        printf("half height: %f\n", elements[iElement][0].hh);
        
        memcpy(&arr[iElement * nSubElementsPerElement * nElePosParams + 2],
               &elements[iElement][0].center[0],
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

      std::vector<std::vector<element_t<T> > >& elements = *m_data->m_elements;

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
      *nElements    = m_data->m_nelements;
      *nSubElements = m_data->m_nsubelements;
      *nParams      = Aperture<T>::nElementPosParameters; // No need to introduce static variable

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

    if (nDim != Aperture<T>::nElementPosParameters) {
      retval = false;
    } else {
      if ((nElements != m_data->m_nelements) || (m_data->m_nsubelements != nSubElementsPerElement)) {

        if (m_data->m_elements)
          delete m_data->m_elements;

        m_data->m_elements     = NULL;
        m_data->m_nelements    = nElements;
        m_data->m_nsubelements = nSubElementsPerElement;

        m_data->m_elements     =
          new std::vector< std::vector<element_t<T> > >(nElements,
              std::vector<element_t<T> >(nSubElementsPerElement));

        if (m_data->m_apodizations) {
          _mm_free(m_data->m_apodizations);
          m_data->m_apodizations = (T*) _mm_malloc(nElements*sizeof(T),16);
        }

        if (m_data->m_pos) {
          _mm_free(m_data->m_pos);
          m_data->m_pos = (sps::point_t<T>*) _mm_malloc(nElements*sizeof(sps::point_t<T>),16);
        }
      }

      std::vector<std::vector<element_t<T> > >& elements = *m_data->m_elements;

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
                 sizeof(sps::point_t<T>));
          elements[iElement][jElement].euler.alpha =
            pos[iElement * nSubElementsPerElement * nElePosParams + jElement*nElePosParams+5];
          elements[iElement][jElement].euler.beta =
            pos[iElement * nSubElementsPerElement * nElePosParams + jElement*nElePosParams+6];
          elements[iElement][jElement].euler.gamma =
            pos[iElement * nSubElementsPerElement * nElePosParams + jElement*nElePosParams+7];
        }
      }
      initElements();
    }
    return retval;
  }

// Calculations

  template <class T>
  void Aperture<T>::CalcCwField(const T* pos, const size_t nPositions, const size_t nDim,
                                std::complex<T>** odata, size_t* nOutPositions)
  {

    assert(nDim == 3);
    if (nDim != 3) {
      *odata = NULL;
      *nOutPositions = 0;
      return;
    }
    size_t nScatterers = nPositions;
    *odata = (std::complex<T>*) malloc(nScatterers*sizeof(std::complex<T>));
    *nOutPositions = nScatterers;

    if (m_data->_focus_validity != m_data->m_focus_type) {
      // TODO: Choose between time-of-flight or phase adjustment
      this->FocusUpdate();
      m_data->_focus_validity = m_data->m_focus_type;
    }

    fnm::CalcCwField<T>(*this->m_data,
                        pos, nPositions,
                        odata);
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

    size_t nScatterers = nPositions;
    *odata = (std::complex<T>*) malloc(nScatterers*sizeof(std::complex<T>));
    *nOutPositions = nScatterers;

    if (m_data->_focus_validity != m_data->m_focus_type) {
      // TODO: Choose between time-of-flight or phase adjustment
      this->FocusUpdate();
      m_data->_focus_validity = m_data->m_focus_type;
    }

    int retval = fnm::CalcCwFieldRef<T>(*this->m_data,
                                        pos, nPositions,
                                        odata);
    return retval;
  }


// Sub-optimal integration range (Old function)
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

    if (m_data->_focus_validity != m_data->m_focus_type) {
      // TODO: Choose between time-of-flight or phase adjustment
      this->FocusUpdate();
      m_data->_focus_validity = m_data->m_focus_type;
    }

    const T lambda = Aperture<T>::_sysparm.c / m_data->m_f0;

    const T k = T(M_2PI)/lambda;

    const size_t nElements = m_data->m_nelements;

    T* apodizations = m_data->m_apodizations;

    T apodization;

    sps::point_t<T> hh_dir, hw_dir, normal;

    // Need weights and abcissa values

    size_t nDivW = Aperture<T>::_sysparm.nDivW;
    size_t nDivH = Aperture<T>::_sysparm.nDivH;

    auto del = [](T* p) {
      _mm_free(p);
    };
    std::unique_ptr<T[], decltype(del)> uweights((T*)_mm_malloc(sizeof(T)*nDivW,16), del);
    std::unique_ptr<T[], decltype(del)> vweights((T*)_mm_malloc(sizeof(T)*nDivH,16), del);
    std::unique_ptr<T[], decltype(del)> uxs((T*)_mm_malloc(sizeof(T)*nDivW,16), del);
    std::unique_ptr<T[], decltype(del)> vxs((T*)_mm_malloc(sizeof(T)*nDivH,16), del);

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

    std::vector<std::vector<element_t<T> > >& elements = *m_data->m_elements;

    sps::point_t<T> projection;

    auto del128 = [](__m128* p) {
      _mm_free(p);
    };
    std::unique_ptr<__m128[],decltype(del128)>  vec_normals((__m128*) _mm_malloc(nElements*sizeof(__m128),16),del128);
    std::unique_ptr<__m128[],decltype(del128)> vec_uvectors((__m128*) _mm_malloc(nElements*sizeof(__m128),16),del128);
    std::unique_ptr<__m128[],decltype(del128)> vec_vvectors((__m128*) _mm_malloc(nElements*sizeof(__m128),16),del128);

    for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
      const element_t<T>& element = elements[iElement][0];

      // Get basis vectors (TODO: Use SSE4)
      basis_vectors(hw_dir, element.euler, 0);
      basis_vectors(hh_dir, element.euler, 1);
      basis_vectors(normal, element.euler, 2);

      vec_normals[iElement]  = _mm_set_ps(0.0f,(float)normal[2],(float)normal[1],(float)normal[0]);
      vec_uvectors[iElement] = _mm_set_ps(0.0f,(float)hw_dir[2], (float)hw_dir[1], (float)hw_dir[0]);
      vec_vvectors[iElement] = _mm_set_ps(0.0f,(float)hh_dir[2], (float)hh_dir[1], (float)hh_dir[0]);
    }

    for (size_t iPoint = 0 ; iPoint < nPositions ; iPoint++) {

      __m128 vec_point = _mm_set_ps(0.0f, float(pos[iPoint*3 + 2]),float(pos[iPoint*3 + 1]),float(pos[iPoint*3]));

      std::complex<T> final = std::complex<T>(T(0.0),T(0.0));

      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {

        apodization = apodizations[iElement];

        if (apodization != 0.0) {

          for (size_t jElement = 0 ; jElement < elements[iElement].size() ; jElement++) {

            const element_t<T>& element = elements[iElement][jElement];

            std::complex<T> result;

#ifdef _WIN32
            __m128 vec_r2p = _mm_sub_ps(vec_point, _mm_loadu_ps((float*)&element.center[0]));
#else
            //_mm_stream_load_si128, _mm_stream_ps(float * p , __m128 a ); stores in *p
            __m128 vec_r2p = _mm_sub_ps(vec_point, _mm_load_ps((float*)&element.center[0]));
#endif
            _mm_store_ss((float*)&projection[0],_mm_dp_ps(vec_uvectors[iElement], vec_r2p,0x71));
            _mm_store_ss((float*)&projection[1],_mm_dp_ps(vec_vvectors[iElement], vec_r2p,0x71));
            _mm_store_ss((float*)&projection[2],_mm_fabs_ps(_mm_dp_ps(vec_normals[iElement], vec_r2p,0x71)));

            result = CalcHzAll<T>(element, projection, k,
                                  uxs.get(), uweights.get(), nDivW,
                                  vxs.get(), vweights.get(), nDivH);
            final = final + result * exp(std::complex<T>(0,m_data->m_phases[iElement]));
#if 0
            T real = result.real();
            T imag = result.imag();

            T carg = cos(m_data->m_phases[iElement]);
            T sarg = sin(m_data->m_phases[iElement]);
            final.real(final.real() + real*carg - imag*sarg);
            final.imag(final.imag() + real*sarg + imag*carg);
#endif
          }
        }
      }
      (*odata)[iPoint].real(final.real());
      (*odata)[iPoint].imag(final.imag());
    }
    return 0;
  }

// Reduced integral and fast (the one to use)
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

    if (m_data->_focus_validity != m_data->m_focus_type) {
      // TODO: Choose between time-of-flight or phase adjustment
      this->FocusUpdate();
      m_data->_focus_validity = m_data->m_focus_type;
    }

    const T lambda = Aperture<T>::_sysparm.c / m_data->m_f0;
    const T k = T(M_2PI)/lambda;

    size_t nDivW = Aperture<T>::_sysparm.nDivW;
    size_t nDivH = Aperture<T>::_sysparm.nDivH;

    auto del = [](T* p) {
      _mm_free(p);
    };
    std::unique_ptr<T[], decltype(del)> uweights((T*)_mm_malloc(sizeof(T)*nDivW,16), del);
    std::unique_ptr<T[], decltype(del)> vweights((T*)_mm_malloc(sizeof(T)*nDivH,16), del);
    std::unique_ptr<T[], decltype(del)> uxs((T*)_mm_malloc(sizeof(T)*nDivW,16), del);
    std::unique_ptr<T[], decltype(del)> vxs((T*)_mm_malloc(sizeof(T)*nDivH,16), del);

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

    // Common basis vectors
    const size_t nElements = m_data->m_nelements;
    std::vector<std::vector<element_t<T> > >& elements = *m_data->m_elements;

    auto del128 = [](__m128* p) {
      _mm_free(p);
    };
    std::unique_ptr<__m128[],decltype(del128)>  vec_normals((__m128*) _mm_malloc(nElements*sizeof(__m128),16),del128);
    std::unique_ptr<__m128[],decltype(del128)> vec_uvectors((__m128*) _mm_malloc(nElements*sizeof(__m128),16),del128);
    std::unique_ptr<__m128[],decltype(del128)> vec_vvectors((__m128*) _mm_malloc(nElements*sizeof(__m128),16),del128);

    for (size_t iElement = 0 ; iElement < nElements ; iElement++) {

      const element_t<T>& element = elements[iElement][0];

      // Get basis vectors
#if ACCURATE_TRIGONOMETRICS
      sps::point_t<T> hh_dir,hw_dir,normal;
      basis_vectors(hw_dir, element.euler, 0);
      basis_vectors(hh_dir, element.euler, 1);
      basis_vectors(normal, element.euler, 2);
      vec_normals[iElement]  = _mm_set_ps(0.0f,(float)normal[2], (float)normal[1], (float)normal[0]);
      vec_uvectors[iElement] = _mm_set_ps(0.0f,(float)hw_dir[2], (float)hw_dir[1], (float)hw_dir[0]);
      vec_vvectors[iElement] = _mm_set_ps(0.0f,(float)hh_dir[2], (float)hh_dir[1], (float)hh_dir[0]);
#else
      // Hack (not working for double anyway)
      sps::euler_t<float> euler;
      memcpy((void*)&euler,(void*)&element.euler,sizeof(sps::euler_t<float>));
      basis_vectors_ps((float*)&vec_uvectors[iElement],
                       (float*)&vec_vvectors[iElement],
                       (float*)&vec_normals[iElement], euler);
#endif
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
    threadarg.nDivU       = nDivW;
    threadarg.vxs         = vxs.get();
    threadarg.vweights    = vweights.get();
    // Unit vectors
    threadarg.normals     = (T*) vec_normals.get();
    threadarg.uvectors    = (T*) vec_uvectors.get();
    threadarg.vvectors    = (T*) vec_vvectors.get();
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
      ApertureThreadArgs<T>::args[i].normals     = (T*) vec_normals.get();
      ApertureThreadArgs<T>::args[i].uvectors    = (T*) vec_uvectors.get();
      ApertureThreadArgs<T>::args[i].vvectors    = (T*) vec_vvectors.get();
      ApertureThreadArgs<T>::args[i].thread_id   = i;
      ApertureThreadArgs<T>::args[i].cpu_id      = ((int) i) % nproc;
      if (i==(nthreads-1))
        ApertureThreadArgs<T>::args[i].iPointEnd = nPositions;
    }

# ifdef HAVE_MQUEUE_H
    // Mesage queues
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

    T* apodizations = m_data->m_apodizations;

    T apodization;

    // Need weights and abcissa values

    size_t nDivW = Aperture<T>::_sysparm.nDivW;
    size_t nDivH = Aperture<T>::_sysparm.nDivH;

    auto del = [](T* p) {
      _mm_free(p);
    };
    std::unique_ptr<T[], decltype(del)> uweights((T*)_mm_malloc(sizeof(T)*nDivW,16), del);
    std::unique_ptr<T[], decltype(del)> vweights((T*)_mm_malloc(sizeof(T)*nDivH,16), del);
    std::unique_ptr<T[], decltype(del)> uxs((T*)_mm_malloc(sizeof(T)*nDivW,16), del);
    std::unique_ptr<T[], decltype(del)> vxs((T*)_mm_malloc(sizeof(T)*nDivH,16), del);

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

    sps::point_t<T> hh_dir, hw_dir, normal;

    std::vector<std::vector<element_t<T> > >& elements = *m_data->m_elements;

    sps::point_t<T> projection;

    // Hack
    auto del128 = [](__m128* p) {
      _mm_free(p);
    };
    std::unique_ptr<__m128[],decltype(del128)>  vec_normals((__m128*) _mm_malloc(nElements*sizeof(__m128),16),del128);
    std::unique_ptr<__m128[],decltype(del128)> vec_uvectors((__m128*) _mm_malloc(nElements*sizeof(__m128),16),del128);
    std::unique_ptr<__m128[],decltype(del128)> vec_vvectors((__m128*) _mm_malloc(nElements*sizeof(__m128),16),del128);

    // TODO: Support 2D
    for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
      const element_t<T>& element = elements[iElement][0];

      basis_vectors(hw_dir, element.euler, 0);
      basis_vectors(hh_dir, element.euler, 1);
      basis_vectors(normal, element.euler, 2);

      vec_normals[iElement]  = _mm_set_ps(0.0f,(float)normal[2],(float)normal[1],(float)normal[0]);
      vec_uvectors[iElement] = _mm_set_ps(0.0f,(float)hw_dir[2], (float)hw_dir[1], (float)hw_dir[0]);
      vec_vvectors[iElement] = _mm_set_ps(0.0f,(float)hh_dir[2], (float)hh_dir[1], (float)hh_dir[0]);
    }

    // vectors, pos, output
    for (size_t iPoint = 0 ; iPoint < nPositions ; iPoint++) {

      __m128 vec_point = _mm_set_ps(0.0f, float(pos[iPoint*3 + 2]),float(pos[iPoint*3 + 1]),float(pos[iPoint*3]));

      std::complex<T> final = std::complex<T>(T(0.0),T(0.0));

      size_t iElement = iPoint % nElements;

      apodization = apodizations[iElement];

      if (apodization != 0.0) {

        for (size_t jElement = 0 ; jElement < elements[iElement].size() ; jElement++) {

          const element_t<T>& element = elements[iElement][jElement];

          std::complex<T> result;
#ifdef _WIN32
          __m128 vec_r2p = _mm_sub_ps(vec_point, _mm_loadu_ps((float*)&element.center[0]));
#else
          __m128 vec_r2p = _mm_sub_ps(vec_point, _mm_load_ps((float*)&element.center[0]));
#endif
          _mm_store_ss((float*)&projection[0],_mm_dp_ps(vec_uvectors[iElement], vec_r2p,0x71));
          _mm_store_ss((float*)&projection[1],_mm_dp_ps(vec_vvectors[iElement], vec_r2p,0x71));
          _mm_store_ss((float*)&projection[2],_mm_fabs_ps(_mm_dp_ps(vec_normals[iElement], vec_r2p,0x71)));

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
    const __m128* vec_uvectors   = (__m128*) pThreadArg->uvectors;
    const __m128* vec_vvectors   = (__m128*) pThreadArg->vvectors;
    const __m128* vec_nvectors   = (__m128*) pThreadArg->normals;
    const size_t nDivW           = pThreadArg->nDivU;
    const size_t nDivH           = pThreadArg->nDivV;
    std::complex<T>* odata       = pThreadArg->field;

#ifndef HAVE_MQUEUE_H
# ifdef HAVE_THREAD
    setcpuid(pThreadArg->cpu_id);
# endif
#endif

    const size_t nElements = m_data->m_nelements;

    std::vector<std::vector<element_t<T> > >& elements = *m_data->m_elements;

    T apodization;
    T k = pThreadArg->k;

    T* apodizations = m_data->m_apodizations;

    sps::point_t<T> projection;
#ifdef _MSC_VER
    debug_print("iPointBegin: %Iu, iPointEnd: %Iu\n", pThreadArg->iPointBegin, pThreadArg->iPointEnd);
#else
    debug_print("iPointBegin: %zu, iPointEnd: %zu\n", pThreadArg->iPointBegin, pThreadArg->iPointEnd);
#endif
    // vectors, pos, output
    for (size_t iPoint = pThreadArg->iPointBegin ; iPoint < pThreadArg->iPointEnd ; iPoint++) {

      __m128 vec_point =
        _mm_set_ps(0.0f,
                   float(pos[iPoint*3 + 2]),
                   float(pos[iPoint*3 + 1]),
                   float(pos[iPoint*3]));

      std::complex<T> final = std::complex<T>(T(0.0),T(0.0));

      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {

        apodization = apodizations[iElement];

        if (apodization != 0.0) {

          for (size_t jElement = 0 ; jElement < elements[iElement].size() ; jElement++) {

            const element_t<T>& element = elements[iElement][jElement];

            std::complex<T> result;
#ifdef _WIN32
            __m128 vec_r2p = _mm_sub_ps(
                               vec_point,
                               _mm_loadu_ps((float*)&element.center[0]));
#else
            __m128 vec_r2p = _mm_sub_ps(
                               vec_point,
                               _mm_load_ps((float*)&element.center[0]));
#endif
            _mm_store_ss((float*)&projection[0],
                         _mm_dp_ps(vec_uvectors[iElement],
                                   vec_r2p,0x71));
            _mm_store_ss((float*)&projection[1],
                         _mm_dp_ps(vec_vvectors[iElement],
                                   vec_r2p,0x71));
            _mm_store_ss((float*)&projection[2],
                         _mm_fabs_ps(_mm_dp_ps(vec_nvectors[iElement],
                                               vec_r2p,0x71)));

            //
            result = CalcHzFast<T>(element, projection, k,
                                   uxs, uweights, nDivW,
                                   vxs, vweights, nDivH);
            T real = result.real();
            T imag = result.imag();

            T carg = cos(m_data->m_phases[iElement]);
            T sarg = sin(m_data->m_phases[iElement]);
            final.real(final.real() + real*carg - imag*sarg);
            final.imag(final.imag() + real*sarg + imag*carg);
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

  template <class T>
  void Aperture<T>::ManagedAllocation(std::complex<T>** outTest, size_t* nOutTest)
  {
    const size_t _nData = 100000000;
    *nOutTest = _nData;
    *outTest = (std::complex<T>*) malloc(_nData*sizeof(std::complex<T>));
  }


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

#ifdef __GNUG__
# if 0
  template <class T>
  struct Aperture<T>::thread_arg Aperture<T>::threadarg[N_MAX_THREADS];
# endif
#elif defined(_WIN32)
  // Explicit instantiations
# if 0
  template struct Aperture<float>::thread_arg;
  template struct Aperture<double>::thread_arg;

  Aperture<float>::thread_arg Aperture<float>::threadarg[N_MAX_THREADS];
  Aperture<double>::thread_arg Aperture<double>::threadarg[N_MAX_THREADS];
# endif
#endif

#ifdef HAVE_MQUEUE_H
  template class ApertureQueue<float>;
  template class ApertureQueue<double>;
#endif
  template class Aperture<float>;
  template class Aperture<double>;
}


// Template instantiation
#ifdef HAVE_MQUEUE_H
template class sps::pthread_launcher<fnm::Aperture<float>,&fnm::Aperture<float>::CalcThreadFunc>;
template class sps::pthread_launcher<fnm::Aperture<double>,&fnm::Aperture<double>::CalcThreadFunc>;
#endif

template class sps::pthread_launcher<fnm::Aperture<float>,&fnm::Aperture<float>::CalcThreaded>;
template class sps::pthread_launcher<fnm::Aperture<double>,&fnm::Aperture<double>::CalcThreaded>;



/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */

