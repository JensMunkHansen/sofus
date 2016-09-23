// TODO: Reuse results,
//       Limits independent of z-coordinate
//       Repeated integral
//       Consider complex f0 or using k

#include <fnm/fnm.hpp>
#include <fnm/config.h>
#include <fnm/FnmSIMD.hpp>
#include <fnm/fnm_data.hpp>
#include <sps/sps_threads.hpp>

#include <sps/mm_malloc.h>
#include <sps/cerr.h>
#include <sps/profiler.h>
#include <sps/extintrin.h>
#include <sps/trigintrin.h>
#include <sps/smath.hpp>

#include <string.h>
#include <memory>
#include <new>
#include <assert.h>

#include "fastgl.hpp"

#ifdef GL_REF
# include "gauss_legendre.h"
#endif

namespace fnm {

#if defined(C99) && !defined(__STRICT_ANSI__)
  template <class T>
  sysparm_t<T> Aperture<T>::_sysparm =
    (sysparm_t<T>){.c = 1500.0, .nDivW = 16, .DivH = 16};
#else
  template <class T>
  sysparm_t<T> Aperture<T>::_sysparm = {T(1500.0), size_t(16), size_t(16)};
#endif

  template <class T>
  size_t Aperture<T>::nthreads = 4;

#ifdef HAVE_PTHREAD_H
  template <class T>
  pthread_t Aperture<T>::threads[N_MAX_THREADS];
  template <class T>
  pthread_attr_t Aperture<T>::attr = {{0}};
#endif

  
  template <class T>
  Aperture<T>::Aperture() : m_data(NULL) {
    m_data = (ApertureData<T>*) _mm_malloc(sizeof(ApertureData<T>),16);
    new (this->m_data) ApertureData<T>();
    initElements();
  }

  template <class T>
  Aperture<T>::Aperture(const size_t nElements, const T width, const T kerf, const T height) : Aperture<T>() {
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

    memset(m_data->m_phases,0,m_data->m_nelements*sizeof(T));

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
      elements[i][0].hw           = width/2;
      elements[i][0].euler.alpha  = T(0.0);
      elements[i][0].euler.beta   = T(0.0);
      elements[i][0].euler.gamma  = T(0.0);
    }
    initElements();
  }
  
  template <class T>
  Aperture<T>::~Aperture() {
    if (m_data) {
      _mm_free(m_data);
      m_data = NULL;
    }
  }

  template <class T>
  const size_t& Aperture<T>::NThreadsGet() const {
    return Aperture<T>::nthreads;
  }

  template <class T>
  void Aperture<T>::NThreadsSet(const size_t &nThreads) {
    Aperture<T>::nthreads = nThreads;
  }

  template <class T>
  void Aperture<T>::PositionsGet(T **pos, size_t* nElements, size_t* nParams) const {

    const size_t _nElements          = m_data->m_nelements;
    const size_t arrSize             = _nElements*3*sizeof(T); 
    const sps::point_t<T>* positions = m_data->m_pos;
    
    // The function allocates
    T* arr = (T*) malloc(arrSize);
    memset(arr,0,arrSize);

    for (size_t i = 0 ; i < _nElements ; i++) {
      arr[i*3 + 0] = positions[i][0];
      arr[i*3 + 1] = positions[i][1];
      arr[i*3 + 2] = positions[i][2];
    }
    // Assign outputs
    *nElements = m_data->m_nelements;
    *nParams   = 3;
    *pos = arr;
  }

  // TODO: Consider storing focus and compute this as part of CalcCW
  
  // Wrong, se FOCUS (Find phase for each element and adjust)
  template <class T>
  void Aperture<T>::FocusSet(const T iFocus[3]) {
    memcpy(&m_data->m_focus[0],&iFocus[0],3*sizeof(T));

    size_t nElements = m_data->m_nelements;

#if 0
    std::unique_ptr<T[]> delays = std::unique_ptr<T[]>(new T[nElements]);
    T maxDelay = T(0.0);
    for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
      delays[iElement] = norm(m_data->m_pos[iElement] - m_data->m_focus) / Aperture<T>::_sysparm.c;
      maxDelay = std::max<T>(maxDelay,delays[iElement]);
    }
    for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
      m_data->m_phases[iElement] = T(M_2PI) * (delays[iElement] - maxDelay) * m_data->m_f0;
    }
#else
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
#endif
  }

  template <class T>
  void Aperture<T>::PhasesGet(T** phases, size_t* nPhases) {
    
    size_t nElements = m_data->m_nelements;
    *phases = (T*) malloc(nElements*sizeof(T));
    if (nElements > 0) {
      memcpy(*phases,m_data->m_phases,nElements*sizeof(T));
      *nPhases = nElements;
    }
    else {
      *nPhases = 0;
    }
  }

  template <class T>
  void Aperture<T>::FocusGet(T oFocus[3]) const {
    memcpy(oFocus,&m_data->m_focus[0],3*sizeof(T));
  }

  template <class T>
  const T& Aperture<T>::F0Get() const {
    return m_data->m_f0;
  }

  template <class T>
  void Aperture<T>::F0Set(const T f0) {
    m_data->m_f0 = f0;
  }

  template <class T>
  const T& Aperture<T>::CGet()  const {
    return Aperture<T>::_sysparm.c;
  }
  
  template <class T>
  void Aperture<T>::CSet (const T c)  {
    Aperture<T>::_sysparm.c = c;
  }
  
  template <class T>
  const size_t& Aperture<T>::NElementsGet() const {
    return m_data->m_nelements;
  }

  template <class T>
  const size_t& Aperture<T>::NSubElementsGet() const {
    return m_data->m_nsubelements;
  }

  template <class T>
  const size_t& Aperture<T>::NDivWGet() const {
    return Aperture<T>::_sysparm.nDivW;
  }

  template <class T>
  void Aperture<T>::NDivWSet(const size_t nDivW) {
    Aperture<T>::_sysparm.nDivW = nDivW;
  }

  template <class T>
  const size_t& Aperture<T>::NDivHGet() const {
    return Aperture<T>::_sysparm.nDivH;
  }

  template <class T>
  void Aperture<T>::NDivHSet(const size_t nDivH) {
    Aperture<T>::_sysparm.nDivH = nDivH;
  }    

  template <class T>
  void Aperture<T>::CalcCwField(const T* pos, const size_t nPositions, const size_t nDim,
                                std::complex<T>** odata, size_t* nOutPositions) {
    
    assert(nDim == 3);
    if (nDim != 3) {
      *odata = NULL;
      *nOutPositions = 0;
      return;
    }
    size_t nScatterers = nPositions;
    *odata = (std::complex<T>*) malloc(nScatterers*sizeof(std::complex<T>));
    *nOutPositions = nScatterers;

    const T lambda = Aperture<T>::_sysparm.c / m_data->m_f0;
  
    const T k = T(M_2PI)/lambda;

    size_t nDivW = Aperture<T>::_sysparm.nDivW;
    size_t nDivH = Aperture<T>::_sysparm.nDivH;
  
    T* uweights = (T*) _mm_malloc(sizeof(T)*nDivW,16);
    T* vweights = (T*) _mm_malloc(sizeof(T)*nDivH,16);
  
    T* uxs = (T*) _mm_malloc(sizeof(T)*nDivW,16);
    T* vxs = (T*) _mm_malloc(sizeof(T)*nDivH,16);
  
    memset(uxs,0,sizeof(T)*nDivW);
    memset(vxs,0,sizeof(T)*nDivH);
  
    memset(uweights,0,sizeof(T)*nDivW);
    memset(vweights,0,sizeof(T)*nDivH);
    
    for (size_t i = 0 ; i < nDivW ; i++) {
      fastgl::QuadPair qp = fastgl::GLPair(nDivW, i+1);
      uxs[i]      = T(qp.x());
      uweights[i] = T(qp.weight);
    }
  
    for (size_t i = 0 ; i < nDivH ; i++) {
      fastgl::QuadPair qp = fastgl::GLPair(nDivH, i+1);
      vxs[i]      = T(qp.x());
      vweights[i] = T(qp.weight);
    }

    sps::point_t<T> point;
  
    size_t nElements = m_data->m_nelements;

    std::vector<std::vector<element_t<T> > >& elements = *m_data->m_elements;
    
    for (size_t iPoint = 0 ; iPoint < nPositions ; iPoint++) {
  
      point[0] = pos[iPoint*3 + 0];
      point[1] = pos[iPoint*3 + 1];
      point[2] = pos[iPoint*3 + 2];

      std::complex<T> field = std::complex<T>(T(0.0),T(0.0));
  
      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {

        std::complex<T> field1 = std::complex<T>(T(0.0),T(0.0));

        for (size_t jElement = 0 ; jElement < elements[iElement].size() ; jElement++) {

          element_t<T> element = elements[iElement][jElement];
    
          // Get basis vectors (can be stored)
          sps::point_t<T> hh_dir, hw_dir, normal;
    
          basis_vectors(hw_dir, element.euler, 0);
          basis_vectors(hh_dir, element.euler, 1);
          basis_vectors(normal, element.euler, 2);
    
          sps::point_t<T> r2p = point - element.center;
  
          // Distance to plane
          T dist2plane = fabs(dot(normal,r2p));
    
          // Projection onto plane
          T u = dot(hw_dir,r2p);
          T v = dot(hh_dir,r2p);
    
          T z  = dist2plane;
  
          T l = fabs(u) + element.hw;
          T s = fabs(v) + element.hh;
    
          field1 += CalcHz<T>(fabs(s),fabs(l),z,k,uxs,uweights,nDivW,
                              vxs,vweights,nDivH)*T(signum<T>(s)*signum<T>(l));
  
          l = element.hw - fabs(u);
  
          field1 += CalcHz<T>(fabs(s),fabs(l),z,k,uxs,uweights,nDivW,
                              vxs,vweights,nDivH)*T(signum<T>(s)*signum<T>(l));
  
          l = fabs(u) + element.hw;
          s = element.hh - fabs(v);
  
          field1 += CalcHz<T>(fabs(s),fabs(l),z,k,uxs,uweights,nDivW,
                              vxs,vweights,nDivH)*T(signum<T>(s)*signum<T>(l));
        
          l = element.hw - fabs(u);
  
          field1 += CalcHz<T>(fabs(s),fabs(l),z,k,uxs,uweights,nDivW,
                              vxs,vweights,nDivH)*T(signum<T>(s)*signum<T>(l));
        }

        T real = field1.real();
        T imag = field1.imag();

        field.real(field.real() + real*cos(m_data->m_phases[iElement]) - imag*sin(m_data->m_phases[iElement]));
        field.imag(field.imag() + real*sin(m_data->m_phases[iElement]) + imag*cos(m_data->m_phases[iElement]));
      }
      (*odata)[iPoint].real(field.real());
      (*odata)[iPoint].imag(field.imag());
    }
    
    _mm_free(uweights);
    _mm_free(vweights);
    _mm_free(uxs);
    _mm_free(vxs);
  }

  template <class T>
  void Aperture<T>::initElements() {

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

      // TODO: Consider keeping position on every sub-element

      // Center element (reference for delays)
      if (nSubElements % 2 == 1) {
        // Odd number of elements
        m_data->m_pos[iElement][0] = elements[iElement][nSubElements/2].center[0];
        m_data->m_pos[iElement][1] = elements[iElement][nSubElements/2].center[1];
        m_data->m_pos[iElement][2] = elements[iElement][nSubElements/2].center[2];
      }
      else {
        // Even number of elements (TODO: On an edge)
        m_data->m_pos[iElement][0] =
          T(0.5) * (elements[iElement][nSubElements-1].center[0] + elements[iElement][0].center[0]);
        m_data->m_pos[iElement][1] =
          T(0.5) * (elements[iElement][nSubElements-1].center[1] + elements[iElement][0].center[1]);
        m_data->m_pos[iElement][2] =
          T(0.5) * (elements[iElement][nSubElements-1].center[2] + elements[iElement][0].center[2]);
      }

      for (jElement=0 ; jElement < nSubElements ; jElement++) {
        element_t<T>& element = elements[iElement][jElement];

        // Hack (not working for double anyway)
        SPS_UNREFERENCED_PARAMETER(normal);
        sps::euler_t<float> euler;
        memcpy((void*)&euler,(void*)&element.euler,sizeof(sps::euler_t<float>));
        __m128 uvector, vvector, normal;
        basis_vectors_ps((float*)&uvector, (float*)&vvector, (float*)&normal, euler);
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
    memset(arr,0,arrSize);

    std::vector<std::vector<element_t<T> > >& elements = *m_data->m_elements;

    for (size_t iElement = 0 ; iElement < _nElements ; iElement++) {
      arr[iElement * nSubElementsPerElement * nElePosParams] = elements[iElement][0].hw;

      arr[iElement * nSubElementsPerElement * nElePosParams + 1] = elements[iElement][0].hh;

      memcpy(&arr[iElement * nSubElementsPerElement * nElePosParams + 2],
             &elements[iElement][0].center[0],
             sizeof(sps::point_t<T>));
      arr[iElement * nSubElementsPerElement * nElePosParams] = elements[iElement][0].euler.alpha;
      arr[iElement * nSubElementsPerElement * nElePosParams] = elements[iElement][0].euler.beta;
      arr[iElement * nSubElementsPerElement * nElePosParams] = elements[iElement][0].euler.gamma;
    }

    *nElements    = m_data->m_nelements;
    *nParams      = Aperture<T>::nElementPosParameters; // No need to introduce static variable

    *out          = arr;

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
  
  template <class T>
  bool Aperture<T>::ElementsSet(const T* pos, const size_t nElements, const size_t nDim) throw (std::runtime_error) {
    bool retval = this->SubElementsSet(pos,nElements,1,nDim);
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
  
  template <class T>
  void Aperture<T>::CalcCwFieldRef(const T* pos, const size_t nPositions, const size_t nDim,
                                   std::complex<T>** odata, size_t* nOutPositions) {

    assert(nDim == 3);
    if (nDim != 3) {
      *odata = NULL;
      *nOutPositions = 0;
      return;
    }
    size_t nScatterers = nPositions;
    *odata = (std::complex<T>*) malloc(nScatterers*sizeof(std::complex<T>));
    *nOutPositions = nScatterers;

    const T lambda = Aperture<T>::_sysparm.c / m_data->m_f0;
  
    const T k = T(M_2PI)/lambda;

    size_t nDivW = Aperture<T>::_sysparm.nDivW;
    size_t nDivH = Aperture<T>::_sysparm.nDivH;
  
    T* uweights = (T*) _mm_malloc(sizeof(T)*nDivW,16);
    T* vweights = (T*) _mm_malloc(sizeof(T)*nDivH,16);
  
    T* uxs = (T*) _mm_malloc(sizeof(T)*nDivW,16);
    T* vxs = (T*) _mm_malloc(sizeof(T)*nDivH,16);
  
    memset(uxs,0,sizeof(T)*nDivW);
    memset(vxs,0,sizeof(T)*nDivH);
  
    memset(uweights,0,sizeof(T)*nDivW);
    memset(vweights,0,sizeof(T)*nDivH);
    
    for (size_t i = 0 ; i < nDivW ; i++) {
      fastgl::QuadPair qp = fastgl::GLPair(nDivW, i+1);
      uxs[i]      = T(qp.x());
      uweights[i] = T(qp.weight);
    }
  
    for (size_t i = 0 ; i < nDivH ; i++) {
      fastgl::QuadPair qp = fastgl::GLPair(nDivH, i+1);
      vxs[i]      = T(qp.x());
      vweights[i] = T(qp.weight);
    }

    sps::point_t<T> point;
  
    size_t nElements = m_data->m_nelements;

    std::vector<std::vector<element_t<T> > >& elements = *m_data->m_elements;
    
    for (size_t iPoint = 0 ; iPoint < nPositions ; iPoint++) {
  
      point[0] = pos[iPoint*3 + 0];
      point[1] = pos[iPoint*3 + 1];
      point[2] = pos[iPoint*3 + 2];

      std::complex<T> field = std::complex<T>(T(0.0),T(0.0));
  
      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {

        std::complex<T> field1 = std::complex<T>(T(0.0),T(0.0));

        for (size_t jElement = 0 ; jElement < elements[iElement].size() ; jElement++) {

          element_t<T> element = elements[iElement][jElement];
    
          // Get basis vectors (can be stored)
          sps::point_t<T> hh_dir, hw_dir, normal;
    
          basis_vectors(hw_dir, element.euler, 0);
          basis_vectors(hh_dir, element.euler, 1);
          basis_vectors(normal, element.euler, 2);
    
          sps::point_t<T> r2p = point - element.center;
  
          // Distance to plane
          T dist2plane = fabs(dot(normal,r2p));
    
          // Projection onto plane
          T u = dot(hw_dir,r2p);
          T v = dot(hh_dir,r2p);
    
          T z  = dist2plane;
  
          T s = fabs(v) + element.hh;
    
          if (fabs(u) > element.hw) {
            field1 += CalcSingle<T>(fabs(u),            fabs(u)+element.hw, s, z, k, uxs,uweights,nDivW);
            field1 += CalcSingle<T>(fabs(u)-element.hw, fabs(u)           , s, z, k, uxs,uweights,nDivW);
            s = fabs(fabs(v) - element.hh);
            field1 += CalcSingle<T>(fabs(u),            fabs(u)+element.hw, s, z, k, uxs,uweights,nDivW);
            field1 += CalcSingle<T>(fabs(u)-element.hw, fabs(u)           , s, z, k, uxs,uweights,nDivW);
          }
          else {
            field1 += CalcSingle<T>(0,                  fabs(u)+element.hw, s, z, k, uxs,uweights,nDivW);
            field1 += CalcSingle<T>(0,                  element.hw-fabs(u), s, z, k, uxs,uweights,nDivW);
            s = element.hh - fabs(v);
            field1 += CalcSingle<T>(0,                  fabs(u)+element.hw, s, z, k, uxs,uweights,nDivW);
            field1 += CalcSingle<T>(0,                  element.hw-fabs(u), s, z, k, uxs,uweights,nDivW);
          }

          s = fabs(u) + element.hw;
          if (fabs(v) > element.hh) {
            field1 += CalcSingle<T>(fabs(v)-element.hh, fabs(v)           , s, z, k, vxs,uweights,nDivH);
            field1 += CalcSingle<T>(fabs(v),            fabs(v)+element.hh, s, z, k, vxs,uweights,nDivH);
            s = element.hw - fabs(u);
            field1 += CalcSingle<T>(fabs(v)-element.hh, fabs(v)           , s, z, k, vxs,uweights,nDivH);
            field1 += CalcSingle<T>(fabs(v),            fabs(v)+element.hh, s, z, k, vxs,uweights,nDivH);
          }
          else {
            field1 += CalcSingle<T>(0,                  fabs(v)+element.hh, s, z, k, vxs,vweights,nDivH);
            s = element.hw - fabs(u);
            field1 += CalcSingle<T>(0,                  fabs(v)+element.hh, s, z, k, vxs,vweights,nDivH);
            s = fabs(u) + element.hw;
            field1 += CalcSingle<T>(0,                  element.hh-fabs(v), s, z, k, vxs,vweights,nDivH);
            s = element.hw - fabs(u);
            field1 += CalcSingle<T>(0,                  element.hh-fabs(v), s, z, k, vxs,vweights,nDivH);
          }
        }
        
        T real = field1.real();
        T imag = field1.imag();

        field.real(field.real() + real*cos(m_data->m_phases[iElement]) - imag*sin(m_data->m_phases[iElement]));
        field.imag(field.imag() + real*sin(m_data->m_phases[iElement]) + imag*cos(m_data->m_phases[iElement]));
      }
      (*odata)[iPoint].real(field.real());
      (*odata)[iPoint].imag(field.imag());
    }

    _mm_free(uweights);
    _mm_free(vweights);
    _mm_free(uxs);
    _mm_free(vxs);
  }
  
  template <class T>
  void Aperture<T>::CalcCwField2(const T* pos, const size_t nPositions, const size_t nDim,
                                 std::complex<T>** odata, size_t* nOutPositions) {
  
    assert(nDim == 3);
    if (nDim != 3) {
      *odata = NULL;
      *nOutPositions = 0;
      return;
    }
    size_t nScatterers = nPositions;
    *odata = (std::complex<T>*) malloc(nScatterers*sizeof(std::complex<T>));
    memset(*odata, 0, nScatterers*sizeof(std::complex<T>));
    *nOutPositions = nScatterers;

    const T lambda = Aperture<T>::_sysparm.c / m_data->m_f0;
  
    const T k = T(M_2PI)/lambda;

    const size_t nElements = m_data->m_nelements;

    T* apodizations = m_data->m_apodizations;
    
    T apodization;

    sps::point_t<T> hh_dir, hw_dir, normal;

    // Need weights and abcissa values
  
    size_t nDivW = Aperture<T>::_sysparm.nDivW;
    size_t nDivH = Aperture<T>::_sysparm.nDivH;
    
    T* uweights = (T*) _mm_malloc(sizeof(T)*nDivW,16);
    T* vweights = (T*) _mm_malloc(sizeof(T)*nDivH,16);
    
    T* uxs = (T*) _mm_malloc(sizeof(T)*nDivW,16);
    T* vxs = (T*) _mm_malloc(sizeof(T)*nDivH,16);
    
    memset(uxs,0,sizeof(T)*nDivW);
    memset(vxs,0,sizeof(T)*nDivH);
    
    memset(uweights,0,sizeof(T)*nDivW);
    memset(vweights,0,sizeof(T)*nDivH);
    
    for (size_t i = 0 ; i < nDivW ; i++) {
      fastgl::QuadPair qp = fastgl::GLPair(nDivW, i+1);
      uxs[i]      = T(qp.x());
      uweights[i] = T(qp.weight);
    }
    
    for (size_t i = 0 ; i < nDivH ; i++) {
      fastgl::QuadPair qp = fastgl::GLPair(nDivH, i+1);
      vxs[i]      = T(qp.x());
      vweights[i] = T(qp.weight);
    }

    std::vector<std::vector<element_t<T> > >& elements = *m_data->m_elements;

    sps::point_t<T> projection;

    // Hack
    __m128* vec_normals  = (__m128*) _mm_malloc(nElements*sizeof(__m128),16);
    __m128* vec_uvectors = (__m128*) _mm_malloc(nElements*sizeof(__m128),16);
    __m128* vec_vvectors = (__m128*) _mm_malloc(nElements*sizeof(__m128),16);

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
                                  uxs, uweights, nDivW,
                                  vxs, vweights, nDivH);
            final = final + result * exp(std::complex<T>(0,m_data->m_phases[iElement]));
          }
        }
      }
      (*odata)[iPoint].real(final.real());
      (*odata)[iPoint].imag(final.imag());
    }

    _mm_free(vec_uvectors);
    _mm_free(vec_vvectors);
    _mm_free(vec_normals);
    _mm_free(uweights);
    _mm_free(vweights);
    _mm_free(uxs);
    _mm_free(vxs);
  }

  template <class T>
  void Aperture<T>::CalcCwFast(const T* pos,
                               const size_t nPositions,
                               const size_t nDim,
                               std::complex<T>** odata,
                               size_t* nOutPositions) {

    assert(nDim == 3);
    if (nDim != 3) {
      *odata = NULL;
      *nOutPositions = 0;
      return;
    }
    *odata = (std::complex<T>*) malloc(nPositions*sizeof(std::complex<T>));
    memset(*odata, 0, nPositions*sizeof(std::complex<T>));
    *nOutPositions = nPositions;

    const T lambda = Aperture<T>::_sysparm.c / m_data->m_f0;
    const T k = T(M_2PI)/lambda;

    size_t nDivW = Aperture<T>::_sysparm.nDivW;
    size_t nDivH = Aperture<T>::_sysparm.nDivH;
    
    T* uweights = (T*) _mm_malloc(sizeof(T)*nDivW,16);
    T* vweights = (T*) _mm_malloc(sizeof(T)*nDivH,16);
    
    T* uxs = (T*) _mm_malloc(sizeof(T)*nDivW,16);
    T* vxs = (T*) _mm_malloc(sizeof(T)*nDivH,16);
    
    memset(uxs,0,sizeof(T)*nDivW);
    memset(vxs,0,sizeof(T)*nDivH);
    
    memset(uweights,0,sizeof(T)*nDivW);
    memset(vweights,0,sizeof(T)*nDivH);

#ifdef GL_REF
    // Hack works for even only
    double* tmp_x = (double*) _mm_malloc(sizeof(double)*nDivW,16);
    double* tmp_w = (double*) _mm_malloc(sizeof(double)*nDivW,16);
    double* pPositiveX = tmp_x + nDivW/2;
    double* pPositiveW = tmp_w + nDivW/2;

    gauss_legendre_tbl(nDivW, pPositiveX, pPositiveW, 1e-10);
    for (size_t i = 0 ; i < nDivW/2 ; i++) {
      *(pPositiveX - 1 - i) = -pPositiveX[i];
      *(pPositiveW - 1 - i) = pPositiveW[i];
    }
    for (size_t i = 0 ; i < nDivW ; i++) {
      uxs[i]      = T(tmp_x[i]);
      uweights[i] = T(tmp_w[i]);
    }
    _mm_free(tmp_x);
    _mm_free(tmp_w);

    tmp_x = (double*) _mm_malloc(sizeof(double)*nDivH,16);
    tmp_w = (double*) _mm_malloc(sizeof(double)*nDivH,16);
    pPositiveX = tmp_x + nDivH/2;
    pPositiveW = tmp_w + nDivH/2;
    gauss_legendre_tbl(nDivH, pPositiveX, pPositiveW, 1e-10);
    for (size_t i = 0 ; i < nDivH/2 ; i++) {
      *(pPositiveX - 1 - i) = -pPositiveX[i];
      *(pPositiveW - 1 - i) = pPositiveW[i];
    }
    for (size_t i = 0 ; i < nDivH ; i++) {
      vxs[i]      = T(tmp_x[i]);
      vweights[i] = T(tmp_w[i]);
    }
    _mm_free(tmp_x);
    _mm_free(tmp_w);
#else
    // Common weights and abcissa
    for (size_t i = 0 ; i < nDivW ; i++) {
      fastgl::QuadPair qp = fastgl::GLPair(nDivW, i+1);
      uxs[i]      = T(qp.x());
      uweights[i] = T(qp.weight);
    }
    
    for (size_t i = 0 ; i < nDivH ; i++) {
      fastgl::QuadPair qp = fastgl::GLPair(nDivH, i+1);
      vxs[i]      = T(qp.x());
      vweights[i] = T(qp.weight);
    }
#endif

    // Common basis vectors
    const size_t nElements = m_data->m_nelements;
    std::vector<std::vector<element_t<T> > >& elements = *m_data->m_elements;

    __m128* vec_normals  = (__m128*) _mm_malloc(nElements*sizeof(__m128),16);
    __m128* vec_uvectors = (__m128*) _mm_malloc(nElements*sizeof(__m128),16);
    __m128* vec_vvectors = (__m128*) _mm_malloc(nElements*sizeof(__m128),16);


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

#ifndef HAVE_THREAD
    thread_arg threadarg;
    threadarg.iPointBegin = 0;
    threadarg.iPointEnd   = nPositions;
    threadarg.pos         = pos;
    threadarg.field       = *odata;
    threadarg.k           = k;
    threadarg.uxs         = uxs;
    threadarg.vweights    = uweights;
    threadarg.nDivU       = nDivW;
    threadarg.uxs         = vxs;
    threadarg.uweights    = vweights;
    // Unit vectors
    threadarg.normals     = (T*) vec_normals;
    threadarg.uvectors    = (T*) vec_uvectors;
    threadarg.vvectors    = (T*) vec_vvectors;
    threadarg.nDivV       = nDivH;
    threadarg.thread_id   = 0;
    threadarg.cpu_id      = 0;
# ifdef _WIN32
    unsigned int retval   = CalcThreaded((void*)&threadarg);
# else
    void* thread_retval   = CalcThreaded((void*)&threadarg);
# endif
    return;
#else

    int nproc;

# if defined(HAVE_PTHREAD_H)
# elif defined(_WIN32)
    unsigned int threadID;
    HANDLE threads[N_MAX_THREADS];
# endif

    nproc = getncpus();

# ifdef HAVE_MQUEUE_H
    sps::pthread_launcher<Aperture<T>,
                     &Aperture<T>::CalcThreaded> launcher[N_MAX_THREADS];
# else
    sps::pthread_launcher<Aperture<T>,
        &Aperture<T>::CalcThreaded> launcher[N_MAX_THREADS];
# endif

    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

# ifdef HAVE_PTHREAD_H
    CallErr(pthread_attr_init,(&Aperture<T>::attr));
#  ifdef HAVE_MQUEUE_H
    CallErr(pthread_attr_setdetachstate, (&Aperture<T>::attr, PTHREAD_CREATE_DETACHED));
#  else
    CallErr(pthread_attr_setdetachstate, (&Aperture<T>::attr, PTHREAD_CREATE_JOINABLE));
#  endif
# endif

    // Populate structs for threads
    for (size_t i=0 ; i < nthreads ; i++) {
      threadarg[i].iPointBegin = 0+i*(nPositions/nthreads);
      threadarg[i].iPointEnd   = (nPositions/nthreads)+i*(nPositions/nthreads);
      threadarg[i].pos         = pos;
      threadarg[i].field       = (*odata);
      threadarg[i].k           = k;
      threadarg[i].uxs         = uxs;
      threadarg[i].uweights    = uweights;
      threadarg[i].nDivU       = nDivW;
      threadarg[i].vxs         = vxs;
      threadarg[i].vweights    = vweights;
      threadarg[i].nDivV       = nDivH;
      threadarg[i].normals     = (T*) vec_normals;
      threadarg[i].uvectors    = (T*) vec_uvectors;
      threadarg[i].vvectors    = (T*) vec_vvectors;
      threadarg[i].thread_id   = i;
      threadarg[i].cpu_id      = ((int) i) % nproc;
      if (i==(nthreads-1))
        threadarg[i].iPointEnd   = nPositions;      
    }
    
# ifdef HAVE_MQUEUE_H
    // Mesage queues
    struct mq_attr qattr;
  
    char buf[MSGMAX];
    unsigned int prio;

    // Create two message queues
    CallErrExit(Aperture<T>::mqd_master = mq_open,
                ("/jmh-master", O_RDWR | O_CREAT,
                 S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
                 NULL), EXIT_FAILURE);
    CallErrExit(mq_close,(Aperture<T>::mqd_master),EXIT_FAILURE);

    CallErrExit(Aperture<T>::mqd_client = mq_open,
                ("/jmh-master", O_RDWR | O_CREAT,
                 S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
                 NULL), EXIT_FAILURE);
    CallErrExit(mq_close,(Aperture<T>::mqd_client),EXIT_FAILURE);

    /* Initiate launcher, note that the threadarg are global and updated
       in between calls to CalcThreaded */
    for (size_t i=0;i<nthreads;i++) {
        launcher[i]               =
            sps::pthread_launcher<Aperture<T>,&Aperture<T>::CalcThreaded>(this,
                                                        &threadarg[i]);
    }

    if ((!Aperture<T>::threads_initialized) && !(nthreads>N_MAX_THREADS)) {

        /* Clear message queues (and create if they don't exist), TODO: Do
           this in two stages */
        CallErrExit(mq_clear,("/jmh-master"),false);
        CallErrExit(mq_clear,("/jmh-client"),false);

        /* Initiate threads */
        for (size_t i=0;i<nthreads;i++) {
            CallErr(pthread_create,
                    (&threads[i], &Aperture<T>::attr,
                     sps::launch_member_function<sps::pthread_launcher<Aperture<T>,
                     &Aperture<T>::CalcThreaded> >, 
                     &launcher[i]));
        }
        Aperture<T>::threads_initialized = true;
    }

    /* Signal idling threads to run */
    CallErrExit(Aperture<T>::mqd_client = mq_open,
                ("/jmh-client", O_WRONLY,
                 S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
                 NULL), false);
    for (size_t i=0;i<nthreads;i++) {
        prio = i;
        if (mq_send(mqd_client,"RUN",3,prio)==-1)
            perror ("mq_send()");
    }
    mq_close(Aperture<T>::mqd_client);
        
    /* Wait for threads to finish */
    CallErrExit(Aperture<T>::mqd_master = mq_open,
                ("/jmh-master", O_RDONLY,
                 S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
                 NULL), false);
        
    mq_getattr(Aperture<T>::mqd_master, &qattr);
        
    size_t n_to_go = nthreads;
    size_t n_idle  = 0;
    /* Wait for threads to finish */
    while (true) {
      mq_receive(Aperture<T>::mqd_master, &buf[0], qattr.mq_msgsize, &prio);
      if (!strncmp(buf,"DONE",4)) {
        n_to_go--;
      }
      else if (!strncmp(buf,"READY",5)) {
        n_idle++;
      }
      else {
        printf("Unexpected message: %s\n",buf);
        retval = false;
        break;
      }
      if (n_to_go==0)
        break;
    }
    mq_close(Aperture<T>::mqd_master);
# else
    /* Without message queues (slower) */
    for (size_t i=0;i<nthreads;i++) {
      launcher[i]               =
        sps::pthread_launcher<Aperture<T>,&Aperture<T>::CalcThreaded>(this,
                                                                      &threadarg[i]);
    }
    for (size_t i=0;i<nthreads;i++) {
#  if defined(HAVE_PTHREAD_H)
      CallErr(pthread_create,
              (&threads[i], &Aperture<T>::attr,
               sps::launch_member_function<sps::pthread_launcher<Aperture<T>,
               &Aperture<T>::CalcThreaded> >, 
               &launcher[i]));
#  elif defined(_WIN32)
      threads[i] =
        (HANDLE)_beginthreadex(NULL, 0,
                               sps::launch_member_function<sps::pthread_launcher<Aperture<T>,
                               &Aperture<T>::CalcThreaded> >,
                               &launcher[i], 0, &threadID );
#  endif
    }
    
    for (size_t i = 0; i < nthreads; i++) {
#  if defined(HAVE_PTHREAD_H)
      CallErr(pthread_join,(threads[i],NULL));
#  elif defined(_WIN32)
      WaitForSingleObject(threads[i], INFINITE );
#  endif
    }
    
#  if defined(HAVE_PTHREAD_H)
    CallErr(pthread_attr_destroy,(&Aperture<T>::attr));
#  endif
# endif
#endif
    
    _mm_free(vec_uvectors);
    _mm_free(vec_vvectors);
    _mm_free(vec_normals);
    _mm_free(uweights);
    _mm_free(vweights);
    _mm_free(uxs);
    _mm_free(vxs);

  }

  template <class T>
  void Aperture<T>::CalcCwFocus(const T* pos, const size_t nPositions, const size_t nDim,
                                std::complex<T>** odata, size_t* nOutPositions) {
  
    assert(nDim == 3);
    if (nDim != 3) {
      *odata = NULL;
      *nOutPositions = 0;
      return;
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
    
    T* uweights = (T*) _mm_malloc(sizeof(T)*nDivW,16);
    T* vweights = (T*) _mm_malloc(sizeof(T)*nDivH,16);
    
    T* uxs = (T*) _mm_malloc(sizeof(T)*nDivW,16);
    T* vxs = (T*) _mm_malloc(sizeof(T)*nDivH,16);
    
    memset(uxs,0,sizeof(T)*nDivW);
    memset(vxs,0,sizeof(T)*nDivH);
    
    memset(uweights,0,sizeof(T)*nDivW);
    memset(vweights,0,sizeof(T)*nDivH);
    
    for (size_t i = 0 ; i < nDivW ; i++) {
      fastgl::QuadPair qp = fastgl::GLPair(nDivW, i+1);
      uxs[i]      = T(qp.x());
      uweights[i] = T(qp.weight);
    }
    
    for (size_t i = 0 ; i < nDivH ; i++) {
      fastgl::QuadPair qp = fastgl::GLPair(nDivH, i+1);
      vxs[i]      = T(qp.x());
      vweights[i] = T(qp.weight);
    }
    sps::point_t<T> hh_dir, hw_dir, normal;

    std::vector<std::vector<element_t<T> > >& elements = *m_data->m_elements;

    sps::point_t<T> projection;

    // Hack
    __m128* vec_normals  = (__m128*) _mm_malloc(nElements*sizeof(__m128),16);
    __m128* vec_uvectors = (__m128*) _mm_malloc(nElements*sizeof(__m128),16);
    __m128* vec_vvectors = (__m128*) _mm_malloc(nElements*sizeof(__m128),16);

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
      (*odata)[iPoint].real(final.real());
      (*odata)[iPoint].imag(final.imag());
    }

    _mm_free(vec_uvectors);
    _mm_free(vec_vvectors);
    _mm_free(vec_normals);
    _mm_free(uweights);
    _mm_free(vweights);
    _mm_free(uxs);
    _mm_free(vxs);
  }

#if defined(HAVE_PTHREAD_H)
  template <class T>
  void* Aperture<T>::CalcThreaded(void* ptarg)
# ifdef HAVE_MQUEUE_H
  template <class T>
  void Aperture<T>::CalcThreaded(void* ptarg)
# endif
#else
  template <class T>
  unsigned int __stdcall Aperture<T>::CalcThreaded(void *ptarg)
#endif
  {
    thread_arg* pThreadArg = reinterpret_cast<thread_arg*>(ptarg);

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

            // CalcHzFast4 works inside only
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
#if HAVE_PTHREAD_H
    pthread_exit(NULL);
#else
    return 0;
#endif
  }

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
                       const size_t nVs) {

  const T carg = cos(-k*z);
  const T sarg = sin(-k*z);
  const T l_2 = T(0.5) * l;
  const T s_2 = T(0.5) * s;

  const T z2 = SQUARE(z);
  const T l2 = SQUARE(l);
  const T s2 = SQUARE(s);
  
  // integral width
  std::complex<T> intW = std::complex<T>(T(0.0),T(0.0));

  for (size_t iu = 0 ; iu < nUs ; iu++) {

    T ls = l_2 * uxs[iu] + l_2;
    T ls2 = SQUARE(ls);

    T argw = -k * sqrt(ls2 + z2 + s2);
    T real = uweights[iu] * (cos(argw) - carg) / (ls2+s2);
    T imag = uweights[iu] * (sin(argw) - sarg) / (ls2+s2);
    intW.real(real + intW.real());
    intW.imag(imag + intW.imag());
  }
  intW *= l_2 * s / (T(M_2PI)*k);
  T realW = intW.real();
  intW.real(intW.imag());
  intW.imag(-realW);
  
  // integral height
  std::complex<T> intH = std::complex<T>(T(0.0),T(0.0));
  
  for(size_t iv = 0 ; iv < nVs ; iv++) {
    T ss = s_2 * vxs[iv] + s_2;
    T ss2 = SQUARE(ss);
    T argh = -k * sqrt(ss2 + z2 + l2);
    T real = vweights[iv] * (cos(argh) - carg) / (ss2+l2);
    T imag = vweights[iv] * (sin(argh) - sarg) / (ss2+l2);
    intH.real(real + intH.real());
    intH.imag(imag + intH.imag());
  }
  intH *= s_2 * l / (T(M_2PI)*k);
  T realH = intH.real();
  intH.real(intH.imag());
  intH.imag(-realH);

  intH = intH + intW;
  return intH;
}

template <class T>
STATIC_INLINE_BEGIN std::complex<T> CalcSingle(const T& s1,
                                               const T& s2,
                                               const T& l,
                                               const T& z,
                                               const T& k,
                                               const T* uxs,
                                               const T* uweights,
                                               const size_t nUs) {

  const T carg = cos(-k*z);
  const T sarg = sin(-k*z);
  const T sm = T(0.5) * (s2 - s1);
  const T sp = T(0.5) * (s2 + s1);

  const T z2 = SQUARE(z);
  const T l2 = SQUARE(l);
  
  // integral width
  std::complex<T> intW = std::complex<T>(T(0.0),T(0.0));

  for (size_t iu = 0 ; iu < nUs ; iu++) {

    T s = sm * uxs[iu] + sp;
    T s2 = SQUARE(s);

    T argw = -k * sqrt(s2 + z2 + l2);
    T real = uweights[iu] * (cos(argw) - carg) / (s2+l2);
    T imag = uweights[iu] * (sin(argw) - sarg) / (s2+l2);
    intW.real(real + intW.real());
    intW.imag(imag + intW.imag());
  }
  intW *= sm * l / (T(M_2PI)*k);
  T realW = intW.real();
  intW.real(intW.imag());
  intW.imag(-realW);
  
  return intW;
}

namespace fnm {

#ifdef __GNUG__
  template <class T>
  struct Aperture<T>::thread_arg Aperture<T>::threadarg[N_MAX_THREADS];
#elif defined(_WIN32)
  // Explicit instantiations
  template struct Aperture<float>::thread_arg;
  template struct Aperture<double>::thread_arg;

  Aperture<float>::thread_arg Aperture<float>::threadarg[N_MAX_THREADS];
  Aperture<double>::thread_arg Aperture<double>::threadarg[N_MAX_THREADS];
#endif

  template class Aperture<float>;
  template class Aperture<double>;

}

  // Template instantiation
#ifdef HAVE_MQUEUE_H
template class sps::pthread_launcher<fnm::Aperture<T>,&fnm::Aperture<T>::thread_func>;
#else
template class sps::pthread_launcher<fnm::Aperture<float>,&fnm::Aperture<float>::CalcThreaded>;
#endif

template std::complex<double> FNM_EXPORT CalcHzVecGL(const double& s,
                                                     const double& l,
                                                     const double& z,
                                                     const double& k,
                                                     const double* uxs,
                                                     const double* uweights,
                                                     const size_t nUs,
                                                     const double* vxs,
                                                     const double* vweights,
                                                     const size_t nVs);

template std::complex<float> FNM_EXPORT CalcHzVecGL(const float& s,
                                                    const float& l,
                                                    const float& z,
                                                    const float& k,
                                                    const float* uxs,
                                                    const float* uweights,
                                                    const size_t nUs,
                                                    const float* vxs,
                                                    const float* vweights,
                                                    const size_t nVs);

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
