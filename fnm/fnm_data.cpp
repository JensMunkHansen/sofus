#include <fnm/config.h>
#include <fnm/fnm_data.hpp>
#include <sps/smath.hpp>
#include <sps/debug.h>
#include <string.h>
#include <assert.h>

namespace fnm {

  template <class T>
  size_t ApertureData<T>::nextID = 0;

  template <class T>
  ApertureData<T>::ApertureData() : m_nelements(0), m_nsubelements(0), m_npos(0)
  {

    m_elements     = element_array(m_nelements,m_nsubelements);

    m_pos          = sps::deleted_aligned_array_create<sps::point_t<T> >(m_npos);
    m_apodizations = sps::deleted_aligned_array_create<T>(m_npos);
    m_phases       = sps::deleted_aligned_array_create<T>(m_nelements*m_nsubelements);
    m_delays       = sps::deleted_aligned_array_create<T>(m_npos);
    m_rectangles   = sps::deleted_aligned_array_create<sps::rect_t<T>>(m_nelements*m_nsubelements);
    m_boxes        = sps::deleted_aligned_array_create<sps::bbox_t<T>>(m_nelements);
    m_h_xyz[0]     = T(0.0);
    m_h_xyz[1]     = T(0.0);
    m_h_xyz[2]     = T(0.0);

    m_f0 = T(1e6);

    m_focus = sps::point_t<T>();
    memset(&this->m_focus[0],0,3*sizeof(T));

    m_center_focus = sps::point_t<T>();
    memset(&this->m_center_focus[0],0,3*sizeof(T));

    m_focus_type = FocusingType::Rayleigh;

    m_focus_valid = FocusingType::FocusingTypeCount;

    // Set unique identifier
    m_id = nextID;
    nextID++;

  }

  template <class T>
  ApertureData<T>::~ApertureData()
  {
  }

  template <class T>
  int ApertureData<T>::DelaysRefGet(size_t* nElements, const T*& delays) const
  {
    *nElements = m_nelements;
    delays = m_delays.get();
    return 0;
  }

  template <class T>
  int ApertureData<T>::ApodizationsRefGet(size_t* nElements, const T*& apodizations) const
  {
    *nElements   = m_nelements;
    apodizations = m_apodizations.get();
    return 0;
  }

  template <class T>
  int ApertureData<T>::ElementsRefGet(size_t* nElements, size_t* nSubElements,
                                      const sps::element_t<T>**& elements) const
  {
    *nElements    = m_nelements;
    *nSubElements = m_nsubelements;
    elements      = m_elements.get();
    return 0;
  }

  template <class T>
  void ApertureData<T>::ElementsSet(element_array &&elements,
                                    const size_t& nRows,
                                    const size_t& nCols)
  {
    if (nRows * nCols > 0) {
      assert(nCols == elements.n);
      this->m_nelements    = nRows;
      this->m_nsubelements = nCols;
      this->m_elements     = std::move(elements);

      initElements();
      initVectors();

      // TODO: Initialize this when needed, i.e. when using "Far field approx" or for display
      initRectangles();

      debug_print("Rectangles of first element\tx: %f %f %f %f\n", m_elements[0][0].vertices[0][0], m_elements[0][0].vertices[0][1], m_elements[0][0].vertices[0][2], m_elements[0][0].vertices[0][3]);
      debug_print("\t\t\t\t\ty: %f %f %f %f\n", m_elements[0][0].vertices[1][0], m_elements[0][0].vertices[1][1], m_elements[0][0].vertices[1][2], m_elements[0][0].vertices[1][3]);
      debug_print("\t\t\t\t\tz: %f %f %f %f\n", m_elements[0][0].vertices[2][0], m_elements[0][0].vertices[2][1], m_elements[0][0].vertices[2][2], m_elements[0][0].vertices[2][3]);
    }
  }

  template <class T>
  void ApertureData<T>::ElementExtentGet(const size_t iElement, sps::bbox_t<T>& bbox) const
  {
    assert(iElement < m_nelements);
    bbox = m_boxes[iElement];
  }

  template <class T>
  void ApertureData<T>::ExtentGet(sps::bbox_t<T>& bbox) const
  {
    // TODO: Consider caching this box
    sps::point_t<T> h_dir,w_dir;

    // Coordinate of lower left position reference
    T pos_ll;

    const size_t nElements       = m_nelements;

    const auto& elements         = m_elements;

    for(size_t iElement = 0 ; iElement < nElements ; iElement++) {
      for(size_t jElement = 0 ; jElement < m_nsubelements ; jElement++) {
        const sps::element_t<T>& element = elements[iElement][jElement];

        sps::basis_vectors<T, sps::EulerIntrinsicYXY>(w_dir,element.euler,0);
        sps::basis_vectors<T, sps::EulerIntrinsicYXY>(h_dir,element.euler,1);

        for(size_t i_xyz=0 ; i_xyz < 3 ; i_xyz++) {

          // Corner position of rectangle (lower left)
          pos_ll =
            element.center[i_xyz]
            - element.hh*h_dir[i_xyz]
            - element.hw*w_dir[i_xyz];

          // Compute aperture bounding box
          if((iElement == 0) && (jElement == 0)) {
            bbox.min[i_xyz] = pos_ll;
            bbox.max[i_xyz] = pos_ll;
          }

          bbox.min[i_xyz] = std::min<T>(bbox.min[i_xyz],pos_ll);
          bbox.max[i_xyz] = std::max<T>(bbox.max[i_xyz],pos_ll);

          T pos_ul = pos_ll+2*element.hh*h_dir[i_xyz];
          bbox.min[i_xyz] = std::min<T>(bbox.min[i_xyz],pos_ul);
          bbox.max[i_xyz] = std::max<T>(bbox.max[i_xyz],pos_ul);

          T pos_lr = pos_ll+2*element.hw*w_dir[i_xyz];
          bbox.min[i_xyz] = std::min<T>(bbox.min[i_xyz],pos_lr);
          bbox.max[i_xyz] = std::max<T>(bbox.max[i_xyz],pos_lr);

          T pos_ur = pos_ll+2*element.hh*h_dir[i_xyz]+2*element.hw*w_dir[i_xyz];
          bbox.min[i_xyz] = std::min<T>(bbox.min[i_xyz],pos_ur);
          bbox.max[i_xyz] = std::max<T>(bbox.max[i_xyz],pos_ur);
        }
      }
    }
  }

  template<class T>
  T ApertureData<T>::AreaGet() const
  {
    T area = T(0.0);
    const size_t nElements       = m_nelements;

    const auto& elements = m_elements;

    for(size_t iElement = 0 ; iElement < nElements ; iElement++) {
      for(size_t jElement = 0 ; jElement < m_nsubelements ; jElement++) {
        const sps::element_t<T>& element = elements[iElement][jElement];
        area = area + T(4.0) * element.hh * element.hw;
      }
    }
    return area;
  }

  template<class T>
  void ApertureData<T>::initVectors()
  {
    // Inputs: elements[][]
    // Output: elements[][]

    const size_t nElements              = m_nelements;
    const size_t nSubElementsPerElement = m_nsubelements;

    element_array& elements = m_elements;

    // Will never fail. Consider removing m_nelements and m_nsubelements
    debug_print("nElements: %zu, nSubElementsPerElement: %zu, m: %zu, n: %zu\n",
                nElements, nSubElementsPerElement, elements.m, elements.n);
    assert(elements.m == nElements);
    assert(elements.n == nSubElementsPerElement);

    for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
      for (size_t jElement = 0 ; jElement < nSubElementsPerElement ; jElement++) {
        auto& element = elements[iElement][jElement];
        // Try replace this
        sps::basis_vectors((T*)&element.uvector[0],
                           (T*)&element.vvector[0],
                           (T*)&element.normal[0],
                           element.euler);
      }
    }
  }

  template <class T>
  void ApertureData<T>::initElements()
  {
    const size_t nElements              = m_nelements;
    const size_t nSubElements           = m_nsubelements;

    auto& elements = m_elements;

    m_apodizations = sps::deleted_aligned_array_create<T>(nElements);
    m_pos          = sps::deleted_aligned_array_create<sps::point_t<T> >(nElements);
    m_phases       = sps::deleted_aligned_array_create<T>(nElements*nSubElements);
    m_delays       = sps::deleted_aligned_array_create<T>(nElements);

    for(size_t iElement=0 ; iElement < nElements ; iElement++) {
      // Set apodizations
      m_apodizations[iElement] = T(1.0);

      // Set delays
      m_delays[iElement] = T(0.0);

      // Center element (reference for delays)
      if(nSubElements % 2 == 1) {
        // Odd number of elements
        m_pos[iElement] = elements[iElement][nSubElements/2].center;
      } else {
        // Even number of elements
        m_pos[iElement] = T(0.5) * (elements[iElement][nSubElements-1].center + elements[iElement][0].center);
      }
      for (size_t jElement = 0 ; jElement < nSubElements ; jElement++) {
        // Set phases - do we need a phase per sub-element?
        m_phases[iElement*nSubElements + jElement] = T(0.0);
      }
    }
  }

  template <class T>
  void ApertureData<T>::initRectangles()
  {
    const size_t nElements              = m_nelements;
    const size_t nSubElements           = m_nsubelements;

    // Not used right now
    const size_t nSubH           = 1;
    const size_t nSubW           = 1;
    const size_t nSubSubElements = nSubH*nSubW;

    auto& elements = m_elements;

    sps::point_t<T> h_dir, w_dir, normal;

    // Delta height and width
    T dh,dw;

    // Lower left position reference
    T pos_ll;

    // Maximum width or height used for memory allocation
    m_h_xyz[0] = m_h_xyz[1] = m_h_xyz[2] = T(0.0);

    m_rectangles   = sps::deleted_aligned_array_create<sps::rect_t<T> >(nElements*nSubElements*nSubSubElements);

    m_boxes        = sps::deleted_aligned_array_create<sps::bbox_t<T> >(nElements);

    for(size_t iElement=0 ; iElement < nElements ; iElement++) {
      for(size_t jElement=0 ; jElement < nSubElements ; jElement++) {
        sps::element_t<T>& element = elements[iElement][jElement];

        sps::euler_t<T> euler;
        memcpy((void*)&euler,(void*)&element.euler,sizeof(sps::euler_t<T>));

        basis_vectors((T*)&w_dir[0],(T*)&h_dir[0],(T*)&normal[0],euler);
        dh = 2*element.hh/nSubH;
        dw = 2*element.hw/nSubW;

        for(size_t i_xyz=0 ; i_xyz < 3 ; i_xyz++) {

          T hh_vec = element.hh*h_dir[i_xyz];
          T hw_vec = element.hw*w_dir[i_xyz];

          // Maximum half width or height used for memory allocation
          m_h_xyz[i_xyz] = std::max<T>(m_h_xyz[i_xyz], fabs(hh_vec));
          m_h_xyz[i_xyz] = std::max<T>(m_h_xyz[i_xyz], fabs(hw_vec));

          // Corner position of rectangle (lower left)
          pos_ll =
            element.center[i_xyz]
            - hh_vec
            - hw_vec;

          T pos_ul = pos_ll + T(2.0)*hh_vec;
          T pos_lr = pos_ll + T(2.0)*hw_vec;
          T pos_ur = pos_ll + T(2.0)*hh_vec + T(2.0)*hw_vec;

          // Vertices of any sub-element are stored.
          element.vertices[i_xyz][0] = pos_ur;
          element.vertices[i_xyz][1] = pos_ul;
          element.vertices[i_xyz][2] = pos_ll;
          element.vertices[i_xyz][3] = pos_lr;

          // Compute coordinates for corners of mathematical elements
          for(size_t i_subh=0 ; i_subh < nSubH ; i_subh++) {
            for(size_t i_subw=0 ; i_subw < nSubW ; i_subw++) {
              m_rectangles[iElement*nSubSubElements*nSubElements +
                           jElement*nSubSubElements +
                           nSubH*i_subw + i_subh][0][i_xyz] =
                             pos_ll + i_subh*dh*h_dir[i_xyz]+i_subw*dw*w_dir[i_xyz];

              m_rectangles[iElement*nSubSubElements*nSubElements +
                           jElement*nSubSubElements +
                           nSubH*i_subw + i_subh][1][i_xyz] =
                             pos_ll + i_subh*dh*h_dir[i_xyz]+(i_subw+1)*dw*w_dir[i_xyz];

              m_rectangles[iElement*nSubSubElements*nSubElements +
                           jElement*nSubSubElements +
                           nSubH*i_subw + i_subh][2][i_xyz] =
                             pos_ll + (i_subh+1)*dh*h_dir[i_xyz]+i_subw*dw*w_dir[i_xyz];

              m_rectangles[iElement*nSubSubElements*nSubElements +
                           jElement*nSubSubElements +
                           nSubH*i_subw + i_subh][3][i_xyz] =
                             pos_ll + (i_subh+1)*dh*h_dir[i_xyz]+(i_subw+1)*dw*w_dir[i_xyz];

            } // for i_subw
          } // for i_subh

          // Test: Caching bounding box for sub-elements
          if (jElement == 0) {
            m_boxes[iElement].min[i_xyz] = pos_ur;
            m_boxes[iElement].max[i_xyz] = pos_ur;
          }

          m_boxes[iElement].min[i_xyz] = std::min<T>(m_boxes[iElement].min[i_xyz],pos_ur);
          m_boxes[iElement].min[i_xyz] = std::min<T>(m_boxes[iElement].min[i_xyz],pos_ul);
          m_boxes[iElement].min[i_xyz] = std::min<T>(m_boxes[iElement].min[i_xyz],pos_ll);
          m_boxes[iElement].min[i_xyz] = std::min<T>(m_boxes[iElement].min[i_xyz],pos_lr);
          m_boxes[iElement].max[i_xyz] = std::max<T>(m_boxes[iElement].max[i_xyz],pos_ur);
          m_boxes[iElement].max[i_xyz] = std::max<T>(m_boxes[iElement].max[i_xyz],pos_ul);
          m_boxes[iElement].max[i_xyz] = std::max<T>(m_boxes[iElement].max[i_xyz],pos_ll);
          m_boxes[iElement].max[i_xyz] = std::max<T>(m_boxes[iElement].max[i_xyz],pos_lr);
        } // for i_xyz
      } // jElement
    } // iElement
  }

  // If this was header-only, we could export it
  template class ApertureData<float>;

#ifdef FNM_DOUBLE_SUPPORT
  template class ApertureData<double>;
#endif
}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
