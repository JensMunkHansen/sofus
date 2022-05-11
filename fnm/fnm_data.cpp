/**
 * @file   fnm_data.cpp
 * @author Jens Munk Hansen <jmh@debian9laptop.parknet.dk>
 * @date   Fri Nov 10 00:18:14 2017
 *
 * @brief
 *
 * Copyright 2017 Jens Munk Hansen
 */
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

#include <fnm/config.h>
#include <sps/debug.h>
#include <string.h>
#include <assert.h>

#include <fnm/fnm_data.hpp>
#include <fnm/fnm_alias.hpp>
#include <sps/smath.hpp>

namespace fnm {

template <class T>
size_t ApertureData<T>::nextID = 0;

template <class T>
ApertureData<T>::ApertureData() : m_nelements(0), m_nsubelements(0), m_npos(0) {
  m_elements     = element_array(m_nelements, m_nsubelements);

  m_pos          = sps::unique_aligned_array_create<sps::point_t<T> >(m_npos);
  m_sensitivities = sps::unique_aligned_array_create<T>(m_npos);
  m_phases       = sps::unique_aligned_array_create<T>(m_nelements*m_nsubelements);
  m_delays       = sps::unique_aligned_array_create<T>(m_npos);

  m_f0 = T(1e6);

  m_focus = sps::point_t<T>();
  memset(&this->m_focus[0], 0, 3*sizeof(T));

  m_focus2 = sps::point_t<T>();
  memset(&this->m_focus2[0], 0, 3*sizeof(T));

  m_center_focus = sps::point_t<T>();
  memset(&this->m_center_focus[0], 0, 3*sizeof(T));

  m_focus_type = FocusingType::Rayleigh;

  m_apodization_type = ApodizationType::ApodizationTypeNonParametric;

  m_fnumber = T(1.0);

  // Internal state variables
  m_rectangles   = sps::unique_aligned_array_create<sps::rect_t<T>>(m_nelements*m_nsubelements);
  m_boxes        = sps::unique_aligned_array_create<sps::bbox_t<T>>(m_nelements);

  memset(m_h_xyz, 0, 3*sizeof(T));

  m_focus_valid = FocusingType::FocusingTypeCount;

  // Set unique identifier
  m_id = nextID;
  nextID++;

#ifdef FNM_CLOSURE_FUNCTIONS
  // Use this instead of thread_local variable
  c_idims[0] = 0;
  c_idims[1] = 0;
  c_idims[2] = 0;
  c_odims[0] = nullptr;
  c_odims[1] = nullptr;
  c_odims[2] = nullptr;
  c_iPtr = nullptr;
  c_ptr = nullptr;
#endif

#if 1
  m_pulses     = new sofus::AperturePulses<T>();
#endif
}

template <class T>
ApertureData<T>::~ApertureData() {
  // For time-domain pulsed-wave simulation, another data object is used
#if 1
  if (m_pulses) {
    delete m_pulses;
    m_pulses = nullptr;
  }
#endif
}

template <class T>
int ApertureData<T>::DelaysRefGet(size_t* nElements,
                                  const T*& delays) const {
  *nElements = m_nelements;
  delays = m_delays.get();
  return 0;
}

template <class T>
int ApertureData<T>::ApodizationsRefGet(size_t* nElements,
                                        const T*& apodizations) const {
  *nElements   = m_nelements;
  apodizations = m_sensitivities.get();
  return 0;
}

template <class T>
int ApertureData<T>::ElementsRefGet(size_t* nElements, size_t* nSubElements,
                                    const sps::element_rect_t<T>**& elements) const {
  *nElements    = m_nelements;
  *nSubElements = m_nsubelements;
  elements      = m_elements.get();
  return 0;
}

// TODO(JMH): Make more general (this is stupid)
template <class T>
int ApertureData<T>::
ApodizationSet(const sps::point_t<T>& direction,
               const T& depth, const ApodizationType& type) {
  const size_t nElements = this->m_nelements;
  const T eps = std::numeric_limits<T>::epsilon();

  int retval = -1;

  sps::point_t<T> center = this->m_center_focus;

  if (this->m_fnumber > eps) {
    retval = 0;
  }

  T apSize = depth / std::max<T>(this->m_fnumber, eps);

  for (size_t iElement = 0; iElement < nElements; iElement++) {
    T val = T(0.0);
    T dist2line = dist_point_to_line(m_pos[iElement], center, direction);

    T indexNorm = dist2line * T(2.0) / std::max<T>(apSize, eps);
    switch (type) {
    case ApodizationType::ApodizationTypeRectangular:
      if (indexNorm < T(1.0)) {
        val = T(1.0);
      }
      break;
    case ApodizationType::ApodizationTypeHamming:
      if (indexNorm < T(1.0)) {
        val = T(0.54) + T(0.46)*cos(T(M_PI)*indexNorm);
      }
      break;
    default:
      break;
    }
    this->m_sensitivities[iElement] = val;
  }
  return retval;
}

template <class T>
int ApertureData<T>::ApodizationSet(sps::unique_aligned_array<T> &&apodization,
                                    const size_t nElements) {
  int retval = -1;

  if (nElements == this->m_nelements) {
    this->m_sensitivities = std::move(apodization);
    retval = 0;
  }
  return retval;
}

template <class T>
int ApertureData<T>::Rotate(const sps::point_t<T>& reference,
                            const sps::euler_t<T>& euler) {
  SPS_UNREFERENCED_PARAMETERS(reference, euler);

  auto& elements = m_elements;

  sps::mat3_t<T> rotm0;
  sps::euler2rot<T, sps::EulerIntrinsicYXY>(euler, &rotm0);

  for (size_t iElement=0 ; iElement < m_nelements ; iElement++) {
    for (size_t jElement=0 ; jElement < m_nsubelements ; jElement++) {
      sps::element_rect_t<T>& element = elements[iElement][jElement];
      sps::point_t<T> p2e = element.center - reference;

      sps::point_t<T> tmp;
      sps::basis_rotate<T, sps::EulerIntrinsicYXY>(
        p2e, euler, &tmp);
      element.center = tmp + reference;

      // Update new Euler values
      sps::mat3_t<T> rotm1;
      sps::euler2rot<T, sps::EulerIntrinsicYXY>(element.euler, &rotm1);

      sps::mat3_t<T> rotm = rotm0*rotm1;

      sps::rot2euler<T, sps::EulerIntrinsicYXY>(rotm, &element.euler);
    }
  }

  initElements();
  initVectors();

  // TODO(Performace): Initialize this when needed only
  initRectangles();
  return 0;
}

template <class T>
void ApertureData<T>::ElementsSet(element_array &&elements,
                                  const size_t& nRows,
                                  const size_t& nCols) {
  if (nRows * nCols > 0) {
    assert(nCols == elements.m_n);
    this->m_npos         = nRows;
    this->m_nelements    = nRows;
    this->m_nsubelements = nCols;
    this->m_elements     = std::move(elements);

    initElements();
    initVectors();

    // TODO(Performace): Initialize this when needed only
    initRectangles();

    debug_print("Rectangles of first element\tx: %f %f %f %f\n",
                m_elements[0][0].vertices[0][0],
                m_elements[0][0].vertices[0][1],
                m_elements[0][0].vertices[0][2],
                m_elements[0][0].vertices[0][3]);
    debug_print("\t\t\t\t\ty: %f %f %f %f\n",
                m_elements[0][0].vertices[1][0],
                m_elements[0][0].vertices[1][1],
                m_elements[0][0].vertices[1][2],
                m_elements[0][0].vertices[1][3]);
    debug_print("\t\t\t\t\tz: %f %f %f %f\n",
                m_elements[0][0].vertices[2][0],
                m_elements[0][0].vertices[2][1],
                m_elements[0][0].vertices[2][2],
                m_elements[0][0].vertices[2][3]);
  }
}

template <class T>
void ApertureData<T>::ElementExtentGet(const size_t iElement,
                                       sps::bbox_t<T>* pBox) const {
  assert(iElement < m_nelements);
  *pBox = m_boxes[iElement];
}

template <class T>
void ApertureData<T>::ExtentGet(sps::bbox_t<T>* pBbox) const {
  sps::point_t<T> h_dir, w_dir;

  sps::bbox_t<T> bbox = sps::bbox_t<T>();
  // Coordinate of lower left position reference
  T pos_ll;

  const size_t nElements       = m_nelements;

  const auto& elements         = m_elements;

  for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
    for (size_t jElement = 0 ; jElement < m_nsubelements ; jElement++) {
      const sps::element_rect_t<T>& element = elements[iElement][jElement];

      sps::basis_vectors<T, sps::EulerIntrinsicYXY>(&w_dir, element.euler, 0);
      sps::basis_vectors<T, sps::EulerIntrinsicYXY>(&h_dir, element.euler, 1);

      for (size_t i_xyz=0 ; i_xyz < 3 ; i_xyz++) {
        // Corner position of rectangle (lower left)
        pos_ll =
          element.center[i_xyz]
          - element.hh*h_dir[i_xyz]
          - element.hw*w_dir[i_xyz];

        // Compute aperture bounding box
        if ((iElement == 0) && (jElement == 0)) {
          bbox.min[i_xyz] = pos_ll;
          bbox.max[i_xyz] = pos_ll;
        }

        bbox.min[i_xyz] = std::min<T>(bbox.min[i_xyz], pos_ll);
        bbox.max[i_xyz] = std::max<T>(bbox.max[i_xyz], pos_ll);

        T pos_ul = pos_ll+2*element.hh*h_dir[i_xyz];
        bbox.min[i_xyz] = std::min<T>(bbox.min[i_xyz], pos_ul);
        bbox.max[i_xyz] = std::max<T>(bbox.max[i_xyz], pos_ul);

        T pos_lr = pos_ll+2*element.hw*w_dir[i_xyz];
        bbox.min[i_xyz] = std::min<T>(bbox.min[i_xyz], pos_lr);
        bbox.max[i_xyz] = std::max<T>(bbox.max[i_xyz], pos_lr);

        T pos_ur = pos_ll+2*element.hh*h_dir[i_xyz]+2*element.hw*w_dir[i_xyz];
        bbox.min[i_xyz] = std::min<T>(bbox.min[i_xyz], pos_ur);
        bbox.max[i_xyz] = std::max<T>(bbox.max[i_xyz], pos_ur);
      }
    }
  }
  *pBbox = bbox;
}

template<class T>
T ApertureData<T>::AreaGet() const {
  T area = T(0.0);
  const size_t nElements       = m_nelements;

  const auto& elements = m_elements;

  for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
    for (size_t jElement = 0 ; jElement < m_nsubelements ; jElement++) {
      const auto& element = elements[iElement][jElement];
      area = area + T(4.0) * element.hh * element.hw;
    }
  }
  return area;
}

template<class T>
void ApertureData<T>::initVectors() {
  // Inputs: elements[][]
  // Output: elements[][]

  const size_t nElements              = m_nelements;
  const size_t nSubElementsPerElement = m_nsubelements;

  element_array& elements = m_elements;

  // Will never fail. Consider removing m_nelements and m_nsubelements
  debug_print("nElements: %zu, nSubElementsPerElement: %zu, m: %zu, n: %zu\n",
              nElements, nSubElementsPerElement, elements.m_m, elements.m_n);
  assert(elements.m_m == nElements);
  assert(elements.m_n == nSubElementsPerElement);

  for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
    for (size_t jElement = 0 ; jElement < nSubElementsPerElement ; jElement++) {
      auto& element = elements[iElement][jElement];
      sps::basis_vectors(element.uvector,
                         element.vvector,
                         element.normal,
                         element.euler);
      // Should not be necessary
      element.uvector[3] = T(0.0);
      element.vvector[3] = T(0.0);
      element.normal[3]  = T(0.0);
    }
  }
}

template <class T>
void ApertureData<T>::initElements() {
  const size_t nElements              = m_nelements;
  const size_t nSubElements           = m_nsubelements;

  auto& elements = m_elements;

  m_sensitivities = sps::unique_aligned_array_create<T>(nElements);
  m_pos          = sps::unique_aligned_array_create<sps::point_t<T> >(nElements);
  m_phases       = sps::unique_aligned_array_create<T>(nElements*nSubElements);
  m_delays       = sps::unique_aligned_array_create<T>(nElements);

  for (size_t iElement=0 ; iElement < nElements ; iElement++) {
    // Set apodizations
    m_sensitivities[iElement] = T(1.0);

    // Set delays
    m_delays[iElement] = T(0.0);

    // Center element (reference for delays)
    if (nSubElements % 2 == 1) {
      // Odd number of elements
      m_pos[iElement] = elements[iElement][nSubElements/2].center;
    } else {
      // Even number of elements
      m_pos[iElement] =
        T(0.5) * (elements[iElement][nSubElements-1].center + elements[iElement][0].center);
    }
    for (size_t jElement = 0 ; jElement < nSubElements ; jElement++) {
      // Set phases - do we need a phase per sub-element?
      m_phases[iElement*nSubElements + jElement] = T(0.0);
    }
  }
}

template <class T>
void ApertureData<T>::initRectangles() {
  const size_t nElements    = m_nelements;
  const size_t nSubElements = m_nsubelements;

  // Not used right now
  const size_t nSubH           = 1;
  const size_t nSubW           = 1;
  const size_t nSubSubElements = nSubH*nSubW;

  auto& elements = m_elements;

  sps::point_t<T> h_dir, w_dir, normal;

  // Delta height and width
  T dh, dw;

  // Lower left position reference
  T pos_ll;

  // Maximum width or height used for memory allocation
  m_h_xyz[0] = m_h_xyz[1] = m_h_xyz[2] = T(0.0);

  m_rectangles =
    sps::unique_aligned_array_create<sps::rect_t<T> >(nElements*nSubElements*nSubSubElements);

  m_boxes      = sps::unique_aligned_array_create<sps::bbox_t<T> >(nElements);

  for (size_t iElement=0 ; iElement < nElements ; iElement++) {
    for (size_t jElement=0 ; jElement < nSubElements ; jElement++) {
      sps::element_rect_t<T>& element = elements[iElement][jElement];

      sps::euler_t<T> euler;
      memcpy(&euler, &element.euler, sizeof(sps::euler_t<T>));

      basis_vectors(static_cast<T*>(&w_dir[0]),
                    static_cast<T*>(&h_dir[0]),
                    static_cast<T*>(&normal[0]), euler);
      dh = 2*element.hh/nSubH;
      dw = 2*element.hw/nSubW;

      for (size_t i_xyz=0 ; i_xyz < 3 ; i_xyz++) {
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
        for (size_t i_subh=0 ; i_subh < nSubH ; i_subh++) {
          for (size_t i_subw=0 ; i_subw < nSubW ; i_subw++) {
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

        m_boxes[iElement].min[i_xyz] =
          std::min<T>(m_boxes[iElement].min[i_xyz], pos_ur);
        m_boxes[iElement].min[i_xyz] =
          std::min<T>(m_boxes[iElement].min[i_xyz], pos_ul);
        m_boxes[iElement].min[i_xyz] =
          std::min<T>(m_boxes[iElement].min[i_xyz], pos_ll);
        m_boxes[iElement].min[i_xyz] =
          std::min<T>(m_boxes[iElement].min[i_xyz], pos_lr);
        m_boxes[iElement].max[i_xyz] =
          std::max<T>(m_boxes[iElement].max[i_xyz], pos_ur);
        m_boxes[iElement].max[i_xyz] =
          std::max<T>(m_boxes[iElement].max[i_xyz], pos_ul);
        m_boxes[iElement].max[i_xyz] =
          std::max<T>(m_boxes[iElement].max[i_xyz], pos_ll);
        m_boxes[iElement].max[i_xyz] =
          std::max<T>(m_boxes[iElement].max[i_xyz], pos_lr);
      }  // for i_xyz
    }  // jElement
  }  // iElement
}

// If this was header-only, we could export it
template class ApertureData<float>;

#ifdef FNM_DOUBLE_SUPPORT
template class ApertureData<double>;
#endif
}  // namespace fnm

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
