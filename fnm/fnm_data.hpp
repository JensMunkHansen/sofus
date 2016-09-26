/**
 * @file   fnm_data.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Thu Jun  9 06:12:31 2016
 *
 * @brief
 *
 *
 */
#pragma once

#include <sps/cenv.h>

#include <fnm/fnm_export.h>
#include <fnm/FnmMath.hpp>

#include <vector>

namespace fnm {

  /*! \brief Aperture data structure
   *
   *
   * A struct containing data of the aperture
   */
#ifdef _MSC_VER
  template <class T>
  struct FNM_EXPORT ApertureData
#else
  template <class T>
  ALIGN16_BEGIN struct ALIGN16_END FNM_EXPORT ApertureData
#endif
  {

    /// Elements
    std::vector< std::vector<element_t<T> > >* m_elements;

    /// Positions
    sps::point_t<T>* m_pos;

    /// Apodization values
    T* m_apodizations;

    /// Phases
    T* m_phases;

    /// Number of elements
    size_t m_nelements;

    /// Number of sub-elements per element
    size_t m_nsubelements;

    /// Rectangles
    sps::rect_t<T>* m_rectangles;

    /// Number of positions
    size_t m_npos;

    /// Center frequency
    T m_f0;

    /// Focus point
    sps::point_t<T> m_focus;

    /// Validity of focus
    bool _focus_valid;

    /// Next unique identifier to use
    static size_t nextID;

    /// Unique identifier
    size_t m_id;

    /// Ctor
    ApertureData();

    /// Dtor
    ~ApertureData();

  };
}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
