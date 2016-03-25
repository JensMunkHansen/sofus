#pragma once

#include <fnm/fnm_export.h>
#include <fnm/FnmMath.hpp>
#include <sps/cenv.h>

#include <vector>

namespace fnm {

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

    sps::point_t<T>* m_pos;

    T* m_apodizations;

    T* m_phases;

    /// Number of elements
    size_t m_nelements;

    size_t m_nsubelements;

    sps::rect_t<T>* m_rectangles;
    
    /// Number of positions
    size_t m_npos;

    T m_f0;
    sps::point_t<T> m_focus;

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
