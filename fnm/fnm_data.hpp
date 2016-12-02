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

#include <fnm/fnm_export.h>
#include <fnm/config.h>

#include <sps/cenv.h>
#include <sps/memory>

#include <fnm/fnm_types.hpp> // FocusingType, element_t

#include <vector>

namespace fnm {

  /*! \brief Aperture data structure
   *
   *
   * A struct containing data of the aperture
   */
#ifdef _MSC_VER
  template <class T>
  class FNM_EXPORT ApertureData
#else
  template <class T>
  ALIGN16_BEGIN struct ALIGN16_END FNM_EXPORT ApertureData
#endif
  {
  public:
    /// Next unique identifier to use
    static size_t nextID;

    /// Number of elements
    size_t m_nelements;

    /// Number of sub-elements per element
    size_t m_nsubelements;

    /// Number of positions, equals number of elements
    size_t m_npos;

    /// Elements
    sps::deleted_aligned_multi_array<sps::element_t<T>, 2U> m_elements;

    /// Positions
    sps::deleted_aligned_array<sps::point_t<T> > m_pos;

    /// Apodization values
    sps::deleted_aligned_array<T> m_apodizations;

    /// Phases, range is (-pi;pi]
    sps::deleted_aligned_array<T> m_phases;

    /// Delays
    sps::deleted_aligned_array<T> m_delays;

    /// Rectangles
    sps::deleted_aligned_array<sps::rect_t<T> > m_rectangles;

    /// Center frequency
    T m_f0;

    /// Focus point
    sps::point_t<T> m_focus;

    /// Focusing Type, Rayleigh or Pythagorean
    int m_focus_type;

    //@{  Internal state variables

    /// Validity of focus, if equal m_focus_type, phases are valid
    int m_focus_valid;

    /// Unique identifier
    size_t m_id;

    //@}
    /// Ctor
    ApertureData();

    void ElementsSet(sps::deleted_aligned_multi_array<sps::element_t<T>,2> &&elements, const size_t& nRows, const size_t& nCols);

    std::vector<std::vector<sps::element_t<T> > > ElementsVectorGet() const;

    void initVectors();

    void initElements();

    void initRectangles();

    void ExtentGet(sps::bbox_t<T>& bbox) const;

    T AreaGet() const;

    /// Dtor
    ~ApertureData();
  private:
    ApertureData(const ApertureData& other) = delete;
    ApertureData& operator=(const ApertureData& other) = delete;
    ApertureData(ApertureData&& other) = delete;
    ApertureData& operator=(ApertureData&& other) = delete;

    // Deep copy is needed, if it is made copyable.
    /*
      ApertureData(const ApertureData& other)
      : m_pos((point_t<T>*)_mm_malloc(other.m_npos*sizeof(point_t<T>),16), [](point_t<T>* f)->void { _mm_free(f);}),
        m_apodizations((T*)_mm_malloc(other.m_npos*sizeof(T),16), [](T* f)->void { _mm_free(f);}),
        m_phases((T*)_mm_malloc(other.m_npos*sizeof(T),16), [](T* f)->void { _mm_free(f);}),
        m_rectangles((rect_t<T>*)_mm_malloc(other.m_npos*sizeof(rect_t<T>),16), [](rect_t<T>* f)->void { _mm_free(f);})
      {
        memcpy(m_pos.get(),other.m_pos.get(),m_npos*sizeof(T));
        // etc
        return *this;
      }

      and std::move is needed, if it is made moveable

      ApertureData(ApertureData&& other) : m_pos(std::move(other.m_pos))
      {
        return *this;
      }
    */
  };
}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
