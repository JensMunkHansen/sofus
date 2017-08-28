/**
 * @file   circular_data.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Mon Jul 24 21:41:54 2017
 *
 * @brief
 *
 *
 */
#pragma once

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

#include <fnm/fnm_export.h>
#include <fnm/config.h>

#include <sps/cenv.h>
#include <sps/memory>
#include <sps/smath_types.hpp> // sps::point_t and sps::euler_t, now element_rect_t

#include <fnm/fnm_types.hpp> // FocusingType, element_rect_t

namespace fnm {

  template <class T>
  struct CircularApertureData : sps::aligned<4*sizeof(T)> {
    /// Next unique identifier to use
    static size_t nextID;

    sps::element_circular_t<T> m_element;

    /// Apodization
    T m_apodization;

    /// Phase, range is (-pi;pi]
    T m_phase;

    /// Delay
    T m_delay;

    /// Center frequency
    T m_f0;

    /// Grid sectoring scala
    T m_gridSectorScale;

    //@{  Internal state variables

    /// Unique identifier
    size_t m_id;

    //@}

    /// Ctor
    CircularApertureData();

    /**
     * Dtor
     *
     *
     * @return
     */
    ~CircularApertureData();

    /**
     * Consider making this private and expose ctor with r as parameter
     *
     */
    void Initialize();
    private:


    // Note: Deep copy is needed, if it is made copyable.
    CircularApertureData(const CircularApertureData& other) = delete;
    CircularApertureData& operator=(const CircularApertureData& other) = delete;
    CircularApertureData(CircularApertureData&& other) = delete;
    CircularApertureData& operator=(CircularApertureData&& other) = delete;
  };
}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
