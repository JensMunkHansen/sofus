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

#include <fnm/circular_data.hpp>
#include <sps/smath.hpp>

namespace fnm {

  template <class T>
  size_t CircularApertureData<T>::nextID = 0;

  template <class T>
  CircularApertureData<T>::CircularApertureData()
  {
    m_element.circle.radius = T(1.0);
    m_element.circle.center[0] = T(0.0);
    m_element.circle.center[1] = T(0.0);
    m_element.circle.center[2] = T(0.0);
    m_element.circle.euler.alpha = T(0.0);
    m_element.circle.euler.beta  = T(0.0);
    m_element.circle.euler.gamma = T(0.0);

    m_apodization = T(1.0);
    m_phase = T(0.0);
    m_delay = T(0.0);
    m_f0 = T(1e6);
    m_gridSectorScale = T(1.0);
    m_id = nextID;
    nextID++;
  }

  template <class T>
  CircularApertureData<T>::~CircularApertureData()
  {
  }

  template <class T>
  void CircularApertureData<T>::Initialize()
  {

    //    ALIGN32_BEGIN sps::point_t<T> h_dir ALIGN32_END;
    //    ALIGN32_BEGIN sps::point_t<T> w_dir ALIGN32_END;

    sps::basis_vectors((T*)&m_element.uvector[0],
                       (T*)&m_element.vvector[0],
                       (T*)&m_element.normal[0],
                       m_element.circle.euler);
  }

  // If this was header-only, we could export it
  template struct CircularApertureData<float>;
#ifdef FNM_DOUBLE_SUPPORT
  template struct CircularApertureData<double>;
#endif
}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
