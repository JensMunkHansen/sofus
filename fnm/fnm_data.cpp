#include <fnm/config.h>
#include <fnm/fnm_data.hpp>
#include <sps/cerr.h>
#include <sps/mm_malloc.h>
#include <string.h>

#ifdef HAVE_MQUEUE_H
# include <mqueue.h>
#endif

namespace fnm {

  template <class T>
  size_t ApertureData<T>::nextID = 0;

  template <class T>
  ApertureData<T>::ApertureData() : m_elements(NULL), m_pos(NULL), m_apodizations(NULL), m_phases(NULL),
    m_nelements(0), m_nsubelements(0), m_rectangles(NULL), m_npos(0)
  {
    m_f0 = T(1e6);

    m_focus = sps::point_t<T>();
    memset(&this->m_focus[0],0,3*sizeof(T));

    _focus_valid = false;

    // Set unique identifier
    m_id = nextID;
    nextID++;

    // Default is one element (remove this)
    m_nelements = 1;
    m_nsubelements = 1;

    m_elements = new std::vector< std::vector<element_t<T> > >(m_nelements, std::vector<element_t<T> >(m_nsubelements));

    m_npos = m_nelements;
    m_pos = (sps::point_t<T>*) _mm_malloc(m_npos*sizeof(sps::point_t<T>),16);
    if (m_pos) {
      memset(&m_pos[0][0],0,3*m_npos*sizeof(T));
    }

    m_apodizations = (T*) _mm_malloc(m_npos*sizeof(T),16);
    if (m_apodizations) {
      m_apodizations[0] = T(1.0);
    }

    m_phases = (T*) _mm_malloc(m_npos*sizeof(T),16);

  }

  template <class T>
  ApertureData<T>::~ApertureData()
  {
    if (this->m_elements) {
      delete this->m_elements;
      this->m_elements = NULL;
    }
    if (this->m_apodizations) {
      _mm_free(this->m_apodizations);
      this->m_apodizations = NULL;
    }
    if (this->m_phases) {
      _mm_free(this->m_phases);
      this->m_phases = NULL;
    }
    if(this->m_pos) {
      _mm_free(this->m_pos);
      this->m_pos = NULL;
    }
  }
  template struct ApertureData<float>;
  template struct ApertureData<double>;
}
/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
