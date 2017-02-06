#include <fnm/config.h>
#include <fnm/fnm_export.h>

#include <cstdlib>

#include <sps/memory>
#include <sps/smath.hpp>

namespace fnm {

  /**
   * Compute focused linear array
   *
   * @param nElements
   * @param nSubH
   * @param nSubW
   * @param pitch
   * @param kerf
   * @param height
   * @param eFocus
   * @param arcPlacement 0 = outside, 1 = inside
   * @param elements output
   */
  template <class T>
  void FNM_EXPORT FocusedLinearArray(const size_t nElements,
                                     const size_t nSubH, // Elevation
                                     const size_t nSubW,
                                     const T pitch,
                                     const T kerf,
                                     const T height,
                                     const T efocus,
                                     const int arcPlacement,
                                     sps::deleted_aligned_multi_array<sps::element_t<T>, 2> &&elements);

  template <class T>
  void FNM_EXPORT FocusedConvexArray(const size_t nElements,
                                     const size_t nSubH, // Elevation
                                     const size_t nSubW,
                                     const T pitch,
                                     const T kerf,
                                     const T height,
                                     const T radius,
                                     const T efocus,
                                     const int arcPlacement,
                                     sps::deleted_aligned_multi_array<sps::element_t<T>, 2> &&elements);

}
/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
