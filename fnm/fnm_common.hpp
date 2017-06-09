/**
 * @file   fnm_common.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Sat Mar 18 13:31:45 2017
 *
 * @brief  Utility functions for fast nearfield method
 *
 *
 */
#include <sps/memory>            // sps::deleted_aligned_array

namespace fnm {

  /// Forward-declare sysparm_t
  template<class T>
  struct sysparm_t;

  /**
   *
   * @param sysparm
   * @param uxs
   * @param uweights
   * @param vxs
   * @param vweights
   */
  template <class T>
  void CalcWeightsAndAbcissae(const sysparm_t<T>* sysparm,
                              sps::deleted_aligned_array<T> &&uxs,
                              sps::deleted_aligned_array<T> &&uweights,
                              sps::deleted_aligned_array<T> &&vxs,
                              sps::deleted_aligned_array<T> &&vweights);
}
