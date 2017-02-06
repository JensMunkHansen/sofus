#include <sps/memory>            // sps::deleted_aligned_array_create

namespace fnm {
  template <class T>
  void CalcWeightsAndAbcissae(sps::deleted_aligned_array<T> &&uxs,
                              sps::deleted_aligned_array<T> &&uweights,
                              sps::deleted_aligned_array<T> &&vxs,
                              sps::deleted_aligned_array<T> &&vweights);
}
