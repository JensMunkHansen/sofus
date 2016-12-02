#include <fnm/fnm_types.hpp>

#include <sps/memory>

namespace fnm {
  //  template struct element_t<float>;
  //  template struct element_t<double>;
  template struct sysparm_t<float>;
  template struct sysparm_t<double>;
}

namespace sps {
  namespace nix {
    template class deleted_aligned_multi_array<sps::element_t<float>, 2U>;
    template class deleted_aligned_multi_array<sps::element_t<double>, 2U>;
  }
}
