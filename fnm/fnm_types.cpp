#include <fnm/fnm_types.hpp>
#include <sps/memory>
#include <sps/smath.hpp> // sps::point_t and sps::euler_t, now element_rect_t

namespace fnm {
  template <class T>
  sysparm_t<T>::sysparm_t()
  {
    /*
      If initialized outside using designated initializer this
      requires C99 and not strict ansi, i.e.

      # if defined(C99) && !defined(__STRICT_ANSI__)
    */

    // Default parameters
    c = T(1500);
    fs = T(100e6);
    att   = T(0.0);
    beta  = T(0.0);
    use_att = false;

    nDivW = 16;
    nDivH = 16;
    nDivA = 16;
	nMaxSectors = 10;
    w     = T(10.0) / T(100e6);
#if FNM_PULSED_WAVE
    /// Time-domain calculation type
    timeDomainCalcType = sofus::TimeDomainCalcType::Sphere;
    timeDomainIntOrder = sofus::TimeDomainIntOrder::Fourth;
#endif
  }
}
namespace fnm {
  template struct sysparm_t<float>;
  template struct sysparm_t<double>;
}

namespace sps {
  namespace nix {
    template class deleted_aligned_multi_array<sps::element_rect_t<float>, 2U>;
    template class deleted_aligned_multi_array<sps::element_rect_t<double>, 2U>;
  }
}

