#include <fnm/fnm_types.hpp>
#include <sps/memory>
#include <sps/smath.hpp>  // sps::point_t and sps::euler_t, now element_rect_t

// Explicit instantiate destructor (otherwise not called)
template int sps::Singleton<fnm::bla<float> >::InstanceDestroy();

namespace fnm {

template <class T>
bla<T>::bla() {
  k = T(1.0);
}

template <class T>
sysparm_t<T>::sysparm_t() {
  /*
    If initialized outside using designated initializer this
    requires C99 and not strict ansi, i.e.

    # if defined(C99) && !defined(__STRICT_ANSI__)
  */

  // Default parameters
  c = T(1540);
  fs = T(100e6);
  f0 = T(7.0e6);
  att   = T(0.0);
  beta  = T(0.0);
  use_att = false;
  rho = T(1000.0);

  nThreads = 1;

  nDivW = 16;
  nDivH = 16;
  nDivA = 16;
  nMaxSectors = 10;
  w     = T(10.0) / T(100e6);
#if FNM_PULSED_WAVE
  /// Time-domain calculation type
  timeDomainCalcType = sofus::TimeDomainCalcType::Sphere;
  timeDomainIntOrder = sofus::TimeDomainIntOrder::Fourth;
  enforcePositiveDelays = false;
  soft_baffle = false;
#endif
}
}  // namespace fnm

namespace fnm {
template struct bla<float>;

template struct sysparm_t<float>;
template struct sysparm_t<double>;
}  // namespace fnm

namespace sps {
namespace nix {
template class unique_aligned_multi_array<sps::element_rect_t<float>, 2U>;
template class unique_aligned_multi_array<sps::element_rect_t<double>, 2U>;
}  // namespace nix
}  // namespace sps
