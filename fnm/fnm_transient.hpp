/**
 * @file   fnm_transient.hpp
 * @author Jens Munk Hansen <jmh@jmhlaptop>
 * @date   Thu May  4 20:25:38 2017
 *
 * @brief
 *
 *
 */
#pragma once
#include <fnm/fnm_types.hpp>

namespace fnm {
#ifdef _MSC_VER
  template<class T>
  class ApertureData;
#else
  template<class T>
  struct ApertureData;
#endif
}

namespace fnm {

  template <class T>
  T CalcFdTransientRef(const sysparm_t<T>* sysparm,
                       const ApertureData<T>* data,
                       const T* pos, const size_t nPositions, const size_t nDim,
                       T** odata, size_t* nSignals, size_t* nSamples);
}
