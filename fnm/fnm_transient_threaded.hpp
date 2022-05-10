/**
 * @file   fnm_transient_threaded.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Mon Sep 25 21:17:39 2017
 *
 * @brief
 *
 *
 */
#pragma once

#include <fnm/config.h>
#include <fnm/fnm_types.hpp>

namespace fnm {
  template<class T>
  class ApertureData;
}

namespace sps {
  class ProgressBarInterface;
}

namespace fnm {
  template <class T,  template <typename> class A>
  T CalcPwFnmThreaded(const sysparm_t<T>* sysparm,
                      const ApertureData<T>* data,
                      const T* pos, const size_t nPositions,
                      T** odata, size_t* nSamples,
                      int mask,
                      sps::ProgressBarInterface* pBar);

#if defined(HAVE_PTHREAD_H)
  template <class T,  template <typename> class A>
  void* CalcPwFnmThreadFunc(void* ptarg);
#else
  template <class T>
  unsigned int __stdcall CalcPwFnmThreadFunc(void *ptarg);
#endif

}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
