/**
 * @file   fnm_response.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Tue Sep 26 19:42:59 2017
 *
 * @brief
 *
 *
 */
#pragma once

#include <sps/smath_types.hpp>

namespace fnm {
template <class T>
class ApertureData;

template <class T>
struct GLQuad2D;

template <class T>
struct sysparm_t;
}  // namespace fnm

namespace sps {
template <class T>
class msignal1D;
}  // namespace sps

namespace fnm {

template <typename T, template <typename> class A>
void FnmResponse(const sysparm_t<T>* pSysparm,
                 const ApertureData<T>* pData,
                 const GLQuad2D<T>* pGL,
                 const T& amplitude,
                 const sps::point_t<T>& point,
                 const int& iSampleSignalStart,
                 const size_t& nSamples,
                 T* pOdata,
                 const int mask = 0x1F);
}  // namespace fnm

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
