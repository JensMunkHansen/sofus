/**
 * @file   fnm_common.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Sat Mar 18 13:31:45 2017
 *
 * @brief  Utility functions for fast nearfield method
 *
 *
 */
#include <sps/memory>            // sps::unique_aligned_array

#include <sps/smath_types.hpp>
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
                            sps::unique_aligned_array<T> &&uxs,
                            sps::unique_aligned_array<T> &&uweights,
                            sps::unique_aligned_array<T> &&vxs,
                            sps::unique_aligned_array<T> &&vweights);

/**
 * Weight and abscissa are scaled, but abcissa are not shifted
 *
 * @param sysparm
 * @param element
 * @param uxs
 * @param uweights
 * @param vxs
 * @param vweights
 */
template <class T>
void CalcWeightsAndAbcissaeScaled(const sysparm_t<T>* sysparm,
                                  const sps::element_rect_t<T>& element,
                                  sps::unique_aligned_array<T> &&uxs,
                                  sps::unique_aligned_array<T> &&uweights,
                                  sps::unique_aligned_array<T> &&vxs,
                                  sps::unique_aligned_array<T> &&vweights);

template <class T>
void CalcWeightsAndAbcissaeSIMD(const sysparm_t<T>* sysparm,
                                T** __restrict uv_xs,
                                T** __restrict uv_ws);

template <class T>
void CalcWeightsAndAbcissaeScaledSIMD(const sysparm_t<T>* sysparm,
                                      const sps::element_rect_t<T>* pElement,
                                      sps::unique_aligned_array<T> &&uv_xs,
                                      sps::unique_aligned_array<T> &&uv_ws);

}
