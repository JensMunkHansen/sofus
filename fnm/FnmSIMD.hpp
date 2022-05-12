/**
 * @file   FnmSIMD.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Thu Jun  9 06:01:12 2016
 *
 * @brief
 *
 *
 */

// TODO: Make HzFast twice as fast by reducing integrals to four

#pragma once

#include <fnm/config.h>
#include <sps/cenv.h>
#include <sps/smath.hpp>
#include <fnm/fnm_types.hpp>
#include <sps/trigintrin.h>

#include <complex>

/**
 *
 *
 * @param s
 * @param l
 * @param z
 * @param k
 * @param us
 * @param uweights
 * @param nUs
 * @param vs
 * @param vweights
 * @param nVs
 *
 * @return
 */
template <class T>
STATIC_INLINE_BEGIN std::complex<T> CalcHzVecGL(const T& s,
    const T& l,
    const T& z,
    const T& k,
    const T* us,
    const T* uweights,
    const size_t nUs,
    const T* vs,
    const T* vweights,
    const size_t nVs) STATIC_INLINE_END;


/**
 *
 *
 * @param element
 * @param projection
 * @param k
 * @param us
 * @param uweights
 * @param nUs
 * @param vs
 * @param vweights
 * @param nVs
 *
 * @return
 */template <class T>
STATIC_INLINE_BEGIN
std::complex<T> CalcHzAll(const sps::element_rect_t<T>& element,
                          const sps::point_t<T>& projection,
                          const T& k,
                          const T* us,
                          const T* uweights,
                          const size_t nUs,
                          const T* vs,
                          const T* vweights,
                          const size_t nVs) STATIC_INLINE_END;


/**
 * Integral with (nUs x nVs) points with projection inside the region of integration
 *
 * @param element
 * @param projection
 * @param k
 * @param us
 * @param uweights
 * @param nUs
 * @param vs
 * @param vweights
 * @param nVs
 *
 * @return
 */
template <class T>
STATIC_INLINE_BEGIN
std::complex<T> CalcHzFast(const sps::element_rect_t<T> &__restrict element,
                           const sps::point_t<T> &__restrict projection,
                           const T &__restrict k,
                           const T* __restrict us,
                           const T* __restrict uweights,
                           const size_t nUs,
                           const T* __restrict vs,
                           const T* __restrict vweights,
                           const size_t nVs) STATIC_INLINE_END;

template <class T>
STATIC_INLINE_BEGIN
std::complex<T> CalcFastFourAny(const T& u,
                                const T& v,
                                const T& hh,
                                const T& hw,
                                const T& z,
                                const T& __restrict k,
                                const T* __restrict s,
                                const T* __restrict weights,
                                const size_t nS) STATIC_INLINE_END;

template <class T>
STATIC_INLINE_BEGIN
std::complex<T> CalcFastFourAny2(const T& u,
                                 const T& v,
                                 const T& hh,
                                 const T& hw,
                                 const T& z,
                                 const T& __restrict k,
                                 const T* __restrict s,
                                 const T* __restrict weights,
                                 const size_t nS) STATIC_INLINE_END;


template <class T>
STATIC_INLINE_BEGIN
std::complex<T> CalcFourFast(const sps::element_rect_t<T> &__restrict element,
                             const sps::point_t<T> &__restrict projection,
                             const T &__restrict k,
                             const T* __restrict uvs,
                             const T* __restrict uvweights,
                             const size_t nUVs) STATIC_INLINE_END;

template <class T>
STATIC_INLINE_BEGIN
std::complex<T> CalcHzFastSingle(const sps::element_rect_t<T>& element,
                                 const sps::point_t<T>& projection,
                                 const T& k,
                                 const T* us,
                                 const T* uweights,
                                 const size_t nUs,
                                 const T* vs,
                                 const T* vweights,
                                 const size_t nVs) STATIC_INLINE_END;

template <>
std::complex<float>
inline CalcHzVecGL(const float& s,
                   const float& l,
                   const float& z,
                   const float& k,
                   const float* us,
                   const float* uweights,
                   const size_t nUs,
                   const float* vs,
                   const float* vweights,
                   const size_t nVs)
{

  const __m128 carg = _mm_set1_ps(cos(-k*z));
  const __m128 sarg = _mm_set1_ps(sin(-k*z));

  const __m128 l_2 = _mm_mul_ps(_mm_set1_ps(l),_m_half_ps);
  const __m128 s_2 = _mm_mul_ps(_mm_set1_ps(s),_m_half_ps);

  const __m128 z2 = _mm_set1_ps(SQUARE(z));
  const __m128 l2 = _mm_set1_ps(SQUARE(l));
  const __m128 s2 = _mm_set1_ps(SQUARE(s));

  __m128 intWreal = _mm_setzero_ps();
  __m128 intWimag = _mm_setzero_ps();

  __m128 real, imag;

  for (size_t iu = 0 ; iu < 4*((nUs+3)/4) ; iu+=4) {

    __m128 ls  = _mm_add_ps(_mm_mul_ps(l_2,_mm_load_ps((float*)&us[iu])),l_2);
    __m128 ls2 = _mm_square_ps(ls);

    __m128 argw = _mm_mul_ps(_mm_set1_ps(-k),_mm_sqrt_ps(_mm_add_ps(_mm_add_ps(ls2,z2),s2)));

    __m128 cargw, sargw;

    _mm_sin_cos_ps(argw, &sargw, &cargw);

    __m128 vec_uweight = _mm_load_ps((float*)&uweights[iu]);

    __m128 denom = _mm_add_ps(ls2,s2);
    __m128 rcp_denom = _mm_rcp_ps(denom);

    real = _mm_mul_ps(_mm_mul_ps(vec_uweight,_mm_sub_ps(cargw, carg)),rcp_denom);
    imag = _mm_mul_ps(_mm_mul_ps(vec_uweight,_mm_sub_ps(sargw, sarg)),rcp_denom);

    // Update integral
    intWreal = _mm_add_ps(intWreal, real);
    intWimag = _mm_add_ps(intWimag, imag);
  }

  __m128 rcp_denom1 = _mm_rcp_ps(_mm_mul_ps(_m_2pi_ps,_mm_set1_ps(k)));

  intWreal = _mm_mul_ps(intWreal, _mm_mul_ps(_mm_mul_ps(l_2,_mm_set1_ps(s)), rcp_denom1));
  intWimag = _mm_mul_ps(intWimag, _mm_mul_ps(_mm_mul_ps(l_2,_mm_set1_ps(s)), rcp_denom1));

  // integral height
  std::complex<float> intH = std::complex<float>(float(0.0),float(0.0));

  __m128 intHreal = _mm_setzero_ps();
  __m128 intHimag = _mm_setzero_ps();

  for(size_t iv = 0 ; iv < 4*((nVs+3)/4) ; iv+=4) {

    __m128 ss  = _mm_add_ps(_mm_mul_ps(s_2,_mm_load_ps((float*)&vs[iv])),s_2);
    __m128 ss2 = _mm_square_ps(ss);

    __m128 argh = _mm_mul_ps(_mm_set1_ps(-k),_mm_sqrt_ps(_mm_add_ps(_mm_add_ps(ss2,z2),l2)));

    __m128 cargh, sargh;

    _mm_sin_cos_ps(argh, &sargh, &cargh);

    __m128 vec_vweight = _mm_load_ps((float*)&vweights[iv]);
    __m128 rcp_denom = _mm_rcp_ps(_mm_add_ps(ss2,l2));

    real = _mm_mul_ps(_mm_mul_ps(vec_vweight,_mm_sub_ps(cargh, carg)),rcp_denom);
    imag = _mm_mul_ps(_mm_mul_ps(vec_vweight,_mm_sub_ps(sargh, sarg)),rcp_denom);
    intHreal = _mm_add_ps(intHreal, real);
    intHimag = _mm_add_ps(intHimag, imag);
  }

  // Divide by denominator
  intHreal = _mm_mul_ps(intHreal, _mm_mul_ps(_mm_mul_ps(s_2,_mm_set1_ps(l)), rcp_denom1));
  intHimag = _mm_mul_ps(intHimag, _mm_mul_ps(_mm_mul_ps(s_2,_mm_set1_ps(l)), rcp_denom1));

  intHreal = _mm_add_ps(intHreal,intWreal);
  intHimag = _mm_add_ps(intHimag,intWimag);

  // Multiply by -i
  __m128 tmp = intHreal;
  intHreal = intHimag;
  intHimag = _mm_neg_ps(tmp);

  // Sum 4 partial integrals
  __m128 result = _mm_dp_ps(_m_one_ps,intHreal,0xF1);
  result = _mm_add_ps(result,_mm_dp_ps(_m_one_ps,intHimag,0xF2));

  ALIGN16_BEGIN float results[4] ALIGN16_END;

  _mm_store_ps(results,result);

  intH.real(results[0]);
  intH.imag(results[1]);

  return intH;
}

template <>
std::complex<float>
inline CalcHzAll(const sps::element_rect_t<float>& element,
                 const sps::point_t<float>& projection, // Consider 4 points
                 const float& k,
                 const float* us,
                 const float* uweights,
                 const size_t nUs,
                 const float* vs,
                 const float* vweights,
                 const size_t nVs)
{

  std::complex<float> retval;

  const float z = projection[2];

  float s0 = fabs(projection[1]) + element.hh;
  float s2 = element.hh - fabs(projection[1]);
  float l0 = fabs(projection[0]) + element.hw;
  float l1 = element.hw - fabs(projection[0]);

  __m128 s = _mm_set_ps(s2,s2,s0,s0);
  __m128 l = _mm_set_ps(l1,l0,l1,l0);

  const __m128 vec_s = _mm_fabs_ps(s);
  const __m128 vec_l = _mm_fabs_ps(l);
  const __m128 cargz = _mm_set1_ps(cos(-k*z));
  const __m128 sargz = _mm_set1_ps(sin(-k*z));

  const __m128 vec_l_2 = _mm_mul_ps(vec_l,_m_half_ps);
  const __m128 vec_s_2 = _mm_mul_ps(vec_s,_m_half_ps);

  const __m128 z2 = _mm_set1_ps(SQUARE(z));
  const __m128 vec_l2 = _mm_square_ps(vec_l);
  const __m128 vec_s2 = _mm_square_ps(vec_s);

  __m128 real, imag;

  __m128 intWreal = _mm_setzero_ps();
  __m128 intWimag = _mm_setzero_ps();

  for (size_t iu = 0 ; iu < nUs ; iu++) {

    __m128 us1       = _mm_load1_ps((float*)&us[iu]);
    __m128 uweights1 = _mm_load1_ps((float*)&uweights[iu]);

    __m128 ls  = _mm_add_ps(_mm_mul_ps(vec_l_2,us1),vec_l_2);
    __m128 ls2 = _mm_square_ps(ls);

    __m128 argw = _mm_mul_ps(
                    _mm_set1_ps(-k),
                    _mm_sqrt_ps(
                      _mm_add_ps(
                        _mm_add_ps(
                          ls2,
                          z2),
                        vec_s2)));

    __m128 cargw, sargw;

    _mm_sin_cos_ps(argw, &sargw, &cargw);

    __m128 rcp_denom = _mm_rcp_ps(_mm_add_ps(ls2,vec_s2));
    real = _mm_mul_ps(_mm_mul_ps(uweights1,_mm_sub_ps(cargw, cargz)),rcp_denom);
    imag = _mm_mul_ps(_mm_mul_ps(uweights1,_mm_sub_ps(sargw, sargz)),rcp_denom);
    intWreal = _mm_add_ps(intWreal, real);
    intWimag = _mm_add_ps(intWimag, imag);
  }

  __m128 rcp_denom1 = _mm_rcp_ps(_mm_mul_ps(_m_2pi_ps,_mm_set1_ps(k)));

  intWreal = _mm_mul_ps(
               intWreal,
               _mm_mul_ps(
                 _mm_mul_ps(
                   vec_l_2,
                   vec_s),
                 rcp_denom1));
  intWimag = _mm_mul_ps(
               intWimag,
               _mm_mul_ps(
                 _mm_mul_ps(
                   vec_l_2,
                   vec_s),
                 rcp_denom1));

  __m128 intHreal = _mm_setzero_ps();
  __m128 intHimag = _mm_setzero_ps();

  for(size_t iv = 0 ; iv < nVs ; iv++) {

    __m128 vs1       = _mm_load1_ps((float*)&vs[iv]);
    __m128 vweights1 = _mm_load1_ps((float*)&vweights[iv]);

    __m128 ss  = _mm_add_ps(_mm_mul_ps(vec_s_2,vs1),vec_s_2);
    __m128 ss2 = _mm_square_ps(ss);

    __m128 argh = _mm_mul_ps(
                    _mm_set1_ps(-k),
                    _mm_sqrt_ps(
                      _mm_add_ps(
                        _mm_add_ps(
                          ss2,
                          z2),
                        vec_l2)));

    __m128 cargh, sargh;

    _mm_sin_cos_ps(argh, &sargh, &cargh);

    __m128 rcp_denom = _mm_rcp_ps(_mm_add_ps(ss2,vec_l2));

    real = _mm_mul_ps(_mm_mul_ps(vweights1,_mm_sub_ps(cargh, cargz)),rcp_denom);
    imag = _mm_mul_ps(_mm_mul_ps(vweights1,_mm_sub_ps(sargh, sargz)),rcp_denom);
    intHreal = _mm_add_ps(intHreal, real);
    intHimag = _mm_add_ps(intHimag, imag);
  }

  // Divide by denominator
  intHreal = _mm_mul_ps(intHreal, _mm_mul_ps(_mm_mul_ps(vec_s_2,vec_l), rcp_denom1));
  intHimag = _mm_mul_ps(intHimag, _mm_mul_ps(_mm_mul_ps(vec_s_2,vec_l), rcp_denom1));

  intHreal = _mm_add_ps(intHreal,intWreal);
  intHimag = _mm_add_ps(intHimag,intWimag);

  // Multiply by -i
  __m128 tmp = intHreal;
  intHreal = intHimag;
  intHimag = _mm_neg_ps(tmp);

  // Filter
  __m128 sign = _mm_mul_ps(s,l);
  intHreal = _mm_mulsign_ps(intHreal,sign);
  intHimag = _mm_mulsign_ps(intHimag,sign);

  // Horizontal sum
  _mm_store_ss(&(reinterpret_cast<float(&)[2]>(retval)[0]),_mm_dp_ps(_m_one_ps,intHreal,0xF1));
  _mm_store_ss(&(reinterpret_cast<float(&)[2]>(retval)[1]),_mm_dp_ps(_m_one_ps,intHimag,0xF1));

  return retval;
}


// Not working accurately!!!

template <>
std::complex<float>
inline CalcHzFastSingle(const sps::element_rect_t<float>& element,
                        const sps::point_t<float>& projection,
                        const float& k,
                        const float* us,
                        const float* uweights,
                        const size_t nUs,
                        const float* vs,
                        const float* vweights,
                        const size_t nVs)
{
  SPS_UNREFERENCED_PARAMETERS(element,projection, k, us, uweights, nUs,
                              vs, vweights, nVs);
  return std::complex<float>();
}

template <>
std::complex<double>
inline CalcHzFastSingle(const sps::element_rect_t<double>& element,
                        const sps::point_t<double>& projection,
                        const double& k,
                        const double* us,
                        const double* uweights,
                        const size_t nUs,
                        const double* vs,
                        const double* vweights,
                        const size_t nVs)
{
  SPS_UNREFERENCED_PARAMETERS(element,projection, k, us, uweights, nUs,
                              vs, vweights, nVs);
  return std::complex<double>();
}

template <>
std::complex<double>
inline CalcHzVecGL(const double& s,
                   const double& l,
                   const double& z,
                   const double& k,
                   const double* us,
                   const double* uweights,
                   const size_t nUs,
                   const double* vs,
                   const double* vweights,
                   const size_t nVs)
{
  SPS_UNREFERENCED_PARAMETERS(s,l,z,k,us,uweights,nUs,vs,vweights,nVs);
  return std::complex<double>();
}

template <>
std::complex<double>
inline CalcHzAll(const sps::element_rect_t<double>& element,
                 const sps::point_t<double>& projection,
                 const double& k,
                 const double* us,
                 const double* uweights,
                 const size_t nUs,
                 const double* vs,
                 const double* vweights,
                 const size_t nVs)
{
  SPS_UNREFERENCED_PARAMETERS(element,projection,k,us,uweights,nUs,vs,vweights,nVs);

  std::complex<double> retval;
  return retval;
}

#include <fnm/fnm_ps.hpp>
#include <fnm/fnm_pd.hpp>

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
