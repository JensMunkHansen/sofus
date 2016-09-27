/**
 * @file   FnmSIMD.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Thu Jun  9 06:01:12 2016
 *
 * @brief
 *
 *
 */

#pragma once

// Example of wrong include
#include <sps/cenv.h>
#include <sps/trigintrin.h>
//#include <sps/smath.hpp>

#include <fnm/FnmMath.hpp>
#include <fnm/config.h>

#include <complex>

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


// Non-reduced integral, but simple
template <class T>
STATIC_INLINE_BEGIN std::complex<T> CalcHzAll(const fnm::element_t<T>& element,
    const sps::point_t<T>& projection,
    const T& k,
    const T* us,
    const T* uweights,
    const size_t nUs,
    const T* vs,
    const T* vweights,
    const size_t nVs) STATIC_INLINE_END;


template <class T>
STATIC_INLINE_BEGIN std::complex<T> CalcHzShort(const fnm::element_t<T>& element,
    const sps::point_t<T>& projection,
    const T& k,
    const T* us,
    const T* uweights,
    const size_t nUs,
    const T* vs,
    const T* vweights,
    const size_t nVs) STATIC_INLINE_END;

/**
 * Integral with (nUs x nVs) points inside the region of integration
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
STATIC_INLINE_BEGIN std::complex<T> CalcHzFast(const fnm::element_t<T>& element,
    const sps::point_t<T>& projection,
    const T& k,
    const T* us,
    const T* uweights,
    const size_t nUs,
    const T* vs,
    const T* vweights,
    const size_t nVs) STATIC_INLINE_END;

template <class T>
STATIC_INLINE_BEGIN std::complex<T> CalcHzFast4(const fnm::element_t<T>& element,
    const sps::point_t<T>& projection,
    const T& k,
    const T* us,
    const T* uweights,
    const size_t nUs,
    const T* vs,
    const T* vweights,
    const size_t nVs) STATIC_INLINE_END;

template <>
std::complex<float> CalcHzVecGL(const float& s,
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

    __m128 rcp_denom = _mm_rcp_ps(_mm_add_ps(ls2,s2));
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
CalcHzAll(const fnm::element_t<float>& element,
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


template <>
std::complex<float>
CalcHzShort(const fnm::element_t<float>& element,
            const sps::point_t<float>& projection,
            const float& k,
            const float* us,
            const float* uweights,
            const size_t nUs,
            const float* vs,
            const float* vweights,
            const size_t nVs)
{

  std::complex<float> retval = std::complex<float>(0.0f,0.0f);

  const float x = projection[0];
  const float y = projection[1];
  const float z = projection[2];

  const float hw = element.hw;
  const float hh = element.hh;

  float s0 = fabs(x) + hw; // [0;|x|+hw] -> [0;hw+|x|]  -> [|x|;|x|+hw]
  float s1 = hw - fabs(x); // [0;hw-|x|] -> -[hw-|x|;0] -> [|x|-hw;|x|]
  float l0 = fabs(y) + hh;
  float l2 = hh - fabs(y);


  // Integral (stop, inside)
  __m128 s_start_in = _mm_setzero_ps();
  __m128 l_start_in = _mm_setzero_ps();

  // Negative
  __m128 s_stop_in = _mm_set_ps(s1,s0,s1,s0);
  __m128 l_stop_in = _mm_set_ps(l2,l2,l0,l0);

  __m128 s_start_out = _mm_set_ps(fabs(x)-hw,fabs(x)   ,fabs(x)-hw,fabs(x));
  __m128 s_stop_out  = _mm_set_ps(fabs(x)   ,fabs(x)+hw,fabs(x)   ,fabs(x)+hw);

  // |x| > hw
  __m128 mask = _mm_cmpgt_ps(_mm_fabs_ps(_mm_set1_ps(x)),_mm_set1_ps(hw));

  __m128 s_start = _mm_sel_ps(s_start_in,s_start_out,mask);
  __m128 s_stop  = _mm_sel_ps(s_stop_in,s_stop_out  ,mask);

  __m128 l_start_out = _mm_set_ps(fabs(y)-hh, fabs(y)-hh, fabs(y)   , fabs(y));
  __m128 l_stop_out  = _mm_set_ps(fabs(y)   , fabs(y)   , fabs(y)+hh, fabs(y)+hh);

  // Okay
  __m128 mask1 = _mm_cmpgt_ps(_mm_fabs_ps(_mm_set1_ps(y)),_mm_set1_ps(hh));

  __m128 l_start = _mm_sel_ps(l_start_in,l_start_out,mask1);
  __m128 l_stop  = _mm_sel_ps(l_stop_in, l_stop_out, mask1);


  // start (abs(x)-a,abs(x),abs(x)-a,abs(x))
  // stop  (abs(x),abs(x)+a,abs(x),abs(x)+a)
  __m128 s = _mm_sub_ps(s_stop,s_start); // Always positive
  __m128 l = _mm_sub_ps(l_stop,l_start);

  const __m128 vec_s = _mm_fabs_ps(s);
  const __m128 vec_l = _mm_fabs_ps(l);

  const __m128 cargz = _mm_set1_ps(cos(-k*z));
  const __m128 sargz = _mm_set1_ps(sin(-k*z));

  const __m128 vec_l_2 = _mm_mul_ps(l,_m_half_ps);

  const __m128 vec_s_2 = _mm_mul_ps(s,_m_half_ps);

  const __m128 z2     = _mm_set1_ps(SQUARE(z));

  __m128 real, imag;

  __m128 intWreal = _mm_setzero_ps();
  __m128 intWimag = _mm_setzero_ps();

  __m128 rcp_denom1 = _mm_rcp_ps(_mm_mul_ps(_m_2pi_ps,_mm_set1_ps(k)));


  // u-integral, s-integral
  for (size_t iu = 0 ; iu < nUs ; iu++) {

    __m128 us1       = _mm_load1_ps(&us[iu]);
    __m128 uweights1 = _mm_load1_ps(&uweights[iu]);

    // [0 ; |x|+hw] , [0 ; hw - |x|];
    __m128 ls  = _mm_add_ps(
                   _mm_mul_ps(
                     _mm_mul_ps(
                       _m_half_ps,
                       _mm_sub_ps(s_stop,s_start)),
                     us1),
                   _mm_mul_ps(
                     _m_half_ps,
                     _mm_add_ps(s_stop,s_start)));

    __m128 ls2 = _mm_square_ps(ls);

    __m128 argw = _mm_mul_ps(
                    _mm_set1_ps(-k),
                    _mm_sqrt_ps(
                      _mm_add_ps(
                        _mm_add_ps(
                          ls2,
                          z2),
                        _mm_square_ps(l_stop_in))));

    __m128 cargw, sargw;

    _mm_sin_cos_ps(argw, &sargw, &cargw);

    __m128 rcp_denom = _mm_rcp_ps(_mm_add_ps(ls2,_mm_square_ps(l_stop_in)));

    real = _mm_mul_ps(
             _mm_mul_ps(
               uweights1,
               _mm_sub_ps(
                 cargw,
                 cargz)),
             rcp_denom);

    imag = _mm_mul_ps(
             _mm_mul_ps(
               uweights1,
               _mm_sub_ps(
                 sargw,
                 sargz)),
             rcp_denom);
    intWreal = _mm_add_ps(intWreal, real);
    intWimag = _mm_add_ps(intWimag, imag);
  }

  // Divide by denominator
  intWreal = _mm_mul_ps(
               intWreal,
               _mm_mul_ps(
                 _mm_mul_ps(
                   _mm_sel_ps(vec_s_2,vec_s,mask), // Half integral
                   l_stop_in),                     // Stop value
                 rcp_denom1));
  intWimag = _mm_mul_ps(
               intWimag,
               _mm_mul_ps(
                 _mm_mul_ps(
                   _mm_sel_ps(vec_s_2,vec_s,mask),
                   l_stop_in),                     // Insert sign
                 rcp_denom1));

  // Filter (consider recover sign earlier)
  __m128 intHreal = _mm_setzero_ps();
  __m128 intHimag = _mm_setzero_ps();

  for(size_t iv = 0 ; iv < nVs ; iv++) {

    __m128 vs1       = _mm_load1_ps(&vs[iv]);
    __m128 vweights1 = _mm_load1_ps(&vweights[iv]);

    __m128 ss  = _mm_add_ps(
                   _mm_mul_ps(
                     _mm_mul_ps(
                       _m_half_ps,
                       _mm_sub_ps(
                         l_stop,
                         l_start)),
                     vs1),
                   _mm_mul_ps(
                     _m_half_ps,
                     _mm_add_ps(
                       l_start,
                       l_stop)));
    __m128 ss2 = _mm_square_ps(ss);

    __m128 argh = _mm_mul_ps(
                    _mm_set1_ps(-k),
                    _mm_sqrt_ps(
                      _mm_add_ps(
                        _mm_add_ps(
                          ss2,
                          z2),
                        _mm_square_ps(s_stop_in))));

    __m128 cargh, sargh;

    _mm_sin_cos_ps(argh, &sargh, &cargh);

    __m128 rcp_denom = _mm_rcp_ps(_mm_add_ps(ss2,_mm_square_ps(s_stop_in)));

    real = _mm_mul_ps(
             _mm_mul_ps(
               vweights1,
               _mm_sub_ps(
                 cargh,
                 cargz)),
             rcp_denom);
    imag = _mm_mul_ps(
             _mm_mul_ps(
               vweights1,
               _mm_sub_ps(
                 sargh,
                 sargz)),
             rcp_denom);
    intHreal = _mm_add_ps(intHreal, real);
    intHimag = _mm_add_ps(intHimag, imag);
  }

  // Divide by denominator
  intHreal = _mm_mul_ps(
               intHreal,
               _mm_mul_ps(
                 _mm_mul_ps(
                   _mm_sel_ps(vec_l_2,vec_l,mask1), // Half integral
                   s_stop_in),                      // Inserts sign
                 rcp_denom1));
  intHimag = _mm_mul_ps(
               intHimag,
               _mm_mul_ps(
                 _mm_mul_ps(
                   _mm_sel_ps(vec_l_2,vec_l,mask1), // Half integral
                   s_stop_in),         // Inserts sign
                 rcp_denom1));

  intHreal = _mm_add_ps(intHreal,intWreal);
  intHimag = _mm_add_ps(intHimag,intWimag);

  // Multiply by -i
  __m128 tmp = intHreal;
  intHreal = intHimag;
  intHimag = _mm_neg_ps(tmp);

  // Horizontal sum
  _mm_store_ss(&(reinterpret_cast<float(&)[2]>(retval)[0]),
               _mm_dp_ps(
                 _m_one_ps,
                 intHreal,0xF1));
  _mm_store_ss(&(reinterpret_cast<float(&)[2]>(retval)[1]),
               _mm_dp_ps(
                 _m_one_ps,
                 intHimag,0xF1));

  return retval;
}

template <>
std::complex<float>
CalcHzFast(const fnm::element_t<float>& element,
           const sps::point_t<float>& projection,
           const float& k,
           const float* us,
           const float* uweights,
           const size_t nUs,
           const float* vs,
           const float* vweights,
           const size_t nVs)
{

  std::complex<float> retval = std::complex<float>(0.0f,0.0f);

  // Temporaries
  const float x  = projection[0];
  const float y  = projection[1];
  const float z  = projection[2];
  const float hw = element.hw;
  const float hh = element.hh;

  const __m128 v_mk = _mm_set1_ps(-k);
  const __m128 v_z  = _mm_set1_ps(z);

  __m128 s_stop_in;
  __m128 l_stop_in;
  __m128 s_start_out;

  __m128 s_stop_out;
  __m128 l_start_out;
  __m128 l_stop_out;

  const __m128 v_absx = _mm_set1_ps(fabs(x));
  const __m128 v_absy = _mm_set1_ps(fabs(y));
  const __m128 v_hw   = _mm_set1_ps(hw);
  const __m128 v_hh   = _mm_set1_ps(hh);

#if 1
  const float absX = fabs(projection[0]);
  const float absY = fabs(projection[1]);

  float s0 = absX + hw; // [0;|x|+hw] -> [0;hw+|x|]  -> [|x|;|x|+hw]
  float s1 = hw - absX; // [0;hw-|x|] -> -[hw-|x|;0] -> [|x|-hw;|x|]
  float l0 = absY + hh;
  float l2 = hh - absY;

  // Negative (maybe)
  // hw + (-1,1,-1,1)*absX
  s_stop_in = _mm_set_ps(s1,s0,s1,s0);
  // hh + (-1,-1,1,1)*absX
  l_stop_in = _mm_set_ps(l2,l2,l0,l0);

  float s1_o = absX+hw;
  float s0_o = absX-hw;
  float l0_o = absY-hh;
  float l2_o = absY+hh;

  // absX + (-1,0,-1,0)*hw
  s_start_out = _mm_set_ps(s0_o, absX, s0_o, absX);
  // absX + (0,1,0,1)*hw
  s_stop_out  = _mm_set_ps(absX, s1_o, absX, s1_o);

  // absY + (-1,-1,0,0)*hh
  l_start_out = _mm_set_ps(l0_o, l0_o, absY, absY);
  // absY + (0,0,1,1)*hh
  l_stop_out  = _mm_set_ps(absY, absY, l2_o, l2_o);
#else

  const __m128 v_pmpm = _mm_set_ps(1.0f,-1.0f,1.0f,-1.0f);
  const __m128 v_ppmm = _mm_set_ps(1.0f,1.0f,-1.0f,-1.0f);
  const __m128 v_p0p0 = _mm_set_ps(1.0f,0.0f,1.0f,0.0f);
  const __m128 v_pp00 = _mm_set_ps(1.0f,1.0f,0.0f,0.0f);
  const __m128 v_absx = _mm_set1_ps(fabs(x));
  const __m128 v_absy = _mm_set1_ps(fabs(y));

  s_stop_in   = _mm_sub_ps(v_hw,_mm_mul_ps(v_pmpm,v_absx));
  l_stop_in   = _mm_sub_ps(v_hh,_mm_mul_ps(v_ppmm,v_absy));

  s_start_out = _mm_sub_ps(v_absx,_mm_mul_ps(v_p0p0,v_hw));
  s_stop_out  = _mm_add_ps(v_absx,_mm_mul_ps(_mm_sub_ps(_m_one_ps,v_p0p0),v_hw));

  l_start_out = _mm_sub_ps(v_absy,_mm_mul_ps(v_pp00,v_hh));
  l_stop_out  = _mm_add_ps(v_absy, _mm_mul_ps(_mm_sub_ps(_m_one_ps,v_pp00),v_hh));

#endif
  // Integral (stop, inside)
  __m128 start_in = _mm_setzero_ps();

  // |x| > hw
  __m128 mask = _mm_cmpgt_ps(v_absx,v_hw);

  __m128 s_start = _mm_sel_ps(start_in, s_start_out, mask);
  __m128 s_stop  = _mm_sel_ps(s_stop_in, s_stop_out, mask);

  // |y| > hh
  __m128 mask1   = _mm_cmpgt_ps(v_absy,v_hh);

  __m128 l_start = _mm_sel_ps(start_in, l_start_out, mask1);
  __m128 l_stop  = _mm_sel_ps(l_stop_in, l_stop_out, mask1);

  __m128 s = _mm_sub_ps(s_stop, s_start); // Always positive
  __m128 l = _mm_sub_ps(l_stop, l_start); // Always positive

  __m128 cargz, sargz;
#if ACCURATE_TRIGONOMETRICS
  sargz = _mm_set1_ps(sin(-z*k));
  cargz = _mm_set1_ps(cos(-z*k));
#else
  _mm_sin_cos_ps(_mm_mul_ps(v_mk,v_z),&sargz,&cargz);
#endif


  const __m128 v_z2     = _mm_square_ps(v_z);
  const __m128 v_l2 = _mm_square_ps(l_stop_in);
  const __m128 v_s2 = _mm_square_ps(s_stop_in);

  __m128 rcp_denom1 = _mm_rcp_ps(_mm_mul_ps(_m_2pi_ps,_mm_neg_ps(v_mk)));


  // u-integral, s-integral
  const __m128 s_offset = _mm_mul_ps(_m_half_ps,_mm_add_ps(s_stop,s_start));
  const __m128 s_scale  = _mm_mul_ps(_m_half_ps,_mm_sub_ps(s_stop,s_start));

  const __m128 l_offset = _mm_mul_ps(_m_half_ps,_mm_add_ps(l_start,l_stop));
  const __m128 l_scale  =_mm_mul_ps(_m_half_ps,_mm_sub_ps(l_stop,l_start));

  __m128 intWreal = _mm_setzero_ps();
  __m128 intWimag = _mm_setzero_ps();

  __m128 real, imag;

  for (size_t iu = 0 ; iu < nUs ; iu++) {

    __m128 us1       = _mm_load1_ps(&us[iu]);
    __m128 uweights1 = _mm_load1_ps(&uweights[iu]);

    // [0 ; |x|+hw] , [0 ; hw - |x|];
    __m128 ls  = _mm_add_ps(_mm_mul_ps(s_scale,us1),s_offset);

    __m128 ls2 = _mm_square_ps(ls);

    __m128 argw = _mm_mul_ps(
                    v_mk,
                    _mm_sqrt_ps(
                      _mm_add_ps(
                        _mm_add_ps(
                          ls2,
                          v_z2),
                        v_l2)));

    __m128 cargw, sargw;

#if ACCURATE_TRIGONOMETRICS
    ALIGN16_BEGIN float in[4] ALIGN16_END;
    ALIGN16_BEGIN float out[4] ALIGN16_END;
    _mm_store_ps(in,argw);
    out[0] = sin(in[0]);
    out[1] = sin(in[1]);
    out[2] = sin(in[2]);
    sargw = _mm_load_ps(out);
    out[0] = cos(in[0]);
    out[1] = cos(in[1]);
    out[2] = cos(in[2]);
    cargw = _mm_load_ps(out);
#else
    _mm_sin_cos_ps(argw, &sargw, &cargw);
#endif

    __m128 rcp_denom = _mm_rcp_ps(_mm_add_ps(ls2,v_l2));

    real = _mm_mul_ps(
             _mm_mul_ps(
               uweights1,
               _mm_sub_ps(
                 cargw,
                 cargz)),
             rcp_denom);

    imag = _mm_mul_ps(
             _mm_mul_ps(
               uweights1,
               _mm_sub_ps(
                 sargw,
                 sargz)),
             rcp_denom);
    intWreal = _mm_add_ps(intWreal, real);
    intWimag = _mm_add_ps(intWimag, imag);
  }

  __m128 intScale = _mm_mul_ps(
                      _mm_mul_ps(
                        _mm_sel_ps(s_scale,_mm_mul_ps(s,_m_half_ps),mask), // Half integral
                        l_stop_in),                 // Stop value
                      rcp_denom1);

  // Divide by denominator
  intWreal = _mm_mul_ps(intWreal,intScale);
  intWimag = _mm_mul_ps(intWimag,intScale);

  // Filter (consider recover sign earlier)
  __m128 intHreal = _mm_setzero_ps();
  __m128 intHimag = _mm_setzero_ps();

  for(size_t iv = 0 ; iv < nVs ; iv++) {

    __m128 vs1       = _mm_load1_ps(&vs[iv]);
    __m128 vweights1 = _mm_load1_ps(&vweights[iv]);

    __m128 ss  = _mm_add_ps(_mm_mul_ps(l_scale,vs1),l_offset);

    __m128 ss2 = _mm_square_ps(ss);

    __m128 argh = _mm_mul_ps(v_mk,
                             _mm_sqrt_ps(
                               _mm_add_ps(
                                 _mm_add_ps(
                                   ss2,
                                   v_z2),
                                 v_s2)));

    __m128 cargh, sargh;

#if ACCURATE_TRIGONOMETRICS
    ALIGN16_BEGIN float in[4] ALIGN16_END;
    ALIGN16_BEGIN float out[4] ALIGN16_END;
    _mm_store_ps(in,argh);
    out[0] = sin(in[0]);
    out[1] = sin(in[1]);
    out[2] = sin(in[2]);
    sargh = _mm_load_ps(out);
    out[0] = cos(in[0]);
    out[1] = cos(in[1]);
    out[2] = cos(in[2]);
    cargh = _mm_load_ps(out);
#else
    _mm_sin_cos_ps(argh, &sargh, &cargh);
#endif

    __m128 rcp_denom = _mm_rcp_ps(_mm_add_ps(ss2,v_s2));

    real = _mm_mul_ps(
             _mm_mul_ps(
               vweights1,
               _mm_sub_ps(
                 cargh,
                 cargz)),
             rcp_denom);
    imag = _mm_mul_ps(
             _mm_mul_ps(
               vweights1,
               _mm_sub_ps(
                 sargh,
                 sargz)),
             rcp_denom);
    intHreal = _mm_add_ps(intHreal, real);
    intHimag = _mm_add_ps(intHimag, imag);
  }

  intScale = _mm_mul_ps(
               _mm_mul_ps(
                 _mm_sel_ps(l_scale,_mm_mul_ps(_m_half_ps,l),mask1), // Half integral
                 s_stop_in),                  // Inserts sign
               rcp_denom1);

  // Divide by denominator
  intHreal = _mm_mul_ps(intHreal,intScale);
  intHimag = _mm_mul_ps(intHimag,intScale);

  intHreal = _mm_add_ps(intHreal,intWreal);
  intHimag = _mm_add_ps(intHimag,intWimag);

  // Multiply by -i
  __m128 tmp = intHreal;
  intHreal = intHimag;
  intHimag = _mm_neg_ps(tmp);

  // Horizontal sums
  _mm_store_ss(
    &(reinterpret_cast<float(&)[2]>(retval)[0]),
    _mm_dp_ps(
      _m_one_ps,
      intHreal,
      0xF1));
  _mm_store_ss(&(reinterpret_cast<float(&)[2]>(retval)[1]),
               _mm_dp_ps(
                 _m_one_ps,
                 intHimag,
                 0xF1));

  return retval;
}

// Only for region |x| > a, |y| > b, forgot the V-part
template <>
std::complex<float>
CalcHzFast4(const fnm::element_t<float>& element,
            const sps::point_t<float>& projection,
            const float& k,
            const float* us,
            const float* uweights,
            const size_t nUs,
            const float* vs,
            const float* vweights,
            const size_t nVs)
{

  SPS_UNREFERENCED_PARAMETER(nVs);
  SPS_UNREFERENCED_PARAMETER(vweights);
  SPS_UNREFERENCED_PARAMETER(vs);

  const float z       = projection[2];

  const __m128 v_pmpm = _mm_set_ps(1.0f,-1.0f,1.0f,-1.0f);
  const __m128 v_xxyy  = _mm_set_ps(projection[0],projection[0],projection[1],projection[1]);
  const __m128 v_aabb  = _mm_set_ps(element.hw,element.hw,element.hh,element.hh);

  const __m128 vec_mk = _mm_set1_ps(-k);
  const __m128 vec_z  = _mm_set1_ps(z);

  const __m128 rcp_denom1 = _mm_rcp_ps(_mm_mul_ps(_m_2pi_ps,_mm_set1_ps(k)));
  const __m128 vec_z2 = _mm_square_ps(vec_z);

  __m128 cargz = _mm_setzero_ps();
  __m128 sargz = _mm_setzero_ps();
  _mm_sin_cos_ps(_mm_mul_ps(vec_mk,vec_z),&sargz,&cargz);

  __m128 v_bbaa  = _mm_permute_ps(v_aabb,   0xB1);
  __m128 v_yyxx  = _mm_permute_ps(v_xxyy,   0xB1);

  __m128 v_stop  = _mm_add_ps(v_xxyy,v_aabb);
  __m128 v_start = _mm_sub_ps(v_xxyy,v_aabb);

  __m128 v_scale  = _mm_mul_ps(_m_half_ps,_mm_sub_ps(v_stop,v_start));
  __m128 v_offset = _mm_mul_ps(_m_half_ps,_mm_add_ps(v_stop,v_start));

  __m128 intReal = _mm_setzero_ps();
  __m128 intImag = _mm_setzero_ps();

  __m128 real, imag;

  __m128 v_adj  = _mm_add_ps(v_yyxx,_mm_mul_ps(v_pmpm,v_bbaa));
  __m128 v_adj2 = _mm_square_ps(v_adj);

  // Both integrals (assumes nUs = nUv)
  for (size_t iu = 0 ; iu < nUs ; iu++) {
    __m128 v_s = _mm_load1_ps(&us[iu]);
    __m128 v_w  = _mm_load1_ps(&uweights[iu]);

    __m128 v_xy = _mm_add_ps(_mm_mul_ps(v_scale,v_s),v_offset);

    __m128 v_xy2 = _mm_square_ps(v_xy);

    __m128 argw = _mm_mul_ps(
                    vec_mk,
                    _mm_sqrt_ps(
                      _mm_add_ps(
                        _mm_add_ps(
                          v_xy2,
                          vec_z2),
                        v_adj2)));

    __m128 cargw, sargw;

    _mm_sin_cos_ps(argw, &sargw, &cargw);

    __m128 rcp_denom = _mm_rcp_ps(_mm_add_ps(v_xy2,v_adj2));

    real = _mm_mul_ps(
             _mm_mul_ps(
               v_w,
               _mm_sub_ps(
                 cargw,
                 cargz)),
             rcp_denom);

    imag = _mm_mul_ps(
             _mm_mul_ps(
               v_w,
               _mm_sub_ps(
                 sargw,
                 sargz)),
             rcp_denom);
    intReal = _mm_add_ps(intReal, real);
    intImag = _mm_add_ps(intImag, imag);
  }

  __m128 v_int_scale = _mm_mul_ps(_mm_mul_ps(
                                    _mm_mul_ps(
                                      v_adj,
                                      v_scale),
                                    rcp_denom1),/*_mm_permute_ps(*/v_pmpm/*,0xB1)*/);

  // Divide by denominator
  intReal = _mm_mul_ps(intReal,v_int_scale);
  intImag = _mm_mul_ps(intImag,v_int_scale);

  // Multiply by -i
  __m128 tmp = intReal;
  intReal = intImag;
  intImag = _mm_neg_ps(tmp);

  std::complex<float> retval = std::complex<float>(0.0f,0.0f);

  // Horizontal sums
  _mm_store_ss(
    &(reinterpret_cast<float(&)[2]>(retval)[0]),
    _mm_dp_ps(
      _m_one_ps,
      intReal,
      0xF1));
  _mm_store_ss(&(reinterpret_cast<float(&)[2]>(retval)[1]),
               _mm_dp_ps(
                 _m_one_ps,
                 intImag,
                 0xF1));

  return retval;

}

template <>
std::complex<double>
CalcHzFast(const fnm::element_t<double>& element,
           const sps::point_t<double>& projection,
           const double& k,
           const double* us,
           const double* uweights,
           const size_t nUs,
           const double* vs,
           const double* vweights,
           const size_t nVs)
{

  std::complex<double> retval = std::complex<double>(0.0f,0.0f);

  const double absX = fabs(projection[0]);
  const double absY = fabs(projection[1]);
  const double z    = projection[2];

  const double hw = element.hw;
  const double hh = element.hh;

  double s0 = absX + hw; // [0;|x|+hw] -> [0;hw+|x|]  -> [|x|;|x|+hw]
  double s1 = hw - absX; // [0;hw-|x|] -> -[hw-|x|;0] -> [|x|-hw;|x|]
  double l0 = absY + hh;
  double l2 = hh - absY;

  // Integral (stop, inside)
  __m128d start_in = _mm_setzero_pd();

  const __m128d vec_mk = _mm_set1_pd(-k);
  const __m128d vec_z  = _mm_set1_pd(z);

  // First half

  // Negative
  __m128d s_stop_in = _mm_set_pd(s1,s0);
  __m128d l_stop_in = _mm_set_pd(l2,l2);

  __m128d s_start_out = _mm_set_pd(absX-hw,absX   );
  __m128d s_stop_out  = _mm_set_pd(absX   ,absX+hw);

  __m128d l_start_out = _mm_set_pd(absY-hh, absY-hh);
  __m128d l_stop_out  = _mm_set_pd(absY   , absY   );

  // |x| > hw
  __m128d mask = _mm_cmpgt_pd(_mm_fabs_pd(_mm_set1_pd(projection[0])),_mm_set1_pd(hw));

  __m128d s_start = _mm_sel_pd(start_in,s_start_out,mask);
  __m128d s_stop  = _mm_sel_pd(s_stop_in,s_stop_out  ,mask);

  // Okay
  __m128d mask1 = _mm_cmpgt_pd(_mm_fabs_pd(_mm_set1_pd(projection[1])),_mm_set1_pd(hh));

  __m128d l_start = _mm_sel_pd(start_in,l_start_out,mask1);
  __m128d l_stop  = _mm_sel_pd(l_stop_in, l_stop_out, mask1);

  __m128d s = _mm_sub_pd(s_stop,s_start); // Always positive
  __m128d l = _mm_sub_pd(l_stop,l_start); // Always positive

  __m128d cargz, sargz;
  _mm_sin_cos_pd(_mm_mul_pd(vec_mk,vec_z),&sargz,&cargz);

  const __m128d z2     = _mm_square_pd(vec_z);
  __m128d vec_l2 = _mm_square_pd(l_stop_in);
  __m128d vec_s2 = _mm_square_pd(s_stop_in);

  __m128d rcp_denom1 = _mm_rcp_pd(_mm_mul_pd(_m_2pi_pd,_mm_set1_pd(k)));


  // u-integral, s-integral
  __m128d s_offset = _mm_mul_pd(_m_half_pd,_mm_add_pd(s_stop,s_start));
  __m128d s_scale  = _mm_mul_pd(_m_half_pd,_mm_sub_pd(s_stop,s_start));

  __m128d l_offset = _mm_mul_pd(_m_half_pd,_mm_add_pd(l_start,l_stop));
  __m128d l_scale  =_mm_mul_pd(_m_half_pd,_mm_sub_pd(l_stop,l_start));

  __m128d intWreal = _mm_setzero_pd();
  __m128d intWimag = _mm_setzero_pd();

  __m128d real, imag;

  for (size_t iu = 0 ; iu < nUs ; iu++) {

    __m128d us1       = _mm_load1_pd(&us[iu]);
    __m128d uweights1 = _mm_load1_pd(&uweights[iu]);

    // [0 ; |x|+hw] , [0 ; hw - |x|];
    __m128d ls  = _mm_add_pd(_mm_mul_pd(s_scale,us1),s_offset);

    __m128d ls2 = _mm_square_pd(ls);

    __m128d argw = _mm_mul_pd(
                     vec_mk,
                     _mm_sqrt_pd(
                       _mm_add_pd(
                         _mm_add_pd(
                           ls2,
                           z2),
                         vec_l2)));

    __m128d cargw, sargw;

    _mm_sin_cos_pd(argw, &sargw, &cargw);

    __m128d rcp_denom = _mm_rcp_pd(_mm_add_pd(ls2,vec_l2));

    real = _mm_mul_pd(
             _mm_mul_pd(
               uweights1,
               _mm_sub_pd(
                 cargw,
                 cargz)),
             rcp_denom);

    imag = _mm_mul_pd(
             _mm_mul_pd(
               uweights1,
               _mm_sub_pd(
                 sargw,
                 sargz)),
             rcp_denom);
    intWreal = _mm_add_pd(intWreal, real);
    intWimag = _mm_add_pd(intWimag, imag);
  }

  __m128d intScale = _mm_mul_pd(
                       _mm_mul_pd(
                         _mm_sel_pd(s_scale,s,mask), // Half integral
                         l_stop_in),                     // Stop value
                       rcp_denom1);

  // Divide by denominator
  intWreal = _mm_mul_pd(intWreal,intScale);
  intWimag = _mm_mul_pd(intWimag,intScale);

  // Filter (consider recover sign earlier)
  __m128d intHreal = _mm_setzero_pd();
  __m128d intHimag = _mm_setzero_pd();

  for(size_t iv = 0 ; iv < nVs ; iv++) {

    __m128d vs1       = _mm_load1_pd(&vs[iv]);
    __m128d vweights1 = _mm_load1_pd(&vweights[iv]);

    __m128d ss  = _mm_add_pd(_mm_mul_pd(l_scale,vs1),l_offset);

    __m128d ss2 = _mm_square_pd(ss);

    __m128d argh = _mm_mul_pd(vec_mk,
                              _mm_sqrt_pd(
                                _mm_add_pd(
                                  _mm_add_pd(
                                    ss2,
                                    z2),
                                  vec_s2)));

    __m128d cargh, sargh;

    _mm_sin_cos_pd(argh, &sargh, &cargh);

    __m128d rcp_denom = _mm_rcp_pd(_mm_add_pd(ss2,vec_s2));

    real = _mm_mul_pd(
             _mm_mul_pd(
               vweights1,
               _mm_sub_pd(
                 cargh,
                 cargz)),
             rcp_denom);
    imag = _mm_mul_pd(
             _mm_mul_pd(
               vweights1,
               _mm_sub_pd(
                 sargh,
                 sargz)),
             rcp_denom);
    intHreal = _mm_add_pd(intHreal, real);
    intHimag = _mm_add_pd(intHimag, imag);
  }

  intScale = _mm_mul_pd(
               _mm_mul_pd(
                 _mm_sel_pd(l_scale,l,mask1), // Half integral
                 s_stop_in),                  // Inserts sign
               rcp_denom1);
  // Divide by denominator
  intHreal = _mm_mul_pd(intHreal,intScale);
  intHimag = _mm_mul_pd(intHimag,intScale);

  intHreal = _mm_add_pd(intHreal,intWreal);
  intHimag = _mm_add_pd(intHimag,intWimag);

  // Multiply by -i
  __m128d tmp = intHreal;
  intHreal = intHimag;
  intHimag = _mm_neg_pd(tmp);

  // Horizontal sum
  ALIGN16_BEGIN double tmpd[2] ALIGN16_END;
  _mm_store_pd(tmpd,_mm_dp_pd(_m_one_pd,intHreal,0xF1));
  retval.real(tmpd[0]);

  _mm_store_pd(tmpd,_mm_dp_pd(_m_one_pd,intHimag,0xF1));
  retval.imag(tmpd[0]);

  // Second half
  // Negative
  s_stop_in = _mm_set_pd(s1,s0);
  l_stop_in = _mm_set_pd(l0,l0);

  s_start_out = _mm_set_pd(absX-hw,absX);
  s_stop_out  = _mm_set_pd(absX   ,absX+hw);

  l_start_out = _mm_set_pd(absY   , absY);
  l_stop_out  = _mm_set_pd(absY+hh, absY+hh);

  // |x| > hw
  mask = _mm_cmpgt_pd(_mm_fabs_pd(_mm_set1_pd(projection[0])),_mm_set1_pd(hw));

  s_start = _mm_sel_pd(start_in,s_start_out,mask);
  s_stop  = _mm_sel_pd(s_stop_in,s_stop_out  ,mask);

  // Okay
  mask1 = _mm_cmpgt_pd(_mm_fabs_pd(_mm_set1_pd(projection[1])),_mm_set1_pd(hh));

  l_start = _mm_sel_pd(start_in,l_start_out,mask1);
  l_stop  = _mm_sel_pd(l_stop_in, l_stop_out, mask1);

  s = _mm_sub_pd(s_stop,s_start); // Always positive
  l = _mm_sub_pd(l_stop,l_start); // Always positive

  _mm_sin_cos_pd(_mm_mul_pd(vec_mk,vec_z),&sargz,&cargz);

  vec_l2 = _mm_square_pd(l_stop_in);
  vec_s2 = _mm_square_pd(s_stop_in);

  rcp_denom1 = _mm_rcp_pd(_mm_mul_pd(_m_2pi_pd,_mm_set1_pd(k)));


  // u-integral, s-integral
  s_offset = _mm_mul_pd(_m_half_pd,_mm_add_pd(s_stop,s_start));
  s_scale  = _mm_mul_pd(_m_half_pd,_mm_sub_pd(s_stop,s_start));

  l_offset = _mm_mul_pd(_m_half_pd,_mm_add_pd(l_start,l_stop));
  l_scale  =_mm_mul_pd(_m_half_pd,_mm_sub_pd(l_stop,l_start));

  intWreal = _mm_setzero_pd();
  intWimag = _mm_setzero_pd();

  for (size_t iu = 0 ; iu < nUs ; iu++) {

    __m128d us1       = _mm_load1_pd(&us[iu]);
    __m128d uweights1 = _mm_load1_pd(&uweights[iu]);

    // [0 ; |x|+hw] , [0 ; hw - |x|];
    __m128d ls  = _mm_add_pd(_mm_mul_pd(s_scale,us1),s_offset);

    __m128d ls2 = _mm_square_pd(ls);

    __m128d argw = _mm_mul_pd(
                     vec_mk,
                     _mm_sqrt_pd(
                       _mm_add_pd(
                         _mm_add_pd(
                           ls2,
                           z2),
                         vec_l2)));

    __m128d cargw, sargw;

    _mm_sin_cos_pd(argw, &sargw, &cargw);

    __m128d rcp_denom = _mm_rcp_pd(_mm_add_pd(ls2,vec_l2));

    real = _mm_mul_pd(
             _mm_mul_pd(
               uweights1,
               _mm_sub_pd(
                 cargw,
                 cargz)),
             rcp_denom);

    imag = _mm_mul_pd(
             _mm_mul_pd(
               uweights1,
               _mm_sub_pd(
                 sargw,
                 sargz)),
             rcp_denom);
    intWreal = _mm_add_pd(intWreal, real);
    intWimag = _mm_add_pd(intWimag, imag);
  }

  intScale = _mm_mul_pd(
               _mm_mul_pd(
                 _mm_sel_pd(s_scale,s,mask), // Half integral
                 l_stop_in),                     // Stop value
               rcp_denom1);

  // Divide by denominator
  intWreal = _mm_mul_pd(intWreal,intScale);
  intWimag = _mm_mul_pd(intWimag,intScale);

  // Filter (consider recover sign earlier)
  intHreal = _mm_setzero_pd();
  intHimag = _mm_setzero_pd();

  for(size_t iv = 0 ; iv < nVs ; iv++) {

    __m128d vs1       = _mm_load1_pd(&vs[iv]);
    __m128d vweights1 = _mm_load1_pd(&vweights[iv]);

    __m128d ss  = _mm_add_pd(_mm_mul_pd(l_scale,vs1),l_offset);

    __m128d ss2 = _mm_square_pd(ss);

    __m128d argh = _mm_mul_pd(vec_mk,
                              _mm_sqrt_pd(
                                _mm_add_pd(
                                  _mm_add_pd(
                                    ss2,
                                    z2),
                                  vec_s2)));

    __m128d cargh, sargh;

    _mm_sin_cos_pd(argh, &sargh, &cargh);

    __m128d rcp_denom = _mm_rcp_pd(_mm_add_pd(ss2,vec_s2));

    real = _mm_mul_pd(
             _mm_mul_pd(
               vweights1,
               _mm_sub_pd(
                 cargh,
                 cargz)),
             rcp_denom);
    imag = _mm_mul_pd(
             _mm_mul_pd(
               vweights1,
               _mm_sub_pd(
                 sargh,
                 sargz)),
             rcp_denom);
    intHreal = _mm_add_pd(intHreal, real);
    intHimag = _mm_add_pd(intHimag, imag);
  }

  intScale = _mm_mul_pd(
               _mm_mul_pd(
                 _mm_sel_pd(l_scale,l,mask1), // Half integral
                 s_stop_in),                  // Inserts sign
               rcp_denom1);
  // Divide by denominator
  intHreal = _mm_mul_pd(intHreal,intScale);
  intHimag = _mm_mul_pd(intHimag,intScale);

  intHreal = _mm_add_pd(intHreal,intWreal);
  intHimag = _mm_add_pd(intHimag,intWimag);

  // Multiply by -i
  tmp = intHreal;
  intHreal = intHimag;
  intHimag = _mm_neg_pd(tmp);

  // Horizontal sum
  _mm_store_pd(tmpd,_mm_dp_pd(_m_one_pd,intHreal,0xF1));
  retval.real(retval.real() + tmpd[0]);

  _mm_store_pd(tmpd,_mm_dp_pd(_m_one_pd,intHimag,0xF1));
  retval.imag(retval.imag() + tmpd[0]);

  return retval;
}

template <>
std::complex<double> CalcHzVecGL(const double& s,
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
  SPS_UNREFERENCED_PARAMETER(s);
  SPS_UNREFERENCED_PARAMETER(l);
  SPS_UNREFERENCED_PARAMETER(z);
  SPS_UNREFERENCED_PARAMETER(k);
  SPS_UNREFERENCED_PARAMETER(us);
  SPS_UNREFERENCED_PARAMETER(uweights);
  SPS_UNREFERENCED_PARAMETER(nUs);
  SPS_UNREFERENCED_PARAMETER(vs);
  SPS_UNREFERENCED_PARAMETER(vweights);
  SPS_UNREFERENCED_PARAMETER(nVs);
  return std::complex<double>();
}

template <>
std::complex<double> CalcHzAll(const fnm::element_t<double>& element,
                               const sps::point_t<double>& projection,
                               const double& k,
                               const double* us,
                               const double* uweights,
                               const size_t nUs,
                               const double* vs,
                               const double* vweights,
                               const size_t nVs)
{
  SPS_UNREFERENCED_PARAMETER(element);
  SPS_UNREFERENCED_PARAMETER(projection);
  SPS_UNREFERENCED_PARAMETER(k);
  SPS_UNREFERENCED_PARAMETER(us);
  SPS_UNREFERENCED_PARAMETER(uweights);
  SPS_UNREFERENCED_PARAMETER(nUs);
  SPS_UNREFERENCED_PARAMETER(vs);
  SPS_UNREFERENCED_PARAMETER(vweights);
  SPS_UNREFERENCED_PARAMETER(nVs);

  std::complex<double> retval;
  return retval;
}

template <>
std::complex<double>
CalcHzShort(const fnm::element_t<double>& element,
            const sps::point_t<double>& projection,
            const double& k,
            const double* us,
            const double* uweights,
            const size_t nUs,
            const double* vs,
            const double* vweights,
            const size_t nVs)
{

  SPS_UNREFERENCED_PARAMETER(element);
  SPS_UNREFERENCED_PARAMETER(projection);
  SPS_UNREFERENCED_PARAMETER(k);
  SPS_UNREFERENCED_PARAMETER(us);
  SPS_UNREFERENCED_PARAMETER(uweights);
  SPS_UNREFERENCED_PARAMETER(nUs);
  SPS_UNREFERENCED_PARAMETER(vs);
  SPS_UNREFERENCED_PARAMETER(vweights);
  SPS_UNREFERENCED_PARAMETER(nVs);
  std::complex<double> retval;
  return retval;
}


template <>
std::complex<double>
CalcHzFast4(const fnm::element_t<double>& element,
            const sps::point_t<double>& projection,
            const double& k,
            const double* us,
            const double* uweights,
            const size_t nUs,
            const double* vs,
            const double* vweights,
            const size_t nVs)
{
  SPS_UNREFERENCED_PARAMETER(element);
  SPS_UNREFERENCED_PARAMETER(projection);
  SPS_UNREFERENCED_PARAMETER(k);
  SPS_UNREFERENCED_PARAMETER(us);
  SPS_UNREFERENCED_PARAMETER(uweights);
  SPS_UNREFERENCED_PARAMETER(nUs);
  SPS_UNREFERENCED_PARAMETER(vs);
  SPS_UNREFERENCED_PARAMETER(vweights);
  SPS_UNREFERENCED_PARAMETER(nVs);
  return std::complex<double>(0.0,0.0);
}
/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */

