/**
 * @file   fnm_ps.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Sat Apr  1 18:50:20 2017
 *
 * @brief
 *
 *
 */

/*
 *  This file is part of SOFUS.
 *
 *  SOFUS is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  SOFUS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with SOFUS.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <fnm/config.h>
#include <sps/cenv.h>
#include <sps/smath.hpp>
#include <fnm/fnm_types.hpp>
#include <sps/trigintrin.h>

#include <complex>

template <>
std::complex<float>
inline CalcHzFast(const sps::element_rect_t<float> &__restrict element,
                  const sps::point_t<float> &__restrict projection,
                  const float &__restrict k,
                  const float* __restrict us,
                  const float* __restrict uweights,
                  const size_t nUs,
                  const float* __restrict vs,
                  const float* __restrict vweights,
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

#if 0
  // Look equally bad
  const float absX = fabs(projection[0]);
  const float absY = fabs(projection[1]);

  float s0 = absX + hw; // [0;|x|+hw] -> [0;hw+|x|]  -> [|x|;|x|+hw]
  float s1 = hw - absX; // [0;hw-|x|] -> -[hw-|x|;0] -> [|x|-hw;|x|]
  float l0 = absY + hh;
  float l2 = hh - absY;

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

  _mm_sin_cos_ps(_mm_mul_ps(v_mk,v_z),&sargz,&cargz);

  const __m128 v_z2 = _mm_square_ps(v_z);
  const __m128 v_l2 = _mm_square_ps(l_stop_in);
  const __m128 v_s2 = _mm_square_ps(s_stop_in);

  __m128 rcp_denom1 = _mm_rcp_ps(_mm_mul_ps(_m_2pi_ps,_mm_neg_ps(v_mk)));

  // u-integral, s-integral
  const __m128 s_offset = _mm_mul_ps(_m_half_ps,_mm_add_ps(s_stop,s_start));
  const __m128 s_scale  = _mm_mul_ps(_m_half_ps,_mm_sub_ps(s_stop,s_start));

  // v-integral, l-integral
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

    _mm_sin_cos_ps(argw, &sargw, &cargw);

    __m128 denom = _mm_add_ps(ls2,v_l2);
    __m128 rcp_denom = _mm_rcp_ps(denom);
    // If the integral has zero length, reciprocal of denominator is set to zero
    __m128 mask_denom = _mm_cmplt_ps(_mm_fabs_ps(s_scale),_m_eps_ps);
    rcp_denom = _mm_sel_ps(rcp_denom,_mm_setzero_ps(),mask_denom);

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
                        _mm_sel_ps(s_scale,_mm_mul_ps(s,_m_half_ps),mask),
                        l_stop_in),
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

    _mm_sin_cos_ps(argh, &sargh, &cargh);

    __m128 denom = _mm_add_ps(ss2,v_s2);
    __m128 rcp_denom = _mm_rcp_ps(denom);
    // If the integral has zero length, reciprocal of denominator is set to zero
    __m128 mask_denom = _mm_cmplt_ps(_mm_fabs_ps(l_scale),_m_eps_ps);
    rcp_denom = _mm_sel_ps(rcp_denom,_mm_setzero_ps(),mask_denom);

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
                 _mm_sel_ps(l_scale,_mm_mul_ps(_m_half_ps,l),mask1),
                 s_stop_in),
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
  retval.real(_mm_cvtss_f32(_mm_dp_ps(_m_one_ps,intHreal,0xF1)));
  retval.imag(_mm_cvtss_f32(_mm_dp_ps(_m_one_ps,intHimag,0xF1)));

  return retval;
}

template <>
std::complex<float>
inline CalcFastFourAny(const float& u,
                       const float& v,
                       const float& hw,
                       const float& hh,
                       const float& z,
                       const float& __restrict k,
                       const float* __restrict s,
                       const float* __restrict weights,
                       const size_t nS)
{
  std::complex<float> retval = std::complex<float>(0.0f,0.0f);

  const __m128 v_mk = _mm_set1_ps(-k);
  const __m128 v_z  = _mm_set1_ps(z);

  const __m128 v_z2 = _mm_square_ps(v_z);

  __m128 adjacent;

  __m128 v_uuvv = _mm_set_ps(v,v,u,u); // 3, 2, 1, 0
  __m128 v_wwhh = _mm_set_ps(hh,hh,hw,hw);

  __m128 low  = _mm_sub_ps(_mm_fabs_ps(v_uuvv),v_wwhh);
  __m128 high = _mm_add_ps(_mm_fabs_ps(v_uuvv),v_wwhh);

  __m128 scale  = _mm_mul_ps(_m_half_ps,_mm_sub_ps(high,low));
  __m128 offset = _mm_mul_ps(_m_half_ps,_mm_add_ps(high,low));

  const __m128 v_pmpm = _mm_set_ps(-1.0f,1.0f,-1.0f,1.0f);

  adjacent = _mm_add_ps(
               _mm_permute_ps(v_wwhh, 0x1B),
               _mm_mul_ps(
                 v_pmpm,
                 _mm_fabs_ps(
                   _mm_permute_ps(v_uuvv, 0x1B))));

  __m128 intReal = _mm_setzero_ps();
  __m128 intImag = _mm_setzero_ps();

  __m128 real, imag;

  __m128 adjacent2 = _mm_square_ps(adjacent);

  __m128 cargz, sargz;

  _mm_sin_cos_ps(_mm_mul_ps(v_mk,v_z),&sargz,&cargz);

  __m128 rcp_denom1 = _mm_rcp_ps(_mm_mul_ps(_m_2pi_ps,_mm_neg_ps(v_mk)));

  for(size_t iuv = 0 ; iuv < nS ; iuv++) {

    __m128 v_uvs   = _mm_load1_ps(&s[iuv]);
    __m128 v_weights = _mm_load1_ps(&weights[iuv]);

    __m128 ss  = _mm_add_ps(_mm_mul_ps(scale,v_uvs),offset);

    __m128 ss2 = _mm_square_ps(ss);

    __m128 argh = _mm_mul_ps(v_mk,
                             _mm_sqrt_ps(
                               _mm_add_ps(
                                 _mm_add_ps(
                                   ss2,
                                   v_z2),
                                 adjacent2)));

    __m128 cargh, sargh;

    _mm_sin_cos_ps(argh, &sargh, &cargh);

    __m128 denom = _mm_add_ps(ss2,adjacent2);
    __m128 rcp_denom = _mm_rcp_ps(denom);

    __m128 mask_denom = _mm_cmplt_ps(_mm_fabs_ps(scale),_m_eps_ps);
    rcp_denom = _mm_sel_ps(rcp_denom,_mm_setzero_ps(),mask_denom);

    real = _mm_mul_ps(
             _mm_mul_ps(
               v_weights,
               _mm_sub_ps(
                 cargh,
                 cargz)),
             rcp_denom);
    imag = _mm_mul_ps(
             _mm_mul_ps(
               v_weights,
               _mm_sub_ps(
                 sargh,
                 sargz)),
             rcp_denom);
    intReal = _mm_add_ps(intReal, real);
    intImag = _mm_add_ps(intImag, imag);
  }

  __m128 intScale = _mm_mul_ps(
                      _mm_mul_ps(
                        scale,
                        adjacent),
                      rcp_denom1);

  // Divide by denominator
  intReal = _mm_mul_ps(intReal,intScale);
  intImag = _mm_mul_ps(intImag,intScale);

  // Multiply by -i
  __m128 tmp = intReal;
  intReal = intImag;
  intImag = _mm_neg_ps(tmp);

  // Horizontal sums (two first integrals)
  retval.real(_mm_cvtss_f32(_mm_dp_ps(_m_one_ps,intReal,0x31)));
  retval.imag(_mm_cvtss_f32(_mm_dp_ps(_m_one_ps,intImag,0x31)));

  return retval;
}



// Integrate two samples at a time (WORKS). TODO: Compare with CalcFastFourAny()
template <>
std::complex<float>
inline CalcFastFourAny2(const float& u,
                        const float& v,
                        const float& hw,
                        const float& hh,
                        const float& z,
                        const float& __restrict k,
                        const float* __restrict s,
                        const float* __restrict weights,
                        const size_t nS)
{
  std::complex<float> retval = std::complex<float>(0.0f,0.0f);

  const __m128 v_mk = _mm_set1_ps(-k);
  const __m128 v_z  = _mm_set1_ps(z);

  const __m128 v_z2 = _mm_square_ps(v_z);

  __m128 adjacent;

  __m128 v_u = _mm_set1_ps(u);
  __m128 v_v = _mm_set1_ps(v);

  __m128 v_hw = _mm_set1_ps(hw);
  __m128 v_hh = _mm_set1_ps(hh);

  __m128 low  = _mm_sub_ps(_mm_fabs_ps(v_u),v_hw);
  __m128 high = _mm_add_ps(_mm_fabs_ps(v_u),v_hw);

  __m128 scale  = _mm_mul_ps(_m_half_ps,_mm_sub_ps(high,low));
  __m128 offset = _mm_mul_ps(_m_half_ps,_mm_add_ps(high,low));

  const __m128 v_ppmm = _mm_set_ps(-1.0f,-1.0f,1.0f,1.0f);

  adjacent = _mm_add_ps(
               v_hh,
               _mm_mul_ps(
                 v_ppmm,
                 _mm_fabs_ps(
                   v_v)));
  // adjacent[0] == adjacent[1], adjacent[2] == adjacent[3]

  __m128 intReal = _mm_setzero_ps();
  __m128 intImag = _mm_setzero_ps();

  __m128 real, imag;

  __m128 adjacent2 = _mm_square_ps(adjacent);

  __m128 cargz, sargz;

  _mm_sin_cos_ps(_mm_mul_ps(v_mk,v_z),&sargz,&cargz);

  __m128 rcp_denom1 = _mm_rcp_ps(_mm_mul_ps(_m_2pi_ps,_mm_neg_ps(v_mk)));


  for(size_t iuv = 0 ; iuv < 2*(nS/2) ; iuv+=2) {
    __m128 v_uvs     = _mm_set_ps(s[iuv+1],s[iuv],s[iuv+1],s[iuv]);
    __m128 v_weights = _mm_set_ps(weights[iuv+1],weights[iuv],weights[iuv+1],weights[iuv]);

    // ss[0]==ss[2], ss[1]==ss[3]
    __m128 ss  = _mm_add_ps(_mm_mul_ps(scale,v_uvs),offset);

    __m128 ss2 = _mm_square_ps(ss);

    __m128 argh = _mm_mul_ps(
                    v_mk,
                    _mm_sqrt_ps(
                      _mm_add_ps(
                        _mm_add_ps(
                          ss2,
                          v_z2),
                        adjacent2)));

    __m128 cargh, sargh;

    _mm_sin_cos_ps(argh, &sargh, &cargh);

    __m128 denom = _mm_add_ps(ss2,adjacent2);
    __m128 rcp_denom = _mm_rcp_ps(denom);

    __m128 mask_denom = _mm_cmplt_ps(_mm_fabs_ps(scale),_m_eps_ps);
    rcp_denom = _mm_sel_ps(rcp_denom,_mm_setzero_ps(),mask_denom);

    real = _mm_mul_ps(
             _mm_mul_ps(
               v_weights,
               _mm_sub_ps(
                 cargh,
                 cargz)),
             rcp_denom);
    imag = _mm_mul_ps(
             _mm_mul_ps(
               v_weights,
               _mm_sub_ps(
                 sargh,
                 sargz)),
             rcp_denom);
    intReal = _mm_add_ps(intReal, real);
    intImag = _mm_add_ps(intImag, imag);
  }

  //  if (__builtin_expect(nS % 2 == 1, 0)) {
  if (SPS_UNLIKELY(nS % 2 == 1)) {
    // Unlikely
    __m128 v_uvs     = _mm_set1_ps(s[nS-1]);
    __m128 v_weights = _mm_set_ps(0,weights[nS-1],0,weights[nS-1]);

    __m128 ss  = _mm_add_ps(_mm_mul_ps(scale,v_uvs),offset);

    __m128 ss2 = _mm_square_ps(ss);

    __m128 argh = _mm_mul_ps(
                    v_mk,
                    _mm_sqrt_ps(
                      _mm_add_ps(
                        _mm_add_ps(
                          ss2,
                          v_z2),
                        adjacent2)));

    __m128 cargh, sargh;

    _mm_sin_cos_ps(argh, &sargh, &cargh);

    __m128 denom = _mm_add_ps(ss2,adjacent2);
    __m128 rcp_denom = _mm_rcp_ps(denom);

    __m128 mask_denom = _mm_cmplt_ps(_mm_fabs_ps(scale),_m_eps_ps);
    rcp_denom = _mm_sel_ps(rcp_denom,_mm_setzero_ps(),mask_denom);

    real = _mm_mul_ps(
             _mm_mul_ps(
               v_weights,
               _mm_sub_ps(
                 cargh,
                 cargz)),
             rcp_denom);
    imag = _mm_mul_ps(
             _mm_mul_ps(
               v_weights,
               _mm_sub_ps(
                 sargh,
                 sargz)),
             rcp_denom);
    intReal = _mm_add_ps(intReal, real);
    intImag = _mm_add_ps(intImag, imag);
  }

  __m128 intScale = _mm_mul_ps(
                      _mm_mul_ps(
                        scale,
                        adjacent),
                      rcp_denom1);

  // Divide by denominator
  intReal = _mm_mul_ps(intReal,intScale);
  intImag = _mm_mul_ps(intImag,intScale);

  // Multiply by -i
  __m128 tmp = intReal;
  intReal = intImag;
  intImag = _mm_neg_ps(tmp);

  // Horizontal sums (two first integrals)
  retval.real(_mm_cvtss_f32(_mm_dp_ps(_m_one_ps,intReal,0xF1)));
  retval.imag(_mm_cvtss_f32(_mm_dp_ps(_m_one_ps,intImag,0xF1)));

  return retval;
}

// TODO: Call CalcFastFourAny2 twice
template <>
std::complex<float>
inline CalcFourFast(const sps::element_rect_t<float> &__restrict element,
                    const sps::point_t<float> &__restrict projection,
                    const float &__restrict k,
                    const float* __restrict uvs,
                    const float* __restrict uvweights,
                    const size_t nUVs)
{

  std::complex<float> retval = std::complex<float>(0.0f,0.0f);

  // Temporaries
  const float u  = projection[0];
  const float v  = projection[1];
  const float z  = projection[2];
  const float hw = element.hw;
  const float hh = element.hh;

  const __m128 v_mk = _mm_set1_ps(-k);
  const __m128 v_z  = _mm_set1_ps(z);

  const __m128 v_z2 = _mm_square_ps(v_z);

  __m128 adjacent;

  __m128 v_uuvv = _mm_set_ps(v,v,u,u);
  __m128 v_wwhh = _mm_set_ps(hh,hh,hw,hw);

  __m128 low  = _mm_sub_ps(_mm_fabs_ps(v_uuvv),v_wwhh);
  __m128 high = _mm_add_ps(_mm_fabs_ps(v_uuvv),v_wwhh);

  __m128 scale  = _mm_mul_ps(_m_half_ps,_mm_sub_ps(high,low));
  __m128 offset = _mm_mul_ps(_m_half_ps,_mm_add_ps(high,low));

  const __m128 v_pmpm = _mm_set_ps(1.0f,-1.0f,1.0f,-1.0f);

  adjacent = _mm_add_ps(
               _mm_permute_ps(v_wwhh, 0x1B),
               _mm_mul_ps(
                 v_pmpm,
                 _mm_fabs_ps(
                   _mm_permute_ps(v_uuvv, 0x1B))));

  __m128 intReal = _mm_setzero_ps();
  __m128 intImag = _mm_setzero_ps();

  __m128 real, imag;

  __m128 adjacent2 = _mm_square_ps(adjacent);

  __m128 cargz, sargz;

  _mm_sin_cos_ps(_mm_mul_ps(v_mk,v_z),&sargz,&cargz);

  __m128 rcp_denom1 = _mm_rcp_ps(_mm_mul_ps(_m_2pi_ps,_mm_neg_ps(v_mk)));

  for(size_t iuv = 0 ; iuv < nUVs ; iuv++) {

    __m128 v_uvs   = _mm_load1_ps(&uvs[iuv]);
    __m128 weights = _mm_load1_ps(&uvweights[iuv]);

    __m128 ss  = _mm_add_ps(_mm_mul_ps(scale,v_uvs),offset);

    __m128 ss2 = _mm_square_ps(ss);

    __m128 argh = _mm_mul_ps(v_mk,
                             _mm_sqrt_ps(
                               _mm_add_ps(
                                 _mm_add_ps(
                                   ss2,
                                   v_z2),
                                 adjacent2)));

    __m128 cargh, sargh;

    _mm_sin_cos_ps(argh, &sargh, &cargh);

    __m128 denom = _mm_add_ps(ss2,adjacent2);
    __m128 rcp_denom = _mm_rcp_ps(denom);

    __m128 mask_denom = _mm_cmplt_ps(_mm_fabs_ps(scale),_m_eps_ps);
    rcp_denom = _mm_sel_ps(rcp_denom,_mm_setzero_ps(),mask_denom);

    real = _mm_mul_ps(
             _mm_mul_ps(
               weights,
               _mm_sub_ps(
                 cargh,
                 cargz)),
             rcp_denom);
    imag = _mm_mul_ps(
             _mm_mul_ps(
               weights,
               _mm_sub_ps(
                 sargh,
                 sargz)),
             rcp_denom);
    intReal = _mm_add_ps(intReal, real);
    intImag = _mm_add_ps(intImag, imag);
  }

  __m128 intScale = _mm_mul_ps(
                      _mm_mul_ps(
                        scale,
                        adjacent),
                      rcp_denom1);

  // Divide by denominator
  intReal = _mm_mul_ps(intReal,intScale);
  intImag = _mm_mul_ps(intImag,intScale);

  // Multiply by -i
  __m128 tmp = intReal;
  intReal = intImag;
  intImag = _mm_neg_ps(tmp);

  // Horizontal sums
  retval.real(_mm_cvtss_f32(_mm_dp_ps(_m_one_ps,intReal,0xF1)));
  retval.imag(_mm_cvtss_f32(_mm_dp_ps(_m_one_ps,intImag,0xF1)));

  return retval;
}

