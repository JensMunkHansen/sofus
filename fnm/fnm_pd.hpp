/**
 * @file   fnm_pd.hpp
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
std::complex<double>
inline CalcFourFast(const sps::element_t<double> &__restrict element,
                    const sps::point_t<double> &__restrict projection,
                    const double &__restrict k,
                    const double* __restrict uvs,
                    const double* __restrict uvweights,
                    const size_t nUVs)
{
  SPS_UNREFERENCED_PARAMETERS(element, projection, k, uvs, uvweights, nUVs);
  std::complex<double> retval = std::complex<double>(0.0f,0.0f);
  return retval;
}

// Works
template <>
std::complex<double>
inline CalcHzFast(const sps::element_t<double> &__restrict element,
                  const sps::point_t<double> &__restrict projection,
                  const double &__restrict k,
                  const double* __restrict us,
                  const double* __restrict uweights,
                  const size_t nUs,
                  const double* __restrict vs,
                  const double* __restrict vweights,
                  const size_t nVs)
{
  std::complex<double> retval = std::complex<double>(0.0f,0.0f);

  // Temporaries
  const double x  = projection[0];
  const double y  = projection[1];
  const double z  = projection[2];
  const double hw = element.hw;
  const double hh = element.hh;

  const __m128d v_mk = _mm_set1_pd(-k);
  const __m128d v_z  = _mm_set1_pd(z);
  const __m128d v_z2 = _mm_square_pd(v_z);

  __m128d s_stop_in;
  __m128d l_stop_in;
  __m128d s_start_out;

  __m128d s_stop_out;
  __m128d l_start_out;
  __m128d l_stop_out;

  const __m128d v_absx = _mm_set1_pd(fabs(x));
  const __m128d v_absy = _mm_set1_pd(fabs(y));
  const __m128d v_hw   = _mm_set1_pd(hw);
  const __m128d v_hh   = _mm_set1_pd(hh);

  // Lower part (split here)
  const __m128d v_pm = _mm_set_pd(1.0f,-1.0f);
  const __m128d v_p0 = _mm_set_pd(1.0f,0.0f);

  s_stop_in   = _mm_sub_pd(v_hw,_mm_mul_pd(v_pm,v_absx));
  l_stop_in   = _mm_add_pd(v_hh,v_absy);

  s_start_out = _mm_sub_pd(v_absx,_mm_mul_pd(v_p0,v_hw));
  s_stop_out  = _mm_add_pd(v_absx,_mm_mul_pd(_mm_sub_pd(_m_one_pd,v_p0),v_hw));

  l_start_out = v_absy;
  l_stop_out  = _mm_add_pd(v_absy, v_hh);

  // Integral (stop, inside)
  __m128d start_in = _mm_setzero_pd();

  // |x| > hw
  __m128d mask = _mm_cmpgt_pd(v_absx,v_hw);

  __m128d s_start = _mm_sel_pd(start_in, s_start_out, mask);
  __m128d s_stop  = _mm_sel_pd(s_stop_in, s_stop_out, mask);

  // |y| > hh
  __m128d mask1   = _mm_cmpgt_pd(v_absy,v_hh);

  __m128d l_start = _mm_sel_pd(start_in, l_start_out, mask1);
  __m128d l_stop  = _mm_sel_pd(l_stop_in, l_stop_out, mask1);

  __m128d s = _mm_sub_pd(s_stop, s_start);
  __m128d l = _mm_sub_pd(l_stop, l_start);

  __m128d cargz, sargz;

  _mm_sin_cos_pd(_mm_mul_pd(v_mk,v_z),&sargz,&cargz);

  __m128d v_l2 = _mm_square_pd(l_stop_in);
  __m128d v_s2 = _mm_square_pd(s_stop_in);

  __m128d rcp_denom1 = _mm_rcp_pd(_mm_mul_pd(_m_2pi_pd,_mm_neg_pd(v_mk)));

  // u-integral, s-integral
  __m128d s_offset = _mm_mul_pd(_m_half_pd,_mm_add_pd(s_stop,s_start));
  __m128d s_scale  = _mm_mul_pd(_m_half_pd,_mm_sub_pd(s_stop,s_start));

  // v-integral, l-integral
  __m128d l_offset = _mm_mul_pd(_m_half_pd,_mm_add_pd(l_start,l_stop));
  __m128d l_scale  = _mm_mul_pd(_m_half_pd,_mm_sub_pd(l_stop,l_start));

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
                     v_mk,
                     _mm_sqrt_pd(
                       _mm_add_pd(
                         _mm_add_pd(
                           ls2,
                           v_z2),
                         v_l2)));

    __m128d cargw, sargw;

    _mm_sin_cos_pd(argw, &sargw, &cargw);

    __m128d denom = _mm_add_pd(ls2,v_l2);
    __m128d rcp_denom = _mm_rcp_pd(denom);
    // If the integral has zero length, reciprocal of denominator is set to zero
    __m128d mask_denom = _mm_cmplt_pd(_mm_fabs_pd(s_scale),_m_eps_pd);
    rcp_denom = _mm_sel_pd(rcp_denom,_mm_setzero_pd(),mask_denom);

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
                         _mm_sel_pd(s_scale,_mm_mul_pd(s,_m_half_pd),mask),
                         l_stop_in),
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

    __m128d argh = _mm_mul_pd(v_mk,
                              _mm_sqrt_pd(
                                _mm_add_pd(
                                  _mm_add_pd(
                                    ss2,
                                    v_z2),
                                  v_s2)));

    __m128d cargh, sargh;

    _mm_sin_cos_pd(argh, &sargh, &cargh);

    __m128d denom = _mm_add_pd(ss2,v_s2);
    __m128d rcp_denom = _mm_rcp_pd(denom);
    // If the integral has zero length, reciprocal of denominator is set to zero
    __m128d mask_denom = _mm_cmplt_pd(_mm_fabs_pd(l_scale),_m_eps_pd);
    rcp_denom = _mm_sel_pd(rcp_denom,_mm_setzero_pd(),mask_denom);

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
                 _mm_sel_pd(l_scale,_mm_mul_pd(_m_half_pd,l),mask1),
                 s_stop_in),
               rcp_denom1);

  // Question: Can we postpone this???

  // Divide by denominator
  intHreal = _mm_mul_pd(intHreal,intScale);
  intHimag = _mm_mul_pd(intHimag,intScale);

  intHreal = _mm_add_pd(intHreal,intWreal);
  intHimag = _mm_add_pd(intHimag,intWimag);

  // Upper part
  // s_stop_in   = _mm_sub_pd(v_hw, _mm_mul_pd(v_pm,v_absx)); /* Unchaged */
  l_stop_in   = _mm_sub_pd(v_hh, v_absy);

  s_start_out = _mm_sub_pd(v_absx,_mm_mul_pd(v_p0,v_hw));
  s_stop_out  = _mm_add_pd(v_absx,_mm_mul_pd(_mm_sub_pd(_m_one_pd,v_p0),v_hw));

  l_start_out = _mm_sub_pd(v_absy, v_hh);
  l_stop_out  = v_absy;

  // |x| > hw
  mask = _mm_cmpgt_pd(v_absx,v_hw);

  s_start = _mm_sel_pd(start_in, s_start_out, mask);
  s_stop  = _mm_sel_pd(s_stop_in, s_stop_out, mask);

  // |y| > hh
  mask1   = _mm_cmpgt_pd(v_absy,v_hh);

  l_start = _mm_sel_pd(start_in, l_start_out, mask1);
  l_stop  = _mm_sel_pd(l_stop_in, l_stop_out, mask1);

  s = _mm_sub_pd(s_stop, s_start); // Always positive
  l = _mm_sub_pd(l_stop, l_start); // Always positive

  _mm_sin_cos_pd(_mm_mul_pd(v_mk,v_z),&sargz,&cargz);

  // HERE
  v_l2 = _mm_square_pd(l_stop_in);
  v_s2 = _mm_square_pd(s_stop_in);

  rcp_denom1 = _mm_rcp_pd(_mm_mul_pd(_m_2pi_pd,_mm_neg_pd(v_mk)));

  // u-integral, s-integral
  s_offset = _mm_mul_pd(_m_half_pd,_mm_add_pd(s_stop,s_start));
  s_scale  = _mm_mul_pd(_m_half_pd,_mm_sub_pd(s_stop,s_start));

  // v-integral, l-integral
  l_offset = _mm_mul_pd(_m_half_pd,_mm_add_pd(l_start,l_stop));
  l_scale  = _mm_mul_pd(_m_half_pd,_mm_sub_pd(l_stop,l_start));

  intWreal = _mm_setzero_pd();
  intWimag = _mm_setzero_pd();

  for (size_t iu = 0 ; iu < nUs ; iu++) {

    __m128d us1       = _mm_load1_pd(&us[iu]);
    __m128d uweights1 = _mm_load1_pd(&uweights[iu]);

    // [0 ; |x|+hw] , [0 ; hw - |x|];
    __m128d ls  = _mm_add_pd(_mm_mul_pd(s_scale,us1),s_offset);

    __m128d ls2 = _mm_square_pd(ls);

    __m128d argw = _mm_mul_pd(
                     v_mk,
                     _mm_sqrt_pd(
                       _mm_add_pd(
                         _mm_add_pd(
                           ls2,
                           v_z2),
                         v_l2)));

    __m128d cargw, sargw;

    _mm_sin_cos_pd(argw, &sargw, &cargw);

    __m128d denom     = _mm_add_pd(ls2,v_l2);
    __m128d rcp_denom = _mm_rcp_pd(denom);
    // If the integral has zero length, reciprocal of denominator is set to zero
    __m128d mask_denom = _mm_cmplt_pd(_mm_fabs_pd(s_scale),_m_eps_pd);
    rcp_denom = _mm_sel_pd(rcp_denom,_mm_setzero_pd(),mask_denom);

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
                 _mm_sel_pd(s_scale,_mm_mul_pd(s,_m_half_pd),mask),
                 l_stop_in),
               rcp_denom1);

  // Divide by denominator
  intWreal = _mm_mul_pd(intWreal,intScale);
  intWimag = _mm_mul_pd(intWimag,intScale);

  // Filter (consider recover sign earlier)
  __m128d _intHreal = _mm_setzero_pd();
  __m128d _intHimag = _mm_setzero_pd();

  for(size_t iv = 0 ; iv < nVs ; iv++) {

    __m128d vs1       = _mm_load1_pd(&vs[iv]);
    __m128d vweights1 = _mm_load1_pd(&vweights[iv]);

    __m128d ss  = _mm_add_pd(_mm_mul_pd(l_scale,vs1),l_offset);

    __m128d ss2 = _mm_square_pd(ss);

    __m128d argh = _mm_mul_pd(v_mk,
                              _mm_sqrt_pd(
                                _mm_add_pd(
                                  _mm_add_pd(
                                    ss2,
                                    v_z2),
                                  v_s2)));

    __m128d cargh, sargh;

    _mm_sin_cos_pd(argh, &sargh, &cargh);

    __m128d denom = _mm_add_pd(ss2,v_s2);
    __m128d rcp_denom = _mm_rcp_pd(denom);
    // If the integral has zero length, reciprocal of denominator is set to zero
    __m128d mask_denom = _mm_cmplt_pd(_mm_fabs_pd(l_scale),_m_eps_pd);
    rcp_denom = _mm_sel_pd(rcp_denom,_mm_setzero_pd(),mask_denom);

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
    _intHreal = _mm_add_pd(_intHreal, real);
    _intHimag = _mm_add_pd(_intHimag, imag);
  }

  intScale = _mm_mul_pd(
               _mm_mul_pd(
                 _mm_sel_pd(l_scale,_mm_mul_pd(_m_half_pd,l),mask1),
                 s_stop_in),
               rcp_denom1);

  // Divide by denominator
  _intHreal = _mm_mul_pd(_intHreal,intScale);
  _intHimag = _mm_mul_pd(_intHimag,intScale);

  _intHreal = _mm_add_pd(_intHreal,intWreal);
  _intHimag = _mm_add_pd(_intHimag,intWimag);

  intHreal = _mm_add_pd(intHreal,_intHreal);
  intHimag = _mm_add_pd(intHimag,_intHimag);

  // Multiply by -i
  __m128d tmp = intHreal;
  intHreal = intHimag;
  intHimag = _mm_neg_pd(tmp);

  // Horizontal sums
#if 0
  _mm_store_sd(
    &(reinterpret_cast<double(&)[2]>(retval)[0]),
    _mm_dp_pd(
      _m_one_pd,
      intHreal,
      0xF1));
  _mm_store_sd(
    &(reinterpret_cast<double(&)[2]>(retval)[1]),
    _mm_dp_pd(
      _m_one_pd,
      intHimag,
      0xF1));
#else
  __m128d tmpd = _mm_dp_pd(_m_one_pd, intHreal, 0xF1);
  tmpd = _mm_add_pd(tmpd, _mm_dp_pd(_m_one_pd, intHimag, 0xF2));
  _mm_store_pd(&(reinterpret_cast<double(&)[2]>(retval)[0]),tmpd);
#endif

  return retval;
}

// Scratch
#if 0
template <>
std::complex<double>
inline CalcHzFast(const sps::element_t<double>& element,
                  const sps::point_t<double>& projection,
                  const double& k,
                  const double* __restrict us,
                  const double* __restrict uweights,
                  const size_t nUs,
                  const double* __restrict vs,
                  const double* __restrict vweights,
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

  const __m128d z2 = _mm_square_pd(vec_z);
  __m128d vec_l2   = _mm_square_pd(l_stop_in);
  __m128d vec_s2   = _mm_square_pd(s_stop_in);

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

    __m128d cargw = _mm_setzero_pd();
    __m128d sargw = _mm_setzero_pd();

    _mm_sin_cos_pd(argw, &sargw, &cargw);

    __m128d rcp_denom = _mm_rcp_pd(_mm_add_pd(ls2,vec_l2));
    __m128d mask_denom = _mm_cmplt_pd(_mm_fabs_pd(s_scale),_m_eps_pd);
    rcp_denom = _mm_sel_pd(rcp_denom,_mm_setzero_pd(),mask_denom);

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

    __m128d mask_denom = _mm_cmplt_pd(_mm_fabs_pd(l_scale),_m_eps_pd);
    rcp_denom = _mm_sel_pd(rcp_denom,_mm_setzero_pd(),mask_denom);

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

  // Horizontal sum (wrong)
  _mm_store_pd(tmpd,_mm_dp_pd(_m_one_pd,intHreal,0xF1));
  retval.real(retval.real() + tmpd[0]);

  _mm_store_pd(tmpd,_mm_dp_pd(_m_one_pd,intHimag,0xF1));
  retval.imag(retval.imag() + tmpd[0]);

  return retval;
}
#endif

template <>
std::complex<double>
inline CalcFastFourAny(const double& u,
                       const double& v,
                       const double& hh,
                       const double& hw,
                       const double& z,
                       const double& __restrict k,
                       const double* __restrict s,
                       const double* __restrict weights,
                       const size_t nS)
{
  SPS_UNREFERENCED_PARAMETERS(u,v,hh,hw,z,k,s,weights,nS);
  std::complex<double> retval = std::complex<double>(0.0f,0.0f);
  return retval;
}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */

