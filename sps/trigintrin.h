/**
 * @file   trigintrin.h
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Thu Feb 12 21:24:24 2015
 *
 * @brief  Trigonometric functions implemented using SSE2/SSE3/SSE4/FMA
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
#include <sps/cenv.h>
#include <sps/math.h>
#include <sps/extintrin.h>

// TODO: Consider making constants extern to use them elsewhere or repeat them

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \defgroup Constants used by functions
 * @{
 */
const __m128 _m_2_pi_ps       = _mm_set1_ps((float)M_2_PI);    // 2 / pi
const __m128 _m_1_pi_ps       = _mm_set1_ps((float)M_1_PI);    // 1 / pi
const __m128 _m_1_2pi_ps      = _mm_set1_ps((float)M_1_2PI);   // 1 / (2pi)
const __m128 _m_2pi_ps        = _mm_set1_ps((float)M_2PI);     // 2 * pi
const __m128 _m_pi_ps         = _mm_set1_ps((float)M_PI);      // pi
const __m128 _m_pi2_ps        = _mm_set1_ps((float)M_PI_2);    // pi/ 2

const __m128 _m_1_4_ps       = _mm_set1_ps(0.25f);
const __m128 _m_1_ps         = _mm_set1_ps(1.0f);
const __m128 _m_2_ps         = _mm_set1_ps(2.0f);

const __m128d _m_2pi_pd      = _mm_set1_pd((double)M_2PI);     // 2 * pi

const __m256d _m256_pi_pd         = _mm256_set1_pd(M_PI);      // pi
const __m256d _m256_1_2pi_pd      = _mm256_set1_pd(M_1_2PI);   // 1 / (2pi)
const __m256d _m256_pi2_pd        = _mm256_set1_pd(M_PI_2);    // pi/ 2

const __m256d _m256_1_pd         = _mm256_set1_pd(1.0f);


/**@}*/

/**
 * Arccos function. Handbook of Mathematical Functions M. Abramowitz
 * and I.A. Stegun, Ed. Absolute error <= 6.7e-5. Execution time is 44
 * clock cycles on Intel(R) Core(TM) i7-2620M (gcc 4.7.2-1)
 *
 * @param x
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128 _mm_arccos_ps(__m128 x) STATIC_INLINE_END;

/**
 * Arccos function.
 *
 * @param x
 *
 * @return
 */
STATIC_INLINE_BEGIN __m256d _mm256_arccos_pd(__m256d x) STATIC_INLINE_END;

/**
 * Arcsine function. Handbook of Mathematical Functions M. Abramowitz
 * and I.A. Stegun, Ed. Absolute error <= 6.7e-5.
 *
 * @param x
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128 _mm_arcsin_ps(__m128 x) STATIC_INLINE_END;

/**
 * Arctan2 function
 *
 *               +
 *              /|
 *             / |
 *            /  | o
 *           /   | p
 *          /    | p
 *   hypotenuse  | o
 *        /      | s
 *       /       | i
 *      /        | t
 *     /         | e
 *    / a        |
 *   +-----------+
 *     adjacent
 *
 * @param y is the opposite
 * @param x is the adjacent
 *
 * @return a
 */
STATIC_INLINE_BEGIN __m128 _mm_arctan2_ps(__m128 y, __m128 x) STATIC_INLINE_END;


/**
 * Exponential function
 *
 * @param d
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128 _mm_exp_ps(__m128 d) STATIC_INLINE_END;


STATIC_INLINE_BEGIN __m256d _mm256_arccos_pd(__m256d x)
{
  // TODO: Increase precision of constants
#if 1
  v4d _x;
  _x.v = x;

  return _mm256_set_pd(acos(_x.f64[3]), acos(_x.f64[2]), acos(_x.f64[1]), acos(_x.f64[0]));
#else
  __m256d mask = _mm256_cmp_pd(x, _mm256_setzero_pd(), _CMP_LT_OQ); // TODO: Manage signals and ordering

  x = _mm256_fabs_pd(x);

  __m256d ret = _mm256_set1_pd(-0.0187293f);
  ret = _mm256_mul_pd(ret,x);
  ret = _mm256_add_pd(ret,_mm256_set1_pd(0.0742610f));
  ret = _mm256_mul_pd(ret,x);
  ret = _mm256_sub_pd(ret,_mm256_set1_pd(0.2121144f));
  ret = _mm256_mul_pd(ret,x);
  ret = _mm256_add_pd(ret, _m256_pi2_pd);
  ret = _mm256_mul_pd(_mm256_sqrt_pd(_mm256_sub_pd(_m256_one_pd,x)),ret);

  return _mm256_sel_pd(ret,_mm256_sub_pd(_m256_pi_pd,ret),mask);
#endif
}

STATIC_INLINE_BEGIN __m128 _mm_arccos_ps(__m128 x)
{
  __m128 mask = _mm_cmplt_ps(x,_mm_setzero_ps());
  x = _mm_fabs_ps(x);

  __m128 ret = _mm_set1_ps(-0.0187293f);
  ret = _mm_mul_ps(ret,x);
  ret = _mm_add_ps(ret,_mm_set1_ps(0.0742610f));
  ret = _mm_mul_ps(ret,x);
  ret = _mm_sub_ps(ret,_mm_set1_ps(0.2121144f));
  ret = _mm_mul_ps(ret,x);
  ret = _mm_add_ps(ret,_m_pi2_ps);
  ret = _mm_mul_ps(_mm_sqrt_ps(_mm_sub_ps(_m_one_ps,x)),ret);

  return _mm_sel_ps(ret,_mm_sub_ps(_m_pi_ps,ret),mask);
}

#if 0
/**
 * Sleefs arctan2 function
 *
 * @param y
 * @param x
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128 atan2kf(__m128 y, __m128 x)
{
  __m128 s, t, u;
  __m128i q;
  __m128 p;

  q = _mm_sel_epi32(_mm_set1_epi32(0),_mm_set1_epi32(-2),_mm_castps_si128(_mm_cmplt_ps(x,_mm_setzero_ps())));

  x = _mm_fabs_ps(x);

  q = _mm_sel_epi32(q,_mm_set1_epi32(1),_mm_castps_si128(_mm_cmplt_ps(x,y)));

  p = _mm_cmplt_ps(x, y);
  s = _mm_sel_ps(y, _mm_neg_ps(x), p);
  t = _mm_max_ps(x, y);

  s = _mm_div_ps(s, t);
  t = _mm_mul_ps(s, s);

  u = _mm_set1_ps(0.00282363896258175373077393f);
  u = _mm_madd_ps(u, t, _mm_set1_ps(-0.0159569028764963150024414f));
  u = _mm_madd_ps(u, t, _mm_set1_ps(0.0425049886107444763183594f));
  u = _mm_madd_ps(u, t, _mm_set1_ps(-0.0748900920152664184570312f));
  u = _mm_madd_ps(u, t, _mm_set1_ps(0.106347933411598205566406f));
  u = _mm_madd_ps(u, t, _mm_set1_ps(-0.142027363181114196777344f));
  u = _mm_madd_ps(u, t, _mm_set1_ps(0.199926957488059997558594f));
  u = _mm_madd_ps(u, t, _mm_set1_ps(-0.333331018686294555664062f));

  t = _mm_madd_ps(s, _mm_mul_ps(t, u), s);
  t = _mm_madd_ps(_mm_cvtepi32_ps(q), _mm_set1_ps((float)(M_PI/2)), t);

  return t;
}

/**
 * Sleefs arcsin function (slower than _mm_arcsin_ps)
 *
 * @param d
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128 xasinf(__m128 d)
{
  __m128 x, y;
  x = _mm_add_ps(_mm_set1_ps(1.0f), d);
  y = _mm_sub_ps(_mm_set1_ps(1.0f), d);
  x = _mm_mul_ps(x, y);
  x = _mm_sqrt_ps(x);
  x = _mm_or_ps(_mm_cmpneq_ps(x,x), atan2kf(_mm_fabs_ps(d), x));
  return _mm_mulsign_ps(x,d);
}

/**
 * Sleefs arccos function (slower than _mm_arccos_ps)
 *
 * @param d
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128 xacosf(__m128 d)
{
  __m128 x, y;
  x = _mm_add_ps(_mm_set1_ps(1.0f), d);
  y = _mm_sub_ps(_mm_set1_ps(1.0f), d);
  x = _mm_mul_ps(x, y);
  x = _mm_sqrt_ps(x);
  x = _mm_mulsign_ps(atan2kf(x, _mm_fabs_ps(d)), d);

  y = _mm_and_ps(_mm_cmplt_ps(d, _mm_setzero_ps()), _mm_set1_ps((float)M_PI));
  x = _mm_add_ps(x, y);
  return x;
}
#endif

STATIC_INLINE_BEGIN __m128 _mm_arcsin_ps(__m128 x)
{
  __m128 mask = _mm_cmplt_ps(x,_mm_setzero_ps());
  x = _mm_fabs_ps(x);
  __m128 ret = _mm_set1_ps(-0.0187293f);
  ret = _mm_mul_ps(ret,x);
  ret = _mm_add_ps(ret,_mm_set1_ps(0.0742610f));
  ret = _mm_mul_ps(ret,x);
  ret = _mm_sub_ps(ret,_mm_set1_ps(0.2121144f));
  ret = _mm_mul_ps(ret,x);
  ret = _mm_add_ps(ret,_m_pi2_ps);
  ret = _mm_sub_ps(_m_pi2_ps,_mm_mul_ps(_mm_sqrt_ps(_mm_sub_ps(_m_one_ps,x)),ret));
  __m128 ret2 = _mm_mul_ps(_m_2_ps,ret);
  ret = _mm_sel_ps(ret,_mm_sub_ps(ret,ret2),mask);
  return ret;
}

STATIC_INLINE_BEGIN __m128 _mm_arctan2_ps(__m128 y, __m128 x)
{
  __m128 t0, t1, t2, t3, t4;

  // Verify this is correct
  SPS_UNREFERENCED_PARAMETER(t2);
  t3 = _mm_fabs_ps(x);
  t1 = _mm_fabs_ps(y);
  __m128 mask = _mm_cmpgt_ps(t1,t3);

  t0 = _mm_max_ps(t3,t1);
  t1 = _mm_min_ps(t3, t1);
  t3 = _mm_rcp_ps(t0);
  t3 = _mm_mul_ps(t1,t3);

  t4 = _mm_mul_ps(t3,t3);
  t0 = _mm_set1_ps(-0.013480470f);
  t0 = _mm_add_ps(_mm_mul_ps(t0,t4),_mm_set1_ps(0.057477314f));

  t0 = _mm_sub_ps(_mm_mul_ps(t0,t4),_mm_set1_ps(0.121239071f));
  t0 = _mm_add_ps(_mm_mul_ps(t0,t4),_mm_set1_ps(0.195635925f));
  t0 = _mm_sub_ps(_mm_mul_ps(t0,t4),_mm_set1_ps(0.332994597f));
  t0 = _mm_add_ps(_mm_mul_ps(t0,t4),_mm_set1_ps(0.999995630f));
  t3 = _mm_mul_ps(t0,t3);

  t3 = _mm_sel_ps(t3,_mm_sub_ps(_m_pi2_ps,t3),mask);
  t3 = _mm_sel_ps(t3,_mm_sub_ps(_m_pi_ps,t3),_mm_cmplt_ps(x,_mm_setzero_ps()));
  t3 = _mm_sel_ps(t3,_mm_neg_ps(t3),_mm_cmplt_ps(y,_mm_setzero_ps()));

  return t3;
}

/**
 * Inaccurate but fast cosine. Error is about 1.0e-3 (slower than _mm_arccos)
 *
 * @param x
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128 _mm_cos_ps_fast(__m128 x)
{
  const __m128 tp = _mm_set1_ps((float) M_1_PI / 2.0f);
  x = _mm_mul_ps(x,tp);
  x = _mm_sub_ps(x,_mm_add_ps(_m_1_4_ps,_mm_floor_ps(_mm_add_ps(x,_m_1_4_ps))));
  x = _mm_mul_ps(x,_mm_mul_ps(_mm_set1_ps(16.0f),_mm_sub_ps(_mm_fabs_ps(x),_m_half_ps)));
  x = _mm_add_ps(x,_mm_mul_ps(_mm_set1_ps(0.225f),_mm_mul_ps(x,_mm_sub_ps(_mm_fabs_ps(x),_m_one_ps))));
  return x;
}

#define PI4_Af 0.78515625f
#define PI4_Bf 0.00024187564849853515625f
#define PI4_Cf 3.7747668102383613586e-08f
#define PI4_Df 1.2816720341285448015e-12f

// How about doxygen

#if ACCURATE_TRIGONOMETRICS

STATIC_INLINE_BEGIN void _mm_sin_cos_ps(__m128 d, __m128* s_, __m128* c_)
{
  ALIGN16_BEGIN float in[4] ALIGN16_END;
  ALIGN16_BEGIN float out[4] ALIGN16_END;
  _mm_store_ps(in,d);
  out[0] = sin(in[0]);
  out[1] = sin(in[1]);
  out[2] = sin(in[2]);
  out[3] = sin(in[3]);
  *s = _mm_load_ps(out);
  out[0] = cos(in[0]);
  out[1] = cos(in[1]);
  out[2] = cos(in[2]);
  out[3] = cos(in[3]);
  *c = _mm_load_ps(out);
}
#else

/**
 * Accurate sine and cosine. Error is less than 6 ulps
 *
 * @param d
 * @param c_
 * @param d_
 */
STATIC_INLINE_BEGIN void _mm_sin_cos_ps(__m128 d, __m128* s_, __m128* c_)
{

  __m128i q,m;
  __m128 u, s, t, rx, ry;

  __m128 r2x, r2y;

  q = _mm_cvtps_epi32(_mm_mul_ps(d, _m_2_pi_ps));

  s = d;

  u = _mm_cvtepi32_ps(q);
  s = _mm_madd_ps(u, _mm_set1_ps(-PI4_Af*2), s);
  s = _mm_madd_ps(u, _mm_set1_ps(-PI4_Bf*2), s);
  s = _mm_madd_ps(u, _mm_set1_ps(-PI4_Cf*2), s);
  s = _mm_madd_ps(u, _mm_set1_ps(-PI4_Df*2), s);

  t = s;

  s = _mm_mul_ps(s, s);

  u = _mm_set1_ps(-0.000195169282960705459117889f);
  u = _mm_madd_ps(u, s, _mm_set1_ps(0.00833215750753879547119141f));
  u = _mm_madd_ps(u, s, _mm_set1_ps(-0.166666537523269653320312f));
  u = _mm_mul_ps(_mm_mul_ps(u, s), t);

  rx = _mm_add_ps(t, u);

  u = _mm_set1_ps(-2.71811842367242206819355e-07f);
  u = _mm_madd_ps(u, s, _mm_set1_ps(2.47990446951007470488548e-05f));
  u = _mm_madd_ps(u, s, _mm_set1_ps(-0.00138888787478208541870117f));
  u = _mm_madd_ps(u, s, _mm_set1_ps(0.0416666641831398010253906f));
  u = _mm_madd_ps(u, s, _mm_set1_ps(-0.5));

  ry = _mm_madd_ps(s, u, _mm_set1_ps(1));

  m = _mm_cmpeq_epi32(_mm_and_si128(q, _mm_set1_epi32(1)), _mm_set1_epi32(0));

  r2x = _mm_sel_ps(ry,rx,_mm_castsi128_ps(m));
  r2y = _mm_sel_ps(rx,ry,_mm_castsi128_ps(m));

  m = _mm_cmpeq_epi32(_mm_and_si128(q, _mm_set1_epi32(2)), _mm_set1_epi32(2));
  r2x = _mm_castsi128_ps(_mm_xor_si128(_mm_and_si128(m, _mm_castps_si128(_mm_set1_ps(-0.0))), _mm_castps_si128(r2x)));

  m = _mm_cmpeq_epi32(_mm_and_si128(_mm_add_epi32(q, _mm_set1_epi32(1)), _mm_set1_epi32(2)), _mm_set1_epi32(2));
  r2y = _mm_castsi128_ps(_mm_xor_si128(_mm_and_si128(m, _mm_castps_si128(_mm_set1_ps(-0.0))), _mm_castps_si128(r2y)));

  __m128 m1 = _mm_is_infinity(d);

  r2x = _mm_or_ps(m1, r2x);
  r2y = _mm_or_ps(m1, r2y);

  *s_ = r2x;
  *c_ = r2y;
}
#endif

STATIC_INLINE_BEGIN void _mm_sin_cos_pd(__m128d d, __m128d* s_, __m128d* c_)
{
  ALIGN16_BEGIN double data[4] ALIGN16_END;
  _mm_store_pd(data, d);
  __m128 d_ps = _mm_set_ps(0.0f,0.0f,(float)data[1],(float)data[0]);
  __m128 c_ps, s_ps;
  _mm_sin_cos_ps(d_ps,&s_ps,&c_ps);
  ALIGN16_BEGIN float fdata[4] ALIGN16_END;
  _mm_store_ps(fdata,s_ps);
  *s_ = _mm_set_pd(fdata[1],fdata[0]);
  _mm_store_ps(fdata,c_ps);
  *c_ = _mm_set_pd(fdata[1],fdata[0]);
}

/**
 * Accurate sine
 *
 * @param d
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128 _mm_sin_ps(__m128 d)
{
  __m128i q;
  __m128 u, s;

  q = _mm_cvtps_epi32(
        _mm_mul_ps(
          d,
          _mm_set1_ps((float)M_1_PI)));
  u = _mm_cvtepi32_ps(q);

  d = _mm_madd_ps(u, _mm_set1_ps(-PI4_Af*4), d);
  d = _mm_madd_ps(u, _mm_set1_ps(-PI4_Bf*4), d);
  d = _mm_madd_ps(u, _mm_set1_ps(-PI4_Cf*4), d);
  d = _mm_madd_ps(u, _mm_set1_ps(-PI4_Df*4), d);

  s = _mm_mul_ps(d, d);

  // Float masking instead
  d = _mm_castsi128_ps(
        _mm_xor_si128(
          _mm_and_si128(
            _mm_cmpeq_epi32(
              _mm_and_si128(
                q,
                _mm_set1_epi32(1)),
              _mm_set1_epi32(1)),
            _mm_castps_si128(_mm_set1_ps(-0.0f))),
          _mm_castps_si128(d)));

  u = _mm_set1_ps(2.6083159809786593541503e-06f);
  u = _mm_madd_ps(u, s, _mm_set1_ps(-0.0001981069071916863322258f));
  u = _mm_madd_ps(u, s, _mm_set1_ps(0.00833307858556509017944336f));
  u = _mm_madd_ps(u, s, _mm_set1_ps(-0.166666597127914428710938f));

  u = _mm_madd_ps(s, _mm_mul_ps(u, d), d);

  u = _mm_or_ps(_mm_is_infinity(d),u);

  return u;
}

STATIC_INLINE_BEGIN __m128 _mm_cos_ps(__m128 d)
{
  __m128i q;
  __m128 u, s;

  q = _mm_cvtps_epi32(
        _mm_sub_ps(
          _mm_mul_ps(
            d,
            _mm_set1_ps((float)M_1_PI)),
          _mm_set1_ps(0.5f)));

  q = _mm_add_epi32(
        _mm_add_epi32(
          q,
          q),
        _mm_set1_epi32(1));

  u = _mm_cvtepi32_ps(q);
  d = _mm_madd_ps(u, _mm_set1_ps(-PI4_Af*2), d);
  d = _mm_madd_ps(u, _mm_set1_ps(-PI4_Bf*2), d);
  d = _mm_madd_ps(u, _mm_set1_ps(-PI4_Cf*2), d);
  d = _mm_madd_ps(u, _mm_set1_ps(-PI4_Df*2), d);

  s = _mm_mul_ps(d, d);

  // TODO: Perform float masking instead
  d =
    _mm_castsi128_ps(
      _mm_xor_si128(
        _mm_and_si128(
          _mm_cmpeq_epi32(
            _mm_and_si128(
              q,
              _mm_set1_epi32(2)),
            _mm_set1_epi32(0)),
          _mm_castps_si128(_mm_set1_ps(-0.0f))),
        _mm_castps_si128(d)));

  u = _mm_set1_ps(2.6083159809786593541503e-06f);
  u = _mm_madd_ps(u, s, _mm_set1_ps(-0.0001981069071916863322258f));
  u = _mm_madd_ps(u, s, _mm_set1_ps(0.00833307858556509017944336f));
  u = _mm_madd_ps(u, s, _mm_set1_ps(-0.166666597127914428710938f));

  u = _mm_madd_ps(s, _mm_mul_ps(u, d), d);

  u = _mm_or_ps(_mm_is_infinity(d), u);

  return u;
}


#define _PS_CONST_TYPE(Name, Type, Val)                                 \
  static const ALIGN16_BEGIN Type _ps_##Name[4] ALIGN16_END = { Val, Val, Val, Val }
#define _PS_CONST(Name, Val)                                            \
  static const ALIGN16_BEGIN float _ps_##Name[4] ALIGN16_END = { Val, Val, Val, Val }
#define _PI32_CONST(Name, Val)                                            \
  static const ALIGN16_BEGIN int _pi32_##Name[4] ALIGN16_END = { Val, Val, Val, Val }

_PS_CONST_TYPE(inv_sign_mask, int, ~0x80000000);
_PS_CONST_TYPE(sign_mask, int, (int)0x80000000);

_PS_CONST(sincof_p0, -1.9515295891E-4f);
_PS_CONST(sincof_p1,  8.3321608736E-3f);
_PS_CONST(sincof_p2, -1.6666654611E-1f);
_PS_CONST(coscof_p0,  2.443315711809948E-005f);
_PS_CONST(coscof_p1, -1.388731625493765E-003f);
_PS_CONST(coscof_p2,  4.166664568298827E-002f);
_PS_CONST(cephes_FOPI, 1.27323954473516f); // 4 / M_PI
_PS_CONST(minus_cephes_DP1, -0.78515625f);
_PS_CONST(minus_cephes_DP2, -2.4187564849853515625e-4f);
_PS_CONST(minus_cephes_DP3, -3.77489497744594108e-8f);

_PI32_CONST(1, 1);
_PI32_CONST(inv1, ~1);
_PI32_CONST(2, 2);
_PI32_CONST(4, 4);
_PS_CONST(1  , 1.0f);
_PS_CONST(0p5, 0.5f);

STATIC_INLINE_BEGIN void _mm_sincos_cephes_ps(__m128 x, __m128 *s, __m128 *c)
{

  __m128 xmm1, xmm2, xmm3 = _mm_setzero_ps(), sign_bit_sin, y;
  __m128i emm0, emm2, emm4;

  sign_bit_sin = x;
  /* take the absolute value */
  x = _mm_and_ps(x, *(__m128*)_ps_inv_sign_mask);
  /* extract the sign bit (upper one) */
  sign_bit_sin = _mm_and_ps(sign_bit_sin, *(__m128*)_ps_sign_mask);

  /* scale by 4/Pi */
  y = _mm_mul_ps(x, *(__m128*)_ps_cephes_FOPI);

  /* store the integer part of y in emm2 */
  emm2 = _mm_cvttps_epi32(y);

  /* j=(j+1) & (~1) (see the cephes sources) */
  emm2 = _mm_add_epi32(emm2, *(__m128i*)_pi32_1);
  emm2 = _mm_and_si128(emm2, *(__m128i*)_pi32_inv1);
  y = _mm_cvtepi32_ps(emm2);

  emm4 = emm2;

  /* get the swap sign flag for the sine */
  emm0 = _mm_and_si128(emm2, *(__m128i*)_pi32_4);
  emm0 = _mm_slli_epi32(emm0, 29);
  __m128 swap_sign_bit_sin = _mm_castsi128_ps(emm0);

  /* get the polynom selection mask for the sine*/
  emm2 = _mm_and_si128(emm2, *(__m128i*)_pi32_2);
  emm2 = _mm_cmpeq_epi32(emm2, _mm_setzero_si128());
  __m128 poly_mask = _mm_castsi128_ps(emm2);

  /* The magic pass: "Extended precision modular arithmetic"
     x = ((x - y * DP1) - y * DP2) - y * DP3; */
  xmm1 = *(__m128*)_ps_minus_cephes_DP1;
  xmm2 = *(__m128*)_ps_minus_cephes_DP2;
  xmm3 = *(__m128*)_ps_minus_cephes_DP3;
  xmm1 = _mm_mul_ps(y, xmm1);
  xmm2 = _mm_mul_ps(y, xmm2);
  xmm3 = _mm_mul_ps(y, xmm3);
  x = _mm_add_ps(x, xmm1);
  x = _mm_add_ps(x, xmm2);
  x = _mm_add_ps(x, xmm3);

  emm4 = _mm_sub_epi32(emm4, *(__m128i*)_pi32_2);
  emm4 = _mm_andnot_si128(emm4, *(__m128i*)_pi32_4);
  emm4 = _mm_slli_epi32(emm4, 29);
  __m128 sign_bit_cos = _mm_castsi128_ps(emm4);

  sign_bit_sin = _mm_xor_ps(sign_bit_sin, swap_sign_bit_sin);


  /* Evaluate the first polynom  (0 <= x <= Pi/4) */
  __m128 z = _mm_mul_ps(x,x);
  y = *(__m128*)_ps_coscof_p0;

  y = _mm_mul_ps(y, z);
  y = _mm_add_ps(y, *(__m128*)_ps_coscof_p1);
  y = _mm_mul_ps(y, z);
  y = _mm_add_ps(y, *(__m128*)_ps_coscof_p2);
  y = _mm_mul_ps(y, z);
  y = _mm_mul_ps(y, z);
  __m128 tmp = _mm_mul_ps(z, *(__m128*)_ps_0p5);
  y = _mm_sub_ps(y, tmp);
  y = _mm_add_ps(y, *(__m128*)_ps_1);

  /* Evaluate the second polynom  (Pi/4 <= x <= 0) */

  __m128 y2 = *(__m128*)_ps_sincof_p0;
  y2 = _mm_mul_ps(y2, z);
  y2 = _mm_add_ps(y2, *(__m128*)_ps_sincof_p1);
  y2 = _mm_mul_ps(y2, z);
  y2 = _mm_add_ps(y2, *(__m128*)_ps_sincof_p2);
  y2 = _mm_mul_ps(y2, z);
  y2 = _mm_mul_ps(y2, x);
  y2 = _mm_add_ps(y2, x);

  // Can be made faster using _mm_blendv_ps

  /* select the correct result from the two polynoms */
  xmm3 = poly_mask;
  __m128 ysin2 = _mm_and_ps(xmm3, y2);
  __m128 ysin1 = _mm_andnot_ps(xmm3, y);
  y2 = _mm_sub_ps(y2,ysin2);
  y = _mm_sub_ps(y, ysin1);

  xmm1 = _mm_add_ps(ysin1,ysin2);
  xmm2 = _mm_add_ps(y,y2);

  /* update the sign */
  *s = _mm_xor_ps(xmm1, sign_bit_sin);
  *c = _mm_xor_ps(xmm2, sign_bit_cos);
}

/**
 * Multiplies the packed singles x by the number 2 raised to the q power
 *
 * @param x
 * @param q
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128 _mm_ldexpf(__m128 x, __m128i q) STATIC_INLINE_END;

STATIC_INLINE_BEGIN __m256d _mm256_ldexpd(__m256d x, __m128i q) STATIC_INLINE_END;


STATIC_INLINE_BEGIN __m128 _mm_exp_approx_ps(__m128 x)
{
  __m128 retval = _mm_setzero_ps();
  __m128 part1  = _mm_mul_ps(x, x);
  __m128 part2  = _mm_mul_ps(x, part1);
  retval =
    _mm_add_ps(_mm_add_ps(_mm_add_ps(_m_one_ps,x),_mm_mul_ps(part1,_m_half_ps)),_mm_mul_ps(part2,_m_sixth_ps));
  return retval;
}

#define L2Uf 0.693145751953125f
#define L2Lf 1.428606765330187045e-06f
#define R_LN2f 1.442695040888963407359924681001892137426645954152985934135449406931f

STATIC_INLINE_BEGIN __m128 _mm_exp_ps(__m128 d)
{
  __m128i q = _mm_cvtps_epi32(_mm_mul_ps(d, _mm_set1_ps(R_LN2f)));
  __m128 s, u;

  s = _mm_madd_ps(_mm_cvtepi32_ps(q), _mm_set1_ps(-L2Uf), d);
  s = _mm_madd_ps(_mm_cvtepi32_ps(q), _mm_set1_ps(-L2Lf), s);

  u = _mm_set1_ps(0.00136324646882712841033936f);
  u = _mm_madd_ps(u, s, _mm_set1_ps(0.00836596917361021041870117f));
  u = _mm_madd_ps(u, s, _mm_set1_ps(0.0416710823774337768554688f));
  u = _mm_madd_ps(u, s, _mm_set1_ps(0.166665524244308471679688f));
  u = _mm_madd_ps(u, s, _mm_set1_ps(0.499999850988388061523438f));

  u = _mm_add_ps(_mm_set1_ps(1.0f), _mm_madd_ps(_mm_mul_ps(s, s), u, s));

  u = _mm_ldexpf(u, q);

  u = _mm_andnot_ps(_mm_is_minfinity(d), u);

  return u;
}

/*
STATIC_INLINE_BEGIN __m256i _mm256_cvtpd_epi64(__m256d d) {

}
*/

STATIC_INLINE_BEGIN __m256d _mm256_madd_pd(__m256d a, __m256d b, __m256d c)
{
  return _mm256_add_pd(_mm256_mul_pd(a,b),c);
}

STATIC_INLINE_BEGIN __m256d _mm256_exp_pd(__m256d d)
{
  __m256d u;
  v4d _d;
  _d.v = d;

  u = _mm256_set_pd(exp(_d.f64[3]),exp(_d.f64[2]),exp(_d.f64[1]),exp(_d.f64[0]));

#if 0
  // Need AVX-512 for _mm256_cvtpd_epi64
  __m128i q = _mm256_cvtpd_epi32(_mm256_mul_pd(d, _mm256_set1_pd(R_LN2f)));
  __m256d s;

  // Need AVX-512 for _mm256_madd_pd
  s = _mm256_madd_pd(_mm256_cvtepi32_pd(q), _mm256_set1_pd(-L2Uf), d);
  s = _mm256_madd_pd(_mm256_cvtepi32_pd(q), _mm256_set1_pd(-L2Lf), s);

  u = _mm256_set1_pd(0.00136324646882712841033936f);
  u = _mm256_madd_pd(u, s, _mm256_set1_pd(0.00836596917361021041870117f));
  u = _mm256_madd_pd(u, s, _mm256_set1_pd(0.0416710823774337768554688f));
  u = _mm256_madd_pd(u, s, _mm256_set1_pd(0.166665524244308471679688f));
  u = _mm256_madd_pd(u, s, _mm256_set1_pd(0.499999850988388061523438f));

  u = _mm256_add_pd(_mm256_set1_pd(1.0f), _mm256_madd_pd(_mm256_mul_pd(s, s), u, s));

  // Not implemented
  u = _mm256_ldexpd(u, q);

  u = _mm256_andnot_pd(_mm256_is_minfinity(d), u);
#endif
  return u;
}




// TODO: Make a faster version
_PI32_CONST(0x7f, 0x7f);

_PS_CONST(exp_hi,   88.3762626647949f);
_PS_CONST(exp_lo,   -88.3762626647949f);

_PS_CONST(cephes_LOG2EF, 1.44269504088896341f);
_PS_CONST(cephes_exp_C1, 0.693359375f);
_PS_CONST(cephes_exp_C2, -2.12194440e-4f);

_PS_CONST(cephes_exp_p0, 1.9875691500E-4f);
_PS_CONST(cephes_exp_p1, 1.3981999507E-3f);
_PS_CONST(cephes_exp_p2, 8.3334519073E-3f);
_PS_CONST(cephes_exp_p3, 4.1665795894E-2f);
_PS_CONST(cephes_exp_p4, 1.6666665459E-1f);
_PS_CONST(cephes_exp_p5, 5.0000001201E-1f);

STATIC_INLINE_BEGIN __m128 _mm_exp_cephes_ps(__m128 x)
{
  __m128 tmp = _mm_setzero_ps(), fx;
  __m128i emm0;

  __m128 one = *(__m128*)_ps_1;

  x = _mm_min_ps(x, *(__m128*)_ps_exp_hi);
  x = _mm_max_ps(x, *(__m128*)_ps_exp_lo);

  /* express exp(x) as exp(g + n*log(2)) */
  fx = _mm_mul_ps(x, *(__m128*)_ps_cephes_LOG2EF);
  fx = _mm_add_ps(fx, *(__m128*)_ps_0p5);

  /* how to perform a floorf with SSE: just below */
  emm0 = _mm_cvttps_epi32(fx);
  tmp  = _mm_cvtepi32_ps(emm0);

  /* if greater, substract 1 */
  __m128 mask = _mm_cmpgt_ps(tmp, fx);
  mask = _mm_and_ps(mask, one);
  fx = _mm_sub_ps(tmp, mask);

  tmp = _mm_mul_ps(fx, *(__m128*)_ps_cephes_exp_C1);
  __m128 z = _mm_mul_ps(fx, *(__m128*)_ps_cephes_exp_C2);
  x = _mm_sub_ps(x, tmp);
  x = _mm_sub_ps(x, z);

  z = _mm_mul_ps(x,x);

  __m128 y = *(__m128*)_ps_cephes_exp_p0;
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, *(__m128*)_ps_cephes_exp_p1);
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, *(__m128*)_ps_cephes_exp_p2);
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, *(__m128*)_ps_cephes_exp_p3);
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, *(__m128*)_ps_cephes_exp_p4);
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, *(__m128*)_ps_cephes_exp_p5);
  y = _mm_mul_ps(y, z);
  y = _mm_add_ps(y, x);
  y = _mm_add_ps(y, one);

  /* build 2^n */
  emm0 = _mm_cvttps_epi32(fx);
  emm0 = _mm_add_epi32(emm0, *(__m128i*)_pi32_0x7f);
  emm0 = _mm_slli_epi32(emm0, 23);
  __m128 pow2n = _mm_castsi128_ps(emm0);

  y = _mm_mul_ps(y, pow2n);
  return y;
}



STATIC_INLINE_BEGIN float sin_fast(float a)
{
  // C simulation gives a max absolute error of less than 1.8e-7
  float c0[4] = { 0.0f,            0.5f,
                  1.0f,            0.0f
                };
  float c1[4] = { 0.25f,          -9.0f,
                  0.75f,           0.159154943091f
                };
  float c2[4] = { 24.9808039603f, -24.9808039603f,
                  -60.1458091736f,  60.1458091736f
                };
  float c3[4] = { 85.4537887573f, -85.4537887573f,
                  -64.9393539429f,  64.9393539429f
                };
  float c4[4] = { 19.7392082214f, -19.7392082214f,
                  -1.0f,            1.0f
                };

  // r0.x = sin(a)
  float r0[3], r1[3], r2[3];

  r1[0] = c1[3] * a - c1[0]; // only difference from cos!
  r1[1] = r1[0] - floor(r1[0]); // and extract fraction

  r2[0] = (float) ( r1[1] < c1[0] );            // range check: 0.0 to 0.25
  r2[1] = (float) ( r1[1] >= c1[1] );    // range check: 0.75 to 1.0
  r2[2] = (float) ( r1[1] >= c1[2] );    // range check: 0.75 to 1.0
  r2[1] = r2[0] * c4[2] + r2[1] * c4[3] + r2[2] * c4[2]; // range check: 0.25 to 0.75


  r0[0]    = c0[0] - r1[1];                // range centering
  r0[1]    = c0[1] - r1[1];                // range centering
  r0[2]    = c0[2] - r1[1];                // range centering

  r0[0]    = r0[0] * r0[0];
  r0[1]    = r0[1] * r0[1];
  r0[2]    = r0[2] * r0[2];

  r1[0]    =     c2[0] * r0[0] + c2[2];           // start power series
  r1[0]    =     r1[0] * r0[0] + c3[0];
  r1[0]    =     r1[0] * r0[0] + c3[2];
  r1[0]    =     r1[0] * r0[0] + c4[0];
  r1[0]    =     r1[0] * r0[0] + c4[2];

  r1[1]    =     c2[1] * r0[1] + c2[3];           // start power series
  r1[1]    =     r1[1] * r0[1] + c3[1];
  r1[1]    =     r1[1] * r0[1] + c3[3];
  r1[1]    =     r1[1] * r0[1] + c4[1];
  r1[1]    =     r1[1] * r0[1] + c4[3];

  r1[2]    =     c2[0] * r0[2] + c2[2];           // start power series
  r1[2]    =     r1[2] * r0[2] + c3[0];
  r1[2]    =     r1[2] * r0[2] + c3[2];
  r1[2]    =     r1[2] * r0[2] + c4[0];
  r1[2]    =     r1[2] * r0[2] + c4[2];

  r0[0] = -r1[0] * r2[0] - r1[1]*r2[1] - r1[2]*r2[2];

  return r0[0];
}
/*

  float sin_fast(float a)
  {
  // C simulation gives a max absolute error of less than 1.8e-7
  float4 c0 = float4( 0.0,            0.5,
  1.0,            0.0            );
  float4 c1 = float4( 0.25,          -9.0,
  0.75,           0.159154943091 );
  float4 c2 = float4( 24.9808039603, -24.9808039603,
  -60.1458091736,  60.1458091736  );
  float4 c3 = float4( 85.4537887573, -85.4537887573,
  -64.9393539429,  64.9393539429  );
  float4 c4 = float4( 19.7392082214, -19.7392082214,
  -1.0,            1.0            );

  // r0.x = sin(a)
  float3 r0, r1, r2;

  r1.x  = c1.w * a - c1.x;                // only difference from cos!
  r1.y  = frac( r1.x );                   // and extract fraction
  r2.x  = (float) ( r1.y < c1.x );        // range check: 0.0 to 0.25
  r2.yz = (float2) ( r1.yy >= c1.yz );    // range check: 0.75 to 1.0
  r2.y  = dot( r2, c4.zwz );              // range check: 0.25 to 0.75
  r0    = c0.xyz - r1.yyy;                // range centering
  r0    = r0 * r0;
  r1    = c2.xyx * r0 + c2.zwz;           // start power series
  r1    =     r1 * r0 + c3.xyx;
  r1    =     r1 * r0 + c3.zwz;
  r1    =     r1 * r0 + c4.xyx;
  r1    =     r1 * r0 + c4.zwz;
  r0.x  = dot( r1, -r2 );                 // range extract

  return r0.x;
  }


  cos(float a)
  {
  // C simulation gives a max absolute error of less than 1.8e-7
  const float4 c0 = float4( 0.0,            0.5,
  1.0,            0.0            );
  const float4 c1 = float4( 0.25,          -9.0,
  0.75,           0.159154943091 );
  const float4 c2 = float4( 24.9808039603, -24.9808039603,
  -60.1458091736,  60.1458091736  );
  const float4 c3 = float4( 85.4537887573, -85.4537887573,
  -64.9393539429,  64.9393539429  );
  const float4 c4 = float4( 19.7392082214, -19.7392082214,
  -1.0,            1.0            );

  // r0.x = cos(a)
  float3 r0, r1, r2;

  r1.x  = c1.w * a;                       // normalize input
  r1.y  = frac( r1.x );                   // and extract fraction
  r2.x  = (float) ( r1.y < c1.x );        // range check: 0.0 to 0.25
  r2.yz = (float2) ( r1.yy >= c1.yz );    // range check: 0.75 to 1.0
  r2.y  = dot( r2, c4.zwz );              // range check: 0.25 to 0.75
  r0    = c0.xyz - r1.yyy;                // range centering
  r0    = r0 * r0;
  r1    = c2.xyx * r0 + c2.zwz;           // start power series
  r1    =     r1 * r0 + c3.xyx;
  r1    =     r1 * r0 + c3.zwz;
  r1    =     r1 * r0 + c4.xyx;
  r1    =     r1 * r0 + c4.zwz;
  r0.x  = dot( r1, -r2 );                 // range extract

  return r0.x;

  tan(float a) {
  return sin(a) / cos(a);
  }


*/

STATIC_INLINE_BEGIN __m128i ilogbp1(__m128 d)
{
  __m128 m = _mm_cmplt_ps(d, _mm_set1_ps(5.421010862427522E-20f));
  d = _mm_sel_ps(d, _mm_mul_ps(_mm_set1_ps(1.8446744073709552E19f), d), m);
  __m128i q = _mm_and_si128(_mm_srli_epi32(_mm_castps_si128(d), 23),
                            _mm_set1_epi32(0xff));
  q = _mm_sub_epi32(q, _mm_sel_epi32(_mm_set1_epi32(0x7e), _mm_set1_epi32(64 + 0x7e), _mm_castps_si128(m)));
  return q;
}

/**
 * Multiplies the packed singles x by the number 2 raised to the q power
 *
 * @param x
 * @param q
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128 _mm_ldexpf(__m128 x, __m128i q)
{
  __m128 u;
  __m128i m = _mm_srai_epi32(q, 31);
  // 31 is size without sign-bit, 23 is size of mantissa

  m = _mm_slli_epi32(_mm_sub_epi32(_mm_srai_epi32(_mm_add_epi32(m, q), 6), m), 4);
  q = _mm_sub_epi32(q, _mm_slli_epi32(m, 2));
  m = _mm_add_epi32(m, _mm_set1_epi32(0x7f));

  m = _mm_and_si128(_mm_cmpgt_epi32(m, _mm_set1_epi32(0)), m);
  __m128i n = _mm_cmpgt_epi32(m, _mm_set1_epi32(0xff));
  m = _mm_or_si128(_mm_andnot_si128(n, m), _mm_and_si128(n, _mm_set1_epi32(0xff)));

  u = _mm_castsi128_ps(_mm_slli_epi32(m, 23));
  x = _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(_mm_mul_ps(x, u), u), u), u);
  u = _mm_castsi128_ps((_mm_slli_epi32(_mm_add_epi32(q, _mm_set1_epi32(0x7f)), 23)));
  return _mm_mul_ps(x, u);
}

// TODO: This is wrong
STATIC_INLINE_BEGIN __m256d _mm256_ldexpd(__m256d x, __m128i q)
{
  SPS_UNREFERENCED_PARAMETERS(x, q);
  // SVML provides _mm256_exp_pd
  return _mm256_setzero_pd();
#if 0
  __m128 u;
  __m128i m = _mm_srai_epi32(q, 31);
  // 31 is size without sign-bit, 23 is size of mantissa

  m = _mm_slli_epi32(_mm_sub_epi32(_mm_srai_epi32(_mm_add_epi32(m, q), 6), m), 4);
  q = _mm_sub_epi32(q, _mm_slli_epi32(m, 2));
  m = _mm_add_epi32(m, _mm_set1_epi32(0x7f));

  m = _mm_and_si128(_mm_cmpgt_epi32(m, _mm_set1_epi32(0)), m);
  __m128i n = _mm_cmpgt_epi32(m, _mm_set1_epi32(0xff));
  m = _mm_or_si128(_mm_andnot_si128(n, m), _mm_and_si128(n, _mm_set1_epi32(0xff)));

  // Sign of mantissa is 23
  u = _mm_castsi128_ps(_mm_slli_epi32(m, 23));
  x = _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(_mm_mul_ps(x, u), u), u), u);
  u = _mm_castsi128_ps((_mm_slli_epi32(_mm_add_epi32(q, _mm_set1_epi32(0x7f)), 23)));
  return _mm_mul_ps(x, u);
#endif
}




/**
 * Cubic root
 *
 * @param d
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128 _mm_cbrtf_ps(__m128 d)
{
  __m128 x, y, q = _mm_set1_ps(1.0), t;
  __m128i e, qu, re;

  e = ilogbp1(_mm_fabs_ps(d));
  d = _mm_ldexpf(d, _mm_neg_epi32(e));

  t = _mm_add_ps(_mm_cvtepi32_ps(e), _mm_set1_ps(6144));

  qu = _mm_cvttps_epi32(_mm_mul_ps(t, _mm_set1_ps(1.0f/3.0f)));
  re = _mm_cvttps_epi32(_mm_sub_ps(t, _mm_mul_ps(_mm_cvtepi32_ps(qu), _mm_set1_ps(3))));

  q = _mm_sel_ps(q, _mm_set1_ps(1.2599210498948731647672106f), _mm_castsi128_ps(_mm_cmpeq_epi32(re,_mm_set1_epi32(1))));
  q = _mm_sel_ps(q, _mm_set1_ps(1.5874010519681994747517056f), _mm_castsi128_ps(_mm_cmpeq_epi32(re,_mm_set1_epi32(2))));

  q = _mm_ldexpf(q, _mm_sub_epi32(qu, _mm_set1_epi32(2048)));

  q = _mm_mulsign_ps(q, d);
  d = _mm_fabs_ps(d);

  x = _mm_set1_ps(-0.601564466953277587890625f);
  x = _mm_madd_ps(x, d, _mm_set1_ps(2.8208892345428466796875f));
  x = _mm_madd_ps(x, d, _mm_set1_ps(-5.532182216644287109375f));
  x = _mm_madd_ps(x, d, _mm_set1_ps(5.898262500762939453125f));
  x = _mm_madd_ps(x, d, _mm_set1_ps(-3.8095417022705078125f));
  x = _mm_madd_ps(x, d, _mm_set1_ps(2.2241256237030029296875f));

  y = _mm_mul_ps(_mm_mul_ps(d, x), x);
  y = _mm_mul_ps(_mm_sub_ps(y, _mm_mul_ps(_mm_mul_ps(_mm_set1_ps(2.0f / 3.0f), y), _mm_madd_ps(y, x, _mm_set1_ps(-1.0f)))), q);

  return y;
}

#ifdef _MSC_VER
# pragma warning( push )
# pragma warning( disable : 4756 )
# pragma warning( disable : 4056 )
#endif
/**
 * Logarithm of packed singles
 *
 * @param d
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128 _mm_log_ps(__m128 d)
{
  __m128 x, x2, t, m;
  __m128i e;

  e = ilogbp1(x = _mm_mul_ps(d, _mm_set1_ps(0.7071f)));
  m = _mm_ldexpf(d, _mm_neg_epi32(e));
  d = x;

  x = _mm_div_ps(_mm_add_ps(_mm_set1_ps(-1.0f), m), _mm_add_ps(_mm_set1_ps(1.0f), m));
  x2 = _mm_mul_ps(x, x);

  t = _mm_set1_ps(0.2371599674224853515625f);
  t = _mm_madd_ps(t, x2, _mm_set1_ps(0.285279005765914916992188f));
  t = _mm_madd_ps(t, x2, _mm_set1_ps(0.400005519390106201171875f));
  t = _mm_madd_ps(t, x2, _mm_set1_ps(0.666666567325592041015625f));
  t = _mm_madd_ps(t, x2, _mm_set1_ps(2.0f));

  x = _mm_madd_ps(x, t, _mm_mul_ps(_mm_set1_ps(0.693147180559945286226764f), _mm_cvtepi32_ps(e)));

  x = _mm_sel_ps(x, _mm_set1_ps(INFINITYf), _mm_is_infinity(d));
  x = _mm_or_ps(_mm_cmpgt_ps(_mm_set1_ps(0),d),x);
  x = _mm_sel_ps(x,_mm_set1_ps(-INFINITYf),_mm_cmpeq_ps(d,_mm_set1_ps(0)));

  return x;
}
#ifdef _MSC_VER
# pragma warning( pop )
#endif

#ifdef __cplusplus
}
#endif

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
