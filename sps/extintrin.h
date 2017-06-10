/**
 * @file   extintrin.h
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Tue May  5 20:13:58 2015
 *
 * @brief  Extensions to Intel's intrinsics
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

// Using AVX2, shifts can be made using _mm_sllv_epi32 (sse_halfs)

#include <sps/cenv.h>
#include <float.h>
#include <stdint.h>

#include <sps/config.h>

#ifdef HAVE_CONFIG_H
#else
// Check is __AVX__ __AVX2__ __SSE4_1__ __SSE3__ __SSE2__ __SSE__
# pragma message("Warning SSE2,SSE3,SSE4(1,2) and AVX is enabled")
# define HAVE_EMMINTRIN_H 1
# define HAVE_PMMINTRIN_H 1
# define HAVE_SMMINTRIN_H 1
# define HAVE_NMMINTRIN_H 1
# define HAVE_IMMINTRIN_H 1
#endif

#ifdef HAVE_EMMINTRIN_H
# include <emmintrin.h> // SSE2
#endif

#ifdef HAVE_PMMINTRIN_H
# include <tmmintrin.h> // SSE3
#endif

#ifdef HAVE_SMMINTRIN_H
# include <smmintrin.h> // SSE4.1
#endif

#ifdef HAVE_NMMINTRIN_H
# include <nmmintrin.h> // SSE4.2
#endif

#ifdef HAVE_IMMINTRIN_H
# include <immintrin.h> // AVX
#endif

#ifndef WIN32
# ifdef HAVE_ZMMINTRIN_H
#  include <zmmintrin.h> // AVX2
# endif
#else
#endif

#ifdef HAVE_STDINT_H
# include <stdint.h>
#endif

#ifdef WIN32
// Bulldozer AMD
// # include <x86intrin.h>
#endif


#ifdef __cplusplus
# include <iostream>
extern "C" {
#endif

/**
 * \defgroup Other constants
 * @{
 */
const __m128 _m_half_ps         = _mm_set1_ps(0.5f);
const __m128 _m_m_half_ps       = _mm_set1_ps(-0.5f);
const __m128 _m_one_ps          = _mm_set1_ps(1.0f);
const __m128 _m_sixth_ps        = _mm_set1_ps(0.166666666f);
// Consider moving these elsewhere (gives warning if non-static)
static const __m128 _m_m_one_ps = _mm_set1_ps(-1.0f);
static const __m128 _m_eps_ps   = _mm_set1_ps(FLT_EPSILON);


const __m128d _m_half_pd         = _mm_set1_pd(0.5);
const __m128d _m_one_pd          = _mm_set1_pd(1.0);

static const __m128d _m_eps_pd   = _mm_set1_pd(DBL_EPSILON);

static const __m256d _m256_m_one_pd = _mm256_set1_pd(-1.0);
static const __m256d _m256_one_pd = _mm256_set1_pd(1.0);
static const __m256d _m256_half_pd = _mm256_set1_pd(0.5);

/**@}*/

/**
 * \defgroup Constants used by functions
 * @{
 */
const ALIGN16_BEGIN int _clear_signmask[4] ALIGN16_END =
{0x7FFFFFFFL,0x7FFFFFFFL,0x7FFFFFFFL,0x7FFFFFFFL};                       ///< clear sign mask for ps

const ALIGN16_BEGIN long long int _clear_signmaskd[2] ALIGN16_END =
{0x7FFFFFFFFFFFFFFF,0x7FFFFFFFFFFFFFFF};                                 ///< clear sign mask for pd

const ALIGN32_BEGIN long long int _clear_signmask_256[4] ALIGN32_END =
{0x7FFFFFFFFFFFFFFF,0x7FFFFFFFFFFFFFFF,0x7FFFFFFFFFFFFFFF,0x7FFFFFFFFFFFFFFF}; ///< clear sign mask for ps


#ifdef _WIN32
const __m128 signmask = _mm_castsi128_ps(_mm_set1_epi32(0x80000000L));     ///< signmask for ps (packed singles)
#elif defined(__GNUC__)
const __m128 signmask = _mm_castsi128_ps(_mm_set1_epi32(-2147483648));     ///< signmask for ps (packed singles)
#endif

const __m128 mantmask = _mm_castsi128_ps(_mm_set1_epi32(0x7FFFFFFFL));     ///< mantissa + exp mask for ps

const ALIGN16_BEGIN int _neg_signmask[4] ALIGN16_END =
{(int)0x80000000L,(int)0x80000000L,(int)0x80000000L,(int)0x80000000L};   ///< negative sign mask for ps

const ALIGN16_BEGIN long long int _neg_signmaskd[2] ALIGN16_END =
{(long long int)0x8000000000000000,(long long int)0x8000000000000000};   ///< negative sign mask for pd

const __m128i negmask = _mm_set1_epi32(-1);                              ///< negate mask 0xFFFFFFFF for epi32


const __m256i negmask_256 = _mm256_set_epi64x(-1,-1,-1,-1);              ///< negate mask 0xFFFFFFFFFFFFFFFF for epi32


/**@}*/

#ifdef _MSC_VER
typedef union __declspec(intrin_type) _CRT_ALIGN(16) v4f {
  __m128 v;
  float               f32[4];
  uint64_t            uint64[2];
  int8_t              int8[16];
  int16_t             int16[8];
  int32_t             int32[4];
  int64_t             int64[2];
  uint8_t             uint8[16];
  uint16_t            uint16[8];
  uint32_t            uint32[4];
  v4f() {}
  v4f(const __m128& _v) {
    v = _v;
  }
} v4f;

typedef union __declspec(intrin_type) _CRT_ALIGN(16) v4i {
  __m128i v;
  float               f32[4];
  uint64_t            uint64[2];
  int8_t              int8[16];
  int16_t             int16[8];
  int32_t             int32[4];
  int64_t             int64[2];
  uint8_t             uint8[16];
  uint16_t            uint16[8];
  uint32_t            uint32[4];
  v4i() {}
  v4i(const __m128i& _v) {
    v = _v;
  }
} v4i;

typedef union __declspec(intrin_type) _CRT_ALIGN(32) v4d {
  __m256d v;
  double               f64[4];
  v4d() {}
  v4d(const __m256d& _v) {
    v = _v;
  }
} v4d;

#else
typedef union v4f {
  __m128 v;
  float               f32[4];
  uint64_t            uint64[2];
  int8_t              int8[16];
  int16_t             int16[8];
  int32_t             int32[4];
  int64_t             int64[2];
  uint8_t             uint8[16];
  uint16_t            uint16[8];
  uint32_t            uint32[4];
  v4f() {}
  v4f(const __m128& _v) {
    v = _v;
  }
} v4f __attribute__ ((aligned (16)));

typedef union v4i {
  __m128i v;
  float               f32[4];
  uint64_t            uint64[2];
  int8_t              int8[16];
  int16_t             int16[8];
  int32_t             int32[4];
  int64_t             int64[2];
  uint8_t             uint8[16];
  uint16_t            uint16[8];
  uint32_t            uint32[4];
  v4i() {}
  v4i(const __m128i& _v) {
    v = _v;
  }
} v4i __attribute__ ((aligned (16)));

typedef union v4d {
  __m256d v;
  double               f64[4];
  v4d() {}
  v4d(const __m256d& _v) {
    v = _v;
  }
}  v4d __attribute__ ((aligned (32)));

#endif

STATIC_INLINE_BEGIN __m256d _mm256_movehl_pd(__m256d a, __m256d b)
{
  return _mm256_castps_pd(_mm256_permute2f128_ps(_mm256_castpd_ps(a), _mm256_castpd_ps(b), 0x11 /*0x21*/));
}


/**
 * Multiply and accumulate
 *
 * @param x
 * @param y
 * @param z
 *
 * @return (x * y) + z
 */
STATIC_INLINE_BEGIN __m128  _mm_madd_ps(__m128 a, __m128 b, __m128 c) STATIC_INLINE_END;

#ifndef _MSC_VER
/**
 * 64-bit full multiplication and shift
 *
 * @param a
 * @param b
 * @param s
 *
 * @return
 */
STATIC_INLINE_BEGIN int64_t mulshift (int64_t a, int64_t b, int s)
{
  int64_t res;
  __asm__ volatile (                // rax = a, rdx = b, ecx = s
    "imulq %%rdx;\n\t"              // rdx:rax = rax * rdx
    "shrdq %%cl, %%rdx, %%rax;\n\t" // rax = int64_t (rdx:rax >> s)
    "movq  %%rax, %0;"              // res = rax
    : "=rm" (res)
    : "a"(a), "d"(b), "c"(s));
  return res;
}
#endif

/**
 * Complex multiply, 2.5 instructions per multiplication instead of 6
 *
 * @param a
 * @param b
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128 _mm_mulcmplx_ps(const __m128 &a, const __m128 &b) STATIC_INLINE_END;

/**
 * Test x==y for each packed single in a and b
 *
 * @param a
 * @param b
 *
 * @return 1 if any of the packed floats are equal, else 0
 */
STATIC_INLINE_BEGIN int _mm_any_eq( __m128 a, __m128 b ) STATIC_INLINE_END;

/**
 * More accurate reciprocal square root using Newton-Rhapson
 *
 * @param a
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128 _mm_rsqrt_nr_ps(const __m128& a) STATIC_INLINE_END;

STATIC_INLINE_BEGIN __m128 _mm_rcp_nr_ps(const __m128& x) STATIC_INLINE_END;

STATIC_INLINE_BEGIN __m128 _mm_rcp_nr_ps(const __m128& x)
{
  __m128 res = _mm_rcp_ps(x);
  __m128 muls = _mm_mul_ps(x, _mm_mul_ps(res, res));
  return res =  _mm_sub_ps(_mm_add_ps(res, res), muls);
}

// AVX-512 (HACK)
#if defined(_MSC_VER) && (_MSC_VER <= 2000)
STATIC_INLINE_BEGIN __m256d _mm256_rcp14_pd(__m256d x)
{
  __m128 x1 = _mm256_cvtpd_ps(x);
  x1 = _mm_rcp_ps(x1);
  return _mm256_cvtps_pd(x1);
}
#endif

// SVML (Intel) provides _mm256_exp_pd in ia32intrin.h

#ifndef INTEL_COMPILER
STATIC_INLINE_BEGIN double
_mm256_cvtsd_f64(__m256d d)
{
  // TODO: Keep input argument in register
  v4d _d;
  _d.v = d;
  return _d.f64[0];
}
#endif

/**
 * Absolute values of packed singles.
 *
 * @param x
 *
 * @return
 */
STATIC_INLINE_BEGIN  __m128 _mm_fabs_ps(__m128 x) STATIC_INLINE_END;

/**
 * Absolute values of packed doubles.
 *
 * @param x
 *
 * @return
 */
STATIC_INLINE_BEGIN  __m128d _mm_fabs_pd(__m128d x) STATIC_INLINE_END;

/**
 * Absolute values of packed doubles.
 *
 * @param x
 *
 * @return
 */
STATIC_INLINE_BEGIN  __m256d _mm256_fabs_pd(__m256d x) STATIC_INLINE_END;

/**
 * Floating point modulo of packed singles
 *
 * @param a
 * @param aDiv
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128 _mm_fmod_ps(const __m128& a, const __m128& aDiv) STATIC_INLINE_END;

STATIC_INLINE_BEGIN __m128 _mm_fmod_ps(const __m128& a, const __m128& aDiv)
{
  __m128 c = _mm_div_ps(a,aDiv);
  __m128i i = _mm_cvttps_epi32(c);
  __m128 cTrunc = _mm_cvtepi32_ps(i);
  __m128 base = _mm_mul_ps(cTrunc, aDiv);
  __m128 r = _mm_sub_ps(a, base);
  return r;
}

#ifndef HAVE_SMMINTRIN_H
//SSE2: multiply 8 16-bit integers and return full 32-bit result in high and low result
STATIC_INLINE_BEGIN void _mm_mul_epi16_full(__m128i &hiResult, __m128i &loResult, const __m128i a, const __m128i b)
{
  __m128i hi  = _mm_mulhi_epi16( a, b );    // (a7*b7[16:31],a6*b6[16:31],a5*b5[16:31],a4*b4[16:31],a3*b3[16:31],a2*b2[16:31],a1*b1[16:31],a0*b0[16:31])
  __m128i low = _mm_mullo_epi16( a, b );    // (a7*b7[0:15] ,a6*b6[0:15] ,a5*b5[0:15] ,a4*b4[0:15],a3*b3[0:15] ,a2*b2[0:15] ,a1*b1[0:15] ,a0*b0[0:15])
  loResult = _mm_unpacklo_epi16( hi, low ); // (a3*b3[0:15],a3*b3[16:31],a2*b2[0:15],a2*b2[16:31],a1*b1[0:15],a1*b1[16:31],a0*b0[0:15],a0*b0[16:31])
  hiResult = _mm_unpackhi_epi16( hi, low ); // (a7*b7[0:15],a7*b7[16:31],a6*b6[0:15],a6*b6[16:31],a5*b5[0:15],a5*b5[16:31],a4*b4[0:15],a4*b4[16:31])
}
#endif

STATIC_INLINE_BEGIN void _mm_mul_epi32_full(__m128i &loResult, __m128i &hiResult, __m128i a, __m128i b)
{
  __m128i _a        = _mm_castps_si128(_mm_permute_ps(_mm_castsi128_ps(a),0xB1));
  __m128i _b        = _mm_castps_si128(_mm_permute_ps(_mm_castsi128_ps(b),0xB1));
  __m128i c1        = _mm_mul_epi32(a,b);   // 0 and 2
  __m128i c2        = _mm_mul_epi32(_a,_b); // 1 and 3
  loResult          = _mm_unpacklo_epi64(c1,c2);
  hiResult          = _mm_unpackhi_epi64(c1,c2);
}


STATIC_INLINE_BEGIN __m128 _mm_mulcmplx_ps(const __m128 &a, const __m128 &b)
{

  // 2.5 instructions instead of 6
  const __m128 xmm4 = _mm_set_ps(1.0f,-1.0f,1.0f,-1.0f);

  __m128 xmm2 = _mm_shuffle_ps(b,b, 0xA0 ); // c3a c3a c4a c4a
  __m128 xmm1 = _mm_shuffle_ps(b,b, 0xF5);
  __m128 xmm3 = _mm_shuffle_ps(a,a, 0xB1);
  __m128 xmm0 = _mm_mul_ps(a,xmm2);
  xmm3 = _mm_mul_ps(xmm1,xmm3);
  xmm3 = _mm_mul_ps(xmm4,xmm3); // Replaced with FMA
  xmm0 = _mm_add_ps(xmm0,xmm3);
  return xmm0;
}

/**
 * Horizontal AND values, 32-bit at a time. Result is placed in the
 * least significant 32-bit integer.
 *
 * @param x
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128i _mm_hand_epi32(__m128i x)
{
  __m128i r = _mm_and_si128(x,_mm_srli_si128(x,8));
  r = _mm_and_si128(r,_mm_srli_si128(r,4));
  return r;
}

/**
 * Horizontal OR values, 32-bit at a time. Result is placed in the
 * least significant 32-bit integer.
 *
 * @param x
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128i _mm_hor_epi32(__m128i x)
{
  __m128i r = _mm_or_si128(x,_mm_srli_si128(x,8));
  r = _mm_or_si128(r,_mm_srli_si128(r,4));
  return r;
}


#if defined(NDEBUG) && !defined(_MSC_VER)
/**
 * Right(logic)-shift packed singles
 *
 * @param x
 * @param imm
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128 _mm_srli_ps(__m128 x, const int imm)
{
  return _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(x),imm));
}
#else
// 8-bit const immediates cannot be considered immediate when debugging
# define _mm_srli_ps(x,imm)                                          \
           _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(x),imm))
#endif

#if defined(NDEBUG) && !defined(_MSC_VER) && !defined(__CYGWIN__)
/**
 * Left-shift of packed singles
 *
 * @param x
 * @param imm Number of bits (a constant immediate)
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128 _mm_slli_ps(__m128 x, const int imm)
{
  return _mm_castsi128_ps(_mm_slli_si128(_mm_castps_si128(x),imm));
}
#else
// 8-bit const immediates cannot be considered immediate when debugging
# define _mm_slli_ps(x,imm)                                          \
           _mm_castsi128_ps(_mm_slli_si128(_mm_castps_si128(x),imm))
#endif


#if HAVE_SMMINTRIN_H
STATIC_INLINE_BEGIN int _mm_dp1_epi32(__m128i x, __m128i y)
{

  __m128i s = _mm_mullo_epi32(x,y);
  s = _mm_add_epi32(s, _mm_srli_si128(s, 8));
  s = _mm_add_epi32(s, _mm_srli_si128(s, 4));
  int sum = _mm_cvtsi128_si32(s);
  return sum;
}

STATIC_INLINE_BEGIN __m128i _mm_dp_epi32(__m128i x, __m128i y)
{

  __m128i s = _mm_mullo_epi32(x,y);
  s = _mm_add_epi32(s, _mm_srli_si128(s, 8));
  s = _mm_add_epi32(s, _mm_srli_si128(s, 4));
  return s;
}
#endif

STATIC_INLINE_BEGIN  __m128d _mm256_hsum_pd(__m256d x1, __m256d x2)
{
  // calculate 4 two-element horizontal sums:
  // lower 64 bits contain x1[0] + x1[1]
  // next 64 bits contain x2[0] + x2[1]
  // next 64 bits contain x1[2] + x1[3]
  // next 64 bits contain x2[2] + x2[3]
  __m256d sum = _mm256_hadd_pd(x1, x2);
  // extract upper 128 bits of result
  __m128d sum_high = _mm256_extractf128_pd(sum, 1);
  // add upper 128 bits of sum to its lower 128 bits
  __m128d result = _mm_add_pd(sum_high, _mm256_castpd256_pd128(sum));
  // lower 64 bits of result contain the sum of x1[0], x1[1], x1[2], x1[3]
  // upper 64 bits of result contain the sum of x2[0], x2[1], x2[2], x2[3]
  return result;
}


STATIC_INLINE_BEGIN __m256d _mm256_dp_pd(const __m256d& x, const __m256d& y, const int mask) STATIC_INLINE_END;

STATIC_INLINE_BEGIN __m256d _mm256_dp_pd(const __m256d& x, const __m256d& y, const int mask)
{
  // Without AVX-512 and __m256_cmp_epu32 a few casts are needed
  const __m256i smask      = _mm256_set_epi64x(0x80,0x40,0x20,0x10);
  const __m256i omask      = _mm256_set_epi64x(0x08,0x04,0x02,0x01);

  const __m256i bum = _mm256_set1_epi64x(mask);

#ifdef HAVE_ZMMINTRIN_H
  const __m256d selectMask = _mm256_cmp_pd(
                               _mm256_castsi256_pd(_mm256_and_si256(smask, bum)),
                               _mm256_castsi256_pd(smask),
                               _CMP_EQ_UQ);
  const __m256d outputMask = _mm256_cmp_pd(
                               _mm256_castsi256_pd(_mm256_and_si256(omask, bum)),
                               _mm256_castsi256_pd(omask),
                               _CMP_EQ_UQ);
#else
  // Without AVX-2 more casts are needed
  const __m256d selectMask = _mm256_cmp_pd(
                               _mm256_castps_pd(_mm256_and_ps(_mm256_castsi256_ps(smask), _mm256_castsi256_ps(bum))),
                               _mm256_castsi256_pd(smask),
                               _CMP_EQ_UQ);
  const __m256d outputMask = _mm256_cmp_pd(
                               _mm256_castps_pd(_mm256_and_ps(_mm256_castsi256_ps(omask), _mm256_castsi256_ps(bum))),
                               _mm256_castsi256_pd(omask),
                               _CMP_EQ_UQ);
#endif

  __m256d xy = _mm256_mul_pd( x, y );
  xy = _mm256_and_pd(xy,selectMask);

  __m256d temp = _mm256_hadd_pd( xy, xy );

  // Without AVX-512 and _mm256_permute2f128_pd, we need a few casts
  __m256d temp_high_low = _mm256_castps_pd(_mm256_permute2f128_ps(_mm256_castpd_ps(temp), _mm256_castpd_ps(temp), 0x21));
  __m256d r3   = _mm256_add_pd(temp, temp_high_low);
  return _mm256_and_pd(r3,outputMask);

}

#if HAVE_SMMINTRIN_H
/*
   SSE4.1 provides _mm_dp_ps.

   Note 4 (4-vector) dot products using _mm_mul_ps and _mm_add_ps
   takes 5 clock cycles (4 clock cycles using FMA). Using
   _mm_dp_ps, it takes 4 clock cycles, 1 clock cycle for each dot
   product.
*/
#else
/**
 * Dot-product of packed singles
 *
 * @param __X
 * @param __Y
 * @param __M
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128 _mm_dp_ps(__m128 __X, __m128 __Y, const int __M)
{

# if HAVE_PMMINTRIN_H

  const __m128 smask      = _mm_set_epi32(0x80,0x40,0x20,0x10);
  const __m128 omask      = _mm_set_epi32(0x08,0x04,0x02,0x01);
  // TODO: Need comparison to form mask
  const __m128 selectMask = _mm_castsi128_ps(_mm_and_si128(smask,_mm_set1_epi32(__M)));
  const __m128 outputMask = _mm_castsi128_ps(_mm_and_si128(omask,_mm_set1_epi32(__M)));

  const r1 = _mm_mul_ps(__X, __Y);
  const r1 = _mm_and_ps(r1,selectMask);
  const r2 = _mm_hadd_ps(r1, r1);
  const r3 = _mm_hadd_ps(r2, r2);
  return _mm_and_ps(r3,outputMask);

# elif HAVE_EMMINTRIN_H

  const __m128 smask      = _mm_set_epi32(0x80,0x40,0x20,0x10);
  const __m128 omask      = _mm_set_epi32(0x08,0x04,0x02,0x01);
  // TODO: Need comparison to form mask
  const __m128 selectMask = _mm_castsi128_ps(_mm_and_si128(smask,_mm_set1_epi32(__M)));
  const __m128 outputMask = _mm_castsi128_ps(_mm_and_si128(omask,_mm_set1_epi32(__M)));

  const __m128 mult  = _mm_mul_ps(a, b);
  const __m128 shuf1 = _mm_shuffle_ps(mult, mult, _MM_SHUFFLE(0, 3, 2, 1));
  const __m128 shuf2 = _mm_shuffle_ps(mult, mult, _MM_SHUFFLE(1, 0, 3, 2));
  const __m128 shuf3 = _mm_shuffle_ps(mult, mult, _MM_SHUFFLE(2, 1, 0, 3));

  return _mm_add_ps(_mm_add_ps(_mm_add_ps(mult, shuf1), shuf2), shuf3);
# else
#  error Enable SSE2, SSE3 or SSE4
# endif
}
#endif

/*

__m256d xy = _mm256_mul_pd( x, y );
__m256d temp = _mm256_hadd_pd( xy, xy );
__m128d hi128 = _mm256_extractf128_pd( temp, 1 );
__m128d dotproduct = _mm_add_pd( (__m128d)temp, hi128 );

Edit:
After an idea of Norbert P. I extended this version to do 4 dot products at one time.

__m256d xy0 = _mm256_mul_pd( x[0], y[0] );
__m256d xy1 = _mm256_mul_pd( x[1], y[1] );
__m256d xy2 = _mm256_mul_pd( x[2], y[2] );
__m256d xy3 = _mm256_mul_pd( x[3], y[3] );

// low to high: xy00+xy01 xy10+xy11 xy02+xy03 xy12+xy13
__m256d temp01 = _mm256_hadd_pd( xy0, xy1 );

// low to high: xy20+xy21 xy30+xy31 xy22+xy23 xy32+xy33
__m256d temp23 = _mm256_hadd_pd( xy2, xy3 );

// low to high: xy02+xy03 xy12+xy13 xy20+xy21 xy30+xy31
__m256d swapped = _mm256_permute2f128_pd( temp01, temp23, 0x21 );

// low to high: xy00+xy01 xy10+xy11 xy22+xy23 xy32+xy33
__m256d blended = _mm256_blend_pd(temp01, temp23, 0b1100);

__m256d dotproduct = _mm256_add_pd( swapped, blended );






 */


/**
 * Absolute value of packed singles
 *
 * @param x
 *
 * @return
 */
STATIC_INLINE_BEGIN  __m128 _mm_fabs_ps(__m128 x)
{
  // Same as _mm_set1_ps(-0.f);
  return _mm_and_ps(x,_mm_load_ps((float *) _clear_signmask));
}

STATIC_INLINE_BEGIN  __m128d _mm_fabs_pd(__m128d x)
{
  return _mm_and_pd(x,_mm_load_pd((double *) _clear_signmaskd));
}

STATIC_INLINE_BEGIN  __m256d _mm256_fabs_pd(__m256d x)
{
  return _mm256_and_pd(x,_mm256_load_pd((double *) _clear_signmask_256));
}

/**
 * Negate packed singles.
 *
 * @param x
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128 _mm_neg_ps(__m128 x)
{
  return _mm_xor_ps(x,_mm_load_ps((float *) _neg_signmask));
}

/**
 * Negate packed doubles.
 *
 * @param x
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128d _mm_neg_pd(__m128d x)
{
  return _mm_xor_pd(x,_mm_load_pd((double *) _neg_signmaskd));
}

/**
 * Negate packed integers.
 *
 * @param x
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128i _mm_neg_epi32(__m128i x)
{
  return _mm_sub_epi32(_mm_setzero_si128(), x);
  // or (x ^ 0xFFFFFFFF) + 1
}

/**
 * Select packed floating point result based on mask.
 *
 * @param a
 * @param b
 * @param mask
 *
 * @return If mask is false, a is returned, otherwise b is returned.
 */
STATIC_INLINE_BEGIN __m128 _mm_sel_ps(__m128 a, __m128 b, __m128 mask )
{
#ifdef HAVE_SMMINTRIN_H
  return _mm_blendv_ps(a,b,mask);
#else
  b = _mm_and_ps( b, mask );
  a = _mm_andnot_ps( mask, a );
  return _mm_or_ps( a, b );
#endif
}

/**
 * Select packed double result based on mask.
 *
 * @param a
 * @param b
 * @param mask
 *
 * @return If mask is false, a is returned, otherwise b is returned.
 */
STATIC_INLINE_BEGIN __m128d _mm_sel_pd(__m128d a, __m128d b, __m128d mask )
{
#ifdef HAVE_SMMINTRIN_H
  return _mm_blendv_pd(a,b,mask);
#else
  b = _mm_and_pd( b, mask );
  a = _mm_andnot_pd( mask, a );
  return _mm_or_pd( a, b );
#endif
}

STATIC_INLINE_BEGIN __m256d _mm256_sel_pd(__m256d a,
    __m256d b,
    __m256d mask )
{
#ifdef HAVE_ZMMINTRIN_H
  return _mm256_blendv_pd(a,b,mask);
#else
  a = _mm256_andnot_pd( mask, a);
  b = _mm256_and_pd( b, mask );
  return _mm256_or_pd( a, b );
#endif
}



#if HAVE_SMMINTRIN_H
// We have _mm_blendv_ps
#else
/**
 * Select packed singles based on a mask
 *
 * @param a
 * @param b
 * @param mask
 *
 * @return If mask is false, a is returned, otherwise b is returned.
 */
STATIC_INLINE_BEGIN __m128 _mm_blendv_ps(__m128 a, __m128 b, __m128 mask )
{
  b = _mm_and_ps( b, mask );
  a = _mm_andnot_ps( mask, a );
  return _mm_or_ps( a, b );
}
#endif

/**
 * Select integer result based on mask.
 *
 * @param a
 * @param b
 * @param mask
 *
 * @return If mask is false, a is returned, otherwise b is returned.
 */
STATIC_INLINE_BEGIN __m128i _mm_sel_epi32(__m128i a, __m128i b, __m128i mask )
{
  // Using _mm_blendv_ps is slower due to register usage
  b = _mm_and_si128( b, mask );
  a = _mm_andnot_si128( mask, a );
  return _mm_or_si128( a, b );
}

#ifdef __GNUC__
# define INFINITYf __builtin_inff()
# define INFINITYd __builtin_inf()
#elif defined(_WIN32)
# define INFINITYf   HUGE_VALF
# define INFINITYd   HUGE_VALD
#endif

#ifdef _MSC_VER
# pragma warning( push )
# pragma warning( disable : 4056 )
# pragma warning( disable : 4756 )
#endif
/**
* Check if value is infinity
*
* @param d
*
* @return
*/
STATIC_INLINE_BEGIN __m128 _mm_is_infinity(__m128 d)
{
  return _mm_cmpeq_ps(_mm_fabs_ps(d), _mm_set1_ps(INFINITYf));
}

STATIC_INLINE_BEGIN __m128 _mm_is_minfinity(__m128 d)
{
  return _mm_cmpeq_ps(_mm_fabs_ps(d), _mm_set1_ps(-INFINITYf));
}

STATIC_INLINE_BEGIN __m256d _mm256_is_minfinity(__m256d d)
{
  return _mm256_cmp_pd(_mm256_fabs_pd(d), _mm256_set1_pd(-INFINITYd), 8);
}

#ifdef _MSC_VER
# pragma warning( pop )
#endif

STATIC_INLINE_BEGIN int _mm_any_eq( __m128 a, __m128 b )
{
  register __m128 mask = _mm_cmpeq_ps( a, b );
  //copy top bit of each result to maskbits
  return _mm_movemask_ps( mask ) != 0;
}

/**
 * Copy sign of sgn to val
 *
 * @param val
 * @param sgn
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128 _mm_copysign_ps(__m128 val, __m128 sgn)
{
  return _mm_or_ps(_mm_and_ps(val,mantmask),_mm_and_ps(signmask,sgn));
}

/**
 * Multiply sign of sgn to val
 *
 * @param val
 * @param sgn
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128 _mm_mulsign_ps(__m128 val, __m128 sgn)
{
  return _mm_xor_ps(val,_mm_and_ps(signmask,sgn));
}

/**
 * Negate packed singles
 *
 * @param a
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128 _mm_not_ps(__m128 a)
{
  return _mm_castsi128_ps(_mm_xor_si128(_mm_castps_si128(a),negmask));
}

STATIC_INLINE_BEGIN __m256d _mm256_not_pd(__m256d a)
{
#ifdef HAVE_ZMMINTRIN_H
  return _mm256_castsi256_pd(_mm256_xor_si256(_mm256_castpd_si256(a), negmask_256));
#else
  // _mm256_xor_si256 not inlined if AVX2 isn't present
  return _mm256_castps_pd(_mm256_xor_ps(_mm256_castpd_ps(a), _mm256_castsi256_ps(negmask_256)));
#endif
}


/**
 * Square packed singles
 *
 * @param a
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128 _mm_square_ps(__m128 a)
{
  return _mm_mul_ps(a,a);
}

STATIC_INLINE_BEGIN __m128d _mm_square_pd(__m128d a)
{
  return _mm_mul_pd(a,a);
}

STATIC_INLINE_BEGIN __m256d _mm256_square_pd(__m256d a)
{
  return _mm256_mul_pd(a,a);
}

STATIC_INLINE_BEGIN __m128d _mm_rcp_pd(__m128d a)
{
  return _mm_div_pd(_mm_set1_pd(1.0),a);
}


/**
 * Square first packed single
 *
 * @param a
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128 _mm_square_ss(__m128 a)
{
  return _mm_mul_ss(a,a);
}

/**
 * Maximum of 16 signed bytes
 *
 * @param buffer
 *
 * @return
 */
STATIC_INLINE_BEGIN int8_t _mm_max1_epi8(__m128i buffer)
{
#ifdef HAVE_SMMINTRIN_H
  __m128i tmp1 = _mm_sub_epi8(_mm_set1_epi8(127), buffer);
  __m128i tmp2 = _mm_min_epu8(tmp1, _mm_srli_epi16(tmp1, 8));

  // Lowest 16-bit are the minimum, next 16-bit are the location
  __m128i tmp3 = _mm_minpos_epu16(tmp2);
  return (int8_t)(127 - _mm_cvtsi128_si32(tmp3));
#else
# error("_mm_minpos_epu16 not available")
#endif
}

STATIC_INLINE_BEGIN __m128  _mm_madd_ps(__m128 a, __m128 b, __m128 c)
{
  /* TODO: Ensure instead that -mfma is passed on to Cygwin */
#if defined(HAVE_FMAINTRIN_H) && !defined(__CYGWIN__)
  return _mm_fmadd_ps(a,b,c);
#else
  return _mm_add_ps(_mm_mul_ps(a, b), c);
#endif
}

#if defined(HAVE_FMAINTRIN_H) || defined(_MSC_VER) || defined(__GNUC__)
/* TODO: Examine for presence of header not instructions */
#else
/**
 * For backwards-compability with the new intrinsics / instructions from Intel
 *
 * @param a
 * @param b
 * @param c
 *
 * @return
 */
STATIC_INLINE_BEGIN __m128  _mm_fmadd_ps(__m128 a, __m128 b, __m128 c)
{
  return _mm_add_ps(_mm_mul_ps(a, b), c);
}
#endif

STATIC_INLINE_BEGIN __m128 _mm_rsqrt_nr_ps( const __m128& a )
{
  __m128 r = _mm_rsqrt_ps(a);
  //  1.5 * r + (-0.5 * a) * r * r * r
  return _mm_add_ps(_mm_mul_ps(_mm_set1_ps(1.5f),r), _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(a, _m_m_half_ps), r), _mm_mul_ps(r, r)));
}

STATIC_INLINE_BEGIN __m128 rsqrt_float4_single(__m128 x)
{
  __m128 three = _mm_set1_ps(3.0f), half = _mm_set1_ps(0.5f);
  __m128 res = _mm_rsqrt_ps(x);
  __m128 muls = _mm_mul_ps(_mm_mul_ps(x, res), res);
  return res = _mm_mul_ps(_mm_mul_ps(half, res), _mm_sub_ps(three, muls));
}



/**
 * Maximum of 4 signed integers
 *
 * @param x
 *
 * @return
 */
STATIC_INLINE_BEGIN int32_t _mm_max1_epi32(__m128i x)
{
  // Move top two floats to lower part of xmm1
  __m128i max1 = _mm_shuffle_epi32(x, _MM_SHUFFLE(0,0,3,2));
  // Get maximum of the two sets of floats
  __m128i max2 = _mm_max_epi32(x,max1);
  // Move second float to lower part of xmm1
  __m128i max3 = _mm_shuffle_epi32(max2, _MM_SHUFFLE(0,0,0,1));
  // Get maximum of the two remaining floats
  __m128i max4 = _mm_max_epi32(max2,max3);

  /*
    movhlps xmm1,xmm0         ; Move top two floats to lower part of xmm1
    maxps   xmm0,xmm1         ; Get maximum of the two sets of floats
    pshufd  xmm1,xmm0,$55     ; Move second float to lower part of xmm1
    maxps   xmm0,xmm1         ; Get minimum of the two remaining floats
  */

  return _mm_cvtsi128_si32(max4);
}

/**
 * Minimum of 4 signed integers
 *
 * @param x
 *
 * @return
 */
STATIC_INLINE_BEGIN int32_t _mm_min1_epi32(__m128i x)
{
  // Move top two floats to lower part of xmm1
  __m128i min1 = _mm_shuffle_epi32(x, _MM_SHUFFLE(0,0,3,2));
  // Get minimum of the two sets of floats
  __m128i min2 = _mm_min_epi32(x,min1);
  // Move second float to lower part of xmm1
  __m128i min3 = _mm_shuffle_epi32(min2, _MM_SHUFFLE(0,0,0,1));
  // Get minimum of the two remaining floats
  __m128i min4 = _mm_min_epi32(min2,min3);

  /*
    movhlps xmm1,xmm0
    minps   xmm0,xmm1
    pshufd  xmm1,xmm0,$55
    minps   xmm0,xmm1
  */

  return _mm_cvtsi128_si32(min4);
}

STATIC_INLINE_BEGIN float _mm_max1_ps(__m128 x)
{
  // Move top two floats to lower part of xmm1
  __m128 max1 = _mm_shuffle_ps(x, x, _MM_SHUFFLE(0,0,3,2));
  // Get maximum of the two sets of floats
  __m128 max2 = _mm_max_ps(x,max1);
  // Move second float to lower part of xmm1
  __m128 max3 = _mm_shuffle_ps(max2, max2, _MM_SHUFFLE(0,0,0,1));
  // Get maximum of the two remaining floats
  __m128 max4 = _mm_max_ps(max2,max3);

  /*
    movhlps xmm1,xmm0         ; Move top two floats to lower part of xmm1
    maxps   xmm0,xmm1         ; Get maximum of the two sets of floats
    pshufd  xmm1,xmm0,$55     ; Move second float to lower part of xmm1
    maxps   xmm0,xmm1         ; Get minimum of the two remaining floats
  */
  float max;
  _mm_store_ss(&max,max4);
  return max;
}

STATIC_INLINE_BEGIN float _mm_min1_ps(__m128 x)
{
  // Move top two floats to lower part of xmm1
  __m128 min1 = _mm_shuffle_ps(x, x, _MM_SHUFFLE(0,0,3,2));
  // Get minimum of the two sets of floats
  __m128 min2 = _mm_min_ps(x,min1);
  // Move second float to lower part of xmm1
  __m128 min3 = _mm_shuffle_ps(min2, min2, _MM_SHUFFLE(0,0,0,1));
  // Get minimum of the two remaining floats
  __m128 min4 = _mm_min_ps(min2,min3);

  float min;
  _mm_store_ss(&min,min4);
  return min;
}

STATIC_INLINE_BEGIN double _mm256_min1_pd(__m256d x)
{
  __m256d ul = _mm256_castps_pd(_mm256_permute2f128_ps(_mm256_castpd_ps(x), _mm256_castpd_ps(x), 0x21));

  __m256d min0 = _mm256_min_pd(x,ul);

  __m256d min1 = _mm256_min_pd(min0, _mm256_permute_pd(min0,   0x05));

  return _mm256_cvtsd_f64(min1);
}

STATIC_INLINE_BEGIN double _mm256_max1_pd(__m256d x)
{
  __m256d ul = _mm256_castps_pd(_mm256_permute2f128_ps(_mm256_castpd_ps(x), _mm256_castpd_ps(x), 0x21));

  __m256d max0 = _mm256_max_pd(x,ul);

  __m256d max1 = _mm256_max_pd(max0, _mm256_permute_pd(max0,   0x05));

  return _mm256_cvtsd_f64(max1);
}


/*
  Transpose can be performed using nxlog2(n) rather than nxn

  n            4(SSE)          8(AVX)    16(AVX512)    32(AVX1024)
  SIMD ops          8              24           64            160
  SIMD +r/w ops    16              40           96            224
  Scalar r/w ops   24             112          480           1984
*/

/*
  __m128 tmp0 ,tmp1, tmp2, tmp3;
  tmp0 = _mm_shuffle_ps(row0, row1, 0x88); // 0 2 4 6
  tmp1 = _mm_shuffle_ps(row0, row1, 0xdd); // 1 3 5 7
  tmp2 = _mm_shuffle_ps(row2, row3, 0x88); // 8 a c e
  tmp3 = _mm_shuffle_ps(row2, row3, 0xdd); // 9 b d f

  row0 = _mm_shuffle_ps(tmp0, tmp2, 0x88); // 0 4 8 c
  row1 = _mm_shuffle_ps(tmp1, tmp3, 0x88); // 1 5 9 d
  row2 = _mm_shuffle_ps(tmp0, tmp2, 0xdd); // 2 6 a e
  row3 = _mm_shuffle_ps(tmp1, tmp3, 0xdd); // 3 7 b f
*/
STATIC_INLINE_BEGIN void _mm_transpose_4x4_ps(const float a[4][4],
    float b[4][4])
{
  __m128 row0, row1, row2, row3;
  row0 = _mm_load_ps((float*)a[0]);
  row1 = _mm_load_ps((float*)a[1]);
  row2 = _mm_load_ps((float*)a[2]);
  row3 = _mm_load_ps((float*)a[3]);

  // Inplace transpose
  _MM_TRANSPOSE4_PS(row0, row1, row2,row3);

  _mm_store_ps(b[0],row0);
  _mm_store_ps(b[1],row1);
  _mm_store_ps(b[2],row2);
  _mm_store_ps(b[3],row3);
}

// Transpose 4x4 blocks within each lane
#define _MM_TRANSPOSE8_LANE4_PS(row0, row1, row2, row3) \
    do { \
      __m256 __t0, __t1, __t2, __t3; \
      __t0 = _mm256_unpacklo_ps(row0, row1); \
      __t1 = _mm256_unpackhi_ps(row0, row1); \
      __t2 = _mm256_unpacklo_ps(row2, row3); \
      __t3 = _mm256_unpackhi_ps(row2, row3); \
      row0 = _mm256_shuffle_ps(__t0, __t2, _MM_SHUFFLE(5, 4, 1, 0)); \
      row1 = _mm256_shuffle_ps(__t0, __t2, _MM_SHUFFLE(7, 6, 3, 2)); \
      row2 = _mm256_shuffle_ps(__t1, __t3, _MM_SHUFFLE(5, 4, 1, 0)); \
      row3 = _mm256_shuffle_ps(__t1, __t3, _MM_SHUFFLE(7, 6, 3, 2)); \
    } while (0)

#ifdef _MSC_VER
# define SPS_IGNORE_CONSTANT_CONDITION_BEGIN \
    __pragma(warning(push))                 \
    __pragma(warning(disable: 4127))
# define SPS_IGNORE_CONSTANT_CONDITION_END \
    __pragma(warning(pop))
#else
# define SPS_IGNORE_CONSTANT_CONDITION_BEGIN
# define SPS_IGNORE_CONSTANT_CONDITION_END
#endif

// http://stackoverflow.com/questions/25622745/transpose-an-8x8-float-using-avx-avx2
#define _MM_TRANSPOSE8_PS(row0, row1, row2, row3, row4, row5, row6, row7) \
SPS_IGNORE_CONSTANT_CONDITION_BEGIN \
    do { \
      __m256 __t0, __t1, __t2, __t3, __t4, __t5, __t6, __t7; \
      __m256 __tt0, __tt1, __tt2, __tt3, __tt4, __tt5, __tt6, __tt7; \
      __t0 = _mm256_unpacklo_ps(row0, row1); \
      __t1 = _mm256_unpackhi_ps(row0, row1); \
      __t2 = _mm256_unpacklo_ps(row2, row3); \
      __t3 = _mm256_unpackhi_ps(row2, row3); \
      __t4 = _mm256_unpacklo_ps(row4, row5); \
      __t5 = _mm256_unpackhi_ps(row4, row5); \
      __t6 = _mm256_unpacklo_ps(row6, row7); \
      __t7 = _mm256_unpackhi_ps(row6, row7); \
      __tt0 = _mm256_shuffle_ps(__t0, __t2, _MM_SHUFFLE(1, 0, 1, 0)); \
      __tt1 = _mm256_shuffle_ps(__t0, __t2, _MM_SHUFFLE(3, 2, 3, 2)); \
      __tt2 = _mm256_shuffle_ps(__t1, __t3, _MM_SHUFFLE(1, 0, 1, 0)); \
      __tt3 = _mm256_shuffle_ps(__t1, __t3, _MM_SHUFFLE(3, 2, 3, 2)); \
      __tt4 = _mm256_shuffle_ps(__t4, __t6, _MM_SHUFFLE(1, 0, 1, 0)); \
      __tt5 = _mm256_shuffle_ps(__t4, __t6, _MM_SHUFFLE(3, 2, 3, 2)); \
      __tt6 = _mm256_shuffle_ps(__t5, __t7, _MM_SHUFFLE(1, 0, 1, 0)); \
      __tt7 = _mm256_shuffle_ps(__t5, __t7, _MM_SHUFFLE(3, 2, 3, 2)); \
      row0 = _mm256_permute2f128_ps(__tt0, __tt4, 0x20); \
      row1 = _mm256_permute2f128_ps(__tt1, __tt5, 0x20); \
      row2 = _mm256_permute2f128_ps(__tt2, __tt6, 0x20); \
      row3 = _mm256_permute2f128_ps(__tt3, __tt7, 0x20); \
      row4 = _mm256_permute2f128_ps(__tt0, __tt4, 0x31); \
      row5 = _mm256_permute2f128_ps(__tt1, __tt5, 0x31); \
      row6 = _mm256_permute2f128_ps(__tt2, __tt6, 0x31); \
      row7 = _mm256_permute2f128_ps(__tt3, __tt7, 0x31); \
    } while (0) \
SPS_IGNORE_CONSTANT_CONDITION_END



STATIC_INLINE_BEGIN void _mm_transpose_8x8_ps(const float a[8][8],
    float b[8][8])
{

  __m256 row0, row1, row2, row3, row4, row5, row6, row7;
  row0 = _mm256_load_ps((float*)a[0]);
  row1 = _mm256_load_ps((float*)a[1]);
  row2 = _mm256_load_ps((float*)a[2]);
  row3 = _mm256_load_ps((float*)a[3]);
  row4 = _mm256_load_ps((float*)a[4]);
  row5 = _mm256_load_ps((float*)a[5]);
  row6 = _mm256_load_ps((float*)a[6]);
  row7 = _mm256_load_ps((float*)a[7]);

  // Inplace transpose
  _MM_TRANSPOSE8_PS(row0, row1, row2,row3,
                    row4, row5, row6, row7);

  _mm256_store_ps(b[0],row0);
  _mm256_store_ps(b[1],row1);
  _mm256_store_ps(b[2],row2);
  _mm256_store_ps(b[3],row3);
  _mm256_store_ps(b[4],row4);
  _mm256_store_ps(b[5],row5);
  _mm256_store_ps(b[6],row6);
  _mm256_store_ps(b[7],row7);
}


#if 0

void _mm_transpose_16x16_epi32(const int a[16][16],
                               int b[16][16])
{

  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, ra, rb, rc, rd, re, rf;

  const int* input = &(a[0][0]);
  int* output = &(b[0][0]);

  r0 = _mm512_load_epi32((void const*)&input[0]  );
  r1 = _mm512_load_epi32((void const*)&input[16] );
  r2 = _mm512_load_epi32((void const*)&input[32] );
  r3 = _mm512_load_epi32((void const*)&input[48] );
  r4 = _mm512_load_epi32((void const*)&input[64] );
  r5 = _mm512_load_epi32((void const*)&input[80] );
  r6 = _mm512_load_epi32((void const*)&input[96] );
  r7 = _mm512_load_epi32((void const*)&input[112]);
  r8 = _mm512_load_epi32((void const*)&input[128]);
  r9 = _mm512_load_epi32((void const*)&input[144]);
  ra = _mm512_load_epi32((void const*)&input[160]);
  rb = _mm512_load_epi32((void const*)&input[176]);
  rc = _mm512_load_epi32((void const*)&input[192]);
  rd = _mm512_load_epi32((void const*)&input[208]);
  re = _mm512_load_epi32((void const*)&input[224]);
  rf = _mm512_load_epi32((void const*)&input[240]);

  __m512i t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, ta, tb, tc, td, te, tf;

  t0 = _mm512_unpacklo_epi32(r0,r1); //   0  16   1  17   4  20   5  21   8  24   9  25  12  28  13  29
  t1 = _mm512_unpackhi_epi32(r0,r1); //   2  18   3  19   6  22   7  23  10  26  11  27  14  30  15  31
  t2 = _mm512_unpacklo_epi32(r2,r3); //  32  48  33  49 ...
  t3 = _mm512_unpackhi_epi32(r2,r3); //  34  50  35  51 ...
  t4 = _mm512_unpacklo_epi32(r4,r5); //  64  80  65  81 ...
  t5 = _mm512_unpackhi_epi32(r4,r5); //  66  82  67  83 ...
  t6 = _mm512_unpacklo_epi32(r6,r7); //  96 112  97 113 ...
  t7 = _mm512_unpackhi_epi32(r6,r7); //  98 114  99 115 ...
  t8 = _mm512_unpacklo_epi32(r8,r9); // 128 ...
  t9 = _mm512_unpackhi_epi32(r8,r9); // 130 ...
  ta = _mm512_unpacklo_epi32(ra,rb); // 160 ...
  tb = _mm512_unpackhi_epi32(ra,rb); // 162 ...
  tc = _mm512_unpacklo_epi32(rc,rd); // 196 ...
  td = _mm512_unpackhi_epi32(rc,rd); // 198 ...
  te = _mm512_unpacklo_epi32(re,rf); // 228 ...
  tf = _mm512_unpackhi_epi32(re,rf); // 230 ...

  r0 = _mm512_unpacklo_epi64(t0,t2); //   0  16  32  48 ...
  r1 = _mm512_unpackhi_epi64(t0,t2); //   1  17  33  49 ...
  r2 = _mm512_unpacklo_epi64(t1,t3); //   2  18  34  49 ...
  r3 = _mm512_unpackhi_epi64(t1,t3); //   3  19  35  51 ...
  r4 = _mm512_unpacklo_epi64(t4,t6); //  64  80  96 112 ...
  r5 = _mm512_unpackhi_epi64(t4,t6); //  65  81  97 114 ...
  r6 = _mm512_unpacklo_epi64(t5,t7); //  66  82  98 113 ...
  r7 = _mm512_unpackhi_epi64(t5,t7); //  67  83  99 115 ...
  r8 = _mm512_unpacklo_epi64(t8,ta); // 128 144 160 176 ...
  r9 = _mm512_unpackhi_epi64(t8,ta); // 129 145 161 178 ...
  ra = _mm512_unpacklo_epi64(t9,tb); // 130 146 162 177 ...
  rb = _mm512_unpackhi_epi64(t9,tb); // 131 147 163 179 ...
  rc = _mm512_unpacklo_epi64(tc,te); // 192 208 228 240 ...
  rd = _mm512_unpackhi_epi64(tc,te); // 193 209 229 241 ...
  re = _mm512_unpacklo_epi64(td,tf); // 194 210 230 242 ...
  rf = _mm512_unpackhi_epi64(td,tf); // 195 211 231 243 ...

  t0 = _mm512_shuffle_i32x4(r0, r4, 0x88); //   0  16  32  48   8  24  40  56  64  80  96  112 ...
  t1 = _mm512_shuffle_i32x4(r1, r5, 0x88); //   1  17  33  49 ...
  t2 = _mm512_shuffle_i32x4(r2, r6, 0x88); //   2  18  34  50 ...
  t3 = _mm512_shuffle_i32x4(r3, r7, 0x88); //   3  19  35  51 ...
  t4 = _mm512_shuffle_i32x4(r0, r4, 0xdd); //   4  20  36  52 ...
  t5 = _mm512_shuffle_i32x4(r1, r5, 0xdd); //   5  21  37  53 ...
  t6 = _mm512_shuffle_i32x4(r2, r6, 0xdd); //   6  22  38  54 ...
  t7 = _mm512_shuffle_i32x4(r3, r7, 0xdd); //   7  23  39  55 ...
  t8 = _mm512_shuffle_i32x4(r8, rc, 0x88); // 128 144 160 176 ...
  t9 = _mm512_shuffle_i32x4(r9, rd, 0x88); // 129 145 161 177 ...
  ta = _mm512_shuffle_i32x4(ra, re, 0x88); // 130 146 162 178 ...
  tb = _mm512_shuffle_i32x4(rb, rf, 0x88); // 131 147 163 179 ...
  tc = _mm512_shuffle_i32x4(r8, rc, 0xdd); // 132 148 164 180 ...
  td = _mm512_shuffle_i32x4(r9, rd, 0xdd); // 133 149 165 181 ...
  te = _mm512_shuffle_i32x4(ra, re, 0xdd); // 134 150 166 182 ...
  tf = _mm512_shuffle_i32x4(rb, rf, 0xdd); // 135 151 167 183 ...

  r0 = _mm512_shuffle_i32x4(t0, t8, 0x88); //   0  16  32  48  64  80  96 112 ... 240
  r1 = _mm512_shuffle_i32x4(t1, t9, 0x88); //   1  17  33  49  66  81  97 113 ... 241
  r2 = _mm512_shuffle_i32x4(t2, ta, 0x88); //   2  18  34  50  67  82  98 114 ... 242
  r3 = _mm512_shuffle_i32x4(t3, tb, 0x88); //   3  19  35  51  68  83  99 115 ... 243
  r4 = _mm512_shuffle_i32x4(t4, tc, 0x88); //   4 ...
  r5 = _mm512_shuffle_i32x4(t5, td, 0x88); //   5 ...
  r6 = _mm512_shuffle_i32x4(t6, te, 0x88); //   6 ...
  r7 = _mm512_shuffle_i32x4(t7, tf, 0x88); //   7 ...
  r8 = _mm512_shuffle_i32x4(t0, t8, 0xdd); //   8 ...
  r9 = _mm512_shuffle_i32x4(t1, t9, 0xdd); //   9 ...
  ra = _mm512_shuffle_i32x4(t2, ta, 0xdd); //  10 ...
  rb = _mm512_shuffle_i32x4(t3, tb, 0xdd); //  11 ...
  rc = _mm512_shuffle_i32x4(t4, tc, 0xdd); //  12 ...
  rd = _mm512_shuffle_i32x4(t5, td, 0xdd); //  13 ...
  re = _mm512_shuffle_i32x4(t6, te, 0xdd); //  14 ...
  rf = _mm512_shuffle_i32x4(t7, tf, 0xdd); //  15  31  47  63  79  96 111 127 ... 255

  _mm512_store_epi32(&output[0]  , t0);
  _mm512_store_epi32(&output[16] , t1);
  _mm512_store_epi32(&output[32] , t2);
  _mm512_store_epi32(&output[48] , t3);
  _mm512_store_epi32(&output[64] , t4);
  _mm512_store_epi32(&output[80] , t5);
  _mm512_store_epi32(&output[96] , t6);
  _mm512_store_epi32(&output[112], t7);
  _mm512_store_epi32(&output[128], t8);
  _mm512_store_epi32(&output[144], t9);
  _mm512_store_epi32(&output[160], ta);
  _mm512_store_epi32(&output[176], tb);
  _mm512_store_epi32(&output[192], tc);
  _mm512_store_epi32(&output[208], td);
  _mm512_store_epi32(&output[224], te);
  _mm512_store_epi32(&output[240], tf);
}
#endif

#ifdef HAVE_IMMINTRIN_H
# if defined(HAVE_FMAINTRIN_H) || defined(_MSC_VER) || defined(__GNUC__)
/* TODO: Examine for presence of header not instructions */
# else
STATIC_INLINE_BEGIN __m256 _mm256_fmadd_ps(__m256 a, __m256 b, __m256 c)
{
  return _mm256_add_ps(_mm256_mul_ps(a,b),c);
}
# endif
#endif

#ifdef __cplusplus
}
#endif

// Ugly stuff (remove)
#ifdef __cplusplus

// #include <type_traits>
// typedef std::aligned_union<sizeof(U),int,char,double>::type U_pod;

template<unsigned i>
float vectorGetByIndex( __m128 V)
{
  v4f converter;
  converter.v = V;
  return converter.f32[i];
}

#include <iostream>
inline std::ostream& operator<<(std::ostream& out, const __m128& vf)
{
#if (defined(_MSC_VER) && defined(_WIN32) && !defined(__ICC))
  out << "__m128[3]: " << (float) vf.m128_f32[3] << " __m128[2]: " << (float) vf.m128_f32[2] \
      << " __m128[1]: " << (float) vf.m128_f32[1] << " __m128[0]: " << (float) vf.m128_f32[0] << std::endl;
#else
  ALIGN16_BEGIN float output[4] ALIGN16_END;
  _mm_store_ps((float*)output,vf);
  out << "__m128[3]: " << (float) output[3] << " __m128[2]: " << (float) output[2] \
      << " __m128[1]: " << (float) output[1] << " __m128[0]: " << (float) output[0] << std::endl;
#endif
  return out;
}

inline std::ostream& operator<<(std::ostream& out, __m128i vi)
{
  ALIGN16_BEGIN int output[4] ALIGN16_END;
  _mm_store_si128((__m128i*)output,vi);
  out << "__m128[3]: " << (int) output[3] << " __m128[2]: " << (int) output[2] \
      << " __m128[1]: " << (int) output[1] << " __m128[0]: " << (int) output[0] << " " << std::endl;
  return out;
}

inline std::ostream& operator<<(std::ostream& out, __m64 vi)
{
#if (defined(_MSC_VER) && defined(_WIN32) && !defined(__ICC))
  out << "__m64[3]: " << vi.m64_i16[3] << " __m64[2]: " << vi.m64_i16[2] \
      << " __m64[1]: " << vi.m64_i16[1] << " __m64[0]: " << vi.m64_i16[0] << std::endl;
#else
  // No need to respect alignment
  ALIGN16_BEGIN short int output[4] ALIGN16_END;
  _mm_stream_pi((__m64*)&output[0],vi);
  out << "__m64[3]: " << output[3] << " __m64[2]: " << output[2]        \
      << " __m64[1]: " << output[1] << " __m64[0]: " << output[0] << std::endl;
#endif
  return out;
}
#endif

#if 1
# define SPS_SPLATS(x) \
  _mm_set1_ps(x)
#else
# define SPS_SPLATS(x) \
  _mm_broadcast_ss((float*)& x)
#endif

#if 0
__m128 rsqrt_float4_single(__m128 x)
{
  __m128 three = _mm_set1_ps(3.0f), half = _mm_set1_ps(0.5f);
  __m128 res = _mm_rsqrt_ps(x);
  __m128 muls = _mm_mul_ps(_mm_mul_ps(x, res), res);
  return res = _mm_mul_ps(_mm_mul_ps(half, res), _mm_sub_ps(three, muls));
}
#endif

/*
_mm_broadcast_ss has weaknesses imposed by the architecture which are
largely hidden by the mm SSE API. The most important difference is as
follows:

_mm_broadcast_ss is limited to loading values from memory only.  What
this means is if you use _mm_broadcast_ss explicitly in a situation
where the source is not in memory then the result will likely be less
efficient than that of using _mm_set1_ps. This sort of situation
typically happens when loading immediate values (constants), or when
using the result of a recent calculation. In those situations the
result will be mapped to a register by the compiler. To use the value
for broadcast, the compiler must dump the value back to
memory. Alternatively, a pshufd could be used to splat directly from
register instead.


_mm_set1_ps is implementation-defined rather than being mapped to a
specific underlying cpu operation (instruction). That means it might
use one of several SSE instructions to perform the splat. A smart
compiler with AVX support enabled should definitely use vbroadcastss
internally when appropriate, but it depends on the AVX implementation
state of the compilers optimizer.

If you're very confident you're loading from memory -- such as
iterating over an array of data -- then direct use of broadcast is
fine. But if there's any doubt at all, I would recommend stick with
_mm_set1_ps.

And in the specific case of a static const float, you absolutely want
to avoid using _mm_broadcast_ss().

 */

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
