/**
 * @file   extintrin.h
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Tue May  5 20:13:58 2015
 * 
 * @brief  Extensions to Intel's intrinsics
 * 
 * 
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

#ifdef HAVE_ZMMINTRIN_H
# include <zmmintrin.h> // AVX2
#endif

#ifdef HAVE_STDINT_H
# include <stdint.h>
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

  /**@}*/

  /**
   * \defgroup Constants used by functions
   * @{
   */
  const ALIGN16_BEGIN int _clear_signmask[4] ALIGN16_END =
    {0x7FFFFFFFL,0x7FFFFFFFL,0x7FFFFFFFL,0x7FFFFFFFL};                       ///< clear sign mask for ps

  const ALIGN16_BEGIN long long int _clear_signmaskd[2] ALIGN16_END =
    {0x7FFFFFFFFFFFFFFF,0x7FFFFFFFFFFFFFFF};                                 ///< clear sign mask for pd

  
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
  
  const __m128i negmask = _mm_set1_epi32(-1);                                ///< negate mask 0xFFFFFFFF for epi32


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
    v4f(const __m128& _v) {v = _v;}
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
    v4i(const __m128i& _v) {v = _v;}
  } v4i;
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
    v4f(const __m128& _v) {v = _v;}
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
    v4i(const __m128i& _v) {v = _v;}
  } v4i __attribute__ ((aligned (16)));
#endif

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
    __asm__ volatile (                                // rax = a, rdx = b, ecx = s
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

  STATIC_INLINE_BEGIN __m128 _mm_rcp_nr_ps(const __m128& x) {
    __m128 res = _mm_rcp_ps(x);
    __m128 muls = _mm_mul_ps(x, _mm_mul_ps(res, res));
    return res =  _mm_sub_ps(_mm_add_ps(res, res), muls);
  }
  
  
  /** 
   * Absolute values of packed singles.
   * 
   * @param x 
   * 
   * @return 
   */
  STATIC_INLINE_BEGIN  __m128 _mm_fabs_ps(__m128 x) STATIC_INLINE_END;

  STATIC_INLINE_BEGIN  __m128d _mm_fabs_pd(__m128d x) STATIC_INLINE_END;

  STATIC_INLINE_BEGIN __m128 _mm_fmod_ps(const __m128& a, const __m128& aDiv) STATIC_INLINE_END;

  STATIC_INLINE_BEGIN __m128 _mm_fmod_ps(const __m128& a, const __m128& aDiv) {
    __m128 c = _mm_div_ps(a,aDiv);
    __m128i i = _mm_cvttps_epi32(c);
    __m128 cTrunc = _mm_cvtepi32_ps(i);
    __m128 base = _mm_mul_ps(cTrunc, aDiv);
    __m128 r = _mm_sub_ps(a, base);
    return r;
  }

#ifndef HAVE_SMMINTRIN_H
  //SSE2: multiply 8 16-bit integers and return full 32-bit result in high and low result
  STATIC_INLINE_BEGIN void _mm_mul_epi16_full(__m128i &hiResult, __m128i &loResult, const __m128i a, const __m128i b) {
    __m128i hi  = _mm_mulhi_epi16( a, b );    // (a7*b7[16:31],a6*b6[16:31],a5*b5[16:31],a4*b4[16:31],a3*b3[16:31],a2*b2[16:31],a1*b1[16:31],a0*b0[16:31])
    __m128i low = _mm_mullo_epi16( a, b );    // (a7*b7[0:15] ,a6*b6[0:15] ,a5*b5[0:15] ,a4*b4[0:15],a3*b3[0:15] ,a2*b2[0:15] ,a1*b1[0:15] ,a0*b0[0:15])
    loResult = _mm_unpacklo_epi16( hi, low ); // (a3*b3[0:15],a3*b3[16:31],a2*b2[0:15],a2*b2[16:31],a1*b1[0:15],a1*b1[16:31],a0*b0[0:15],a0*b0[16:31])
    hiResult = _mm_unpackhi_epi16( hi, low ); // (a7*b7[0:15],a7*b7[16:31],a6*b6[0:15],a6*b6[16:31],a5*b5[0:15],a5*b5[16:31],a4*b4[0:15],a4*b4[16:31])
  }
#endif

  STATIC_INLINE_BEGIN void _mm_mul_epi32_full(__m128i &loResult, __m128i &hiResult, __m128i a, __m128i b) {
    __m128i _a        = _mm_castps_si128(_mm_permute_ps(_mm_castsi128_ps(a),0xB1));
    __m128i _b        = _mm_castps_si128(_mm_permute_ps(_mm_castsi128_ps(b),0xB1));
    __m128i c1        = _mm_mul_epi32(a,b);   // 0 and 2
    __m128i c2        = _mm_mul_epi32(_a,_b); // 1 and 3
    loResult          = _mm_unpacklo_epi64(c1,c2);
    hiResult          = _mm_unpackhi_epi64(c1,c2);
  }


  STATIC_INLINE_BEGIN __m128 _mm_mulcmplx_ps(const __m128 &a, const __m128 &b) {

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

  // It is actually vertical

  /** 
   * Horizontal AND values, 32-bit at a time. Result is placed in the
   * least significant 32-bit integer.
   * 
   * @param x 
   * 
   * @return 
   */
  STATIC_INLINE_BEGIN __m128i _mm_hand_epi32(__m128i x) {
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
  STATIC_INLINE_BEGIN __m128i _mm_hor_epi32(__m128i x) {
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
  STATIC_INLINE_BEGIN __m128 _mm_srli_ps(__m128 x, const int imm) {
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
  STATIC_INLINE_BEGIN __m128 _mm_slli_ps(__m128 x, const int imm) {
    return _mm_castsi128_ps(_mm_slli_si128(_mm_castps_si128(x),imm));
  }
#else
  // 8-bit const immediates cannot be considered immediate when debugging
# define _mm_slli_ps(x,imm)                                          \
           _mm_castsi128_ps(_mm_slli_si128(_mm_castps_si128(x),imm))
#endif


#if HAVE_SMMINTRIN_H
  STATIC_INLINE_BEGIN int _mm_dp1_epi32(__m128i x, __m128i y) {
  
    __m128i s = _mm_mullo_epi32(x,y);
    s = _mm_add_epi32(s, _mm_srli_si128(s, 8));
    s = _mm_add_epi32(s, _mm_srli_si128(s, 4));
    int sum = _mm_cvtsi128_si32(s);
    return sum;
  }

  STATIC_INLINE_BEGIN __m128i _mm_dp_epi32(__m128i x, __m128i y) {
  
    __m128i s = _mm_mullo_epi32(x,y);
    s = _mm_add_epi32(s, _mm_srli_si128(s, 8));
    s = _mm_add_epi32(s, _mm_srli_si128(s, 4));
    return s;
  }
#endif

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
  STATIC_INLINE_BEGIN __m128 _mm_dp_ps(__m128 __X, __m128 __Y, const int __M) {

# if HAVE_PMMINTRIN_H

    const __m128 smask      = _mm_set_epi32(0x80,0x40,0x20,0x10);
    const __m128 omask      = _mm_set_epi32(0x08,0x04,0x02,0x01);
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

  /** 
   * Absolute value of packed singles 
   * 
   * @param x 
   * 
   * @return 
   */
  STATIC_INLINE_BEGIN  __m128 _mm_fabs_ps(__m128 x) {
    return _mm_and_ps(x,_mm_load_ps((float *) _clear_signmask));
  }

  STATIC_INLINE_BEGIN  __m128d _mm_fabs_pd(__m128d x) {
    return _mm_and_pd(x,_mm_load_pd((double *) _clear_signmask));
  }
  
  /** 
   * Negate packed singles.
   * 
   * @param x 
   * 
   * @return 
   */
  STATIC_INLINE_BEGIN __m128 _mm_neg_ps(__m128 x) {
    return _mm_xor_ps(x,_mm_load_ps((float *) _neg_signmask));
  }

  STATIC_INLINE_BEGIN __m128d _mm_neg_pd(__m128d x) {
    return _mm_xor_pd(x,_mm_load_pd((double *) _neg_signmaskd));
  }
  
  STATIC_INLINE_BEGIN __m128i _mm_neg_epi32(__m128i x) {
    return _mm_sub_epi32(_mm_setzero_si128(), x);
    // or (x ^ 0xFFFFFFFF) + 1
  }
  
  /** 
   * Select floating point result based on mask.
   * 
   * @param a 
   * @param b 
   * @param mask 
   * 
   * @return If mask is false, a is returned, otherwise b is returned.
   */
  STATIC_INLINE_BEGIN __m128 _mm_sel_ps(__m128 a, __m128 b, __m128 mask ) {
#ifdef HAVE_SMMINTRIN_H
    return _mm_blendv_ps(a,b,mask);
#else
    b = _mm_and_ps( b, mask );
    a = _mm_andnot_ps( mask, a );
    return _mm_or_ps( a, b );
#endif
  }

  STATIC_INLINE_BEGIN __m128d _mm_sel_pd(__m128d a, __m128d b, __m128d mask ) {
#ifdef HAVE_SMMINTRIN_H
    return _mm_blendv_pd(a,b,mask);
#else
    b = _mm_and_pd( b, mask );
    a = _mm_andnot_pd( mask, a );
    return _mm_or_pd( a, b );
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
  STATIC_INLINE_BEGIN __m128 _mm_blendv_ps(__m128 a, __m128 b, __m128 mask ) {
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
  STATIC_INLINE_BEGIN __m128i _mm_sel_epi32(__m128i a, __m128i b, __m128i mask ) {
    // Using _mm_blendv_ps is slower due to register usage
    b = _mm_and_si128( b, mask );
    a = _mm_andnot_si128( mask, a );
    return _mm_or_si128( a, b );
  }

#ifdef __GNUC__
# define INFINITYf __builtin_inff()
#elif defined(_WIN32)
# define INFINITYf   HUGE_VALF
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
   STATIC_INLINE_BEGIN __m128 _mm_is_infinity(__m128 d) {
     return _mm_cmpeq_ps(_mm_fabs_ps(d), _mm_set1_ps(INFINITYf));
   }

   STATIC_INLINE_BEGIN __m128 _mm_is_minfinity(__m128 d) {
     return _mm_cmpeq_ps(_mm_fabs_ps(d), _mm_set1_ps(-INFINITYf));
   }
#ifdef _MSC_VER
# pragma warning( pop )
#endif

   STATIC_INLINE_BEGIN int _mm_any_eq( __m128 a, __m128 b ) {
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
  STATIC_INLINE_BEGIN __m128 _mm_copysign_ps(__m128 val, __m128 sgn) {
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
  STATIC_INLINE_BEGIN __m128 _mm_mulsign_ps(__m128 val, __m128 sgn) {
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
  STATIC_INLINE_BEGIN __m128 _mm_square_ss(__m128 a) {
    return _mm_mul_ss(a,a);
  }

  /** 
   * Maximum of 16 signed bytes
   * 
   * @param buffer 
   * 
   * @return 
   */  
  STATIC_INLINE_BEGIN int8_t _mm_max1_epi8(__m128i buffer) {
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

  STATIC_INLINE_BEGIN __m128  _mm_madd_ps(__m128 a, __m128 b, __m128 c) {
#ifdef HAVE_FMAINTRIN_H
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
  STATIC_INLINE_BEGIN __m128  _mm_fmadd_ps(__m128 a, __m128 b, __m128 c) {
    return _mm_add_ps(_mm_mul_ps(a, b), c);
  }
#endif

  STATIC_INLINE_BEGIN __m128 _mm_rsqrt_nr_ps( const __m128& a ) {
    __m128 r = _mm_rsqrt_ps(a);
    //  1.5 * r + (-0.5 * a) * r * r * r
    return _mm_add_ps(_mm_mul_ps(_mm_set1_ps(1.5f),r), _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(a, _m_m_half_ps), r), _mm_mul_ps(r, r)));
  }

  STATIC_INLINE_BEGIN __m128 rsqrt_float4_single(__m128 x) {
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
  STATIC_INLINE_BEGIN int32_t _mm_max1_epi32(__m128i x) {
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
  STATIC_INLINE_BEGIN int32_t _mm_min1_epi32(__m128i x) {
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

  STATIC_INLINE_BEGIN float _mm_max1_ps(__m128 x) {
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

  STATIC_INLINE_BEGIN float _mm_min1_ps(__m128 x) {
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

#ifdef HAVE_IMMINTRIN_H
# if defined(HAVE_FMAINTRIN_H) || defined(_MSC_VER) || defined(__GNUC__)
  /* TODO: Examine for presence of header not instructions */  
# else
  STATIC_INLINE_BEGIN __m256 _mm256_fmadd_ps(__m256 a, __m256 b, __m256 c) {
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
  float vectorGetByIndex( __m128 V) {
    v4f converter;
    converter.v = V;
    return converter.f32[i];
  }

#include <iostream>
inline std::ostream& operator<<(std::ostream& out, const __m128& vf) {
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

inline std::ostream& operator<<(std::ostream& out, __m128i vi) {
  ALIGN16_BEGIN int output[4] ALIGN16_END;
  _mm_store_si128((__m128i*)output,vi);
  out << "__m128[3]: " << (int) output[3] << " __m128[2]: " << (int) output[2] \
      << " __m128[1]: " << (int) output[1] << " __m128[0]: " << (int) output[0] << " " << std::endl;  
  return out;
}

inline std::ostream& operator<<(std::ostream& out, __m64 vi) {
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
__m128 rsqrt_float4_single(__m128 x) {
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
