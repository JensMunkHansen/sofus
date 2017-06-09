/**
 * @file   math.h
 * @author Jens Munk Hansen <jens.munk.hansen\@gmail.com>
 * @date   Sat Jun 20 20:01:01 2015
 *
 * @brief  Interface to math.h. Takes care of missing C99 functions
 *
 *
 */
#pragma once

#include <sps/cenv.h>

// TODO: This is not C
#include <cstddef>

//#ifdef __GNUC__
# include<limits.h>
//#endif

#ifndef SQUARE
# define SQUARE(z) ((z) * (z))
#endif

#ifndef CUBIC
# define CUBIC(z) ((z) * (z) * (z))
#endif

#if defined(_MSC_VER)
# define _USE_MATH_DEFINES 1
# define NOMINMAX 1
#endif

#include <math.h>

/**
 * Next power of two. Works until 2^30
 *
 * @param x
 *
 * @return
 */
STATIC_INLINE_BEGIN int next_power_two (int x) STATIC_INLINE_END;

#if defined(__cplusplus)

# if defined(CXX11) || defined(C11)
#  ifdef __GNUG__
/**
 * Is power of two
 *
 * @param x
 *
 * @return
 */
constexpr bool is_power_of_two(size_t x)
{
  return x && ((x & (x-1)) == 0);
}
#  endif
# endif

/**
 * Next power of two inlined and using __builtin_clz
 *
 * Works for signed and unsigned until 2 ^ (sizeof(T)*8-1)
 *
 * @param value
 *
 * @return
 */
template<typename T>
STATIC_INLINE_BEGIN T next_power_two(T value) STATIC_INLINE_END;

# ifdef __GNUG__

template<typename T> T next_power_two(T value)
{
  --value;
  for(size_t i = 1; i < sizeof(T) * CHAR_BIT; i*=2)
    value |= value >> i;
  return value+1;
}

template <> int
next_power_two(int value)
{
  return 1 << ((sizeof(int) * CHAR_BIT) - __builtin_clz(value-1));
}

template <> long int
next_power_two(long int value)
{
  return 1 << ((sizeof(int) * CHAR_BIT) - __builtin_clzl(value-1));
}

# elif defined(_MSC_VER)

template<typename T> T next_power_two(T value)
{
  --value;
  for(size_t i = 1; i < sizeof(T) * CHAR_BIT ; i*=2)
    value |= value >> i;
  return value+1;
}

# endif
#endif

// C++11, C++0x enforces strict ANSI and M_PI etc. are undefined
#ifdef __STRICT_ANSI__
# ifndef M_PI
constexpr auto M_PI     = 3.14159265358979323846;
# endif
# ifndef M_PI_2
constexpr auto M_PI_2   = 1.57079632679489661923;
# endif
constexpr auto M_2PI    = 6.28318530717958623200;
# ifndef M_3PI_2
constexpr auto M_3PI_2  = 4.71238898038468967400;
# endif
# ifndef M_1_PI
constexpr auto M_1_PI   = 0.31830988618379067154;
# endif
# ifndef M_2_PI
constexpr auto M_2_PI   = 0.63661977236758134308; // 2 / PI
# endif
# ifndef M_1_2PI
constexpr auto M_1_2PI   = 0.15915494309189535; // 1 / (2 PI)
# endif
#else
# ifndef M_2PI
#  define M_2PI   6.28318530717958623200
# endif
# ifndef M_3PI_2
#  define M_3PI_2 4.71238898038468967400
# endif
# ifndef M_1_PI
#  define M_1_PI  0.31830988618379067154
# endif
# ifndef M_2_PI
#  define M_2_PI  0.63661977236758134308
# endif
# ifndef M_1_2PI
#  define M_1_2PI 0.15915494309189535 // 1 / (2 PI)
# endif
#endif

STATIC_INLINE_BEGIN int next_power_two (int x)
{
  if (x < 0)
    return 0;
  --x;
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  return x+1;
}

// Missing features from C99 (TODO: Support C)
#ifndef C99
# include <stdint.h>
# include <assert.h>
# if defined(__cplusplus)
#  ifdef __GNUC__
static inline int log2i(int x)
{
  assert(x > 0);
  return sizeof(int) * CHAR_BIT - __builtin_clz(x) - 1;
}
#  elif defined(_MSC_VER)
#   if  (_MSC_VER < 1800)
static inline double log2( double n )
{
  return log( n ) / log( 2.0 );
}
#   endif
#  else
// Compilers other than Microsoft all support inline assembly (TODO: Implement log2i for Microsoft)
static inline int log2i(uint32_t x)
{
  uint32_t y;
  asm ( "\tbsr %1, %0\n"
        : "=r"(y)
        : "r" (x)
      );
  return y;
}
#  endif
# endif
#endif

/* Local variables: */
/* indent-tab-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* indent-tabs-mode: nil */
/* End: */
