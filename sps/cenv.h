/**
 * @file   cenv.h
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Sat Nov 29 17:48:45 2014
 * 
 * @brief  Environment macros introduced for portability
 *         This file must be kept C compliant
 * 
 * 
 */
#pragma once

// TODO: Add other OS'es
#if defined(_WIN32) || defined(WIN32) || defined(__CYGWIN__) || defined(__MINGW32__) || defined(__BORLANDC__)
# define SPS_OS_WIN
#endif

#ifndef SPS_UNREFERENCED_PARAMETER
# define SPS_UNREFERENCED_PARAMETER(x) ((void)(x))
#endif

#ifdef __GNUC__
# include <features.h>
#endif

#if defined(__STDC__)
# define C89
# if defined(__STDC_VERSION__)
#  define C90
#  if (__STDC_VERSION__ >= 199409L)
#   define C94
#  endif
#  if (__STDC_VERSION__ >= 199901L)
#   define C99
#  endif
#  if (__STDC_VERSION__ >= 201112L)
#   define C11
#  endif
# elif defined(__cplusplus)
#  ifdef __USE_ISOC99
#   define C99
#  endif
# endif
#endif

// Alignment
#if (defined(_MSC_VER) && defined(_WIN32))
# define ALIGN16_BEGIN __declspec(align(16))
# define ALIGN16_END
#elif defined(__GNUC__)
# define ALIGN16_BEGIN
# define ALIGN16_END __attribute__((aligned(16)))
#endif

// Static inlines
#if (defined(_MSC_VER) && defined(_WIN32))
// Note when used inside a namespace, the static is superfluous
# define STATIC_INLINE_BEGIN static inline //__forceinline
# define STATIC_INLINE_END
#elif (defined(__GNUC__))
# define STATIC_INLINE_BEGIN static inline
# if defined(__CYGWIN__)
#  define STATIC_INLINE_END
# else
#  define STATIC_INLINE_END __attribute__ ((always_inline))
# endif
#endif

#define __BEGIN__       {
#define __END__         goto exit; exit: ; }

// Thread-local storage
#if defined(_WIN32) && defined(_MSC_VER)
# define __THREAD __declspec(thread)
#elif defined(__GNUC__)
# define __THREAD __thread
#endif

#ifdef __GNUG__
# define GCC_VERSION (__GNUC__ * 10000 \
		      + __GNUC_MINOR__ * 100		\
		      + __GNUC_PATCHLEVEL__)
#endif

// TODO: Make file with helper macros
#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

