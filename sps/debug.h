/**
 * @file   debug.h
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Sun Jul 12 16:16:54 2015
 *
 * @brief  Debug macros
 *
 *
 */
#pragma once

#include <sps/stdio.h>
#include <sps/cenv.h>

#ifndef SPS_DEBUG
# define SPS_DEBUG 0
#endif

#ifndef __SPS_FUNCTION__
# ifdef _MSC_VER
#  define __SPS_FUNCTION__ __FUNCTION__
# else
#  define __SPS_FUNCTION__ __FUNCTION__
# endif
#endif

#ifndef __SPS_DEBUG_FILE__
# define __SPS_DEBUG_FILE__ 0
#endif

#ifndef __SPS_DEBUG_LINE__
# define __SPS_DEBUG_LINE__ 1
#endif

#ifndef __SPS_DEBUG_FUNCTION__
# define __SPS_DEBUG_FUNCTION__ 1
#endif

#define MULTI_LINE_MACRO_BEGIN do {
#define MULTI_LINE_MACRO_END \
    __pragma(warning(push)) \
    __pragma(warning(disable: 4127)) \
    } while(0) \
    __pragma(warning(pop))

// Compiler bug in MSVC. Disabling C4127 is not working for a do { } while (0)
#if defined(_WIN32)
# if (__SPS_DEBUG_FUNCTION__) && (__SPS_DEBUG_LINE__) && (__SPS_DEBUG_FILE__)
#  define debug_print(fmt, ...)
__pragma(warning(push)) \
__pragma(warning(disable: 4127)) \
for(;;)
{
  if (SPS_DEBUG) fprintf(stderr, "%s:%s():%d: " fmt, __FILE__, __SPS_FUNCTION__, __LINE__,  __VA_ARGS__);
}
break;
} \
__pragma(warning(pop))
# else
#  define debug_print(fmt,...)                                                               \
   __pragma(warning(push)) \
   __pragma(warning(disable: 4127)) \
   __pragma(warning(disable: 6271)) \
   for(;;) { if (SPS_DEBUG) { fprintf(stderr, "%s(): " fmt, __SPS_FUNCTION__,  __VA_ARGS__);} break; } \
   __pragma(warning(pop))
# endif
#elif defined(C99)
# if (__SPS_DEBUG_FUNCTION__) && (__SPS_DEBUG_LINE__) && (__SPS_DEBUG_FILE__)
#  define debug_print(...)                                                        \
  do { if (SPS_DEBUG) fprintf(stderr, "%s:%d:%s(): " FIRST(__VA_ARGS__), __FILE__,        \
                          __SPS_FUNCTION__, __LINE__  REST_(__VA_ARGS__)); } while (0)
# else
#  define debug_print(...) \
   do { if (SPS_DEBUG) fprintf(stderr, "%s(): " FIRST(__VA_ARGS__), __SPS_FUNCTION__ REST(__VA_ARGS__)); } while (0)
# endif
#else
# if (__SPS_DEBUG_FUNCTION__) && (__SPS_DEBUG_LINE__) && (__SPS_DEBUG_FILE__)
#  define debug_print(fmt, ...)                                                   \
  do { if (SPS_DEBUG) fprintf(stderr, "%s:%s():%d: " fmt, __FILE__,                   \
                          __SPS_FUNCTION__, __LINE__,  ## __VA_ARGS__); } while (0)
# else
#  define debug_print(fmt, ...)                                                                \
  do { if (SPS_DEBUG) fprintf(stderr, "%s(): " fmt, __SPS_FUNCTION__, ## __VA_ARGS__); } while (0)
# endif
#endif

// MULTI_LINE_MACRO_BEGIN if (DEBUG) fprintf(stderr, "%s(): " fmt, __SPS_FUNCTION__, __VA_ARGS__); MULTI_LINE_MACRO_END

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
