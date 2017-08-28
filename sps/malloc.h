/**
 * @file   malloc.h
 * @author Jens Munk Hansen <jmh@jmhlaptop>
 * @date   Thu Aug  3 20:11:01 2017
 *
 * @brief  wrapper for malloc.h
 *
 *
 */

#pragma once

#include <sps/cenv.h>
#include <malloc.h>

/**
 * Return size of allocated memory
 *
 * @param data
 *
 * @return
 */
STATIC_INLINE_BEGIN size_t msize(void* data)
{
#ifdef _WIN32
  return _msize(data);
#elif __APPLE__
  return malloc_size(data);
#elif defined(__GNUG__)
  return malloc_usable_size(data);
#endif
}

