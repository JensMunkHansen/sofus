/**
 * @file   stdlib.h
 * @author Jens Munk Hansen <jmh@jmhlaptop.parknet.dk>
 * @date   Mon Sep 14 19:13:35 2015
 *
 * @brief  Wrapper for stdlib
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

#ifdef HAVE_CONFIG
# include <sps/config.h>
#endif

#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include <sps/cenv.h>

namespace sps {
// Consider adding a control size argument instead
  const int xtoa_max_length = 32;

#ifndef HAVE_ITOA
  STATIC_INLINE_BEGIN char *itoa(int val, char* memory, int base = 10)
  {
    char* res = memory;
    int i;
    bool negative = val < 0;
    res[xtoa_max_length - 1] = '\0';
    for (i = xtoa_max_length - 2; val != 0 && i != 0; i--, val /= base) {
      res[i] = "0123456789ABCDEF"[abs(val) % base];
    }
    if (negative) {
      res[i--] = '-';
    }
    return &res[i + 1];
  }
#endif

#ifndef HAVE_UTOA
  STATIC_INLINE_BEGIN char *utoa(uintptr_t val, char* memory, int base = 10)
  {
    char* res = memory;
    int i;
    res[xtoa_max_length - 1] = '\0';
    for (i = xtoa_max_length - 2; val != 0 && i != 0; i--, val /= base) {
      res[i] = "0123456789ABCDEF"[val % base];
    }
    return &res[i + 1]; // Start of string
  }
#endif

#ifndef HAVE_PTOA
  STATIC_INLINE_BEGIN char *ptoa(const void *val, char* memory)
  {
    const int res_max_length = 32;
    // New area, utoa take at most 32 chars (incl null termination)
    char* buf = sps::utoa(reinterpret_cast<uintptr_t>(val), memory + res_max_length, 16);
    char* result = memory;
    strcpy(result + 2, buf);
    result[0] = '0';
    result[1] = 'x';
    return result;
  }
#endif
}
