/**
 * @file   stdio.h
 * @author Jens Munk Hansen <jmh@jmhlaptop.parknet.dk>
 * @date   Tue Jun 23 22:40:52 2015
 *
 * @brief  Interface to stdio.h. Takes care of missing C99 functionality on Windows (snprintf)
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

#include <stdio.h>

#ifndef C99
#include <stdarg.h>

#define snprintf c99_snprintf

#ifdef _MSC_VER
inline int c99_vsnprintf(char* str, size_t size, const char* format, va_list ap)
{
  int count = -1;

  if (size != 0)
    count = _vsnprintf_s(str, size, _TRUNCATE, format, ap);
  if (count == -1)
    count = _vscprintf(format, ap);

  return count;
}

inline int c99_snprintf(char* str, size_t size, const char* format, ...)
{
  int count;
  va_list ap;

  va_start(ap, format);
  count = c99_vsnprintf(str, size, format, ap);
  va_end(ap);

  return count;
}

#elif defined(__CYGWIN__)
// TODO: Handle C99 on cygwin
#else
# error No compatible snprintf function
#endif

#endif // C99
