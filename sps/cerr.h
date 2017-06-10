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

#include <sps/stdio.h>
#include <sps/cenv.h>
#include <string.h> // strerror
#include <errno.h>


/**
 * Macro for calling functions.
 *
 * @param fun
 * @param arg
 *
 */
#define CallErr(fun, arg)  {                                  \
  if ((fun arg)<0)                                            \
    FailErr(#fun)                                             \
}

/**
 * Macro for calling functions. Functions registered using atexit will
 * be called before the process exists
 *
 * @param fun
 * @param arg
 *
 * @return
 */
#define CallErrExit(fun, arg) {                               \
  if ((fun arg)<0) {                                          \
    FailErr(#fun);                                            \
    exit(EXIT_FAILURE);                                       \
  }                                                           \
}

/**
 * Macro for calling functions. If function fails (retval < 0), the
 * process will be terminated without calling functions registered
 * using atexit.
 *
 * @param fun
 * @param arg
 *
 * @return
 */
#define CallErr_Exit(fun, arg) {                              \
  if ((fun arg)<0) {                                          \
    FailErr(#fun);                                            \
    _exit(EXIT_FAILURE);                                      \
  }                                                           \
}

#define FailErr(msg) {                                        \
    (void)fprintf(stderr, "%s, %s(%d)\n",                     \
                msg, strerror(errno), errno);                 \
    (void)fflush(stderr);}

#define CallErrReturn(fun, arg, ret)  {                       \
  if ((fun arg)<0) {                                          \
    FailErr(#fun);                                            \
    return ret;                                               \
  }                                                           \
}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
