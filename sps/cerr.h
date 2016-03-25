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
