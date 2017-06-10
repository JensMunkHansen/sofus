/**
 * @file   strace.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Sun Jul 20 04:43:06 2014
 *
 * @brief  Stack-trace functionality (C++-interface)
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

#ifndef _STRACE_HPP_
#define _STRACE_HPP_

// TODO: Use fcntl(fd, F_GETFD) for testing file descriptor

#include <cstdio>
#include <cstddef>

// Forward declarations
#if 0
typedef struct siginfo siginfo_t;
#else
# include <sps/signal.h>
#endif

typedef struct sigaction sigaction_t;

#ifndef _STRACE_H_
namespace sps {

  class STrace {
  public:
#endif

    /**
     * Options which are all configurable using integers
     *
     */
    typedef enum {
      STRACE_GEN_CORE_DUMP ,  ///< Generate core dump (default=1)
      STRACE_EXEC_EXIT_FUNCS, ///< Call function registered using atexit
      STRACE_FAST_EXIT,       ///< Exit fast. If this is set, clean up and core dump are ignored (default=0)
      STRACE_STACK_DEPTH,     ///< Depth of stack (frames) (default=16)
      STRACE_COMMON_PATH,     ///< Include the common path for executables
      STRACE_REL_PATH,        ///< Relative path
      STRACE_PID,             ///< Include PID
      STRACE_COLOR_OUTPUT,    ///< Enable color output
      STRACE_THREAD_SAFETY,   ///< Enable thread safety
      STRACE_FILEDESCRIPTOR,  ///< File descriptor (default=STDERR_FILENO)
      STRACE_OPTION_COUNT     ///< Number of options
    } straceOption;

    typedef enum {
      STRACE_ERR_OK                  = 0,
      STRACE_ERR_INVALID_STACK_DEPTH,    ///<
      STRACE_ERR_INVALID_FD,             ///<  Invalid file descriptor
      STRACE_ERR_COUNT
    } straceErrorCodes;

#ifndef _STRACE_H_

    /**
     * Print to file descriptor
     *
     * @param msg
     * @param len
     */
    static void print2fd(const char *msg, size_t len = 0);

    /**
     * Acquire singleton instance
     *
     *
     * @return
     */
    static STrace& Instance();

    /**
     * Enable stack trace (current handler are stored)
     *
     *
     * @return
     */
    sps::STrace::straceErrorCodes enable();

    /**
     * Disable stack trace (old handlers are restored)
     *
     *
     * @return
     */
    sps::STrace::straceErrorCodes disable();

    /**
     * Set options
     *
     * @param opt
     * @param val
     *
     * @return
     */
    sps::STrace::straceErrorCodes set_opt(straceOption opt, int val);

  private:
    STrace();                                   // Private so that it can  not be called
    ~STrace();
    STrace(STrace const&) {} //           = default;  // copy constructor is private
    STrace& operator=(STrace const&);//= default;  // assignment operator is private

    int init();

    static bool m_bInitialized; ///< Is it initialized

    static struct sigaction m_sa_abrt; ///< Handler for abortion
    static struct sigaction m_sa_segv; ///< Handler for segmentation fault

    /**
     * The signal handler. It is very limited stuff which can be used
     * inside this function, e.g. no malloc's
     *
     * @param sig
     * @param info
     * @param secret
     */
    static void signalHandler(int sig, ::siginfo_t* info, void* secret);

    /**
     * Malloc hook used for preventing backtrace for allocating during
     * the handling of the signals
     *
     * @param size
     * @param caller
     *
     * @return
     */
    static void* mallocHook(size_t size, const void* caller);

    /**
     * Wrapper around the program addr2line (Debian).
     *
     * @param image
     * @param addr
     * @param color_output
     * @param memory
     *
     * @return
     */
    static char* addr2line(const char *image, void *addr,
                           bool color_output,
                           char** memory);

    static bool m_bGenerateCoreDump; ///< Generate a core dump on exit
    static bool m_bCallAtExit;       ///< Call functions registered using atexit
    static bool m_bQuickExit;        ///< Call quickexit()
    static size_t m_nFrames;         ///< Number of stack frames
    static bool m_bIncCommonPath;    ///< Include common path
    static bool m_bIncRelativePaths; ///< Include relative paths, e.g ../
    static bool m_bAppendPID;        ///< Append process id's
    static bool m_bEnableColorOutput;///< Use color for output
    static bool m_bThreadSafe;       ///< Thread-safe handling of signals (verify SIGCONT)
    static int  m_fdOutput;          ///< File descriptor used for output (default=STDOUT_FILENO)
    static char* memory_;            ///< Heap used for operation

    static const size_t kNeededMemory;
  };
}
#endif

#endif /* _STRACE_HPP_ */

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
