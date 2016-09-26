/**
 * @file   strace.cpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Sun Jul 20 03:19:30 2014
 *
 * @brief
 *
 *
 */

/* TODO: Study strace, which intercepts system calls, ltrace, which
   intercepts library calls and ptrace. ltrace only suports ELF32
   binaries and fails to tracks some children, when using the -f
   option. Also libraries opened using dlopen is not traced. */

// Capture all system calls using ptrace

#include <sps/strace.hpp>

#include <sps/cerr.h>
#include <sps/stdlib.h>

#include <cstring> // strlen, strcat, strcmp, strcpy
#include <stdint.h>
#include <malloc.h> // __malloc_hook

#include <cxxabi.h>
#include <execinfo.h> // backtrace

#include <signal.h>
#include <wait.h>
#include <unistd.h>
#include <fcntl.h>
#include <pthread.h>

#ifndef _GNU_SOURCE
# define _GNU_SOURCE
#endif
#include <dlfcn.h>

// We are not allowed to do malloc inside signal handlers
#pragma GCC poison malloc realloc free backtrace_symbols        \
  printf fprintf sprintf snprintf scanf sscanf

#ifndef STATIC_INLINE
# define STATIC_INLINE __attribute__((always_inline)) static inline
#endif

static inline void safe_abort()
{
  // Read signal handler for abort
  struct sigaction sa;
  sigaction(SIGABRT, NULL, &sa);

  // Insert default handler
  sa.sa_handler = SIG_DFL;

  // Clear SIG_INFO flag (calls sa_sigaction)
  sa.sa_flags = sa.sa_flags & ~SA_SIGINFO;

  // Tell process to continue in case children are hanging
  kill(getppid(), SIGCONT);

  // Install handler
  sigaction(SIGABRT, &sa, NULL);

  // Abort process
  abort();
}

#if __GNUC__ >= 4 && __GNUC_MINOR__ >= 6
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif

namespace sps {

//  STrace* STrace::m_pInstance = NULL; // Old style
  bool    STrace::m_bInitialized = false;

  bool   STrace::m_bGenerateCoreDump = false;
  bool   STrace::m_bCallAtExit       = true;
  bool   STrace::m_bQuickExit        = false;
  size_t STrace::m_nFrames           = 16U;
  bool   STrace::m_bIncCommonPath    = false;
  bool   STrace::m_bIncRelativePaths = false;
  bool   STrace::m_bAppendPID        = false;
  bool   STrace::m_bEnableColorOutput= true;
  bool   STrace::m_bThreadSafe       = false;
  int    STrace::m_fdOutput          = STDERR_FILENO;

// REMOVE
  char*  STrace::memory_             = NULL;

  const size_t STrace::kNeededMemory = 12288;

  struct sigaction STrace::m_sa_segv;
  struct sigaction STrace::m_sa_abrt;

  void STrace::print2fd(const char *msg, size_t len)
  {
    if (len > 0) {
      CallErr_Exit(write, (m_fdOutput, msg, len));
    } else {
      CallErr_Exit(write, (m_fdOutput, msg, strlen(msg)));
    }
  }

  STrace::STrace()
  {
    if (memory_ == NULL) {
      memory_ = new char[kNeededMemory];
    }
  }

  STrace::~STrace()
  {
    if (memory_) {
      delete [] memory_;
    }
  }

/// @brief Used to workaround backtrace() usage of malloc().
  void* STrace::mallocHook(size_t size,
                           const void* caller )
  {
    char* malloc_buffer = memory_ + kNeededMemory - 512;
    if (size > 512U) {
      const char* msg = "malloc() replacement function should not return "
                        "a memory block larger than 512 bytes\n";
      print2fd(msg, strlen(msg) + 1);
      _exit(EXIT_FAILURE);
    }
    return malloc_buffer;
  }

  STrace& STrace::Instance()
  {
    static STrace singleton;
    return singleton;
  }

  sps::STrace::straceErrorCodes STrace::set_opt(straceOption opt, int val)
  {

    sps::STrace::straceErrorCodes err = sps::STrace::straceErrorCodes(0);
    switch (opt) {
    case (STRACE_FILEDESCRIPTOR):
      if (fcntl(val,F_GETFL,0)) {
        err = STRACE_ERR_INVALID_FD;
        break;
      }
      m_fdOutput = val;
      break;
    case (STRACE_GEN_CORE_DUMP):
      m_bGenerateCoreDump = (bool) val;
      break;
    case (STRACE_EXEC_EXIT_FUNCS):
      m_bCallAtExit = (bool) val;
      break;
    case (STRACE_FAST_EXIT):
      m_bQuickExit = (bool) val;
      break;
    case (STRACE_STACK_DEPTH):
      m_nFrames = (size_t) val;
      break;
    case (STRACE_COMMON_PATH):
      m_bIncCommonPath = (bool) val;
      break;
    case (STRACE_REL_PATH):
      m_bIncRelativePaths = (bool) val;
      break;
    case (STRACE_PID):
      m_bAppendPID = (bool) val;
      break;
    case (STRACE_COLOR_OUTPUT):
      m_bEnableColorOutput = (bool) val;
      break;
    default:
      break;
    }
    return err;
  }

  void STrace::signalHandler(int sig, ::siginfo_t* info, void* secret)
  {

    // Stop all other running threads by forking
    pid_t forkedPid = fork();
    if (forkedPid != 0) {
      // Main process
      int status;
      if (m_bThreadSafe) {
        // Freeze the original process, until it's child prints the stack trace
        kill(getpid(), SIGSTOP);
        // Wait for the child without blocking and exit as soon as possible,
        // so that no zombies are left.
        waitpid(forkedPid, &status, WNOHANG);
      } else {
        // Wait for the child, blocking only the current thread.
        // All other threads will continue to run, potentially crashing the parent.
        waitpid(forkedPid, &status, 0);
      }
      if (m_bQuickExit) {
        ::quick_exit(EXIT_FAILURE);
      }
      if (m_bGenerateCoreDump) {
        struct sigaction sa;
        sigaction(SIGABRT, NULL, &sa);
        sa.sa_handler = SIG_DFL;
        sigaction(SIGABRT, &sa, NULL);
        abort();
      } else {
        if (m_bCallAtExit) {
          exit(EXIT_FAILURE);
        } else {
          _exit(EXIT_FAILURE);
        }
      }
    } else {
      // Child process - the else clause is for clarity

      // TODO: Consider using atexit for closing fd

      ucontext_t *uc = reinterpret_cast<ucontext_t *>(secret);

      if (dup2(STDERR_FILENO, STDOUT_FILENO) == -1) {  // redirect stdout to stderr
        print2fd("Failed to redirect stdout to stderr\n");
      }
      char* memory = memory_;
      {
        char* msg = memory;
        const int msg_max_length = 128;
        if (m_bEnableColorOutput) {
          // Bold red
          strcpy(msg, "\033[31;1m");
        } else {
          msg[0] = '\0';
        }
        switch (sig) {
        case SIGSEGV:
          strcat(msg, "Segmentation fault");
          break;
        case SIGABRT:
          strcat(msg, "Aborted");
          break;
        default:
          strcat(msg, "Caught signal ");
          strcat(msg, sps::itoa(sig, msg + msg_max_length));
          break;
        }
        if (m_bEnableColorOutput)
          strcat(msg, "\033[0m");

        strcat(msg, " (thread");

        if (m_bEnableColorOutput)
          strcat(msg, "\033[33;1m"); // yellow bold

        // If a pthread
        strcat(msg, sps::utoa(pthread_self(), msg + msg_max_length));

        if (m_bEnableColorOutput)
          strcat(msg, "\033[0m");

        strcat(msg, ", pid ");

        if (m_bEnableColorOutput)
          strcat(msg, "\033[33;1m");  // yellow bold

        strcat(msg, sps::itoa(getppid(), msg + msg_max_length));

        if (m_bEnableColorOutput)
          strcat(msg, "\033[0m");

        strcat(msg, ")");
        print2fd(msg);
      }

      print2fd("\n\nStack trace:\n");
      void **trace = reinterpret_cast<void**>(memory);
      memory += (m_nFrames + 2) * sizeof(void*);

      // Workaround malloc() inside backtrace()
      void* (*oldMallocHook)(size_t, const void*) = __malloc_hook;
      void (*oldFreeHook)(void *, const void *)   = __free_hook;
      __malloc_hook = mallocHook;
      __free_hook = NULL;
      int trace_size = backtrace(trace, m_nFrames + 2);
      __malloc_hook = oldMallocHook;
      __free_hook = oldFreeHook;
      if (trace_size <= 2) {
        close(m_fdOutput);
        safe_abort();
      }

      // Overwrite sigaction with caller's address
#if defined(__arm__)
      trace[1] = reinterpret_cast<void *>(uc->uc_mcontext.arm_pc);
#elif defined(__x86_64__)
      trace[1] = reinterpret_cast<void *>(uc->uc_mcontext.gregs[REG_RIP]);
#elif defined(__i386__)
      trace[1] = reinterpret_cast<void *>(uc->uc_mcontext.gregs[REG_EIP]);
#else
# error Only ARM, x86 and x86-64 are supported
#endif

      const int path_max_length = 2048;
      char* name_buf = memory;
      ssize_t name_buf_length = readlink("/proc/self/exe", name_buf,
                                         path_max_length - 1);
      if (name_buf_length < 1) {
        close(m_fdOutput);
        safe_abort();
      }
      name_buf[name_buf_length] = 0;
      memory += name_buf_length + 1;
      char* cwd = memory;
      if (getcwd(cwd, path_max_length) == NULL) {
        close(m_fdOutput);
        safe_abort();
      }
      strcat(cwd, "/");
      memory += strlen(cwd) + 1;
      char* prev_memory = memory;

      int stackOffset = trace[2] == trace[1]? 2 : 1;
      for (int i = stackOffset; i < trace_size; i++) {
        memory = prev_memory;
        char *line;
        Dl_info dlinf;
        if (dladdr(trace[i], &dlinf) == 0 || dlinf.dli_fname[0] != '/' ||
            !strcmp(name_buf, dlinf.dli_fname)) {
          line = addr2line(name_buf, trace[i], m_bEnableColorOutput, &memory);
        } else {
          line = addr2line(dlinf.dli_fname,
                           reinterpret_cast<void *>(reinterpret_cast<char *>(trace[i]) -
                               reinterpret_cast<char *>(dlinf.dli_fbase)),
                           m_bEnableColorOutput, &memory);
        }

        char *function_name_end = strstr(line, "\n");
        if (function_name_end != NULL) {
          *function_name_end = 0;
          {
            // "\033[34;1m[%s]\033[0m \033[33;1m(%i)\033[0m\n
            char* msg = memory;
            const int msg_max_length = 512;
            if (m_bEnableColorOutput) {
              strcpy(msg, "\033[34;1m"); // bold blue
            } else {
              msg[0] = 0;
            }
            strcat(msg, "[");
            strcat(msg, line);
            strcat(msg, "]");
            if (m_bAppendPID) {
              if (m_bEnableColorOutput)
                strcat(msg, "\033[33;1m"); // bold yellow
              //                strcat(msg, "\033[0m\033[33;1m"); // normal bold yellow

              strcat(msg, " (");
              strcat(msg, sps::itoa(getppid(), msg + msg_max_length));
              strcat(msg, ")");
              if (m_bEnableColorOutput)
                strcat(msg, "\033[0m");

              strcat(msg, "\n");
            } else {
              if (m_bEnableColorOutput)
                strcat(msg, "\033[0m");

              strcat(msg, "\n");
            }
            print2fd(msg);
          }
          line = function_name_end + 1;

          // Remove the common path root
          if (!m_bIncCommonPath) {
            int cpi;
            for (cpi = 0; cwd[cpi] == line[cpi]; cpi++) {};
            if (line[cpi - 1] != '/') {
              for (; line[cpi - 1] != '/'; cpi--) {};
            }
            if (cpi > 1) {
              line = line + cpi;
            }
          }

          // Remove relative path root
          if (!m_bIncRelativePaths) {
            char *path_cut_pos = strstr(line, "../");
            if (path_cut_pos != NULL) {
              path_cut_pos += 3;
              while (!strncmp(path_cut_pos, "../", 3)) {
                path_cut_pos += 3;
              }
              line = path_cut_pos;
            }
          }

          // Mark line number
          if (m_bEnableColorOutput) {
            char* number_pos = strstr(line, ":");
            if (number_pos != NULL) {
              char* line_number = memory;  // 128
              strcpy(line_number, number_pos);
              // Overwrite the new line char
              line_number[strlen(line_number) - 1] = 0;
              // \033[32;1m%s\033[0m\n
              strcpy(number_pos, "\033[32;1m"); // bold green
              strcat(line, line_number);
              strcat(line, "\033[0m\n");
            }
          }
        }

        // Overwrite the new line char
        line[strlen(line) - 1] = 0;

        // Append pid
        if (m_bAppendPID) {
          // %s\033[33;1m(%i)\033[0m\n
          strcat(line, " ");
          if (m_bEnableColorOutput)
            strcat(line, "\033[33;1m"); // bold yellow

          strcat(line, "(");
          strcat(line, sps::itoa(getppid(), memory));
          strcat(line, ")");
          if (m_bEnableColorOutput)
            strcat(line, "\033[0m");
        }

        strcat(line, "\n");
        print2fd(line);
      }

      // Write '\0' to indicate the end of the output
      char end = '\0';
      write(STDERR_FILENO, &end, 1);

      if (m_bThreadSafe) {
        // Resume the parent process
        kill(getppid(), SIGCONT);
      }

      close(m_fdOutput);

      // This is called in the child process
      _exit(EXIT_SUCCESS);
    }
  }

  sps::STrace::straceErrorCodes STrace::enable()
  {

    if (!m_bInitialized) {
      this->init();
    }
    struct sigaction sa;
    //  sa.sa_handler = (__sighandler_t)signalHandler;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_RESTART | SA_SIGINFO | SA_ONSTACK;

    sa.sa_sigaction = signalHandler;
    sigaction(SIGSEGV, &sa, NULL);
    sigaction(SIGABRT, &sa, NULL);
    return sps::STrace::straceErrorCodes(EXIT_SUCCESS);
  }

  sps::STrace::straceErrorCodes STrace::disable()
  {
    if (m_bInitialized) {
      sigaction(SIGSEGV, &m_sa_segv, NULL);
      sigaction(SIGABRT, &m_sa_abrt, NULL);
    }
    return sps::STrace::straceErrorCodes(EXIT_SUCCESS);
  }

  int STrace::init()
  {

    if (!m_bInitialized) {
      // Store old state, once
      sigaction(SIGABRT, NULL, &m_sa_abrt);
      sigaction(SIGSEGV, NULL, &m_sa_abrt);
      m_bInitialized = true;
    }
    return EXIT_SUCCESS;
  }

// modify memory
  char* STrace::addr2line(const char *image /* input exe name */,
                          void *addr /* input */,
                          bool color_output /* input */,
                          char** memory)
  {

    // Pipe from child to parent process
    int pipefd[2];
    if (pipe(pipefd) != 0) {
      safe_abort();
    }
    pid_t pid = fork();

    if (pid == 0) {
      // Child
      close(pipefd[0]);
      dup2(pipefd[1], STDOUT_FILENO);
      dup2(pipefd[1], STDERR_FILENO);
      if (execlp("addr2line", "addr2line",
                 sps::ptoa(addr, *memory), "-f", "-C", "-e", image,
                 reinterpret_cast<void*>(NULL)) == -1) {
        // Take care if execlp fails
        close(pipefd[1]);
        safe_abort();
      }
    }

    // Only parent
    close(pipefd[1]);

    const int line_max_length = 4096;
    char* line = *memory;
    *memory += line_max_length;
    ssize_t len = read(pipefd[0], line, line_max_length);
    close(pipefd[0]);
    if (len == 0) {
      safe_abort();
    }
    line[len] = 0;

    // Wait for child
    if (waitpid(pid, NULL, 0) != pid) {
      safe_abort();
    }
    if (line[0] == '?') {
      char* straddr = sps::ptoa(addr, *memory);

      if (color_output)
        strcpy(line, "\033[32;1m");

      strcat(line, straddr);

      if (color_output)
        strcat(line, "\033[0m");

      strcat(line, " at ");
      strcat(line, image);
      strcat(line, " ");
    } else {
      if (*(strstr(line, "\n") + 1) == '?') {
        char* straddr = sps::ptoa(addr, *memory);
        strcpy(strstr(line, "\n") + 1, image);
        strcat(line, ":");
        strcat(line, straddr);
        strcat(line, "\n");
      }
    }
    return line;
  }

}


#if __GNUC__ >= 4 && __GNUC_MINOR__ >= 6
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif

// See also sysdig, systemtap(kernel level) and dtrace4linux

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
