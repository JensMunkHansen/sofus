#pragma once

// Assume pthreads available
#include <sps/cerr.h>
#include <sps/config.h>

#ifdef _WIN32
# ifndef NOMINMAX
#  define NOMINMAX
# endif
# include <windows.h>
# include <process.h>
# include <io.h>
# if defined(HAVE_PTHREAD_H)
#  include <pthread.h>
# endif
#else
# if defined(HAVE_SIGNAL_H)
#  include <csignal>
# endif
# if defined(HAVE_PTHREAD_H)
#  include <pthread.h>
# endif
# include <unistd.h>
# include <sys/syscall.h>
#endif

#ifdef __linux__
static pid_t gettid( void ) {
    pid_t pid;
    CallErr(pid = syscall, (__NR_gettid));
    return pid;
}
#endif

static int setcpuid(int cpu_id) {
#ifdef __GNUG__
  cpu_set_t set;
  CPU_ZERO( &set );
  CPU_SET( cpu_id, &set );
  CallErrReturn(sched_setaffinity,
                (gettid(), sizeof(cpu_set_t), &set), -1);
#elif defined(_WIN32)
  HANDLE hThread = GetCurrentThread();
  SetThreadIdealProcessor(hThread,cpu_id);
#endif
  return 0;
}

static int getncpus() {
  int nproc;
#if defined(__linux__)
  cpu_set_t cpuset;
  CPU_ZERO( &cpuset );
  CallErr(sched_getaffinity,
          (gettid(), sizeof( cpu_set_t ), &cpuset));
  
  nproc = CPU_COUNT(&cpuset);
#elif defined(_WIN32)
  SYSTEM_INFO info;
  GetSystemInfo(&info);
  nproc = info.dwNumberOfProcessors;
#endif
  return nproc;
}

namespace sps {

#if defined(HAVE_PTHREAD_H)
  template<class T, void*(T::*thread_func)(void*)>
#elif defined(_WIN32)
  template<class T, unsigned int(__stdcall T::*thread_func)(void*)>
#endif
  
  class pthread_launcher {
  public:
    pthread_launcher(T* obj=NULL, void* arg=NULL) : _obj(obj), _arg(arg) {}
#if defined(HAVE_PTHREAD_H)
    void *launch() { return (_obj->*thread_func)(_arg);}
#elif defined(_WIN32)
    unsigned int launch() { return (_obj->*thread_func)(_arg);}
#endif
  private:
    /// Object pointer
    T* _obj;
    /// Command argument for member function
    void *_arg;
  };
  
  // Launch thread function
  template<class T>
#if defined(HAVE_PTHREAD_H)
  void *launch_member_function(void *obj)
#elif defined(_WIN32)
    unsigned int __stdcall launch_member_function(void *obj)
#endif
  {
    T* launcher = reinterpret_cast<T*>(obj);
    return launcher->launch();
  }
}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
