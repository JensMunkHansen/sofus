#ifndef SPS_PROFILER_H
#define SPS_PROFILER_H

#include <sps/cenv.h>

#ifdef _MSC_VER
# include <windows.h>
#else
# include <sys/time.h>
# include <sys/resource.h>
#endif

namespace sps {
  class profiler {
  public:
#ifdef _MSC_VER
    static double time(void)
    {
      static __THREAD double freq = 0.0;
      if (freq == 0.0) {
        LARGE_INTEGER llDiv;
        QueryPerformanceFrequency(&llDiv);
        freq = double(llDiv.QuadPart);
      }
      LARGE_INTEGER ll;
      QueryPerformanceCounter(&ll);
      return double(ll.QuadPart) / freq;
    }
#else
    static double time(void)
    {

      register double sec;
      register double usec;

#if 1
      struct timespec tspec;
      clock_gettime(CLOCK_MONOTONIC, &tspec);
      sec = (double) tspec.tv_sec;
      usec = ((double) tspec.tv_nsec) / 1000000000.0;
#else
      static __THREAD struct rusage foo;
      getrusage (RUSAGE_SELF, &foo); // process values
      sec = (double) foo.ru_utime.tv_sec;
      usec = ((double) foo.ru_utime.tv_usec) / 1000000.0;
#endif
      return (sec + usec);
    }
#endif
  };
}



#endif // SPS_PROFILER_H

/*

clock() from <time.h> (20ms or 10ms resolution?)

gettimeofday() from Posix <sys/time.h> (microseconds)

clock_gettime() on Posix (nanoseconds?) // CLOCK_PROCESS_CPUTIME_ID, CLOCK_THREAD_CPUTIME_ID,
CLOCK_MONOTONIC, CLOCK_REALTIME
*/
