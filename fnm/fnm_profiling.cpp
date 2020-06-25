/**
 * @file   fnm_profiling.cpp
 * @author Jens Munk Hansen <jmh@debian9laptop.parknet.dk>
 * @date   Tue Feb 20 21:17:45 2018
 *
 * @brief
 *
 *
 */

#include <sps/profiler.h>
#include <sps/malloc.h>
#include <sps/mm_malloc.h>

#include <cstdio>
#include <cstring>
#include <fnm/fnm_profiling.hpp>

// TODO(JEM): Consider moving them into namespace + static
char PerformanceInfo[32] = "\n";
double g_tStart = 0.0;
volatile bool g_bExit(false);

void Cancel() {
  g_bExit = true;
}

void ProfilerStart() {
  g_tStart = sps::profiler::time();
  g_bExit = false;
}

void ProfilerStop() {
  double duration = sps::profiler::time() - g_tStart;
  sprintf(PerformanceInfo, "%3.2f seconds", duration);
}

void ProfileInfoGet(char** ostring) {
  *ostring = static_cast<char*>(SPS_MALLOC(32));
  strncpy(*ostring, PerformanceInfo, 32);
}
