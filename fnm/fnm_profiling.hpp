/**
 * @file   fnm_profiling.hpp
 * @author Jens Munk Hansen <jmh@debian9laptop.parknet.dk>
 * @date   Tue Feb 20 21:18:08 2018
 *
 * @brief
 *
 *
 */

#pragma once

#include <fnm/fnm_export.h>

FNM_EXPORT void ProfilerStart();
FNM_EXPORT void ProfilerStop();
FNM_EXPORT void ProfileInfoGet(char** ostring);
FNM_EXPORT void Cancel();

extern char PerformanceInfo[32];  // = "\n";
extern double g_tStart;  // = 0.0;

#ifndef SWIG_VERSION
// TODO(JEM): SWIG on Windows does not like when this is an std::atomic<bool>
// extern std::atomic<bool> g_bExit;
extern volatile bool g_bExit;
#endif

namespace fnm {

}  // namespace fnm
