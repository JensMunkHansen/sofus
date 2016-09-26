/**
 * @file   strace.h
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Fri Jul 18 14:44:43 2014
 *
 * @brief  Stack-trace functionality
 *
 *
 */

#ifndef _STRACE_H_
#define _STRACE_H_

#include <sps/strace.hpp>
#include <sps/sps_export.h> // Not used (only POSIX)

#ifdef __cplusplus
extern "C" {
#endif

// This can be avoided if only singletons are supported
typedef struct STrace STrace;

/**
 * Create singleton
 *
 *
 * @return
 */
STrace* SPS_DLL_EXPORT strace_create();

/**
 * Enable stack trace
 *
 * @param obj
 *
 * @return
 */
straceErrorCodes SPS_DLL_EXPORT strace_enable(STrace* obj);

/**
 * Disable stack trace
 *
 * @param obj
 *
 * @return
 */
straceErrorCodes SPS_DLL_EXPORT strace_disable(STrace* obj);

/**
 * Set options
 *
 * @param obj
 *
 * @return
 */
straceErrorCodes SPS_DLL_EXPORT strace_set_opt(STrace* obj, straceOption opt, int val);

#ifdef __cplusplus
}
#endif

#endif /* _STRACE_H_ */
