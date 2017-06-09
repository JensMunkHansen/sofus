/**
 * @file   strace.h
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Fri Jul 18 14:44:43 2014
 *
 * @brief  Stack-trace functionality (C-interface)
 *
 *
 */

#ifndef _STRACE_H_
#define _STRACE_H_

#ifndef SPS_EXPORT
# include <sps/sps_export.h>
#endif

#include <sps/strace.hpp>

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
STrace* SPS_EXPORT strace_create();

/**
 * Enable stack trace
 *
 * @param obj
 *
 * @return
 */
straceErrorCodes SPS_EXPORT strace_enable(STrace* obj);

/**
 * Disable stack trace
 *
 * @param obj
 *
 * @return
 */
straceErrorCodes SPS_EXPORT strace_disable(STrace* obj);

/**
 * Set options
 *
 * @param obj
 *
 * @return
 */
straceErrorCodes SPS_EXPORT strace_set_opt(STrace* obj, straceOption opt, int val);

#ifdef __cplusplus
}
#endif

#endif /* _STRACE_H_ */
