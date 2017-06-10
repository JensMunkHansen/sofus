/**
 * @file   strace.h
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Fri Jul 18 14:44:43 2014
 *
 * @brief  Stack-trace functionality (C-interface)
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
