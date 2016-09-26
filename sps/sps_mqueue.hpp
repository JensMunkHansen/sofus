/**
 * @file   sps_mqueue.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Thu Jun 16 21:21:58 2016
 *
 * @brief
 *
 *
 */

#pragma once

#include <sps/config.h>
#include <sps/sps_export.h>
#include <sps/cenv.h>

#ifdef HAVE_MQUEUE_H
# include <mqueue.h>
#endif

/**
 * Clear POSIX message queue
 *
 * @param qname name of the queue
 *
 * @return
 */
SPS_EXPORT int mq_clear(const char* qname);

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
