#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <fnm/fnm_export.h>
#include <fnm/fnm_types.h>
#include <stddef.h>

#define FNM_EXTERNAL_API FNM_EXPORT

// Enough to write aperture instead of struct aperture
typedef struct Aperture Aperture;

/**
 * Constructor
 *
 *
 * @return
 */
FNM_EXTERNAL_API int ApertureCreate(Aperture** obj);
