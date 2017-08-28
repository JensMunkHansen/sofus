#ifndef __IF_FNM_H
#define __IF_FNM_H

#include <fnm/config.h>
#include <fnm/fnm_export.h>
#include <fnm/fnm_types.h>
#include <stddef.h>

// TODO: Issue that complex.h includes sps/memory rather than memory!!!
#include <complex.h>

#define FNM_EXTERNAL_API FNM_EXPORT

#ifdef __cplusplus
extern "C" {
#endif

// Enough to write aperture instead of struct aperture
typedef struct Aperture Aperture;

/**
 * Constructor
 *
 *
 * @return error code
 */
FNM_EXTERNAL_API int ApertureCreate(Aperture** obj);

/**
 *
 *
 * @param obj
 * @param nElements
 * @param width
 * @param kerf
 * @param height
 *
 * @return error code
 */
FNM_EXTERNAL_API int ApertureLinearCreate(Aperture** obj, size_t nElements, float width, float kerf, float height);

/**
 * Destructor
 *
 * @param obj object
 */
FNM_EXTERNAL_API int ApertureDestroy(Aperture* obj);

#ifdef FNM_CLOSURE_FUNCTIONS
/**
 * ApertureRwFloatParamSet
 *
 * @param obj Aperture object
 * @param fsel float parameter selector
 * @param f Input data
 * @param nDim Number of dimensions
 *
 * @return error code
 *
 * Variable arguments are used for specifying each dimensions, e.g.
 * @code
 * Aperture* obj;
 * ApertureLinearCreate(&obj,10,1.0,0.0,1.0);
 * int fsel = ElementDelays;
 * float* f = (float*) malloc(10*sizeof(float));
 * ApertureRwFloatParamSet(obj, fsel, f, 1, 10);
 * @code
 */
FNM_EXTERNAL_API int ApertureRwFloatParamSet(Aperture* obj, int fsel, const float* f, size_t nDim, ...);

/**
 * ApertureRwFloatParamGet
 *
 * @param obj Aperture object
 * @param fsel float parameter selector
 * @param f Output data
 * @param nDim
 *
 * @return error code
 *
 * When getting parameters, memory is allocated and data are copied.
 */
FNM_EXTERNAL_API int ApertureRwFloatParamGet(Aperture* obj, int fsel, float** f, size_t nDim, ...);
#endif

FNM_EXTERNAL_API int ApertureNDivWSet(Aperture* obj, size_t nDiv);

FNM_EXTERNAL_API int ApertureNDivHSet(Aperture* obj, size_t nDiv);

FNM_EXTERNAL_API int ApertureNThreadsSet(Aperture* obj, size_t nthreads);

FNM_EXTERNAL_API int ApertureFocusingTypeSet(Aperture* obj, int ftype);

#ifdef _MSC_VER
FNM_EXTERNAL_API int ApertureCalcCwFieldRef(Aperture* obj, float* pos, size_t nPositions, size_t nDim,
																						_Fcomplex** odata, size_t* nOutPositions);

FNM_EXTERNAL_API int ApertureCalcCwFieldFast(Aperture* obj, float* pos, size_t nPositions, size_t nDim,
																						 _Fcomplex** odata, size_t* nOutPositions);
#else
FNM_EXTERNAL_API int ApertureCalcCwFieldRef(Aperture* obj, float* pos, size_t nPositions, size_t nDim,
																						float _Complex** odata, size_t* nOutPositions);

FNM_EXTERNAL_API int ApertureCalcCwFieldFast(Aperture* obj, float* pos, size_t nPositions, size_t nDim,
																						 float _Complex** odata, size_t* nOutPositions);
#endif
	
#ifdef __cplusplus
}
#endif

#endif // #ifndef __IF_FNM_H

/* Local variables: */
/* indent-tab-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
