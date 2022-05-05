/**
 * @file   if_fnm.h
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Fri Sep  8 21:14:16 2017
 *
 * @brief  External C interface for creation of apertures and computation
 *         of continous-wave fields.
 *
 */

#ifndef __IF_FNM_H
#define __IF_FNM_H

#include <fnm/config.h>
#include <fnm/fnm_export.h>

#include <stddef.h>
#include <complex.h>
#include <stdarg.h>

#ifndef SPS_FCOMPLEX
# include <sps/cenv.h>
#endif

#ifndef FNM_EXTERNAL_API
# define FNM_EXTERNAL_API FNM_EXPORT
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct Aperture Aperture;

/**
 * Constructor
 *
 *
 * @return error code
 */
FNM_EXTERNAL_API int ApertureCreate(Aperture** obj);

/**
 * Create Linear Array. @see fnm::FocusedLinearArray
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


// FNM_EXTERNAL_API int ApertureRwFloatParamSet1D(Aperture* obj, int fsel, const float* f, size_t nData);


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


FNM_EXTERNAL_API int FreeCArray(void* pData);

/**
 * ApertureNDivWSet
 *
 * See @ref fnm::Aperture<T>::ApertureNDivWSet
 *
 * @param obj
 * @param nDiv
 *
 * @return
 */
FNM_EXTERNAL_API int ApertureNDivWSet(Aperture* obj, size_t nDiv);

/**
 * ApertureNDivWSet
 *
 * See @ref fnm::Aperture<T>::ApertureNDivWSet
 *
 * @param obj
 *
 * @return
 */
FNM_EXTERNAL_API size_t ApertureNDivWGet(Aperture* obj);

/**
 *
 *
 * @param obj
 * @param nDiv
 *
 * @return
 */
FNM_EXTERNAL_API int ApertureNDivHSet(Aperture* obj, size_t nDiv);

/**
 *
 *
 * @param obj
 *
 * @return
 */
FNM_EXTERNAL_API size_t ApertureNDivHGet(Aperture* obj);

/**
 *
 *
 * @param obj
 *
 * @return
 */
FNM_EXTERNAL_API size_t ApertureNThreadsGet(Aperture* obj);

/**
 *
 *
 * @param[in] obj
 * @param[in] nThreads
 *
 * @return
 */
FNM_EXTERNAL_API int ApertureNThreadsSet(Aperture* obj, size_t nThreads);

/**
 *
 *
 * @param obj
 * @param ftype
 *
 * @return
 */
FNM_EXTERNAL_API int ApertureFocusingTypeSet(Aperture* obj, int ftype);

/**
 *
 *
 * @param obj
 * @param pos
 * @param nPositions
 * @param nDim
 * @param odata
 * @param nOutPositions
 *
 * @return
 */
FNM_EXTERNAL_API int ApertureCalcCwFieldRef(Aperture* obj, float* pos, size_t nPositions, size_t nDim,
    SPS_FCOMPLEX** odata, size_t* nOutPositions);

/**
 *
 *
 * @param obj
 * @param pos
 * @param nPositions
 * @param nDim
 * @param odata
 * @param nOutPositions
 *
 * @return
 */
FNM_EXTERNAL_API int ApertureCalcCwFieldFast(Aperture* obj, float* pos, size_t nPositions, size_t nDim,
    SPS_FCOMPLEX** odata, size_t* nOutPositions);

#ifdef __cplusplus
}
#endif

#endif // #ifndef __IF_FNM_H

/* Local variables: */
/* indent-tab-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */

