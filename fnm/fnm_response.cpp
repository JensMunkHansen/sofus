/**
 * @file   fnm_response.cpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Tue Sep 26 19:42:59 2017
 *
 * @brief
 *
 *
 */

#include <fnm/config.h>
#include <fnm/fnm_response.hpp>
#include <fnm/fnm_basis.hpp>
#include <fnm/fnm_data.hpp>

// TODO: Move stuff to fnm_response.hpp
#include <fnm/fnm_transient.hpp>

#if 1 //def FNM_PULSED_WAVE
# include <sofus/rect_int_limits.hpp> // calcProjectionAndLimits
# include <sps/msignals.hpp>
#endif

namespace fnm {
/*
  T offset = tStart (sampled) - delay;
  iSample = int(fSampleStart - offset*fs) + 1; // first non-zero sample
  t = invfs * iSample + offset; // fractional time to use

  DirectWave needs to know: iSample, tStart, delay

  and increment

  invfs * iSample + offset

  by

  invfs
 */


template <typename T, template <typename> class A>
void FnmResponse(const sysparm_t<T>* pSysparm,
                 const ApertureData<T>* pData,
                 const GLQuad2D<T>* pGL,
                 const T& amplitude,
                 const sps::point_t<T>& point,
                 const int& iSampleSignalStart,
                 const size_t& nSamples,
                 T* pOdata,
                 const int mask) {
  const sps::element_rect_t<T>** ppElements = NULL;
  const T* pApodizations = nullptr;
  const T* pDelays = nullptr;

  size_t nElements = 0, nSubElements = 0;
  pData->ElementsRefGet(&nElements, &nSubElements, ppElements);
  pData->ApodizationsRefGet(&nElements, pApodizations);
  pData->DelaysRefGet(&nElements, pDelays);

  sofus::sysparm_t<T> timeDomainParm;
  timeDomainParm.fs = pSysparm->fs;
  timeDomainParm.c  = pSysparm->c;

  sofus::proj_limit_dist_t<T> pld;

  T delay;
  T apodization;

  for (size_t iElement = 0 ; iElement < nElements ; iElement++) {
    apodization = pApodizations[iElement];
    delay       = pDelays[iElement];
    if (fabs(apodization) > T(0.0)) {
      for (size_t jElement = 0 ; jElement <  nSubElements ; jElement++) {
        const sps::element_rect_t<T>& element = ppElements[iElement][jElement];

        assert(sps::is_aligned<4*sizeof(T)>::value(&element.center[0]));
        assert(sps::is_aligned<4*sizeof(T)>::value(&element.uvector[0]));
        assert(sps::is_aligned<4*sizeof(T)>::value(&element.vvector[0]));

        // TODO(Double): Fix for double precision
        sofus::calcProjectionAndLimits(timeDomainParm, element,
                                       point, delay,
                                       &pld);
        debug_print("pld.u: %f, pld.v: %f, pld.dist2plane: %f,"
                    " fSampleStart: %f, fSampleStop: %f\n",
                    pld.u, pld.v,
                    pld.dist2plane, pld.fSampleStart, pld.fSampleStop);

        T scale = amplitude * apodization;

        // Compute response
        if (mask & 0x01) {
          DirectWaveSingle<T, A>(
            pSysparm, &element, pGL, scale, &pld,
            delay,
            iSampleSignalStart, nSamples, pOdata);
        }
        if (mask & 0x02) {
          EdgeResponse<T, A>(
            pSysparm, &element, scale, 0, &pld,
            &pGL->v,
            delay,
            iSampleSignalStart, nSamples, pOdata);
        }
        if (mask & 0x04) {
          EdgeResponse<T, A>(
            pSysparm, &element, scale, 1, &pld,
            &pGL->u,
            delay,
            iSampleSignalStart, nSamples, pOdata);
        }
        if (mask & 0x08) {
          EdgeResponse<T, A>(
            pSysparm, &element, scale, 2, &pld,
            &pGL->v,
            delay,
            iSampleSignalStart, nSamples, pOdata);
        }
        if (mask & 0x10) {
          EdgeResponse<T, A>(
            pSysparm, &element, scale, 3, &pld,
            &pGL->u,
            delay,
            iSampleSignalStart, nSamples, pOdata);
        }
      }
    }
  }
}

template void
FnmResponse<float, ToneBurst>(
    const sysparm_t<float>* pSysparm,
    const ApertureData<float>* pData,
    const GLQuad2D<float>* pGL,
    const float& amplitude,
    const sps::point_t<float>& point,
    const int& iSampleSignalStart,
    const size_t& nSamples,
    float* pOdata,
    const int mask);

template void
FnmResponse<float, HanningWeightedPulse>(
  const sysparm_t<float>* pSysparm,
  const ApertureData<float>* pData,
  const GLQuad2D<float>* pGL,
  const float& amplitude,
  const sps::point_t<float>& point,
  const int& iSampleSignalStart,
  const size_t& nSamples,
  float* pOdata,
  const int mask);

template void
FnmResponse<float, Identity>(
  const sysparm_t<float>* pSysparm,
  const ApertureData<float>* pData,
  const GLQuad2D<float>* pGL,
  const float& amplitude,
  const sps::point_t<float>& point,
  const int& iSampleSignalStart,
  const size_t& nSamples,
  float* pOdata,
  const int mask);

#ifdef FNM_DOUBLE_SUPPORT
template void FnmResponse<double, ToneBurst>(
  const sysparm_t<double>* pSysparm,
  const ApertureData<double>* pData,
  const GLQuad2D<double>* pGL,
  const double& amplitude,
  const sps::point_t<double>& point,
  const int& iSampleSignalStart,
  const size_t& nSamples,
  double* pOdata,
  const int mask);

template void FnmResponse<double, HanningWeightedPulse>(
  const sysparm_t<double>* pSysparm,
  const ApertureData<double>* pData,
  const GLQuad2D<double>* pGL,
  const double& amplitude,
  const sps::point_t<double>& point,
  const int& iSampleSignalStart,
  const size_t& nSamples,
  double* pOdata,
  const int mask);

template void FnmResponse<double, Identity>(
  const sysparm_t<double>* pSysparm,
  const ApertureData<double>* pData,
  const GLQuad2D<double>* pGL,
  const double& amplitude,
  const sps::point_t<double>& point,
  const int& iSampleSignalStart,
  const size_t& nSamples,
  double* pOdata,
  const int mask);
#endif
}  // namespace fnm


/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
