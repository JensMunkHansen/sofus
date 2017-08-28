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

#include <fnm/fnm_transient.hpp>
#include <fnm/fnm_data.hpp>
#include <fnm/fnm.hpp>            // Aperture<T>::fs
#include <fnm/fnm_common.hpp>

#include <sps/debug.h>

namespace fnm {

  template <class T>
  T (*const ToneBurst<T>::SpatialBasisFunction[ToneBurst<T>::nTerms])(T,T,T) = {
    ToneBurst<T>::sbf0,
    ToneBurst<T>::sbf1,
  };

  template <class T>
  T (*const ToneBurst<T>::TemporalBasisFunction[ToneBurst<T>::nTerms])(T,T,T) = {
    ToneBurst<T>::tbf0,
    ToneBurst<T>::tbf1,
  };

  template <class T>
  T (*const HanningWeightedPulse<T>::SpatialBasisFunction[HanningWeightedPulse<T>::nTerms])(T,T,T) = {
    HanningWeightedPulse<T>::sbf0,
    HanningWeightedPulse<T>::sbf1,
    HanningWeightedPulse<T>::sbf2,
    HanningWeightedPulse<T>::sbf3,
    HanningWeightedPulse<T>::sbf4,
    HanningWeightedPulse<T>::sbf5
  };

  template <class T>
  T (*const HanningWeightedPulse<T>::TemporalBasisFunction[HanningWeightedPulse<T>::nTerms])(T,T,T) = {
    HanningWeightedPulse<T>::tbf0,
    HanningWeightedPulse<T>::tbf1,
    HanningWeightedPulse<T>::tbf2,
    HanningWeightedPulse<T>::tbf3,
    HanningWeightedPulse<T>::tbf4,
    HanningWeightedPulse<T>::tbf5
  };

#ifdef FNM_PULSED_WAVE
  
  template <class T>
  T TransientSingleRect(const sysparm_t<T>* sysparm,
                        const ApertureData<T>* data,
                        const T* pos, const size_t nPositions,
                        T** odata, size_t* nSamples)
  {
    // Reset output
    *odata = nullptr;
    *nSamples = 0;

    const T c  = sysparm->c;
    const T fs = sysparm->fs;
    const T W  = sysparm->w;

    const T f0 = data->m_f0;

    const T invfs = T(1.0) / fs;

    //const T eps = std::numeric_limits<T>::epsilon();

    sps::bbox_t<T> box;
    data->ExtentGet(box);

    debug_print("nPositions: %zu\n",nPositions);

    size_t nElements,nSubElements;

    const T* delays;
    const T* apodizations;
    const sps::element_rect_t<T>** elementsArray = nullptr;

    data->ApodizationsRefGet(&nElements, apodizations);
    data->DelaysRefGet(&nElements, delays);
    data->ElementsRefGet(&nElements, &nSubElements, elementsArray);

    SPS_UNREFERENCED_PARAMETERS(apodizations, delays);

    // TODO: Support more than one element with zero delay
    const sps::element_rect_t<T>& element = elementsArray[0][0];

    // Simple box surrounding the scatters
    sps::bbox_t<T> scatter_box = sps::bbox_t<T>();
    sps::compute_bounding_box3(pos, nPositions, &scatter_box);

    // Box surrounding single element
    sps::bbox_t<T> element_box = sps::bbox_t<T>();
    data->ElementExtentGet(0, element_box);

    T distNear;
    T distFar;

    // Uses the 8 corner points
    sps::dists_most_distant_and_closest(scatter_box,
                                        element_box,
                                        &distNear,
                                        &distFar);

    T tStart = (distNear / c);
    T tEnd   = (distFar / c) + W;

    T fSampleSignalStart = fs * tStart;

    int iSampleSignalStart = (int) floor(fSampleSignalStart);

    // Allocate output
    int _nSamples = (int) ceil(T(2.0) + fs*(tEnd - tStart));

    *odata = (T*) malloc(nPositions*_nSamples*sizeof(T));
    memset(*odata, 0, nPositions*_nSamples*sizeof(T));
    *nSamples = _nSamples;

    sofus::sysparm_t<T> timeDomainParm;
    timeDomainParm.fs = sysparm->fs;
    timeDomainParm.c  = sysparm->c;

    auto uweights = sps::deleted_aligned_array<T>();
    auto vweights = sps::deleted_aligned_array<T>();
    auto uxs      = sps::deleted_aligned_array<T>();
    auto vxs      = sps::deleted_aligned_array<T>();

    CalcWeightsAndAbcissae(sysparm, std::move(uxs), std::move(uweights),
                           std::move(vxs),std::move(vweights));

    GLQuad2D<T> uv;
    uv.u.xs     = uxs.get();
    uv.u.ws     = uweights.get();
    uv.u.nx     = sysparm->nDivW;
    uv.v.xs     = vxs.get();
    uv.v.ws     = vweights.get();
    uv.v.nx     = sysparm->nDivH;

    for (size_t iPosition = 0 ; iPosition < nPositions ; iPosition++) {
      sofus::proj_limit_dist_t<T> pld;

      sps::point_t<T> point;
      point[0] = pos[iPosition * 3];
      point[1] = pos[iPosition * 3 + 1];
      point[2] = pos[iPosition * 3 + 2];

      T delay = T(0.0);
      sofus::calcProjectionAndLimits(timeDomainParm, elementsArray[0][0],
                                     point, delay,
                                     &pld);

      T z = pld.dist2plane;

      int iSampleStart = (int)pld.fSampleStart + 1; // First non-zero value

      int iSampleStop  = (int) (pld.fSampleStop + W*fs) + 1; // Ends with a non-zero value

      T t1;

      // Todo add delay (or subtract t_offset)
      t1 = z / c; // = tau_2

      // Solve integral using Gauss-Legendre
      T spatial = FourDirect(element, pld, &uv);

      // Missing factor (rho_0 c a)/pi
      for (int iSample = iSampleStart ; iSample < iSampleStop ; iSample++) {

        int iSampleSignal = iSample - iSampleSignalStart;
        T t = invfs * iSample;

        T direct = - ToneBurst<T>::eval(t-t1, W, f0) * spatial; // -v(t-tau_2)

        //TODO: Compute edge terms

        // update output
        T value = direct;
        (*odata)[iSampleSignal + iPosition*_nSamples] = value;
      }
    }

    return T(0.0);
  }

  template <class T>
  T FourEdge(const sysparm_t<T>* sysparm,
             const sps::element_rect_t<T>& element,
             const sofus::proj_limit_dist_t<T>& limit,
             const T& tau,
             const GLQuad2D<T>* uv);

  template <class T>
  T FourDirect(const sps::element_rect_t<T>& element,
               const sofus::proj_limit_dist_t<T>& limit,
               const GLQuad2D<T>* uv)
  {

    T result = T(0.0);

    // Integration limits
    T uLow  = fabs(limit.u) - element.hw;
    T uHigh = fabs(limit.u) + element.hw;

    T vLow  = fabs(limit.v) - element.hh;
    T vHigh = fabs(limit.v) + element.hh;

    // Adjacent edges
    T v[2];
    v[0] = element.hh - fabs(limit.v); // -vLow
    v[1] = element.hh + fabs(limit.v); //  vHigh

    T u[2];
    u[0] = element.hw - fabs(limit.u); // -uLow
    u[1] = element.hw + fabs(limit.u); //  uHigh

    T sm = T(0.5) * (uHigh - uLow);
    T sp = T(0.5) * (uHigh + uLow);

    T l2[2];
    l2[0] = SQUARE(v[0]);
    l2[1] = SQUARE(v[1]);

    T sInt[2];
    memset(sInt,0,2*sizeof(T));

    for (size_t iU = 0 ; iU < uv->u.nx ; iU++) {
      T s = sm * uv->u.ws[iU] + sp;
      T s2 = SQUARE(s);
      sInt[0] += T(1.0) / (l2[0]+s2);
      sInt[1] += T(1.0) / (l2[1]+s2);
    }

    result += sInt[0]*v[0];
    result += sInt[1]*v[1];

    sm = T(0.5) * (vHigh - vLow);
    sp = T(0.5) * (vHigh + vLow);

    l2[0] = SQUARE(u[0]);
    l2[1] = SQUARE(u[1]);

    memset(sInt,0,2*sizeof(T));

    for (size_t iV = 0 ; iV < uv->v.nx ; iV++) {
      T s = sm * uv->v.ws[iV] + sp;
      T s2 = SQUARE(s);
      sInt[0] += T(1.0) / (l2[0]+s2);
      sInt[1] += T(1.0) / (l2[1]+s2);
    }

    result += sInt[0]*u[0];
    result += sInt[1]*u[1];

    return result;
  }

  template <class T>
  void DirectWave(const sysparm_t<T>* sysparm,
                  const sps::element_rect_t<T>& element,
                  const sofus::proj_limit_dist_t<T>& limit,
                  const T f0,
                  const T delay,
                  const T apodization,
                  const GLQuad2D<T>* uv,
                  sps::msignal1D<T>* output)
  {
    SPS_UNREFERENCED_PARAMETER(delay);
    const T fs    = sysparm->fs;
    const T c     = sysparm->c;
    const T W     = sysparm->w;
    const T _invfs = T(1.0) / fs;

    int iSampleStart = ((int)limit.fSampleStart) + 1;

    int iSampleSignalStart = iSampleStart - output->offset;

    int nSamples = (int)floor(W*fs) + 1U;

    int iSampleSignalStop = iSampleSignalStart + nSamples;

    // Spatial part using 4 integrals. Missing a factor (\rho c) / (2\pi)

    T spatial = FourDirect(element, limit, uv);

    for (int i = iSampleSignalStart ; i < iSampleSignalStop ; i++) {
      int iSample = i + output->offset;
      T t = _invfs * iSample;
      (*output)[i] += apodization * spatial * HanningWeightedPulse<T>::eval(t - limit.dist2plane / c, W, f0);
    }
  }

  /**
   * Compute 6 edge-terms
   *
   * @param s0           sigma start
   * @param s1           sigma stop
   * @param l            adjacent
   * @param z            z-coordinate
   * @param ss           abcissas
   * @param sweights     weights
   * @param nSs          # of abcissas
   * @param f0           f0
   * @param W            width of pulse [s]
   * @param output       edge terms
   */
  template <class T>
  STATIC_INLINE_BEGIN
  void EdgeTerms(const T& s0, const T& s1,
                 const T& l, const T& z,
                 const T* ss,
                 const T* sweights,
                 const size_t nSs,
                 const T& f0,
                 const T& W,
                 T output[6])
  {

    const T sm = T(0.5) * (s1 - s0);
    const T sp = T(0.5) * (s1 + s0);

    const T z2 = SQUARE(z);
    const T l2 = SQUARE(l);

    // tau, W, f0
    memset(output, 0, 6*sizeof(T));

    for (size_t iS = 0 ; iS < nSs ; iS++) {

      T s = sm * ss[iS] + sp;
      T s2 = SQUARE(s);

      T tau = sqrt(z2+s2+l2);
      T denom = std::max<T>((s2+l2),std::numeric_limits<T>::epsilon());
      for (size_t iTerm = 0 ; iTerm < 6 ; iTerm++) {
        output[iTerm] += HanningWeightedPulse<T>::SpatialBasisFunction[iTerm](tau, W, f0)/denom;
      }
    }

    // Scale: rho*l / (2pi)
    for (size_t iTerm = 0 ; iTerm < 6 ; iTerm++) {
      output[iTerm] = output[iTerm] * l / T(T(M_2PI));
    }
  }

  template <class T>
  void EdgeWaves(const sysparm_t<T>* sysparm,
                 const sps::element_rect_t<T>& element,
                 const sofus::proj_limit_dist_t<T>& limit,
                 const T f0,
                 const T delay,
                 const T apodization,
                 const GLQuad2D<T>* uv,
                 sps::msignal1D<T>* output)
  {
    const T c     = sysparm->c;
    const T fs    = Aperture<T>::fs;
    const T _invfs = T(1.0) / fs;

    // First non-zero sample
    int iSampleStart = ((int)limit.fSampleStart) + 1;

    // First non-zero sample in signal
    int iSampleSignalStart = iSampleStart - output->offset;

    // Number of extra samples due to pulse
    int nTemporalSamples = ((int) (limit.fSampleStart + sysparm->W*fs)) + 1 - iSampleStart;

    // Number of spatio-temporal samples
    int nSpatioTemporalSamples =  ((int) (limit.fSampleStart + sysparm->W*fs + (limit.fSampleStop - limit.fSampleStart))) + 1 - iSampleStart;

    int iSampleSignalStop = iSampleSignalStart + nSpatioTemporalSamples;

    T spatial = T(0.0);

    // realdist = ((sample * invfs) - delay)*c
    // nTemporalSamples (*) when to start changing sigma
    // nSpatialSamples

    for (int i = iSampleSignalStart ; i < iSampleSignalStop ; i++) {
      int iSample = i + output->offset;

      T tau = iSample * _invfs - delay;

      // tau_0 = fSampleStart * invfs - delay
      // tau_d2p = limit.dist2plane / c

      // T t = _invfs * iSample; // Time beginning of pulse reaches point (if no delays)

      // Find limits of integration and compute spatial part
      spatial = FourEdge<T>(sysparm, element, limit, tau, uv);

      //(*output)[i] += spatial * HanningWeightedPulse<T>::eval(t - limit.dist2plane / c, sysparm->W, f0);
    }
  }

  template <class T>
  T FourEdge(const sysparm_t<T>* sysparm,
             const sps::element_rect_t<T>& element,
             const sofus::proj_limit_dist_t<T>& limit,
             const T& tau,
             const GLQuad2D<T>* uv)
  {
    SPS_UNREFERENCED_PARAMETERS(element,limit,tau,uv);

    //

    T tau_d2p = limit.dist2plane / sysparm->c;
    //
    // tau_0, tau_1, tau_2, tau_3

    T tau_1;
    T W;

    //if (t-(tau_1+0.5*W))

    // sigma_0 = sqrt(SQUARE(t-W) - SQUARE(l))

    T result = T(0.0);
    return result;
  }

  template <class T>
  T CalcFdTransientRef(const sysparm_t<T>* sysparm,
                       const ApertureData<T>* data,
                       const T* pos, const size_t nPositions, const size_t nDim,
                       T** odata, size_t* nSignals, size_t* nSamples)
  {
    SPS_UNREFERENCED_PARAMETERS(nDim);
    const T eps    = std::numeric_limits<T>::epsilon();
    const T fs = sysparm->fs;
    const T _invfs = T(1.0) / fs;
    SPS_UNREFERENCED_PARAMETER(_invfs);

    sofus::sysparm_t<T> timeDomainParm;
    timeDomainParm.fs = sysparm->fs;
    timeDomainParm.c  = sysparm->c;

    // Pulse is 10 samples long
    const T f0 =  data->m_f0;
    sps::bbox_t<T> box;
    data->ExtentGet(box);

    debug_print("nPositions: %zu\n",nPositions);

    int iSampleStart = 0;
    T tStart = T(0.0);
    size_t nData;
    size_t nElements,nSubElements;

    const T* delays;
    const T* apodizations;
    const sps::element_rect_t<T>** elementsArray = nullptr;

    data->ApodizationsRefGet(&nElements, apodizations);
    data->DelaysRefGet(&nElements, delays);
    data->ElementsRefGet(&nElements, &nSubElements, elementsArray);

    ComputeBoxTimes(timeDomainParm,
                    box,
                    pos, nPositions,
                    delays, apodizations, nElements,
                    &tStart,
                    &iSampleStart,
                    &nData);

    size_t nTemporalData = nData + (size_t) floor(sysparm->w*fs) + 1U;

    // Allocate and zero output
    T* _data = (T*) malloc(nPositions*nTemporalData*sizeof(T));
    memset(_data,0,nPositions*nTemporalData*sizeof(T));

    auto plds = sps::deleted_aligned_multi_array_create<sofus::proj_limit_dist_t<T>,2>(nElements,nSubElements);

    T apodization;
    T delay;

    auto uweights = sps::deleted_aligned_array<T>();
    auto vweights = sps::deleted_aligned_array<T>();
    auto uxs      = sps::deleted_aligned_array<T>();
    auto vxs      = sps::deleted_aligned_array<T>();

    CalcWeightsAndAbcissae(sysparm, std::move(uxs), std::move(uweights),
                           std::move(vxs),std::move(vweights));

    GLQuad2D<T> uv;
    uv.u.xs     = uxs.get();
    uv.u.ws     = uweights.get();
    uv.u.nx     = sysparm->nDivW;
    uv.v.xs     = vxs.get();
    uv.v.ws     = vweights.get();
    uv.v.nx     = sysparm->nDivH;

    sps::msignal1D<T> signal(nTemporalData, (nData+3)*sizeof(T));
    signal.offset = iSampleStart;

    // Loop over positions
    for (size_t iPoint = 0 ; iPoint < nPositions ; iPoint++) {

      sps::point_t<T> point;
      point[0] = pos[iPoint*3 + 0];
      point[1] = pos[iPoint*3 + 1];
      point[2] = pos[iPoint*3 + 2];

      T fSampleStartMin = std::numeric_limits<T>::max();
      T fSampleStopMax  = std::numeric_limits<T>::min();

      memset((void*)signal.get(),0,signal.nbytes);

      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {

        apodization = apodizations[iElement];

        if (fabs(apodization) > eps) {
          for (size_t iSubElement = 0 ; iSubElement < nSubElements ; iSubElement++) {
            sofus::calcProjectionAndLimits(timeDomainParm, elementsArray[iElement][iSubElement],
                                           point, delays[iElement],
                                           &(plds[iElement][iSubElement]));
            fSampleStopMax  = std::max<T>(fSampleStopMax, plds[iElement][iSubElement].fSampleStop);
            fSampleStartMin = std::min<T>(fSampleStartMin, plds[iElement][iSubElement].fSampleStart);
          }
        }
      }

      int iSampleStartMin = ((int)fSampleStartMin) + 1; // First non-zero value
      int iSampleStopMax  = ((int)fSampleStopMax) + 1;  // Ends with a non-zero value

      int nXmtSamples = iSampleStopMax - iSampleStartMin;

      SPS_UNREFERENCED_PARAMETERS(nXmtSamples, iSampleStartMin, iSampleStopMax);
      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {

        apodization = apodizations[iElement];
        delay = delays[iElement];
        SPS_UNREFERENCED_PARAMETERS(apodization, delay);
        if (fabs(apodization) > eps) {
          for (size_t iSubElement = 0 ; iSubElement < nSubElements ; iSubElement++) {
            const sps::element_rect_t<T>& element = elementsArray[iElement][iSubElement];
            const sofus::proj_limit_dist_t<T>& limit = plds[iElement][iSubElement];
            SPS_UNREFERENCED_PARAMETERS(element,limit);

            DirectWave(sysparm, element, limit,
                       f0,
                       delay, apodization,
                       &uv,
                       &signal);

            // Compute loop (8 integrals, update 5 terms)

            // Make single integral function

            // Find du and dv

            // Integrate one of four corners (spatial term)

            // Multi with temporal

            // Save
          }
        }
      }
      memcpy(&_data[iPoint*nTemporalData],
             signal.get(), signal.ndata*sizeof(T));

    } // positions
    // Find limits

    // Find sigma limits

    // Perform integrals (4)

    // Multiply with temporal functions

    // spatial
    // cos(2 pi f0 tau)
    // sin(2 pi f0 tau)
    // cos(2 pi tau / W) cos(2 pi f0 tau)
    // cos               sin
    // sin               cos
    // sin               sin

    // temporal
    // 0.5 sin(2 pi f0 t)
    // -0.5 cos(2 pi f0 t)
    // -0.5 cos(2 pi t/W) sin(2 pi f0 t)
    //  0.5 cos           cos
    // -0.5 sin           sin
    //  0.5 sin           cos

    // (exp(-jk sqrt(z**2+sigma**2+s**2)) - exp(-jkz)) / (sigma**2 + s**2)
    //
    // -> (v(t-1/c sqrt(z**2+sigma**2 + s**2)) - v(t-z/c))  / (sigma**2 + s**2)
    //
    // v(t-tau) = rect((t-tau)/W) sum f_n(tau) g_n(t)

    *nSamples = nTemporalData;
    *nSignals = nPositions;
    *odata    = _data;

    return 0;
  }
#endif

  template struct ToneBurst<float>;

#ifdef FNM_PULSED_WAVE
  template float
  FNM_EXPORT CalcFdTransientRef(const fnm::sysparm_t<float>* sysparm,
                                const ApertureData<float>* data,
                                const float* pos, const size_t nPositions, const size_t nDim,
                                float** odata, size_t* nSignals, size_t* nSamples);
  
  template float TransientSingleRect(const sysparm_t<float>* sysparm,
                                     const ApertureData<float>* data,
                                     const float* pos, const size_t nPositions,
                                     float** odata, size_t* nSamples);
#endif



#if FNM_DOUBLE_SUPPORT
  template struct ToneBurst<double>;

#ifdef FNM_PULSED_WAVE
  template double
  FNM_EXPORT CalcFdTransientRef(const fnm::sysparm_t<double>* sysparm,
                                const ApertureData<double>* data,
                                const double* pos, const size_t nPositions, const size_t nDim,
                                double** odata, size_t* nSignals, size_t* nSamples);

  template double TransientSingleRect(const sysparm_t<double>* sysparm,
                                      const ApertureData<double>* data,
                                      const double* pos, const size_t nPositions,
                                      double** odata, size_t* nSamples);
#endif

#endif
}
