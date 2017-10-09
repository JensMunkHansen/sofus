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

// TODO: Make this work for non-sampled delays, Evaluate()

#include <fnm/fnm_transient.hpp>
#include <fnm/fnm_data.hpp>
#include <fnm/fnm_common.hpp>
#include <fnm/fnm_basis.hpp>

#include <sps/debug.h>

#include <sofus/rect_int_limits.hpp> // calcProjectionAndLimits

#include <sps/algorithm>

#include <utility>

namespace fnm {

  template <class T>
  T FourDirect(const sps::element_rect_t<T>* element,
               const sofus::proj_limit_dist_t<T>* limit,
               const GLQuad2D<T>* uv)
  {
    const T eps = std::numeric_limits<T>::epsilon();

    T result = T(0.0);

    // Integration limits
    T uLow  = fabs(limit->u) - element->hw;
    T uHigh = fabs(limit->u) + element->hw;

    T vLow  = fabs(limit->v) - element->hh;
    T vHigh = fabs(limit->v) + element->hh;

    // Adjacent
    T v[2];
    v[0] = element->hh - fabs(limit->v);
    v[1] = element->hh + fabs(limit->v);

    T u[2];
    u[0] = element->hw - fabs(limit->u);
    u[1] = element->hw + fabs(limit->u);

    T sp;

    // Offset of abcissa
    sp = T(0.5) * (uHigh + uLow);

    T l2[2];
    l2[0] = SQUARE(v[0]);
    l2[1] = SQUARE(v[1]);

    T sInt[2];
    memset(sInt,0,2*sizeof(T));

    for (size_t iU = 0 ; iU < uv->u.nx ; iU++) {
      T s2 = SQUARE(uv->u.xs[iU] + sp);
      sInt[0] += uv->u.ws[iU] * (T(1.0) / std::max<T>(s2+l2[0],eps));
      sInt[1] += uv->u.ws[iU] * (T(1.0) / std::max<T>(s2+l2[1],eps));
    }
    result += sInt[0]*v[0];
    result += sInt[1]*v[1];

    // Offset of abcissa
    sp = T(0.5) * (vHigh + vLow);

    l2[0] = SQUARE(u[0]);
    l2[1] = SQUARE(u[1]);

    memset(sInt,0,2*sizeof(T));

    // Without this, 0x01 is not confined to inside the element
    for (size_t iV = 0 ; iV < uv->v.nx ; iV++) {
      T s2 = SQUARE(uv->v.xs[iV] + sp);
      sInt[0] += uv->v.ws[iV] * (T(1.0) / std::max<T>(s2+l2[0],eps));
      sInt[1] += uv->v.ws[iV] * (T(1.0) / std::max<T>(s2+l2[1],eps));
    }
    result += sInt[0]*u[0];
    result += sInt[1]*u[1];

    return result;
  }

  template <class T,  template <typename> class A>
  void DirectWaveSingle(const sysparm_t<T>* sysparm,
                        const sps::element_rect_t<T>* element,
                        const T& scale,
                        const T& f0,
                        const sofus::proj_limit_dist_t<T>* pld,
                        const GLQuad2D<T>* uv,
                        const T delay,
                        const int iSampleSignalStart,
                        const size_t nSamples,
                        T* odata)
  {
    const T W = sysparm->w;
    const T c = sysparm->c;
    const T fs = sysparm->fs;

    const T invfs = T(1.0) / fs;

    const T factor = - scale * sysparm->rho * c / T(M_2PI);

    SPS_UNREFERENCED_PARAMETER(nSamples);

    T z = pld->dist2plane;

    T t1;

    t1 = z / c; // tau_2

    int iSampleStart = (int) (pld->fSampleStart + delay*fs) + 1; // First non-zero value
    int iSampleStop  = (int) (pld->fSampleStop + delay*fs + W*fs) + 1; // Ends with a non-zero value

    T spatial = FourDirect(element, pld, uv);

    assert(iSampleStart > (iSampleSignalStart - 1) && "Access before start");
    assert((iSampleStop - 1 - iSampleSignalStart) < (int)nSamples && "Access after stop");

    // Vectorize four sample at a time + (nSamples % 4)
    for (int iSample = iSampleStart ; iSample < iSampleStop ; iSample++) {

      int iSampleSignal = iSample - iSampleSignalStart;
      T t = invfs * iSample - delay; // TEST: If we remove this, it gets crazy

      // TODO: Make this work with non-sampled delays
      T direct = - A<T>::Evaluate(t-t1, W, f0) * spatial;

      T value = direct;
      // HERE SEGV
      odata[iSampleSignal] += factor * value;
    }
  }

  template <typename T, template <typename> class A>
  void UpdateSIR(const sysparm_t<T>* sysparm,
                 const T f0,
                 const T ssigma, const T sweight,
                 const T adjacent, const T z2,
                 A<T>* pulse)
  {
    const T eps = std::numeric_limits<T>::epsilon();

    // Uses scaled weight and sigma
    T denom = std::max<T>(eps, SQUARE(adjacent) + SQUARE(ssigma));
    T tau = sqrt(std::max<T>(T(0.0),z2 + denom)) / sysparm->c;
    T scale = sweight / denom;

    pulse->UpdateSpatial(scale, tau, sysparm->w, f0);
  }

  // TODO: Vectorize - compute 4 edges at a time

  // Safe implementation using two integration ranges for each edge
  template <typename T, template <typename> class A>
  void EdgeResponse(const sysparm_t<T>* sysparm,
                    const sps::element_rect_t<T>* element,
                    const T& scale,
                    const size_t iEdge,
                    const sofus::proj_limit_dist_t<T>* pld,
                    const GLQuad1D<T>* pGL,
                    const T& f0,
                    const T& delay,
                    const int& iSampleSignalStart,
                    const size_t& nSamples,
                    T* odata)
  {
    //    const T eps = std::numeric_limits<T>::epsilon();

    const T fs = sysparm->fs;
    const T invfs = T(1.0) / fs;
    const T W = sysparm->w;
    const T c = sysparm->c;

    const T factor = - scale * sysparm->rho * c / T(M_2PI);

    const T z2 = SQUARE(pld->dist2plane);

    T hl = T(0.0); // Half long
    T pl = T(0.0); // Point long

    T hs = T(0.0); // Half short
    T ps = T(0.0); // Point short

    T adjacent = T(0.0); // Orthogonal to integration range

    size_t iVertices[2] = {0,0}; // Vertices for active edge

    size_t nl = 0; // Number of abcissas for integration along l

    switch (iEdge) {
    case 0:
      iVertices[0] = 0;
      iVertices[1] = 3;
      break;
    case 1:
      iVertices[0] = 0;
      iVertices[1] = 1;
      break;
    case 2:
      iVertices[0] = 1;
      iVertices[1] = 2;
      break;
    case 3:
      iVertices[0] = 3;
      iVertices[1] = 2;
      break;
    };

    switch (iEdge) {
    case 1:
    case 3:
      hl = element->hw;
      hs = element->hh;
      pl = pld->u;
      ps = pld->v;
      nl = sysparm->nDivH;
      break;
    case 0:
    case 2:
      hl = element->hh;
      hs = element->hw;
      pl = pld->v;
      ps = pld->u;
      nl = sysparm->nDivW;
      break;
    };

    switch (iEdge) {
    case 0:
    case 1:
      if (ps > T(0.0)) {
        adjacent = hs - fabs(ps);
      } else {
        adjacent = hs + fabs(ps);
      }
      break;
    case 2:
    case 3:
      if (ps < T(0.0)) {
        adjacent = hs - fabs(ps);
      } else {
        adjacent = hs + fabs(ps);
      }
      break;
    };

    T t1 = pld->vdists[iVertices[0]] / c;
    T t2 = pld->vdists[iVertices[1]] / c;

    T t12min, t12max;

    enum Segment {
      Lower = 0,
      Upper = 1,
    };

    // Limits - used for disabling upper or lower
    int iSigmaMin[2];
    int iSigmaMax[2];

    // Nothing is disabled for now
    iSigmaMin[Lower] = 0;
    iSigmaMax[Lower] = (int)nl;
    iSigmaMin[Upper] = 0;
    iSigmaMax[Upper] = (int)nl;

    // Dynamic limits
    int iSigmaLow[2];
    int iSigmaUp[2];

    // Reset to zero
    memset(iSigmaLow, 0, 2*sizeof(int));
    memset(iSigmaUp, 0, 2*sizeof(int));

    if (fabs(pl) < hl) {
      t12min = sqrt(std::max<T>(T(0.0), z2 + SQUARE(adjacent))) / c;

      auto end = pGL->xs + nl;

      // Greater or equal (We should end here, verify)
      auto it = std::find_if(pGL->xs, end,[&](T a)->bool {return !(pl > a);});

      if (it != end) {
        int iProj = (int) std::distance(pGL->xs, it);
        iSigmaLow[Lower] = iProj;
        iSigmaUp[Lower]  = iProj;
        iSigmaLow[Upper] = iProj;
        iSigmaUp[Upper]  = iProj;
      } else {
        iSigmaLow[Lower] = (int)nl;
        iSigmaUp[Lower] = (int)nl;
        // Upper integral is disabled usig iSigmaMax
        iSigmaLow[Upper] = (int)nl;
        iSigmaUp[Upper]  = (int)nl;
      }
    } else {
      t12min = std::min<T>(t1,t2);
      if (pl > hl) {
        iSigmaLow[Lower] = (int)nl;
        iSigmaUp[Lower]  = (int)nl;
        // Upper integral is disabled
        iSigmaMax[Upper] = 0;
      } else {
        iSigmaLow[Upper] = 0;
        iSigmaUp[Upper]  = 0;
        // Lower integral is disabled
        iSigmaMin[Lower] = (int) nl;
      }
    }

    t12max = std::max<T>(t1,t2);

    T fSampleStart = fs * t12min;
    T fSampleStop = fs * t12max;

    int iSampleStart = (int)(fSampleStart + delay*fs) + 1;
    int iSampleStop  = (int) (fSampleStop + delay*fs + W*fs) + 1;

    // Lower values for ranges: [-hl ; pl] and [pl, hl]
    T fSigmaMin[2];

    // Upper values for ranges: [-hl ; pl] and [pl, hl]
    T fSigmaMax[2];

    T deltaSigma;

    T sigma;

    A<T> ImpulseResponse = A<T>();

#if SPS_DEBUG && 0
    typedef std::vector<std::pair<int,int> > limitVector;
    limitVector upperLimits = limitVector();
    limitVector lowerLimits = limitVector();
#endif

#ifdef SPS_DEBUG
    if (iSampleStart < iSampleSignalStart) {
      debug_print("\033[31;1mSamples too early: u: %f, v: %f\033[0m\n", pld->u, pld->v);
    }
    if (iSampleStop >= int(nSamples + iSampleSignalStart)) {
      debug_print("\033[31;1mSamples too late: u: %f, v: %f\033[0m\n", pld->u, pld->v);
    }
#endif

    assert(iSampleStart > (iSampleSignalStart - 1) && "Access before start");
    assert((iSampleStop - 1 - iSampleSignalStart) < (int)nSamples && "Access after stop");

    for (int iSample = iSampleStart ; iSample < iSampleStop ; iSample++) {

      int iSampleSignal = iSample - iSampleSignalStart;
      T t = invfs * iSample - delay;

      // Compute two sigma segments (asymmetry could be if we start to early)
      T sqrtArg = SQUARE(c*t) - z2 - SQUARE(adjacent);
      sqrtArg = std::max<T>(T(0.0), sqrtArg);
      deltaSigma = sqrt(sqrtArg);

      // Lower sigma range
      fSigmaMax[Lower] = std::clamp(pl, -hl, hl);
      fSigmaMin[Lower] = std::clamp(pl - deltaSigma, -hl, hl);

      // Upper sigma range
      fSigmaMax[Upper] = std::clamp(pl + deltaSigma, -hl, hl);
      fSigmaMin[Upper] = std::clamp(pl, -hl, hl);

      if ((t > t12min + W) && (t < t12max + W)) {
        // Shorten ranges
        deltaSigma = sqrt(std::max<T>(T(0.0), SQUARE(c*(t-W)) - z2 - SQUARE(adjacent)));
        fSigmaMin[Upper] = std::clamp(pl + deltaSigma, -hl, hl);
        fSigmaMax[Lower] = std::clamp(pl - deltaSigma, -hl, hl);
      }
      // Remove sigma from upper sigmas
      for (int iSigma = iSigmaLow[Upper] ; (iSigma < iSigmaMax[Upper] && iSigma < iSigmaUp[Upper] ) ; iSigma++) {
        sigma = pGL->xs[iSigma];
        if (sigma > fSigmaMin[Upper]) {
          break;
        }
        iSigmaLow[Upper] = iSigma + 1;
        UpdateSIR<T,A>(sysparm, f0, pl - sigma, -pGL->ws[iSigma], adjacent, z2, &ImpulseResponse);
      }

      // Add sigma to upper sigmas
      for (int iSigma = iSigmaUp[Upper] ; iSigma < iSigmaMax[Upper] ; iSigma++) {
        sigma = pGL->xs[iSigma];
        if (sigma > fSigmaMax[Upper]) {
          break;
        }
        iSigmaUp[Upper] = iSigma + 1;
        UpdateSIR<T,A>(sysparm, f0, pl - sigma, pGL->ws[iSigma], adjacent, z2, &ImpulseResponse);
      }

      // Remove sigma from lower sigmas
      for (int iSigma = iSigmaUp[Lower] ; (iSigma > iSigmaMin[Lower] && iSigma > iSigmaLow[Lower]) ; iSigma--) {
        assert(iSigma > 0 && "Lower sigma out of bounds");
        sigma = pGL->xs[iSigma-1];
        if (sigma < fSigmaMax[Lower]) {
          break;
        }
        iSigmaUp[Lower] = iSigma - 1;
        UpdateSIR<T,A>(sysparm, f0, pl - sigma, -pGL->ws[iSigma-1], adjacent, z2, &ImpulseResponse);
      }

      for (int iSigma = iSigmaLow[Lower] ; (iSigma > iSigmaMin[Lower]) ; iSigma--) {

        assert(iSigma > 0 && "Lower sigma out of bounds");
        sigma = pGL->xs[iSigma-1];
        if ((sigma < fSigmaMin[Lower])) {
          break;
        }
        iSigmaLow[Lower] = iSigma - 1;
        UpdateSIR<T,A>(sysparm, f0, pl - sigma, pGL->ws[iSigma-1], adjacent, z2, &ImpulseResponse);
      }

#if SPS_DEBUG && 0
      upperLimits.push_back(std::make_pair(iSigmaLow[Upper], iSigmaUp[Upper]));
      lowerLimits.push_back(std::make_pair(iSigmaLow[Lower], iSigmaUp[Lower]));
#endif

#if 0
      debug_print("fSigmaMin[Lower]: %1.5f, fSigmaMax[Lower]: %1.5f\n", fSigmaMin[Lower], fSigmaMax[Lower]);
#elif 0
      debug_print("fSigmaMin[Upper]: %1.5f, fSigmaMax[Upper]: %1.5f\n", fSigmaMin[Upper], fSigmaMax[Upper]);
#endif

      // TODO: Make this work for non-sampled delays
      T edge = ImpulseResponse.EvaluateTSD(t,W,f0);

      assert(iSampleSignal < (int)nSamples);

      // update output
      odata[iSampleSignal] += factor * adjacent * edge;
    }

#if SPS_DEBUG && 0
    std::cout << "[";
    for (auto &elem : lowerLimits) {
      std::cout << elem.second - elem.first << " ";
    }
    std::cout << "]" << std::endl;

    std::cout << "[";
    for (auto &elem : upperLimits) {
      std::cout << elem.second - elem.first << " ";
    }
    std::cout << "]" << std::endl;
#endif
  }


  template <class T, template <typename> class A>
  T TransientSingleRect(const sysparm_t<T>* sysparm,
                        const ApertureData<T>* data,
                        const T* pos, const size_t nPositions,
                        T** odata, size_t* nSamples, int mask)
  {
    const T c  = sysparm->c;
    const T fs = sysparm->fs;
    const T W  = sysparm->w;

    const T f0 = data->m_f0;

    // Reset output
    *odata = nullptr;
    *nSamples = 0;

    sps::bbox_t<T> box;
    data->ExtentGet(box);

    debug_print("nPositions: %zu\n",nPositions);

    size_t nElements = 0, nSubElements = 0;

    const T* delays = nullptr;
    const T* apodizations = nullptr;
    const sps::element_rect_t<T>** elementsArray = nullptr;

    data->ApodizationsRefGet(&nElements, apodizations);
    data->DelaysRefGet(&nElements, delays);
    data->ElementsRefGet(&nElements, &nSubElements, elementsArray);

    SPS_UNREFERENCED_PARAMETERS(apodizations, delays);

    // Simple box surrounding the scatters
    sps::bbox_t<T> scatter_box = sps::bbox_t<T>();
    sps::compute_bounding_box3(pos, nPositions, &scatter_box);

    auto uws = sps::deleted_aligned_array<T>();
    auto vws = sps::deleted_aligned_array<T>();
    auto uxs = sps::deleted_aligned_array<T>();
    auto vxs = sps::deleted_aligned_array<T>();

    const size_t iElement = 0;
    const size_t jElement = 0;
    const sps::element_rect_t<T>& element = elementsArray[iElement][jElement];

    /********************************************
     * Compute abcissas and weights             *
     ********************************************/
    CalcWeightsAndAbcissaeScaled(sysparm, element, std::move(uxs), std::move(uws),
                                 std::move(vxs),std::move(vws));

    GLQuad2D<T> uv;
    uv.u.xs     = uxs.get();
    uv.u.ws     = uws.get();
    uv.u.nx     = sysparm->nDivW;
    uv.v.xs     = vxs.get();
    uv.v.ws     = vws.get();
    uv.v.nx     = sysparm->nDivH;

    GLQuad1D<T> u;
    u.xs = uxs.get();
    u.ws = uws.get();
    u.nx = sysparm->nDivW;

    GLQuad1D<T> v;
    v.xs = vxs.get();
    v.ws = vws.get();
    v.nx = sysparm->nDivH;

    sps::bbox_t<T> element_box = sps::bbox_t<T>();
    data->ElementExtentGet(iElement, element_box);

    T distNear;
    T distFar;

    // Uses the 8 corner points
    sps::dists_most_distant_and_closest(scatter_box, element_box, &distNear, &distFar);

    T tStart = (distNear / c);
    T tEnd   = (distFar / c) + W;

    T fSampleSignalStart = fs * tStart;

    int iSampleSignalStart = (int) floor(fSampleSignalStart);

    // Allocate output
    int _nSamples = (int) ceil(T(2.0) + fs*(tEnd - tStart));

    *odata = (T*) malloc(nPositions*_nSamples*sizeof(T));
    memset(*odata, 0, nPositions*_nSamples*sizeof(T));
    *nSamples = _nSamples;

    // Ugly to use sofus::sysparm_t<T>
    sofus::sysparm_t<T> timeDomainParm;
    timeDomainParm.fs = sysparm->fs;
    timeDomainParm.c  = sysparm->c;

    const T scale = T(1.0);

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

      if (mask & 0x01) {
        DirectWaveSingle<T,A>(sysparm,
                              &element, scale,
                              f0,
                              &pld,
                              &uv,
                              delay,
                              iSampleSignalStart,
                              _nSamples,
                              &((*odata)[iPosition*_nSamples]));
      }

      if (mask & 0x02) {
        EdgeResponse<T, A>(sysparm, &element, scale,
                           0, &pld,
                           &v,
                           f0,
                           delay,
                           iSampleSignalStart,
                           _nSamples,
                           &((*odata)[iPosition*_nSamples]));
      }
      if (mask & 0x04) {
        EdgeResponse<T, A>(sysparm, &element, scale,
                           1, &pld,
                           &u,
                           f0,
                           delay,
                           iSampleSignalStart,
                           _nSamples,
                           &((*odata)[iPosition*_nSamples]));
      }
      if (mask & 0x08) {
        EdgeResponse<T, A>(sysparm, &element, scale,
                           2, &pld,
                           &v,
                           f0,
                           delay,
                           iSampleSignalStart,
                           _nSamples,
                           &((*odata)[iPosition*_nSamples]));
      }
      if (mask & 0x10) {
        EdgeResponse<T, A>(sysparm, &element, scale,
                           3, &pld,
                           &u,
                           f0,
                           delay,
                           iSampleSignalStart,
                           _nSamples,
                           &((*odata)[iPosition*_nSamples]));
      }
    }
    return T(0.0);
  }

  template void
  EdgeResponse<float, ToneBurst>(const sysparm_t<float>* sysparm,
                                 const sps::element_rect_t<float>* element,
                                 const float& scale,
                                 const size_t iEdge,
                                 const sofus::proj_limit_dist_t<float>* pld,
                                 const GLQuad1D<float>* pGL,
                                 const float& f0,
                                 const float& delay,
                                 const int& iSampleSignalStart,
                                 const size_t& nSamples,
                                 float* odata);

  template void
  EdgeResponse<float, HanningWeightedPulse>(const sysparm_t<float>* sysparm,
      const sps::element_rect_t<float>* element,
      const float& scale,
      const size_t iEdge,
      const sofus::proj_limit_dist_t<float>* pld,
      const GLQuad1D<float>* pGL,
      const float& f0,
      const float& delay,
      const int& iSampleSignalStart,
      const size_t& nSamples,
      float* odata);

  template void
  DirectWaveSingle<float, ToneBurst>(const sysparm_t<float>* sysparm,
                                     const sps::element_rect_t<float>* element,
                                     const float& scale,
                                     const float& f0,
                                     const sofus::proj_limit_dist_t<float>* pld,
                                     const GLQuad2D<float>* uv,
                                     const float delay,
                                     const int iSampleSignalStart,
                                     const size_t nSamples,
                                     float* odata);

  template void
  DirectWaveSingle<float, HanningWeightedPulse>(const sysparm_t<float>* sysparm,
      const sps::element_rect_t<float>* element,
      const float& scale,
      const float& f0,
      const sofus::proj_limit_dist_t<float>* pld,
      const GLQuad2D<float>* uv,
      const float delay,
      const int iSampleSignalStart,
      const size_t nSamples,
      float* odata);

  template float
  TransientSingleRect<float, ToneBurst>(const sysparm_t<float>* sysparm,
                                        const ApertureData<float>* data,
                                        const float* pos, const size_t nPositions,
                                        float** odata, size_t* nSamples, int mask);

  template float
  TransientSingleRect<float, HanningWeightedPulse>(const sysparm_t<float>* sysparm,
      const ApertureData<float>* data,
      const float* pos, const size_t nPositions,
      float** odata, size_t* nSamples, int mask);

#if FNM_DOUBLE_SUPPORT

  template void
  DirectWaveSingle(const sysparm_t<double>* sysparm,
                   const sps::element_rect_t<double>* element,
                   const double& scale,
                   const double& f0,
                   const sofus::proj_limit_dist_t<double>* pld,
                   const GLQuad2D<double>* uv,
                   const double delay,
                   const int iSampleSignalStart,
                   const size_t nSamples,
                   double* odata);

  template void
  EdgeResponse<double, ToneBurst>(const sysparm_t<double>* sysparm,
                                  const sps::element_rect_t<double>* element,
                                  const double& scale,
                                  const size_t iEdge,
                                  const sofus::proj_limit_dist_t<double>* pld,
                                  const GLQuad1D<double>* pGL,
                                  const double& f0,
                                  const double& delay,
                                  const int& iSampleSignalStart,
                                  const size_t& nSamples,
                                  double* odata);

  template double
  TransientSingleRect<double, ToneBurst>(const sysparm_t<double>* sysparm,
                                         const ApertureData<double>* data,
                                         const double* pos, const size_t nPositions,
                                         double** odata, size_t* nSamples, int mask);
#endif

}
