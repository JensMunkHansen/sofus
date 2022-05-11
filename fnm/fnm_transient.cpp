/**
 * @file   fnm_transient.cpp
 * @author  <jens.munk.hansen@gmail.com>
 * @date   Tue Nov  7 14:57:46 2017
 *
 * @brief
 *
 * Copyright 2017 Jens Munk Hansen
 *
 */

/*
 *
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
 *
 */

#include <sps/debug.h>
#include <fnm/config.h>

#include <utility>
#include <fnm/fnm_transient.hpp>
#include <fnm/fnm_data.hpp>
#include <fnm/fnm_common.hpp>
#include <fnm/fnm_basis.hpp>

#include <sofus/rect_int_limits.hpp>  // calcProjectionAndLimits
#include <sofus/sofus_calc.hpp>       // ComputeBoxTimes
#include <sps/msignals.hpp>

#include <sps/algorithm>

#include <sps/sse2-linear-search.hpp>

#ifndef SHOW_LIMITS
# define SHOW_LIMITS 0
#endif

#ifndef DST
# define DST 0
#endif

namespace fnm {

#if 1 //def FNM_PULSED_WAVE

// TODO(JMH): Support odd number of abcissas.
template <class T>
T FourDirect(const GLQuad2D<T>* pUV,
             const sps::element_rect_t<T>* pElement,
             const sofus::proj_limit_dist_t<T>* pLimit) {
  T result = T(0.0);

  T uLow  = fabs(pLimit->u) - pElement->hw;
  T uHigh = fabs(pLimit->u) + pElement->hw;
  T vLow  = fabs(pLimit->v) - pElement->hh;
  T vHigh = fabs(pLimit->v) + pElement->hh;

  // Adjacent
  T v[2];
  v[0] = pElement->hh - fabs(pLimit->v);
  v[1] = pElement->hh + fabs(pLimit->v);

  T u[2];
  u[0] = pElement->hw - fabs(pLimit->u);
  u[1] = pElement->hw + fabs(pLimit->u);

  // Offset of abcissa
  T sp = T(0.5) * (uHigh + uLow);

  T l2[2];
  l2[0] = SQUARE(v[0]);
  l2[1] = SQUARE(v[1]);

  T sInt[2] = {T(0.0), T(0.0)};

  for (size_t iU = 0 ; iU < pUV->u.nx ; iU++) {
    T s2 = SQUARE(pUV->u.xs[iU] + sp);
    // Issue is here for odd number of abcissas
    if (fabs(s2+l2[0]) > T(0.0))
      sInt[0] += pUV->u.ws[iU] * (T(1.0) / (s2+l2[0]));
    if (fabs(s2+l2[1]) > T(0.0))
      sInt[1] += pUV->u.ws[iU] * (T(1.0) / (s2+l2[1]));
  }
  // Multiplied in the loop to prevent near infinities
  result += sInt[0]*v[0];
  result += sInt[1]*v[1];

  // Offset of abcissa
  sp = T(0.5) * (vHigh + vLow);

  l2[0] = SQUARE(u[0]);
  l2[1] = SQUARE(u[1]);

  sInt[0] = T(0.0);
  sInt[1] = T(0.0);

  for (size_t iV = 0 ; iV < pUV->v.nx ; iV++) {
    T s2 = SQUARE(pUV->v.xs[iV] + sp);
    if (fabs(s2+l2[0]) > T(0.0))
      sInt[0] += pUV->v.ws[iV] * (T(1.0) / (s2+l2[0]));
    if (fabs(s2+l2[1]) > T(0.0))
      sInt[1] += pUV->v.ws[iV] * (T(1.0) / (s2+l2[1]));
  }
  // Stripes if a limit on size of u[0] is introduced
  result += sInt[0]*u[0];
  result += sInt[1]*u[1];

  return result;
}


template <class T,  template <typename> class A>
void DirectWaveSingle(
  const sysparm_t<T>* pSysparm,
  const sps::element_rect_t<T>* pElement,
  const GLQuad2D<T>* pUV,
  const T& scale,
  const sofus::proj_limit_dist_t<T>* pld, const T delay,
  const int iSampleSignalStart, const size_t nSamples, T* odata) {
  // Constants
  const T W = pSysparm->w;
  const T c = pSysparm->c;
  const T fs = pSysparm->fs;
  const T f0 = pSysparm->f0;

  const T invfs = T(1.0) / fs;

  // M_1_2PI
  const T factor = - scale * pSysparm->rho * c / T(M_2PI);

  debug_cond_print(0, "factor: %f\n", factor);

  T z = pld->dist2plane;

  T t1;

  t1 = z / c;  // tau_2

  const float fSampleStart = pld->fSampleStart;
  const float fSampleStop  = pld->fSampleStop;

  // Delay is already included in fSampleStart - this is the right solution.
  int iSampleStart = static_cast<int>(fSampleStart) + 1;
  int iSampleStop  = static_cast<int>(fSampleStop + W*fs) + 1;

  T spatial = FourDirect(pUV, pElement, pld);

  debug_cond_print(0, "spatial: %f\n", spatial);

  assert(iSampleStart > (iSampleSignalStart - 1) && "Access before start");
  assert((iSampleStop - 1 - iSampleSignalStart) <
         static_cast<int>(nSamples) && "Access after stop");

  if ((iSampleStart < iSampleSignalStart) ||
      (iSampleStop > static_cast<int>(nSamples) + iSampleSignalStart)) {
    return;
  }

  debug_cond_print(0, "iSampleStart: %d, iSampleStop: %d\n", iSampleStart, iSampleStop);

  // Vectorize four sample at a time + (nSamples % 4)
  for (int iSample = iSampleStart ; iSample < iSampleStop ; iSample++) {
    int iSampleSignal = iSample - iSampleSignalStart;
    T t = invfs * iSample - delay;

    T direct = - A<T>::Evaluate(t - t1, W, f0);

    T value = direct * spatial;

#if SOFUS_ENABLE_ATTENUATION
    // TODO(JMH): Support frequency dependent attenuation
    if (pSysparm->use_att) {
      value = value * exp(- t * pSysparm->c * pSysparm->att);
    }
#endif
    odata[iSampleSignal] += factor * value;
  }
}

template <typename T, template <typename> class A>
void UpdateSIR(const sysparm_t<T>* pSysparm,
               const T ssigma, const T sweight,
               const T adjacent, const T z2,
               A<T>* pulse) {
  const T f0 = pSysparm->f0;
  T denom = SQUARE(adjacent) + SQUARE(ssigma);
  T tau = sqrt(std::max<T>(T(0.0), z2 + denom)) / pSysparm->c;
  T scale = denom > T(0.0) ? sweight / denom : T(0.0);

  // Could skip the call if denom == 0.0
  pulse->UpdateSpatial(scale, tau, pSysparm->w, f0);
}

// TODO(JMH):
//  - Support odd number of abcissas.
//  - Issue when nDivH is odd, iEdge=0,2 and limit->v = 0.0
template <typename T, template <typename> class A>
void EdgeResponse(
  const sysparm_t<T>* pSysparm,
  const sps::element_rect_t<T>* pElement,
  const T& scale,
  const size_t iEdge,
  const sofus::proj_limit_dist_t<T>* pld,
  const GLQuad1D<T>* pGL,
  const T& delay,
  const int& iSampleSignalStart, const size_t& nSamples, T* odata) {
  const T fs = pSysparm->fs;
  const T f0 = pSysparm->f0;
  const T invfs = T(1.0) / fs;
  const T W = pSysparm->w;
  const T c = pSysparm->c;

  const T factor = - scale * pSysparm->rho * c / T(M_2PI);

  const T z2 = SQUARE(pld->dist2plane);

  T hl = T(0.0);  // Half long
  T pl = T(0.0);  // Point long

  T hs = T(0.0);  // Half short
  T ps = T(0.0);  // Point short

  T adjacent = T(0.0);  // Orthogonal to integration range

  size_t iVertices[2] = {0, 0};  // Vertices for active edge

  size_t nl = 0;  // Number of abcissas for integration along l

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
  }

  switch (iEdge) {
  case 1:
  case 3:
    hl = pElement->hw;
    hs = pElement->hh;
    pl = pld->u;
    ps = pld->v;
    nl = pSysparm->nDivW;
    break;
  case 0:
  case 2:
    hl = pElement->hh;
    hs = pElement->hw;
    pl = pld->v;
    ps = pld->u;
    nl = pSysparm->nDivH;
    break;
  }

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
  }

  debug_print("ps: %f, pl: %f, hs: %f, hl: %f, adjacent: %f, nl: %zu\n",
              ps, pl, hs, hl, adjacent, nl);

  T t1 = pld->vdists[iVertices[0]] / c;
  T t2 = pld->vdists[iVertices[1]] / c;

  T t12min, t12max;

  enum Segment {
    Lower = 0,
    Upper = 1,
  };

  // Limits - used for disabling upper or lower
  int iSigmaMin[2] = {0, 0};
  int iSigmaMax[2] = {0, 0};

  // Nothing is disabled by setting this
  iSigmaMin[Lower] = 0;
  iSigmaMax[Lower] = static_cast<int>(nl);
  iSigmaMin[Upper] = 0;
  iSigmaMax[Upper] = static_cast<int>(nl);

  // Dynamic limits
  int iSigmaLow[2] = {0, 0};
  int iSigmaUp[2] = {0, 0};

  if (SPS_UNLIKELY(fabs(pl) < hl)) {
    // Time to reach edge
    t12min = sqrt(std::max<T>(T(0.0), z2 + SQUARE(adjacent))) / c;

    auto end = pGL->xs + nl;

    // Greater or equal
    auto it = std::find_if(pGL->xs, end, [&](T a)->bool {return !(pl > a);});
    if (it != end) {
      int iProj = static_cast<int>(std::distance(pGL->xs, it));
      iSigmaLow[Lower] = iProj;
      iSigmaUp[Lower]  = iProj;
      iSigmaLow[Upper] = iProj;
      iSigmaUp[Upper]  = iProj;
    } else {
      iSigmaLow[Lower] = static_cast<int>(nl);
      iSigmaUp[Lower] = static_cast<int>(nl);
      // Upper integral is disabled usig iSigmaMax
      iSigmaLow[Upper] = static_cast<int>(nl);
      iSigmaUp[Upper]  = static_cast<int>(nl);
    }
  } else {
    /*
      TODO(JMH): Consider avoiding fabs(pl) < hl and
      work instead with [-oo,-hl[, [-hl;hl[ ,[hl;oo]
    */

    // Time to nearest vertex
    t12min = std::min<T>(t1, t2);
    if (pl >= hl) {  // pl == hl => Lower integral is active.
      iSigmaLow[Lower] = static_cast<int>(nl);
      iSigmaUp[Lower]  = static_cast<int>(nl);
      // Upper integral is disabled
      iSigmaMax[Upper] = 0;
    } else {
      iSigmaLow[Upper] = 0;
      iSigmaUp[Upper]  = 0;
      // Lower integral is disabled
      iSigmaMin[Lower] = static_cast<int>(nl);
    }
  }

  // Time to furthest vertex
  t12max = std::max<T>(t1, t2);

  debug_print("t12min: %f\n", t12min);

  T fSampleStart = fs * t12min;
  T fSampleStop = fs * t12max;

  debug_print("fSampleStart: %f, fSampleStop: %f\n", fSampleStart, fSampleStop);

  int iSampleStart = static_cast<int>(fSampleStart + delay*fs) + 1;
  int iSampleStop  = static_cast<int>(fSampleStop + delay*fs + W*fs) + 1;

  if ((iSampleStart < iSampleSignalStart) ||
      (iSampleStop > static_cast<int>(nSamples) + iSampleSignalStart)) {
    return;
  }

  debug_print("iSampleStart: %d, iSampleStop: %d\n", iSampleStart, iSampleStop);

  // Lower values for ranges: [-hl ; pl] and [pl, hl]
  T fSigmaMin[2] = {0.0f, 0.0f};

  // Upper values for ranges: [-hl ; pl] and [pl, hl]
  T fSigmaMax[2] = {0.0f, 0.0f};

  A<T> ImpulseResponse = A<T>();

  ImpulseResponse.ResetSpatial();

#if SPS_DEBUG && SHOW_LIMITS
  typedef std::vector<std::pair<int, int> > limitVector;
  limitVector upperLimits = limitVector();
  limitVector lowerLimits = limitVector();

  // iEdge, iSegment, iSample, breakmask, iSigma
  std::vector<std::array<int, 5> > stateVector =
    std::vector<std::array<int, 5> >();

  std::vector<std::array<float, 4> > flimits =
    std::vector<std::array<float, 4> >();
#endif

#ifdef SPS_DEBUG
  if (iSampleStart < iSampleSignalStart) {
    debug_print("\033[31;1mSamples too early: u: %f, v: %f\033[0m\n",
                pld->u, pld->v);
  }
  if (iSampleStop >= static_cast<int>(nSamples) + iSampleSignalStart) {
    debug_print("\033[31;1mSamples too late: u: %f, v: %f\033[0m\n",
                pld->u, pld->v);
  }
#endif

  // Issue
  assert(iSampleStart > (iSampleSignalStart - 1) && "Access before start");
  assert((iSampleStop - 1 - iSampleSignalStart) <
         static_cast<int>(nSamples) && "Access after stop");

  if ((iSampleStart < iSampleSignalStart) ||
      (iSampleStop > iSampleSignalStart + static_cast<int>(nSamples))) {
    return;
  }

  debug_cond_print(DST && 1, "SigmaUpUpper: %d\n", iSigmaUp[Upper]);
  debug_cond_print(DST && 1, "adjacent: %f\n", adjacent);
  for (int iSample = iSampleStart ; iSample < iSampleStop ; iSample++) {
    int iSampleSignal = iSample - iSampleSignalStart;
    T t = invfs * iSample - delay;

    // Compute two sigma segments (asymmetry could be if we start to early)
    T sqrtArg = SQUARE(c*t) - z2 - SQUARE(adjacent);
    sqrtArg = std::max<T>(T(0.0), sqrtArg);
    T deltaSigma = sqrt(sqrtArg);

    // Lower sigma range
    fSigmaMax[Lower] = sps::clamp(pl /*  Subtract shortening here */, -hl, hl);
    fSigmaMin[Lower] = sps::clamp(pl - deltaSigma, -hl, hl);

    // Upper sigma range
    fSigmaMax[Upper] = sps::clamp(pl + deltaSigma, -hl, hl);
    fSigmaMin[Upper] = sps::clamp(pl /* Add shortening here */, -hl, hl);

    if (SPS_UNLIKELY((t > t12min + W) && (t < t12max + W))) {
      // Shorten ranges
      deltaSigma =
        sqrt(std::max<T>(T(0.0), SQUARE(c*(t-W)) - z2 - SQUARE(adjacent)));
      fSigmaMin[Upper] = sps::clamp(pl + deltaSigma, -hl, hl);
      fSigmaMax[Lower] = sps::clamp(pl - deltaSigma, -hl, hl);
    }

#if SPS_DEBUG && SHOW_LIMITS
    std::array<float, 4> state = {{
        fSigmaMax[Upper],
        fSigmaMin[Upper],
        fSigmaMax[Lower],
        fSigmaMin[Lower]
      }
    };
    flimits.push_back(state);
#endif

    debug_cond_print(0, "fSigmaMinUpper: %f, fSigmaMaxUpper: %f, \
                     fSigmaMinLower: %f, fSigmaMaxLower: %f\n", fSigmaMin[Upper], fSigmaMax[Upper],
                     fSigmaMin[Lower], fSigmaMax[Lower]);

    /*****************************************
     * Remove sigma from upper sigmas
     *****************************************/
    for (int iSigma = iSigmaLow[Upper] ;
         (iSigma < iSigmaMax[Upper] && iSigma < iSigmaUp[Upper]);
         iSigma++) {
      T sigma = pGL->xs[iSigma];
      debug_cond_print(DST &&iSample == iSampleStart + 30,
                       "RU, iSigma: %d, fSigmaMinUpper: %f\n",
                       iSigma, fSigmaMin[Upper]);

      if (SPS_UNLIKELY(sigma > fSigmaMin[Upper])) {
        break;
      }
      iSigmaLow[Upper] = iSigma + 1;
      debug_cond_print(DST &&iSample == iSampleStart + 30,
                       "RU, iSigma: %d, sigma: %f\n",
                       iSigma, sigma);

#if SPS_DEBUG && SHOW_LIMITS
      std::array<int, 5> state;
      state[0] = (int) iEdge;
      state[1] = 0;
      state[2] = (int) iSample;
      state[3] = 1;
      state[4] = (int) iSigma;
      stateVector.push_back(state);
#endif

      UpdateSIR<T, A>(pSysparm, pl - sigma,
                      -pGL->ws[iSigma], adjacent, z2, &ImpulseResponse);
    }

    /*****************************************
     * Add sigma from upper sigmas
     *****************************************/
    for (int iSigma = iSigmaUp[Upper] ; iSigma < iSigmaMax[Upper] ; iSigma++) {
      T sigma = pGL->xs[iSigma];
      debug_cond_print(DST &&iSample == iSampleStart + 30,
                       "AU, iSigma: %d, fSigmaMaxUpper: %f\n",
                       iSigma, fSigmaMax[Upper]);
      if (SPS_UNLIKELY(sigma > fSigmaMax[Upper])) {
        break;
      }
      debug_cond_print(DST &&iSample == iSampleStart + 30,
                       "AU, iSigma: %d, sigma: %f\n",
                       iSigma, sigma);

      iSigmaUp[Upper] = iSigma + 1;
      debug_cond_print(DST &&iSample == iSampleStart + 30,
                       "ws: %f\n",
                       pGL->ws[iSigma]);

      UpdateSIR<T, A>(pSysparm, pl - sigma,
                      pGL->ws[iSigma], adjacent, z2, &ImpulseResponse);

#if SPS_DEBUG && SHOW_LIMITS
      std::array<int, 5> state;
      state[0] = (int) iEdge;
      state[1] = 1;
      state[2] = (int) iSample;
      state[3] = 1;
      state[4] = (int) iSigma;
      stateVector.push_back(state);
#endif

      debug_cond_print(DST &&iSample == iSampleStart + 30,
                       "AU: terms: %f %f\n",
                       ImpulseResponse.m_fTerms[0],
                       ImpulseResponse.m_fTerms[1]);
    }

#if 1
    /*****************************************
     * Remove sigma from lower sigmas
     *****************************************/
    for (int iSigma = iSigmaUp[Lower] ;
         (iSigma > iSigmaMin[Lower] && iSigma > iSigmaLow[Lower]);
         iSigma--) {
      assert(iSigma > 0 && "Lower sigma out of bounds");
      T sigma = pGL->xs[iSigma-1];
      if (SPS_UNLIKELY(sigma < fSigmaMax[Lower])) {
        break;
      }
      iSigmaUp[Lower] = iSigma - 1;
      UpdateSIR<T, A>(pSysparm, pl - sigma,
                      -pGL->ws[iSigma-1], adjacent, z2, &ImpulseResponse);

#if SPS_DEBUG && SHOW_LIMITS
      std::array<int, 5> state;
      state[0] = (int) iEdge;
      state[1] = 2;
      state[2] = (int) iSample;
      state[3] = 1;
      state[4] = (int) iSigma - 1;
      stateVector.push_back(state);
#endif

      debug_cond_print(DST &&iSample == iSampleStart + 30,
                       "RL: terms: %f %f\n",
                       ImpulseResponse.m_fTerms[0],
                       ImpulseResponse.m_fTerms[1]);
    }

    /*****************************************
     * Add sigma from lower sigmas
     *****************************************/
    for (int iSigma = iSigmaLow[Lower] ;
         (iSigma > iSigmaMin[Lower]);
         iSigma--) {
      assert(iSigma > 0 && "Lower sigma out of bounds");
      T sigma = pGL->xs[iSigma-1];
      if (SPS_UNLIKELY(sigma < fSigmaMin[Lower])) {
        break;
      }
      iSigmaLow[Lower] = iSigma - 1;
      UpdateSIR<T, A>(pSysparm, pl - sigma,
                      pGL->ws[iSigma-1], adjacent, z2, &ImpulseResponse);
#if SPS_DEBUG && SHOW_LIMITS
      std::array<int, 5> state;
      state[0] = (int) iEdge;
      state[1] = 3;
      state[2] = (int) iSample;
      state[3] = 1;
      state[4] = (int) iSigma - 1;
      stateVector.push_back(state);
#endif
      debug_cond_print(DST && iSample == iSampleStart + 30,
                       "AL: terms: %f %f\n",
                       ImpulseResponse.m_fTerms[0],
                       ImpulseResponse.m_fTerms[1]);
    }

#endif

#if SPS_DEBUG && SHOW_LIMITS
    upperLimits.push_back(std::make_pair(iSigmaLow[Upper], iSigmaUp[Upper]));
    lowerLimits.push_back(std::make_pair(iSigmaLow[Lower], iSigmaUp[Lower]));
#endif

    T edge = ImpulseResponse.EvaluateTSD(t, W, f0);

    assert(iSampleSignal < static_cast<int>(nSamples));

#if SOFUS_ENABLE_ATTENUATION
    if (pSysparm->use_att) {
      edge = edge * exp(- t * pSysparm->c * pSysparm->att);
    }
#endif

    edge = edge * factor * adjacent;

    debug_cond_print(0, "edge_term: %f\n", edge);
    // update output
    odata[iSampleSignal] += edge;
  }

#if SPS_DEBUG && SHOW_LIMITS
#if 0
  std::cout << "[";
  for (auto &elem : lowerLimits) {
    std::cout << elem.second - elem.first
              << "[" << elem.first << ";" << elem.second << "]" << " ";
  }
  std::cout << "]" << std::endl;

  std::cout << "[";
  for (auto &elem : upperLimits) {
    std::cout << elem.second - elem.first
              << "[" << elem.first << ";" << elem.second << "]" << " ";
  }
  std::cout << "]" << std::endl;
#endif
  std::cout << "Scalar" << std::endl;
  for (auto &elem : stateVector) {
    for (auto s : elem) {
      std::cout << s << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "fLimits" << std::endl;
  for (auto &elem : flimits) {
    for (auto s : elem) {
      std::cout << s << " ";
    }
    std::cout << std::endl;
  }
#endif
}


template <class T, template <typename> class A>
T TransientSingleRect(const sysparm_t<T>* pSysparm,
                      const ApertureData<T>* data,
                      const T* pos, const size_t nPositions,
                      T** odata, size_t* nSamples, int mask) {
  const T c  = pSysparm->c;
  const T fs = pSysparm->fs;
  const T W  = pSysparm->w;

  // Reset output
  *odata = nullptr;
  *nSamples = 0;

  sps::bbox_t<T> box;
  data->ExtentGet(&box);

  debug_print("nPositions: %zu\n", nPositions);

  size_t nElements = 0, nSubElements = 0;

  const T* delays = nullptr;
  const T* apodizations = nullptr;
  const sps::element_rect_t<T>** ppElements = nullptr;

  data->ApodizationsRefGet(&nElements, apodizations);
  data->DelaysRefGet(&nElements, delays);
  data->ElementsRefGet(&nElements, &nSubElements, ppElements);

  SPS_UNREFERENCED_PARAMETERS(apodizations, delays);

  // Simple box surrounding the scatters
  sps::bbox_t<T> scatter_box = sps::bbox_t<T>();
  sps::compute_bounding_box3(pos, nPositions, &scatter_box);

  auto uws = sps::unique_aligned_array<T>();
  auto vws = sps::unique_aligned_array<T>();
  auto uxs = sps::unique_aligned_array<T>();
  auto vxs = sps::unique_aligned_array<T>();

  const size_t iElement = 0;
  const size_t jElement = 0;
  const sps::element_rect_t<T>& element = ppElements[iElement][jElement];

  /********************************************
   * Compute abcissas and weights             *
   ********************************************/
  CalcWeightsAndAbcissaeScaled(pSysparm, element,
                               std::move(uxs), std::move(uws),
                               std::move(vxs), std::move(vws));

  GLQuad2D<T> uv;
  uv.u.xs     = uxs.get();
  uv.u.ws     = uws.get();
  uv.u.nx     = pSysparm->nDivW;
  uv.v.xs     = vxs.get();
  uv.v.ws     = vws.get();
  uv.v.nx     = pSysparm->nDivH;

  GLQuad1D<T> u;
  u.xs = uxs.get();
  u.ws = uws.get();
  u.nx = pSysparm->nDivW;

  GLQuad1D<T> v;
  v.xs = vxs.get();
  v.ws = vws.get();
  v.nx = pSysparm->nDivH;

  sps::bbox_t<T> element_box = sps::bbox_t<T>();
  data->ElementExtentGet(iElement, &element_box);

  T distNear = T(0.0);
  T distFar  = T(0.0);

  // Uses the 8 corner points
  sps::dists_most_distant_and_closest(scatter_box, element_box,
                                      &distNear, &distFar);

  T tStart = (distNear / c);
  T tEnd   = (distFar / c) + W;

  T fSampleSignalStart = fs * tStart;

  int iSampleSignalStart = static_cast<int>(floor(fSampleSignalStart));

  // Allocate output
  int _nSamples = static_cast<int>(ceil(T(2.0) + fs*(tEnd - tStart)));

  _nSamples = _nSamples + static_cast<int>(W*fs);

  *odata = static_cast<T*>(SPS_MALLOC(nPositions*_nSamples*sizeof(T)));
  memset(*odata, 0, nPositions*_nSamples*sizeof(T));
  *nSamples = _nSamples;

  // Ugly to use sofus::sysparm_t<T>
  sofus::sysparm_t<T> timeDomainParm;
  timeDomainParm.fs = pSysparm->fs;
  timeDomainParm.c  = pSysparm->c;

  const T scale = T(1.0);

  for (size_t iPosition = 0 ; iPosition < nPositions ; iPosition++) {
    sofus::proj_limit_dist_t<T> pld;

    sps::point_t<T> point;
    point[0] = pos[iPosition * 3];
    point[1] = pos[iPosition * 3 + 1];
    point[2] = pos[iPosition * 3 + 2];

    T delay = T(0.0);
    sofus::calcProjectionAndLimits(timeDomainParm, ppElements[0][0],
                                   point, delay,
                                   &pld);
    debug_print("pld.u: %f, pld.v: %f, pld.dist2plane: %f, "
                "fSampleStart: %f, fSampleStop: %f\n",
                pld.u, pld.v, pld.dist2plane,
                pld.fSampleStart, pld.fSampleStop);

    // SGV: Double
    if (mask & 0x01) {
      DirectWaveSingle<T, A>(
        pSysparm, &element, &uv, scale,
        &pld, delay,
        iSampleSignalStart, _nSamples, &((*odata)[iPosition*_nSamples]));
    }
    if (mask & 0x02) {
      EdgeResponse<T, A>(
        pSysparm, &element, scale,
        0, &pld,
        &v,
        delay,
        iSampleSignalStart, _nSamples, &((*odata)[iPosition*_nSamples]));
    }
    if (mask & 0x04) {
      EdgeResponse<T, A>(
        pSysparm, &element, scale,
        1, &pld,
        &u,
        delay,
        iSampleSignalStart, _nSamples, &((*odata)[iPosition*_nSamples]));
    }
    if (mask & 0x08) {
      EdgeResponse<T, A>(
        pSysparm, &element, scale,
        2, &pld,
        &v,
        delay,
        iSampleSignalStart, _nSamples, &((*odata)[iPosition*_nSamples]));
    }
    if (mask & 0x10) {
      EdgeResponse<T, A>(
        pSysparm, &element, scale,
        3, &pld,
        &u,
        delay,
        iSampleSignalStart, _nSamples, &((*odata)[iPosition*_nSamples]));
    }
  }
  return tStart;
}

#endif

template void
EdgeResponse<float, ToneBurst>(
  const sysparm_t<float>* pSysparm,
  const sps::element_rect_t<float>* pElement,
  const float& scale,
  const size_t iEdge,
  const sofus::proj_limit_dist_t<float>* pld,
  const GLQuad1D<float>* pGL,
  const float& delay,
  const int& iSampleSignalStart, const size_t& nSamples, float* odata);

template void
EdgeResponse<float, HanningWeightedPulse>(
  const sysparm_t<float>* pSysparm,
  const sps::element_rect_t<float>* pElement,
  const float& scale,
  const size_t iEdge,
  const sofus::proj_limit_dist_t<float>* pld,
  const GLQuad1D<float>* pGL,
  const float& delay,
  const int& iSampleSignalStart, const size_t& nSamples, float* odata);

template void
EdgeResponse<float, Identity>(
  const sysparm_t<float>* pSysparm,
  const sps::element_rect_t<float>* pElement,
  const float& scale,
  const size_t iEdge,
  const sofus::proj_limit_dist_t<float>* pld,
  const GLQuad1D<float>* pGL,
  const float& delay,
  const int& iSampleSignalStart, const size_t& nSamples, float* odata);

template void
DirectWaveSingle<float, ToneBurst>(
  const sysparm_t<float>* pSysparm,
  const sps::element_rect_t<float>* pElement,
  const GLQuad2D<float>* uv,
  const float& scale,
  const sofus::proj_limit_dist_t<float>* pld,
  const float delay,
  const int iSampleSignalStart, const size_t nSamples, float* odata);

template void
DirectWaveSingle<float, HanningWeightedPulse>(
  const sysparm_t<float>* pSysparm,
  const sps::element_rect_t<float>* pElement,
  const GLQuad2D<float>* uv,
  const float& scale,
  const sofus::proj_limit_dist_t<float>* pld,
  const float delay,
  const int iSampleSignalStart, const size_t nSamples, float* odata);

template void
DirectWaveSingle<float, Identity>(
  const sysparm_t<float>* pSysparm,
  const sps::element_rect_t<float>* pElement,
  const GLQuad2D<float>* uv,
  const float& scale,
  const sofus::proj_limit_dist_t<float>* pld,
  const float delay,
  const int iSampleSignalStart, const size_t nSamples, float* odata);

template float
TransientSingleRect<float, ToneBurst>(
  const sysparm_t<float>* pSysparm,
  const ApertureData<float>* data,
  const float* pos, const size_t nPositions,
  float** odata, size_t* nSamples, int mask);

template float
TransientSingleRect<float, HanningWeightedPulse>(
  const sysparm_t<float>* pSysparm,
  const ApertureData<float>* data,
  const float* pos, const size_t nPositions,
  float** odata, size_t* nSamples, int mask);

template float
TransientSingleRect<float, Identity>(
  const sysparm_t<float>* pSysparm,
  const ApertureData<float>* data,
  const float* pos, const size_t nPositions,
  float** odata, size_t* nSamples, int mask);

#if FNM_DOUBLE_SUPPORT

template void
DirectWaveSingle<double, ToneBurst>(
  const sysparm_t<double>* pSysparm,
  const sps::element_rect_t<double>* pElement,
  const GLQuad2D<double>* uv,
  const double& scale,
  const sofus::proj_limit_dist_t<double>* pld,
  const double delay,
  const int iSampleSignalStart,
  const size_t nSamples,
  double* odata);

template void
DirectWaveSingle<double, HanningWeightedPulse>(
  const sysparm_t<double>* pSysparm,
  const sps::element_rect_t<double>* pElement,
  const GLQuad2D<double>* uv,
  const double& scale,
  const sofus::proj_limit_dist_t<double>* pld,
  const double delay,
  const int iSampleSignalStart, const size_t nSamples, double* odata);

template void
EdgeResponse<double, ToneBurst>(
  const sysparm_t<double>* pSysparm,
  const sps::element_rect_t<double>* pElement,
  const double& scale,
  const size_t iEdge,
  const sofus::proj_limit_dist_t<double>* pld,
  const GLQuad1D<double>* pGL,
  const double& delay,
  const int& iSampleSignalStart, const size_t& nSamples, double* odata);

template void
EdgeResponse<double, HanningWeightedPulse>(const sysparm_t<double>* pSysparm,
    const sps::element_rect_t<double>* pElement,
    const double& scale,
    const size_t iEdge,
    const sofus::proj_limit_dist_t<double>* pld,
    const GLQuad1D<double>* pGL,
    const double& delay,
    const int& iSampleSignalStart, const size_t& nSamples, double* odata);

template double
TransientSingleRect<double, ToneBurst>(
  const sysparm_t<double>* pSysparm,
  const ApertureData<double>* data,
  const double* pos, const size_t nPositions,
  double** odata, size_t* nSamples, int mask);

template double
TransientSingleRect<double, HanningWeightedPulse>(
  const sysparm_t<double>* pSysparm,
  const ApertureData<double>* data,
  const double* pos, const size_t nPositions,
  double** odata, size_t* nSamples, int mask);
#endif

}  // namespace fnm

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
