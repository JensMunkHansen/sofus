/**
 * @file   circular_calc.cpp
 * @author Jens Munk Hansen <jmh@debian9laptop.parknet.dk>
 * @date   Thu Apr  5 21:55:02 2018
 *
 * @brief Calculation routines for harmonic and transient for a circular disc
 *
 * Copyright 2017 Jens Munk Hansen
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

#include <sps/debug.h>
#include <sps/profiler.h>
#include <sps/algorithm>

#include <vector>
#include <sps/queue.hpp>
#include <sps/sps_threads.hpp>

#include <fnm/fnm.hpp>
#include <fnm/fnm_profiling.hpp>
#include <fnm/circular_data.hpp>
#include <fnm/circular_calc.hpp>
#include <fnm/fnm_basis.hpp>  // Basis functions

#include <gl/gl.hpp>

#if defined(__GNUG__) || (defined(_MSC_VER) && (_MSC_VER >= 1900))
# include <sps/threadpool.hpp>
#endif

namespace fnm {
template <class T>
struct PwFnmThreadArg {
  /// Pointer to system parameters
  const fnm::sysparm_t<T>* pSysparm;
  /// Pointer to aperture data
  const CircularApertureData<T>* pApertureData;
  /// Gauss-Legendre coefficienter
  const GLQuad1D<T>* pGL;
  /// Position data
  const T* pPos;

  /// Start index of output
  size_t iOutputBegin;
  /// Stop index of output
  size_t iOutputEnd;
  /// Number of positions (all threads)
  size_t nPositions;
  /// Output
  T* pField;
  /// Size of output data (without temporal - there isn't any)
  size_t nSamples;
  /// Integration start
  int iSampleStart;
  /// Pointer to queue
  sps::queue<float>* pQueue;
  /// Thread ID
  size_t threadId;
  /// CPU ID
  int cpuId;
};

template <class T>
void CalcOneDimensionalWeightsAndAbcissae(
  const size_t N,
  const T begin,
  const T end,
  sps::unique_aligned_array<T> &&uxs,
  sps::unique_aligned_array<T> &&uweights) {
  // Abcissa and weights
  uxs      = sps::unique_aligned_array_create<T>(N);
  uweights = sps::unique_aligned_array_create<T>(N);

  const T sm = T(0.5) * (end - begin);
  const T sp = T(0.5) * (end + begin);

  // Common weights and abcissa
  for (size_t i = 0 ; i < N ; i++) {
    gl::GLNode qp = gl::GL(N, i);
    // Conversion from double to T
    uxs[i]        = sm * T(qp.value) + sp;
    uweights[i]   = T(qp.weight);
  }
}

template <class T,  template <typename> class A>
#if defined(HAVE_PTHREAD_H)
void* CalcCircularPwThreadFunc(void* ptarg)
#else
unsigned int __stdcall CalcCircularPwThreadFunc(void *ptarg)
#endif
{
  fnm::PwFnmThreadArg<T>* pThreadArg =
    reinterpret_cast<fnm::PwFnmThreadArg<T>*>(ptarg);

  const fnm::sysparm_t<T>* pSysparm = pThreadArg->pSysparm;

  const CircularApertureData<T>* pApertureData = pThreadArg->pApertureData;
  const GLQuad1D<T>* pGL               = pThreadArg->pGL;
  const T* pos                         = pThreadArg->pPos;
  T* odata                             = pThreadArg->pField;
  size_t nSamples                      = pThreadArg->nSamples;
  int iSampleSignalStart               = pThreadArg->iSampleStart;

#ifdef HAVE_THREAD
  setcpuid(pThreadArg->cpuId);
#endif
  const size_t nPositions          = pThreadArg->nPositions;
  sps::queue<float>* pQueue        = pThreadArg->pQueue;

  // Initial time
  double tProfStart = 0.0;

  if (pThreadArg->threadId == 0) {
    tProfStart = sps::profiler::time();
  }

  size_t nCallbackPeriod = 1;
  double duration = 0.0;

  const T eps = std::numeric_limits<T>::epsilon();

  const sps::element_circular_t<T>& element = pApertureData->m_element;
  const T a = element.circle.radius;
  const T a2 = SQUARE(a);

  const T fs = pSysparm->fs;
  const T invfs = T(1.0) / fs;
  const T c = pSysparm->c;
  const T W  = pSysparm->w;
  const T f0 = pApertureData->m_f0;

  const size_t nMaxAbcissa = pSysparm->nDivA;

  const T* axs      = pGL->xs;
  const T* aweights = pGL->ws;

  T scale = pSysparm->rho * c * a / T(M_PI);

  scale *= T(M_PI_2);

  T distNear;
  T distFar;

  T fSampleStart, fSampleStop;

  for (size_t iPoint = pThreadArg->iOutputBegin ;
       iPoint < pThreadArg->iOutputEnd ; iPoint++) {
    // Should be aligned as type also
    alignas(sizeof(T)*4) sps::point_t<T> point;

    point[0] = pos[iPoint*3 + 0];
    point[1] = pos[iPoint*3 + 1];
    point[2] = pos[iPoint*3 + 2];

    // Compute response - use a function instead
    T r, z;
    sps::dist_point_to_circle_local(point, element.circle,
                                    &r, &z, &distNear, &distFar);

    T z2 = SQUARE(z);
    T r2 = SQUARE(r);

    T raz2 = r2 + a2 + z2;

    fSampleStart = fs * (distNear / c);
    fSampleStop  = fs * (distFar / c);

    int iSampleStart = (int)fSampleStart + 1;

    int iSampleStop  = (int) (fSampleStop + W*fs) + 1;

    // Solve integral using Gauss-Legendre
    T spatial = T(0.0);
    for (size_t ia = 0 ; ia < nMaxAbcissa ; ia++) {
      // Can be done once and for all
      T psi = axs[ia];
      T weight = aweights[ia];

      T sqarg = r2 + a2 - T(2.0)*a*r*cos(psi);

      T factor = (r*cos(psi) - a) / sqarg;
      spatial += weight*factor;
    }

    T t1, t2, t3;

    t1 = z / c;
    t2 = sqrt(z2 + SQUARE((r-a))) / c;
    t3 = sqrt(z2 + SQUARE((r+a))) / c;

    // Two terms for tone burst
    T E[A<T>::kTerms];
    for (size_t iTerm = 0 ; iTerm < A<T>::kTerms ; iTerm++) {
      E[iTerm] = T(0.0);
    }

    T psi_1 = T(0.0), psi_2 = T(0.0);

    T psi, weight;

    bool bedge = true;

    size_t i_1 = 0, i_2 = 0;

    // Unit is Pa

    for (int iSample = iSampleStart ; iSample < iSampleStop ; iSample++) {
      bedge = true;
      int iSampleSignal = iSample - iSampleSignalStart;
      T t = invfs * iSample;

      T direct = - A<T>::Evaluate(t-t1, W, f0) * spatial;  // -v(t-tau_2)

      if (t < t2 + W) {
        psi_1 = T(0.0);
      } else if (t < t3 + W) {
        T arg = (raz2 - SQUARE(c)*SQUARE(t-W)) / (T(2.0)*a*std::max<T>(r, eps));
        arg = sps::clamp<T>(arg, T(-1.0), T(1.0));
        psi_1 = acos(arg);
      } else {
        psi_1 = T(0.0);
        bedge = false;
      }

      if ((t >= t3) && (t < (t3+W))) {
        psi_2 = T(M_PI);
      } else if (t > t2) {
        T arg = (raz2 - SQUARE(c*t)) / (T(2.0)*a*std::max<T>(r, eps));
        arg = sps::clamp<T>(arg, T(-1.0), T(1.0));
        psi_2 = acos(arg);
      } else {
        psi_2 = T(0.0);
        bedge = false;
      }

      T edge = T(0.0);
      if (bedge) {
        // Remove any psi
        for (size_t ia = i_1 ; (ia < i_2 && ia < nMaxAbcissa) ; ia++) {
          psi = axs[ia];
          if (psi > psi_1) {
            break;
          }
          i_1 = ia + 1;

          weight = aweights[ia];

          T sqarg = r2+a2-2*a*r*cos(psi);
          T factor = (r*cos(psi) - a) / sqarg;

          T tau = sqrt(std::max<T>(z2 + sqarg, T(0.0))) / c;  // tau_1

          for (size_t iTerm = 0 ; iTerm < A<T>::kTerms ; iTerm++) {
            E[iTerm] -=
              (weight * factor * A<T>::SpatialBasisFunction[iTerm](tau, W, f0));
          }
        }
        // Add new psi
        for (size_t ia = i_2 ; ia < nMaxAbcissa ; ia++) {
          psi = axs[ia];
          if (psi > psi_2) {
            break;
          }
          i_2 = ia + 1;
          weight = aweights[ia];

          T sqarg = r2+a2-T(2.0)*a*r*cos(psi);
          T factor = (r*cos(psi) - a) / sqarg;

          T tau = sqrt(std::max<T>(z2+sqarg, T(0.0))) / c;  // tau_1

          for (size_t iTerm = 0 ; iTerm < A<T>::kTerms ; iTerm++) {
            E[iTerm] +=
              (weight * factor * A<T>::SpatialBasisFunction[iTerm](tau, W, f0));
          }
        }

        for (size_t iTerm = 0 ; iTerm < A<T>::kTerms ; iTerm++) {
          edge +=
            E[iTerm] * A<T>::TemporalBasisFunction[iTerm](t, W, f0);
        }
      }

      T value = scale * (direct + edge);
      odata[iSampleSignal + iPoint * nSamples] = value;
    }

    if (pThreadArg->threadId == 0) {
      // Check performance after 9 points
      if (iPoint == 9) {
        duration = sps::profiler::time() - tProfStart;
        nCallbackPeriod = (size_t) (10.0 / duration);
        nCallbackPeriod = std::max<size_t>(nCallbackPeriod, 1);
      }

      if ( (iPoint > 9) && (iPoint % nCallbackPeriod == 0) ) {
        float val = 100.0f * static_cast<float>(iPoint) / nPositions;
        pQueue->push(val);
# ifdef __GNUC__
        sched_yield();
# elif defined(_WIN32)
        SwitchToThread();
# endif
      }
    }
    if (g_bExit) {
      break;
    }
  }
  debug_print("Thread done\n");

  // Will always happen, even for nPoints < 9
  if (pThreadArg->threadId == 0) {
    pQueue->push(100.0f);
  }
#if HAVE_PTHREAD_H
  return nullptr;
#else
  return 0;
#endif
}

template <class T, template <typename> class A>
T CalcCircularTransientFieldRefZeroDelay(
  const fnm::sysparm_t<T>* sysparm,
  const fnm::CircularApertureData<T>& data,
  const T* pos, const size_t nPositions,
  T** odata, size_t* nSamples) {
  // Reset output
  *odata = nullptr;
  *nSamples = 0;

  const T c  = sysparm->c;
  const T fs = sysparm->fs;
  const T W  = sysparm->w;

  const T f0 = data.m_f0;

  const T invfs = T(1.0) / fs;

  const T eps = std::numeric_limits<T>::epsilon();

  const sps::element_circular_t<T>& element = data.m_element;

  // Simple box surrounding the scatters
  sps::bbox_t<T> scatter_box = sps::bbox_t<T>();
  sps::compute_bounding_box3(pos, nPositions, &scatter_box);

  // Box surrounding circular aperture
  sps::bbox_t<T> circle_box = sps::bbox_t<T>();
  sps::compute_bounding_box_circle<T>(element.circle, &circle_box);

  T distNear;
  T distFar;
  debug_print("scatter_box: %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f\n", scatter_box.min[0], scatter_box.min[1], scatter_box.min[2],
              scatter_box.max[0], scatter_box.max[1], scatter_box.max[2]);
  debug_print("circle_box: %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f\n", circle_box.min[0], circle_box.min[1], circle_box.min[2],
              circle_box.max[0], circle_box.max[1], circle_box.max[2]);
  dists_most_distant_and_closest(scatter_box, circle_box, &distNear, &distFar);

  T tStart = (distNear / c);
  T tEnd   = (distFar / c) + W;

  T fSampleSignalStart = fs * tStart;

  int iSampleSignalStart = (int) floor(fSampleSignalStart);

  // Allocate output
  int _nSamples = (int) ceil(T(2.0) + fs*(tEnd - tStart));

  *odata = static_cast<T*>(malloc(nPositions*_nSamples*sizeof(T)));
  memset(*odata, 0, nPositions*_nSamples*sizeof(T));
  *nSamples = _nSamples;

  const T a = element.circle.radius;
  const T a2 = SQUARE(a);

  const T s1 = 0;
  const T s2 = T(M_PI);
  const size_t nMaxAbcissa = sysparm->nDivA;

  sps::unique_aligned_array<T> aweights;
  sps::unique_aligned_array<T> axs;

  CalcOneDimensionalWeightsAndAbcissae(nMaxAbcissa, s1, s2, std::move(axs), std::move(aweights));

  T fSampleStart, fSampleStop;

  T scale = sysparm->rho * c * a / T(M_PI);

  scale *= T(M_PI_2);

  for (size_t iPosition = 0; iPosition < nPositions; iPosition++) {
    sps::point_t<T> point;
    point[0] = pos[iPosition * 3];
    point[1] = pos[iPosition * 3 + 1];
    point[2] = pos[iPosition * 3 + 2];

    T r, z;
    sps::dist_point_to_circle_local(point, element.circle,
                                    &r, &z, &distNear, &distFar);

    T z2 = SQUARE(z);
    T r2 = SQUARE(r);

    T raz2 = r2 + a2 + z2;

    fSampleStart = fs * (distNear / c);
    fSampleStop  = fs * (distFar / c);

    int iSampleStart = (int)fSampleStart + 1;

    int iSampleStop  = (int) (fSampleStop + W*fs) + 1;

    // Solve integral using Gauss-Legendre
    T spatial = T(0.0);
    for (size_t ia = 0 ; ia < nMaxAbcissa ; ia++) {
      // Can be done once and for all
      T psi = axs[ia];
      T weight = aweights[ia];

      T sqarg = r2 + a2 - T(2.0)*a*r*cos(psi);

      T factor = (r*cos(psi) - a) / sqarg;
      spatial += weight*factor;
    }

    T t1, t2, t3;

    t1 = z / c;
    t2 = sqrt(z2 + SQUARE((r-a))) / c;
    t3 = sqrt(z2 + SQUARE((r+a))) / c;

    // Two terms for tone burst
    T E[A<T>::kTerms];
    for (size_t iTerm = 0 ; iTerm < A<T>::kTerms ; iTerm++) {
      E[iTerm] = T(0.0);
    }

    T psi_1 = T(0.0), psi_2 = T(0.0);

    T psi, weight;

    bool bedge = true;

    size_t i_1 = 0, i_2 = 0;

    // Unit is Pa

    for (int iSample = iSampleStart ; iSample < iSampleStop ; iSample++) {
      bedge = true;
      int iSampleSignal = iSample - iSampleSignalStart;
      T t = invfs * iSample;

      T direct = - A<T>::Evaluate(t-t1, W, f0) * spatial;  // -v(t-tau_2)

      if (t < t2 + W) {
        psi_1 = T(0.0);
      } else if (t < t3 + W) {
        T arg = (raz2 - SQUARE(c)*SQUARE(t-W)) / (T(2.0)*a*std::max<T>(r, eps));
        arg = sps::clamp<T>(arg, T(-1.0), T(1.0));
        psi_1 = acos(arg);
      } else {
        psi_1 = T(0.0);
        bedge = false;
      }

      if ((t >= t3) && (t < (t3+W))) {
        psi_2 = T(M_PI);
      } else if (t > t2) {
        T arg = (raz2 - SQUARE(c*t)) / (T(2.0)*a*std::max<T>(r, eps));
        arg = sps::clamp<T>(arg, T(-1.0), T(1.0));
        psi_2 = acos(arg);
      } else {
        psi_2 = T(0.0);
        bedge = false;
      }

      T edge = T(0.0);
      if (bedge) {
        // Remove any psi
        for (size_t ia = i_1 ; (ia < i_2 && ia < nMaxAbcissa) ; ia++) {
          psi = axs[ia];
          if (psi > psi_1) {
            break;
          }
          i_1 = ia + 1;

          weight = aweights[ia];

          T sqarg = r2+a2-2*a*r*cos(psi);
          T factor = (r*cos(psi) - a) / sqarg;

          T tau = sqrt(std::max<T>(z2 + sqarg, T(0.0))) / c;  // tau_1

          for (size_t iTerm = 0 ; iTerm < A<T>::kTerms ; iTerm++) {
            E[iTerm] -=
              (weight * factor * A<T>::SpatialBasisFunction[iTerm](tau, W, f0));
          }
        }
        // Add new psi
        for (size_t ia = i_2 ; ia < nMaxAbcissa ; ia++) {
          psi = axs[ia];
          if (psi > psi_2) {
            break;
          }
          i_2 = ia + 1;
          weight = aweights[ia];

          T sqarg = r2+a2-T(2.0)*a*r*cos(psi);
          T factor = (r*cos(psi) - a) / sqarg;

          T tau = sqrt(std::max<T>(z2+sqarg, T(0.0))) / c;  // tau_1

          for (size_t iTerm = 0 ; iTerm < A<T>::kTerms ; iTerm++) {
            E[iTerm] +=
              (weight * factor * A<T>::SpatialBasisFunction[iTerm](tau, W, f0));
          }
        }

        for (size_t iTerm = 0 ; iTerm < A<T>::kTerms ; iTerm++) {
          edge +=
            E[iTerm] * A<T>::TemporalBasisFunction[iTerm](t, W, f0);
        }
      }

      T value = scale * (direct + edge);
      (*odata)[iSampleSignal + iPosition*_nSamples] = value;
    }
  }

  return 0;
}

template <class T>
int CalcCircularCwFieldRef(const fnm::sysparm_t<T>* sysparm,
                           const fnm::CircularApertureData<T>& data,
                           const T* pos, const size_t nPositions,
                           std::complex<T>** odata) {
  int retval = 0;

  // Allocate output
  *odata =
    static_cast<std::complex<T>*>(SPS_MALLOC(nPositions*sizeof(std::complex<T>)));
  memset(*odata, 0, nPositions*sizeof(std::complex<T>));

  const T eps = std::numeric_limits<T>::epsilon();

  const T a = data.m_element.circle.radius;
  const T a2 = SQUARE(a);

  const T lambda = sysparm->c / data.m_f0;
  const T k = T(M_2PI)/lambda;

  const T scale = (a/(k*T(M_PI)));

  const size_t nMaxSectors = sysparm->nMaxSectors;

  size_t nMaxAbcissa = sysparm->nDivA;
  T S = data.m_gridSectorScale;

  size_t nSectors = nMaxSectors;

  if (S == T(1.0)) {
    nSectors = 1;
  }

  // Number of abcissas, M should be chosen according to which sector a point belong to
  //
  // r/z = tan(arcsin((i-1)/ (n_l-1))) for i=1 to n_l
  //
  // There is n_s = n_l-1 sectors (n_l is number of lines dividing sectors).

  // M = N_max(1-S) / (n_s ) * (i-n_s) + N_max for i=1 to n_s

  std::vector<sps::unique_aligned_array<T> > aweights(nSectors);
  std::vector<sps::unique_aligned_array<T> > axs(nSectors);
  std::vector<size_t> nAbcissas(nSectors);

  const T s1 = 0;
  const T s2 = T(M_PI);

  for (size_t iSector = 0; iSector < nSectors; iSector++) {
    T fAbcissas = nMaxAbcissa * (T(1.0) - S) / nSectors * (T(iSector) + 1 - nSectors) + nMaxAbcissa;
    size_t nAbcissa = (size_t)fAbcissas;
    CalcOneDimensionalWeightsAndAbcissae(
      nAbcissa,
      s1, s2,
      std::move(axs[iSector]), std::move(aweights[iSector]));
    nAbcissas[iSector] = nAbcissa;
  }

  for (size_t iPosition = 0 ; iPosition < nPositions ; iPosition++) {
    sps::point_t<T> point;
    point[0] = pos[iPosition*3];
    point[1] = pos[iPosition*3+1];
    point[2] = pos[iPosition*3+2];

    T field_real = T(0.0);
    T field_imag = T(0.0);

    // No loop, there is only a single element
    const sps::element_circular_t<T>& element = data.m_element;

    sps::point_t<T> r2p = point - element.circle.center;

    sps::point_t<T> normal = sps::point_t<T>();
    sps::basis_vectors<T, sps::EulerIntrinsicYXY>(&normal, element.circle.euler, 2);

    T z = dot(normal, r2p);
    T z2 = SQUARE(z);
    T r = sqrt(dot(r2p, r2p) - z2);
    T r2 = SQUARE(r);

    // Compute number of relevant abcissas
    size_t iSector =
      (size_t)floor(std::max<T>(sin(atan2(r, z))*nSectors - eps, T(0.0)));

    // If S==1.0, we only use a single sector, iSector = 0
    if (S == T(1.0)) {
      iSector = 0;
    }

    T carg = cos(-k*z);  // k = \omega / c
    T sarg = sin(-k*z);

    T real, imag;

    // Solve integral using Gauss-Legendre
    for (size_t ia = 0 ; ia < nAbcissas[iSector] ; ia++) {
      // Can be done once and for all
      T psi = axs[iSector][ia];
      T weight = aweights[iSector][ia];

      T sqarg = r2 + a2 - 2*a*r*cos(psi);
      T factor = (r*cos(psi) - a) / sqarg;

      T sqr = sqrt(sqarg+z2);

      real = weight * factor * (cos(-k*sqr) - carg);
      imag = weight * factor * (sin(-k*sqr) - sarg);

      // No phase supported yet
      field_real += real;
      field_imag += imag;
    }

    // Multiply with -i
    std::swap(field_real, field_imag);
    field_imag = -field_imag;

    // Scale with radius and k
    field_real = field_real * scale;
    field_imag = field_imag * scale;

    (*odata)[iPosition] = std::complex<T>(field_real, field_imag);
  }
  return retval;
}

#if __GNUG__ || defined(_MSC_VER) && (_MSC_VER >= 1900)
template <class T, template <typename> class A>
T CalcCircularPwThreaded(
  const fnm::sysparm_t<T>* pSysparm,
  const fnm::CircularApertureData<T>& data,
  const T* pos, const size_t nPositions,
  T** odata, size_t* nSamples) {
  // Reset output
  *odata = nullptr;
  *nSamples = 0;

  const T c  = pSysparm->c;
  const T fs = pSysparm->fs;
  const T W  = pSysparm->w;

  const sps::element_circular_t<T>& element = data.m_element;

  // Simple box surrounding the scatters
  sps::bbox_t<T> scatter_box = sps::bbox_t<T>();
  sps::compute_bounding_box3(pos, nPositions, &scatter_box);

  // Box surrounding circular aperture
  sps::bbox_t<T> circle_box = sps::bbox_t<T>();
  sps::compute_bounding_box_circle<T>(element.circle, &circle_box);

  T distNear;
  T distFar;
  debug_print("scatter_box: %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f\n", scatter_box.min[0], scatter_box.min[1], scatter_box.min[2],
              scatter_box.max[0], scatter_box.max[1], scatter_box.max[2]);
  debug_print("circle_box: %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f\n", circle_box.min[0], circle_box.min[1], circle_box.min[2],
              circle_box.max[0], circle_box.max[1], circle_box.max[2]);
  dists_most_distant_and_closest(scatter_box, circle_box, &distNear, &distFar);

  T tStart = (distNear / c);
  T tEnd   = (distFar / c) + W;

  T fSampleSignalStart = fs * tStart;

  int iSampleSignalStart = (int) floor(fSampleSignalStart);

  // Allocate output
  int _nSamples = (int) ceil(T(2.0) + fs*(tEnd - tStart));

  *odata = static_cast<T*>(SPS_MALLOC(nPositions*_nSamples*sizeof(T)));
  memset(*odata, 0, nPositions*_nSamples*sizeof(T));
  *nSamples = _nSamples;

  const T s1 = 0;
  const T s2 = T(M_PI);
  const size_t nMaxAbcissa = pSysparm->nDivA;

  sps::unique_aligned_array<T> aweights;
  sps::unique_aligned_array<T> axs;

  CalcOneDimensionalWeightsAndAbcissae(
    nMaxAbcissa, s1, s2, std::move(axs), std::move(aweights));

  GLQuad1D<T> GL;
  GL.xs = axs.get();
  GL.ws = aweights.get();
  GL.nx = nMaxAbcissa;

  int nproc = getncpus();

  sps::queue<float>* progressQueue = new sps::queue<float>();

  PwFnmThreadArg<T> threadArg[N_MAX_THREADS];

  size_t nThreads = pSysparm->nThreads;

  for (size_t iThread = 0; iThread < nThreads; ++iThread) {
    threadArg[iThread].pSysparm         = pSysparm;
    threadArg[iThread].pApertureData    = &data;
    threadArg[iThread].pGL              = &GL;
    threadArg[iThread].pPos             = pos;
    threadArg[iThread].iOutputBegin = 0 + iThread*(nPositions/nThreads);
    threadArg[iThread].iOutputEnd   =
      (nPositions/nThreads) + iThread*(nPositions/nThreads);
    threadArg[iThread].nPositions = nPositions;
    threadArg[iThread].pField       = *odata;

    threadArg[iThread].nSamples = _nSamples;
    threadArg[iThread].iSampleStart = iSampleSignalStart;

    if (iThread == 0) {
      threadArg[iThread].pQueue    = progressQueue;
    }

    threadArg[iThread].threadId    = iThread;
    threadArg[iThread].cpuId       = ((int) iThread) % nproc;
    if (iThread == nThreads-1) {
      threadArg[iThread].iOutputEnd = nPositions;
    }
  }
#ifdef _MSC_VER
  std::vector<sps::ThreadPool::TaskFuture<unsigned int> > futures;
  for (size_t iThread = 0 ; iThread < nThreads ; iThread++) {
    sps::ThreadPool::TaskFuture<unsigned int> taskFuture =
    sps::thread::defaultpool::submitJob([&](void* pThreadArg) mutable {
      return CalcCircularPwThreadFunc<T, A>(pThreadArg);
    }, static_cast<void*>(&threadArg[iThread]));
    futures.push_back(std::move(taskFuture));
  }
#else
  std::vector<sps::ThreadPool::TaskFuture<void*> > futures;
  for (size_t iThread = 0 ; iThread < nThreads ; iThread++) {
    sps::ThreadPool::TaskFuture<void*> taskFuture =
    sps::thread::defaultpool::submitJob([&](void* pThreadArg) mutable {
      return CalcCircularPwThreadFunc<T, A>(pThreadArg);
    }, static_cast<void*>(&threadArg[iThread]));
    futures.push_back(std::move(taskFuture));
  }
#endif
  for  (size_t iThread = 0 ; iThread < nThreads ; iThread++) {
    auto _dummy = futures[iThread].Get();
    SPS_UNREFERENCED_PARAMETER(_dummy);
  }

  delete progressQueue;

  // Assign output
  *nSamples = _nSamples;
  return tStart;
}
#endif

template
int CalcCircularCwFieldRef(const sysparm_t<float>* sysparm,
                           const CircularApertureData<float>& data,
                           const float* pos, const size_t nPositions,
                           std::complex<float>** odata);

template
float CalcCircularTransientFieldRefZeroDelay<float, ToneBurst>(
  const sysparm_t<float>* sysparm,
  const CircularApertureData<float>& data,
  const float* pos, const size_t nPositions,
  float** odata, size_t* nSamples);

template
float CalcCircularTransientFieldRefZeroDelay<float, HanningWeightedPulse>(
  const sysparm_t<float>* sysparm,
  const CircularApertureData<float>& data,
  const float* pos, const size_t nPositions,
  float** odata, size_t* nSamples);

template float CalcCircularPwThreaded<float, ToneBurst>(
  const fnm::sysparm_t<float>* pSysparm,
  const fnm::CircularApertureData<float>& data,
  const float* pos, const size_t nPositions,
  float** odata, size_t* nSamples);

template
float CalcCircularPwThreaded<float, HanningWeightedPulse>(
  const fnm::sysparm_t<float>* pSysparm,
  const fnm::CircularApertureData<float>& data,
  const float* pos, const size_t nPositions,
  float** odata, size_t* nSamples);

#ifdef FNM_DOUBLE_SUPPORT
template
int CalcCircularCwFieldRef(const sysparm_t<double>* sysparm,
                           const CircularApertureData<double>& data,
                           const double* pos, const size_t nPositions,
                           std::complex<double>** odata);

template
double CalcCircularTransientFieldRefZeroDelay<double, ToneBurst>(
  const sysparm_t<double>* sysparm,
  const CircularApertureData<double>& data,
  const double* pos, const size_t nPositions,
  double** odata, size_t* nSamples);

template
double CalcCircularTransientFieldRefZeroDelay<double, HanningWeightedPulse>(
  const sysparm_t<double>* sysparm,
  const CircularApertureData<double>& data,
  const double* pos, const size_t nPositions,
  double** odata, size_t* nSamples);

template
double CalcCircularPwThreaded<double, ToneBurst>(
  const fnm::sysparm_t<double>* pSysparm,
  const fnm::CircularApertureData<double>& data,
  const double* pos, const size_t nPositions,
  double** odata, size_t* nSamples);

template
double CalcCircularPwThreaded<double, HanningWeightedPulse>(
  const fnm::sysparm_t<double>* pSysparm,
  const fnm::CircularApertureData<double>& data,
  const double* pos, const size_t nPositions,
  double** odata, size_t* nSamples);

#endif
}  // namespace fnm
