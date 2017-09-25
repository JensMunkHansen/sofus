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

#include <fnm/circular_data.hpp>
#include <fnm/circular_calc.hpp>
#include <fnm/fnm_basis.hpp> // Basis functions

#include <gl/gl.hpp>

#include <sps/algorithm>
#include <sps/debug.h>

#include <vector>

// TODO: Work on version supporting delay: all t -> t - t_offset

#define CIRCULAR_FAST_EN_INTEGRALS 1

namespace fnm {
  template <class T>
  void CalcOneDimensionalWeightsAndAbcissae(const size_t N,
      const T begin,
      const T end,
      sps::deleted_aligned_array<T> &&uxs,
      sps::deleted_aligned_array<T> &&uweights)
  {

    uxs      = sps::deleted_aligned_array_create<T>(N);
    uweights = sps::deleted_aligned_array_create<T>(N);

    const T sm = T(0.5) * (end - begin);
    const T sp = T(0.5) * (end + begin);

    // Common weights and abcissa
    for (size_t i = 0 ; i < N ; i++) {
      gl::GLNode qp = gl::GL(N,i);
      // Conversion from double to T
      uxs[i]        = sm * T(qp.value) + sp;
      uweights[i]   = T(qp.weight);
    }
  }

  template <class T, template <typename> class A>
  T CalcCircularTransientFieldRefZeroDelay(const fnm::sysparm_t<T>* sysparm,
      const fnm::CircularApertureData<T>& data,
      const T* pos, const size_t nPositions,
      T** odata, size_t* nSamples)
  {
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

    *odata = (T*) malloc(nPositions*_nSamples*sizeof(T));
    memset(*odata, 0, nPositions*_nSamples*sizeof(T));
    *nSamples = _nSamples;

    const T a = element.circle.radius;
    const T a2 = SQUARE(a);

    const T s1 = 0;
    const T s2 = T(M_PI);
    const size_t nMaxAbcissa = sysparm->nDivA;

    sps::deleted_aligned_array<T> aweights;
    sps::deleted_aligned_array<T> axs;

    CalcOneDimensionalWeightsAndAbcissae(nMaxAbcissa, s1, s2, std::move(axs), std::move(aweights));

    T fSampleStart, fSampleStop;

    const T rho = T(1000.0); // [kg/m3]

    T scale = rho * c * a / T(M_PI);

    // If we multiply with an additional T(M_PI_2) then it matches the results produced using SOFUS
    scale *= T(M_PI_2);

    for (size_t iPosition = 0; iPosition < nPositions; iPosition++) {
      sps::point_t<T> point;
      point[0] = pos[iPosition * 3];
      point[1] = pos[iPosition * 3 + 1];
      point[2] = pos[iPosition * 3 + 2];

      T r, z;
      sps::dist_point_to_circle_local(point, element.circle, &r, &z, &distNear, &distFar);

      T z2 = SQUARE(z);
      T r2 = SQUARE(r);

      T raz2 = r2 + a2 + z2;

      fSampleStart = fs * (distNear / c);
      fSampleStop  = fs * (distFar / c);

      int iSampleStart = (int)fSampleStart + 1; // First non-zero value

      int iSampleStop  = (int) (fSampleStop + W*fs) + 1; // Ends with a non-zero value

      // Missing factor (rho_0 c a)/pi

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

      t1 = z / c; // (=tau_2)
      t2 = sqrt(z2 + SQUARE((r-a))) / c;
      t3 = sqrt(z2 + SQUARE((r+a))) / c;

      // Two terms for tone burst
      T E[A<T>::nTerms];
      for (size_t iTerm = 0 ; iTerm < A<T>::nTerms ; iTerm++) {
        E[iTerm] = T(0.0);
      }

      T psi_1 = T(0.0), psi_2=T(0.0);

      T psi, weight;

      bool bedge = true;

      size_t i_1 = 0, i_2 = 0;

      // Unit is Pa

      for (int iSample = iSampleStart ; iSample < iSampleStop ; iSample++) {

        bedge = true;
        int iSampleSignal = iSample - iSampleSignalStart;
        T t = invfs * iSample;

        T direct = - A<T>::Evaluate(t-t1, W, f0) * spatial; // -v(t-tau_2)

        if (t < t2 + W) {
          psi_1 = T(0.0);
        } else if (t < t3 + W) {
          T arg = (raz2 - SQUARE(c)*SQUARE(t-W)) / (T(2.0)*a*std::max<T>(r,eps));
          arg = std::clamp<T>(arg,T(-1.0),T(1.0));
          psi_1 = acos(arg);
        } else {
          psi_1 = T(0.0);
          bedge = false;
        }

        if ( (t >= t3) && (t < (t3+W))) {
          psi_2 = T(M_PI);
        } else if ( (t > t2)) {
          T arg = (raz2 - SQUARE(c*t)) / (T(2.0)*a*std::max<T>(r,eps));
          arg = std::clamp<T>(arg,T(-1.0),T(1.0));
          psi_2 = acos(arg);
        } else {
          psi_2 = T(0.0);
          bedge = false;
        }

#if CIRCULAR_FAST_EN_INTEGRALS

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

            T tau = sqrt(std::max<T>(z2 + sqarg,T(0.0))) / c; // tau_1

            for (size_t iTerm = 0 ; iTerm < A<T>::nTerms ; iTerm++) {
              E[iTerm] -= (weight * factor * A<T>::SpatialBasisFunction[iTerm](tau,W,f0));
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

            T tau = sqrt( std::max<T>(z2+sqarg, T(0.0))) / c; // tau_1

            for (size_t iTerm = 0 ; iTerm < A<T>::nTerms ; iTerm++) {
              E[iTerm] += (weight * factor * A<T>::SpatialBasisFunction[iTerm](tau,W,f0));
            }
          }

          for (size_t iTerm = 0 ; iTerm < A<T>::nTerms ; iTerm++) {
            edge += E[iTerm] * A<T>::TemporalBasisFunction[iTerm](t,W,f0);
          }
        }
#else
        // Works but is slower than reusing terms
        T edge = T(0.0);
        E[0] = T(0.0);
        E[1] = T(0.0);
        if (bedge) {

          for (size_t ia = 0 ; ia < nMaxAbcissa ; ia++) {
            psi = axs[ia];
            weight = aweights[ia];
            if ( (psi < psi_1) || (psi > psi_2) ) {
              continue;
            }

            T sqarg = SQUARE(r)+SQUARE(a)-T(2.0)*a*r*cos(psi);
            T factor = (r*cos(psi) - a) / sqarg;

            T tau = sqrt( std::max<T>(z2+sqarg, T(0.0))) / c; // tau_1
            for (size_t iTerm = 0 ; iTerm < A<T>::nTerms ; iTerm++) {
              E[iTerm] += (weight * factor * A<T>::SpatialBasisFunction[iTerm](tau,W,f0));
            }
          }
          for (size_t iTerm = 0 ; iTerm < A<T>::nTerms ; iTerm++) {
            edge += E[iTerm] * A<T>::TemporalBasisFunction[iTerm](t,W,f0);
          }
        }
#endif
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
                             std::complex<T>** odata)
  {
    int retval = 0;

    // Allocate output
    *odata = (std::complex<T>*) malloc(nPositions*sizeof(std::complex<T>));
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

    std::vector<sps::deleted_aligned_array<T> > aweights(nSectors);
    std::vector<sps::deleted_aligned_array<T> > axs(nSectors);
    std::vector<size_t> nAbcissas(nSectors);

    const T s1 = 0;
    const T s2 = T(M_PI);

    for (size_t iSector = 0; iSector < nSectors; iSector++) {
      T fAbcissas = nMaxAbcissa * (T(1.0) - S) / nSectors * (T(iSector) + 1 - nSectors) + nMaxAbcissa;
      size_t nAbcissa = (size_t)fAbcissas;
      CalcOneDimensionalWeightsAndAbcissae(nAbcissa,
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
      sps::basis_vectors<T, sps::EulerIntrinsicYXY>(normal, element.circle.euler, 2);

      T z = dot(normal,r2p);
      T z2 = SQUARE(z);
      T r = sqrt(dot(r2p,r2p) - z2);
      T r2 = SQUARE(r);

      // Compute number of relevant abcissas
      size_t iSector = (size_t)floor(std::max<T>(sin(atan2(r, z))*nSectors - eps, T(0.0)));

      // If S==1.0, we only use a single sector, iSector = 0
      if (S == T(1.0)) {
        iSector = 0;
      }

      T carg = cos(-k*z); // k = \omega / c
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
      std::swap(field_real,field_imag);
      field_imag = -field_imag;

      // Scale with radius and k
      field_real = field_real * scale;
      field_imag = field_imag * scale;

      (*odata)[iPosition] = std::complex<T>(field_real,field_imag);

    }
    return retval;
  }

  template int CalcCircularCwFieldRef(const sysparm_t<float>* sysparm,
                                      const CircularApertureData<float>& data,
                                      const float* pos, const size_t nPositions,
                                      std::complex<float>** odata);

  template float CalcCircularTransientFieldRefZeroDelay<float, ToneBurst>(const sysparm_t<float>* sysparm,
      const CircularApertureData<float>& data,
      const float* pos, const size_t nPositions,
      float** odata, size_t* nSamples);

  template float CalcCircularTransientFieldRefZeroDelay<float, HanningWeightedPulse>(const sysparm_t<float>* sysparm,
      const CircularApertureData<float>& data,
      const float* pos, const size_t nPositions,
      float** odata, size_t* nSamples);

#ifdef FNM_DOUBLE_SUPPORT
  template int CalcCircularCwFieldRef(const sysparm_t<double>* sysparm,
                                      const CircularApertureData<double>& data,
                                      const double* pos, const size_t nPositions,
                                      std::complex<double>** odata);

  template double CalcCircularTransientFieldRefZeroDelay<double, ToneBurst>(const sysparm_t<double>* sysparm,
      const CircularApertureData<double>& data,
      const double* pos, const size_t nPositions,
      double** odata, size_t* nSamples);

  template double CalcCircularTransientFieldRefZeroDelay<double, HanningWeightedPulse>(const sysparm_t<double>* sysparm,
      const CircularApertureData<double>& data,
      const double* pos, const size_t nPositions,
      double** odata, size_t* nSamples);

#endif
}
