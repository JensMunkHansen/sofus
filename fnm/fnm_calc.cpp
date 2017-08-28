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

#include <fnm/config.h>
#include <fnm/fnm.hpp>
#include <fnm/fnm_common.hpp>
#include <fnm/fnm_calc.hpp>
#include <fnm/FnmSIMD.hpp> // CalcHzFast

#include <gl/gl.hpp>

#include <sps/smath.hpp>
#include <sps/mm_malloc.h>
#include <sps/memory>
#include <sps/sps_threads.hpp>
#include <sps/profiler.h>
#include <sps/progress.hpp>
#include <sps/queue.hpp>
#include <sps/debug.h>

#include <memory>

#if 1
# define CALC_SELECT CalcSingle
#endif

namespace fnm {

#ifdef HAVE_PTHREAD_H
  // External
  pthread_t threads[N_MAX_THREADS] = {};
  pthread_attr_t attr;
#endif

  /////////////////////////////////////////////////
  // Types visible for this compilation unit only
  /////////////////////////////////////////////////
  template <class T>
  struct CwFieldThreadArg {
    const sysparm_t<T>* sysparm;
    const ApertureData<T>* data;
    size_t iPointBegin;
    size_t iPointEnd;
    const T* pos;
    std::complex<T>* field;
    T k;
    T* uxs;
    T* uweights;
    size_t nDivU;
    T* vxs;
    T* vweights;
    size_t nDivV;
    /// Number of positions (all threads)
    size_t nPositions;
    sps::queue<float>* pQueue;
    size_t threadId;
    int cpu_id;
  };

  template <class T>
  void CalcWeightsAndAbcissae(const sysparm_t<T>* sysparm,
                              sps::deleted_aligned_array<T> &&uxs,
                              sps::deleted_aligned_array<T> &&uweights,
                              sps::deleted_aligned_array<T> &&vxs,
                              sps::deleted_aligned_array<T> &&vweights)
  {

    size_t nDivW = sysparm->nDivW;
    size_t nDivH = sysparm->nDivH;

    uxs      = sps::deleted_aligned_array_create<T>(nDivW);
    uweights = sps::deleted_aligned_array_create<T>(nDivW);
    vxs      = sps::deleted_aligned_array_create<T>(nDivH);
    vweights = sps::deleted_aligned_array_create<T>(nDivH);

    // Common weights and abcissa
    for (size_t i = 0 ; i < nDivW ; i++) {
      gl::GLNode qp = gl::GL(nDivW,i);
      // Conversion from double to T
      uxs[i]        = T(qp.value);
      uweights[i]   = T(qp.weight);
    }

    for (size_t i = 0 ; i < nDivH ; i++) {
      gl::GLNode qp = gl::GL(nDivH, i);
      // Conversion from double to T
      vxs[i]      = T(qp.value);
      vweights[i] = T(qp.weight);
    }
  }

#if defined(HAVE_PTHREAD_H)
  template <class T>
  void* CalcCwThreadFunc(void* ptarg)
#else
  template <class T>
  unsigned int __stdcall CalcCwThreadFunc(void *ptarg)
#endif
  {
    fnm::CwFieldThreadArg<T>* pThreadArg = reinterpret_cast<fnm::CwFieldThreadArg<T>*>(ptarg);

#if FNM_ENABLE_ATTENUATION
    const fnm::sysparm_t<T>* sysparm = pThreadArg->sysparm;
#endif
    const ApertureData<T>* m_data    = pThreadArg->data;
    const T k                        = pThreadArg->k;
    const T* uxs                     = pThreadArg->uxs;
    const T* vxs                     = pThreadArg->vxs;
    const T* __restrict uweights     = pThreadArg->uweights;
    const T* __restrict vweights     = pThreadArg->vweights;
    const T* pos                     = pThreadArg->pos;
    const size_t nDivW               = pThreadArg->nDivU;
    const size_t nDivH               = pThreadArg->nDivV;
    std::complex<T>* odata           = pThreadArg->field;

#ifdef HAVE_THREAD
    setcpuid(pThreadArg->cpu_id);
#endif
    const size_t nPositions          = pThreadArg->nPositions;
    sps::queue<float>* pQueue        = pThreadArg->pQueue;

    const size_t nElements    = m_data->m_nelements;
    const size_t nSubElements = m_data->m_nsubelements;

    size_t _nElements, _nSubElements;
    const sps::element_rect_t<T>** elements = NULL;
    m_data->ElementsRefGet(&_nElements, &_nSubElements, elements);

    const T* apodizations = m_data->m_apodizations.get();

    T apodization;

    ALIGN16_BEGIN sps::point_t<T> projection ALIGN16_END;

#ifdef _MSC_VER
    debug_print("iPointBegin: %Iu, iPointEnd: %Iu\n", pThreadArg->iPointBegin, pThreadArg->iPointEnd);
#else
    debug_print("iPointBegin: %zu, iPointEnd: %zu\n", pThreadArg->iPointBegin, pThreadArg->iPointEnd);
#endif

#if FNM_ENABLE_ATTENUATION
    T alpha = sysparm->att;
    if (!(sysparm->use_att)) {
      alpha = T(0.0);
    }
#endif

    // Initial time
    double tProfStart = 0.0;

    if (pThreadArg->threadId == 0) {
      tProfStart = sps::profiler::time();
    }

    size_t nCallbackPeriod = 1;
    double duration = 0.0;

    for (size_t iPoint = pThreadArg->iPointBegin ; iPoint < pThreadArg->iPointEnd ; iPoint++) {

#ifdef FNM_DOUBLE_SUPPORT
      sps::point_t<T> point;
      point[0] = pos[iPoint*3];
      point[1] = pos[iPoint*3+1];
      point[2] = pos[iPoint*3+2];
#else
      __m128 vec_point =
        _mm_set_ps(0.0f,
                   float(pos[iPoint*3 + 2]),
                   float(pos[iPoint*3 + 1]),
                   float(pos[iPoint*3 + 0]));
#endif

      std::complex<T> final = std::complex<T>(T(0.0),T(0.0));

      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {

        apodization = apodizations[iElement];

        if (apodization != 0.0) {

          for (size_t jElement = 0 ; jElement <  m_data->m_nsubelements ; jElement++) {

            const auto& element = elements[iElement][jElement];
            std::complex<T> result;
            assert( ((uintptr_t)&element.center[0] & 0xF) == 0);
            assert( ((uintptr_t)&element.uvector[0] & 0xF) == 0);
            assert( ((uintptr_t)&element.vvector[0] & 0xF) == 0);

#ifdef FNM_DOUBLE_SUPPORT
            // Scalar implementation
            sps::point_t<T> r2p = point - element.center;
            sps::point_t<T> hh_dir,hw_dir,normal;

            sps::basis_vectors<T, sps::EulerIntrinsicYXY>(hw_dir, element.euler, 0);
            sps::basis_vectors<T, sps::EulerIntrinsicYXY>(hh_dir, element.euler, 1);
            sps::basis_vectors<T, sps::EulerIntrinsicYXY>(normal, element.euler, 2);

            // Projection onto plane
            projection[0] = dot(hw_dir,r2p);
            projection[1] = dot(hh_dir,r2p);
            projection[2] = dot(normal,r2p);
#else
            __m128 vec_r2p = _mm_sub_ps(
                               vec_point,
                               _mm_load_ps((float*)&element.center[0]));
            _mm_store_ss((float*)&projection[0],
                         _mm_dp_ps(_mm_load_ps(&element.uvector[0]),
                                   vec_r2p,0x71));
            _mm_store_ss((float*)&projection[1],
                         _mm_dp_ps(_mm_load_ps(&element.vvector[0]),
                                   vec_r2p,0x71));
            _mm_store_ss((float*)&projection[2],
                         _mm_fabs_ps(_mm_dp_ps(_mm_load_ps(&element.normal[0]),
                                               vec_r2p,0x71)));
#endif

            // TODO: Use CalcHzFast for comparison with reference
#if 1
            // This is worse not reducing the number of integrals
            result = CalcHzFast<T>(element, projection, k,
                                   uxs, uweights, nDivW,
                                   vxs, vweights, nDivH);
#elif 1
            // This is more correct - number of integral are reduced to 4 (but do not match Python reference)
            // Any works for all nDivH, nDivW. Any2 works for nDivH/nDivW even

            // TODO: Make CalcFastFour call these two and use GLQuad2D
            //       Make CalcCwFourRef match this
            result = CalcFastFourAny2<T>(projection[0],
                                         projection[1],
                                         element.hw,
                                         element.hh,
                                         projection[2],
                                         k,
                                         uxs,
                                         uweights, nDivW);
            result += CalcFastFourAny2<T>(projection[1],
                                          projection[0],
                                          element.hh,
                                          element.hw,
                                          projection[2],
                                          k,
                                          vxs,
                                          vweights, nDivH);
#endif
            T real = apodization * result.real();
            T imag = apodization * result.imag();

            T arg = m_data->m_phases[iElement*nSubElements+jElement];
            T carg = cos(arg);
            T sarg = sin(arg);

#if FNM_ENABLE_ATTENUATION

# ifdef FNM_DOUBLE_SUPPORT
            T dist = norm(r2p);
# else
            T dist = (T) _mm_cvtss_f32(_mm_sqrt_ps(_mm_dp_ps(vec_r2p,vec_r2p,0x71)));
# endif

            // Valgrind reports uninitialized variable (must be a false positive)
            T factor = exp(-dist*alpha);

            final.real(final.real() + factor*(real*carg - imag*sarg));
            final.imag(final.imag() + factor*(real*sarg + imag*carg));
#else
            final.real(final.real() + real*carg - imag*sarg);
            final.imag(final.imag() + real*sarg + imag*carg);
#endif
          }
        }
      }
      odata[iPoint] = final;

      if (pThreadArg->threadId==0) {
        // Check time
        if (iPoint == 9) {
          duration = sps::profiler::time() - tProfStart;
          nCallbackPeriod = (size_t) (10.0 / duration);
        }

        if ( (iPoint > 9) && (iPoint % nCallbackPeriod == 0) ) {
          float val = 100.0f * float(iPoint) / nPositions;
          pQueue->push(val);
# ifdef __GNUC__
          sched_yield();
# elif defined(_WIN32)
          SwitchToThread();
# endif
        }
      }
    }
    debug_print("Thread done\n");

    // Will always happen, even for nPoints < 9
    if (pThreadArg->threadId == 0) {
      pQueue->push(100.0f);
    }
#if HAVE_PTHREAD_H
# ifdef HAVE_MQUEUE_H
    return NULL;
# else
    pthread_exit(NULL);
# endif
#else
    return 0;
#endif
  }

  /**
   * Computation of the individual integrals in the reference implementation @ref CalcCwFieldRef
   *
   * @param s1
   * @param s2
   * @param l
   * @param z
   * @param k
   * @param uxs
   * @param uweights
   * @param nUs
   *
   * @return
   */
  template <class T>
  STATIC_INLINE_BEGIN
  std::complex<T> CalcSingle(const T& s1,
                             const T& s2,
                             const T& l,
                             const T& z,
                             const T& k,
                             const T* uxs,
                             const T* uweights,
                             const size_t nUs)
  {

    const T carg = cos(-k*z);
    const T sarg = sin(-k*z);
    const T sm = T(0.5) * (s2 - s1);
    const T sp = T(0.5) * (s2 + s1);

    const T z2 = SQUARE(z);
    const T l2 = SQUARE(l);

    // Integral
    T intWreal = T(0.0), intWimag = T(0.0);

    for (size_t iu = 0 ; iu < nUs ; iu++) {

      T s = sm * uxs[iu] + sp;
      T s22 = SQUARE(s);

      T argw = -k * sqrt(s22 + z2 + l2);
      T real = uweights[iu] * (cos(argw) - carg) / std::max<T>((s22+l2),std::numeric_limits<T>::epsilon());
      T imag = uweights[iu] * (sin(argw) - sarg) / std::max<T>((s22+l2),std::numeric_limits<T>::epsilon());
      intWreal += real;
      intWimag += imag;
    }

    intWreal *= sm * l;
    intWimag *= sm * l;

    return std::complex<T>(intWreal,intWimag);
  }

  template <class T>
  STATIC_INLINE_BEGIN
  std::complex<T> CalcSingleFast(const T& s1,
                                 const T& s2,
                                 const T& l,
                                 const T& z,
                                 const T& k,
                                 const T* uxs,
                                 const T* uweights,
                                 const size_t nUs)
  {
    SPS_UNREFERENCED_PARAMETERS(s1, s2, l, z, k, uxs, uweights, nUs);
    return std::complex<T>();
  }

  template <>
  inline
  std::complex<float> CalcSingleFast(const float& s1,
                                     const float& s2,
                                     const float& l,
                                     const float& z,
                                     const float& k,
                                     const float* uxs, // Use restrict
                                     const float* uweights,
                                     const size_t nUs)
  {

    assert(((uintptr_t)&uxs[0] & 0x0F) == 0      && "Data must be aligned");
    assert(((uintptr_t)&uweights[0] & 0x0F) == 0 && "Data must be aligned");

    const __m128 v_s1 = _mm_set1_ps(s1);
    const __m128 v_s2 = _mm_set1_ps(s2);

    const __m128 sm = _mm_mul_ps(_m_half_ps,_mm_sub_ps(v_s2,v_s1));
    const __m128 sp = _mm_mul_ps(_m_half_ps,_mm_add_ps(v_s2,v_s1));

    __m128 carg = _mm_setzero_ps();
    __m128 sarg = _mm_setzero_ps();

    _mm_sin_cos_ps(_mm_set1_ps(-k*z), &sarg, &carg);

    const __m128 v_z2 = _mm_square_ps(_mm_set1_ps(z));
    const __m128 v_l2 = _mm_square_ps(_mm_set1_ps(l));

    __m128 intWreal = _mm_setzero_ps();
    __m128 intWimag = _mm_setzero_ps();

    __m128 cargw = _mm_setzero_ps();
    __m128 sargw = _mm_setzero_ps();

    for (size_t iu = 0 ; iu < nUs ; iu += 4) {

      __m128 s   = _mm_add_ps(
                     _mm_mul_ps(
                       sm,
                       _mm_load_ps(&uxs[iu])),
                     sp);
      __m128 s22 = _mm_square_ps(s);

      __m128 argw = _mm_mul_ps(
                      _mm_set1_ps(-k),
                      _mm_sqrt_ps(
                        _mm_add_ps(
                          v_l2,
                          _mm_add_ps(
                            v_z2,
                            s22))));

      _mm_sin_cos_ps(argw,&sargw,&cargw);

      __m128 denom = _mm_add_ps(s22,
                                v_l2);
      denom = _mm_rcp_ps(denom);
      __m128 uweight = _mm_load_ps((float*)&uweights[iu]);

      __m128 real = _mm_mul_ps(
                      _mm_mul_ps(
                        uweight,
                        _mm_sub_ps(
                          cargw,
                          carg)),
                      denom);

      __m128 imag = _mm_mul_ps(
                      _mm_mul_ps(
                        uweight,
                        _mm_sub_ps(
                          sargw,
                          sarg)),
                      denom);

      intWreal = _mm_add_ps(real,intWreal);
      intWimag = _mm_add_ps(imag,intWimag);
    }

    __m128 scale = _mm_mul_ps(sm,_mm_set1_ps(l));
    intWreal = _mm_mul_ps(scale,intWreal);
    intWimag = _mm_mul_ps(scale,intWimag);

    __m128 tmp = _mm_dp_ps(intWreal, _m_one_ps, 0xF1);

    float real = 0.0f;
    real = _mm_cvtss_f32(tmp);

    tmp = _mm_dp_ps(intWimag, _m_one_ps, 0xF1);
    float imag = 0.0f;
    imag = _mm_cvtss_f32(tmp);

    std::complex<float> c;
    c.real(real);
    c.imag(imag);
    return c;
  }

  template <class T>
  std::complex<T> CalcHz(const T& s,
                         const T& l,
                         const T& z,
                         const T& k,
                         const GLQuad2D<T>* gl)
  {

    const T carg = cos(-k*z);
    const T sarg = sin(-k*z);
    const T l_2 = T(0.5) * l;
    const T s_2 = T(0.5) * s;

    const T z2 = SQUARE(z);
    const T l2 = SQUARE(l);
    const T s2 = SQUARE(s);

    // integral width
    std::complex<T> intW = std::complex<T>(T(0.0),T(0.0));

    for (size_t iu = 0 ; iu < gl->u.nx ; iu++) {

      // Is this right?
      T ls = l_2 * gl->u.xs[iu] + l_2;
      T ls2 = SQUARE(ls);

      T argw = -k * sqrt(ls2 + z2 + s2);
      T real = gl->u.ws[iu] * (cos(argw) - carg) / (ls2+s2);
      T imag = gl->u.ws[iu] * (sin(argw) - sarg) / (ls2+s2);
      intW.real(real + intW.real());
      intW.imag(imag + intW.imag());
    }
    intW *= l_2 * s / (T(M_2PI)*k);

    // mult by i
    T realW = intW.real();
    intW.real(intW.imag());
    intW.imag(-realW);

    // integral height
    std::complex<T> intH = std::complex<T>(T(0.0),T(0.0));

    for(size_t iv = 0 ; iv < gl->v.nx ; iv++) {
      T ss = s_2 * gl->u.xs[iv] + s_2;
      T ss2 = SQUARE(ss);
      T argh = -k * sqrt(ss2 + z2 + l2);
      T real = gl->v.ws[iv] * (cos(argh) - carg) / (ss2+l2);
      T imag = gl->v.ws[iv] * (sin(argh) - sarg) / (ss2+l2);
      intH.real(real + intH.real());
      intH.imag(imag + intH.imag());
    }
    intH *= s_2 * l / (T(M_2PI)*k);
    T realH = intH.real();
    intH.real(intH.imag());
    intH.imag(-realH);

    intH = intH + intW;
    return intH;
  }

  template <class T>
  void CalcCwField(const sysparm_t<T>* sysparm,
                   const ApertureData<T>& data,
                   const T* pos, const size_t nPositions,
                   std::complex<T>** odata)
  {

    const T lambda = sysparm->c / data.m_f0;

    const T k = T(M_2PI)/lambda;

    const size_t nDivW = sysparm->nDivW;
    const size_t nDivH = sysparm->nDivH;

    auto uweights = sps::deleted_aligned_array<T>();
    auto vweights = sps::deleted_aligned_array<T>();
    auto uxs      = sps::deleted_aligned_array<T>();
    auto vxs      = sps::deleted_aligned_array<T>();

    CalcWeightsAndAbcissae(sysparm, std::move(uxs), std::move(uweights),
                           std::move(vxs),std::move(vweights));

    GLQuad2D<T> uv;
    uv.u.xs = uxs.get();
    uv.u.ws = uweights.get();
    uv.u.nx = nDivW;

    uv.v.xs = vxs.get();
    uv.v.ws = vweights.get();
    uv.v.nx = nDivH;

    sps::point_t<T> point;

    const size_t nElements = data.m_nelements;

    const size_t nSubElements = data.m_nsubelements;

    // const auto& elements = data.m_elements;
    size_t _nElements, _nSubElements;
    const sps::element_rect_t<T>** elements = NULL;
    data.ElementsRefGet(&_nElements, &_nSubElements, elements);

    for (size_t iPoint = 0 ; iPoint < nPositions ; iPoint++) {

      point[0] = pos[iPoint*3 + 0];
      point[1] = pos[iPoint*3 + 1];
      point[2] = pos[iPoint*3 + 2];

      std::complex<T> field = std::complex<T>(T(0.0),T(0.0));

      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {

        T apodization = data.m_apodizations[iElement];
        if (apodization != 0.0) {
          std::complex<T> field1 = std::complex<T>(T(0.0),T(0.0));

          for (size_t jElement = 0 ; jElement < nSubElements ; jElement++) {

            const sps::element_rect_t<T>& element = elements[iElement][jElement];

            // Get basis vectors (can be stored)
            sps::point_t<T> hh_dir, hw_dir, normal;

            sps::basis_vectors<T, sps::EulerIntrinsicYXY>(hw_dir, element.euler, 0);
            sps::basis_vectors<T, sps::EulerIntrinsicYXY>(hh_dir, element.euler, 1);
            sps::basis_vectors<T, sps::EulerIntrinsicYXY>(normal, element.euler, 2);

            sps::point_t<T> r2p = point - element.center;

            // Distance to plane
            T dist2plane = fabs(dot(normal,r2p));

            // Projection onto plane
            T u = dot(hw_dir,r2p);
            T v = dot(hh_dir,r2p);

            T z  = dist2plane;

            T l = fabs(u) + element.hw;
            T s = fabs(v) + element.hh;

            field1 += CalcHz<T>(fabs(s), fabs(l), z, k,
                                &uv)*T(signum<T>(s)*signum<T>(l));

            l = element.hw - fabs(u);

            field1 += CalcHz<T>(fabs(s),fabs(l),z,k,
                                &uv)*T(signum<T>(s)*signum<T>(l));

            l = fabs(u) + element.hw;
            s = element.hh - fabs(v);

            field1 += CalcHz<T>(fabs(s),fabs(l),z,k,
                                &uv)*T(signum<T>(s)*signum<T>(l));

            l = element.hw - fabs(u);

            field1 += CalcHz<T>(fabs(s),fabs(l),z,k,
                                &uv)*T(signum<T>(s)*signum<T>(l));
          }

          T real = apodization * field1.real();
          T imag = apodization * field1.imag();

          // Using common phase
          T arg = data.m_phases[iElement*nSubElements];
          T carg = cos(arg);
          T sarg = sin(arg);
          field.real(field.real() + real*carg - imag*sarg);
          field.imag(field.imag() + real*sarg + imag*carg);
        }
      }
      (*odata)[iPoint] = field;
    }
  }



  template <class T>
  int CalcCwThreaded(const fnm::sysparm_t<T>* sysparm,
                     const ApertureData<T>* data,
                     const T* pos, const size_t nPositions,
                     std::complex<T>** odata,
                     sps::ProgressBarInterface* pBar)
  {
    SPS_UNREFERENCED_PARAMETER(pBar);

    ProfilerStart();

    int retval = 0;

    const T lambda = sysparm->c / data->m_f0;
    const T k = T(M_2PI)/lambda;

    size_t nDivW = sysparm->nDivW;
    size_t nDivH = sysparm->nDivH;

    auto uweights = sps::deleted_aligned_array<T>();
    auto vweights = sps::deleted_aligned_array<T>();
    auto uxs      = sps::deleted_aligned_array<T>();
    auto vxs      = sps::deleted_aligned_array<T>();

    CalcWeightsAndAbcissae(sysparm,std::move(uxs),std::move(uweights),
                           std::move(vxs),std::move(vweights));

#ifndef HAVE_THREAD
    fnm::CwFieldThreadArg<T> threadarg;
    threadarg.data        = data;
    threadarg.iPointBegin = 0;
    threadarg.iPointEnd   = nPositions;
    threadarg.pos         = pos;
    threadarg.field       = *odata;
    threadarg.k           = k;
    threadarg.uxs         = uxs.get();
    threadarg.uweights    = uweights.get();
    threadarg.nDivU       = nDivW; // Works for threaded
    threadarg.vxs         = vxs.get();
    threadarg.vweights    = vweights.get();
    threadarg.nDivV       = nDivH;
    threadarg.threadId    = 0;
    threadarg.cpu_id      = 0;

# ifdef _WIN32
    retval                = (unsigned int)CalcCwThreadFunc((void*)&threadarg);
# else
    void* thread_retval   = CalcCwThreadFunc((void*)&threadarg);
    SPS_UNREFERENCED_PARAMETER(thread_retval);
# endif
#else

    int nproc = getncpus();

# if defined(_WIN32)
    unsigned int threadID;
    uintptr_t threads[N_MAX_THREADS];
# endif

    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

# ifdef HAVE_PTHREAD_H
    CallErr(pthread_attr_init,(&attr));
    CallErr(pthread_attr_setdetachstate, (&attr, PTHREAD_CREATE_JOINABLE));
# endif

    debug_print("nthreads: %zu\n", Aperture<T>::nthreads);

    sps::queue<float>* progressQueue = new sps::queue<float>();

    CwFieldThreadArg<T> threadarg[N_MAX_THREADS];

    // Populate structs for threads
    for (size_t i=0 ; i < Aperture<T>::nthreads ; i++) {
      threadarg[i].sysparm     = sysparm;
      threadarg[i].data        = data;
      threadarg[i].iPointBegin = 0+i*(nPositions/Aperture<T>::nthreads);
      threadarg[i].iPointEnd   = (nPositions/Aperture<T>::nthreads)+i*(nPositions/Aperture<T>::nthreads);
      threadarg[i].pos         = pos;
      threadarg[i].field       = (*odata);
      threadarg[i].k           = k;
      threadarg[i].uxs         = uxs.get();
      threadarg[i].uweights    = uweights.get();
      threadarg[i].nDivU       = nDivW;
      threadarg[i].vxs         = vxs.get();
      threadarg[i].vweights    = vweights.get();
      threadarg[i].nDivV       = nDivH;
      threadarg[i].threadId   = i;
      threadarg[i].cpu_id      = ((int) i) % nproc;

      threadarg[i].nPositions    = nPositions / Aperture<T>::nthreads; /* Used for profiling only */

      if (i==0)
        threadarg[i].pQueue    = progressQueue;

      if (i==(Aperture<T>::nthreads-1))
        threadarg[i].iPointEnd = nPositions;
    }

    /* Without message queues (slower) */
    for (size_t i=0; i<Aperture<T>::nthreads; i++) {
# if defined(HAVE_PTHREAD_H)
      CallErr(pthread_create,
              (&fnm::threads[i],
               &fnm::attr,
               &CalcCwThreadFunc<T>,
               &threadarg[i]));
# elif defined(_WIN32)
      threads[i] =
        _beginthreadex(NULL, 0,
                       &CalcCwThreadFunc<T>,
                       &threadarg[i], 0, &threadID );
# endif
    }

    // Monitor progress
    float fProgress = 0.0f;

#ifdef _MSC_VER
# pragma warning(push)
# pragma warning(disable:4127)
#endif
    while(true) {
      fProgress = progressQueue->pop();
      if (fProgress == 100.0f) {
        break;
      }
#ifdef USE_PROGRESS_BAR
      if (pBar) {
        pBar->show(fProgress); // Causes Segfault (consider flushing)
      }
#endif
#ifdef __GNUC__
      sched_yield();
#elif defined(_WIN32)
      SwitchToThread();
#endif
    }
#ifdef _MSC_VER
# pragma warning(pop)
#endif

    for (size_t i = 0; i < Aperture<T>::nthreads; i++) {
#  if defined(HAVE_PTHREAD_H)
      CallErr(pthread_join,(threads[i],NULL));
#  elif defined(_WIN32)
      WaitForSingleObject((HANDLE) threads[i], INFINITE );
#  endif
    }

    delete progressQueue;

    ProfilerStop();

    // Without message queues we destroy attributes
#  if defined(HAVE_PTHREAD_H)
    CallErr(pthread_attr_destroy,(&attr));
#  endif
#endif
    return retval;
  }

  template <class T>
  int CalcCwFieldFourRef(const sysparm_t<T>* sysparm,
                         const ApertureData<T>* data,
                         const T* pos, const size_t nPositions,
                         std::complex<T>** odata)
  {
    const T lambda = sysparm->c / data->m_f0;

    const T k = T(M_2PI) / lambda;

    size_t nDivW = sysparm->nDivW;
    size_t nDivH = sysparm->nDivH;

    auto uweights = sps::deleted_aligned_array<T>();
    auto vweights = sps::deleted_aligned_array<T>();
    auto uxs      = sps::deleted_aligned_array<T>();
    auto vxs      = sps::deleted_aligned_array<T>();

    CalcWeightsAndAbcissae(sysparm, std::move(uxs), std::move(uweights),
                           std::move(vxs), std::move(vweights));

    sps::point_t<T> point;

    const size_t nElements = data->m_nelements;

    const size_t nSubElements = data->m_nsubelements;

    size_t _nElements, _nSubElements;
    const sps::element_rect_t<T>** elements = NULL;
    data->ElementsRefGet(&_nElements, &_nSubElements, elements);

    const T* apodizations = data->m_apodizations.get();

    T apodization;

    for (size_t iPoint = 0 ; iPoint < nPositions ; iPoint++) {

      point[0] = pos[iPoint*3 + 0];
      point[1] = pos[iPoint*3 + 1];
      point[2] = pos[iPoint*3 + 2];

      std::complex<T> field = std::complex<T>(T(0.0),T(0.0));

      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {

        apodization = apodizations[iElement];

        // Output of individual integrals are complex
        std::complex<T> field1 = std::complex<T>(T(0.0),T(0.0));

        if (apodization != 0.0) {

          for (size_t jElement = 0 ; jElement < data->m_nsubelements ; jElement++) {

            const sps::element_rect_t<T>& element = elements[iElement][jElement];

            // Get basis vectors (can be stored - are they initialized)
            sps::point_t<T> hh_dir = sps::point_t<T>();
            sps::point_t<T> hw_dir = sps::point_t<T>();
            sps::point_t<T> normal = sps::point_t<T>();

            sps::basis_vectors<T, sps::EulerIntrinsicYXY>(hw_dir, element.euler, 0);
            sps::basis_vectors<T, sps::EulerIntrinsicYXY>(hh_dir, element.euler, 1);
            sps::basis_vectors<T, sps::EulerIntrinsicYXY>(normal, element.euler, 2);

            debug_print("hw_dir: %f %f %f\n", hw_dir[0], hw_dir[1], hw_dir[2]);
            debug_print("hh_dir: %f %f %f\n", hh_dir[0], hh_dir[1], hh_dir[2]);
            debug_print("normal: %f %f %f\n", normal[0], normal[1], normal[2]);
            sps::point_t<T> r2p = point - element.center;

            // Distance to plane
            T dist2plane = dot(normal,r2p);

            // Projection onto plane
            T u = dot(hw_dir,r2p);
            T v = dot(hh_dir,r2p);
            T z  = dist2plane;
            debug_print("u: %f, v: %f, z: %f\n",u,v,z);

            T s = fabs(v) + element.hh;
            std::complex<T> tmp = std::complex<T>(T(0.0),T(0.0));

            // u-integral  x (Python), hw is a (Python)
            if (fabs(u) > element.hw) {
              debug_print("outside\n");
              // Outside
              tmp = CALC_SELECT<T>(fabs(u)-element.hw, fabs(u)+element.hw, s, z, k, uxs.get(),uweights.get(),nDivW);
              field1 += tmp;
              s = element.hh - fabs(v);
              tmp = CALC_SELECT<T>(fabs(u)-element.hw, fabs(u)+element.hw, s, z, k, uxs.get(),uweights.get(),nDivW);
              field1 += tmp;
            } else {

              debug_print("inside\n");
              debug_print("low: %f, high: %f, l: %f, z: %f, k: %f\n",
                          T(0.0), fabs(u) + element.hw, s, z, k);
              tmp = CALC_SELECT<T>(fabs(u)-element.hw, element.hw+fabs(u), s, z, k, uxs.get(),uweights.get(),nDivW);
              field1 += tmp;

              s = element.hh - fabs(v);
              tmp = CALC_SELECT<T>(fabs(u)-element.hw, element.hw+fabs(u), s, z, k, uxs.get(),uweights.get(),nDivW);
              field1 += tmp;
            }

            s = fabs(u) + element.hw;
            if (fabs(v) > element.hh) {
              debug_print("outside\n");
              // Outside
              tmp = CALC_SELECT<T>(fabs(v)-element.hh, fabs(v)+element.hh, s, z, k, vxs.get(),vweights.get(),nDivH);
              field1 += tmp;
              s = element.hw - fabs(u);
              tmp = CALC_SELECT<T>(fabs(v)-element.hh, fabs(v)+element.hh, s, z, k, vxs.get(),vweights.get(),nDivH);
              field1 += tmp;
            } else {
              // Inside (changes result if CalcSingleFast)
              s = fabs(u) + element.hw;
              debug_print("inside\n");
              tmp = CALC_SELECT<T>(fabs(v)-element.hh, fabs(v)+element.hh, s, z, k, vxs.get(),vweights.get(),nDivH);
              debug_print("int_v: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
              s = element.hw - fabs(u);
              tmp = CALC_SELECT<T>(fabs(v)-element.hh, fabs(v)+element.hh, s, z, k, vxs.get(),vweights.get(),nDivH);
              debug_print("int_v: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
            }
          }

          T real = field1.real();
          T imag = field1.imag();

          // Multiply with i
          std::swap(real,imag);
          imag = -imag;

          // Divide with 2*pi*k
          real = real / (T(M_2PI)*k);
          imag = imag / (T(M_2PI)*k);

          // Phases (common, we ignore index jElement)
          T arg = data->m_phases[iElement*nSubElements];
          debug_print("arg: %f\n",arg);
          T carg = cos(arg);
          T sarg = sin(arg);
          field.real(field.real() + real*carg - imag*sarg);
          field.imag(field.imag() + real*sarg + imag*carg);
        } /* if (apodization != 0.0) */
      } /* for (size_t iElement = 0 ; iElement < nElements ; iElement++) */

      (*odata)[iPoint] = field;
    }
    return 0;
  }

  // TODO: FIX ME
#ifdef USE_PROGRESS_BAR
  template <class T>
  int CalcCwFieldRef(const sysparm_t<T>* sysparm,
                     const ApertureData<T>* data,
                     const T* pos, const size_t nPositions,
                     std::complex<T>** odata,
                     sps::ProgressBarInterface* pBar)
#else
  template <class T>
  int CalcCwFieldRef(const sysparm_t<T>* sysparm,
                     const ApertureData<T>* data,
                     const T* pos, const size_t nPositions,
                     std::complex<T>** odata,
                     void* pBar)
#endif
  {

    const T lambda = sysparm->c / data->m_f0;

    const T k = T(M_2PI) / lambda;

    size_t nDivW = sysparm->nDivW;
    size_t nDivH = sysparm->nDivH;

    auto uweights = sps::deleted_aligned_array<T>();
    auto vweights = sps::deleted_aligned_array<T>();
    auto uxs      = sps::deleted_aligned_array<T>();
    auto vxs      = sps::deleted_aligned_array<T>();

    CalcWeightsAndAbcissae(sysparm, std::move(uxs), std::move(uweights),
                           std::move(vxs), std::move(vweights));

    sps::point_t<T> point;

    const size_t nElements = data->m_nelements;

    const size_t nSubElements = data->m_nsubelements;

    size_t _nElements, _nSubElements;
    const sps::element_rect_t<T>** elements = NULL;
    data->ElementsRefGet(&_nElements, &_nSubElements, elements);

    const T* apodizations = data->m_apodizations.get();

    T apodization;

    // Initial time
    double tProfStart = sps::profiler::time();

    // Every n'th point, we issue a callback
    size_t nCallbackPeriod = 0;
    double duration = 0;

    for (size_t iPoint = 0 ; iPoint < nPositions ; iPoint++) {

      point[0] = pos[iPoint*3 + 0];
      point[1] = pos[iPoint*3 + 1];
      point[2] = pos[iPoint*3 + 2];

      std::complex<T> field = std::complex<T>(T(0.0),T(0.0));

      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {

        apodization = apodizations[iElement];

        // Output of individual integrals are complex
        std::complex<T> field1 = std::complex<T>(T(0.0),T(0.0));

        if (apodization != 0.0) {

          for (size_t jElement = 0 ; jElement < data->m_nsubelements ; jElement++) {

            const sps::element_rect_t<T>& element = elements[iElement][jElement];

            // Get basis vectors (can be stored - are they initialized)
            sps::point_t<T> hh_dir = sps::point_t<T>();
            sps::point_t<T> hw_dir = sps::point_t<T>();
            sps::point_t<T> normal = sps::point_t<T>();

            sps::basis_vectors<T, sps::EulerIntrinsicYXY>(hw_dir, element.euler, 0);
            sps::basis_vectors<T, sps::EulerIntrinsicYXY>(hh_dir, element.euler, 1);
            sps::basis_vectors<T, sps::EulerIntrinsicYXY>(normal, element.euler, 2);

            debug_print("hw_dir: %f %f %f\n", hw_dir[0], hw_dir[1], hw_dir[2]);
            debug_print("hh_dir: %f %f %f\n", hh_dir[0], hh_dir[1], hh_dir[2]);
            debug_print("normal: %f %f %f\n", normal[0], normal[1], normal[2]);
            sps::point_t<T> r2p = point - element.center;

            // Distance to plane
            T dist2plane = dot(normal,r2p);

            // Projection onto plane
            T u = dot(hw_dir,r2p);
            T v = dot(hh_dir,r2p);
            T z  = dist2plane;
            debug_print("u: %f, v: %f, z: %f\n",u,v,z);

            // We could wait multiplying by -i and dividing with (2*pi*k) till the end

            T s = fabs(v) + element.hh;
            std::complex<T> tmp = std::complex<T>(T(0.0),T(0.0));

            // u-integral  x (Python), hw is a (Python)
            if (fabs(u) > element.hw) {
              debug_print("outside\n");
              // Outside
              tmp = CALC_SELECT<T>(fabs(u),            fabs(u)+element.hw, s, z, k, uxs.get(),uweights.get(),nDivW);
              debug_print("int_u: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
              tmp = CALC_SELECT<T>(fabs(u)-element.hw, fabs(u), s, z, k, uxs.get(),uweights.get(),nDivW);
              debug_print("int_u: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
              s = element.hh - fabs(v);
              tmp = CALC_SELECT<T>(fabs(u),            fabs(u)+element.hw, s, z, k, uxs.get(),uweights.get(),nDivW);
              debug_print("int_u: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
              tmp = CALC_SELECT<T>(fabs(u)-element.hw, fabs(u), s, z, k, uxs.get(),uweights.get(),nDivW);
              debug_print("int_u: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
            } else {
              debug_print("inside\n");
              debug_print("low: %f, high: %f, l: %f, z: %f, k: %f\n",
                          T(0.0), fabs(u) + element.hw, s, z, k);
              tmp = CALC_SELECT<T>(0,                  fabs(u)+element.hw, s, z, k, uxs.get(),uweights.get(),nDivW);
              debug_print("int_u: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
              tmp = CALC_SELECT<T>(0,                  element.hw-fabs(u), s, z, k, uxs.get(),uweights.get(),nDivW);
              debug_print("int_u: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
              s = element.hh - fabs(v);
              tmp = CALC_SELECT<T>(0,                  fabs(u)+element.hw, s, z, k, uxs.get(),uweights.get(),nDivW);
              debug_print("int_u: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
              tmp = CALC_SELECT<T>(0,                  element.hw-fabs(u), s, z, k, uxs.get(),uweights.get(),nDivW);
              debug_print("int_u: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
            }

            // v-integral
            s = fabs(u) + element.hw;
            if (fabs(v) > element.hh) {
              debug_print("outside\n");
              // Outside
              tmp = CALC_SELECT<T>(fabs(v)-element.hh, fabs(v), s, z, k, vxs.get(),vweights.get(),nDivH);
              field1 += tmp;
              debug_print("int_v: %f %f\n",tmp.real(),tmp.imag());
              tmp = CALC_SELECT<T>(fabs(v),            fabs(v)+element.hh, s, z, k, vxs.get(),vweights.get(),nDivH);
              field1 += tmp;
              debug_print("int_v: %f %f\n",tmp.real(),tmp.imag());
              s = element.hw - fabs(u);
              tmp = CALC_SELECT<T>(fabs(v)-element.hh, fabs(v), s, z, k, vxs.get(),vweights.get(),nDivH);
              debug_print("int_v: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
              tmp = CALC_SELECT<T>(fabs(v),            fabs(v)+element.hh, s, z, k, vxs.get(),vweights.get(),nDivH);
              debug_print("int_v: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
            } else {
              // Inside (changes result if CalcSingleFast)
              debug_print("inside\n");
              tmp = CALC_SELECT<T>(0,                  fabs(v)+element.hh, s, z, k, vxs.get(),vweights.get(),nDivH);
              debug_print("int_v: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
              s = element.hw - fabs(u);
              tmp = CALC_SELECT<T>(0,                  fabs(v)+element.hh, s, z, k, vxs.get(),vweights.get(),nDivH);
              debug_print("int_v: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
              s = fabs(u) + element.hw;
              tmp = CALC_SELECT<T>(0,                  element.hh-fabs(v), s, z, k, vxs.get(),vweights.get(),nDivH);
              debug_print("int_v: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
              s = element.hw - fabs(u);
              tmp = CALC_SELECT<T>(0,                  element.hh-fabs(v), s, z, k, vxs.get(),vweights.get(),nDivH);
              debug_print("int_v: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
            }
          }

          T real = field1.real();
          T imag = field1.imag();

          // Multiply with i
          std::swap(real,imag);
          imag = -imag;

          // Divide with 2*pi*k
          real = real / (T(M_2PI)*k);
          imag = imag / (T(M_2PI)*k);

          // Phases (common, we ignore index jElement)
          T arg = data->m_phases[iElement*nSubElements];
          debug_print("arg: %f\n",arg);
          T carg = cos(arg);
          T sarg = sin(arg);
          field.real(field.real() + real*carg - imag*sarg);
          field.imag(field.imag() + real*sarg + imag*carg);
        } /* if (apodization != 0.0) */
      } /* for (size_t iElement = 0 ; iElement < nElements ; iElement++) */

      (*odata)[iPoint] = field;

      // Check time
      if (iPoint == 9) {
        duration = sps::profiler::time() - tProfStart;
        nCallbackPeriod = (size_t) (10.0 / duration);
        nCallbackPeriod = std::max<size_t>(1, nCallbackPeriod);
      }

      if ( (iPoint > 9) && (iPoint % nCallbackPeriod == 0) && (pBar)) {
        float val = 100.0f * float(iPoint) / std::max<size_t>(1,nPositions);

        SPS_UNREFERENCED_PARAMETER(val);
#if USE_PROGRESS_BAR
        pBar->show(val);
#endif

#ifdef __GNUC__
        sched_yield();
#elif defined(_WIN32)
        SwitchToThread();
#endif
      }
    }
    return 0;
  }

  template <class T>
  int CalcCwFocusNaiveFast(const sysparm_t<T>* sysparm,
                           const ApertureData<T>& data,
                           const T* pos, const size_t nPositions,
                           std::complex<T>** odata)
  {

    const T lambda = sysparm->c / data.m_f0;

    const T k = T(M_2PI)/lambda;

    const size_t nElements = data.m_nelements;

    const size_t nSubElements = data.m_nsubelements;

    const T* apodizations = data.m_apodizations.get();

    T apodization;


    // Need weights and abcissa values
    size_t nDivW = sysparm->nDivW;
    size_t nDivH = sysparm->nDivH;

    auto uweights = sps::deleted_aligned_array<T>();
    auto vweights = sps::deleted_aligned_array<T>();
    auto uxs      = sps::deleted_aligned_array<T>();
    auto vxs      = sps::deleted_aligned_array<T>();

    CalcWeightsAndAbcissae(sysparm, std::move(uxs), std::move(uweights),
                           std::move(vxs), std::move(vweights));

    // sps::deleted_aligned_multi_array<sps::element_rect_t<T>,2>& elements = data.m_elements;
    size_t _nElements, _nSubElements;
    const sps::element_rect_t<T>** elements = NULL;
    data.ElementsRefGet(&_nElements, &_nSubElements, elements);

    sps::point_t<T> projection;

    for (size_t iPoint = 0 ; iPoint < nPositions ; iPoint++) {

#ifdef FNM_DOUBLE_SUPPORT
      sps::point_t<T> point;
      point[0] = pos[iPoint*3];
      point[1] = pos[iPoint*3 + 1];
      point[2] = pos[iPoint*3 + 1];
#else
      __m128 vec_point = _mm_set_ps(0.0f, float(pos[iPoint*3 + 2]),float(pos[iPoint*3 + 1]),float(pos[iPoint*3]));
#endif

      std::complex<T> final = std::complex<T>(T(0.0),T(0.0));

      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {

        apodization = apodizations[iElement];

        if (apodization != 0.0) {

          for (size_t jElement = 0 ; jElement < nSubElements ; jElement++) {

            const sps::element_rect_t<T>& element = elements[iElement][jElement];

            std::complex<T> result;

#ifdef FNM_DOUBLE_SUPPORT
            // Scalar implementation
            sps::point_t<T> r2p = point - element.center;
            sps::point_t<T> hh_dir,hw_dir,normal;

            sps::basis_vectors<T, sps::EulerIntrinsicYXY>(hw_dir, element.euler, 0);
            sps::basis_vectors<T, sps::EulerIntrinsicYXY>(hh_dir, element.euler, 1);
            sps::basis_vectors<T, sps::EulerIntrinsicYXY>(normal, element.euler, 2);

            // Projection onto plane
            projection[0] = dot(hw_dir,r2p);
            projection[1] = dot(hh_dir,r2p);
            projection[2] = dot(normal,r2p);
#else
            assert(((uintptr_t)&element.center[0] & 0x0F) == 0 && "Data must be aligned");
            __m128 vec_r2p = _mm_sub_ps(vec_point, _mm_load_ps((float*)&element.center[0]));
            _mm_store_ss((float*)&projection[0],_mm_dp_ps(_mm_load_ps(&element.uvector[0]), vec_r2p,0x71));
            _mm_store_ss((float*)&projection[1],_mm_dp_ps(_mm_load_ps(&element.vvector[0]), vec_r2p,0x71));
            _mm_store_ss((float*)&projection[2],_mm_fabs_ps(_mm_dp_ps(_mm_load_ps(&element.normal[0]), vec_r2p,0x71)));
#endif
            result = CalcHzAll<T>(element, projection, k,
                                  uxs.get(), uweights.get(), nDivW,
                                  vxs.get(), vweights.get(), nDivH);
            final = final + apodization *
                    result * exp(std::complex<T>(0,data.m_phases[iElement*nSubElements+jElement]));
          }
        }
      }
      (*odata)[iPoint] = final;
    }
    return 0;
  }

  template <class T>
  int CalcCwFocusRef(const sysparm_t<T>* sysparm,
                     const ApertureData<T>& data,
                     const T* pos, const size_t nPositions,
                     std::complex<T>** odata)
  {
    const T lambda = sysparm->c / data.m_f0;

    const T k = T(M_2PI)/lambda;

    const size_t nElements    = data.m_nelements;

    const size_t nSubElements = data.m_nsubelements;

    const T* apodizations     = data.m_apodizations.get();

    T apodization;

    size_t nDivW = sysparm->nDivW;
    size_t nDivH = sysparm->nDivH;

    auto uweights = sps::deleted_aligned_array<T>();
    auto vweights = sps::deleted_aligned_array<T>();
    auto uxs      = sps::deleted_aligned_array<T>();
    auto vxs      = sps::deleted_aligned_array<T>();

    CalcWeightsAndAbcissae(sysparm,std::move(uxs),std::move(uweights),
                           std::move(vxs),std::move(vweights));

    //const auto& elements = data.m_elements;
    size_t _nElements, _nSubElements;
    const sps::element_rect_t<T>** elements = NULL;
    data.ElementsRefGet(&_nElements, &_nSubElements, elements);

    for (size_t iPoint = 0 ; iPoint < nPositions ; iPoint++) {

      sps::point_t<T> point;
      memcpy(&point[0],&pos[iPoint*3],3*sizeof(T));

      size_t iElement = iPoint % nElements;

      apodization = apodizations[iElement];

      std::complex<T> final = std::complex<T>(T(0.0),T(0.0));

      if (apodization != 0.0) {

        // Average over sub-elements
        for (size_t jElement = 0 ; jElement < nSubElements ; jElement++) {

          const sps::element_rect_t<T>& element = elements[iElement][jElement];

          // Get basis vectors (can be stored)
          sps::point_t<T> hh_dir, hw_dir, normal;

          sps::basis_vectors<T, sps::EulerIntrinsicYXY>(hw_dir, element.euler, 0);
          sps::basis_vectors<T, sps::EulerIntrinsicYXY>(hh_dir, element.euler, 1);
          sps::basis_vectors<T, sps::EulerIntrinsicYXY>(normal, element.euler, 2);

          sps::point_t<T> r2p = point - element.center;

          // Distance to plane
          T dist2plane = dot(normal,r2p);

          // Projection onto plane
          T u = dot(hw_dir,r2p);
          T v = dot(hh_dir,r2p);
          T z = dist2plane;

          T s = fabs(v) + element.hh;
          std::complex<T> tmp;
          // u-integral  x (Python), hw is a (Python)
          if (fabs(u) > element.hw) {
            // Outside
            tmp = CalcSingle<T>(fabs(u),            fabs(u)+element.hw, s, z, k, uxs.get(),uweights.get(),nDivW);
            debug_print("int_u: %f %f\n",tmp.real(),tmp.imag());
            final += tmp;
            tmp = CalcSingle<T>(fabs(u)-element.hw, fabs(u), s, z, k, uxs.get(),uweights.get(),nDivW);
            debug_print("int_u: %f %f\n",tmp.real(),tmp.imag());
            final += tmp;
            s = element.hh - fabs(v);
            tmp = CalcSingle<T>(fabs(u),            fabs(u)+element.hw, s, z, k, uxs.get(),uweights.get(),nDivW);
            debug_print("int_u: %f %f\n",tmp.real(),tmp.imag());
            final += tmp;
            tmp = CalcSingle<T>(fabs(u)-element.hw, fabs(u), s, z, k, uxs.get(),uweights.get(),nDivW);
            debug_print("int_u: %f %f\n",tmp.real(),tmp.imag());
            final += tmp;
          } else {
            debug_print("inside\n");
            // Inside
            debug_print("low: %f, high: %f, l: %f, z: %f, k: %f\n",
                        T(0.0), fabs(u) + element.hw, s, z, k);
            tmp = CalcSingle<T>(0,                  fabs(u)+element.hw, s, z, k, uxs.get(),uweights.get(),nDivW);
            debug_print("int_u: %f %f\n",tmp.real(),tmp.imag());
            final += tmp;
            tmp = CalcSingle<T>(0,                  element.hw-fabs(u), s, z, k, uxs.get(),uweights.get(),nDivW);
            debug_print("int_u: %f %f\n",tmp.real(),tmp.imag());
            final += tmp;
            s = element.hh - fabs(v);
            tmp = CalcSingle<T>(0,                  fabs(u)+element.hw, s, z, k, uxs.get(),uweights.get(),nDivW);
            debug_print("int_u: %f %f\n",tmp.real(),tmp.imag());
            final += tmp;
            tmp = CalcSingle<T>(0,                  element.hw-fabs(u), s, z, k, uxs.get(),uweights.get(),nDivW);
            debug_print("int_u: %f %f\n",tmp.real(),tmp.imag());
            final += tmp;
          }

          // v-integral
          s = fabs(u) + element.hw;
          if (fabs(v) > element.hh) {
            debug_print("outside\n");
            // Outside
            tmp = CalcSingle<T>(fabs(v)-element.hh, fabs(v), s, z, k, vxs.get(),vweights.get(),nDivH);
            final += tmp;
            debug_print("int_v: %f %f\n",tmp.real(),tmp.imag());
            tmp = CalcSingle<T>(fabs(v),            fabs(v)+element.hh, s, z, k, vxs.get(),vweights.get(),nDivH);
            final += tmp;
            debug_print("int_v: %f %f\n",tmp.real(),tmp.imag());
            s = element.hw - fabs(u);
            tmp = CalcSingle<T>(fabs(v)-element.hh, fabs(v), s, z, k, vxs.get(),vweights.get(),nDivH);
            debug_print("int_v: %f %f\n",tmp.real(),tmp.imag());
            final += tmp;
            tmp = CalcSingle<T>(fabs(v),            fabs(v)+element.hh, s, z, k, vxs.get(),vweights.get(),nDivH);
            debug_print("int_v: %f %f\n",tmp.real(),tmp.imag());
            final += tmp;
          } else {
            // Inside
            debug_print("inside\n");
            tmp = CalcSingle<T>(0,                  fabs(v)+element.hh, s, z, k, vxs.get(),vweights.get(),nDivH);
            debug_print("int_v: %f %f\n",tmp.real(),tmp.imag());
            final += tmp;
            s = element.hw - fabs(u);
            tmp = CalcSingle<T>(0,                  fabs(v)+element.hh, s, z, k, vxs.get(),vweights.get(),nDivH);
            debug_print("int_v: %f %f\n",tmp.real(),tmp.imag());
            final += tmp;
            s = fabs(u) + element.hw;
            tmp = CalcSingle<T>(0,                  element.hh-fabs(v), s, z, k, vxs.get(),vweights.get(),nDivH);
            debug_print("int_v: %f %f\n",tmp.real(),tmp.imag());
            final += tmp;
            s = element.hw - fabs(u);
            tmp = CalcSingle<T>(0,                  element.hh-fabs(v), s, z, k, vxs.get(),vweights.get(),nDivH);
            debug_print("int_v: %f %f\n",tmp.real(),tmp.imag());
            final += tmp;
          }
        }

        T real = final.real();
        T imag = final.imag();

        // Multiply with i
        std::swap(real,imag);
        imag = -imag;

        // Divide with 2*pi*k
        real = real / (T(M_2PI)*k);
        imag = imag / (T(M_2PI)*k);

        final.real(final.real() + real);
        final.imag(final.imag() + imag);
      } /* if (apodization != 0.0) */
      (*odata)[iPoint] = final;
    } /* for (size_t iPoint = 0 ; iPoint < nPositions ; iPoint++) */
    return 0;
  }

  template <class T>
  int CalcCwFocus(const sysparm_t<T>* sysparm,
                  const ApertureData<T>& data,
                  const T* pos, const size_t nPositions,
                  std::complex<T>** odata)
  {
    const T lambda = sysparm->c / data.m_f0;

    const T k = T(M_2PI)/lambda;

    const size_t nElements    = data.m_nelements;

    const size_t nSubElements = data.m_nsubelements;

    const T* apodizations     = data.m_apodizations.get();

    T apodization;

    size_t nDivW = sysparm->nDivW;
    size_t nDivH = sysparm->nDivH;

    auto uweights = sps::deleted_aligned_array<T>();
    auto vweights = sps::deleted_aligned_array<T>();
    auto uxs      = sps::deleted_aligned_array<T>();
    auto vxs      = sps::deleted_aligned_array<T>();

    CalcWeightsAndAbcissae(sysparm, std::move(uxs),std::move(uweights),
                           std::move(vxs),std::move(vweights));

    size_t _nElements, _nSubElements;
    const sps::element_rect_t<T>** elements = NULL;
    data.ElementsRefGet(&_nElements, &_nSubElements, elements);

    sps::point_t<T> projection;

    for (size_t iPoint = 0 ; iPoint < nPositions ; iPoint++) {

#ifdef FNM_DOUBLE_SUPPORT
      sps::point_t<T> point;
      memcpy(&point[0],&pos[iPoint*3],3*sizeof(T));
#else
      __m128 vec_point = _mm_set_ps(0.0f, float(pos[iPoint*3 + 2]),float(pos[iPoint*3 + 1]),float(pos[iPoint*3]));
#endif

      size_t iElement;
#if FNM_PHASED_FOCUS
      iElement = iPoint % (nElements*nSubElements);
#else
      iElement = iPoint % nElements;
#endif

      apodization = apodizations[iElement];

      std::complex<T> final = std::complex<T>(T(0.0),T(0.0));

      if (apodization != 0.0) {

        // Average over sub-elements
        for (size_t jElement = 0 ; jElement < nSubElements ; jElement++) {

          const sps::element_rect_t<T>& element = elements[iElement][jElement];

          std::complex<T> result;

#ifdef FNM_DOUBLE_SUPPORT
          // Scalar implementation
          sps::point_t<T> r2p = point - element.center;
          sps::point_t<T> hh_dir,hw_dir,normal;

          sps::basis_vectors<T, sps::EulerIntrinsicYXY>(hw_dir, element.euler, 0);
          sps::basis_vectors<T, sps::EulerIntrinsicYXY>(hh_dir, element.euler, 1);
          sps::basis_vectors<T, sps::EulerIntrinsicYXY>(normal, element.euler, 2);

          // Projection onto plane
          projection[0] = dot(hw_dir,r2p);
          projection[1] = dot(hh_dir,r2p);
          projection[2] = dot(normal,r2p);
#else
          // Vector implementation
          assert(((uintptr_t)&element.center[0] & 0x0F) == 0 && "Data must be aligned");
          assert(((uintptr_t)&element.uvector[0] & 0x0F) == 0 && "Data must be aligned");
          assert(((uintptr_t)&element.vvector[0] & 0x0F) == 0 && "Data must be aligned");

          __m128 vec_r2p = _mm_sub_ps(vec_point, _mm_load_ps((float*)&element.center[0]));

          _mm_store_ss((float*)&projection[0],_mm_dp_ps(_mm_load_ps((float*)&element.uvector[0]), vec_r2p,0x71));
          _mm_store_ss((float*)&projection[1],_mm_dp_ps(_mm_load_ps((float*)&element.vvector[0]), vec_r2p,0x71));
          _mm_store_ss((float*)&projection[2],_mm_dp_ps(_mm_load_ps((float*)&element.normal[0]), vec_r2p,0x71));

#endif

          result = CalcHzFast<T>(element, projection, k,
                                 uxs.get(), uweights.get(), nDivW,
                                 vxs.get(), vweights.get(), nDivH);

          T real = result.real();
          T imag = result.imag();

          final.real(final.real() + real);
          final.imag(final.imag() + imag);
        }
      } /* if (apodization != 0.0) */
      (*odata)[iPoint] = final;
    }
    return 0;
  }


  template std::complex<float>
  CalcSingle(const float& s1,
             const float& s2,
             const float& l,
             const float& z,
             const float& k,
             const float* uxs,
             const float* uweights,
             const size_t nUs);

  template std::complex<float>
  CalcSingleFast(const float& s1,
                 const float& s2,
                 const float& l,
                 const float& z,
                 const float& k,
                 const float* uxs,
                 const float* uweights,
                 const size_t nUs);

  template void FNM_EXPORT CalcCwField(const sysparm_t<float>* sysparm,
                                       const ApertureData<float>& data,
                                       const float* pos, const size_t nPositions,
                                       std::complex<float>** odata);

# ifdef USE_PROGRESS_BAR
  template int FNM_EXPORT CalcCwFieldRef(const sysparm_t<float>* sysparm,
                                         const ApertureData<float>* data,
                                         const float* pos, const size_t nPositions,
                                         std::complex<float>** odata,
                                         sps::ProgressBarInterface* pbar);
# else
  template int FNM_EXPORT CalcCwFieldRef(const sysparm_t<float>* sysparm,
                                         const ApertureData<float>* data,
                                         const float* pos, const size_t nPositions,
                                         std::complex<float>** odata,
                                         void* pbar);
# endif

  template int FNM_EXPORT CalcCwFieldFourRef(const sysparm_t<float>* sysparm,
      const ApertureData<float>* data,
      const float* pos, const size_t nPositions,
      std::complex<float>** odata);

  template int FNM_EXPORT CalcCwFocus(const sysparm_t<float>* sysparm,
                                      const ApertureData<float>& data,
                                      const float* pos, const size_t nPositions,
                                      std::complex<float>** odata);

  template int FNM_EXPORT CalcCwFocusRef(const sysparm_t<float>* sysparm,
                                         const ApertureData<float>& data,
                                         const float* pos, const size_t nPositions,
                                         std::complex<float>** odata);

  template int FNM_EXPORT CalcCwFocusNaiveFast(const sysparm_t<float>* sysparm,
      const ApertureData<float>& data,
      const float* pos, const size_t nPositions,
      std::complex<float>** odata);


  template std::complex<float> FNM_EXPORT CalcSingleFast(const float& s1,
      const float& s2,
      const float& l,
      const float& z,
      const float& k,
      const float* uxs,
      const float* uweights,
      const size_t nUs);

  template int FNM_EXPORT CalcCwThreaded(const fnm::sysparm_t<float>* sysparm,
                                         const ApertureData<float>* data,
                                         const float* pos, const size_t nPositions,
                                         std::complex<float>** odata,
                                         sps::ProgressBarInterface* pBar);

#ifdef FNM_DOUBLE_SUPPORT

  template std::complex<double>
  CalcSingle(const double& s1,
             const double& s2,
             const double& l,
             const double& z,
             const double& k,
             const double* uxs,
             const double* uweights,
             const size_t nUs);

  template std::complex<double>
  CalcSingleFast(const double& s1,
                 const double& s2,
                 const double& l,
                 const double& z,
                 const double& k,
                 const double* uxs,
                 const double* uweights,
                 const size_t nUs);

  template void FNM_EXPORT CalcCwField(const sysparm_t<double>* sysparm,
                                       const ApertureData<double>& data,
                                       const double* pos, const size_t nPositions,
                                       std::complex<double>** odata);

# ifdef USE_PROGRESS_BAR
  template int FNM_EXPORT CalcCwFieldRef(const sysparm_t<double>* sysparm,
                                         const ApertureData<double>* data,
                                         const double* pos, const size_t nPositions,
                                         std::complex<double>** odata,
                                         sps::ProgressBarInterface* pbar);
# else
  template int FNM_EXPORT CalcCwFieldRef(const sysparm_t<double>* sysparm,
                                         const ApertureData<double>* data,
                                         const double* pos, const size_t nPositions,
                                         std::complex<double>** odata,
                                         void* pbar);
# endif

  template int FNM_EXPORT CalcCwFieldFourRef(const sysparm_t<double>* sysparm,
      const ApertureData<double>* data,
      const double* pos, const size_t nPositions,
      std::complex<double>** odata);

  template int FNM_EXPORT CalcCwFocus(const sysparm_t<double>* sysparm,
                                      const ApertureData<double>& data,
                                      const double* pos, const size_t nPositions,
                                      std::complex<double>** odata);

  template int FNM_EXPORT CalcCwFocusRef(const sysparm_t<double>* sysparm,
                                         const ApertureData<double>& data,
                                         const double* pos, const size_t nPositions,
                                         std::complex<double>** odata);

  template int FNM_EXPORT CalcCwFocusNaiveFast(const sysparm_t<double>* sysparm,
      const ApertureData<double>& data,
      const double* pos, const size_t nPositions,
      std::complex<double>** odata);

  template int FNM_EXPORT CalcCwThreaded(const fnm::sysparm_t<double>* sysparm,
                                         const ApertureData<double>* data,
                                         const double* pos, const size_t nPositions,
                                         std::complex<double>** odata,
                                         sps::ProgressBarInterface* pBar);


#endif
}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */

