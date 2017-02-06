
#include <fnm/fnm.hpp>
#include <fnm/fnm_common.hpp>
#include <fnm/fnm_calc.hpp>
#include <fnm/FnmSIMD.hpp>

#include <gl/gl.hpp>

#include <sps/smath.hpp>
#include <sps/mm_malloc.h>
#include <sps/memory>
#include <sps/sps_threads.hpp>
#include <sps/debug.h>

#include <memory>

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
    size_t thread_id;
    int cpu_id;
  };

  template <class T>
  void CalcWeightsAndAbcissae(sps::deleted_aligned_array<T> &&uxs,
                              sps::deleted_aligned_array<T> &&uweights,
                              sps::deleted_aligned_array<T> &&vxs,
                              sps::deleted_aligned_array<T> &&vweights)
  {

    size_t nDivW = Aperture<T>::_sysparm.nDivW;
    size_t nDivH = Aperture<T>::_sysparm.nDivH;

    uxs      = sps::deleted_aligned_array_create<T>(nDivW);
    uweights = sps::deleted_aligned_array_create<T>(nDivW);
    vxs      = sps::deleted_aligned_array_create<T>(nDivH);
    vweights = sps::deleted_aligned_array_create<T>(nDivH);

    memset(uxs.get(),0,sizeof(T)*nDivW);
    memset(uweights.get(),0,sizeof(T)*nDivW);
    memset(vxs.get(),0,sizeof(T)*nDivH);
    memset(vweights.get(),0,sizeof(T)*nDivH);

    // Common weights and abcissa
    for (size_t i = 0 ; i < nDivW ; i++) {
      gl::GLNode qp = gl::GL(nDivW,i);
      uxs[i]        = T(qp.value);
      uweights[i]   = T(qp.weight);
    }

    for (size_t i = 0 ; i < nDivH ; i++) {
      gl::GLNode qp = gl::GL(nDivH, i);
      vxs[i]      = T(qp.value);
      vweights[i] = T(qp.weight);
    }
  }

#ifndef FNM_DOUBLE_SUPPORT
#if defined(HAVE_PTHREAD_H)
  template <class T>
  void* CalcCwThreadFunc(void* ptarg)
#else
  template <class T>
  unsigned int __stdcall CalcCwThreadFunc(void *ptarg)
#endif
  {
    fnm::CwFieldThreadArg<T>* pThreadArg = reinterpret_cast<fnm::CwFieldThreadArg<T>*>(ptarg);

    const ApertureData<T>* m_data = pThreadArg->data;
    const T k                     = pThreadArg->k;
    const T* uxs                  = pThreadArg->uxs;
    const T* vxs                  = pThreadArg->vxs;
    const T* uweights             = pThreadArg->uweights;
    const T* vweights             = pThreadArg->vweights;
    const T* pos                  = pThreadArg->pos;
    const size_t nDivW            = pThreadArg->nDivU;
    const size_t nDivH            = pThreadArg->nDivV;
    std::complex<T>* odata        = pThreadArg->field;

#ifdef HAVE_THREAD
    setcpuid(pThreadArg->cpu_id);
#endif

    const size_t nElements    = m_data->m_nelements;
    const size_t nSubElements = m_data->m_nsubelements;
    const auto& elements      = m_data->m_elements;

    const T* apodizations = m_data->m_apodizations.get();

    T apodization;

    ALIGN16_BEGIN sps::point_t<T> projection ALIGN16_END;

#ifdef _MSC_VER
    debug_print("iPointBegin: %Iu, iPointEnd: %Iu\n", pThreadArg->iPointBegin, pThreadArg->iPointEnd);
#else
    debug_print("iPointBegin: %zu, iPointEnd: %zu\n", pThreadArg->iPointBegin, pThreadArg->iPointEnd);
#endif

#if FNM_ENABLE_ATTENUATION
    //T alpha = Aperture<T>::_sysparm.att * T(100) / T(1e6) * m_data->m_f0 * T(0.1151);
    // SI units
    T alpha = Aperture<T>::_sysparm.att;
    if (!(Aperture<T>::_sysparm.use_att)) {
      alpha = T(0.0);
    }
#endif

    // vectors, pos, output
    for (size_t iPoint = pThreadArg->iPointBegin ; iPoint < pThreadArg->iPointEnd ; iPoint++) {

      __m128 vec_point =
        _mm_set_ps(0.0f,
                   float(pos[iPoint*3 + 2]),
                   float(pos[iPoint*3 + 1]),
                   float(pos[iPoint*3 + 0]));

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

            result = CalcHzFast<T>(element, projection, k,
                                   uxs, uweights, nDivW,
                                   vxs, vweights, nDivH);

            T real = apodization * result.real();
            T imag = apodization * result.imag();

            T arg = m_data->m_phases[iElement*nSubElements+jElement];
            T carg = cos(arg);
            T sarg = sin(arg);

            // TODO: Fix attenuation (HERE). It is not working
#if FNM_ENABLE_ATTENUATION
            T dist = (T) _mm_cvtss_f32(_mm_sqrt_ps(_mm_dp_ps(vec_r2p,vec_r2p,0x71)));

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
    }
    debug_print("Thread done\n");
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
#endif

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
  std::complex<T> CalcHz(const T& s,
                         const T& l,
                         const T& z,
                         const T& k,
                         const T* uxs,
                         const T* uweights,
                         const size_t nUs,
                         const T* vxs,
                         const T* vweights,
                         const size_t nVs)
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

    for (size_t iu = 0 ; iu < nUs ; iu++) {

      T ls = l_2 * uxs[iu] + l_2;
      T ls2 = SQUARE(ls);

      T argw = -k * sqrt(ls2 + z2 + s2);
      T real = uweights[iu] * (cos(argw) - carg) / (ls2+s2);
      T imag = uweights[iu] * (sin(argw) - sarg) / (ls2+s2);
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

    for(size_t iv = 0 ; iv < nVs ; iv++) {
      T ss = s_2 * vxs[iv] + s_2;
      T ss2 = SQUARE(ss);
      T argh = -k * sqrt(ss2 + z2 + l2);
      T real = vweights[iv] * (cos(argh) - carg) / (ss2+l2);
      T imag = vweights[iv] * (sin(argh) - sarg) / (ss2+l2);
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
  void CalcCwField(const ApertureData<T>& data,
                   const T* pos, const size_t nPositions,
                   std::complex<T>** odata)
  {

    const T lambda = Aperture<T>::_sysparm.c / data.m_f0;

    const T k = T(M_2PI)/lambda;

    const size_t nDivW = Aperture<T>::_sysparm.nDivW;
    const size_t nDivH = Aperture<T>::_sysparm.nDivH;

    auto uweights = sps::deleted_aligned_array<T>();
    auto vweights = sps::deleted_aligned_array<T>();
    auto uxs      = sps::deleted_aligned_array<T>();
    auto vxs      = sps::deleted_aligned_array<T>();

    CalcWeightsAndAbcissae(std::move(uxs),std::move(uweights),
                           std::move(vxs),std::move(vweights));

    sps::point_t<T> point;

    const size_t nElements = data.m_nelements;

    const size_t nSubElements = data.m_nsubelements;

    const auto& elements = data.m_elements;

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

            const sps::element_t<T>& element = elements[iElement][jElement];

            // Get basis vectors (can be stored)
            sps::point_t<T> hh_dir, hw_dir, normal;

            basis_vectors(hw_dir, element.euler, 0);
            basis_vectors(hh_dir, element.euler, 1);
            basis_vectors(normal, element.euler, 2);

            sps::point_t<T> r2p = point - element.center;

            // Distance to plane
            T dist2plane = fabs(dot(normal,r2p));

            // Projection onto plane
            T u = dot(hw_dir,r2p);
            T v = dot(hh_dir,r2p);

            T z  = dist2plane;

            T l = fabs(u) + element.hw;
            T s = fabs(v) + element.hh;

            field1 += CalcHz<T>(fabs(s),fabs(l),z,k,uxs.get(),uweights.get(),nDivW,
                                vxs.get(),vweights.get(),nDivH)*T(signum<T>(s)*signum<T>(l));

            l = element.hw - fabs(u);

            field1 += CalcHz<T>(fabs(s),fabs(l),z,k,uxs.get(),uweights.get(),nDivW,
                                vxs.get(),vweights.get(),nDivH)*T(signum<T>(s)*signum<T>(l));

            l = fabs(u) + element.hw;
            s = element.hh - fabs(v);

            field1 += CalcHz<T>(fabs(s),fabs(l),z,k,uxs.get(),uweights.get(),nDivW,
                                vxs.get(),vweights.get(),nDivH)*T(signum<T>(s)*signum<T>(l));

            l = element.hw - fabs(u);

            field1 += CalcHz<T>(fabs(s),fabs(l),z,k,uxs.get(),uweights.get(),nDivW,
                                vxs.get(),vweights.get(),nDivH)*T(signum<T>(s)*signum<T>(l));
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



  // Works outside class (without message queues)
  template <class T>
  int CalcCwThreaded(const ApertureData<T>* data,
                     const T* pos, const size_t nPositions,
                     std::complex<T>** odata)
  {

    int retval = 0;

    const T lambda = Aperture<T>::_sysparm.c / data->m_f0;
    const T k = T(M_2PI)/lambda;

    size_t nDivW = Aperture<T>::_sysparm.nDivW;
    size_t nDivH = Aperture<T>::_sysparm.nDivH;

    auto uweights = sps::deleted_aligned_array<T>();
    auto vweights = sps::deleted_aligned_array<T>();
    auto uxs      = sps::deleted_aligned_array<T>();
    auto vxs      = sps::deleted_aligned_array<T>();

    CalcWeightsAndAbcissae(std::move(uxs),std::move(uweights),
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
    threadarg.thread_id   = 0;
    threadarg.cpu_id      = 0;

# ifdef _WIN32
    retval                = (unsigned int)CalcCwThreaded((void*)&threadarg);
# else
    void* thread_retval   = CalcCwThreaded((void*)&threadarg);
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

    CwFieldThreadArg<T> threadarg[N_MAX_THREADS];

    // Populate structs for threads
    for (size_t i=0 ; i < Aperture<T>::nthreads ; i++) {
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
      threadarg[i].thread_id   = i;
      threadarg[i].cpu_id      = ((int) i) % nproc;
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

    for (size_t i = 0; i < Aperture<T>::nthreads; i++) {
#  if defined(HAVE_PTHREAD_H)
      CallErr(pthread_join,(threads[i],NULL));
#  elif defined(_WIN32)
      WaitForSingleObject((HANDLE) threads[i], INFINITE );
#  endif
    }

    // Without message queues we destroy attributes
#  if defined(HAVE_PTHREAD_H)
    CallErr(pthread_attr_destroy,(&attr));
#  endif
#endif
    return retval;
  }

  template <class T>
  int CalcCwFieldRef(const ApertureData<T>& data,
                     const T* pos, const size_t nPositions,
                     std::complex<T>** odata)
  {
    const T lambda = Aperture<T>::_sysparm.c / data.m_f0;

    const T k = T(M_2PI)/lambda;

    size_t nDivW = Aperture<T>::_sysparm.nDivW;
    size_t nDivH = Aperture<T>::_sysparm.nDivH;

    auto uweights = sps::deleted_aligned_array<T>();
    auto vweights = sps::deleted_aligned_array<T>();
    auto uxs      = sps::deleted_aligned_array<T>();
    auto vxs      = sps::deleted_aligned_array<T>();

    CalcWeightsAndAbcissae(std::move(uxs),std::move(uweights),
                           std::move(vxs),std::move(vweights));

    sps::point_t<T> point;

    const size_t nElements = data.m_nelements;

    const size_t nSubElements = data.m_nsubelements;

    const auto& elements = data.m_elements;

    const T* apodizations = data.m_apodizations.get();

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

          for (size_t jElement = 0 ; jElement < data.m_nsubelements ; jElement++) {

            const sps::element_t<T>& element = elements[iElement][jElement];

            // Get basis vectors (can be stored)
            sps::point_t<T> hh_dir, hw_dir, normal;

            basis_vectors(hw_dir, element.euler, 0);
            basis_vectors(hh_dir, element.euler, 1);
            basis_vectors(normal, element.euler, 2);

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
            std::complex<T> tmp;
            // u-integral  x (Python), hw is a (Python)
            if (fabs(u) > element.hw) {
              // Outside
              tmp = CalcSingle<T>(fabs(u),            fabs(u)+element.hw, s, z, k, uxs.get(),uweights.get(),nDivW);
              debug_print("int_u: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
              tmp = CalcSingle<T>(fabs(u)-element.hw, fabs(u)           , s, z, k, uxs.get(),uweights.get(),nDivW);
              debug_print("int_u: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
              s = element.hh - fabs(v);
              tmp = CalcSingle<T>(fabs(u),            fabs(u)+element.hw, s, z, k, uxs.get(),uweights.get(),nDivW);
              debug_print("int_u: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
              tmp = CalcSingle<T>(fabs(u)-element.hw, fabs(u)           , s, z, k, uxs.get(),uweights.get(),nDivW);
              debug_print("int_u: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
            } else {
              // Inside
              debug_print("low: %f, high: %f, l: %f, z: %f, k: %f\n",
                          T(0.0), fabs(u) + element.hw, s, z, k);
              tmp = CalcSingle<T>(0,                  fabs(u)+element.hw, s, z, k, uxs.get(),uweights.get(),nDivW);
              debug_print("int_u: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
              tmp = CalcSingle<T>(0,                  element.hw-fabs(u), s, z, k, uxs.get(),uweights.get(),nDivW);
              debug_print("int_u: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
              s = element.hh - fabs(v);
              tmp = CalcSingle<T>(0,                  fabs(u)+element.hw, s, z, k, uxs.get(),uweights.get(),nDivW);
              debug_print("int_u: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
              tmp = CalcSingle<T>(0,                  element.hw-fabs(u), s, z, k, uxs.get(),uweights.get(),nDivW);
              debug_print("int_u: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
            }

            // v-integral
            s = fabs(u) + element.hw;
            if (fabs(v) > element.hh) {
              // Outside
              tmp = CalcSingle<T>(fabs(v)-element.hh, fabs(v)           , s, z, k, vxs.get(),vweights.get(),nDivH);
              field1 += tmp;
              debug_print("int_v: %f %f\n",tmp.real(),tmp.imag());
              tmp = CalcSingle<T>(fabs(v),            fabs(v)+element.hh, s, z, k, vxs.get(),vweights.get(),nDivH);
              field1 += tmp;
              debug_print("int_v: %f %f\n",tmp.real(),tmp.imag());
              s = element.hw - fabs(u);
              tmp = CalcSingle<T>(fabs(v)-element.hh, fabs(v)           , s, z, k, vxs.get(),vweights.get(),nDivH);
              debug_print("int_v: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
              tmp = CalcSingle<T>(fabs(v),            fabs(v)+element.hh, s, z, k, vxs.get(),vweights.get(),nDivH);
              debug_print("int_v: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
            } else {
              // Inside
              tmp = CalcSingle<T>(0,                  fabs(v)+element.hh, s, z, k, vxs.get(),vweights.get(),nDivH);
              debug_print("int_v: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
              s = element.hw - fabs(u);
              tmp = CalcSingle<T>(0,                  fabs(v)+element.hh, s, z, k, vxs.get(),vweights.get(),nDivH);
              debug_print("int_v: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
              s = fabs(u) + element.hw;
              tmp = CalcSingle<T>(0,                  element.hh-fabs(v), s, z, k, vxs.get(),vweights.get(),nDivH);
              debug_print("int_v: %f %f\n",tmp.real(),tmp.imag());
              field1 += tmp;
              s = element.hw - fabs(u);
              tmp = CalcSingle<T>(0,                  element.hh-fabs(v), s, z, k, vxs.get(),vweights.get(),nDivH);
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

          // Phases (common)
          T arg = data.m_phases[iElement*nSubElements];
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

  template <class T>
  int CalcCwFocus(const ApertureData<T>& data,
                  const T* pos, const size_t nPositions,
                  std::complex<T>** odata)
  {

    const T lambda = Aperture<T>::_sysparm.c / data.m_f0;

    const T k = T(M_2PI)/lambda;

    const size_t nElements    = data.m_nelements;

    const size_t nSubElements = data.m_nsubelements;

    const T* apodizations     = data.m_apodizations.get();

    T apodization;

    // Need weights and abcissa values

    size_t nDivW = Aperture<T>::_sysparm.nDivW;
    size_t nDivH = Aperture<T>::_sysparm.nDivH;

    auto uweights = sps::deleted_aligned_array<T>();
    auto vweights = sps::deleted_aligned_array<T>();
    auto uxs      = sps::deleted_aligned_array<T>();
    auto vxs      = sps::deleted_aligned_array<T>();

    CalcWeightsAndAbcissae(std::move(uxs),std::move(uweights),
                           std::move(vxs),std::move(vweights));

    const auto& elements = data.m_elements;

    sps::point_t<T> projection;

    // TODO: Update to work for double

    // vectors, pos, output
    for (size_t iPoint = 0 ; iPoint < nPositions ; iPoint++) {

#ifdef FNM_DOUBLE_SUPPORT
      sps::point_t<T> point;
      memcpy(&point[0],&pos[iPoint*3],3*sizeof(T));
#else
      __m128 vec_point = _mm_set_ps(0.0f, float(pos[iPoint*3 + 2]),float(pos[iPoint*3 + 1]),float(pos[iPoint*3]));
#endif

      size_t iElement = iPoint % nElements;

      apodization = apodizations[iElement];

      std::complex<T> final = std::complex<T>(T(0.0),T(0.0));

      if (apodization != 0.0) {

        // Only the phase multiplication can be vectorized if unrolled
        // and if we have sub-elements!!!
        for (size_t jElement = 0 ; jElement < nSubElements ; jElement++) {

          const sps::element_t<T>& element = elements[iElement][jElement];

          std::complex<T> result;

#ifdef FNM_DOUBLE_SUPPORT
          // Scalar implementation
          sps::point_t<T> r2p = point - element.center;
          sps::point_t<T> hh_dir,hw_dir,normal;

          // Projection onto plane
          projection[0] = dot(hw_dir,r2p);
          projection[1] = dot(hh_dir,r2p);
          projection[2] = fabs(dot(normal,r2p));
#else
          // Vector implementation, TODO: Specialize for double
          assert(((uintptr_t)&element.center[0] & 0x0F) == 0 && "Data must be aligned");
          assert(((uintptr_t)&element.uvector[0] & 0x0F) == 0 && "Data must be aligned");
          assert(((uintptr_t)&element.vvector[0] & 0x0F) == 0 && "Data must be aligned");
          __m128 vec_r2p = _mm_sub_ps(vec_point, _mm_load_ps((float*)&element.center[0]));

          // Use vectors stored with elements (not working, if not set using ElementsSet)
          _mm_store_ss((float*)&projection[0],_mm_dp_ps(_mm_load_ps((float*)&element.uvector[0]), vec_r2p,0x71));
          _mm_store_ss((float*)&projection[1],_mm_dp_ps(_mm_load_ps((float*)&element.vvector[0]), vec_r2p,0x71));
          _mm_store_ss((float*)&projection[2],_mm_fabs_ps(_mm_dp_ps(_mm_load_ps((float*)&element.normal[0]), vec_r2p,0x71)));
#endif
          result = CalcHzFast<T>(element, projection, k,
                                 uxs.get(), uweights.get(), nDivW,
                                 vxs.get(), vweights.get(), nDivH);
          T real = result.real();
          T imag = result.imag();

          T arg = data.m_phases[iElement*nSubElements + jElement];
          T carg = cos(arg);
          T sarg = sin(arg);
          final.real(final.real() + real*carg - imag*sarg);
          final.imag(final.imag() + real*sarg + imag*carg);
        }
      }
      (*odata)[iPoint] = final;
    }
    return 0;
  }


  template int FNM_EXPORT CalcCwFieldRef(const ApertureData<float>& data,
                                         const float* pos, const size_t nPositions,
                                         std::complex<float>** odata);

  template void FNM_EXPORT CalcCwField(const ApertureData<float>& data,
                                       const float* pos, const size_t nPositions,
                                       std::complex<float>** odata);

  template int FNM_EXPORT CalcCwThreaded(const ApertureData<float>* data,
                                         const float* pos, const size_t nPositions,
                                         std::complex<float>** odata);

  template int FNM_EXPORT CalcCwFocus(const ApertureData<float>& data,
                                      const float* pos, const size_t nPositions,
                                      std::complex<float>** odata);


#ifdef FNM_DOUBLE_SUPPORT
  template int FNM_EXPORT CalcCwFieldRef(const ApertureData<double>& data,
                                         const double* pos, const size_t nPositions,
                                         std::complex<double>** odata);

  template void FNM_EXPORT CalcCwField(const ApertureData<double>& data,
                                       const double* pos, const size_t nPositions,
                                       std::complex<double>** odata);
#endif
}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */

