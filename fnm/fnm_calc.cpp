#include <fnm/fnm.hpp>
#include <fnm/fnm_calc.hpp>
#include <fnm/FnmSIMD.hpp>

#include <gl/gl.hpp>

#include <sps/smath.hpp>
#include <sps/mm_malloc.h>
#include <sps/memory>
#include <memory>
#include <sps/debug.h>

namespace fnm {

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

  // Cannot pass a unique_ptr by reference
  template <class T>
  void CalcCwField(const ApertureData<T>& data,
                   const T* pos, const size_t nPositions,
                   std::complex<T>** odata)
  {

    const T lambda = Aperture<T>::_sysparm.c / data.m_f0;

    const T k = T(M_2PI)/lambda;

    const size_t nDivW = Aperture<T>::_sysparm.nDivW;
    const size_t nDivH = Aperture<T>::_sysparm.nDivH;

    auto uweights = sps::deleted_aligned_array_create<T>(nDivW);
    auto vweights = sps::deleted_aligned_array_create<T>(nDivH);
    auto uxs      = sps::deleted_aligned_array_create<T>(nDivW);
    auto vxs      = sps::deleted_aligned_array_create<T>(nDivH);

    memset(uxs.get(),0,sizeof(T)*nDivW);
    memset(vxs.get(),0,sizeof(T)*nDivH);
    memset(uweights.get(),0,sizeof(T)*nDivW);
    memset(vweights.get(),0,sizeof(T)*nDivH);

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

    sps::point_t<T> point;

    size_t nElements = data.m_nelements;

    auto& elements = data.m_elements;

    for (size_t iPoint = 0 ; iPoint < nPositions ; iPoint++) {

      point[0] = pos[iPoint*3 + 0];
      point[1] = pos[iPoint*3 + 1];
      point[2] = pos[iPoint*3 + 2];

      std::complex<T> field = std::complex<T>(T(0.0),T(0.0));

      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {

        T apodization = data.m_apodizations[iElement];
        if (apodization != 0.0) {
          std::complex<T> field1 = std::complex<T>(T(0.0),T(0.0));

          for (size_t jElement = 0 ; jElement < data.m_nsubelements ; jElement++) {

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

          field.real(field.real() + real*cos(data.m_phases[iElement]) - imag*sin(data.m_phases[iElement]));
          field.imag(field.imag() + real*sin(data.m_phases[iElement]) + imag*cos(data.m_phases[iElement]));
        }
      }
      (*odata)[iPoint].real(field.real());
      (*odata)[iPoint].imag(field.imag());
    }
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

    auto uweights = sps::deleted_aligned_array_create<T>(nDivW);
    auto vweights = sps::deleted_aligned_array_create<T>(nDivH);
    auto uxs      = sps::deleted_aligned_array_create<T>(nDivW);
    auto vxs      = sps::deleted_aligned_array_create<T>(nDivH);

    memset(uxs.get(),0,sizeof(T)*nDivW);
    memset(vxs.get(),0,sizeof(T)*nDivH);

    memset(uweights.get(),0,sizeof(T)*nDivW);
    memset(vweights.get(),0,sizeof(T)*nDivH);

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

    sps::point_t<T> point;

    size_t nElements = data.m_nelements;

    auto& elements = data.m_elements;

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

          // Phases
          field.real(field.real() + real*cos(data.m_phases[iElement]) - imag*sin(data.m_phases[iElement]));
          field.imag(field.imag() + real*sin(data.m_phases[iElement]) + imag*cos(data.m_phases[iElement]));

        } /* if (apodization != 0.0) */
      } /* for (size_t iElement = 0 ; iElement < nElements ; iElement++) */

      (*odata)[iPoint].real(field.real());
      (*odata)[iPoint].imag(field.imag());
    }
    return 0;
  }

  template int FNM_EXPORT CalcCwFieldRef(const ApertureData<float>& data,
                                         const float* pos, const size_t nPositions,
                                         std::complex<float>** odata);

  template void FNM_EXPORT CalcCwField(const ApertureData<float>& data,
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

