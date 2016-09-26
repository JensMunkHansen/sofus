#include <fnm/fnm.hpp>
#include <fnm/fnm_calc.hpp>

#include <gl/gl.hpp>

#include <sps/mm_malloc.h>
#include <memory>

namespace fnm {

  template <class T>
  STATIC_INLINE_BEGIN std::complex<T> CalcSingle(const T& s1,
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

    // integral width
    std::complex<T> intW = std::complex<T>(T(0.0),T(0.0));

    for (size_t iu = 0 ; iu < nUs ; iu++) {

      T s = sm * uxs[iu] + sp;
      T s22 = SQUARE(s);

      T argw = -k * sqrt(s2 + z2 + l2);
      T real = uweights[iu] * (cos(argw) - carg) / (s22+l2);
      T imag = uweights[iu] * (sin(argw) - sarg) / (s22+l2);
      intW.real(real + intW.real());
      intW.imag(imag + intW.imag());
    }
    intW *= sm * l / (T(M_2PI)*k);
    T realW = intW.real();
    intW.real(intW.imag());
    intW.imag(-realW);

    return intW;
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

  template <class T>
  void CalcCwField(const ApertureData<T>& data,
                   const T* pos, const size_t nPositions,
                   std::complex<T>** odata)
  {

    const T lambda = Aperture<T>::_sysparm.c / data.m_f0;

    const T k = T(M_2PI)/lambda;

    size_t nDivW = Aperture<T>::_sysparm.nDivW;
    size_t nDivH = Aperture<T>::_sysparm.nDivH;

    auto del = [](T* p) {
      _mm_free(p);
    };
    std::unique_ptr<T[], decltype(del)> uweights((T*)_mm_malloc(sizeof(T)*nDivW,16), del);
    std::unique_ptr<T[], decltype(del)> vweights((T*)_mm_malloc(sizeof(T)*nDivH,16), del);
    std::unique_ptr<T[], decltype(del)> uxs((T*)_mm_malloc(sizeof(T)*nDivW,16), del);
    std::unique_ptr<T[], decltype(del)> vxs((T*)_mm_malloc(sizeof(T)*nDivH,16), del);

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

    std::vector<std::vector<element_t<T> > >& elements = *data.m_elements;

    for (size_t iPoint = 0 ; iPoint < nPositions ; iPoint++) {

      point[0] = pos[iPoint*3 + 0];
      point[1] = pos[iPoint*3 + 1];
      point[2] = pos[iPoint*3 + 2];

      std::complex<T> field = std::complex<T>(T(0.0),T(0.0));

      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {

        std::complex<T> field1 = std::complex<T>(T(0.0),T(0.0));

        for (size_t jElement = 0 ; jElement < elements[iElement].size() ; jElement++) {

          element_t<T> element = elements[iElement][jElement];

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

        T real = field1.real();
        T imag = field1.imag();

        field.real(field.real() + real*cos(data.m_phases[iElement]) - imag*sin(data.m_phases[iElement]));
        field.imag(field.imag() + real*sin(data.m_phases[iElement]) + imag*cos(data.m_phases[iElement]));
      }
      (*odata)[iPoint].real(field.real());
      (*odata)[iPoint].imag(field.imag());
    }
  }

  template <class T>
  void CalcCwFieldRef(const ApertureData<T>& data,
                      const T* pos, const size_t nPositions,
                      std::complex<T>** odata)
  {
    const T lambda = Aperture<T>::_sysparm.c / data.m_f0;

    const T k = T(M_2PI)/lambda;

    size_t nDivW = Aperture<T>::_sysparm.nDivW;
    size_t nDivH = Aperture<T>::_sysparm.nDivH;

    auto del = [](T* p) {
      _mm_free(p);
    };
    std::unique_ptr<T[], decltype(del)> uweights((T*)_mm_malloc(sizeof(T)*nDivW,16), del);
    std::unique_ptr<T[], decltype(del)> vweights((T*)_mm_malloc(sizeof(T)*nDivH,16), del);
    std::unique_ptr<T[], decltype(del)> uxs((T*)_mm_malloc(sizeof(T)*nDivW,16), del);
    std::unique_ptr<T[], decltype(del)> vxs((T*)_mm_malloc(sizeof(T)*nDivH,16), del);

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

    std::vector<std::vector<element_t<T> > >& elements = *data.m_elements;

    for (size_t iPoint = 0 ; iPoint < nPositions ; iPoint++) {

      point[0] = pos[iPoint*3 + 0];
      point[1] = pos[iPoint*3 + 1];
      point[2] = pos[iPoint*3 + 2];

      std::complex<T> field = std::complex<T>(T(0.0),T(0.0));

      for (size_t iElement = 0 ; iElement < nElements ; iElement++) {

        std::complex<T> field1 = std::complex<T>(T(0.0),T(0.0));

        for (size_t jElement = 0 ; jElement < elements[iElement].size() ; jElement++) {

          element_t<T> element = elements[iElement][jElement];

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

          T s = fabs(v) + element.hh;

          if (fabs(u) > element.hw) {
            field1 += CalcSingle<T>(fabs(u)-element.hw, fabs(u)           , s, z, k, uxs.get(),uweights.get(),nDivW);
            s = fabs(fabs(v) - element.hh);
            field1 += CalcSingle<T>(fabs(u),            fabs(u)+element.hw, s, z, k, uxs.get(),uweights.get(),nDivW);
            field1 += CalcSingle<T>(fabs(u)-element.hw, fabs(u)           , s, z, k, uxs.get(),uweights.get(),nDivW);
          } else {
            field1 += CalcSingle<T>(0,                  fabs(u)+element.hw, s, z, k, uxs.get(),uweights.get(),nDivW);
            field1 += CalcSingle<T>(0,                  element.hw-fabs(u), s, z, k, uxs.get(),uweights.get(),nDivW);
            s = element.hh - fabs(v);
            field1 += CalcSingle<T>(0,                  fabs(u)+element.hw, s, z, k, uxs.get(),uweights.get(),nDivW);
            field1 += CalcSingle<T>(0,                  element.hw-fabs(u), s, z, k, uxs.get(),uweights.get(),nDivW);
          }

          s = fabs(u) + element.hw;
          if (fabs(v) > element.hh) {
            field1 += CalcSingle<T>(fabs(v)-element.hh, fabs(v)           , s, z, k, vxs.get(),uweights.get(),nDivH);
            field1 += CalcSingle<T>(fabs(v),            fabs(v)+element.hh, s, z, k, vxs.get(),uweights.get(),nDivH);
            s = element.hw - fabs(u);
            field1 += CalcSingle<T>(fabs(v)-element.hh, fabs(v)           , s, z, k, vxs.get(),uweights.get(),nDivH);
            field1 += CalcSingle<T>(fabs(v),            fabs(v)+element.hh, s, z, k, vxs.get(),uweights.get(),nDivH);
          } else {
            field1 += CalcSingle<T>(0,                  fabs(v)+element.hh, s, z, k, vxs.get(),vweights.get(),nDivH);
            s = element.hw - fabs(u);
            field1 += CalcSingle<T>(0,                  fabs(v)+element.hh, s, z, k, vxs.get(),vweights.get(),nDivH);
            s = fabs(u) + element.hw;
            field1 += CalcSingle<T>(0,                  element.hh-fabs(v), s, z, k, vxs.get(),vweights.get(),nDivH);
            s = element.hw - fabs(u);
            field1 += CalcSingle<T>(0,                  element.hh-fabs(v), s, z, k, vxs.get(),vweights.get(),nDivH);
          }
        }

        T real = field1.real();
        T imag = field1.imag();

        field.real(field.real() + real*cos(data.m_phases[iElement]) - imag*sin(data.m_phases[iElement]));
        field.imag(field.imag() + real*sin(data.m_phases[iElement]) + imag*cos(data.m_phases[iElement]));
      }
      (*odata)[iPoint].real(field.real());
      (*odata)[iPoint].imag(field.imag());
    }
  }

  template void FNM_EXPORT CalcCwFieldRef(const ApertureData<float>& data,
                                          const float* pos, const size_t nPositions,
                                          std::complex<float>** odata);

  template void FNM_EXPORT CalcCwFieldRef(const ApertureData<double>& data,
                                          const double* pos, const size_t nPositions,
                                          std::complex<double>** odata);

  template void FNM_EXPORT CalcCwField(const ApertureData<float>& data,
                                       const float* pos, const size_t nPositions,
                                       std::complex<float>** odata);

  template void FNM_EXPORT CalcCwField(const ApertureData<double>& data,
                                       const double* pos, const size_t nPositions,
                                       std::complex<double>** odata);

}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */

