#include <fnm/fnm_common.hpp>
#include <fnm/fnm_types.hpp>
#include <gl/gl.hpp>

namespace fnm {
  template <>
  void CalcWeightsAndAbcissaeSIMD<float>(const sysparm_t<float>* sysparm,
                                         float** __restrict uv_xs,
                                         float** __restrict uv_ws)
  {

    size_t nDiv = std::max<size_t>(sysparm->nDivW,sysparm->nDivH);

    *uv_xs    = (float*) _mm_malloc(4*nDiv*sizeof(float),16);
    *uv_ws    = (float*) _mm_malloc(4*nDiv*sizeof(float),16);

    // Common weights and abcissa
    for (size_t i = 0 ; i < nDiv ; i++) {
      gl::GLNode qp = gl::GL(nDiv,i);
      float value  = float(qp.value);
      float weight = float(qp.weight);
      // us (could be scaled by hw later)
      (*uv_xs)[i*4+0] = value;
      (*uv_xs)[i*4+1] = value;
      // vs (the same here, but could be scaled using hh)
      (*uv_xs)[i*4+2] = value;
      (*uv_xs)[i*4+3] = value;
      // uws
      (*uv_ws)[i*4+0]     = weight;
      (*uv_ws)[i*4+1]     = weight;
      // vws
      (*uv_ws)[i*4+2]     = weight;
      (*uv_ws)[i*4+3]     = weight;
    }
  }

  template <>
  void CalcWeightsAndAbcissaeScaledSIMD<float>(const sysparm_t<float>* sysparm,
      const sps::element_rect_t<float>* pElement,
      sps::deleted_aligned_array<float> &&uv_xs,
      sps::deleted_aligned_array<float> &&uv_ws)
  {

    size_t nDiv = std::max<size_t>(sysparm->nDivW,sysparm->nDivH);

    uv_xs = sps::deleted_aligned_array_create<float>(4*nDiv);
    uv_ws = sps::deleted_aligned_array_create<float>(4*nDiv);

    // Common weights and abcissa
    for (size_t i = 0 ; i < nDiv ; i++) {
      gl::GLNode qp = gl::GL(nDiv,i);
      float value  = float(qp.value);
      float weight = float(qp.weight);
      // values
      float tmp = pElement->hw * value;
      uv_xs[i*4+0] = tmp;
      uv_xs[i*4+1] = tmp;
      tmp = pElement->hh * value;
      uv_xs[i*4+2] = tmp;
      uv_xs[i*4+3] = tmp;

      // weights
      tmp = pElement->hw * weight;
      uv_ws[i*4+0] = tmp;
      uv_ws[i*4+1] = tmp;
      tmp = pElement->hh * weight;
      uv_ws[i*4+2] = tmp;
      uv_ws[i*4+3] = tmp;
    }
  }

  template <>
  void CalcWeightsAndAbcissaeSIMD<double>(const sysparm_t<double>* sysparm,
                                          double** __restrict uv_xs,
                                          double** __restrict uv_ws)
  {
    SPS_UNREFERENCED_PARAMETERS(sysparm,uv_xs,uv_ws);
  }

  template <>
  void CalcWeightsAndAbcissaeScaledSIMD<double>(const sysparm_t<double>* sysparm,
      const sps::element_rect_t<double>* pElement,
      sps::deleted_aligned_array<double> &&uv_xs,
      sps::deleted_aligned_array<double> &&uv_ws)
  {
    SPS_UNREFERENCED_PARAMETERS(sysparm,pElement,uv_xs,uv_ws);
  }

  template <class T>
  void CalcWeightsAndAbcissaeScaled(const sysparm_t<T>* sysparm,
                                    const sps::element_rect_t<T>& element,
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
      uxs[i]        = element.hw * T(qp.value);
      uweights[i]   = element.hw * T(qp.weight);
    }

    for (size_t i = 0 ; i < nDivH ; i++) {
      gl::GLNode qp = gl::GL(nDivH, i);
      // Conversion from double to T
      vxs[i]      = element.hh * T(qp.value);
      vweights[i] = element.hh * T(qp.weight);
    }
  }

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

  template void CalcWeightsAndAbcissaeScaled(const sysparm_t<float>* sysparm,
      const sps::element_rect_t<float>& element,
      sps::deleted_aligned_array<float> &&uxs,
      sps::deleted_aligned_array<float> &&uweights,
      sps::deleted_aligned_array<float> &&vxs,
      sps::deleted_aligned_array<float> &&vweights);

  template void CalcWeightsAndAbcissae(const sysparm_t<float>* sysparm,
                                       sps::deleted_aligned_array<float> &&uxs,
                                       sps::deleted_aligned_array<float> &&uweights,
                                       sps::deleted_aligned_array<float> &&vxs,
                                       sps::deleted_aligned_array<float> &&vweights);

  template void CalcWeightsAndAbcissaeSIMD<float>(const sysparm_t<float>* sysparm,
      float** __restrict uv_xs,
      float** __restrict uv_ws);


#ifdef FNM_DOUBLE_SUPPORT
  template void CalcWeightsAndAbcissaeScaled(const sysparm_t<double>* sysparm,
      const sps::element_rect_t<double>& element,
      sps::deleted_aligned_array<double> &&uxs,
      sps::deleted_aligned_array<double> &&uweights,
      sps::deleted_aligned_array<double> &&vxs,
      sps::deleted_aligned_array<double> &&vweights);

  template void CalcWeightsAndAbcissae(const sysparm_t<double>* sysparm,
                                       sps::deleted_aligned_array<double> &&uxs,
                                       sps::deleted_aligned_array<double> &&uweights,
                                       sps::deleted_aligned_array<double> &&vxs,
                                       sps::deleted_aligned_array<double> &&vweights);

  template void CalcWeightsAndAbcissaeSIMD<double>(const sysparm_t<double>* sysparm,
      double** __restrict uv_xs,
      double** __restrict uv_ws);
#endif
}
