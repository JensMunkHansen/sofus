// TODO: Clean up - it is obsolete

#include <sofus/config.h>
#include <sofus/sofus_export.h>

#include <sofus/ApertureData.hpp>
#include <sps/math.h>
#include <sps/smath.hpp>
#include <sps/trigintrin.h>
#include <iostream>
#include <algorithm>

namespace sps {
  template <class T>
  void element_t<T>::edges(__m128* dx, __m128* dy, __m128* dz)
  {
    __m128 x,y,z;
    _m_vertices(&x,&y,&z);

    __m128 x1 = _mm_shuffle_ps(x, x, _MM_SHUFFLE(0,3,2,1));
    *dx = _mm_sub_ps(x,x1);

    __m128 y1 = _mm_shuffle_ps(y, y, _MM_SHUFFLE(0,3,2,1));
    *dy = _mm_sub_ps(y,y1);

    __m128 z1 = _mm_shuffle_ps(z, z, _MM_SHUFFLE(0,3,2,1));
    *dz = _mm_sub_ps(z,z1);
  }

#if INLINE_VERTICES_CALCULATOR
#else
  template <class T>
  void element_t<T>::_m_vertices(__m128* x0, __m128* y0, __m128* z0)
  {

    const float alpha = (float) this->euler.alpha;
    const float beta  = (float) this->euler.beta;
    const float gamma = (float) this->euler.gamma;

    // Sine and cosines
    __m128 v_cabc = _m_one_ps;
    __m128 v_sabc = _mm_setzero_ps();

    __m128 euler = _mm_set_ps(0.0,gamma,beta,alpha);

#pragma message("No inline vertice calculations")

    // Faster than 2 calls
    _mm_sin_cos_ps(euler,&v_sabc,&v_cabc);
    //_sabc = _mm_sin_ps(euler);
    //_cabc = _mm_cos_ps(euler);

    ALIGN16_BEGIN float cabc[4] ALIGN16_END;
    ALIGN16_BEGIN float sabc[4] ALIGN16_END;

    _mm_store_ps((float*)&cabc[0],v_cabc);
    _mm_store_ps((float*)&sabc[0],v_sabc);

    // Spread results out, we need 4 vertices anyway
#if HAVE_IMMINTRIN_H
    __m128 sa = _mm_broadcast_ss((float*)&sabc[0]);
    __m128 sb = _mm_broadcast_ss((float*)&sabc[1]);
    __m128 sc = _mm_broadcast_ss((float*)&sabc[2]);

    __m128 ca = _mm_broadcast_ss((float*)&cabc[0]);
    __m128 cb = _mm_broadcast_ss((float*)&cabc[1]);
    __m128 cc = _mm_broadcast_ss((float*)&cabc[2]);
#else
    __m128 sa = _mm_set1_ps(sabc[0]);
    __m128 sb = _mm_set1_ps(sabc[1]);
    __m128 sc = _mm_set1_ps(sabc[2]);

    __m128 ca = _mm_set1_ps(cabc[0]);
    __m128 cb = _mm_set1_ps(cabc[1]);
    __m128 cc = _mm_set1_ps(cabc[2]);
#endif

#if HAVE_IMMINTRIN_H
    // Original vertices (note no z-component)
    __m128 v_vertice_x = _mm_mul_ps(_mm_broadcast_ss((float*)&this->hw),_mm_set_ps(1.0,-1.0,-1.0,1.0f));
    __m128 v_vertice_y = _mm_mul_ps(_mm_broadcast_ss((float*)&this->hh),_mm_set_ps(-1.0,-1.0,1.0,1.0f));
#else
    __m128 v_vertice_x = _mm_mul_ps(_mm_set1_ps((float)this->hw),_mm_set_ps(1.0,-1.0,-1.0,1.0f));
    __m128 v_vertice_y = _mm_mul_ps(_mm_set1_ps((float)this->hh),_mm_set_ps(-1.0,-1.0,1.0,1.0f));
#endif

#if ROTATION_CONVENTION == ROTATION_CONVENTION_EULER_YXY
    // Rotation (only yxy is supported)
    __m128 v_vertice_xx =
      _mm_add_ps(
        _mm_mul_ps(
          _mm_mul_ps(sa,sb),
          v_vertice_y),
        _mm_mul_ps(
          _mm_sub_ps(
            _mm_mul_ps(ca,cc),
            _mm_mul_ps(
              cb,
              _mm_mul_ps(
                sa,
                sc))),
          v_vertice_x));

    __m128 v_vertice_yy =
      _mm_add_ps(
        _mm_mul_ps(
          cb,
          v_vertice_y),
        _mm_mul_ps(
          _mm_mul_ps(
            sb,
            sc),
          v_vertice_x));

    __m128 v_vertice_zz =
      _mm_sub_ps(
        _mm_mul_ps(
          _mm_mul_ps(
            ca,
            sb),
          v_vertice_y),
        _mm_mul_ps(
          _mm_add_ps(
            _mm_mul_ps(
              cc,
              sa),
            _mm_mul_ps(
              ca,
              _mm_mul_ps(
                cb,
                sc))),
          v_vertice_x));
#else
#error Unsupported rotation convention
#endif

    // Add origin
#if HAVE_IMMINTRIN_H
    *x0 = _mm_add_ps(v_vertice_xx,_mm_broadcast_ss((float*)&this->center[0]));
    *y0 = _mm_add_ps(v_vertice_yy,_mm_broadcast_ss((float*)&this->center[1]));
    *z0 = _mm_add_ps(v_vertice_zz,_mm_broadcast_ss((float*)&this->center[2]));
#else
    *x0 = _mm_add_ps(v_vertice_xx,_mm_set1_ps((float)this->center[0]));
    *y0 = _mm_add_ps(v_vertice_yy,_mm_set1_ps((float)this->center[1]));
    *z0 = _mm_add_ps(v_vertice_zz,_mm_set1_ps((float)this->center[2]));
#endif

  }
#endif

  template <class T>
  void element_t<T>::calc_vertices(sps::point_t<T> (*vertices)[4]) const
  {

    sps::point_t<T> hh_dir, hw_dir, normal;

    basis_vectors(hw_dir, this->euler, 0);
    basis_vectors(hh_dir, this->euler, 1);
    basis_vectors(normal, this->euler, 2);

    (*vertices)[0] = this->center + this->hw * hw_dir + this->hh * hh_dir;
    (*vertices)[1] = this->center - this->hw * hw_dir + this->hh * hh_dir;
    (*vertices)[2] = this->center - this->hw * hw_dir - this->hh * hh_dir;
    (*vertices)[3] = this->center + this->hw * hw_dir - this->hh * hh_dir;
  }

  inline std::ostream& operator<<(std::ostream& out, const float _array[8])
  {
    out << _array[0] << ", " << _array[1] << ", " << _array[2] << ", " << _array[3] << ", " << _array[4] << ", "
        << _array[5] << ", " << _array[6] << ", " << _array[7] << std::endl;
    return out;
  }

  template struct point_t<float>;
  template struct element_t<float>;
  
  template struct point_t<double>;
  template struct element_t<double>;
}
