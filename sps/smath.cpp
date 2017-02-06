// TODO: Support multiple configuration files or make it header only
#include <sps/config.h>

#include <sps/smath.hpp>

#include <sps/trigintrin.h>

namespace sps {

  template <typename T>
  void basis_rotate(const sps::point_t<T>& input, const euler_t<T>& euler, sps::point_t<T>& output)
  {
    sps::point_t<T> vec0, vec1, vec2;
    basis_vectors<T>(vec0, euler, 0);
    basis_vectors<T>(vec1, euler, 1);
    basis_vectors<T>(vec2, euler, 2);

    output[0] = vec0[0]*input[0] + vec1[0]*input[1] + vec2[0]*input[2];
    output[1] = vec0[1]*input[0] + vec1[1]*input[1] + vec2[1]*input[2];
    output[2] = vec0[2]*input[0] + vec1[2]*input[1] + vec2[2]*input[2];
  }

  template <typename T>
  void basis_vectors(sps::point_t<T>& output, const euler_t<T>& euler, size_t index)
  {

    const T alpha = euler.alpha;
    const T beta  = euler.beta;
    const T gamma = euler.gamma;

    const T sa = sin(alpha);
    const T ca = cos(alpha);
    const T sb = sin(beta);
    const T cb = cos(beta);
    const T sc = sin(gamma);
    const T cc = cos(gamma);

#if ROTATION_CONVENTION == ROTATION_CONVENTION_EULER_ZYZ
#pragma message("Rotation used is zyz")
    switch (index) {
    case 0:
      output[0]= ca*cc-cb*sa*sc;
      output[1]=-cb*cc*sa-ca*sc;
      output[2]= sa*sb;
      break;
    case 1:
      output[0]= cc*sa+ca*cb*sc;
      output[1]= ca*cb*cc-sa*sc;
      output[2]=-ca*sb;
      break;
    case 2:
      output[0]= sb*sc;
      output[1]= cc*sb;
      output[2]= cb;
      break;
    }
#elif ROTATION_CONVENTION == ROTATION_CONVENTION_EULER_YXY
#pragma message("Rotation used is yxy")
    // Intrinsic rotations, y, x and y'
    switch (index) {
    case 0:
      // Mult by (a,b,c) to get x component
      output[0]= ca*cc - cb*sa*sc;
      output[1]= sb*sc;
      output[2]=-ca*cb*sc - cc*sa;
      break;
    case 1:
      output[0]=sa*sb;
      output[1]=cb;
      output[2]=ca*sb;
      break;
    case 2:
      output[0]= ca*sc + cb*cc*sa;
      output[1]=-cc*sb;
      output[2]= ca*cb*cc - sa*sc;
      break;
    }
#elif ROTATION_CONVENTION == ROTATION_CONVENTION_EXTRINSIC_XYZ
#pragma message("Rotation used is extrinsic xyz")
    switch (index) {
    case 0:
      output[0]= cb*cc;
      output[1]= cb*sc;
      output[2]=-sb;
      break;
    case 1:
      output[0]=-ca*sc + cc*sa*sb;
      output[1]= ca*cc + sa*sb*sc;
      output[2]= cb*sa;
      break;
    case 2:
      output[0]=ca*cc*sb + sa*sc;
      output[1]=ca*sb*sc - cc*sa;
      output[2]=ca*cb;
      break;
    }
#endif
  }

  template <>
  void basis_vectors(float* vec0, float* vec1, float* vec2, const sps::euler_t<float>& euler)
  {

    __m128 a;
#ifdef _WIN32
    // Euler is not aligned on Windows
    a = _mm_loadu_ps((float*)&euler);
#else
    a = _mm_load_ps((float*)&euler);
#endif
    v4f c,s;
    _mm_sin_cos_ps(a,&s.v,&c.v);

    v4f* p0 = (v4f*)vec0;
    v4f* p1 = (v4f*)vec1;
    v4f* p2 = (v4f*)vec2;

#if ROTATION_CONVENTION == ROTATION_CONVENTION_EULER_ZYZ
    p0->f32[0] =  c.f32[0]*c.f32[2]          - c.f32[2]*s.f32[0]*s.f32[2];
    p0->f32[1] = -c.f32[1]*c.f32[2]*s.f32[0] - c.f32[0]*s.f32[2];
    p0->f32[2] =  s.f32[0]*s.f32[1];

    p1->f32[0] =  c.f32[2]*s.f32[0]          + c.f32[2]*c.f32[1]*s.f32[2];
    p1->f32[1] =  c.f32[0]*c.f32[1]*c.f32[2] - s.f32[0]*s.f32[2];
    p1->f32[2] = -c.f32[0]*s.f32[1];

    p2->f32[0] = s.f32[1]*s.f32[2];
    p2->f32[1] = c.f32[2]*s.f32[1];
    p2->f32[2] = c.f32[1];
#elif ROTATION_CONVENTION == ROTATION_CONVENTION_EULER_YXY
    p0->f32[0] =  c.f32[0]*c.f32[2]          - c.f32[2]*s.f32[0]*s.f32[2];
    p0->f32[1] =  s.f32[1]*s.f32[2];
    p0->f32[2] = -c.f32[0]*c.f32[1]*s.f32[2] - c.f32[2]*s.f32[0];

    p1->f32[0] =  s.f32[0]*s.f32[1];
    p1->f32[1] =  c.f32[1];
    p1->f32[2] =  c.f32[0]*s.f32[1];

    p2->f32[0] = c.f32[0]*s.f32[2] + c.f32[1]*c.f32[2]*s.f32[0];
    p2->f32[1] = -c.f32[2]*s.f32[1];
    p2->f32[2] = c.f32[0]*c.f32[1]*c.f32[2] - s.f32[0]*s.f32[2];
#endif
  }

  template <>
  void basis_vectors(double* vec0, double* vec1, double* vec2, const sps::euler_t<double>& euler)
  {

    // TODO: TEST ME
    const double alpha = euler.alpha;
    const double beta  = euler.beta;
    const double gamma = euler.gamma;

    const double sa = sin(alpha);
    const double ca = cos(alpha);
    const double sb = sin(beta);
    const double cb = cos(beta);
    const double sc = sin(gamma);
    const double cc = cos(gamma);

    // Intrinsic rotations, y, x and y'
#if ROTATION_CONVENTION == ROTATION_CONVENTION_EULER_ZYZ
    vec0[0]= ca*cc-cb*sa*sc;
    vec0[1]=-cb*cc*sa-ca*sc;
    vec0[2]= sa*sb;

    vec1[0]= cc*sa+ca*cb*sc;
    vec1[1]= ca*cb*cc-sa*sc;
    vec2[2]=-ca*sb;

    vec2[0]= sb*sc;
    vec2[1]= cc*sb;
    vec2[2]= cb;
#elif ROTATION_CONVENTION == ROTATION_CONVENTION_EULER_YXY
    vec0[0]= ca*cc - cb*sa*sc;
    vec0[1]= sb*sc;
    vec0[2]=-ca*cb*sc - cc*sa;

    vec1[0]=sa*sb;
    vec1[1]=cb;
    vec1[2]=ca*sb;

    vec2[0]= ca*sc + cb*cc*sa;
    vec2[1]=-cc*sb;
    vec2[2]= ca*cb*cc - sa*sc;
#endif
  }

#ifdef _WIN32
  template class std::aligned_array<float,4U>;
  template class std::aligned_array<double,4U>;
#endif

  template struct euler_t<float>;
  template struct point_t<float>;

#ifdef _WIN32
// Not possible to move to fnm library
  template class std::aligned_array<sps::point_t<float>,4U>;
  template class std::aligned_array<sps::point_t<double>,4U>;
#endif

  template point_t<float> operator-(const point_t<float> &a, const point_t<float> &b);
  template point_t<float> operator+(const point_t<float> &a, const point_t<float> &b);
  template point_t<float> operator*(const float &a, const point_t<float> &b);
  template float dot(const point_t<float> &a, const point_t<float> &b);
  template point_t<float> cross(const point_t<float> &a, const point_t<float> &b);
  template float norm(const point_t<float> &a);
  template float dist_to_line(const point_t<float>& point, const point_t<float>& pointOnLine,
                              const point_t<float>& direction);

  template void SPS_EXPORT basis_vectors(sps::point_t<float>& output, const sps::euler_t<float>& euler, size_t index);

  template std::ostream& operator<<(std::ostream& out, const point_t<float>& point);
  template void SPS_EXPORT basis_vectors(float* vec0, float* vec1, float* vec2, const sps::euler_t<float>& euler);

  template void SPS_EXPORT basis_rotate(const sps::point_t<float>& input, const sps::euler_t<float>& euler,
                                        sps::point_t<float>& output);

  template void SPS_EXPORT dists_most_distant_and_closest(const sps::bbox_t<float> &box0,
      const sps::bbox_t<float> &box1,
      float* distNear,
      float* distFar);


  template struct euler_t<double>;
  template struct point_t<double>;

  template struct element_t<float>;
  template struct element_t<double>;


  template point_t<double> operator-(const point_t<double> &a, const point_t<double> &b);
  template point_t<double> operator+(const point_t<double> &a, const point_t<double> &b);
  template point_t<double> operator*(const double &a, const point_t<double> &b);
  template double dot(const point_t<double> &a, const point_t<double> &b);
  template point_t<double> cross(const point_t<double> &a, const point_t<double> &b);
  template double norm(const point_t<double> &a);
  template double dist_to_line(const point_t<double>& point, const point_t<double>& pointOnLine,
                               const point_t<double>& direction);
  template void SPS_EXPORT basis_vectors(sps::point_t<double>& output, const sps::euler_t<double>& euler, size_t index);
  template std::ostream& operator<<(std::ostream& out, const point_t<double>& point);
  template void SPS_EXPORT basis_vectors(double* vec0, double* vec1, double* vec2, const sps::euler_t<double>& euler);

  template void SPS_EXPORT basis_rotate(const sps::point_t<double>& input, const sps::euler_t<double>& euler,
                                        sps::point_t<double>& output);


  template void SPS_EXPORT dists_most_distant_and_closest(const sps::bbox_t<double> &box0,
      const sps::bbox_t<double> &box1,
      double* distNear,
      double* distFar);
}



/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
