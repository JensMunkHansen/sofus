/**
 * @file   smath.hpp
 * @author Jens Munk Hansen <jmh@jmhlaptop>
 * @date   Sat Oct 10 18:41:43 2015
 * 
 * @brief  Simple math
 * 
 * 
 */

#pragma once
#include <sps/cenv.h>
#include <sps/sps_export.h>
#include <sps/math.h>

#ifdef _MSC_VER
# include <sps/aligned_array.hpp>
#else
# include <array>
#endif

#include <sps/trigintrin.h>

// Remove eventually
#include <iostream>


#ifdef _WIN32
template <class T>
T signum(const T& x) {
  return T((T(0) < x) - (x < T(0)));
}
#else
template <typename T> inline constexpr
int signum(T x, std::false_type is_signed) {
  return T(0) < x;
}

template <typename T> inline constexpr
int signum(T x, std::true_type is_signed) {
  return (T(0) < x) - (x < T(0));
}

template <typename T> inline constexpr
int signum(T x) {
  return signum(x, std::is_signed<T>());
}
#endif

namespace sps {

#ifdef _MSC_VER
  template <typename T>
  struct SPS_EXPORT point_t : public std::aligned_array<T,4> {};
#else
  template <typename T>
  ALIGN16_BEGIN struct ALIGN16_END SPS_EXPORT point_t : public std::array<T,4> {};
#endif

#ifdef _MSC_VER
  template <typename T>
  struct SPS_EXPORT rect_t : public std::aligned_array<point_t<T> ,4> {};
#else
  template <typename T>
  struct rect_t : public std::array<point_t<T> ,4> {};
#endif
  
	
  /*! \brief Bounding box
   *
   *
   *  The bounding box is used for fast parallel computation and to
   *  optimize the cache usage.
   */
  template <typename T>
  struct SPS_EXPORT bbox_t {
    point_t<T> min;
    point_t<T> max;
  };

  /*! \brief Euler angles
   *
   *
   *  Euler angles, the z-x-z' convention is used.
   */
#ifdef _MSC_VER
  template <typename T>
  struct SPS_EXPORT euler_t {
    T alpha;
    T beta;
    T gamma;
    T dummy;
  };
#else
  template <typename T>
  ALIGN16_BEGIN struct ALIGN16_END euler_t {
		T alpha;
		T beta;
		T gamma;
		T dummy;
  };
#endif
	
  template <class T>
  inline T dist_point_to_point(const point_t<T>& a, const point_t<T>& b) {
    return sqrt(SQUARE(a[0]-b[0])+SQUARE(a[1]-b[1])+SQUARE(a[2]-b[2]));
  }
  
  template <typename T>
  inline T dot(const point_t<T> &a, const point_t<T> &b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  }

  template <typename T>
  inline point_t<T> operator-(const point_t<T> &a, const point_t<T> &b)  {
    point_t<T> c;
    c[0] = a[0]-b[0];
    c[1] = a[1]-b[1];
    c[2] = a[2]-b[2];
    
    return c;
  }
  
  template <typename T>
  inline point_t<T> operator+(const point_t<T> &a, const point_t<T> &b) {
    point_t<T> c;
    c[0] = a[0] + b[0];
    c[1] = a[1] + b[1];
    c[2] = a[2] + b[2];
    return c;
  }
  
  template <typename T>
  inline point_t<T> cross(const point_t<T> &a, const point_t<T> &b) {
    point_t<T> c;
    c[0] = a[1]*b[2] - b[1]*a[2];
    c[1] = a[2]*b[0] - b[2]*a[0];
    c[2] = a[0]*b[1] - b[0]*a[1];
    return c;
  }

  template <typename T>
  inline point_t<T> operator*(const T &a, const point_t<T> &b) {
    point_t<T> c;
    c[0] = a*b[0];
    c[1] = a*b[1];
    c[2] = a*b[2];
    return c;
  }

  template <typename T>
  inline T norm(const point_t<T> &a) {
    return sqrt(dot(a,a));
  }

  template <typename T>
  inline T dist_to_line(const point_t<T>& point, const point_t<T>& pointOnLine,
                        const point_t<T>& direction) {
    return norm(cross(direction, point - pointOnLine));
  }

  template <typename T>
	inline T sgn_dist_to_plane(const point_t<T>& point, const point_t<T>& pointOnPlane,
														 const point_t<T>& unitNormal)
	{
		return dot((point - pointOnPlane),unitNormal);
	}
	
  template <typename T>
  inline sps::point_t<T> clamp_vector(const sps::point_t<T> &point, const sps::bbox_t<T> &box) {
		sps::point_t<T> clamped;
    clamped[0] = (point[0] < box.min[0]) ? box.min[0] : (point[0] > box.max[0]) ? box.max[0] : point[0];
    clamped[1] = (point[1] < box.min[1]) ? box.min[1] : (point[1] > box.max[1]) ? box.max[1] : point[1];
    clamped[2] = (point[2] < box.min[2]) ? box.min[2] : (point[2] > box.max[2]) ? box.max[2] : point[2];
    return clamped;
  }
  
  template <typename T>
  inline sps::point_t<T> nearest_point_on_bbox(const sps::point_t<T> &point, const sps::bbox_t<T> &box) {
    return clamp_vector(point,box);
  }
  
  template <typename T>
  inline sps::point_t<T> farthest_point_on_bbox(const sps::point_t<T> &point, const sps::bbox_t<T> &box) {
		sps::point_t<T> farthest;
    
    farthest[0] = (point[0] < box.min[0]) ? box.max[0] :
      (point[0] > box.max[0]) ? box.min[0] : 
      (point[0] - box.min[0]) > (box.max[0] - point[0]) ? box.min[0] : box.max[0];
    farthest[1] = (point[1] < box.min[1]) ? box.max[1] :
      (point[1] > box.max[1]) ? box.min[1] :
      (point[1] - box.min[1]) > (box.max[1] - point[1]) ? box.min[1] : box.max[1];
    farthest[2] = (point[2] < box.min[2]) ? box.max[2] :
      (point[2] > box.max[2]) ? box.min[2] :
      (point[2] - box.min[2]) > (box.max[2] - point[2]) ? box.min[2] : box.max[2];
    return farthest;
  }

  /** 
   * Function for returning the basis vector for a given coordinate
   * system defined using 3 Euler angles according to the the z-x-z'
   * convention.
   * 
   * @param output 
   * @param euler 
   * @param index 
   */
  template <typename T>
  void basis_vectors(sps::point_t<T>& output, const sps::euler_t<T>& euler, size_t index);

  template <typename T>
  void basis_vectors_ps(T* vec0, T* vec1, T* vec2, const sps::euler_t<T>& euler);

  template <typename T>
  inline std::ostream& operator<<(std::ostream& out, const point_t<T>& point) {
    out << "x: " << point[0] << " y: " << point[1] << " z: " << point[2] << std::endl;
    return out;
  }	
}

/* Local variables: */
/* indent-tab-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
