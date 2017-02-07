/**
 * @file   smath.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
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

#include <sps/smath_types.hpp>

// Remove eventually
#include <iostream>

#ifdef _WIN32
/**
 * Signum function
 *
 * @param x
 *
 * @return
 */
template <class T>
T signum(const T& x)
{
  return T((T(0) < x) - (x < T(0)));
}
#else
/**
 * Signum function
 *
 * @param x
 * @param is_signed
 *
 * @return
 */
template <typename T> inline constexpr
int signum(T x, std::false_type is_signed)
{
  return T(0) < x;
}

/**
 * Signum function
 *
 * @param x
 * @param is_signed
 *
 * @return
 */
template <typename T> inline constexpr
int signum(T x, std::true_type is_signed)
{
  return (T(0) < x) - (x < T(0));
}

/**
 * Signum function
 *
 * @param x
 *
 * @return
 */
template <typename T> inline constexpr
int signum(T x)
{
  return signum(x, std::is_signed<T>());
}
#endif

/** @addtogroup SPS */
namespace sps {


  /**
   * Distance from point to point
   *
   * @param a
   * @param b
   *
   * @return
   */
  template <class T>
  inline T dist_point_to_point(const point_t<T>& a, const point_t<T>& b)
  {
    return sqrt(SQUARE(a[0]-b[0])+SQUARE(a[1]-b[1])+SQUARE(a[2]-b[2]));
  }

  /**
   * Dot product of vectors or points
   *
   * @param a
   * @param b
   *
   * @return
   */
  template <typename T>
  inline T dot(const point_t<T> &a, const point_t<T> &b)
  {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  }

  /**
   * Subtract two points
   *
   * @param a
   * @param b
   *
   * @return a - b
   */
  template <typename T>
  inline point_t<T> operator-(const point_t<T> &a, const point_t<T> &b)
  {
    point_t<T> c;
    c[0] = a[0]-b[0];
    c[1] = a[1]-b[1];
    c[2] = a[2]-b[2];

    return c;
  }

  /**
   * Add two points
   *
   * @param a
   * @param b
   *
   * @return a + b
   */
  template <typename T>
  inline point_t<T> operator+(const point_t<T> &a, const point_t<T> &b)
  {
    point_t<T> c;
    c[0] = a[0] + b[0];
    c[1] = a[1] + b[1];
    c[2] = a[2] + b[2];
    return c;
  }

  /**
   * Cross product of two points or vectors
   *
   * @param a
   * @param b
   *
   * @return a x b
   */
  template <typename T>
  inline point_t<T> cross(const point_t<T> &a, const point_t<T> &b)
  {
    point_t<T> c;
    c[0] = a[1]*b[2] - b[1]*a[2];
    c[1] = a[2]*b[0] - b[2]*a[0];
    c[2] = a[0]*b[1] - b[0]*a[1];
    return c;
  }

  /**
   * Scalar multiplication of a vector or point
   *
   * @param a scalar
   * @param b point
   *
   * @return
   */
  template <typename T>
  inline point_t<T> operator*(const T &a, const point_t<T> &b)
  {
    point_t<T> c;
    c[0] = a*b[0];
    c[1] = a*b[1];
    c[2] = a*b[2];
    return c;
  }

  /**
   * Norm of a vector or point
   *
   * @param a
   *
   * @return |a|
   */
  template <typename T>
  inline T norm(const point_t<T> &a)
  {
    return sqrt(dot(a,a));
  }

  /**
   * Distance from point to line
   *
   * @param point
   * @param pointOnLine
   * @param direction
   *
   * @return
   */
  template <typename T>
  inline T dist_to_line(const point_t<T>& point, const point_t<T>& pointOnLine,
                        const point_t<T>& direction)
  {
    return norm(cross(direction, point - pointOnLine));
  }

  /**
   * Signed distance from point to plane
   *
   * @param point
   * @param pointOnPlane
   * @param unitNormal
   *
   * @return
   */
  template <typename T>
  inline T sgn_dist_to_plane(const point_t<T>& point, const point_t<T>& pointOnPlane,
                             const point_t<T>& unitNormal)
  {
    return dot((point - pointOnPlane),unitNormal);
  }

  /**
   * Clamp a vector inside a box
   *
   * @param point
   * @param box
   *
   * @return
   */
  template <typename T>
  inline sps::point_t<T> clamp_vector(const sps::point_t<T> &point, const sps::bbox_t<T> &box)
  {
    sps::point_t<T> clamped;
    clamped[0] = (point[0] < box.min[0]) ? box.min[0] : (point[0] > box.max[0]) ? box.max[0] : point[0];
    clamped[1] = (point[1] < box.min[1]) ? box.min[1] : (point[1] > box.max[1]) ? box.max[1] : point[1];
    clamped[2] = (point[2] < box.min[2]) ? box.min[2] : (point[2] > box.max[2]) ? box.max[2] : point[2];
    return clamped;
  }

  /**
   * Nearest point on a box (from a point)
   *
   * @param point
   * @param box
   *
   * @return
   */
  template <typename T>
  inline sps::point_t<T> nearest_point_on_bbox(const sps::point_t<T> &point, const sps::bbox_t<T> &box)
  {
    return clamp_vector(point,box);
  }

  /**
   * Farthest point on a box (from a point)
   *
   * @param point
   * @param box
   *
   * @return
   */
  template <typename T>
  inline sps::point_t<T> farthest_point_on_bbox(const sps::point_t<T> &point, const sps::bbox_t<T> &box)
  {
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

  template <typename T>
  inline void dists_most_distant_and_closest(const sps::bbox_t<T> &box0,
      const sps::bbox_t<T> &box1,
      T* distNear,
      T* distFar)
  {
    // Corners
    T boundaries[6];
    sps::point_t<T> border_points0[8];
    sps::point_t<T> border_points1[8];

    memcpy(&boundaries[0],&box0.min[0],3*sizeof(T));
    memcpy(&boundaries[3],&box0.max[0],3*sizeof(T));
    for (size_t i = 0 ; i < 2 ; i++) {
      for (size_t j = 0 ; j < 2 ; j++) {
        for (size_t k = 0 ; k < 2 ; k++) {
          border_points0[i*4+j*2+k][0] = boundaries[i*3];   // 0,3
          border_points0[i*4+j*2+k][1] = boundaries[1+j*3]; // 1,4
          border_points0[i*4+j*2+k][2] = boundaries[2+k*3]; // 2,5
        }
      }
    }

    memcpy(&boundaries[0],&box1.min[0],3*sizeof(T));
    memcpy(&boundaries[3],&box1.max[0],3*sizeof(T));
    for (size_t i = 0 ; i < 2 ; i++) {
      for (size_t j = 0 ; j < 2 ; j++) {
        for (size_t k = 0 ; k < 2 ; k++) {
          border_points1[i*4+j*2+k][0] = boundaries[i*3];   // 0,3
          border_points1[i*4+j*2+k][1] = boundaries[1+j*3]; // 1,4
          border_points1[i*4+j*2+k][2] = boundaries[2+k*3]; // 2,5
        }
      }
    }

    *distNear = std::numeric_limits<T>::max();
    *distFar  = std::numeric_limits<T>::min();

    for (size_t iBorderPoint = 0 ; iBorderPoint < 8 ; iBorderPoint++) {
      sps::point_t<T> corner0 = border_points0[iBorderPoint];
      sps::point_t<T> near0 = nearest_point_on_bbox(corner0,
                              box1);
      for (size_t jBorderPoint = 0 ; jBorderPoint < 8 ; jBorderPoint++) {
        sps::point_t<T> corner1 = border_points1[jBorderPoint];

        *distFar = std::max<T>(*distFar, dist_point_to_point<T>(corner0,corner1));

        sps::point_t<T> near1 = nearest_point_on_bbox(near0,box0);
        *distNear = std::min<T>(*distNear, dist_point_to_point<T>(near0,near1));
      }
    }
  }

#ifdef _WIN32
  template struct SPS_EXPORT point_t<float>;
  template struct SPS_EXPORT point_t<double>;
  template struct SPS_EXPORT euler_t<float>;
  template struct SPS_EXPORT euler_t<double>;
#endif

  template <typename T>
  void SPS_EXPORT basis_rotate(const sps::point_t<T>& input, const euler_t<T>& euler, sps::point_t<T>& output);

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
  void SPS_EXPORT basis_vectors(sps::point_t<T>& output, const sps::euler_t<T>& euler, size_t index);

  /**
  * Function for returning the basis vector for a given coordinate
  * system defined using 3 Euler angles according to the the z-x-z'
  * convention. Using SIMD for packed singles
   *
   * @param vec0
   * @param vec1
   * @param vec2
   * @param euler
   */
  template <typename T>
  void SPS_EXPORT basis_vectors(T* vec0, T* vec1, T* vec2, const sps::euler_t<T>& euler);

  /**
   * Operator for printing points to a stream
   *
   * @param out
   * @param point
   *
   * @return
   */
  template <typename T>
  inline std::ostream& operator<<(std::ostream& out, const point_t<T>& point)
  {
    out << "x: " << point[0] << " y: " << point[1] << " z: " << point[2] << std::endl;
    return out;
  }

}

/*@}*/

/* Local variables: */
/* indent-tab-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */