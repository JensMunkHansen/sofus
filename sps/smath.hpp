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

#ifdef _MSC_VER
// TODO: Try out: alignas (16) std::array<char, 10>
# include <sps/aligned_array.hpp>
#else
# include <array>
#endif


#ifdef _WIN32
namespace std {
  template class SPS_EXPORT aligned_array<float,4>;
  template class SPS_EXPORT aligned_array<double,4>;
}
#endif

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
/*@{*/

namespace sps {

  /*! \brief Point type (aligned)
   *
   *
   *  Point as an aligned array
   */
#ifdef _MSC_VER
//    template <typename T>
//    struct SPS_EXPORT point_t: public std::aligned_array<T,4> {};
    template <typename T>
    struct SPS_EXPORT point_t : public std::aligned_array<T,4> {};
#else
  template <typename T>
  ALIGN16_BEGIN struct ALIGN16_END SPS_EXPORT point_t : public std::array<T,4> {};
#endif

  /*! \brief Rectangle
   *
   *
   *  Rectangle structure as an array of \ref sps::point_t
   */
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
    /// Boundary minimum
    point_t<T> min;
    /// Boundary maximum
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
    /// alpha
    T alpha;
    /// beta
    T beta;
    /// gamma
    T gamma;
    /// Dummy used by alignment
    T dummy;
  };
#else
  template <typename T>
  ALIGN16_BEGIN struct ALIGN16_END euler_t {
    /// alpha
    T alpha;
    /// beta
    T beta;
    /// gamma
    T gamma;
    /// Dummy used by alignment
    T dummy;
  };
#endif

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

#ifdef _WIN32
  template struct SPS_EXPORT point_t<float>;
  template struct SPS_EXPORT point_t<double>;
  template struct SPS_EXPORT euler_t<float>;
  template struct SPS_EXPORT euler_t<double>;
#endif

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
  void SPS_EXPORT basis_vectors_ps(T* vec0, T* vec1, T* vec2, const sps::euler_t<T>& euler);

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
