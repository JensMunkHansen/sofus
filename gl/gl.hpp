/**
 * @file   gl.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Mon Jun 20 21:11:40 2016
 *
 * @brief  Gauss-Legendre integration library
 *
 *
 */

#pragma once

#include <gl/gl_export.h>
#include <gl/config.h>
#include <stddef.h>

/** @defgroup GL Gauss-Legendre library
 *  @brief Nodes and weight for quadrature
 *
 * @section gl_abstract Abstract
 *
 * This module is used for performing integration using Gauss-Legendre
 * quadrature. Nodes and weights are computed using an iteration-free
 * concept adapted from the paper.
 *
 * I. Bogaert, "Iteration-Free Computation of Gauss-Legendre Quadrature Nodes and Weights",
 * SIAM Journal of Scientific Computing.
 *
 * The first number of nodes and weights are precomputed using a
 * iterative approach and tabulated. The remaining nodes and weights
 * are computed using the above mentioned iteration-free concept.
 *
 * @section gl_introduction Introduction
 *
 * A very brief introduction to quadrature rules is given here
 * together with instruction of how to apply the library for numerical
 * integration.
 *
 * In numerical analysis, a quadrature rule is an approximation of the
 * definite integral of a function, usually stated as a weighted sum
 * of function values at specified points within the domain of
 * integration. An n-point Gaussian quadrature rule, named after Carl
 * Friedrich Gauss, is a quadrature rule constructed to yield an exact
 * result for polynomials of degree 2n − 1 or less by a suitable
 * choice of the points xi and weights wi for i = 1, ..., n. The
 * domain of integration for such a rule is conventionally taken as
 * [−1, 1], so the rule is stated as
 *
 * \f[
 * \int_{-1}^{1} f(x) dx = \sum_{i=1}^{n} \omega_i f(x_i)
 * \f]
 *
 * Gaussian quadrature as above will only produce good results if the
 * function \f$f(x)\f$ is well approximated by a polynomial function
 * within the range \f$[−1, 1]\f$. The method is not, for example,
 * suitable for functions with singularities. However, if the
 * integrated function can be written as \f$ f(x)=\omega (x)g(x)\,\f$,
 * where \f$g(x)\f$ is approximately polynomial and \f$\omega(x)\f$ is
 * known, then alternative weights \f$\omega_{i}'\f$ and points \f$
 * x_{i}'\f$ that depend on the weighting function \f$\omega(x)\f$ may
 * give better results, where
 *
 * \f[
 * \int_{-1}^{1}f(x)\,dx=\int_{-1}^{1}\omega (x)g(x)\,dx\approx \sum_{i=1}^{n}\omega_{i}'g(x_{i}').
 * \f]
 *
 * Common weighting functions include \f$ \omega (x)=1/{\sqrt
 * {1-x^{2}}}\,\f$ (Chebyshev–Gauss) and \f$\omega (x)=e^{-x^{2}}\f$
 * (Gauss–Hermite). It can be shown (see Press, et al., or Stoer and
 * Bulirsch) that the evaluation points xi are just the roots of a
 * polynomial belonging to a class of orthogonal polynomials.
 *
 * @subsection gl_gl_quadrature Gauss-Legendre quadrature
 *
 * For the simplest integration problem stated above, i.e. with
 * \f$\omega (x)=1\f$, the associated polynomials are Legendre
 * polynomials, \f$P_n(x)\f$, and the method is usually known as
 * Gauss–Legendre quadrature. With the n-th polynomial normalized to
 * give \f$P_n(1) = 1\f$, the i-th Gauss node, \f$x_i\f$, is the i-th
 * root of \f$P_n\f$; its weight is given by (Abramowitz & Stegun
 * 1972, p. 887)
 *
 * \f[
 * \omega_{i}={\frac {2}{\left(1-x_{i}^{2}\right)[P'_{n}(x_{i})]^{2}}}.
 * \f]
 *
 * @subsection gl_change_of_interval Change of interval
 *
 * An integral over [a, b] must be changed into an integral over \f$[−1,
 * 1]\f$ before applying the Gaussian quadrature rule. This change of
 * interval can be done in the following way:
 *
 * \f[
 * \int_{a}^{b}f(x)\,dx={\frac {b-a}{2}}\int_{-1}^{1}f\left({\frac {b-a}{2}}x+{\frac {a+b}{2}}\right)\,dx.
 * \f]
 *
 * Applying the Gaussian quadrature rule then results in the following approximation:
 *
 * \f[
 * \int_{a}^{b}f(x)\,dx={\frac {b-a}{2}}\sum _{i=1}^{n}\omega_{i}f\left({\frac {b-a}{2}}x_{i}+{\frac {a+b}{2}}\right).
 * \f]
 *
 * @section gl_usage Using the library
 */

/** @addtogroup GL */
/*@{*/

//! Gauss-Legendre interfaces and implementations
namespace gl {

  /*! @struct GLNode
    @brief Container for weight and nodes
   */
  struct GL_EXPORT GLNode {
    /// Nodes
    double value;
    /// Weight
    double weight;
  };

  /**
   * Gauss-Legendre quadrature using an n-point rule.
   *
   * @param n Quadrature order
   * @param f Integrand
   * @param data Pointer to user-defined data which will
                be passed to \p f every time it is called (as second parameter).
   * @param a Lower integration limit
   * @param b Upper integration limit
   *
   * @return
   */
  double GL_EXPORT GLQuad(size_t n, double (*f)(double,void*), void* data, double a, double b);

  /**
   *
   *  Purpose:
   *
   *    GL computes the kth GL pair of an n-point rule. It uses look-up tables for the n < 101.
   *
   *  Licensing:
   *
   *    This code is distributed under the BSD license.
   *
   *  Modified:
   *
   *    22 December 2015
   *
   *  Author:
   *
   *    Ignace Bogaert
   *
   *  Reference:
   *
   *    Ignace Bogaert,
   *    Iteration-free computation of Gauss-Legendre quadrature nodes and weights,
   *    SIAM Journal on Scientific Computing,
   *    Volume 36, Number 3, 2014, pages A1008-1026.
   *
   * The only function that needs to be public
   *
   * @param l The number of points in the given rule
   * @param k The index of the point to be returned
   *
   * @return The location and weight of the point
   */
  gl::GLNode GL_EXPORT GL(size_t l, size_t k);

  /**
   *  GLS computes the kth GL pair of an n-point rule using the
   *  formula rather than look-up tables.
   *
   *  Reference:
   *
   *    Ignace Bogaert,
   *    Iteration-free computation of Gauss-Legendre quadrature nodes and weights,
   *    SIAM Journal on Scientific Computing,
   *    Volume 36, Number 3, 2014, pages A1008-1026.
   *
   * @param n The number of points in the given rule
   * @param k The index of the point to be returned
   *
   * @return
   */
  gl::GLNode GL_EXPORT GLS(size_t n, size_t k);

  /**
   * This function computes the k'th zero of the BesselJ(0,x)
   *
   * @param k
   *
   * @return
   */
  double GL_EXPORT besseljzero(int k);

  /**
   * This function computes the square of BesselJ(1, BesselZero(0,k))
   *
   * @param k
   *
   * @return
   */
  double GL_EXPORT besselj1squared(int k);
}
/*@}*/

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
