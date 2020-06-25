/**
 * @file   fastgl.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Thu Jun  9 05:23:58 2016
 *
 * @brief Fast Gauss-Legendre integration
 *
 * *******************************************
 *    Copyright (C) 2014 by Ignace Bogaert   *
 * *******************************************
 *
 *  This software package is based on the paper
 *     I. Bogaert, "Iteration-Free Computation of Gauss-Legendre Quadrature Nodes and Weights",
 *     to be published in the SIAM Journal of Scientific Computing.
 *
 *  The main features of this software are:
 *  - Speed: due to the simple formulas and the O(1) complexity computation of individual Gauss-Legendre
 *    quadrature nodes and weights. This makes it compatible with parallel computing paradigms.
 *  - Accuracy: the error on the nodes and weights is within a few ulps (see the paper for details).
 *
 *  Disclaimer:
 *  THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 *  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 *  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 *  OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef FASTGL_H
#define FASTGL_H

#include <gl/gl_export.h>
#include <stddef.h>
#include <cmath>

/** @defgroup FastGL Fast Gauss-Legendre library
 *
 * This module is used for performing Gauss-Legendre integration
 *
 */

/** @addtogroup FastGL */
/*@{*/

// Functions for fastgl in double precision
namespace fastgl {

  /*! \brief Node-weigth pair
   *
   *
   * A struct for containing a Node-weight pair
   */
  struct GL_EXPORT QuadPair {
    double theta, weight;

    /**
     *
     * A function for getting the node in x-space
     *
     * @return
     */
    double x()
    {
      return cos(theta);
    }

    /**
     * A constructor
     *
     * @param t
     * @param w
     *
     * @return
     */
    QuadPair(double t, double w) : theta(t), weight(w) {}

    /**
     * Default ctor
     *
     *
     * @return
     */
    QuadPair() {}
  };

  /**
   *
   * Function for getting Gauss-Legendre nodes and weights
   * Theta values of the zeros are in [0,pi], and monotonically increasing.
   * The index of the zero k should always be in [1,n].
   * Compute a node-weight pair:
   *
   * @param n
   * @param k
   *
   * @return
   */
  GL_EXPORT QuadPair GLPair(size_t n, size_t k);
}

/*@}*/

#endif

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
