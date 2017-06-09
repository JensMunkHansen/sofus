/**
 * @file   gl_lut.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Wed Jun 22 13:22:00 2016
 *
 * @brief Look up tables LUTs with weight and values for Gauss-Legendre knots
 *
 *
 */
#pragma once

#include <gl/config.h>
#include <gl/gl_export.h>

/** @addtogroup GL */
/*@{*/

namespace gl {
  /// Weight used for Gauss-Legendre integration
  extern const double* weights[_GL_LUT_TABLE_SIZE];
  /// Abcissae or coordinates used for Gauss-Legendre integration
  extern const double* abcissas[_GL_LUT_TABLE_SIZE];
}
/*@}*/

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
