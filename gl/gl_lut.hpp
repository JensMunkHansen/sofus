/**
 * @file   gl_lut.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Wed Jun 22 13:22:00 2016
 *
 * @brief Look up tables LUTs with weight and values for Gauss-Legendre knots
 *
 *
 */

/*
 *  This file is part of SOFUS.
 *
 *  SOFUS is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  SOFUS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with SOFUS.  If not, see <http://www.gnu.org/licenses/>.
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
