/**
 * @file   bessel_lut.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Tue Jun 21 21:11:51 2016
 *
 * @brief  Look up tables LUTs for Bessel functions
 *
 *
 */
#pragma once

#include <gl/config.h>
#include <gl/gl_export.h>

namespace gl {
  extern double JZ[_BESSELJZERO_LUT_TABLE_SIZE];
  extern double J1[_BESSELJ_1_SQUARED_LUT_TABLE_SIZE];
}
