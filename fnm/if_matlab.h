/**
 * @file   if_matlab.h
 * @author Jens Munk Hansen <jmh@debian9laptop.parknet.dk>
 * @date   Mon Jun 18 21:42:28 2018
 *
 * @brief  Interface header used for saving Matlab for parsing e.g. __attribute__ and
 *         other keywords, which it cannot handle
 *
 *
 */

#ifndef __IF_MATLAB_H
#define __IF_MATLAB_H

# define FNM_STATIC_DEFINE
#ifdef __GNUC__
# define FNM_EXTERNAL_API
#elif defined(_MSC_VER)
# define FNM_EXTERNAL_API __declspec(dllexport)
#endif
# ifndef SPS_FCOMPLEX
#  ifdef __GNUC__
#   define SPS_FCOMPLEX float _Complex
#  elif defined(_MSC_VER)
#   define SPS_FCOMPLEX float  // Windows does not handle C99 complex numbers
#  endif
# endif
# include <fnm/if_fnm.h>
#endif
