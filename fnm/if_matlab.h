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
# define FNM_EXPORT
#elif defined(_MSC_VER)
# define FNM_EXTERNAL_API __declspec(dllexport)
# define FNM_EXPORT __declspec(dllexport)
#endif

enum MyMatlabEnum {
  First = 0,
  Last = 1
};

# include <fnm/fnm_types.h>
# include <fnm/if_fnm.h>

#endif
