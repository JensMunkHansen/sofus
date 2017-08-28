/**
 * @file   multi_malloc.h
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Fri May 5 19:00:21 2007
 *
 * @brief  Allocation of multi-dimensional C arrays using row pointers
 *
 *
 */
#pragma once

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

#include <sps/cenv.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>
#include <string.h>

STATIC_INLINE_BEGIN void* multi_malloc(size_t s, size_t d, ...) STATIC_INLINE_END;

STATIC_INLINE_BEGIN void multi_free(void *r, size_t d) STATIC_INLINE_END;


/**
 * Allocate multi-dimensional array and establish row pointers
 *
 * @param s size of each element
 * @param d number of dimension
 *
 * @return
 */
STATIC_INLINE_BEGIN void* multi_malloc(size_t s, size_t d, ...)
{

  char* tree;

  va_list ap;             /* varargs list traverser */
  size_t max,             /* size of array to be declared */
         *q;                   /* pointer to dimension list */
  char **r,               /* pointer to beginning of the array of the
                           * pointers for a dimension */
       **s1, *t;             /* base pointer to beginning of first array */
  size_t i, j;            /* loop counters */
  size_t *d1;             /* dimension list */

  va_start(ap,d);
  d1 = (size_t *) malloc(d*sizeof(size_t));

  for(i=0; i<d; i++)
    d1[i] = va_arg(ap,size_t);

  r = &tree;
  q = d1;                 /* first dimension */

  if (d==1) {
    max = s * (*q);
    va_end(ap);
    free(d1);
    return malloc(max);
  }

  max = 1;
  for (i = 0; i < d - 1; i++, q++) {
    /* for each of the dimensions
                                             * but the last */
    max *= (*q);
    r[0]=(char *)malloc(max * sizeof(char **));
    r = (char **) r[0];     /* step through to beginning of next
                             * dimension array */
  }
  max *= s * (size_t) (*q);        /* grab actual array memory */
  r[0] = (char *)malloc(max * sizeof(char));

  /*
   * r is now set to point to the beginning of each array so that we can
   * use it to scan down each array rather than having to go across and
   * then down
   */
  r = (char **) tree;     /* back to the beginning of list of arrays */
  q = d1;                 /* back to the first dimension */
  max = 1;
  for (i = 0; i < d - 2; i++, q++) {
    /* we deal with the last
                                           * array of pointers later on */
    max *= (*q);                        /* number of elements in this dimension */
    for (j=1, s1=r+1, t=r[0]; j<max; j++) {
      /* scans down array for
                                               * first and subsequent
                                               * elements */

      /*  modify each of the pointers so that it points to
       * the correct position (sub-array) of the next
       * dimension array. s1 is the current position in the
       * current array. t is the current position in the
       * next array. t is incremented before s1 is, but it
       * starts off one behind. *(q+1) is the dimension of
       * the next array. */

      *s1 = (t += sizeof (char **) **(q + 1));
      s1++;
    }
    r = (char **) r[0];     /* step through to begining of next
                             * dimension array */
  }
  max *= (*q);              /* max is total number of elements in the
                             * last pointer array */

  /* same as previous loop, but different size factor */
  for (j = 1, s1 = r + 1, t = r[0]; j < max; j++)
    *s1++ = (t += s **(q + 1));

  va_end(ap);
  free(d1);

  return((void *)tree);              /* return base pointer */
}

/**
 * Free multi-dimensional array and corresponding row pointers
 *
 * @param r data
 * @param d number of dimensions
 */
STATIC_INLINE_BEGIN void multi_free(void *r, size_t d)
{

  void **p;
  void *next=NULL;
  size_t i;

  for (p = (void **)r, i = 0; i < d; p = (void **) next,i++)
    if (p != NULL) {
      next = *p;
      free(p);
      p = NULL;
    }
}

/* Local variables: */
/* indent-tab-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* indent-tabs-mode: nil */
/* End: */
