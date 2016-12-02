/**
 * @file   multi_malloc.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Fri May 23 11:42:25 2014
 *
 * @brief  Allocation routine for multi-dimensional arrays
 *
 */
#pragma once

#include <cstddef>
#include <cstdlib>
#include <cstdarg>
#include <cassert>

#include <type_traits>

#include <sps/math.h>      // is_power_of_two

// TEST
#include <sps/mm_malloc.h>
#include <sps/memory>

#ifndef DOXYGEN_SHOULD_SKIP_THIS

template<std::size_t Size, std::size_t Alignment> struct check_size_divides_alignment {

  static_assert(Alignment % Size == 0, "Size must be a divisor of alignment");
#ifdef __GNUG__
  static_assert(is_power_of_two(Alignment), "Alignment must be power of two");
#endif
};

#endif

#if defined(_MSC_VER) && (_MSC_VER < 1800)
template<typename T, size_t Alignment>
#else
template<typename T, size_t Alignment=alignof(max_align_t)>
#endif
void* __mm_multi_malloc(size_t d, va_list ap)
{
  char* tree;

  size_t max,             /* size of array to be declared */
         *q;                   /* pointer to dimension list */
  char **r,               /* pointer to beginning of the array of the
                           * pointers for a dimension */
       **s1, *t;             /* base pointer to beginning of first array */
  size_t i, j;            /* loop counters */
  size_t *d1;             /* dimension list */

  size_t sa;              /* elements per alignment */

  //va_start(ap,d);
  d1 = (size_t *) malloc(d*sizeof(size_t));

  for(i=0; i<d; i++)
    d1[i] = va_arg(ap,size_t);

  sa = Alignment / sizeof(T);

  r = &tree;
  q = d1;                 /* first dimension */

  if (d==1) {
    max = sa * (((*q)+(sa-1)) / sa);
    free(d1);
    //va_end(ap);
    return _mm_malloc(max*sizeof(T), Alignment);
  }

  max = 1;
  for (i = 0; i < d - 1; i++, q++) {
    /* for each of the dimensions
                                             * but the last */
    max *= (*q);
    r[0]=(char *)_mm_malloc(max * sizeof(char **),Alignment);
    r = (char **) r[0];     /* step through to beginning of next
                             * dimension array */
  }

  /* grab actual array memory */
  max *= sizeof(T) * (size_t) (sa * ( ( (*q) + (sa-1) ) / sa ) );
  r[0] = (char *)_mm_malloc(max * sizeof(char),Alignment);

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
    max *= (*q);    /* number of elements in this dimension */
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
    *s1++ = (t += sizeof(T) **(q + 1));

  //va_end(ap);
  free(d1);

  return((void *)tree);              /* return base pointer */
}

/**
 * _mm_multi_malloc
 *
 * Allocate multi-dimensional array, where the inner-most dimension
 * is aligned by Alignment bytes (default=16). The size of the
 * individual dimensions are given as variadic arguments.
 *
 * Example:
 *
 * float** f = _mm_multi_malloc<float>(2,100,100);
 * f[99][99] = 1.0f;
 * _mm_multi_free<float>(f,2);
 *
 * @param d Number of dimensions
 *
 * @return
 */
#if defined(_MSC_VER) && (_MSC_VER < 1800)
template<typename T, size_t Alignment>
#else
template<typename T, size_t Alignment=alignof(max_align_t)>
#endif
void* _mm_multi_malloc(size_t d, ...)
{

  check_size_divides_alignment<sizeof(T), Alignment>();

#ifndef __STRICT_ANSI__
  static_assert(std::is_scalar<T>::value,"Type T must be scalar");
#endif

  va_list ap;             /* varargs list traverser */
  va_start(ap, d);
  void* retval = __mm_multi_malloc<T,Alignment>(d, ap);
  va_end(ap);
  return retval;
}



/**
 *  _mm_multi_free
 *
 * Free multi-dimensional arrays allocated using _mm_multi_malloc
 *
 * @param r The pointer
 * @param d Number of dimensions
 */
#if defined(_MSC_VER) && (_MSC_VER < 1800)
template<typename T, size_t Alignment>
#else
template<typename T, size_t Alignment=alignof(max_align_t)>//std::alignment_of<max_align_t>::value>
#endif
void _mm_multi_free(void *r, size_t d)
{

  check_size_divides_alignment<sizeof(T), Alignment>();

  void **p;
  void *next=NULL;
  size_t i;

  for (p = (void **)r, i = 0; i < d; p = (void **) next,i++)
    if (p != NULL) {
      next = *p;
      _mm_free(p); // Removed a & before p
      p = NULL;
    }
}

#ifndef _MSC_VER
/**
 * _mm_multi_malloc_nc
 *
 * Allocate multi-dimensional arrays for SIMD usage. All sub-
 * dimensions are aligned at the cost of an additional memory
 * usage. Be careful, the actual content of a matrix may not
 * be contigously accessible using the address of the first
 * element
 *
 * Example:
 *
 * float** f = _mm_multi_malloc_nc<float>(2,100,100);
 * f[99][99] = 1.0f;
 * _mm_multi_free_nc<float>(f,2);
 *
 * @param d Number of dimensions
 *
 * @return
 *
 * TODO: Change to only inner-most dimension is padded
 */
template<typename T, size_t Alignment=16>
void* _mm_multi_malloc_nc(size_t d, ...)
{

  //  static_assert(Alignment % sizeof(T) == 0, "size must be a divisor of alignment");
  check_size_divides_alignment<sizeof(T), Alignment>();

  char* tree;

  va_list ap;             /* varargs list traverser */
  size_t max,             /* size of array to be declared */
         *q;                   /* pointer to dimension list */
  char **r,               /* pointer to beginning of the array of the
                           * pointers for a dimension */
       **s1, *t;             /* base pointer to beginning of first array */
  size_t i, j;            /* loop counters */
  size_t *d1;             /* dimension list */

  size_t sa;              /* elements per alignment */

  va_start(ap,d);
  d1 = (size_t *) malloc(d*sizeof(size_t));

  for(i=0; i<d; i++)
    d1[i] = va_arg(ap,size_t);

  sa = Alignment / sizeof(T);

  r = &tree;
  q = d1;                 /* first dimension */

  if (d==1) {
    max = sa * (( sizeof(T)*(*q) + (sa-1))/sa);
    va_end(ap);
    return _mm_malloc(max, Alignment);
  }

  max = 1;
  for (i = 0; i < d - 1; i++, q++) {
    /* for each of the dimensions
     * but the last */
    max *= sa * (((*q)+(sa-1))/sa);
    //
    r[0]=(char *)_mm_malloc(max * sizeof(char **),Alignment);
    r = (char **) r[0];     /* step through to beginning of next
                             * dimension array */
  }
  max *= sizeof(T) * (size_t) (sa * (((*q)+(sa-1))/sa));        /* grab actual array memory */
  r[0] = (char *)_mm_malloc(max * sizeof(char),Alignment);

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
    max *= sa * (((*q)+(sa-1))/sa);    /* number of elements in this dimension */
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

      *s1 = (t += sizeof (char **) *  (sa * (((*q+1)+(sa-1))/sa)) );
      s1++;
    }
    r = (char **) r[0];     /* step through to begining of next
                             * dimension array */
  }
  max *= sa * (((*q)+(sa-1))/sa);              /* max is total number of elements in the
                             * last pointer array */

  /* same as previous loop, but different size factor */
  for (j = 1, s1 = r + 1, t = r[0]; j < max; j++)
    *s1++ = (t += sizeof(T) * (sa * (((*q+1)+(sa-1))/sa)) );

  va_end(ap);
  free(d1);

  return((void *)tree);              /* return base pointer */
}

/**
 *  _mm_multi_free_nc
 *
 * Free multi-dimensional arrays allocated using _mm_multi_malloc_nc
 *
 * @param r The pointer
 * @param d Number of dimensions
 */
template<typename T, size_t Alignment=16>
void _mm_multi_free_nc(void *r, size_t d)
{

  check_size_divides_alignment<sizeof(T), Alignment>();

  void **p;
  void *next=NULL;
  size_t i;

  for (p = (void **)r, i = 0; i < d; p = (void **) next,i++)
    if (p != NULL) {
      next = *p;
      _mm_free(p); // Remove a & before p
    }
}
#endif

// Deleted multi_arrays

namespace sps {
#if (defined(__GNUG__) && __cplusplus >= 201103L)

  namespace test {
    // Works
    template<typename T>
    using deleted_unique_multi_array = std::unique_ptr<T[], std::function<void(T*)> >;

    template<typename T, size_t d>
    sps::test::deleted_unique_multi_array<T> deleted_aligned_multi_array(size_t _d, ...)
    {
      va_list ap;
      va_start(ap, _d);
      auto retval = deleted_unique_multi_array<T>((T*)__mm_multi_malloc<T>(d, ap), [](T* f)->void {_mm_multi_free<T>((void*)f, d);});
      va_end(ap);
      return retval;
    }

    // This works
#if 0
    std::unique_ptr<float[][10][10], std::function< void(float (*)[10][10])> >  aa =
      std::unique_ptr<float[][10][10], std::function< void(float(*)[10][10])> >(new (float[2][10][10]), [](float (*f)[10][10])->void {delete[] f;});
#endif

    // Not working
    template<typename T, size_t a, size_t b>
    using deleted_unique_multi_array2 = std::unique_ptr<T[][b], std::function< void(T(*)[b])> >;

    template<typename T, size_t a, size_t b> // size_t... and using const size_t ndims = sizeof...(size_t)
    sps::test::deleted_unique_multi_array2<T,a,b> deleted_aligned_multi_array2()
    {
      return deleted_unique_multi_array2<T,a,b>( (T (*)[b]) _mm_multi_malloc<T>(2, a, b),
             [](T (*f)[b])->void {_mm_multi_free<T>( (void*)&(f[0][0]), 2);});
    }

    template<typename T, size_t a, size_t b, size_t c>
    using deleted_unique_multi_array3 = std::unique_ptr<T[][b][c], std::function< void(T (*)[b][c])> >;

    template<typename T, size_t a, size_t b,size_t c> // size_t... and using const size_t ndims = sizeof...(size_t)
    sps::test::deleted_unique_multi_array3<T,a,b,c> deleted_aligned_multi_array3()
    {
      return deleted_unique_multi_array3<T,a,b,c>( (T (*)[b][c]) _mm_multi_malloc<T>(3, a, b, c),
             [](T (*f)[b][c])->void {_mm_multi_free<T>(&(f[0][0][0]), 3);});
    }
  }
#endif
}

/*
#define RETURNS(x) ->decltype(x) { return (x); }

template<typename ...Args>
auto get_last( Args&&... args )
  RETURNS( std::get< sizeof...(Args)-1 >( std::tie(std::forward<Args>(args)...) ) )
we can then use this in another function:

template<typename ...Args>
void foo( Args&&... args ) {
  auto&& last = get_last(std::forward<Args>(args)...);
}
 */


// Wrappers needed because Microsoft do not support default template arguments for functions
#if defined(_MSC_VER) && (_MSC_VER < 1800)
template <typename T>
void* _mm_multi_malloc(size_t d, size_t dim0)
{
  return _mm_multi_malloc<T,16>(d, dim0);
}
template <typename T>
void* _mm_multi_malloc(size_t d, size_t dim0, size_t dim1)
{
  return _mm_multi_malloc<T,16>(d, dim0, dim1);
}
template <typename T>
void _mm_multi_free(void* r, size_t d)
{
  return _mm_multi_free<T,16>(r,d);
}
#endif

/*
You could std::bind your deleter's second argument before passing it as the deleter:

auto deleter = std::bind(myDeleter, std::placeholders::_1, 5);
std::shared_ptr<A> myA(a, deleter);
Alternatively, your deleter could be a functor that takes the int through its constructor:

struct myDeleter
{
  myDeleter(int);
  void operator()(A*);
};

myDeleter deleter(5);
std::shared_ptr<A> myA(a, deleter);
Alternatively you could use a lambda expression:

std::shared_ptr<A> myA(a, [](A* a){ myDeleter(a, 5); });

std::unique_ptr<int*, std::function<void(int**)>> x(
    new int*[10](),
    [](int** x) {
        std::for_each(x, x + 10, std::default_delete<int[]>());
        delete[] x;
    }
);

The unique_ptr declaration takes care of allocating the row dimension
of the array. The trailing () in new int*[10]() ensures that each
column pointer is initialized to nullptr.

A for loop then allocates the column arrays:

for (size_t row = 0; row < 10; ++row) {
    (x.get())[row] = new int[5];
}


Anyway, what you ask can be done, more or less:

int (*c)[2] = (int(*)[2])new int[2];
But a typedef will make it easier:

typedef int ai[2];
ai *c = (ai*)new int[2];
And to be safe, the delete should be done using the original type:

delete [](int*)c;



*/


/* Local variables: */
/* indent-tab-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* indent-tabs-mode: nil */
/* End: */
