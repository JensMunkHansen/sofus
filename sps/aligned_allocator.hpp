/**
 * @file   aligned_allocator.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Wed Jul 22 23:18:58 2011
 *
 * @brief  Aligned allocator for STL types
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

#include <cstddef>           // Required for size_t and ptrdiff_t and NULL
#include <new>               // Required for placement new and std::bad_alloc
#include <stdexcept>         // Required for std::length_error

#include <cstdlib>           // Required for malloc() and free()

#ifdef _WIN32
# include <malloc.h>         // Required for _mm_malloc() and _mm_free()
#else
# include <mm_malloc.h>      // Required for _mm_malloc() and _mm_free()
#endif

#ifndef UNUSED
# define UNUSED(p) ((void)(p))
#endif

/*! Aligned allocator for STL containers. */
template <typename T, std::size_t Alignment = 16> class aligned_allocator {
public:

  /// <summary>   STL standard aliases. . </summary>
  typedef T * pointer;
  typedef const T * const_pointer;
  typedef T& reference;
  typedef const T& const_reference;
  typedef T value_type;
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;

  /**
   * Address of type T.
   *
   * @param r [in] The object of type T.
   *
   * @return the address
   */
  T* address(T& r) const
  {
    return &r;
  }

  /**
   * Const address of type T.
   *
   * @param s The const object of type T.
   *
   * @return the address
   */
  const T* address(const T& s) const
  {
    return &s;
  }

  /**
   * Gets the maximum size.
   *
   *
   * @return the size
   */
  size_t max_size() const
  {
    return (static_cast<size_t>(0) - static_cast<size_t>(1)) / sizeof(T);
  }

  /*! Structure for rebinding */
  template <typename U> struct rebind {
    typedef aligned_allocator<U,Alignment> other;
  };

  /**
   * Not equal to
   *
   * @param other
   *
   * @return The other
   */
  bool operator!=(const aligned_allocator& other) const
  {
    return !(*this == other);
  }

  /**
   * Construct
   *
   * @param p [in] If non-null, the T* p is the address used for construction.
   * @param t The object to displace.
   */
  void construct(T* const p, const T& t) const
  {
    void* const pv = static_cast<void*>(p);
    new (pv) T(t);
  }

  /**
   * Destroys the given p.
   *
   * @param p [in] If non-null, the T* const to destroy.
   */
  void destroy(T* const p) const;

  /**
   * Equality operator.
   *
   * @param other
   *
   * @return true for stateless allocators.
   */
  bool operator==(const aligned_allocator& other) const
  {
    return true;
  }

  /**
   * Default constructor - empty for stateless allocators.
   *
   */
  aligned_allocator() { }

  /**
   * Default copy constructor - empty for stateless allocators.
   *
   */
  aligned_allocator(const aligned_allocator&) { }

  /**
   * Default rebinding constructor - empty for stateless allocators.
   *
   * @param other
   */
  template <typename U> aligned_allocator(const aligned_allocator<U,Alignment>& other) { }

  /**
   * Destructor
   *
   */
  ~aligned_allocator() { }

  /// <summary>   Allocates memory. </summary>
  /// <param name="n">    The. </param>
  /// <returns>

  /**
   * Allocates memory
   *
   * @exception std::length_error      Thrown when length error.
   * @exception std::bad_alloc         Thrown when bad allocate.
   * @param n The size
   *
   * @return null if it fails, else a reference pointer.
   */
  T* allocate(const size_t n) const
  {
    // The return value of allocate(0) is unspecified.
    // aligned_allocator returns NULL in order to avoid depending
    // on malloc(0)'s implementation-defined behavior
    // (the implementation can define malloc(0) to return NULL,
    // in which case the bad_alloc check below would fire).
    // All allocators can return NULL in this case.
    if (n == 0) {
      return NULL;
    }

    // All allocators should contain an integer overflow check.
    // The Standardization Committee recommends that std::length_error
    // be thrown in the case of integer overflow.
    if (n > max_size()) {
      throw std::length_error("aligned_allocator<T>::allocate() - Integer overflow.");
    }

    // aligned_allocator wraps _mm_malloc().
    void* const pv = _mm_malloc(n * sizeof(T),Alignment);

    // Allocators should throw std::bad_alloc in the case of memory allocation failure.
    if (pv == NULL) {
      throw std::bad_alloc();
    }
    return static_cast<T*>(pv);
  }

  /**
   * Deallocates the memory
   *
   * @param p [in] If non-null, the T* p is the address used for construction.
   * @param n The length of the buffer.
   */
  void deallocate(T * const p, const size_t n) const
  {
    UNUSED(n);
    _mm_free(p);
  }

  /**
   * The allocator ignores hints, so the same as allocate.
   *
   * @param n
   *
   * @return null if it fails, else.
   */
  template <typename U> T * allocate(const size_t n, const U * /* const hint */) const
  {
    return allocate(n);
  }

  // Allocators are not required to be assignable, so
  // all allocators should have a private unimplemented
  // assignment operator. Note that this will trigger the
  // off-by-default (enabled under /Wall) warning C4626
  // "assignment operator could not be generated because a
  // base class assignment operator is inaccessible" within
  // the STL headers, but that warning is useless.

private:

  /**
   * Private unimplemented assignment operator.
   *
   * @param other
   *
   * @return A shallow copy of this object.
   */
  aligned_allocator& operator=(const aligned_allocator& other);

};

// A compiler bug causes it to believe that p->~T() doesn't reference p.
#ifdef _MSC_VER
# pragma warning(push)
# pragma warning(disable: 4100) // Unreferenced formal parameter
#endif

/**
 * Destroys the given p. The definition of destroy() must be the same for all allocators.
 *
 * @param p [in] If non-null, the T * const to destroy.
 */
template <typename T, std::size_t Alignment>
void aligned_allocator<T,Alignment>::destroy(T * const p) const
{
  p->~T();
}

#ifdef _MSC_VER
# pragma warning(pop)
#endif

/* Local variables: */
/* indent-tab-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* indent-tabs-mode: nil */
/* End: */
