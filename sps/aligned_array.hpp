/**
 * @file   aligned_array.hpp
 * @author  <jens.munk.hansen@gmail.com>
 * @date   Tue Jul 14 14:43:20 2015
 *
 * @brief The following code declares class aligned_array, an STL-like
 * container (as wrapper) for arrays of constant size.  The memory is
 * aligned on a 16-byte boundary
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

#include <sps/cenv.h>
#include <sps/sps_export.h>

#include <cstddef>
#include <stdexcept>

#include <iterator>
#include <algorithm>

namespace std {

  template <typename T, std::size_t N>
  class SPS_EXPORT aligned_array {
  public:
    ALIGN16_BEGIN T elems[N] ALIGN16_END;    /// fixed-size array of elements of type T
    //__declspec( align(4*sizeof(T)) ) T elems[N];    /// fixed-size array of elements of type T

  public:
    /** @name Type definitions
     *
     */
    ///@{
    typedef T              value_type;
    typedef T*             iterator;
    typedef const T*       const_iterator;
    typedef T&             reference;
    typedef const T&       const_reference;
    typedef std::size_t    size_type;
    typedef std::ptrdiff_t difference_type;
    ///@}

    /** @name Iterator support
     *
     */
    ///@{

    iterator begin()
    {
      return elems;
    }
    const_iterator begin() const
    {
      return elems;
    }
    iterator end()
    {
      return elems+N;
    }
    const_iterator end() const
    {
      return elems+N;
    }

    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    reverse_iterator rbegin()
    {
      return reverse_iterator(end());
    }
    const_reverse_iterator rbegin() const
    {
      return const_reverse_iterator(end());
    }
    reverse_iterator rend()
    {
      return reverse_iterator(begin());
    }
    const_reverse_iterator rend() const
    {
      return const_reverse_iterator(begin());
    }
    ///@}

    reference operator[](size_type i)
    {
      return elems[i];
    }
    const_reference operator[](size_type i) const
    {
      return elems[i];
    }

    /** @name at() with range check
     *
     */
    ///@{
    reference at(size_type i)
    {
      rangecheck(i);
      return elems[i];
    }
    const_reference at(size_type i) const
    {
      rangecheck(i);
      return elems[i];
    }
    ///@}

    /** @name front() and back()
     *
     */
    ///@{
    reference front()
    {
      return elems[0];
    }
    const_reference front() const
    {
      return elems[0];
    }
    reference back()
    {
      return elems[N-1];
    }
    const_reference back() const
    {
      return elems[N-1];
    }
    ///@}

    // size is constant
    static size_type size()
    {
      return N;
    }
    static bool empty()
    {
      return false;
    }
    static size_type max_size()
    {
      return N;
    }
    enum { static_size = N };

    /**
     * swap (note: linear complexity in N, constant for given instantiation)
     *
     * @param T
     * @param y
     */
    void swap (aligned_array<T,N>& y)
    {
#if defined(_MSC_VER) && defined(NDEBUG)
      // Gives no warnings, but linker errors in debug mode
      std::swap_ranges(begin(),end(), stdext::make_checked_array_iterator(y.begin(), std::distance(begin(),end())));
#else
      std::swap_ranges(begin(),end(), y.begin());
#endif
    }

    /**
     * Direct access to data (read-only)
     *
     *
     * @return
     */
    const T* data() const
    {
      return elems;
    }

    /**
     * Use array as C array (direct read/write access to data)
     *
     *
     * @return
     */
    T* data()
    {
      return elems;
    }

    /**
     * assignment with type conversion
     *
     * @tparam T2
     * @param rhs
     *
     * @return
     */
    template <typename T2>
    aligned_array<T,N>& operator= (const aligned_array<T2,N>& rhs)
    {
      std::copy(rhs.begin(),rhs.end(), begin());
      return *this;
    }

    /**
     * assign one value to all elements
     *
     * @param value
     */
    void assign (const T& value)
    {
#if defined(_MSC_VER) && defined(NDEBUG)
      std::fill_n(stdext::make_checked_array_iterator(begin(),std::distance(begin(),end())),size(),value);
#else
      std::fill_n(begin(),size(),value);
#endif
    }

  private:

    /**
     * check range (may be private because it is static)
     *
     * @param i
     */
    static void rangecheck (size_type i)
    {
      if (i >= size()) {
        throw std::out_of_range("aligned_array<>: index out of range");
      }
    }

  };

  template <typename T>
  class aligned_array<T,0> {
  public:
    char c;  // to ensure different array intances return unique values for begin/end

  public:
    // type definitions
    typedef T              value_type;
    typedef T*             iterator;
    typedef const T*       const_iterator;
    typedef T&             reference;
    typedef const T&       const_reference;
    typedef std::size_t    size_type;
    typedef std::ptrdiff_t difference_type;

    // iterator support
    iterator begin()
    {
      return reinterpret_cast< iterator >( &c );
    }
    const_iterator begin() const
    {
      return reinterpret_cast< const_iterator >( &c );
    }
    iterator end()
    {
      return reinterpret_cast< iterator >( &c );
    }
    const_iterator end() const
    {
      return reinterpret_cast< const_iterator >( &c );
    }

    // reverse iterator support
#if 1
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
#else
    // workaround for broken reverse_iterator implementations
    typedef std::reverse_iterator<iterator,T> reverse_iterator;
    typedef std::reverse_iterator<const_iterator,T> const_reverse_iterator;
#endif

    reverse_iterator rbegin()
    {
      return reverse_iterator(end());
    }
    const_reverse_iterator rbegin() const
    {
      return const_reverse_iterator(end());
    }
    reverse_iterator rend()
    {
      return reverse_iterator(begin());
    }
    const_reverse_iterator rend() const
    {
      return const_reverse_iterator(begin());
    }

    // at() with range check
    reference at(size_type i)
    {
      throw std::out_of_range("aligned_array<0>: index out of range");
    }
    const_reference at(size_type i) const
    {
      throw std::out_of_range("<0>: index out of range");
    }

    // size is constant
    static size_type size()
    {
      return 0;
    }
    static bool empty()
    {
      return true;
    }
    static size_type max_size()
    {
      return 0;
    }
    enum { static_size = 0 };

    // swap
    void swap (aligned_array<T,0>& y)
    {
      //  could swap value of c, but value is not part of documented array state
    }

    // direct access to data
    const T* data() const
    {
      return NULL;
    }
    T* data()
    {
      return NULL;
    }

    // assignment with type conversion
    template < typename T2 >
    aligned_array< T,0 >& operator= (const aligned_array< T2, 0>& rhs)
    {
      return *this;
    }

    //  Calling these operations are undefined behaviour for 0-size arrays,
    //  but Library TR1 requires their presence.
    // operator[]
    reference operator[](size_type i)
    {
      makes_no_sense();
    }
    const_reference operator[](size_type i) const
    {
      makes_no_sense();
    }

    // front() and back()
    reference front()
    {
      makes_no_sense();
    }
    const_reference front() const
    {
      makes_no_sense();
    }
    reference back()
    {
      makes_no_sense();
    }
    const_reference back() const
    {
      makes_no_sense();
    }

  private:

    // helper for operations that have undefined behaviour for 0-size arrays,
    //  assert( false ); added to make lack of support clear
    static void makes_no_sense ()
    {
      assert(true);
      throw std::out_of_range("aligned_array<0>: index out of range");
    }
  };

// comparisons
  template<class T, std::size_t N>
  bool operator== (const aligned_array<T,N>& x, const aligned_array<T,N>& y)
  {
    return std::equal(x.begin(), x.end(), y.begin());
  }
  template<class T, std::size_t N>
  bool operator< (const aligned_array<T,N>& x, const aligned_array<T,N>& y)
  {
    return std::lexicographical_compare(x.begin(),x.end(),y.begin(),y.end());
  }
  template<class T, std::size_t N>
  bool operator!= (const aligned_array<T,N>& x, const aligned_array<T,N>& y)
  {
    return !(x==y);
  }
  template<class T, std::size_t N>
  bool operator> (const aligned_array<T,N>& x, const aligned_array<T,N>& y)
  {
    return y<x;
  }
  template<class T, std::size_t N>
  bool operator<= (const aligned_array<T,N>& x, const aligned_array<T,N>& y)
  {
    return !(y<x);
  }
  template<class T, std::size_t N>
  bool operator>= (const aligned_array<T,N>& x, const aligned_array<T,N>& y)
  {
    return !(x<y);
  }

// global swap()
  template<class T, std::size_t N>
  inline void swap (aligned_array<T,N>& x, aligned_array<T,N>& y)
  {
    x.swap(y);
  }

  template class SPS_EXPORT aligned_array<float,4>;
  template class SPS_EXPORT aligned_array<double,4>;
}

/* Local variables: */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* indent-tabs-mode: nil */
/* End: */

