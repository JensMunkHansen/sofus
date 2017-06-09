/**
 * @file   globals.hpp
 * @author Jens Munk Hansen <jmh@jmhlaptop>
 * @date   Fri Mar 31 20:45:27 2017
 *
 * @brief  Templated singletons
 *
 *
 */
#pragma once

namespace sps {
  /*
    Template <struct V> is not allowed so this gives as warning, when V is a struct
  */
#ifdef _MSC_VER
# pragma warning(push)
# pragma warning(disable : 4099)
#endif

  template <typename T, template <typename> class V> class globalstruct
  {
  public:

    template <typename U> using MyTemplate = V<U>;

    static class V<T> *pVar;
  private:
    globalstruct() = delete;
    globalstruct(const globalstruct& rhs) = delete;
    void operator=(const globalstruct& rhs) = delete;
  };

#ifdef _MSC_VER
# pragma warning(pop)
#endif

  template <typename T> class global {
  public:
    static T *pVar;

  private:
    global() = delete;
    global(const global& rhs) = delete;
    void operator=(const global& rhs) = delete;
  };
}

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
