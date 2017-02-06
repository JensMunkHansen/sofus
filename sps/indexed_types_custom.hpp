#pragma once

namespace sps {
  class bad_indexed_auto_ptr : public std::exception {
  public:
    char const*
    what() const noexcept
    {
      return "sps::bad_indexed_auto_ptr";
    }
  };
// Substitute for bad_weak_ptr object in the case of -fno-exceptions.
  inline void
  __throw_bad_indexed_auto_ptr()
  {
#if __EXCEPTIONS
    throw bad_indexed_auto_ptr();
#else
    __builtin_abort();
#endif
  }

};
