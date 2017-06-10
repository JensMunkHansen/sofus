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
