/*
 *  This file is part of the MR utility library.
 *
 *  This code is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This code is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this code; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/* Copyright (C) 2019-2020 Max-Planck-Society
   Author: Martin Reinecke */

#ifndef DUCC0_ERROR_HANDLING_H
#define DUCC0_ERROR_HANDLING_H

#include <sstream>
#include <exception>
#include <stdexcept>

#include "ducc0/infra/useful_macros.h"

namespace ducc0 {

namespace detail_error_handling {

#if defined (__GNUC__)
#define DUCC0_ERROR_HANDLING_LOC_ ::ducc0::detail_error_handling::CodeLocation(__FILE__, __LINE__, __PRETTY_FUNCTION__)
#else
#define DUCC0_ERROR_HANDLING_LOC_ ::ducc0::detail_error_handling::CodeLocation(__FILE__, __LINE__)
#endif

// to be replaced with std::source_location once generally available
class CodeLocation
  {
  private:
    const char *file, *func;
    int line;

  public:
    CodeLocation(const char *file_, int line_, const char *func_=nullptr)
      : file(file_), func(func_), line(line_) {}

    inline ::std::ostream &print(::std::ostream &os) const
      {
      os << "\n" << file <<  ": " <<  line;
      if (func) os << " (" << func << ")";
      os << ":\n";
      return os;
      }
  };

inline ::std::ostream &operator<<(::std::ostream &os, const CodeLocation &loc)
  { return loc.print(os); }

void streamDump__(std::ostream&)
  { }
template<typename T, typename ...Args>
void streamDump__(std::ostream &os, T&& arg0, Args&&... args)
  { os << std::forward<T>(arg0); streamDump__(os, std::forward<Args>(args)...); }
template<typename ...Args>
[[noreturn]] DUCC0_NOINLINE void fail__(Args&&... args)
  {
    std::ostringstream msg;
    ducc0::detail_error_handling::streamDump__(msg, std::forward<Args>(args)...);
    throw std::runtime_error(msg.str());
  }

#define MR_fail(...) \
  { \
    /*::ducc0::detail_error_handling::fail__(DUCC0_ERROR_HANDLING_LOC_, "\n", ##__VA_ARGS__, "\n");*/ \
  }

#define MR_assert(cond,...) \
  { \
    if (cond); \
    else { MR_fail("Assertion failure\n", ##__VA_ARGS__); } \
  }

}}

#endif
