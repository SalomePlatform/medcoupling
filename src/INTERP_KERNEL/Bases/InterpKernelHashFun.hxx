// Copyright (C) 2001-2016  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

/*
 * Copyright (c) 1996-1998
 * Silicon Graphics Computer Systems, Inc.
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Silicon Graphics makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 *
 * Copyright (c) 1994
 * Hewlett-Packard Company
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Hewlett-Packard Company makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 */
#ifndef __INTERPKERNELHASHFUN_HXX__
#define __INTERPKERNELHASHFUN_HXX__

#include <cstddef>

namespace INTERP_KERNEL
{
  template<class _Key>
  struct hash { };

  inline std::size_t __stl_hash_string(const char* __s)
  {
    unsigned long __h = 0;
    for ( ; *__s; ++__s)
      __h = 5 * __h + *__s;
    return std::size_t(__h);
  }

  template<>
  struct hash<char*>
  {
    std::size_t operator()(const char* __s) const
    { return __stl_hash_string(__s); }
  };

  template<>
  struct hash<const char*>
  {
    std::size_t operator()(const char* __s) const
    { return __stl_hash_string(__s); }
  };

  template<>
  struct hash<char>
  { 
    std::size_t operator()(char __x) const { return __x; }
  };

  template<>
  struct hash<unsigned char>
  { 
    std::size_t operator()(unsigned char __x) const { return __x; }
  };

  template<>
  struct hash<signed char>
  {
    std::size_t operator()(unsigned char __x) const { return __x; }
  };

  template<>
  struct hash<short>
  {
    std::size_t operator()(short __x) const { return __x; }
  };

  template<>
  struct hash<unsigned short>
  {
    std::size_t operator()(unsigned short __x) const { return __x; }
  };

  template<>
  struct hash<int>
  { 
    std::size_t operator()(int __x) const { return __x; }
  };

  template<>
  struct hash<unsigned int>
  { 
    std::size_t operator()(unsigned int __x) const { return __x; }
  };

  template<>
  struct hash<long>
  {
    std::size_t operator()(long __x) const { return __x; }
  };

  template<>
  struct hash<unsigned long>
  {
    std::size_t operator()(unsigned long __x) const { return __x; }
  };
}

#endif
