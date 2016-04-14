// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
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

#ifndef __INTERPKERNELSTLEXT_HXX__
#define __INTERPKERNELSTLEXT_HXX__

#include <functional>

namespace INTERP_KERNEL
{
  namespace STLEXT
  {
    template<typename _Pair>
    struct Select1st : public std::unary_function<_Pair, typename _Pair::first_type>
    {
      typename _Pair::first_type& operator()(_Pair& __x) const { return __x.first; }
      const typename _Pair::first_type&operator()(const _Pair& __x) const { return __x.first; }
    };

    template<typename _T1, typename _T2>
    inline void Construct(_T1* __p, const _T2& __value__) { ::new(static_cast<void*>(__p)) _T1(__value__); }

    template<typename _Tp> inline void Destroy(_Tp* __pointer) { __pointer->~_Tp(); }
  }
}

#endif
