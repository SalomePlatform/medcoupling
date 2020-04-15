// Copyright (C) 2017-2020  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (EDF R&D)

#ifndef __IK_MCIDTYPE_HXX__
#define __IK_MCIDTYPE_HXX__

#include <cstdint>
#include <stddef.h>
#include <cstddef>

#ifndef MEDCOUPLING_USE_64BIT_IDS

typedef std::int32_t mcIdType;

#else

typedef std::int64_t mcIdType;

#endif

template <class T> inline mcIdType ToIdType(T val)
{
  return static_cast<mcIdType>(val);
}
template <class T> inline T FromIdType(mcIdType val)
{
  return static_cast<T>(val);
}


#endif
