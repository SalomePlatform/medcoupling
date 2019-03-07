// Copyright (C) 2017-2019  CEA/DEN, EDF R&D
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

#ifndef __MEDCOUPLING_MCTYPE_HXX__
#define __MEDCOUPLING_MCTYPE_HXX__

#include <cstdint>
#include <stddef.h>
#include <cstddef>

namespace MEDCoupling
{
  using Int64 = std::int64_t;
  using Int32 = std::int32_t;
#ifndef MEDCOUPLING_USE_64BIT_IDS
  using mcIdType = std::int32_t;
#else
  using mcIdType = std::int64_t;
#endif
  inline mcIdType ToIdType(std::size_t val) { return mcIdType(val); }
}

#define DataArrayInt DataArrayInt32
#define DataArrayIdType DataArrayInt32

#define DataArrayIntIterator DataArrayInt32Iterator

#endif
