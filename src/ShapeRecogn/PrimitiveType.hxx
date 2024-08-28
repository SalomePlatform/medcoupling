// Copyright (C) 2024  CEA, EDF
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

#pragma once

#include <vector>
#include <string>
#include <cstdint>

namespace MEDCoupling
{
  enum class PrimitiveType : std::uint8_t
  {
    Plane = 0,
    Sphere = 1,
    Cylinder = 2,
    Cone = 3,
    Torus = 4,
    Unknown = 5
  };

  std::vector<PrimitiveType> AllManagedPrimitives();

  std::vector<std::string> AllManagedPrimitivesStr();

  std::string ConvertPrimitiveToString(PrimitiveType type);
  
  PrimitiveType ConvertStringToPrimitive(const std::string& type);

  int ConvertPrimitiveToInt(PrimitiveType type);
};
