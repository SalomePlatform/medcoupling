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
#include "ShapeRecognDefines.hxx"
namespace MEDCoupling
{
enum class SHAPE_RECOGNITION_EXPORT PrimitiveType : unsigned char
{
    Plane = 0,
    Sphere = 1,
    Cylinder = 2,
    Cone = 3,
    Torus = 4,
    Unknown = 5
};

SHAPE_RECOGNITION_EXPORT std::vector<PrimitiveType>
AllManagedPrimitives();

SHAPE_RECOGNITION_EXPORT std::vector<std::string>
AllManagedPrimitivesStr();

SHAPE_RECOGNITION_EXPORT std::string
ConvertPrimitiveToString(PrimitiveType type);

SHAPE_RECOGNITION_EXPORT PrimitiveType
ConvertStringToPrimitive(const std::string &type);

SHAPE_RECOGNITION_EXPORT int
ConvertPrimitiveToInt(PrimitiveType type);
};  // namespace MEDCoupling
