// Copyright (C) 2007-2024  CEA, EDF
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

#ifndef __PRIMITIVETYPE_HXX__
#define __PRIMITIVETYPE_HXX__

#include <string>
namespace MEDCoupling
{
    enum PrimitiveType
    {
        Plane = 0,
        Sphere = 1,
        Cylinder = 2,
        Cone = 3,
        Torus = 4,
        Unknown = 5
    };

    inline std::string convertPrimitiveToString(PrimitiveType type)
    {
        std::string typeName = "";
        switch (type)
        {
        case PrimitiveType::Plane:
            typeName = "Plane";
            break;
        case PrimitiveType::Sphere:
            typeName = "Sphere";
            break;
        case PrimitiveType::Cylinder:
            typeName = "Cylinder";
            break;
        case PrimitiveType::Cone:
            typeName = "Cone";
            break;
        case PrimitiveType::Torus:
            typeName = "Torus";
            break;
        case PrimitiveType::Unknown:
            typeName = "Unknown";
            break;
        default:
            break;
        }
        return typeName;
    };

    inline int convertPrimitiveToInt(PrimitiveType type)
    {
        int typeInt = 5;
        switch (type)
        {
        case PrimitiveType::Plane:
            typeInt = 0;
            break;
        case PrimitiveType::Sphere:
            typeInt = 1;
            break;
        case PrimitiveType::Cylinder:
            typeInt = 2;
            break;
        case PrimitiveType::Cone:
            typeInt = 3;
            break;
        case PrimitiveType::Torus:
            typeInt = 4;
            break;
        case PrimitiveType::Unknown:
        default:
            break;
        }
        return typeInt;
    };
};

#endif // __PRIMITIVETYPE_HXX__