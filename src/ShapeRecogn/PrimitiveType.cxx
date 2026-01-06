// Copyright (C) 2024-2026  CEA, EDF
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

#include "PrimitiveType.hxx"

#include <iterator>
#include <algorithm>

constexpr char PLANE_STR[] = "Plane";
constexpr char SPHERE_STR[] = "Sphere";
constexpr char CYLINDER_STR[] = "Cylinder";
constexpr char CONE_STR[] = "Cone";
constexpr char TORUS_STR[] = "Torus";
constexpr char UNKNOWN_STR[] = "Unknown";

std::vector<MEDCoupling::PrimitiveType>
MEDCoupling::AllManagedPrimitives()
{
    return {
        PrimitiveType::Plane,
        PrimitiveType::Sphere,
        PrimitiveType::Cylinder,
        PrimitiveType::Cone,
        PrimitiveType::Torus,
        PrimitiveType::Unknown
    };
}

std::vector<std::string>
MEDCoupling::AllManagedPrimitivesStr()
{
    std::vector<PrimitiveType> apt(AllManagedPrimitives());
    std::vector<std::string> ret;
    std::transform(
        apt.begin(),
        apt.end(),
        std::back_inserter(ret),
        [](PrimitiveType type) { return MEDCoupling::ConvertPrimitiveToString(type); }
    );
    return ret;
}

std::string
MEDCoupling::ConvertPrimitiveToString(MEDCoupling::PrimitiveType type)
{
    std::string typeName;
    switch (type)
    {
        case PrimitiveType::Plane:
            typeName = PLANE_STR;
            break;
        case PrimitiveType::Sphere:
            typeName = SPHERE_STR;
            break;
        case PrimitiveType::Cylinder:
            typeName = CYLINDER_STR;
            break;
        case PrimitiveType::Cone:
            typeName = CONE_STR;
            break;
        case PrimitiveType::Torus:
            typeName = TORUS_STR;
            break;
        case PrimitiveType::Unknown:
            typeName = UNKNOWN_STR;
            break;
        default:
            break;
    }
    return typeName;
};

MEDCoupling::PrimitiveType
MEDCoupling::ConvertStringToPrimitive(const std::string &type)
{
    if (type == PLANE_STR)
        return PrimitiveType::Plane;
    if (type == SPHERE_STR)
        return PrimitiveType::Sphere;
    if (type == CYLINDER_STR)
        return PrimitiveType::Cylinder;
    if (type == CONE_STR)
        return PrimitiveType::Cone;
    if (type == TORUS_STR)
        return PrimitiveType::Torus;
    if (type == UNKNOWN_STR)
        return PrimitiveType::Unknown;
    return PrimitiveType::Unknown;
}

int
MEDCoupling::ConvertPrimitiveToInt(MEDCoupling::PrimitiveType type)
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
}
