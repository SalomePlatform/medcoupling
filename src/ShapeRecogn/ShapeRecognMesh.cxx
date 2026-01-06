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

#include "ShapeRecognMesh.hxx"

using namespace MEDCoupling;

ShapeRecognMesh::ShapeRecognMesh()
    : nodeK1(nullptr),
      nodeK2(nullptr),
      nodePrimitiveType(nullptr),
      nodeNormal(nullptr),
      areaId(nullptr),
      areaPrimitiveType(nullptr),
      areaNormal(nullptr),
      minorRadius(nullptr),
      radius(nullptr),
      angle(nullptr),
      center(nullptr),
      axis(nullptr),
      apex(nullptr)
{
}

std::size_t
ShapeRecognMesh::getHeapMemorySizeWithoutChildren() const
{
    return 0;
}

std::vector<const BigMemoryObject *>
ShapeRecognMesh::getDirectChildrenWithNull() const
{
    std::vector<const BigMemoryObject *> ret;
    ret.push_back(nodeK1);
    ret.push_back(nodeK2);
    ret.push_back(nodePrimitiveType);
    ret.push_back(nodeNormal);
    ret.push_back(areaId);
    ret.push_back(areaPrimitiveType);
    ret.push_back(areaNormal);
    ret.push_back(minorRadius);
    ret.push_back(radius);
    ret.push_back(angle);
    ret.push_back(center);
    ret.push_back(axis);
    ret.push_back(apex);
    return ret;
}

ShapeRecognMesh *
ShapeRecognMesh::New()
{
    return new ShapeRecognMesh;
}

const MEDCouplingFieldDouble *
ShapeRecognMesh::getNodeK1() const
{
    return nodeK1;
}

const MEDCouplingFieldDouble *
ShapeRecognMesh::getNodeK2() const
{
    return nodeK2;
}

const MEDCouplingFieldInt32 *
ShapeRecognMesh::getNodePrimitiveType() const
{
    return nodePrimitiveType;
}

const MEDCouplingFieldDouble *
ShapeRecognMesh::getNodeNormal() const
{
    return nodeNormal;
}

const MEDCouplingFieldInt32 *
ShapeRecognMesh::getAreaId() const
{
    return areaId;
}

const MEDCouplingFieldInt32 *
ShapeRecognMesh::getAreaPrimitiveType() const
{
    return areaPrimitiveType;
}

const MEDCouplingFieldDouble *
ShapeRecognMesh::getAreaNormal() const
{
    return areaNormal;
}

const MEDCouplingFieldDouble *
ShapeRecognMesh::getMinorRadius() const
{
    return minorRadius;
}

const MEDCouplingFieldDouble *
ShapeRecognMesh::getRadius() const
{
    return radius;
}

const MEDCouplingFieldDouble *
ShapeRecognMesh::getAngle() const
{
    return angle;
}

const MEDCouplingFieldDouble *
ShapeRecognMesh::getCenter() const
{
    return center;
}

const MEDCouplingFieldDouble *
ShapeRecognMesh::getAxis() const
{
    return axis;
}

const MEDCouplingFieldDouble *
ShapeRecognMesh::getApex() const
{
    return apex;
}
