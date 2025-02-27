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

#include <string>

#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldInt32.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "ShapeRecognDefines.hxx"
namespace MEDCoupling
{
    class SHAPE_RECOGNITION_EXPORT ShapeRecognMesh : public RefCountObject
    {
        friend class ShapeRecognMeshBuilder;

    public:
        static ShapeRecognMesh *New();
        std::size_t getHeapMemorySizeWithoutChildren() const;
        std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;

        // Node properties
        const MEDCouplingFieldDouble *getNodeK1() const;
        const MEDCouplingFieldDouble *getNodeK2() const;
        const MEDCouplingFieldInt32  *getNodePrimitiveType() const;
        const MEDCouplingFieldDouble *getNodeNormal() const;
        // see ShapeRecognMeshBuilder::buildNodeWeakDirections
        // see ShapeRecognMeshBuilder::buildNodeMainDirections

        // Area properties
        const MEDCouplingFieldInt32  *getAreaId() const;
        const MEDCouplingFieldInt32  *getAreaPrimitiveType() const;
        const MEDCouplingFieldDouble *getAreaNormal() const;
        const MEDCouplingFieldDouble *getMinorRadius() const;
        const MEDCouplingFieldDouble *getRadius() const;
        const MEDCouplingFieldDouble *getAngle() const;
        const MEDCouplingFieldDouble *getCenter() const;
        const MEDCouplingFieldDouble *getAxis() const;
        const MEDCouplingFieldDouble *getApex() const;

        // see ShapeRecognMeshBuilder::buildAreaAxisPoint
        // see ShapeRecognMeshBuilder::buildAreaAffinePoint

    protected:
        ShapeRecognMesh();

    private:
        MCAuto<MEDCouplingFieldDouble> nodeK1;
        MCAuto<MEDCouplingFieldDouble> nodeK2;
        MCAuto<MEDCouplingFieldInt32>  nodePrimitiveType;
        MCAuto<MEDCouplingFieldDouble> nodeNormal;
        MCAuto<MEDCouplingFieldInt32>  areaId;
        MCAuto<MEDCouplingFieldInt32>  areaPrimitiveType;
        MCAuto<MEDCouplingFieldDouble> areaNormal;
        MCAuto<MEDCouplingFieldDouble> minorRadius;
        MCAuto<MEDCouplingFieldDouble> radius;
        MCAuto<MEDCouplingFieldDouble> angle;
        MCAuto<MEDCouplingFieldDouble> center;
        MCAuto<MEDCouplingFieldDouble> axis;
        MCAuto<MEDCouplingFieldDouble> apex;
    };
}
