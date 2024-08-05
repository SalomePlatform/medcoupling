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

#ifndef __SHAPERECOGNMESH_HXX__
#define __SHAPERECOGNMESH_HXX__

#include <string>

#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingRefCountObject.hxx"

namespace MEDCoupling
{
    class ShapeRecognMesh : public RefCountObject
    {
        friend class ShapeRecognMeshBuilder;

    public:
        static ShapeRecognMesh *New();
        std::size_t getHeapMemorySizeWithoutChildren() const;
        std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;

        void save(const std::string &outputFile, bool writeFromScratch = true) const;

        // Node properties
        MEDCoupling::MEDCouplingFieldDouble *getNodeK1() const;
        MEDCoupling::MEDCouplingFieldDouble *getNodeK2() const;
        MEDCoupling::MEDCouplingFieldDouble *getNodePrimitiveType() const;
        MEDCoupling::MEDCouplingFieldDouble *getNodeNormal() const;

        // Area properties
        MEDCoupling::MEDCouplingFieldDouble *getAreaId() const;
        MEDCoupling::MEDCouplingFieldDouble *getAreaPrimitiveType() const;
        MEDCoupling::MEDCouplingFieldDouble *getAreaNormal() const;
        MEDCoupling::MEDCouplingFieldDouble *getMinorRadius() const;
        MEDCoupling::MEDCouplingFieldDouble *getRadius() const;
        MEDCoupling::MEDCouplingFieldDouble *getAngle() const;
        MEDCoupling::MEDCouplingFieldDouble *getCenter() const;
        MEDCoupling::MEDCouplingFieldDouble *getAxis() const;
        MEDCoupling::MEDCouplingFieldDouble *getApex() const;

    protected:
        ShapeRecognMesh();
        ~ShapeRecognMesh();

    private:
        MEDCoupling::MEDCouplingFieldDouble *nodeK1;
        MEDCoupling::MEDCouplingFieldDouble *nodeK2;
        MEDCoupling::MEDCouplingFieldDouble *nodePrimitiveType;
        MEDCoupling::MEDCouplingFieldDouble *nodeNormal;
        MEDCoupling::MEDCouplingFieldDouble *areaId;
        MEDCoupling::MEDCouplingFieldDouble *areaPrimitiveType;
        MEDCoupling::MEDCouplingFieldDouble *areaNormal;
        MEDCoupling::MEDCouplingFieldDouble *minorRadius;
        MEDCoupling::MEDCouplingFieldDouble *radius;
        MEDCoupling::MEDCouplingFieldDouble *angle;
        MEDCoupling::MEDCouplingFieldDouble *center;
        MEDCoupling::MEDCouplingFieldDouble *axis;
        MEDCoupling::MEDCouplingFieldDouble *apex;
    };
};

#endif // __SHAPERECOGNMESH_HXX__
