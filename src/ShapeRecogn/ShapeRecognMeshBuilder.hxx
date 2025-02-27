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

#include "Nodes.hxx"
#include "Areas.hxx"

#include <memory>
#include <functional>
#include "ShapeRecognDefines.hxx"
namespace MEDCoupling
{
    class MEDCouplingFieldInt32;
    class MEDCouplingFieldInt64;
    class MEDCouplingFieldDouble;

    class ShapeRecognMesh;

    class SHAPE_RECOGNITION_EXPORT ShapeRecognMeshBuilder
    {
    public:
        ShapeRecognMeshBuilder(MCAuto< MEDCouplingUMesh > mesh);
        ShapeRecognMeshBuilder(MEDCouplingUMesh *mesh);

        const Nodes *getNodes() const;
        const Areas *getAreas() const;

        MCAuto<ShapeRecognMesh> recognize();

        // Node properties
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> buildNodeWeakDirections() const;
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> buildNodeMainDirections() const;

        //Area properties
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> buildAreaAxisPoint() const;
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> buildAreaAffinePoint() const;
    private:
        // Node properties
        MEDCoupling::MEDCouplingFieldDouble *buildNodeK1() const;
        MEDCoupling::MEDCouplingFieldDouble *buildNodeK2() const;
        MEDCoupling::MEDCouplingFieldInt32  *buildNodePrimitiveType() const;
        MEDCoupling::MEDCouplingFieldDouble *buildNodeNormal() const;

        // Area properties
        MEDCoupling::MEDCouplingFieldInt32  *buildAreaId() const;
        MEDCoupling::MEDCouplingFieldInt32  *buildAreaPrimitiveType() const;
        MEDCoupling::MEDCouplingFieldDouble *buildAreaNormal() const;
        MEDCoupling::MEDCouplingFieldDouble *buildMinorRadius() const;
        MEDCoupling::MEDCouplingFieldDouble *buildRadius() const;
        MEDCoupling::MEDCouplingFieldDouble *buildAngle() const;
        MEDCoupling::MEDCouplingFieldDouble *buildCenter() const;
        MEDCoupling::MEDCouplingFieldDouble *buildAxis() const;
        MEDCoupling::MEDCouplingFieldDouble *buildApex() const;

        double *buildArea3DArray(std::function<std::array<double, 3>(Areas *, Int32)> areaFunc) const;
        double *buildAreaArray(std::function<double(Areas *, Int32)> areaFunc) const;
        void assign(MCAuto< MEDCouplingUMesh > mesh);
        void checkNodesBeforeBuildingField() const;
        void checkAreasBeforeBuildingField() const;
    private:
        MCConstAuto< MEDCouplingUMesh > mesh;
        std::unique_ptr<Nodes> nodes;
        std::unique_ptr<Areas> areas;
    };
}
