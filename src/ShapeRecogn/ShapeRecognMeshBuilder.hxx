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

#ifndef __SHAPERECOGNMESHBUILDER_HXX__
#define __SHAPERECOGNMESHBUILDER_HXX__

#include <string>

#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"

namespace MEDCoupling
{
    class Nodes;
    class Areas;
    class ShapeRecognMesh;

    class ShapeRecognMeshBuilder
    {
    public:
        ShapeRecognMeshBuilder(const std::string &fileName, int meshDimRelToMax = 0);
        ~ShapeRecognMeshBuilder();

        const Nodes *getNodes() const;
        const Areas *getAreas() const;

        ShapeRecognMesh *recognize();

    private:
        // Node properties
        MEDCoupling::MEDCouplingFieldDouble *buildNodeK1() const;
        MEDCoupling::MEDCouplingFieldDouble *buildNodeK2() const;
        MEDCoupling::MEDCouplingFieldDouble *buildNodePrimitiveType() const;
        MEDCoupling::MEDCouplingFieldDouble *buildNodeNormal() const;

        // Area properties
        MEDCoupling::MEDCouplingFieldDouble *buildAreaId() const;
        MEDCoupling::MEDCouplingFieldDouble *buildAreaPrimitiveType() const;
        MEDCoupling::MEDCouplingFieldDouble *buildAreaNormal() const;
        MEDCoupling::MEDCouplingFieldDouble *buildMinorRadius() const;
        MEDCoupling::MEDCouplingFieldDouble *buildRadius() const;
        MEDCoupling::MEDCouplingFieldDouble *buildAngle() const;
        MEDCoupling::MEDCouplingFieldDouble *buildCenter() const;
        MEDCoupling::MEDCouplingFieldDouble *buildAxis() const;
        MEDCoupling::MEDCouplingFieldDouble *buildApex() const;

        template <typename T>
        MEDCouplingFieldDouble *buildField(
            const std::string &name,
            size_t nbOfCompo,
            const std::vector<T> &values) const;
        MEDCouplingFieldDouble *buildField(
            const std::string &name,
            size_t nbOfCompo,
            double *values) const;
        MEDCouplingFieldDouble *buildField(
            const std::string &name,
            size_t nbOfCompo,
            DataArrayDouble *values) const;
        double *buildArea3DArray(const std::array<double, 3> &(*areaFunc)(Areas *, mcIdType)) const;
        double *buildAreaArray(double (*areaFunc)(Areas *, mcIdType)) const;

        const MEDCouplingUMesh *mesh;
        Nodes *nodes = nullptr;
        Areas *areas = nullptr;
    };
};

#endif // __SHAPERECOGNMESHBUILDER_HXX__
