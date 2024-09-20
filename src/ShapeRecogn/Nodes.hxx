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
#include "MEDCouplingUMesh.hxx"
#include "PrimitiveType.hxx"
#include <array>
#include "ShapeRecognDefines.hxx"

namespace MEDCoupling
{
    class SHAPE_RECOGNITION_EXPORT Nodes
    {
    public:
        friend class NodesBuilder;
        Nodes(const MEDCouplingUMesh *mesh,
              const DataArrayInt64 *neighbors,
              const DataArrayInt64 *neighborsIdx);

        mcIdType getNbNodes() const;
        const std::vector<double> &getK1() const;
        const std::vector<double> &getK2() const;
        const std::vector<PrimitiveType> &getPrimitiveType() const;
        const std::vector<double> &getNormals() const;
        const std::vector<double> &getWeakDirections() const;
        const std::vector<double> &getMainDirections() const;

        std::array<double, 3> getNormal(mcIdType nodeId) const;
        double getK1(mcIdType nodeId) const;
        double getK2(mcIdType nodeId) const;
        double getKdiff0(mcIdType nodeId) const;
        double getAdimK1(mcIdType nodeId) const;
        double getAdimK2(mcIdType nodeId) const;
        double getAdimKdiff0(mcIdType nodeId) const;
        const std::vector<mcIdType> getNeighbors(mcIdType nodeId) const;
        PrimitiveType getPrimitiveType(mcIdType nodeId) const;
        std::array<double, 3> getWeakDirection(mcIdType nodeId) const;
        std::array<double, 3> getMainDirection(mcIdType nodeId) const;
        std::array<double, 3> getCoordinates(mcIdType nodeId) const;

    private:
        MCConstAuto<MEDCouplingUMesh> mesh;
        const DataArrayDouble *coords;
        // normals 3 * nbNodes
        std::vector<double> normals;
        // neighbors
        MCConstAuto<DataArrayInt64> neighbors;
        MCConstAuto<DataArrayInt64> neighborsIdx;
        // curvature
        std::vector<double> k1;
        std::vector<double> k2;
        std::vector<double> adimK1;
        std::vector<double> adimK2;
        std::vector<double> weakDirections;
        std::vector<double> mainDirections;
        std::vector<PrimitiveType> primitives;
    };
}
