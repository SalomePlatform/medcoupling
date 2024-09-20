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
#include <vector>
#include <array>
#include "MEDCouplingUMesh.hxx"
#include "PrimitiveType.hxx"
#include "ShapeRecognDefines.hxx"
namespace MEDCoupling
{
    class Nodes;

    class SHAPE_RECOGNITION_EXPORT NodesBuilder
    {
    public:
        NodesBuilder(const MEDCouplingUMesh *);

        Nodes *build();

    private:
        void computeNormals();
        void computeCurvatures(double tol = 0.000001);
        void computeCurvatures(mcIdType nodeId, double tol);
        PrimitiveType findPrimitiveType(double k1, double k2, double kdiff0, double kis0) const;
        PrimitiveType findPrimitiveType2(double k1, double k2, double kdiff0, double kis0) const;
        std::vector<double> computeNormalCurvatureCoefficients(
            const std::vector<double> &discreteCurvatures,
            const std::vector<double> &tangents,
            const std::array<double, 3> &normal,
            const std::array<double, 3> &e1) const;
        void computeCellNormal(const std::vector<mcIdType> &nodeIds, std::array<double, 3> &cellNormal) const;
        double computeAverageDistance(mcIdType nodeId, const std::vector<mcIdType> &neighborIds) const;
        std::vector<double> computeDiscreteCurvatures(mcIdType nodeId, const std::vector<mcIdType> &neighborIds) const;
        double computeDiscreteCurvature(mcIdType nodeId, mcIdType neighborId) const;
        std::vector<double> computeTangentDirections(mcIdType nodeId, const std::vector<mcIdType> &neighborIds) const;

        const MEDCouplingUMesh *mesh;
        Nodes *nodes;
    };
}
