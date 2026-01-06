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

#pragma once
#include "ShapeRecognDefines.hxx"
#include "PrimitiveType.hxx"
#include "Nodes.hxx"

#include <array>

namespace MEDCoupling
{

struct Area
{
    PrimitiveType primitive = PrimitiveType::Unknown;
    double k1 = 0.0;
    double k2 = 0.0;
    double adimK1 = 0.0;
    double adimK2 = 0.0;
    double adimKdiff0 = 0.0;
    double minorRadius = 0.0;
    double radius = 0.0;
    double angle = 0.0;
    std::array<double, 3> normal{0.0, 0.0, 0.0};
    std::array<double, 3> center{0.0, 0.0, 0.0};
    std::array<double, 3> axis{0.0, 0.0, 0.0};
    std::array<double, 3> axisPoint{0.0, 0.0, 0.0};
    std::array<double, 3> apex{0.0, 0.0, 0.0};
    std::vector<mcIdType> nodeIds;
};

class SHAPE_RECOGNITION_EXPORT Areas
{
   public:
    Areas(const Nodes *nodes);

    mcIdType addArea(PrimitiveType primitive = PrimitiveType::Unknown);
    void cleanArea(mcIdType areaId);
    void cancelArea(mcIdType areaId, PrimitiveType primitive = PrimitiveType::Unknown);
    void removeArea(mcIdType areaId);
    void addNode(mcIdType areaId, mcIdType nodeId);
    void cleanInvalidNodeAreas();

    mcIdType getAreaId(mcIdType nodeId) const;
    const std::vector<Int32> &getAreaIdByNodes() const;

    bool isEmpty(mcIdType areaId) const;
    size_t getNumberOfAreas() const;
    size_t getNumberOfNodes(mcIdType areaId) const;

    bool isNodeCompatible(mcIdType areaId, mcIdType nodeId) const;

    PrimitiveType getPrimitiveType(mcIdType areaId) const;
    std::string getPrimitiveTypeName(mcIdType areaId) const;
    int getPrimitiveTypeInt(mcIdType areaId) const;

    const std::vector<mcIdType> &getNodeIds(mcIdType areaId) const;

    double getAdimK1(mcIdType areaId) const;
    double getAdimK2(mcIdType areaId) const;
    double getAdimKdiff0(mcIdType areaId) const;

    void computeProperties(mcIdType areaId);
    double getMinorRadius(mcIdType areaId) const;
    double getRadius(mcIdType areaId) const;
    double getAngle(mcIdType areaId) const;
    const std::array<double, 3> &getNormal(mcIdType areaId) const;
    const std::array<double, 3> &getCenter(mcIdType areaId) const;
    const std::array<double, 3> &getAxis(mcIdType areaId) const;
    const std::array<double, 3> &getAxisPoint(mcIdType areaId) const;
    const std::array<double, 3> &getApex(mcIdType areaId) const;
    std::array<double, 3> getAffinePoint(mcIdType areaId) const;

   private:
    void cleanArea(mcIdType areaId, mcIdType newAreaId);
    void removeNode(mcIdType nodeId);

    void computePlaneProperties(mcIdType areaId);
    void computeSphereProperties(mcIdType areaId);
    void computeCylinderProperties(mcIdType areaId);
    void computeConeProperties(mcIdType areaId);
    void computeTorusProperties(mcIdType areaId);

    std::vector<Area> areas;
    const Nodes *nodes;
    std::vector<Int32> areaIdByNodes;
};
}  // namespace MEDCoupling
