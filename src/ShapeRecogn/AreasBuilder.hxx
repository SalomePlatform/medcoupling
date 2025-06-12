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

#include "Nodes.hxx"
#include "Areas.hxx"
#include "ShapeRecognDefines.hxx"
namespace MEDCoupling
{
class SHAPE_RECOGNITION_EXPORT AreasBuilder
{
   public:
    AreasBuilder(const Nodes *nodes);

    void build();

    Areas *getAreas() const;

   private:
    void explore();
    void expand();
    void rebuild();
    void exploreAreas();
    void expandAreas();
    void expandAreasByType(PrimitiveType primitive);
    void rebuildInvalidAreas();
    void filterHighPass();
    bool doesItMatch(mcIdType areaId, mcIdType nodeId) const;
    bool doesItBelong(mcIdType areaId, mcIdType nodeId) const;
    bool isInvalidCylinderNode(mcIdType nodeId) const;

    static double distanceToPlane(
        const std::array<double, 3> &nodeCoords, const std::array<double, 3> &point, const std::array<double, 3> &normal
    );
    static double distanceToSphere(const std::array<double, 3> &nodeCoords, const std::array<double, 3> &center);
    static double distanceToCylinder(
        const std::array<double, 3> &nodeCoords,
        const std::array<double, 3> &axis,
        const std::array<double, 3> &axisPoint
    );
    static double distanceToCone(
        const std::array<double, 3> &nodeCoords,
        const std::array<double, 3> &axis,
        const std::array<double, 3> &apex,
        double angle
    );

    const Nodes *nodes;
    Areas *areas;

    size_t threshold = 5;
};
}  // namespace MEDCoupling
