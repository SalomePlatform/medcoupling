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

#include "NodesBuilder.hxx"
#include "Nodes.hxx"
#include "MathOps.hxx"
#include "ShapeRecongConstants.hxx"

using namespace MEDCoupling;

NodesBuilder::NodesBuilder(const MEDCouplingUMesh *mesh) : mesh(mesh) {}

Nodes *
NodesBuilder::build()
{
    DataArrayInt64 *neighbors;
    DataArrayInt64 *neighborsIdx;
    mesh->computeNeighborsOfNodes(neighbors, neighborsIdx);
    nodes = new Nodes(mesh, neighbors, neighborsIdx);
    computeNormals();
    computeCurvatures();
    return nodes;
}

void
NodesBuilder::computeNormals()
{
    mcIdType nbNodes = mesh->getNumberOfNodes();
    mcIdType nbCells = mesh->getNumberOfCells();
    std::array<double, 3> cellNormal;
    cellNormal.fill(0.0);
    // product of normal AreaByCell
    std::vector<double> prodNormalAreaByCell(3 * nbCells, 0.0);
    for (int cellId = 0; cellId < nbCells; cellId++)
    {
        std::vector<mcIdType> nodeIds;
        mesh->getNodeIdsOfCell(cellId, nodeIds);
        computeCellNormal(nodeIds, cellNormal);
        prodNormalAreaByCell[3 * cellId + 0] = cellNormal[0] / 2.0;
        prodNormalAreaByCell[3 * cellId + 1] = cellNormal[1] / 2.0;
        prodNormalAreaByCell[3 * cellId + 2] = cellNormal[2] / 2.0;
    }
    //
    nodes->normals.resize(3 * nbNodes, 1.0);
    MCAuto<DataArrayInt64> revNodal = DataArrayInt64::New();
    MCAuto<DataArrayInt64> revNodalIdx = DataArrayInt64::New();
    mesh->getReverseNodalConnectivity(revNodal, revNodalIdx);
    for (size_t nodeId = 0; nodeId < (size_t)nbNodes; nodeId++)
    {
        mcIdType nbCells = revNodalIdx->getIJ(nodeId + 1, 0) - revNodalIdx->getIJ(nodeId, 0);
        std::vector<mcIdType> cellIds(nbCells, 0);
        mcIdType start = revNodalIdx->getIJ(nodeId, 0);
        for (size_t i = 0; i < cellIds.size(); ++i) cellIds[i] = revNodal->getIJ(start + i, 0);
        double normal = 0.0;
        for (size_t i = 0; i < 3; i++)
        {
            nodes->normals[3 * nodeId + i] = 0.0;
            for (mcIdType j = 0; j < nbCells; j++)
            {
                nodes->normals[3 * nodeId + i] += prodNormalAreaByCell[3 * cellIds[j] + i];
            }
            nodes->normals[3 * nodeId + i] /= (double)nbCells;
            normal += pow(nodes->normals[3 * nodeId + i], 2);
        }
        for (size_t i = 0; i < 3; i++) nodes->normals[3 * nodeId + i] /= sqrt(normal);
    }
}

void
NodesBuilder::computeCurvatures(double tol)
{
    mcIdType nbNodes = mesh->getNumberOfNodes();
    nodes->k1.resize(nbNodes);
    nodes->k2.resize(nbNodes);
    nodes->adimK1.resize(nbNodes);
    nodes->adimK2.resize(nbNodes);
    nodes->primitives.resize(nbNodes);
    nodes->weakDirections.resize(3 * nbNodes);
    nodes->mainDirections.resize(3 * nbNodes);
    for (mcIdType nodeId = 0; nodeId < nbNodes; nodeId++) computeCurvatures(nodeId, tol);
}

void
NodesBuilder::computeCurvatures(mcIdType nodeId, double tol)
{
    std::array<double, 3> normal = nodes->getNormal(nodeId);
    std::vector<mcIdType> neighborIds = nodes->getNeighbors(nodeId);
    double theta0 = 0.0;
    double k1 = 0.0;
    double k2 = 0.0;
    double adimK1 = 0.0;
    double adimK2 = 0.0;
    PrimitiveType primitive = PrimitiveType::Unknown;
    std::array<double, 3> mainDirection{0.0, 0.0, 0.0};
    std::array<double, 3> weakDirection{0.0, 0.0, 0.0};
    if (neighborIds.size() > 0)
    {
        std::vector<double> discreteCurvatures = computeDiscreteCurvatures(nodeId, neighborIds);
        std::vector<double> tangents = computeTangentDirections(nodeId, neighborIds);
        mcIdType maxCurvatureId = std::distance(
            discreteCurvatures.begin(), std::max_element(discreteCurvatures.begin(), discreteCurvatures.end())
        );
        std::array<double, 3> e1{
            tangents[3 * maxCurvatureId], tangents[3 * maxCurvatureId + 1], tangents[3 * maxCurvatureId + 2]
        };
        std::array<double, 3> e2 = MathOps::normalize(MathOps::cross(e1, normal));
        std::vector<double> coeffs = computeNormalCurvatureCoefficients(discreteCurvatures, tangents, normal, e1);
        double a = coeffs[0], b = coeffs[1], c = coeffs[2];
        double h = (a + c) / 2.0;
        double kg = a * c - pow(b, 2) / 4.0;
        k1 = h + sqrt(fabs(pow(h, 2) - kg));
        k2 = h - sqrt(fabs(pow(h, 2) - kg));
        if (fabs(k1 - k2) >= tol)
            theta0 = 0.5 * asin(b / (k1 - k2));
        std::array<double, 3> direction1{0.0, 0.0, 0.0};
        std::array<double, 3> direction2{0.0, 0.0, 0.0};
        for (size_t i = 0; i < 3; ++i)
        {
            direction1[i] = cos(theta0) * e1[i] + sin(theta0) * e2[i];
            direction2[i] = -sin(theta0) * e1[i] + cos(theta0) * e2[i];
        }
        double averageDistance = computeAverageDistance(nodeId, neighborIds);
        double adimK1 = k1 * averageDistance;
        double adimK2 = k2 * averageDistance;
        double adimKdiff0, adimKis0;
        if (fabs(k1) > fabs(k2))
        {
            adimKdiff0 = adimK1;
            adimKis0 = adimK2;
            mainDirection = direction1;
            weakDirection = direction2;
        }
        else
        {
            adimKdiff0 = adimK2;
            adimKis0 = adimK1;
            mainDirection = direction2;
            weakDirection = direction1;
        }
        primitive = findPrimitiveType(adimK1, adimK2, adimKdiff0, adimKis0);
    }
    // Populate nodes
    nodes->k1[nodeId] = k1;
    nodes->k2[nodeId] = k2;
    nodes->adimK1[nodeId] = adimK1;
    nodes->adimK2[nodeId] = adimK2;
    for (size_t i = 0; i < 3; ++i)
    {
        nodes->weakDirections[3 * nodeId + i] = weakDirection[i];
        nodes->mainDirections[3 * nodeId + i] = mainDirection[i];
    }
    nodes->primitives[nodeId] = primitive;
}

PrimitiveType
NodesBuilder::findPrimitiveType(double k1, double k2, double kdiff0, double kis0) const
{
    if ((fabs(k1 - k2) < EPSILON_PRIMITIVE) && (fabs((k1 + k2) / 2) < EPSILON_PRIMITIVE))
        return PrimitiveType::Plane;
    else if ((fabs(k1 - k2) < EPSILON_PRIMITIVE) && (fabs((k1 + k2) / 2) > EPSILON_PRIMITIVE))
        return PrimitiveType::Sphere;
    else if ((fabs(kdiff0) > EPSILON_PRIMITIVE) && (fabs(kis0) < EPSILON_PRIMITIVE))
        return PrimitiveType::Cylinder;
    else if ((fabs(k1) > EPSILON_PRIMITIVE) && (fabs(k2) > EPSILON_PRIMITIVE))
        return PrimitiveType::Torus;
    else
        return PrimitiveType::Unknown;
}

PrimitiveType
NodesBuilder::findPrimitiveType2(double k1, double k2, double kdiff0, double kis0) const
{
    double epsilon2 = pow(EPSILON_PRIMITIVE, 2);
    double diffCurvature = fabs(k1 - k2);
    double gaussianCurvature = k1 * k2;
    double meanCurvature = (k1 + k2) / 2.0;
    if (fabs(k1) < EPSILON_PRIMITIVE && fabs(k2) < EPSILON_PRIMITIVE && gaussianCurvature < epsilon2 &&
        meanCurvature < EPSILON_PRIMITIVE)
        return PrimitiveType::Plane;
    else if (diffCurvature < EPSILON_PRIMITIVE && k1 > EPSILON_PRIMITIVE && k2 > EPSILON_PRIMITIVE)
        return PrimitiveType::Sphere;
    else if ((fabs(k1) > EPSILON_PRIMITIVE && fabs(k2) < EPSILON_PRIMITIVE) ||
             (fabs(k1) < EPSILON_PRIMITIVE && fabs(k2) > EPSILON_PRIMITIVE))
        return PrimitiveType::Cylinder;
    else if (std::signbit(k1) != std::signbit(k2) || (fabs(k1) < EPSILON_PRIMITIVE && fabs(k2) < EPSILON_PRIMITIVE))
        return PrimitiveType::Torus;
    else
        return PrimitiveType::Unknown;
}

std::vector<double>
NodesBuilder::computeNormalCurvatureCoefficients(
    const std::vector<double> &discreteCurvatures,
    const std::vector<double> &tangents,
    const std::array<double, 3> &normal,
    const std::array<double, 3> &e1
) const
{
    size_t nbNeighbors = discreteCurvatures.size();
    std::vector<double> a(3 * nbNeighbors, 0.0);
    for (size_t i = 0; i < nbNeighbors; ++i)
    {
        std::array<double, 3> tangent{tangents[3 * i], tangents[3 * i + 1], tangents[3 * i + 2]};
        double theta = MathOps::computeOrientedAngle(normal, tangent, e1);
        double cos_theta = cos(theta);
        double sin_theta = sin(theta);
        a[i] = pow(cos_theta, 2);
        a[nbNeighbors + i] = cos_theta * sin_theta;
        a[2 * nbNeighbors + i] = pow(sin_theta, 2);
    }
    return MathOps::lstsq(a, discreteCurvatures);
}

void
NodesBuilder::computeCellNormal(const std::vector<mcIdType> &nodeIds, std::array<double, 3> &cellNormal) const
{
    std::vector<double> point1;
    std::vector<double> point2;
    std::vector<double> point3;
    mesh->getCoordinatesOfNode(nodeIds[0], point1);
    mesh->getCoordinatesOfNode(nodeIds[1], point2);
    mesh->getCoordinatesOfNode(nodeIds[2], point3);
    std::array<double, 3> a;
    a.fill(3);
    std::array<double, 3> b;
    b.fill(3);
    for (int i = 0; i < 3; i++)
    {
        a[i] = point2[i] - point1[i];
        b[i] = point3[i] - point1[i];
    }
    cellNormal[0] = a[1] * b[2] - a[2] * b[1];
    cellNormal[1] = a[2] * b[0] - a[0] * b[2];
    cellNormal[2] = a[0] * b[1] - a[1] * b[0];
}

double
NodesBuilder::computeAverageDistance(mcIdType nodeId, const std::vector<mcIdType> &neighborIds) const
{
    double distance = 0.0;
    std::array<double, 3> nodeCoords = nodes->getCoordinates(nodeId);
    for (size_t i = 0; i < neighborIds.size(); ++i)
    {
        std::array<double, 3> neighborCoords = nodes->getCoordinates(neighborIds[i]);
        double distanceToNeighbor = 0.0;
        for (size_t j = 0; j < 3; ++j) distanceToNeighbor += pow(neighborCoords[j] - nodeCoords[j], 2);
        distance += sqrt(distanceToNeighbor);
    }
    return distance / (double)neighborIds.size();
}

std::vector<double>
NodesBuilder::computeDiscreteCurvatures(mcIdType nodeId, const std::vector<mcIdType> &neighborIds) const
{
    std::vector<double> discreteCurvatures(neighborIds.size(), 0.0);
    for (size_t i = 0; i < neighborIds.size(); ++i)
        discreteCurvatures[i] = computeDiscreteCurvature(nodeId, neighborIds[i]);
    return discreteCurvatures;
}

double
NodesBuilder::computeDiscreteCurvature(mcIdType nodeId, mcIdType neighborId) const
{
    double curvature = 0.0;
    double n = 0.0;
    for (size_t i = 0; i < 3; i++)
    {
        double ni = nodes->coords->getIJ(neighborId, i) - nodes->coords->getIJ(nodeId, i);
        curvature += ni * (nodes->normals[3 * neighborId + i] - nodes->normals[3 * nodeId + i]);
        n += ni * ni;
    }
    return curvature / n;
}

std::vector<double>
NodesBuilder::computeTangentDirections(mcIdType nodeId, const std::vector<mcIdType> &neighborIds) const
{
    size_t nbNeighbors = neighborIds.size();
    std::vector<double> tangent(3 * nbNeighbors, 0.0);
    for (size_t i = 0; i < nbNeighbors; ++i)
    {
        mcIdType neighborId = neighborIds[i];
        double s = 0.0;
        for (size_t j = 0; j < 3; j++)
        {
            tangent[3 * i + j] = nodes->coords->getIJ(neighborId, j) - nodes->coords->getIJ(nodeId, j);
            s += tangent[3 * i + j] * nodes->normals[3 * nodeId + j];
        }
        double n = 0.0;
        for (size_t j = 0; j < 3; j++)
        {
            tangent[3 * i + j] -= s * nodes->normals[3 * nodeId + j];
            n += tangent[3 * i + j] * tangent[3 * i + j];
        }
        for (size_t j = 0; j < 3; j++) tangent[3 * i + j] /= sqrt(n);
    }
    return tangent;
}
