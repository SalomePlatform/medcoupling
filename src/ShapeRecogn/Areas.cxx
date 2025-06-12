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

#include "Areas.hxx"
#include "MathOps.hxx"

#include <random>
#include <cblas.h>

using namespace MEDCoupling;

Areas::Areas(const Nodes *nodes) : nodes(nodes), areaIdByNodes(nodes->getNbNodes(), -1) {}

mcIdType
Areas::addArea(PrimitiveType primitive)
{
    Area area;
    area.primitive = primitive;
    areas.push_back(area);
    return areas.size() - 1;
}

void
Areas::cleanArea(mcIdType areaId)
{
    cleanArea(areaId, -1);
}

void
Areas::cancelArea(mcIdType areaId, PrimitiveType primitive)
{
    cleanArea(areaId, -2);
    Area &area = areas[areaId];
    area.primitive = primitive;
}

void
Areas::removeArea(mcIdType areaId)
{
    for (size_t nodeId = 0; nodeId < areaIdByNodes.size(); ++nodeId)
    {
        if (areaIdByNodes[nodeId] == areaId)
            areaIdByNodes[nodeId] = -1;
        else if (areaIdByNodes[nodeId] > areaId)
            areaIdByNodes[nodeId] -= 1;
    }
    areas.erase(areas.begin() + areaId);
}

void
Areas::addNode(mcIdType areaId, mcIdType nodeId)
{
    removeNode(nodeId);
    areaIdByNodes[nodeId] = FromIdType<Int32>(areaId);
    Area &area = areas[areaId];
    area.nodeIds.push_back(nodeId);
    size_t nbNodes = area.nodeIds.size();
    area.k1 = ((double)(nbNodes - 1) * area.k1 + nodes->getK1(nodeId)) / (double)nbNodes;
    area.k2 = ((double)(nbNodes - 1) * area.k2 + nodes->getK2(nodeId)) / (double)nbNodes;
    area.adimK1 = ((double)(nbNodes - 1) * area.adimK1 + nodes->getAdimK1(nodeId)) / (double)nbNodes;
    area.adimK2 = ((double)(nbNodes - 1) * area.adimK2 + nodes->getAdimK2(nodeId)) / (double)nbNodes;
    area.adimKdiff0 = ((double)(nbNodes - 1) * area.adimKdiff0 + nodes->getAdimKdiff0(nodeId)) / (double)nbNodes;
}

void
Areas::cleanInvalidNodeAreas()
{
    for (mcIdType nodeId = 0; nodeId < (mcIdType)areaIdByNodes.size(); ++nodeId)
    {
        if (areaIdByNodes[nodeId] <= -2)
            areaIdByNodes[nodeId] = -1;
    }
}

mcIdType
Areas::getAreaId(mcIdType nodeId) const
{
    return areaIdByNodes[nodeId];
}

bool
Areas::isEmpty(mcIdType areaId) const
{
    return areas[areaId].nodeIds.empty();
}

size_t
Areas::getNumberOfAreas() const
{
    return areas.size();
}

size_t
Areas::getNumberOfNodes(mcIdType areaId) const
{
    return areas[areaId].nodeIds.size();
}

bool
Areas::isNodeCompatible(mcIdType areaId, mcIdType nodeId) const
{
    PrimitiveType areaType = getPrimitiveType(areaId);
    PrimitiveType nodeType = nodes->getPrimitiveType(nodeId);
    return (
        (areaType != PrimitiveType::Cylinder || nodeType != PrimitiveType::Plane) &&
        (areaType != PrimitiveType::Sphere || nodeType != PrimitiveType::Plane)
    );
}

PrimitiveType
Areas::getPrimitiveType(mcIdType areaId) const
{
    const Area &area = areas[areaId];
    if (area.primitive != PrimitiveType::Unknown)
        return area.primitive;
    else if (area.nodeIds.empty())
        return PrimitiveType::Unknown;
    else
        return nodes->getPrimitiveType(area.nodeIds[0]);
}

std::string
Areas::getPrimitiveTypeName(mcIdType areaId) const
{
    return ConvertPrimitiveToString(getPrimitiveType(areaId));
}

int
Areas::getPrimitiveTypeInt(mcIdType areaId) const
{
    return ConvertPrimitiveToInt(getPrimitiveType(areaId));
}

const std::vector<mcIdType> &
Areas::getNodeIds(mcIdType areaId) const
{
    return areas[areaId].nodeIds;
}

double
Areas::getAdimK1(mcIdType areaId) const
{
    return areas[areaId].adimK1;
}

double
Areas::getAdimK2(mcIdType areaId) const
{
    return areas[areaId].adimK2;
}

double
Areas::getAdimKdiff0(mcIdType areaId) const
{
    return areas[areaId].adimKdiff0;
}

void
Areas::computeProperties(mcIdType areaId)
{
    switch (getPrimitiveType(areaId))
    {
        case PrimitiveType::Plane:
            computePlaneProperties(areaId);
            break;
        case PrimitiveType::Sphere:
            computeSphereProperties(areaId);
            break;
        case PrimitiveType::Cylinder:
            computeCylinderProperties(areaId);
            break;
        case PrimitiveType::Cone:
            computeConeProperties(areaId);
            break;
        case PrimitiveType::Torus:
            computeTorusProperties(areaId);
            break;
        case PrimitiveType::Unknown:
        default:
            break;
    }
}

double
Areas::getMinorRadius(mcIdType areaId) const
{
    const Area &area = areas[areaId];
    return area.minorRadius;
}

double
Areas::getRadius(mcIdType areaId) const
{
    const Area &area = areas[areaId];
    return area.radius;
}

double
Areas::getAngle(mcIdType areaId) const
{
    return areas[areaId].angle;
}

const std::array<double, 3> &
Areas::getNormal(mcIdType areaId) const
{
    return areas[areaId].normal;
}

const std::array<double, 3> &
Areas::getCenter(mcIdType areaId) const
{
    const Area &area = areas[areaId];
    return area.center;
}

const std::array<double, 3> &
Areas::getAxis(mcIdType areaId) const
{
    return areas[areaId].axis;
}

const std::array<double, 3> &
Areas::getAxisPoint(mcIdType areaId) const
{
    return areas[areaId].axisPoint;
}

const std::array<double, 3> &
Areas::getApex(mcIdType areaId) const
{
    return areas[areaId].apex;
}

std::array<double, 3>
Areas::getAffinePoint(mcIdType areaId) const
{
    const Area &area = areas[areaId];
    if (area.nodeIds.empty())
        return {0.0, 0.0, 0.0};
    else
        return nodes->getCoordinates(area.nodeIds[0]);
}

void
Areas::cleanArea(mcIdType areaId, mcIdType newAreaId = -1)
{
    Area &area = areas[areaId];
    for (mcIdType nodeId : area.nodeIds) areaIdByNodes[nodeId] = FromIdType<Int32>(newAreaId);
    area.primitive = PrimitiveType::Unknown;
    area.k1 = 0.0;
    area.k2 = 0.0;
    area.adimK1 = 0.0;
    area.adimK2 = 0.0;
    area.adimKdiff0 = 0.0;
    area.minorRadius = 0.0;
    area.radius = 0.0;
    area.angle = 0.0;
    area.normal = {0.0, 0.0, 0.0};
    area.center = {0.0, 0.0, 0.0};
    area.axis = {0.0, 0.0, 0.0};
    area.axisPoint = {0.0, 0.0, 0.0};
    area.apex = {0.0, 0.0, 0.0};
    area.nodeIds.clear();
}

void
Areas::removeNode(mcIdType nodeId)
{
    mcIdType areaId = areaIdByNodes[nodeId];
    if (areaId > -1)
    {
        Area &area = areas[areaId];
        area.nodeIds.erase(std::remove(area.nodeIds.begin(), area.nodeIds.end(), nodeId), area.nodeIds.end());
        areaIdByNodes[nodeId] = -1;
        // TODO: Update the parameters of the area ?
    }
}

void
Areas::computePlaneProperties(mcIdType areaId)
{
    Area &area = areas[areaId];
    const std::vector<double> &normals = nodes->getNormals();
    mcIdType nbNodes = area.nodeIds.size();
    area.normal = {0.0, 0.0, 0.0};
    for (mcIdType nodeId : area.nodeIds)
    {
        for (size_t i = 0; i < 3; ++i) area.normal[i] += normals[3 * nodeId + i];
    }
    for (size_t i = 0; i < 3; ++i) area.normal[i] /= (double)nbNodes;
}

void
Areas::computeSphereProperties(mcIdType areaId)
{
    Area &area = areas[areaId];
    const std::vector<double> &normals = nodes->getNormals();
    area.radius = (2 / (area.k1 + area.k2));
    std::array<double, 3> center{0.0, 0.0, 0.0};
    if (!area.nodeIds.empty())
    {
        size_t nbNodes = area.nodeIds.size();
        for (mcIdType nodeId : area.nodeIds)
        {
            std::array<double, 3> nodeCoords = nodes->getCoordinates(nodeId);
            for (size_t i = 0; i < 3; ++i) center[i] += nodeCoords[i] - area.radius * normals[3 * nodeId + i];
        }
        for (size_t i = 0; i < 3; ++i) center[i] /= (double)nbNodes;
    }
    area.center = center;
    area.radius = fabs(area.radius);
}

void
Areas::computeCylinderProperties(mcIdType areaId)
{
    Area &area = areas[areaId];
    size_t nbNodes = area.nodeIds.size();
    // Project the nodes to the central axis of the cylinder
    std::vector<double> projectedNodes(3 * nbNodes, 0.0);
    area.radius = 0;
    area.axisPoint.fill(0.0);
    for (size_t i = 0; i < nbNodes; ++i)
    {
        mcIdType nodeId = area.nodeIds[i];
        std::array<double, 3> nodeCoords = nodes->getCoordinates(nodeId);
        double approxRadius = 1.0 / nodes->getKdiff0(nodeId);
        for (size_t j = 0; j < 3; ++j)
        {
            projectedNodes[3 * i + j] = nodeCoords[j] - approxRadius * nodes->getNormals()[3 * nodeId + j];
            area.axisPoint[j] += projectedNodes[3 * i + j];
        }
        area.radius += approxRadius;
    }
    // Axis point is the mean of the projected nodes
    for (size_t i = 0; i < 3; ++i) area.axisPoint[i] /= (double)nbNodes;
    // Radius of the cylinder is the mean of the approximate radius of each node
    area.radius = fabs(area.radius / (double)nbNodes);
    // Compute the axis of the cylinder
    area.axis = MathOps::computePCAFirstAxis(projectedNodes);
}

void
Areas::computeConeProperties(mcIdType areaId)
{
    Area &area = areas[areaId];
    size_t nbNodes = area.nodeIds.size();
    // Project the nodes to the central axis of the cone
    std::vector<double> projectedNodes(3 * nbNodes, 0.0);
    std::vector<double> radiusNodes(nbNodes, 0.0);
    area.axisPoint.fill(0.0);
    for (size_t i = 0; i < nbNodes; ++i)
    {
        mcIdType nodeId = area.nodeIds[i];
        std::array<double, 3> nodeCoords = nodes->getCoordinates(nodeId);
        radiusNodes[i] = 1.0 / nodes->getKdiff0(nodeId);
        for (size_t j = 0; j < 3; ++j)
        {
            projectedNodes[3 * i + j] = nodeCoords[j] - radiusNodes[i] * nodes->getNormals()[3 * nodeId + j];
            area.axisPoint[j] += projectedNodes[3 * i + j];
        }
    }
    // Axis point is the mean of the projected nodes
    for (size_t i = 0; i < 3; ++i) area.axisPoint[i] /= (double)nbNodes;
    // Compute the axis of the cone
    area.axis = MathOps::computePCAFirstAxis(projectedNodes);
    double normAxis = MathOps::computeNorm(area.axis);
    for (size_t i = 0; i < 3; ++i) area.axis[i] /= normAxis;
    // Compute the angle of the cone
    const std::vector<double> &weakDirections = nodes->getWeakDirections();
    std::vector<double> weakDirectionNodes(3 * nbNodes, 0.0);
    for (size_t i = 0; i < nbNodes; ++i)
    {
        mcIdType nodeId = area.nodeIds[i];
        for (size_t j = 0; j < 3; ++j) weakDirectionNodes[3 * i + j] = weakDirections[3 * nodeId + j];
    }
    std::vector<double> angles = MathOps::computeAngles(weakDirectionNodes, area.axis);
    // Correct the angles > pi/2 with the supplementary angle
    for (size_t i = 0; i < angles.size(); ++i)
    {
        if (angles[i] > M_PI_2)
            angles[i] = M_PI - angles[i];
    }
    area.angle = MathOps::mean(angles);
    // Compute the radius of the cone which is the mean of the approximate radius of each node
    area.radius = MathOps::mean(radiusNodes);
    // Select extrem nodes
    double q1 = MathOps::computeQuantile(radiusNodes, 0.1);
    double q2 = MathOps::computeQuantile(radiusNodes, 0.9);
    std::vector<mcIdType> q1_indices;
    for (mcIdType idx = 0; idx < (mcIdType)radiusNodes.size(); ++idx)
    {
        if (radiusNodes[idx] < q1)
            q1_indices.push_back(idx);
    }
    std::vector<mcIdType> q2_indices;
    for (mcIdType idx = 0; idx < (mcIdType)radiusNodes.size(); ++idx)
    {
        if (radiusNodes[idx] > q2)
            q2_indices.push_back(idx);
    }
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(q2_indices.begin(), q2_indices.end(), g);
    // Compute the height of the cone
    // std::vector<double> heights(q1_indices.size(), 0.0);
    // std::vector<double> distancesToApex(q1_indices.size(), 0.0);
    std::array<double, 3> p{0.0, 0.0, 0.0};
    size_t min_q_size = std::min<size_t>(q1_indices.size(), q2_indices.size());
    for (size_t i = 0; i < min_q_size; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
            p[j] = (projectedNodes[3 * q1_indices[i] + j] - projectedNodes[3 * q2_indices[i] + j]);
        double height = MathOps::computeNorm(p);
        double distanceToApex =
            height * radiusNodes[q1_indices[i]] / (radiusNodes[q2_indices[i]] - radiusNodes[q1_indices[i]]);
        double orientationAxis = MathOps::dot(p, area.axis);
        for (size_t j = 0; j < 3; ++j)
        {
            area.apex[j] += projectedNodes[3 * q1_indices[i] + j];
            if (orientationAxis >= 0)
                area.apex[j] += distanceToApex * area.axis[j];
            else
                area.apex[j] -= distanceToApex * area.axis[j];
        }
    }
    for (size_t j = 0; j < 3; ++j) area.apex[j] /= (double)min_q_size;
}

void
Areas::computeTorusProperties(mcIdType areaId)
{
    Area &area = areas[areaId];
    size_t n = area.nodeIds.size();
    if (n == 0)
        return;
    std::vector<double> areaNodesK1(n, 0.0);
    std::vector<double> areaNodesK2(n, 0.0);
    for (size_t i = 0; i < n; ++i)
    {
        areaNodesK1[i] = nodes->getK1(area.nodeIds[i]);
        areaNodesK2[i] = nodes->getK2(area.nodeIds[i]);
    }
    double var1 = MathOps::computeVariance(areaNodesK1);
    double var2 = MathOps::computeVariance(areaNodesK2);
    double minorCurvature;
    if (var1 > var2)
        minorCurvature = MathOps::mean(areaNodesK2);
    else
        minorCurvature = MathOps::mean(areaNodesK1);
    area.minorRadius = 1.0 / minorCurvature;
    std::vector<double> majorRadiusNodes(3 * n, 0.0);
    for (size_t i = 0; i < n; ++i)
    {
        std::array<double, 3> coords = nodes->getCoordinates(area.nodeIds[i]);
        std::array<double, 3> normal = nodes->getNormal(area.nodeIds[i]);
        for (size_t j = 0; j < 3; ++j) majorRadiusNodes[3 * i + j] = coords[j] - normal[j] * area.minorRadius;
    }
    std::array<double, 3> meanMajorRadiusNodes = MathOps::meanCoordinates(majorRadiusNodes);
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < 3; ++j) majorRadiusNodes[3 * i + j] -= meanMajorRadiusNodes[j];
    }
    std::array<double, 3> normal = MathOps::computePCAThirdAxis(majorRadiusNodes);
    std::array<double, 6> base2d = MathOps::computeBaseFromNormal(normal);
    std::vector<double> projectedMajorRadiusNodes(2 * n, 0.0);
    cblas_dgemm(
        CBLAS_LAYOUT::CblasRowMajor,
        CBLAS_TRANSPOSE::CblasNoTrans,
        CBLAS_TRANSPOSE::CblasTrans,
        (int)n,
        2,
        3,
        1.0,
        majorRadiusNodes.data(),
        3,
        base2d.data(),
        3,
        0.0,
        projectedMajorRadiusNodes.data(),
        2
    );
    std::vector<double> A(3 * n, 1.0);
    std::vector<double> B(n, 0.0);
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < 2; ++j) A[3 * i + j] = projectedMajorRadiusNodes[2 * i + j];
        B[i] = pow(projectedMajorRadiusNodes[2 * i], 2) + pow(projectedMajorRadiusNodes[2 * i + 1], 2);
    }
    std::vector<double> fit = MathOps::lstsqRow(A, B);
    double a = fit[0];
    double b = fit[1];
    double c = fit[2];
    double xc = a / 2.0;
    double yc = b / 2.0;
    area.radius = sqrt(4.0 * c + pow(a, 2) + pow(b, 2)) / 2.0;
    for (size_t i = 0; i < 3; ++i) area.center[i] = xc * base2d[i] + yc * base2d[3 + i] + meanMajorRadiusNodes[i];
}

const std::vector<Int32> &
Areas::getAreaIdByNodes() const
{
    return areaIdByNodes;
}
