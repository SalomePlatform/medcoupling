#include "AreasBuilder.hxx"
#include "MathOps.hxx"
#include "ShapeRecongConstants.hxx"

#include <unordered_set>

using namespace MEDCoupling;

AreasBuilder::AreasBuilder(const Nodes *nodes) : nodes(nodes), areas(new Areas(nodes))
{
    size_t nbNodes = nodes->getNbNodes();
    threshold = std::max<size_t>(THRESHOLD_MIN_NB_NODES, nbNodes / THRESHOLD_MAX_NB_AREAS);
}

void AreasBuilder::build()
{
    explore();
    expand();
    rebuild();
}

void AreasBuilder::explore()
{
    exploreAreas();
    areas->cleanInvalidNodeAreas();
    filterHighPass();
}

void AreasBuilder::expand()
{
    expandAreas();
    filterHighPass();
}

void AreasBuilder::rebuild()
{
    rebuildInvalidAreas();
    filterHighPass();
    expandAreasByType(PrimitiveType::Cone);
    filterHighPass();
}

Areas *AreasBuilder::getAreas() const
{
    return areas;
}

void AreasBuilder::exploreAreas()
{
    // int nbNodesExplored = 0;
    std::vector<bool> exploredNodeIds(nodes->getNbNodes(), false);
    std::unordered_set<mcIdType> nodesToExplore;
    // Reserve a set with the size of nodes to avoid reallocation for each insert/erase
    // TODO: Improve the size ? Nb of Nodes is too much ?
    nodesToExplore.reserve(nodes->getNbNodes());
    mcIdType areaId = -1;
    for (mcIdType nodeId = 0; nodeId < nodes->getNbNodes(); ++nodeId)
    {
        if (!exploredNodeIds[nodeId] &&
            nodes->getPrimitiveType(nodeId) != PrimitiveType::Unknown)
        {
            exploredNodeIds[nodeId] = true;
            if (areaId != -1 && areas->getNumberOfNodes(areaId) < threshold)
                areas->cancelArea(areaId, nodes->getPrimitiveType(nodeId));
            else
                areaId = areas->addArea(nodes->getPrimitiveType(nodeId));
            areas->addNode(areaId, nodeId);
            const std::vector<mcIdType> neighbors = nodes->getNeighbors(nodeId);
            for (mcIdType neighborId : neighbors)
            {
                if (nodes->getPrimitiveType(neighborId) == areas->getPrimitiveType(areaId) &&
                    areas->getAreaId(neighborId) <= -1)
                    nodesToExplore.insert(neighborId);
            }
            // Explore all the neighbors matching the area
            while (!nodesToExplore.empty())
            {
                mcIdType neighborId = *nodesToExplore.begin();
                nodesToExplore.erase(neighborId);
                if (doesItMatch(areaId, neighborId))
                {
                    exploredNodeIds[neighborId] = true;
                    areas->addNode(areaId, neighborId);
                    const std::vector<mcIdType> neighborsOfNeighbor = nodes->getNeighbors(neighborId);
                    for (mcIdType neighborIdOfNeighbor : neighborsOfNeighbor)
                    {
                        if (!exploredNodeIds[neighborIdOfNeighbor] &&
                            areas->getAreaId(neighborIdOfNeighbor) <= -1 &&
                            // Already in doesItMatch but avoid useless insertion
                            nodes->getPrimitiveType(neighborIdOfNeighbor) == areas->getPrimitiveType(areaId))
                            nodesToExplore.insert(neighborIdOfNeighbor);
                    }
                }
            }
        }
        // if (!exploredNodeIds[nodeId])
        //     nbNodesExplored += 1;
        exploredNodeIds[nodeId] = true;
    }
}

void AreasBuilder::expandAreas()
{
    // Expand by topological order
    expandAreasByType(PrimitiveType::Plane);
    expandAreasByType(PrimitiveType::Sphere);
    expandAreasByType(PrimitiveType::Cylinder);
    expandAreasByType(PrimitiveType::Cone);
    expandAreasByType(PrimitiveType::Torus);
}

void AreasBuilder::expandAreasByType(PrimitiveType primitive)
{
    std::unordered_set<mcIdType> nodesToExplore;
    // Reserve a set with the size of nodes to avoid reallocation for each insert/erase
    // TODO: Improve the size ? Nb of Nodes is too much ?
    nodesToExplore.reserve(nodes->getNbNodes());
    for (mcIdType areaId = 0; areaId < (mcIdType)areas->getNumberOfAreas(); ++areaId)
    {
        if (areas->getPrimitiveType(areaId) == primitive)
        {
            std::vector<bool> exploredNodeIds(nodes->getNbNodes(), false);
            areas->computeProperties(areaId);
            const std::vector<mcIdType> &nodeIds = areas->getNodeIds(areaId);
            for (mcIdType nodeId : nodeIds)
            {
                exploredNodeIds[nodeId] = true;
                nodesToExplore.insert(nodeId);
            }
            while (!nodesToExplore.empty())
            {
                mcIdType nodeId = *nodesToExplore.begin();
                nodesToExplore.erase(nodeId);
                if (doesItBelong(areaId, nodeId))
                {
                    // TODO: Is the properties need to be updated after adding a node ?
                    // It gives bad results for the cone and the cylinder
                    // mcIdType oldAreaId = areas->getAreaId(nodeId);
                    areas->addNode(areaId, nodeId);
                    // areas->computeProperties(areaId);
                    // areas->computeProperties(oldAreaId);
                    const std::vector<mcIdType> neighborIds = nodes->getNeighbors(nodeId);
                    for (mcIdType neighborId : neighborIds)
                    {
                        if (!exploredNodeIds[neighborId])
                            nodesToExplore.insert(neighborId);
                    }
                }
                exploredNodeIds[nodeId] = true;
            }
        }
    }
}

void AreasBuilder::rebuildInvalidAreas()
{
    std::vector<mcIdType> exploredNodeIds(nodes->getNbNodes(), false);
    std::vector<bool> isIinvalidNodes(nodes->getNbNodes(), false);
    std::unordered_set<mcIdType> nodesToExplore;
    // Reserve a set with the size of nodes to avoid reallocation for each insert/erase
    // TODO: Improve the size ? Nb of Nodes is too much ?
    nodesToExplore.reserve(nodes->getNbNodes());
    for (mcIdType nodeId = 0; nodeId < nodes->getNbNodes(); ++nodeId)
        isIinvalidNodes[nodeId] = isInvalidCylinderNode(nodeId);
    for (mcIdType nodeId = 0; nodeId < nodes->getNbNodes(); ++nodeId)
    {
        if (isIinvalidNodes[nodeId] && !exploredNodeIds[nodeId])
        {
            mcIdType areaId = areas->addArea(PrimitiveType::Cone);
            areas->addNode(areaId, nodeId);
            const std::vector<mcIdType> neighbors = nodes->getNeighbors(nodeId);
            for (mcIdType neighborId : neighbors)
                nodesToExplore.insert(neighborId);
            while (!nodesToExplore.empty())
            {
                mcIdType neighborId = *nodesToExplore.begin();
                nodesToExplore.erase(neighborId);
                if (
                    (areas->getAreaId(neighborId) == -1 ||
                     isIinvalidNodes[neighborId]))
                {
                    exploredNodeIds[neighborId] = true;
                    areas->addNode(areaId, neighborId);
                    const std::vector<mcIdType> neighborsOfNeighbor = nodes->getNeighbors(neighborId);
                    for (mcIdType neighborIdOfNeighbor : neighborsOfNeighbor)
                    {
                        if (!exploredNodeIds[neighborIdOfNeighbor])
                            nodesToExplore.insert(neighborIdOfNeighbor);
                    }
                }
            }
        }
    }
}

void AreasBuilder::filterHighPass()
{
    mcIdType nbAreas = areas->getNumberOfAreas();
    for (mcIdType areaId = (nbAreas - 1); areaId >= 0; --areaId)
    {
        if (areas->getNumberOfNodes(areaId) < threshold)
            areas->removeArea(areaId);
    }
}

bool AreasBuilder::doesItMatch(mcIdType areaId, mcIdType nodeId) const
{
    PrimitiveType areaType = areas->getPrimitiveType(areaId);
    PrimitiveType neighborType = nodes->getPrimitiveType(nodeId);
    bool isMatching = false;
    if (areaType == neighborType)
    {
        switch (areaType)
        {
        case PrimitiveType::Plane:
        case PrimitiveType::Torus:
            isMatching = true;
            break;
        case PrimitiveType::Sphere:
        {
            double kmoy = fabs(
                (areas->getAdimK1(areaId) + areas->getAdimK2(areaId)) / 2.0);
            double nodeKmoy = fabs(
                (nodes->getAdimK1(nodeId) + nodes->getAdimK2(nodeId)) / 2.0);
            isMatching = fabs((nodeKmoy - kmoy)) < TOL_MATCH_SPHERE;
        }
        break;
        case PrimitiveType::Cylinder:
            isMatching = fabs(areas->getAdimKdiff0(areaId) -
                              nodes->getAdimKdiff0(nodeId)) < TOL_MATCH_CYLINDER;
            break;
        case PrimitiveType::Cone:
        case PrimitiveType::Unknown:
        default:
            break;
        }
    }
    return isMatching;
}

bool AreasBuilder::doesItBelong(mcIdType areaId, mcIdType nodeId) const
{
    bool isClose = false;
    if (areas->isNodeCompatible(areaId, nodeId))
    {
        PrimitiveType areaType = areas->getPrimitiveType(areaId);
        switch (areaType)
        {
        case PrimitiveType::Plane:
        {
            isClose = distanceToPlane(
                          nodes->getCoordinates(nodeId),
                          areas->getAffinePoint(areaId),
                          areas->getNormal(areaId)) < DELTA_PLANE;
        }
        break;
        case PrimitiveType::Sphere:
        {
            double distanceToCenter = distanceToSphere(
                nodes->getCoordinates(nodeId),
                areas->getCenter(areaId));
            double radius = areas->getRadius(areaId);
            isClose = fabs((distanceToCenter - radius) / radius) < DELTA_SPHERE;
        }
        break;
        case PrimitiveType::Cylinder:
        {
            double distance = distanceToCylinder(
                nodes->getCoordinates(nodeId),
                areas->getAxis(areaId),
                areas->getAxisPoint(areaId));
            double radius = areas->getRadius(areaId);
            isClose = fabs((distance - radius) / radius) < DELTA_CYLINDER;
        }
        break;
        case PrimitiveType::Cone:
        {
            double radius = areas->getRadius(areaId);
            isClose = distanceToCone(
                          nodes->getCoordinates(nodeId),
                          areas->getAxis(areaId),
                          areas->getApex(areaId),
                          areas->getAngle(areaId)) /
                          fabs(radius) <
                      DELTA_CONE;
        }
        break;
        case PrimitiveType::Torus:
        case PrimitiveType::Unknown:
        default:
            break;
        }
    }
    return isClose;
}

bool AreasBuilder::isInvalidCylinderNode(mcIdType nodeId) const
{
    mcIdType areaId = areas->getAreaId(nodeId);
    if (areaId != -1 &&
        nodes->getPrimitiveType(nodeId) == PrimitiveType::Cylinder &&
        areas->getPrimitiveType(areaId) == PrimitiveType::Cylinder)
    {
        double angle = MathOps::computeAngle(
            nodes->getWeakDirection(nodeId),
            areas->getAxis(areaId));
        return angle >= THETA_MAX_CYLINDER && angle <= (M_PI - THETA_MAX_CYLINDER);
    }
    return false;
}

double AreasBuilder::distanceToPlane(
    const std::array<double, 3> &nodeCoords,
    const std::array<double, 3> &point,
    const std::array<double, 3> &normal)
{
    std::array<double, 3> vec{0.0, 0.0, 0.0};
    for (size_t i = 0; i < nodeCoords.size(); ++i)
        vec[i] = nodeCoords[i] - point[i];
    return fabs(MathOps::dot(vec, normal)) / MathOps::computeNorm(normal);
}

double AreasBuilder::distanceToSphere(
    const std::array<double, 3> &nodeCoords,
    const std::array<double, 3> &center)
{
    return MathOps::computeNorm(
        std::array<double, 3>{
            nodeCoords[0] - center[0],
            nodeCoords[1] - center[1],
            nodeCoords[2] - center[2]});
}

double AreasBuilder::distanceToCylinder(
    const std::array<double, 3> &nodeCoords,
    const std::array<double, 3> &axis,
    const std::array<double, 3> &axisPoint)
{

    std::array<double, 3> pa = {
        axisPoint[0] - nodeCoords[0],
        axisPoint[1] - nodeCoords[1],
        axisPoint[2] - nodeCoords[2]};
    double innerProduct = MathOps::dot(pa, axis);
    return MathOps::computeNorm(
        std::array<double, 3>({pa[0] - innerProduct * axis[0],
                               pa[1] - innerProduct * axis[1],
                               pa[2] - innerProduct * axis[2]}));
}

double AreasBuilder::distanceToCone(
    const std::array<double, 3> &nodeCoords,
    const std::array<double, 3> &axis,
    const std::array<double, 3> &apex,
    double angle)
{
    std::array<double, 3> ps{
        apex[0] - nodeCoords[0],
        apex[1] - nodeCoords[1],
        apex[2] - nodeCoords[2]};
    std::array<double, 3> v(axis);
    if (MathOps::dot(axis, ps) <= 0)
    {
        for (size_t i = 0; i < 3; ++i)
            v[i] *= -1;
    }
    double a = MathOps::computeNorm(
        MathOps::cross(ps, v));
    double b = MathOps::dot(ps, v);
    return fabs(a * cos(angle) - b * sin(angle));
}
