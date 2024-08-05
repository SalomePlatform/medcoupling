#include "ShapeRecognMesh.hxx"

#include "MEDLoader.hxx"

using namespace MEDCoupling;

ShapeRecognMesh::ShapeRecognMesh()
    : nodeK1(0), nodeK2(0), nodePrimitiveType(0),
      nodeNormal(0), areaId(0), areaPrimitiveType(0),
      areaNormal(0), minorRadius(0), radius(0),
      angle(0), center(0), axis(0), apex(0)
{
}
ShapeRecognMesh::~ShapeRecognMesh()
{
    nodeK1->decrRef();
    nodeK2->decrRef();
    nodePrimitiveType->decrRef();
    nodeNormal->decrRef();
    areaId->decrRef();
    areaPrimitiveType->decrRef();
    areaNormal->decrRef();
    minorRadius->decrRef();
    radius->decrRef();
    angle->decrRef();
    center->decrRef();
    axis->decrRef();
    apex->decrRef();
}

std::size_t ShapeRecognMesh::getHeapMemorySizeWithoutChildren() const
{
    return 0;
}

std::vector<const BigMemoryObject *> ShapeRecognMesh::getDirectChildrenWithNull() const
{
    std::vector<const BigMemoryObject *> ret;
    ret.push_back(nodeK1);
    ret.push_back(nodeK2);
    ret.push_back(nodePrimitiveType);
    ret.push_back(nodeNormal);
    ret.push_back(areaId);
    ret.push_back(areaPrimitiveType);
    ret.push_back(areaNormal);
    ret.push_back(minorRadius);
    ret.push_back(radius);
    ret.push_back(angle);
    ret.push_back(center);
    ret.push_back(axis);
    ret.push_back(apex);
    return ret;
}

ShapeRecognMesh *ShapeRecognMesh::New()
{
    return new ShapeRecognMesh;
}

void ShapeRecognMesh::save(const std::string &outputFile, bool writeFromScratch) const
{
    // Nodes
    // - k1
    WriteField(outputFile, nodeK1, writeFromScratch);
    // - k2
    WriteField(outputFile, nodeK2, false);
    // - primitive types
    WriteField(outputFile, nodePrimitiveType, false);
    // - Normal
    WriteField(outputFile, nodeNormal, false);
    // Areas
    // - Area Id
    WriteField(outputFile, areaId, false);
    // - Primitive Types
    WriteField(outputFile, areaPrimitiveType, false);
    // - Normal
    WriteField(outputFile, areaNormal, false);
    // - Minor Radius
    WriteField(outputFile, minorRadius, false);
    // - Radius
    WriteField(outputFile, radius, false);
    // - Angle
    WriteField(outputFile, angle, false);
    // - Center
    WriteField(outputFile, center, false);
    // - Axis
    WriteField(outputFile, axis, false);
    // - Apex
    WriteField(outputFile, apex, false);
}

MEDCouplingFieldDouble *ShapeRecognMesh::getNodeK1() const
{
    return nodeK1;
}

MEDCouplingFieldDouble *ShapeRecognMesh::getNodeK2() const
{
    return nodeK2;
}

MEDCouplingFieldDouble *ShapeRecognMesh::getNodePrimitiveType() const
{
    return nodePrimitiveType;
}

MEDCouplingFieldDouble *ShapeRecognMesh::getNodeNormal() const
{
    return nodeNormal;
}

MEDCouplingFieldDouble *ShapeRecognMesh::getAreaId() const
{
    return areaId;
}

MEDCouplingFieldDouble *ShapeRecognMesh::getAreaPrimitiveType() const
{
    return areaPrimitiveType;
}

MEDCouplingFieldDouble *ShapeRecognMesh::getAreaNormal() const
{
    return areaNormal;
}

MEDCouplingFieldDouble *ShapeRecognMesh::getMinorRadius() const
{
    return minorRadius;
}

MEDCouplingFieldDouble *ShapeRecognMesh::getRadius() const
{
    return radius;
}

MEDCouplingFieldDouble *ShapeRecognMesh::getAngle() const
{
    return angle;
}

MEDCouplingFieldDouble *ShapeRecognMesh::getCenter() const
{
    return center;
}

MEDCouplingFieldDouble *ShapeRecognMesh::getAxis() const
{
    return axis;
}

MEDCouplingFieldDouble *ShapeRecognMesh::getApex() const
{
    return apex;
}
