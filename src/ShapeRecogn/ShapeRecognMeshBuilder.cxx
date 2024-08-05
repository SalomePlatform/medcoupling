#include "ShapeRecognMeshBuilder.hxx"

#include "NodesBuilder.hxx"
#include "AreasBuilder.hxx"
#include "MEDLoader.hxx"
#include "ShapeRecognMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"

using namespace MEDCoupling;

ShapeRecognMeshBuilder::ShapeRecognMeshBuilder(const std::string &fileName, int meshDimRelToMax)
{
    mesh = ReadUMeshFromFile(fileName, meshDimRelToMax);
    if (mesh->getMeshDimension() != 2)
        throw INTERP_KERNEL::Exception("Expect a mesh with a dimension equal to 2");
    if (mesh->getNumberOfCellsWithType(INTERP_KERNEL::NORM_TRI3) != mesh->getNumberOfCells())
        throw INTERP_KERNEL::Exception("Expect a mesh containing exclusively triangular cells");
}

ShapeRecognMeshBuilder::~ShapeRecognMeshBuilder()
{
    if (areas != nullptr)
        delete areas;
    if (nodes != nullptr)
        delete nodes;
    mesh->decrRef();
}

ShapeRecognMesh *ShapeRecognMeshBuilder::recognize()
{
    mesh->incrRef();
    NodesBuilder nodesBuilder(mesh);
    nodes = nodesBuilder.build();
    AreasBuilder areasBuilder(nodes);
    areasBuilder.build();
    areas = areasBuilder.getAreas();
    MCAuto<ShapeRecognMesh> recognMesh = ShapeRecognMesh::New();
    recognMesh->nodeK1 = buildNodeK1();
    recognMesh->nodeK2 = buildNodeK2();
    recognMesh->nodePrimitiveType = buildNodePrimitiveType();
    recognMesh->nodeNormal = buildNodeNormal();
    recognMesh->areaId = buildAreaId();
    recognMesh->areaPrimitiveType = buildAreaPrimitiveType();
    recognMesh->areaNormal = buildAreaNormal();
    recognMesh->minorRadius = buildMinorRadius();
    recognMesh->radius = buildRadius();
    recognMesh->angle = buildAngle();
    recognMesh->center = buildCenter();
    recognMesh->axis = buildAxis();
    recognMesh->apex = buildApex();
    return recognMesh.retn();
}

const Nodes *ShapeRecognMeshBuilder::getNodes() const
{
    return nodes;
}

const Areas *ShapeRecognMeshBuilder::getAreas() const
{
    return areas;
}

MEDCouplingFieldDouble *ShapeRecognMeshBuilder::buildNodeK1() const
{
    if (nodes == nullptr)
        throw INTERP_KERNEL::Exception("recognize must be called before building any fields");
    return buildField("K1 (Node)", 1, nodes->getK1());
}

MEDCouplingFieldDouble *ShapeRecognMeshBuilder::buildNodeK2() const
{
    if (nodes == nullptr)
        throw INTERP_KERNEL::Exception("recognize must be called before building any fields");
    return buildField("K2 (Node)", 1, nodes->getK2());
}

MEDCouplingFieldDouble *ShapeRecognMeshBuilder::buildNodePrimitiveType() const
{
    if (nodes == nullptr)
        throw INTERP_KERNEL::Exception("recognize must be called before building any fields");
    return buildField("Primitive Type (Node)", 1, nodes->getPrimitiveType());
}

MEDCouplingFieldDouble *ShapeRecognMeshBuilder::buildNodeNormal() const
{
    if (nodes == nullptr)
        throw INTERP_KERNEL::Exception("recognize must be called before building any fields");
    return buildField("Normal (Node)", 3, nodes->getNormals());
}

MEDCouplingFieldDouble *ShapeRecognMeshBuilder::buildAreaId() const
{
    if (areas == nullptr)
        throw INTERP_KERNEL::Exception("recognize must be called before building any fields");
    return buildField("Area Id", 1, areas->getAreaIdByNodes());
}

MEDCouplingFieldDouble *ShapeRecognMeshBuilder::buildAreaPrimitiveType() const
{
    if (areas == nullptr)
        throw INTERP_KERNEL::Exception("recognize must be called before building any fields");
    double *values = buildAreaArray([](Areas *areas, mcIdType areaId) -> double
                                    { return (double)areas->getPrimitiveType(areaId); });
    return buildField("Primitive Type (Area)", 1, values);
}

MEDCouplingFieldDouble *ShapeRecognMeshBuilder::buildAreaNormal() const
{
    if (areas == nullptr)
        throw INTERP_KERNEL::Exception("recognize must be called before building any fields");
    double *values = buildArea3DArray([](Areas *areas, mcIdType areaId) -> const std::array<double, 3> &
                                      { return areas->getNormal(areaId); });
    return buildField("Normal (Area)", 3, values);
}

MEDCouplingFieldDouble *ShapeRecognMeshBuilder::buildMinorRadius() const
{
    if (areas == nullptr)
        throw INTERP_KERNEL::Exception("recognize must be called before building any fields");
    double *values = buildAreaArray([](Areas *areas, mcIdType areaId) -> double
                                    { return areas->getMinorRadius(areaId); });
    return buildField("Minor Radius (Area)", 1, values);
}

MEDCouplingFieldDouble *ShapeRecognMeshBuilder::buildRadius() const
{
    if (areas == nullptr)
        throw INTERP_KERNEL::Exception("recognize must be called before building any fields");
    double *values = buildAreaArray([](Areas *areas, mcIdType areaId) -> double
                                    { return areas->getRadius(areaId); });
    return buildField("Radius (Area)", 1, values);
}

MEDCouplingFieldDouble *ShapeRecognMeshBuilder::buildAngle() const
{
    if (areas == nullptr)
        throw INTERP_KERNEL::Exception("recognize must be called before building any fields");
    double *values = buildAreaArray([](Areas *areas, mcIdType areaId) -> double
                                    { return areas->getAngle(areaId); });
    return buildField("Angle (Area)", 1, values);
}

MEDCouplingFieldDouble *ShapeRecognMeshBuilder::buildCenter() const
{
    if (areas == nullptr)
        throw INTERP_KERNEL::Exception("recognize must be called before building any fields");
    double *values = buildArea3DArray([](Areas *areas, mcIdType areaId) -> const std::array<double, 3> &
                                      { return areas->getCenter(areaId); });
    return buildField("Center (Area)", 3, values);
}

MEDCouplingFieldDouble *ShapeRecognMeshBuilder::buildAxis() const
{
    if (areas == nullptr)
        throw INTERP_KERNEL::Exception("recognize must be called before building any fields");
    double *values = buildArea3DArray([](Areas *areas, mcIdType areaId) -> const std::array<double, 3> &
                                      { return areas->getAxis(areaId); });
    return buildField("Axis (Area)", 3, values);
}

MEDCouplingFieldDouble *ShapeRecognMeshBuilder::buildApex() const
{
    if (areas == nullptr)
        throw INTERP_KERNEL::Exception("recognize must be called before building any fields");
    double *values = buildArea3DArray([](Areas *areas, mcIdType areaId) -> const std::array<double, 3> &
                                      { return areas->getApex(areaId); });
    return buildField("Apex (Area)", 3, values);
}

template <typename T>
MEDCouplingFieldDouble *ShapeRecognMeshBuilder::buildField(
    const std::string &name,
    size_t nbOfCompo,
    const std::vector<T> &values) const
{
    DataArrayDouble *data = DataArrayDouble::New();
    data->setName(name);
    data->alloc(nodes->getNbNodes(), nbOfCompo);
    std::copy(values.begin(), values.end(), data->getPointer());
    data->declareAsNew();
    return buildField(name, nbOfCompo, data);
}

MEDCouplingFieldDouble *ShapeRecognMeshBuilder::buildField(
    const std::string &name,
    size_t nbOfCompo,
    double *values) const
{
    DataArrayDouble *data = DataArrayDouble::New();
    data->setName(name);
    data->useArray(
        values,
        true,
        MEDCoupling::DeallocType::CPP_DEALLOC,
        nodes->getNbNodes(),
        nbOfCompo);
    return buildField(name, nbOfCompo, data);
}

MEDCouplingFieldDouble *ShapeRecognMeshBuilder::buildField(
    const std::string &name,
    size_t nbOfCompo,
    DataArrayDouble *data) const
{
    MEDCouplingFieldDouble *field = MEDCouplingFieldDouble::New(ON_NODES);
    field->setName(name);
    field->setMesh(mesh);
    if (nbOfCompo == 3)
    {
        data->setInfoOnComponent(0, "X");
        data->setInfoOnComponent(1, "Y");
        data->setInfoOnComponent(2, "Z");
    }
    field->setArray(data);
    data->decrRef();
    return field;
}

double *ShapeRecognMeshBuilder::buildArea3DArray(
    const std::array<double, 3> &(*areaFunc)(Areas *, mcIdType)) const
{
    double *values = new double[3 * nodes->getNbNodes()];
    const std::vector<mcIdType> &areaIdByNodes = areas->getAreaIdByNodes();
    for (size_t nodeId = 0; nodeId < areaIdByNodes.size(); ++nodeId)
    {
        mcIdType areaId = areaIdByNodes[nodeId];
        if (areaId != -1)
        {
            const std::array<double, 3> &areaValues = areaFunc(areas, areaId);
            values[3 * nodeId] = areaValues[0];
            values[3 * nodeId + 1] = areaValues[1];
            values[3 * nodeId + 2] = areaValues[2];
        }
    }
    return values;
}

double *ShapeRecognMeshBuilder::buildAreaArray(double (*areaFunc)(Areas *, mcIdType)) const
{
    const std::vector<mcIdType> &areaIdByNodes = areas->getAreaIdByNodes();
    double *values = new double[nodes->getNbNodes()];
    for (size_t nodeId = 0; nodeId < areaIdByNodes.size(); ++nodeId)
    {
        mcIdType areaId = areaIdByNodes[nodeId];
        if (areaId != -1)
            values[nodeId] = areaFunc(areas, areaId);
    }
    return values;
}
