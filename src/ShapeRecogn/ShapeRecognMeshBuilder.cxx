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

#include "ShapeRecognMeshBuilder.hxx"

#include "NodesBuilder.hxx"
#include "AreasBuilder.hxx"
#include "ShapeRecognMesh.hxx"

#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldInt64.hxx"
#include "MEDCouplingFieldInt32.hxx"
#include "MEDCouplingMemArray.txx"

#include <algorithm>
#include <iterator>
#include <limits>

using namespace MEDCoupling;

ShapeRecognMeshBuilder::ShapeRecognMeshBuilder(MCAuto<MEDCouplingUMesh> mesh) { assign(mesh); }

ShapeRecognMeshBuilder::ShapeRecognMeshBuilder(MEDCouplingUMesh *mesh)
{
    MCAuto<MEDCouplingUMesh> meshIn = MCAuto<MEDCouplingUMesh>::TakeRef(mesh);
    assign(meshIn);
}

void
ShapeRecognMeshBuilder::assign(MCAuto<MEDCouplingUMesh> mesh)
{
    this->mesh = mesh;
    if (mesh->getMeshDimension() != 2)
        throw INTERP_KERNEL::Exception("Expect a mesh with a dimension equal to 2");
    if (mesh->getNumberOfCellsWithType(INTERP_KERNEL::NORM_TRI3) != mesh->getNumberOfCells())
        throw INTERP_KERNEL::Exception("Expect a mesh containing exclusively triangular cells");
}

MCAuto<ShapeRecognMesh>
ShapeRecognMeshBuilder::recognize()
{
    mesh->incrRef();
    NodesBuilder nodesBuilder(mesh);
    nodes.reset(nodesBuilder.build());
    AreasBuilder areasBuilder(nodes.get());
    areasBuilder.build();
    areas.reset(areasBuilder.getAreas());
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
    return recognMesh;
}

const Nodes *
ShapeRecognMeshBuilder::getNodes() const
{
    return nodes.get();
}

const Areas *
ShapeRecognMeshBuilder::getAreas() const
{
    return areas.get();
}

void
ShapeRecognMeshBuilder::checkNodesBeforeBuildingField() const
{
    if (!nodes.get())
        throw INTERP_KERNEL::Exception("recognize must be called before building any fields");
}

void
ShapeRecognMeshBuilder::checkAreasBeforeBuildingField() const
{
    if (!areas.get())
        throw INTERP_KERNEL::Exception("recognize must be called before building any fields");
}

template <typename T>
typename Traits<T>::FieldType *
buildField(const std::string &name, MCAuto<typename Traits<T>::ArrayType> data, const MEDCouplingUMesh *mesh)
{
    using ZeField = typename Traits<T>::FieldType;
    ZeField *field = ZeField::New(ON_NODES);
    mcIdType nbNodes = data->getNumberOfTuples();
    size_t nbOfCompo = data->getNumberOfComponents();
    field->setName(name);
    field->setMesh(mesh);
    if (nbOfCompo == 3)
    {
        data->setInfoOnComponent(0, "X");
        data->setInfoOnComponent(1, "Y");
        data->setInfoOnComponent(2, "Z");
    }
    field->setArray(data);
    return field;
}

template <typename T>
typename Traits<T>::FieldType *
buildField(const std::string &name, size_t nbOfCompo, const std::vector<T> &values, const MEDCouplingUMesh *mesh)
{
    using ZeArray = typename Traits<T>::ArrayType;
    mcIdType nbNodes = mesh->getNumberOfNodes();
    MCAuto<ZeArray> data(ZeArray::New());
    data->setName(name);
    data->alloc(nbNodes, nbOfCompo);
    std::copy(values.begin(), values.end(), data->getPointer());
    data->declareAsNew();
    return buildField<T>(name, data, mesh);
}

template <typename T>
typename Traits<T>::FieldType *
buildField(const std::string &name, size_t nbOfCompo, T *values, const MEDCouplingUMesh *mesh)
{
    using ZeArray = typename Traits<T>::ArrayType;
    mcIdType nbNodes = mesh->getNumberOfNodes();
    ZeArray *data = ZeArray::New();
    data->setName(name);
    data->useArray(values, true, MEDCoupling::DeallocType::CPP_DEALLOC, nbNodes, nbOfCompo);
    return buildField<T>(name, data, mesh);
}

template <class T>
T *
buildAreaArrayT(Areas *areas, Nodes *nodes, std::function<T(Areas *, mcIdType)> areaFunc)
{
    const std::vector<Int32> &areaIdByNodes = areas->getAreaIdByNodes();
    T *values = new T[nodes->getNbNodes()];
    for (size_t nodeId = 0; nodeId < areaIdByNodes.size(); ++nodeId)
    {
        Int32 areaId = areaIdByNodes[nodeId];
        if (areaId != -1)
            values[nodeId] = areaFunc(areas, areaId);
        else
            values[nodeId] = std::numeric_limits<T>::max();
    }
    return values;
}

MEDCouplingFieldDouble *
ShapeRecognMeshBuilder::buildNodeK1() const
{
    checkNodesBeforeBuildingField();
    return buildField<double>("K1 (Node)", 1, nodes->getK1(), mesh);
}

MEDCouplingFieldDouble *
ShapeRecognMeshBuilder::buildNodeK2() const
{
    checkNodesBeforeBuildingField();
    return buildField<double>("K2 (Node)", 1, nodes->getK2(), mesh);
}

MEDCouplingFieldInt32 *
ShapeRecognMeshBuilder::buildNodePrimitiveType() const
{
    checkNodesBeforeBuildingField();
    std::vector<Int32> tmp;
    std::transform(
        nodes->getPrimitiveType().begin(),
        nodes->getPrimitiveType().end(),
        std::back_inserter(tmp),
        [](const PrimitiveType &elt) { return Int32(elt); }
    );
    return buildField<Int32>("Primitive Type (Node)", 1, tmp, mesh);
}

MEDCouplingFieldDouble *
ShapeRecognMeshBuilder::buildNodeNormal() const
{
    checkNodesBeforeBuildingField();
    return buildField<double>("Normal (Node)", 3, nodes->getNormals(), mesh);
}

MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble>
ShapeRecognMeshBuilder::buildNodeWeakDirections() const
{
    checkNodesBeforeBuildingField();
    return MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble>(
        buildField<double>("WeakDirection (Node)", 3, nodes->getWeakDirections(), mesh)
    );
}

MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble>
ShapeRecognMeshBuilder::buildNodeMainDirections() const
{
    checkNodesBeforeBuildingField();
    return MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble>(
        buildField<double>("MainDirection (Node)", 3, nodes->getMainDirections(), mesh)
    );
}

MEDCouplingFieldInt32 *
ShapeRecognMeshBuilder::buildAreaId() const
{
    checkAreasBeforeBuildingField();
    return buildField<Int32>("Area Id", 1, areas->getAreaIdByNodes(), mesh);
}

MEDCouplingFieldInt32 *
ShapeRecognMeshBuilder::buildAreaPrimitiveType() const
{
    checkAreasBeforeBuildingField();
    Int32 *values = buildAreaArrayT<Int32>(
        areas.get(),
        nodes.get(),
        [](Areas *areas, mcIdType areaId) -> Int32 { return (Int32)areas->getPrimitiveType(areaId); }
    );
    return buildField<Int32>("Primitive Type (Area)", 1, values, mesh);
}

MEDCouplingFieldDouble *
ShapeRecognMeshBuilder::buildAreaNormal() const
{
    checkAreasBeforeBuildingField();
    double *values = buildArea3DArray(
        [](Areas *areas, mcIdType areaId) -> const std::array<double, 3> & { return areas->getNormal(areaId); }
    );
    return buildField<double>("Normal (Area)", 3, values, mesh);
}

MEDCouplingFieldDouble *
ShapeRecognMeshBuilder::buildMinorRadius() const
{
    checkAreasBeforeBuildingField();
    double *values =
        buildAreaArray([](Areas *areas, mcIdType areaId) -> double { return areas->getMinorRadius(areaId); });
    return buildField<double>("Minor Radius (Area)", 1, values, mesh);
}

MEDCouplingFieldDouble *
ShapeRecognMeshBuilder::buildRadius() const
{
    checkAreasBeforeBuildingField();
    double *values = buildAreaArray([](Areas *areas, mcIdType areaId) -> double { return areas->getRadius(areaId); });
    return buildField<double>("Radius (Area)", 1, values, mesh);
}

MEDCouplingFieldDouble *
ShapeRecognMeshBuilder::buildAngle() const
{
    checkAreasBeforeBuildingField();
    double *values = buildAreaArray([](Areas *areas, mcIdType areaId) -> double { return areas->getAngle(areaId); });
    return buildField<double>("Angle (Area)", 1, values, mesh);
}

MEDCouplingFieldDouble *
ShapeRecognMeshBuilder::buildCenter() const
{
    checkAreasBeforeBuildingField();
    double *values = buildArea3DArray(
        [](Areas *areas, mcIdType areaId) -> std::array<double, 3> { return areas->getCenter(areaId); }
    );
    return buildField<double>("Center (Area)", 3, values, mesh);
}

MEDCouplingFieldDouble *
ShapeRecognMeshBuilder::buildAxis() const
{
    checkAreasBeforeBuildingField();
    double *values =
        buildArea3DArray([](Areas *areas, mcIdType areaId) -> std::array<double, 3> { return areas->getAxis(areaId); });
    return buildField<double>("Axis (Area)", 3, values, mesh);
}

MEDCouplingFieldDouble *
ShapeRecognMeshBuilder::buildApex() const
{
    checkAreasBeforeBuildingField();
    double *values =
        buildArea3DArray([](Areas *areas, mcIdType areaId) -> std::array<double, 3> { return areas->getApex(areaId); });
    return buildField<double>("Apex (Area)", 3, values, mesh);
}

MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble>
ShapeRecognMeshBuilder::buildAreaAxisPoint() const
{
    checkAreasBeforeBuildingField();
    double *values = buildArea3DArray(
        [](Areas *areas, mcIdType areaId) -> std::array<double, 3> { return areas->getAxisPoint(areaId); }
    );
    return MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble>(
        buildField<double>("AxisPoint (Area)", 3, values, mesh)
    );
}

MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble>
ShapeRecognMeshBuilder::buildAreaAffinePoint() const
{
    checkAreasBeforeBuildingField();
    double *values = buildArea3DArray(
        [](Areas *areas, mcIdType areaId) -> std::array<double, 3> { return areas->getAffinePoint(areaId); }
    );
    return MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble>(
        buildField<double>("AffinePoint (Area)", 3, values, mesh)
    );
}

double *
ShapeRecognMeshBuilder::buildArea3DArray(std::function<std::array<double, 3>(Areas *, Int32)> areaFunc) const
{
    double *values = new double[3 * nodes->getNbNodes()];
    const std::vector<Int32> &areaIdByNodes = areas->getAreaIdByNodes();
    for (size_t nodeId = 0; nodeId < areaIdByNodes.size(); ++nodeId)
    {
        Int32 areaId = areaIdByNodes[nodeId];
        if (areaId != -1)
        {
            const std::array<double, 3> areaValues(areaFunc(areas.get(), areaId));
            values[3 * nodeId] = areaValues[0];
            values[3 * nodeId + 1] = areaValues[1];
            values[3 * nodeId + 2] = areaValues[2];
        }
        else
        {
            std::for_each(
                values + 3 * nodeId,
                values + 3 * (nodeId + 1),
                [](double &val) { val = std::numeric_limits<double>::max(); }
            );
        }
    }
    return values;
}

double *
ShapeRecognMeshBuilder::buildAreaArray(std::function<double(Areas *, Int32)> areaFunc) const
{
    return buildAreaArrayT<double>(areas.get(), nodes.get(), areaFunc);
}
