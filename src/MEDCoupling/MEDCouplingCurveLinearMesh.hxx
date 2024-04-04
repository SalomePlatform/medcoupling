// Copyright (C) 2007-2024  CEA, EDF
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
// Author : Anthony Geay (CEA/DEN)

#ifndef __PARAMEDMEM_MEDCOUPLINGCURVELINEARMESH_HXX__
#define __PARAMEDMEM_MEDCOUPLINGCURVELINEARMESH_HXX__

#include "MCType.hxx"
#include "MEDCoupling.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "MEDCouplingMesh.hxx"
#include "MEDCouplingStructuredMesh.hxx"
#include "MCAuto.hxx"
#include <string>
#include <cstddef>
#include <vector>
#include <utility>
#include <ostream>

namespace MEDCoupling
{
  class MEDCouplingCurveLinearMesh : public MEDCouplingStructuredMesh
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingCurveLinearMesh *New();
    MEDCOUPLING_EXPORT static MEDCouplingCurveLinearMesh *New(const std::string& meshName);
    MEDCOUPLING_EXPORT std::string getClassName() const override { return std::string("MEDCouplingCurveLinearMesh"); }
    MEDCOUPLING_EXPORT MEDCouplingCurveLinearMesh *deepCopy() const override;
    MEDCOUPLING_EXPORT MEDCouplingCurveLinearMesh *clone(bool recDeepCpy) const override;
    MEDCOUPLING_EXPORT void updateTime() const override;
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const override;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const override;
    MEDCOUPLING_EXPORT MEDCouplingMeshType getType() const override { return CURVE_LINEAR; }
    MEDCOUPLING_EXPORT void copyTinyStringsFrom(const MEDCouplingMesh *other) override;
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const override;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const override;
    MEDCOUPLING_EXPORT void checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                 DataArrayIdType *&cellCor, DataArrayIdType *&nodeCor) const override;
    MEDCOUPLING_EXPORT void checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                            DataArrayIdType *&cellCor) const override;
    MEDCOUPLING_EXPORT void checkConsistencyLight() const override;
    MEDCOUPLING_EXPORT void checkConsistency(double eps=1e-12) const override;
    MEDCOUPLING_EXPORT mcIdType getNumberOfCells() const override;
    MEDCOUPLING_EXPORT mcIdType getNumberOfNodes() const override;
    MEDCOUPLING_EXPORT int getSpaceDimension() const override;
    MEDCOUPLING_EXPORT void getCoordinatesOfNode(mcIdType nodeId, std::vector<double>& coo) const override;
    MEDCOUPLING_EXPORT std::string simpleRepr() const override;
    MEDCOUPLING_EXPORT std::string advancedRepr() const override;
    MEDCOUPLING_EXPORT const DataArrayDouble *getDirectAccessOfCoordsArrIfInStructure() const override;
    MEDCOUPLING_EXPORT DataArrayDouble *getCoords();
    MEDCOUPLING_EXPORT const DataArrayDouble *getCoords() const;
    MEDCOUPLING_EXPORT void setCoords(const DataArrayDouble *coords);
    MEDCOUPLING_EXPORT void setNodeGridStructure(const mcIdType *gridStructBg, const mcIdType *gridStructEnd);
    MEDCOUPLING_EXPORT std::vector<mcIdType> getNodeGridStructure() const override;
    MEDCOUPLING_EXPORT MEDCouplingStructuredMesh *buildStructuredSubPart(const std::vector< std::pair<mcIdType,mcIdType> >& cellPart) const override;
    // tools
    MEDCOUPLING_EXPORT void getBoundingBox(double *bbox) const override;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureField(bool isAbs) const override;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureFieldOnNode(bool isAbs) const override;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildOrthogonalField() const override;
    MEDCOUPLING_EXPORT mcIdType getCellContainingPoint(const double *pos, double eps) const override;
    MEDCOUPLING_EXPORT void getCellsContainingPoint(const double *pos, double eps, std::vector<mcIdType>& elts) const override;
    MEDCOUPLING_EXPORT void rotate(const double *center, const double *vector, double angle) override;
    MEDCOUPLING_EXPORT void translate(const double *vector) override;
    MEDCOUPLING_EXPORT void scale(const double *point, double factor) override;
    MEDCOUPLING_EXPORT MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const override;
    MEDCOUPLING_EXPORT DataArrayDouble *getCoordinatesAndOwner() const override;
    MEDCOUPLING_EXPORT DataArrayDouble *computeCellCenterOfMass() const override;
    MEDCOUPLING_EXPORT DataArrayDouble *computeIsoBarycenterOfNodesPerCell() const override;
    MEDCOUPLING_EXPORT void renumberCells(const mcIdType *old2NewBg, bool check=true) override;
    //some useful methods
    MEDCOUPLING_EXPORT void getNodeGridStructure(mcIdType *res) const override;
    //serialisation-unserialization
    MEDCOUPLING_EXPORT void getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<mcIdType>& tinyInfo, std::vector<std::string>& littleStrings) const override;
    MEDCOUPLING_EXPORT void resizeForUnserialization(const std::vector<mcIdType>& tinyInfo, DataArrayIdType *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const override;
    MEDCOUPLING_EXPORT void serialize(DataArrayIdType *&a1, DataArrayDouble *&a2) const override;
    MEDCOUPLING_EXPORT void unserialization(const std::vector<double>& tinyInfoD, const std::vector<mcIdType>& tinyInfo, const DataArrayIdType *a1, DataArrayDouble *a2,
                                            const std::vector<std::string>& littleStrings) override;
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const override;
    MEDCOUPLING_EXPORT std::string getVTKFileExtension() const override;
  private:
    void getMeasureFieldMeshDim1(bool isAbs, MEDCouplingFieldDouble *field) const;
    void getMeasureFieldMeshDim2(bool isAbs, MEDCouplingFieldDouble *field) const;
    void getMeasureFieldMeshDim3(bool isAbs, MEDCouplingFieldDouble *field) const;
    void getBarycenterAndOwnerMeshDim3(DataArrayDouble *bary) const;
    void getBarycenterAndOwnerMeshDim2(DataArrayDouble *bary) const;
    void getBarycenterAndOwnerMeshDim1(DataArrayDouble *bary) const;
  private:
    MEDCouplingCurveLinearMesh();
    MEDCouplingCurveLinearMesh(const MEDCouplingCurveLinearMesh& other, bool deepCpy);
    ~MEDCouplingCurveLinearMesh() override;
    void writeVTKLL(std::ostream& ofs, const std::string& cellData, const std::string& pointData, DataArrayByte *byteData) const override;
    std::string getVTKDataSetType() const override;
  private:
    MCAuto<DataArrayDouble> _coords;
    std::vector<mcIdType> _structure;
  };
}

#endif
