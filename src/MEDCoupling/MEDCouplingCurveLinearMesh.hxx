// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
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

#include "MEDCoupling.hxx"
#include "MEDCouplingStructuredMesh.hxx"
#include "MCAuto.hxx"

namespace MEDCoupling
{
  class MEDCouplingCurveLinearMesh : public MEDCouplingStructuredMesh
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingCurveLinearMesh *New();
    MEDCOUPLING_EXPORT static MEDCouplingCurveLinearMesh *New(const std::string& meshName);
    MEDCOUPLING_EXPORT MEDCouplingCurveLinearMesh *deepCopy() const;
    MEDCOUPLING_EXPORT MEDCouplingCurveLinearMesh *clone(bool recDeepCpy) const;
    MEDCOUPLING_EXPORT void updateTime() const;
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDCOUPLING_EXPORT MEDCouplingMeshType getType() const { return CURVE_LINEAR; }
    MEDCOUPLING_EXPORT void copyTinyStringsFrom(const MEDCouplingMesh *other);
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const;
    MEDCOUPLING_EXPORT void checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                 DataArrayInt *&cellCor, DataArrayInt *&nodeCor) const;
    MEDCOUPLING_EXPORT void checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                            DataArrayInt *&cellCor) const;
    MEDCOUPLING_EXPORT void checkConsistencyLight() const;
    MEDCOUPLING_EXPORT void checkConsistency(double eps=1e-12) const;
    MEDCOUPLING_EXPORT std::size_t getNumberOfCells() const;
    MEDCOUPLING_EXPORT int getNumberOfNodes() const;
    MEDCOUPLING_EXPORT int getSpaceDimension() const;
    MEDCOUPLING_EXPORT void getCoordinatesOfNode(int nodeId, std::vector<double>& coo) const;
    MEDCOUPLING_EXPORT std::string simpleRepr() const;
    MEDCOUPLING_EXPORT std::string advancedRepr() const;
    MEDCOUPLING_EXPORT const DataArrayDouble *getDirectAccessOfCoordsArrIfInStructure() const;
    MEDCOUPLING_EXPORT DataArrayDouble *getCoords();
    MEDCOUPLING_EXPORT const DataArrayDouble *getCoords() const;
    MEDCOUPLING_EXPORT void setCoords(const DataArrayDouble *coords);
    MEDCOUPLING_EXPORT void setNodeGridStructure(const int *gridStructBg, const int *gridStructEnd);
    MEDCOUPLING_EXPORT std::vector<int> getNodeGridStructure() const;
    MEDCOUPLING_EXPORT MEDCouplingStructuredMesh *buildStructuredSubPart(const std::vector< std::pair<int,int> >& cellPart) const;
    // tools
    MEDCOUPLING_EXPORT void getBoundingBox(double *bbox) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureField(bool isAbs) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureFieldOnNode(bool isAbs) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildOrthogonalField() const;
    MEDCOUPLING_EXPORT int getCellContainingPoint(const double *pos, double eps) const;
    MEDCOUPLING_EXPORT void getCellsContainingPoint(const double *pos, double eps, std::vector<int>& elts) const;
    MEDCOUPLING_EXPORT void rotate(const double *center, const double *vector, double angle);
    MEDCOUPLING_EXPORT void translate(const double *vector);
    MEDCOUPLING_EXPORT void scale(const double *point, double factor);
    MEDCOUPLING_EXPORT MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const;
    MEDCOUPLING_EXPORT DataArrayDouble *getCoordinatesAndOwner() const;
    MEDCOUPLING_EXPORT DataArrayDouble *computeCellCenterOfMass() const;
    MEDCOUPLING_EXPORT DataArrayDouble *computeIsoBarycenterOfNodesPerCell() const;
    MEDCOUPLING_EXPORT void renumberCells(const int *old2NewBg, bool check=true);
    //some useful methods
    MEDCOUPLING_EXPORT void getNodeGridStructure(int *res) const;
    //serialisation-unserialization
    MEDCOUPLING_EXPORT void getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const;
    MEDCOUPLING_EXPORT void resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const;
    MEDCOUPLING_EXPORT void serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const;
    MEDCOUPLING_EXPORT void unserialization(const std::vector<double>& tinyInfoD, const std::vector<int>& tinyInfo, const DataArrayInt *a1, DataArrayDouble *a2,
                                            const std::vector<std::string>& littleStrings);
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
    MEDCOUPLING_EXPORT std::string getVTKFileExtension() const;
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
    ~MEDCouplingCurveLinearMesh();
    void writeVTKLL(std::ostream& ofs, const std::string& cellData, const std::string& pointData, DataArrayByte *byteData) const;
    std::string getVTKDataSetType() const;
  private:
    MCAuto<DataArrayDouble> _coords;
    std::vector<int> _structure;
  };
}

#endif
