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

#ifndef __PARAMEDMEM_MEDCOUPLINGIMESH_HXX__
#define __PARAMEDMEM_MEDCOUPLINGIMESH_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingStructuredMesh.hxx"

namespace MEDCoupling
{
  class MEDCouplingCMesh;

  class MEDCouplingIMesh : public MEDCouplingStructuredMesh
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingIMesh *New();
    MEDCOUPLING_EXPORT static MEDCouplingIMesh *New(const std::string& meshName, int spaceDim, const mcIdType *nodeStrctStart, const mcIdType *nodeStrctStop,
                                                    const double *originStart, const double *originStop, const double *dxyzStart, const double *dxyzStop);
    MEDCOUPLING_EXPORT std::string getClassName() const override { return std::string("MEDCouplingIMesh"); }
    //
    MEDCOUPLING_EXPORT void setSpaceDimension(int spaceDim);
    MEDCOUPLING_EXPORT void setNodeStruct(const mcIdType *nodeStrctStart, const mcIdType *nodeStrctStop);
    MEDCOUPLING_EXPORT std::vector<mcIdType> getNodeStruct() const;
    MEDCOUPLING_EXPORT void setOrigin(const double *originStart, const double *originStop);
    MEDCOUPLING_EXPORT std::vector<double> getOrigin() const;
    MEDCOUPLING_EXPORT void setDXYZ(const double *dxyzStart, const double *dxyzStop);
    MEDCOUPLING_EXPORT std::vector<double> getDXYZ() const;
    MEDCOUPLING_EXPORT void setAxisUnit(const std::string& unitName);
    MEDCOUPLING_EXPORT std::string getAxisUnit() const;
    MEDCOUPLING_EXPORT double getMeasureOfAnyCell() const;
    MEDCOUPLING_EXPORT MEDCouplingCMesh *convertToCartesian() const;
    MEDCOUPLING_EXPORT void refineWithFactor(const std::vector<mcIdType>& factors);
    MEDCOUPLING_EXPORT MEDCouplingIMesh *asSingleCell() const;
    MEDCOUPLING_EXPORT static void CondenseFineToCoarse(const std::vector<mcIdType>& coarseSt, const DataArrayDouble *fineDA, const std::vector< std::pair<mcIdType,mcIdType> >& fineLocInCoarse, const std::vector<mcIdType>& facts, DataArrayDouble *coarseDA);
    MEDCOUPLING_EXPORT static void CondenseFineToCoarseGhost(const std::vector<mcIdType>& coarseSt, const DataArrayDouble *fineDA, const std::vector< std::pair<mcIdType,mcIdType> >& fineLocInCoarse, const std::vector<mcIdType>& facts, DataArrayDouble *coarseDA, mcIdType ghostSize);
    MEDCOUPLING_EXPORT static void SpreadCoarseToFine(const DataArrayDouble *coarseDA, const std::vector<mcIdType>& coarseSt, DataArrayDouble *fineDA, const std::vector< std::pair<mcIdType,mcIdType> >& fineLocInCoarse, const std::vector<mcIdType>& facts);
    MEDCOUPLING_EXPORT static void SpreadCoarseToFineGhost(const DataArrayDouble *coarseDA, const std::vector<mcIdType>& coarseSt, DataArrayDouble *fineDA, const std::vector< std::pair<mcIdType,mcIdType> >& fineLocInCoarse, const std::vector<mcIdType>& facts, mcIdType ghostSize);
    MEDCOUPLING_EXPORT static void SpreadCoarseToFineGhostZone(const DataArrayDouble *coarseDA, const std::vector<mcIdType>& coarseSt, DataArrayDouble *fineDA, const std::vector< std::pair<mcIdType,mcIdType> >& fineLocInCoarse, const std::vector<mcIdType>& facts, mcIdType ghostSize);
    //
    MEDCOUPLING_EXPORT MEDCouplingIMesh *deepCopy() const;
    MEDCOUPLING_EXPORT MEDCouplingIMesh *clone(bool recDeepCpy) const;
    MEDCOUPLING_EXPORT const DataArrayDouble *getDirectAccessOfCoordsArrIfInStructure() const;
    MEDCOUPLING_EXPORT MEDCouplingIMesh *buildWithGhost(mcIdType ghostLev) const;
    MEDCOUPLING_EXPORT void updateTime() const;
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDCOUPLING_EXPORT MEDCouplingMeshType getType() const { return IMAGE_GRID; }
    MEDCOUPLING_EXPORT void copyTinyStringsFrom(const MEDCouplingMesh *other);
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const;
    MEDCOUPLING_EXPORT void checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                 DataArrayIdType *&cellCor, DataArrayIdType *&nodeCor) const;
    MEDCOUPLING_EXPORT void checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                            DataArrayIdType *&cellCor) const;
    MEDCOUPLING_EXPORT void checkConsistencyLight() const;
    MEDCOUPLING_EXPORT void checkConsistency(double eps=1e-12) const;
    MEDCOUPLING_EXPORT int getSpaceDimension() const;
    MEDCOUPLING_EXPORT void getCoordinatesOfNode(mcIdType nodeId, std::vector<double>& coo) const;
    MEDCOUPLING_EXPORT std::string simpleRepr() const;
    MEDCOUPLING_EXPORT std::string advancedRepr() const;
    // tools
    MEDCOUPLING_EXPORT void getBoundingBox(double *bbox) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureField(bool isAbs) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureFieldOnNode(bool isAbs) const;
    MEDCOUPLING_EXPORT mcIdType getCellContainingPoint(const double *pos, double eps) const;
    MEDCOUPLING_EXPORT void getCellsContainingPoint(const double *pos, double eps, std::vector<mcIdType>& elts) const;
    MEDCOUPLING_EXPORT void rotate(const double *center, const double *vector, double angle);
    MEDCOUPLING_EXPORT void translate(const double *vector);
    MEDCOUPLING_EXPORT void scale(const double *point, double factor);
    MEDCOUPLING_EXPORT MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const;
    MEDCOUPLING_EXPORT DataArrayDouble *getCoordinatesAndOwner() const;
    MEDCOUPLING_EXPORT DataArrayDouble *computeCellCenterOfMass() const;
    MEDCOUPLING_EXPORT DataArrayDouble *computeIsoBarycenterOfNodesPerCell() const;
    MEDCOUPLING_EXPORT void renumberCells(const mcIdType *old2NewBg, bool check=true);
    //some useful methods
    MEDCOUPLING_EXPORT void getNodeGridStructure(mcIdType *res) const;
    MEDCOUPLING_EXPORT std::vector<mcIdType> getNodeGridStructure() const;
    MEDCouplingStructuredMesh *buildStructuredSubPart(const std::vector< std::pair<mcIdType,mcIdType> >& cellPart) const;
    //serialisation-unserialization
    MEDCOUPLING_EXPORT void getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<mcIdType>& tinyInfo, std::vector<std::string>& littleStrings) const;
    MEDCOUPLING_EXPORT void resizeForUnserialization(const std::vector<mcIdType>& tinyInfo, DataArrayIdType *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const;
    MEDCOUPLING_EXPORT void serialize(DataArrayIdType *&a1, DataArrayDouble *&a2) const;
    MEDCOUPLING_EXPORT void unserialization(const std::vector<double>& tinyInfoD, const std::vector<mcIdType>& tinyInfo, const DataArrayIdType *a1, DataArrayDouble *a2,
                                            const std::vector<std::string>& littleStrings);
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
    MEDCOUPLING_EXPORT std::string getVTKFileExtension() const;
  private:
    MEDCouplingIMesh();
    MEDCouplingIMesh(const MEDCouplingIMesh& other, bool deepCopy);
    ~MEDCouplingIMesh();
    void writeVTKLL(std::ostream& ofs, const std::string& cellData, const std::string& pointData, DataArrayByte *byteData) const;
    std::string getVTKDataSetType() const;
    bool isEqualWithoutConsideringStrInternal(const MEDCouplingMesh *other, double prec, std::string& reason) const;
    std::vector<std::string> buildInfoOnComponents() const;
    void checkSpaceDimension() const;
    static void CheckSpaceDimension(int spaceDim);
    static int FindIntRoot(int val, int order);
    static void SpreadCoarseToFineGhost2D(const double *inPtr, double *outPtr, std::size_t nbCompo, const std::vector<mcIdType>& coarseSt, const std::vector< std::pair<mcIdType,mcIdType> >& fineLocInCoarse, const std::vector<mcIdType>& facts, mcIdType ghostSize);
    static void SpreadCoarseToFineGhostZone2D(const double *inPtr, double *outPtr, std::size_t nbCompo, const std::vector<mcIdType>& coarseSt, const std::vector< std::pair<mcIdType,mcIdType> >& fineLocInCoarse, const std::vector<mcIdType>& facts, mcIdType ghostSize);
  private:
    int _space_dim;
    double _origin[3];
    double _dxyz[3];
    mcIdType _structure[3];
    std::string _axis_unit;
  };
}

#endif
