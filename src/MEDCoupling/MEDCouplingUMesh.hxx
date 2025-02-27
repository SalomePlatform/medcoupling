// Copyright (C) 2007-2025  CEA, EDF
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

#pragma once

#include "MEDCoupling.hxx"
#include "MEDCouplingPointSet.hxx"
#include "MEDCouplingMemArray.hxx"

#include "CellModel.hxx"

#include <set>

namespace MEDCoupling
{
  class MEDCouplingUMeshCellByTypeEntry;
  class MEDCouplingUMeshCellIterator;
  class MEDCoupling1SGTUMesh;
  class MEDCoupling1GTUMesh;
  class MEDCouplingSkyLineArray;

  class MEDCouplingUMesh : public MEDCouplingPointSet
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *New();
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *New(const std::string& meshName, int meshDim);
    std::string getClassName() const override { return std::string("MEDCouplingUMesh"); }
    // Copy methods
    MEDCOUPLING_EXPORT MEDCouplingUMesh *deepCopy() const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *clone(bool recDeepCpy) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *deepCopyConnectivityOnly() const;

    MEDCOUPLING_EXPORT void shallowCopyConnectivityFrom(const MEDCouplingPointSet *other);
    MEDCOUPLING_EXPORT void updateTime() const;
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDCouplingMeshType getType() const { return UNSTRUCTURED; }
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const;
    MEDCOUPLING_EXPORT void checkFastEquivalWith(const MEDCouplingMesh *other, double prec) const;
    MEDCOUPLING_EXPORT void checkConsistencyLight() const;
    MEDCOUPLING_EXPORT void checkConsistency(double eps=1e-12) const;
    MEDCOUPLING_EXPORT void checkGeomConsistency(double eps=1e-12) const;
    MEDCOUPLING_EXPORT void setMeshDimension(int meshDim);
    MEDCOUPLING_EXPORT void allocateCells(mcIdType nbOfCells=0);
    MEDCOUPLING_EXPORT void insertNextCell(INTERP_KERNEL::NormalizedCellType type, mcIdType size, const mcIdType *nodalConnOfCell);
    MEDCOUPLING_EXPORT void finishInsertingCells();
    MEDCOUPLING_EXPORT MEDCouplingUMeshCellIterator *cellIterator();
    MEDCOUPLING_EXPORT MEDCouplingUMeshCellByTypeEntry *cellsByType();
    MEDCOUPLING_EXPORT std::set<INTERP_KERNEL::NormalizedCellType> getAllGeoTypes() const;
    MEDCOUPLING_EXPORT std::vector<INTERP_KERNEL::NormalizedCellType> getAllGeoTypesSorted() const;
    MEDCOUPLING_EXPORT std::set<INTERP_KERNEL::NormalizedCellType> getTypesOfPart(const mcIdType *begin, const mcIdType *end) const;
    MEDCOUPLING_EXPORT void setConnectivity(DataArrayIdType *conn, DataArrayIdType *connIndex, bool isComputingTypes=true);
    const DataArrayIdType *getNodalConnectivity() const { return _nodal_connec; }
    const DataArrayIdType *getNodalConnectivityIndex() const { return _nodal_connec_index; }
    DataArrayIdType *getNodalConnectivity() { return _nodal_connec; }
    DataArrayIdType *getNodalConnectivityIndex() { return _nodal_connec_index; }
    MEDCOUPLING_EXPORT INTERP_KERNEL::NormalizedCellType getTypeOfCell(mcIdType cellId) const;
    MEDCOUPLING_EXPORT DataArrayIdType *giveCellsWithType(INTERP_KERNEL::NormalizedCellType type) const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const;
    MEDCOUPLING_EXPORT void getNodeIdsOfCell(mcIdType cellId, std::vector<mcIdType>& conn) const;
    MEDCOUPLING_EXPORT std::string simpleRepr() const;
    MEDCOUPLING_EXPORT std::string advancedRepr() const;
    MEDCOUPLING_EXPORT std::string cppRepr() const;
    MEDCOUPLING_EXPORT std::string reprConnectivityOfThis() const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildSetInstanceFromThis(std::size_t spaceDim) const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfNodesInCell(mcIdType cellId) const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfCells() const;
    MEDCOUPLING_EXPORT int getMeshDimension() const;
    MEDCOUPLING_EXPORT mcIdType getNodalConnectivityArrayLen() const;
    MEDCOUPLING_EXPORT void computeTypes();
    //! size of returned tinyInfo must be always the same.
    MEDCOUPLING_EXPORT void getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<mcIdType>& tinyInfo, std::vector<std::string>& littleStrings) const;
    MEDCOUPLING_EXPORT bool isEmptyMesh(const std::vector<mcIdType>& tinyInfo) const;
    MEDCOUPLING_EXPORT void resizeForUnserialization(const std::vector<mcIdType>& tinyInfo, DataArrayIdType *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const;
    MEDCOUPLING_EXPORT void serialize(DataArrayIdType *&a1, DataArrayDouble *&a2) const;
    MEDCOUPLING_EXPORT void unserialization(const std::vector<double>& tinyInfoD, const std::vector<mcIdType>& tinyInfo, const DataArrayIdType *a1, DataArrayDouble *a2, const std::vector<std::string>& littleStrings);
    MEDCOUPLING_EXPORT std::string getVTKDataSetType() const;
    MEDCOUPLING_EXPORT std::string getVTKFileExtension() const;
    MEDCOUPLING_EXPORT void writeVTKLL(std::ostream& ofs, const std::string& cellData, const std::string& pointData, DataArrayByte *byteData) const;
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
    //tools
    MEDCOUPLING_EXPORT static int AreCellsEqual(const mcIdType *conn, const mcIdType *connI, mcIdType cell1, mcIdType cell2, int compType);
    MEDCOUPLING_EXPORT static int AreCellsEqualPolicy0(const mcIdType *conn, const mcIdType *connI, mcIdType cell1, mcIdType cell2);
    MEDCOUPLING_EXPORT static int AreCellsEqualPolicy1(const mcIdType *conn, const mcIdType *connI, mcIdType cell1, mcIdType cell2);
    MEDCOUPLING_EXPORT static int AreCellsEqualPolicy2(const mcIdType *conn, const mcIdType *connI, mcIdType cell1, mcIdType cell2);
    MEDCOUPLING_EXPORT static int AreCellsEqualPolicy2NoType(const mcIdType *conn, const mcIdType *connI, mcIdType cell1, mcIdType cell2);
    MEDCOUPLING_EXPORT static int AreCellsEqualPolicy7(const mcIdType *conn, const mcIdType *connI, mcIdType cell1, mcIdType cell2);
    MEDCOUPLING_EXPORT void convertToPolyTypes(const mcIdType *cellIdsToConvertBg, const mcIdType *cellIdsToConvertEnd);
    MEDCOUPLING_EXPORT void convertAllToPoly();
    MEDCOUPLING_EXPORT void convertExtrudedPolyhedra();
    MEDCOUPLING_EXPORT bool unPolyze();
    MEDCOUPLING_EXPORT void simplifyPolyhedra(double eps);
    MEDCOUPLING_EXPORT void colinearizeEdges(double eps);
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildSpreadZonesWithPoly() const;
    MEDCOUPLING_EXPORT std::vector<DataArrayIdType *> partitionBySpreadZone() const;
    MEDCOUPLING_EXPORT DataArrayIdType *computeFetchedNodeIds() const;
    MEDCOUPLING_EXPORT DataArrayIdType *getNodeIdsInUse(mcIdType& nbrOfNodesInUse) const;
    MEDCOUPLING_EXPORT void computeNodeIdsAlg(std::vector<bool>& nodeIdsInUse) const;
    MEDCOUPLING_EXPORT DataArrayIdType *computeNbOfNodesPerCell() const;
    MEDCOUPLING_EXPORT DataArrayIdType *computeNbOfFacesPerCell() const;
    MEDCOUPLING_EXPORT DataArrayIdType *computeEffectiveNbOfNodesPerCell() const;
    MEDCOUPLING_EXPORT DataArrayIdType *zipCoordsTraducer();
    MEDCOUPLING_EXPORT void findCommonCells(int compType, mcIdType startCellId, DataArrayIdType *& commonCellsArr, DataArrayIdType *& commonCellsIArr) const;
    MEDCOUPLING_EXPORT bool areCellsIncludedIn(const MEDCouplingUMesh *other, int compType, DataArrayIdType *& arr) const;
    MEDCOUPLING_EXPORT bool areCellsIncludedInPolicy7(const MEDCouplingUMesh *other, DataArrayIdType *& arr) const;
    MEDCOUPLING_EXPORT void getReverseNodalConnectivity(DataArrayIdType *revNodal, DataArrayIdType *revNodalIndx) const;
    MEDCOUPLING_EXPORT MCAuto<MEDCouplingUMesh> explodeIntoEdges(MCAuto<DataArrayIdType>& desc, MCAuto<DataArrayIdType>& descIndex, MCAuto<DataArrayIdType>& revDesc, MCAuto<DataArrayIdType>& revDescIndx) const;
    MEDCOUPLING_EXPORT MCAuto<MEDCouplingUMesh> explodeMeshTo(int targetDeltaLevel, MCAuto<DataArrayIdType>& desc, MCAuto<DataArrayIdType>& descIndx, MCAuto<DataArrayIdType>& revDesc, MCAuto<DataArrayIdType>& revDescIndx) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *explode3DMeshTo1D(DataArrayIdType *desc, DataArrayIdType *descIndx, DataArrayIdType *revDesc, DataArrayIdType *revDescIndx) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildDescendingConnectivity(DataArrayIdType *desc, DataArrayIdType *descIndx, DataArrayIdType *revDesc, DataArrayIdType *revDescIndx) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildDescendingConnectivity2(DataArrayIdType *desc, DataArrayIdType *descIndx, DataArrayIdType *revDesc, DataArrayIdType *revDescIndx) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *explodeMeshIntoMicroEdges(DataArrayIdType *desc, DataArrayIdType *descIndx, DataArrayIdType *revDesc, DataArrayIdType *revDescIndx) const;
    MEDCOUPLING_EXPORT void computeNeighborsOfCells(DataArrayIdType *&neighbors, DataArrayIdType *&neighborsIdx) const;
    MEDCOUPLING_EXPORT void computeCellNeighborhoodFromNodesOne(const DataArrayIdType *nodeNeigh, const DataArrayIdType *nodeNeighI, MCAuto<DataArrayIdType>& cellNeigh, MCAuto<DataArrayIdType>& cellNeighIndex) const;
    MEDCOUPLING_EXPORT static void ComputeNeighborsOfCellsAdv(const DataArrayIdType *desc, const DataArrayIdType *descI, const DataArrayIdType *revDesc, const DataArrayIdType *revDescI,
                                                              DataArrayIdType *&neighbors, DataArrayIdType *&neighborsIdx);
    MEDCOUPLING_EXPORT void computeNeighborsOfNodes(DataArrayIdType *&neighbors, DataArrayIdType *&neighborsIdx) const;
    MEDCOUPLING_EXPORT void computeEnlargedNeighborsOfNodes(MCAuto<DataArrayIdType> &neighbors, MCAuto<DataArrayIdType>& neighborsIdx) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *mergeMyselfWithOnSameCoords(const MEDCouplingPointSet *other) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildPartOfMySelf(const mcIdType *begin, const mcIdType *end, bool keepCoords=true) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildPartOfMySelfSlice(mcIdType start, mcIdType end, mcIdType step, bool keepCoords=true) const;
    MEDCOUPLING_EXPORT void setPartOfMySelf(const mcIdType *cellIdsBg, const mcIdType *cellIdsEnd, const MEDCouplingUMesh& otherOnSameCoordsThanThis);
    MEDCOUPLING_EXPORT void setPartOfMySelfSlice(mcIdType start, mcIdType end, mcIdType step, const MEDCouplingUMesh& otherOnSameCoordsThanThis);
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildFacePartOfMySelfNode(const mcIdType *begin, const mcIdType *end, bool fullyIn) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildUnstructured() const;
    MEDCOUPLING_EXPORT DataArrayIdType *findBoundaryNodes() const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildBoundaryMesh(bool keepCoords) const;
    MEDCOUPLING_EXPORT DataArrayIdType *findCellIdsOnBoundary() const;
    MEDCOUPLING_EXPORT void findCellIdsLyingOn(const MEDCouplingUMesh& otherDimM1OnSameCoords, DataArrayIdType *&cellIdsRk0, DataArrayIdType *&cellIdsRk1) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *computeSkin() const;
    MEDCOUPLING_EXPORT DataArrayIdType *findNodesToDuplicate(const MEDCouplingUMesh& otherDimM1OnSameCoords) const;
    MEDCOUPLING_EXPORT void findCellsToRenumber(const MEDCouplingUMesh& otherDimM1OnSameCoords, const mcIdType *nodeIdsToDuplicateBg, const mcIdType *nodeIdsToDuplicateEnd,
                                                 DataArrayIdType *& cellIdsNeededToBeRenum, DataArrayIdType *& cellIdsNotModified) const;
    MEDCOUPLING_EXPORT void duplicateNodes(const mcIdType *nodeIdsToDuplicateBg, const mcIdType *nodeIdsToDuplicateEnd);
    MEDCOUPLING_EXPORT void renumberNodesWithOffsetInConn(mcIdType offset);
    MEDCOUPLING_EXPORT void renumberNodesInConn(const INTERP_KERNEL::HashMap<mcIdType,mcIdType>& newNodeNumbersO2N);
    MEDCOUPLING_EXPORT void renumberNodesInConn(const mcIdType *newNodeNumbersO2N);
    MEDCOUPLING_EXPORT void renumberNodesInConn(const std::map<mcIdType,mcIdType>& newNodeNumbersO2N) override;
    MEDCOUPLING_EXPORT void shiftNodeNumbersInConn(mcIdType delta);
    MEDCOUPLING_EXPORT void duplicateNodesInConn(const mcIdType *nodeIdsToDuplicateBg, const mcIdType *nodeIdsToDuplicateEnd, mcIdType offset);
    MEDCOUPLING_EXPORT void renumberCells(const mcIdType *old2NewBg, bool check=true);
    MEDCOUPLING_EXPORT DataArrayIdType *getCellsInBoundingBox(const double *bbox, double eps) const;
    MEDCOUPLING_EXPORT DataArrayIdType *getCellsInBoundingBox(const INTERP_KERNEL::DirectedBoundingBox& bbox, double eps);
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureField(bool isAbs) const;
    MEDCOUPLING_EXPORT DataArrayDouble *getPartMeasureField(bool isAbs, const mcIdType *begin, const mcIdType *end) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureFieldOnNode(bool isAbs) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildOrthogonalField() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildPartOrthogonalField(const mcIdType *begin, const mcIdType *end) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildDirectionVectorField() const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildSlice3D(const double *origin, const double *vec, double eps, DataArrayIdType *&cellIds) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildSlice3DSurf(const double *origin, const double *vec, double eps, DataArrayIdType *&cellIds) const;
    MEDCOUPLING_EXPORT MCAuto<MEDCouplingUMesh> clipSingle3DCellByPlane(const double origin[3], const double vec[3], double eps) const;
    MEDCOUPLING_EXPORT DataArrayIdType *getCellIdsCrossingPlane(const double *origin, const double *vec, double eps) const;
    MEDCOUPLING_EXPORT bool isContiguous1D() const;
    MEDCOUPLING_EXPORT void project1D(const double *pt, const double *v, double eps, double *res) const;
    MEDCOUPLING_EXPORT double distanceToPoint(const double *ptBg, const double *ptEnd, mcIdType& cellId) const;
    MEDCOUPLING_EXPORT DataArrayDouble *distanceToPoints(const DataArrayDouble *pts, DataArrayIdType *& cellIds) const;
    MEDCOUPLING_EXPORT mcIdType getCellContainingPoint(const double *pos, double eps) const;
    MEDCOUPLING_EXPORT void getCellsContainingPoint(const double *pos, double eps, std::vector<mcIdType>& elts) const;
    MEDCOUPLING_EXPORT void getCellsContainingPoints(const double *pos, mcIdType nbOfPoints, double eps, MCAuto<DataArrayIdType>& elts, MCAuto<DataArrayIdType>& eltsIndex) const override;
    MEDCOUPLING_EXPORT void getCellsContainingPointsLinearPartOnlyOnNonDynType(const double *pos, mcIdType nbOfPoints, double eps, MCAuto<DataArrayIdType>& elts, MCAuto<DataArrayIdType>& eltsIndex) const override;
    MEDCOUPLING_EXPORT void checkButterflyCells(std::vector<mcIdType>& cells, double eps=1e-12) const;
    MEDCOUPLING_EXPORT DataArrayIdType *convexEnvelop2D();
    MEDCOUPLING_EXPORT DataArrayIdType *findAndCorrectBadOriented3DExtrudedCells();
    MEDCOUPLING_EXPORT DataArrayIdType *findAndCorrectBadOriented3DCells();
    MEDCOUPLING_EXPORT DataArrayDouble *getBoundingBoxForBBTree(double arcDetEps=1e-12) const;
    MEDCOUPLING_EXPORT DataArrayDouble *getBoundingBoxForBBTreeFast() const;
    MEDCOUPLING_EXPORT DataArrayDouble *getBoundingBoxForBBTree2DQuadratic(double arcDetEps=1e-12) const;
    MEDCOUPLING_EXPORT DataArrayDouble *getBoundingBoxForBBTree1DQuadratic(double arcDetEps=1e-12) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildExtrudedMesh(const MEDCouplingUMesh *mesh1D, int policy);
    MEDCOUPLING_EXPORT bool isFullyQuadratic() const;
    MEDCOUPLING_EXPORT bool isPresenceOfQuadratic() const;
    MEDCOUPLING_EXPORT void convertQuadraticCellsToLinear();
    MEDCOUPLING_EXPORT DataArrayIdType *convertLinearCellsToQuadratic(int conversionType=0);
    MEDCOUPLING_EXPORT MCAuto<MEDCouplingUMesh> convertToQuadraticBasedOnSeg3(const MEDCoupling1SGTUMesh *seg3) const;
    MEDCOUPLING_EXPORT void tessellate2D(double eps);
    MEDCOUPLING_EXPORT MEDCoupling1SGTUMesh *tetrahedrize(int policy, DataArrayIdType *& n2oCells, mcIdType& nbOfAdditionalPoints) const;
    MEDCOUPLING_EXPORT DataArrayIdType *simplexize(int policy);
    MEDCOUPLING_EXPORT bool areOnlySimplexCells() const;
    MEDCOUPLING_EXPORT void convertDegeneratedCells();
    MEDCOUPLING_EXPORT DataArrayIdType *convertDegeneratedCellsAndRemoveFlatOnes();
    MEDCOUPLING_EXPORT bool removeDegenerated1DCells();
    MEDCOUPLING_EXPORT void are2DCellsNotCorrectlyOriented(const double *vec, bool polyOnly, std::vector<mcIdType>& cells) const;
    MEDCOUPLING_EXPORT void orientCorrectly2DCells(const double *vec, bool polyOnly);
    MEDCOUPLING_EXPORT void orientCorrectly2DCells(const MEDCouplingUMesh* refFaces = nullptr);
    MEDCOUPLING_EXPORT void changeOrientationOfCells();
    MEDCOUPLING_EXPORT void arePolyhedronsNotCorrectlyOriented(std::vector<mcIdType>& cells) const;
    MEDCOUPLING_EXPORT void orientCorrectlyPolyhedrons();
    MEDCOUPLING_EXPORT void invertOrientationOfAllCells();
    MEDCOUPLING_EXPORT void getFastAveragePlaneOfThis(double *vec, double *pos) const;
    MEDCOUPLING_EXPORT void attractSeg3MidPtsAroundNodes(double ratio, const mcIdType *nodeIdsBg, const mcIdType *nodeIdsEnd);
    //Mesh quality
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getEdgeRatioField() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getAspectRatioField() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getWarpField() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getSkewField() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *computeDiameterField() const;
    //utilities for MED File RW
    MEDCOUPLING_EXPORT std::vector<mcIdType> getDistributionOfTypes() const;
    MEDCOUPLING_EXPORT DataArrayIdType *checkTypeConsistencyAndContig(const std::vector<mcIdType>& code, const std::vector<const DataArrayIdType *>& idsPerType) const;
    MEDCOUPLING_EXPORT void splitProfilePerType(const DataArrayIdType *profile, std::vector<mcIdType>& code, std::vector<DataArrayIdType *>& idsInPflPerType, std::vector<DataArrayIdType *>& idsPerType, bool smartPflKiller=true) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *emulateMEDMEMBDC(const MEDCouplingUMesh *nM1LevMesh, DataArrayIdType *desc, DataArrayIdType *descIndx, DataArrayIdType *&revDesc, DataArrayIdType *&revDescIndx, DataArrayIdType *& nM1LevMeshIds, DataArrayIdType *&meshnM1Old2New) const;
    MEDCOUPLING_EXPORT DataArrayIdType *sortCellsInMEDFileFrmt();
    MEDCOUPLING_EXPORT bool checkConsecutiveCellTypes() const;
    MEDCOUPLING_EXPORT bool checkConsecutiveCellTypesForMEDFileFrmt() const;
    MEDCOUPLING_EXPORT bool checkConsecutiveCellTypesAndOrder(const INTERP_KERNEL::NormalizedCellType *orderBg, const INTERP_KERNEL::NormalizedCellType *orderEnd) const;
    MEDCOUPLING_EXPORT DataArrayIdType *getLevArrPerCellTypes(const INTERP_KERNEL::NormalizedCellType *orderBg, const INTERP_KERNEL::NormalizedCellType *orderEnd, DataArrayIdType *&nbPerType) const;
    MEDCOUPLING_EXPORT DataArrayIdType *getRenumArrForMEDFileFrmt() const;
    MEDCOUPLING_EXPORT DataArrayIdType *getRenumArrForConsecutiveCellTypesSpec(const INTERP_KERNEL::NormalizedCellType *orderBg, const INTERP_KERNEL::NormalizedCellType *orderEnd) const;
    MEDCOUPLING_EXPORT DataArrayIdType *rearrange2ConsecutiveCellTypes();
    MEDCOUPLING_EXPORT std::vector<MEDCouplingUMesh *> splitByType() const;
    MEDCOUPLING_EXPORT MEDCoupling1GTUMesh *convertIntoSingleGeoTypeMesh() const;
    MEDCOUPLING_EXPORT DataArrayIdType *convertNodalConnectivityToStaticGeoTypeMesh() const;
    MEDCOUPLING_EXPORT void convertNodalConnectivityToDynamicGeoTypeMesh(DataArrayIdType *&nodalConn, DataArrayIdType *&nodalConnIndex) const;
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *AggregateSortedByTypeMeshesOnSameCoords(const std::vector<const MEDCouplingUMesh *>& ms,
                                                                                        DataArrayIdType *&szOfCellGrpOfSameType,
                                                                                        DataArrayIdType *&idInMsOfCellGrpOfSameType);
    MEDCOUPLING_EXPORT DataArrayIdType *keepCellIdsByType(INTERP_KERNEL::NormalizedCellType type, const mcIdType *begin, const mcIdType *end) const;
    MEDCOUPLING_EXPORT DataArrayIdType *convertCellArrayPerGeoType(const DataArrayIdType *da) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *keepSpecifiedCells(INTERP_KERNEL::NormalizedCellType type, const mcIdType *idsPerGeoTypeBg, const mcIdType *idsPerGeoTypeEnd) const;
    MEDCOUPLING_EXPORT std::vector<bool> getQuadraticStatus() const;
    //
    MEDCOUPLING_EXPORT MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const;
    MEDCOUPLING_EXPORT DataArrayDouble *computeCellCenterOfMass() const;
    MEDCOUPLING_EXPORT DataArrayDouble *computeCellCenterOfMassWithPrecision(double eps) const;
    MEDCOUPLING_EXPORT DataArrayDouble *computeIsoBarycenterOfNodesPerCell() const;
    MEDCOUPLING_EXPORT DataArrayDouble *getPartBarycenterAndOwner(const mcIdType *begin, const mcIdType *end) const;
    MEDCOUPLING_EXPORT DataArrayDouble *computePlaneEquationOf3DFaces() const;
    MEDCOUPLING_EXPORT DataArrayIdType *conformize2D(double eps);
    MEDCOUPLING_EXPORT DataArrayIdType *colinearize2D(double eps);
    MEDCOUPLING_EXPORT DataArrayIdType *colinearizeKeepingConform2D(double eps);
    MEDCOUPLING_EXPORT DataArrayIdType *conformize3D(double eps);
    MEDCOUPLING_EXPORT mcIdType split2DCells(const DataArrayIdType *desc, const DataArrayIdType *descI, const DataArrayIdType *subNodesInSeg, const DataArrayIdType *subNodesInSegI, const DataArrayIdType *midOpt=0, const DataArrayIdType *midOptI=0);
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *Build0DMeshFromCoords(DataArrayDouble *da);
    MEDCOUPLING_EXPORT static MCAuto<MEDCouplingUMesh> Build1DMeshFromCoords(DataArrayDouble *da);
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *MergeUMeshes(const MEDCouplingUMesh *mesh1, const MEDCouplingUMesh *mesh2);
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *MergeUMeshes(const std::vector<const MEDCouplingUMesh *>& a);
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *MergeUMeshesOnSameCoords(const MEDCouplingUMesh *mesh1, const MEDCouplingUMesh *mesh2);
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *MergeUMeshesOnSameCoords(const std::vector<const MEDCouplingUMesh *>& meshes);
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *FuseUMeshesOnSameCoords(const std::vector<const MEDCouplingUMesh *>& meshes, int compType, std::vector<DataArrayIdType *>& corr);
    MEDCOUPLING_EXPORT static void PutUMeshesOnSameAggregatedCoords(const std::vector<MEDCouplingUMesh *>& meshes);
    MEDCOUPLING_EXPORT static void MergeNodesOnUMeshesSharingSameCoords(const std::vector<MEDCouplingUMesh *>& meshes, double eps);
    MEDCOUPLING_EXPORT static bool IsPolygonWellOriented(bool isQuadratic, const double *vec, const mcIdType *begin, const mcIdType *end, const double *coords);
    MEDCOUPLING_EXPORT static bool IsPolyhedronWellOriented(const mcIdType *begin, const mcIdType *end, const double *coords);
    MEDCOUPLING_EXPORT static bool Is3DExtrudedStaticCellWellOriented(const mcIdType *begin, const mcIdType *end, const double *coords);
    MEDCOUPLING_EXPORT static void CorrectExtrudedStaticCell(mcIdType *begin, mcIdType *end);
    MEDCOUPLING_EXPORT static bool IsTetra4WellOriented(const mcIdType *begin, const mcIdType *end, const double *coords);
    MEDCOUPLING_EXPORT static bool IsPyra5WellOriented(const mcIdType *begin, const mcIdType *end, const double *coords);
    MEDCOUPLING_EXPORT static void SimplifyPolyhedronCell(double eps, const DataArrayDouble *coords, mcIdType index, DataArrayIdType *res, MEDCouplingUMesh *faces,
                                                          DataArrayIdType *E_Fi, DataArrayIdType *E_F, DataArrayIdType *F_Ei, DataArrayIdType *F_E);
    MEDCOUPLING_EXPORT static void ComputeVecAndPtOfFace(double eps, const double *coords, const mcIdType *begin, const mcIdType *end, double *v, double *p);
    MEDCOUPLING_EXPORT static void TryToCorrectPolyhedronOrientation(mcIdType *begin, mcIdType *end, const double *coords);
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *Intersect2DMeshes(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2, double eps, DataArrayIdType *&cellNb1, DataArrayIdType *&cellNb2);
    MEDCOUPLING_EXPORT static void Intersect2DMeshWith1DLine(const MEDCouplingUMesh *mesh2D, const MEDCouplingUMesh *mesh1D,
                                                             double eps, MEDCouplingUMesh *&splitMesh2D, MEDCouplingUMesh *&splitMesh1D, DataArrayIdType *&cellIdInMesh2D, DataArrayIdType *&cellIdInMesh1D);
    MEDCOUPLING_EXPORT static bool BuildConvexEnvelopOf2DCellJarvis(const double *coords, const mcIdType *nodalConnBg, const mcIdType *nodalConnEnd, DataArrayIdType *nodalConnecOut);
    MEDCOUPLING_EXPORT static std::vector<DataArrayIdType *> PartitionBySpreadZone(const DataArrayIdType *arrIn, const DataArrayIdType *arrIndxIn);
    MEDCOUPLING_EXPORT static DataArrayIdType *ComputeSpreadZoneGradually(const DataArrayIdType *arrIn, const DataArrayIdType *arrIndxIn);
    MEDCOUPLING_EXPORT static DataArrayIdType *ComputeSpreadZoneGraduallyFromSeed(const mcIdType *seedBg, const mcIdType *seedEnd, const DataArrayIdType *arrIn, const DataArrayIdType *arrIndxIn, mcIdType nbOfDepthPeeling, mcIdType& nbOfDepthPeelingPerformed);
    MEDCOUPLING_EXPORT static void FindCommonCellsAlg(int compType, mcIdType startCellId, const DataArrayIdType *nodal, const DataArrayIdType *nodalI, const DataArrayIdType *revNodal, const DataArrayIdType *revNodalI,
                                                      DataArrayIdType *& commonCellsArr, DataArrayIdType *& commonCellsIArr);
    MEDCOUPLING_EXPORT DataArrayIdType *buildUnionOf2DMesh() const;
    MEDCOUPLING_EXPORT DataArrayIdType *buildUnionOf3DMesh() const;
    MEDCOUPLING_EXPORT DataArrayIdType *orderConsecutiveCells1D() const;
    MEDCOUPLING_EXPORT MEDCouplingSkyLineArray* generateGraph() const;
  private: // all private methods are impl in MEDCouplingUMesh_internal.cxx

    MEDCouplingUMesh();
    MEDCouplingUMesh(const MEDCouplingUMesh& other, bool deepCpy);
    ~MEDCouplingUMesh();
    void checkFullyDefined() const;
    void checkConnectivityFullyDefined() const;
    void reprConnectivityOfThisLL(std::ostringstream& stream) const;
    //tools
    DataArrayIdType *simplexizePol0();
    DataArrayIdType *simplexizePol1();
    DataArrayIdType *simplexizePlanarFace5();
    DataArrayIdType *simplexizePlanarFace6();
    void tessellate2DInternal(double eps);
    void tessellate2DCurveInternal(double eps);
    void subDivide2DMesh(const mcIdType *nodeSubdived, const mcIdType *nodeIndxSubdived, const mcIdType *desc, const mcIdType *descIndex);
    void fillCellIdsToKeepFromNodeIds(const mcIdType *begin, const mcIdType *end, bool fullyIn, DataArrayIdType *&cellIdsKeptArr) const;
    void split3DCurveWithPlane(const double *origin, const double *vec, double eps, std::vector<mcIdType>& cut3DCurve);
    MEDCouplingUMesh *buildExtrudedMeshFromThisLowLev(mcIdType nbOfNodesOf1Lev, bool isQuad) const;
    DataArrayDouble *fillExtCoordsUsingTranslation(const MEDCouplingUMesh *mesh1D, bool isQuad) const;
    DataArrayDouble *fillExtCoordsUsingTranslAndAutoRotation(const MEDCouplingUMesh *mesh1D, bool isQuad) const;
    DataArrayDouble *fillExtCoordsUsingTranslAndAutoRotation2D(const MEDCouplingUMesh *mesh1D, bool isQuad) const;
    DataArrayDouble *fillExtCoordsUsingTranslAndAutoRotation3D(const MEDCouplingUMesh *mesh1D, bool isQuad) const;
    static bool AreCellsEqualInPool(const std::vector<mcIdType>& candidates, int compType, const mcIdType *conn, const mcIdType *connI, DataArrayIdType *result) ;
    MEDCouplingUMesh *buildPartOfMySelfKeepCoords(const mcIdType *begin, const mcIdType *end) const;
    MEDCouplingUMesh *buildPartOfMySelfKeepCoordsSlice(mcIdType start, mcIdType end, mcIdType step) const;
    DataArrayIdType *convertLinearCellsToQuadratic1D0(DataArrayIdType *&conn, DataArrayIdType *&connI, DataArrayDouble *& coords, std::set<INTERP_KERNEL::NormalizedCellType>& types) const;
    DataArrayIdType *convertLinearCellsToQuadratic2DAnd3D0(const MEDCouplingUMesh *m1D, const DataArrayIdType *desc, const DataArrayIdType *descI, DataArrayIdType *&conn, DataArrayIdType *&connI, DataArrayDouble *& coords, std::set<INTERP_KERNEL::NormalizedCellType>& types) const;
    DataArrayIdType *convertLinearCellsToQuadratic2D0(DataArrayIdType *&conn, DataArrayIdType *&connI, DataArrayDouble *& coords, std::set<INTERP_KERNEL::NormalizedCellType>& types) const;
    DataArrayIdType *convertLinearCellsToQuadratic2D1(DataArrayIdType *&conn, DataArrayIdType *&connI, DataArrayDouble *& coords, std::set<INTERP_KERNEL::NormalizedCellType>& types) const;
    DataArrayIdType *convertLinearCellsToQuadratic3D0(DataArrayIdType *&conn, DataArrayIdType *&connI, DataArrayDouble *& coords, std::set<INTERP_KERNEL::NormalizedCellType>& types) const;
    DataArrayIdType *convertLinearCellsToQuadratic3D1(DataArrayIdType *&conn, DataArrayIdType *&connI, DataArrayDouble *& coords, std::set<INTERP_KERNEL::NormalizedCellType>& types) const;
    DataArrayIdType *buildUnionOf2DMeshLinear(const MEDCouplingUMesh *skin, const DataArrayIdType *n2o) const;
    DataArrayIdType *buildUnionOf2DMeshQuadratic(const MEDCouplingUMesh *skin, const DataArrayIdType *n2o) const;
    template<int SPACEDIM>
    void getCellsContainingPointsAlg(const double *coords, const double *pos, mcIdType nbOfPoints,
                                     double eps, MCAuto<DataArrayIdType>& elts, MCAuto<DataArrayIdType>& eltsIndex, std::function<bool(INTERP_KERNEL::NormalizedCellType,int)> sensibilityTo2DQuadraticLinearCellsFunc) const;
    void getCellsContainingPointsZeAlg(const double *pos, mcIdType nbOfPoints, double eps,
                                       MCAuto<DataArrayIdType>& elts, MCAuto<DataArrayIdType>& eltsIndex,
                                       std::function<bool(INTERP_KERNEL::NormalizedCellType,mcIdType)> sensibilityTo2DQuadraticLinearCellsFunc) const;
/// @cond INTERNAL
    static void DeleteCellTypeInIndexedArray(const DataArrayIdType *arrIn, const DataArrayIdType *arrIndxIn, MCAuto<DataArrayIdType>& arrOut, MCAuto<DataArrayIdType>& arrIndxOut);
    static MEDCouplingUMesh *MergeUMeshesLL(const std::vector<const MEDCouplingUMesh *>& a);
    typedef mcIdType (*DimM1DescNbrer)(mcIdType id, mcIdType nb, const INTERP_KERNEL::CellModel& cm, bool compute, const mcIdType *conn1, const mcIdType *conn2);
    template<class SonsGenerator>
    MEDCouplingUMesh *buildDescendingConnectivityGen(DataArrayIdType *desc, DataArrayIdType *descIndx, DataArrayIdType *revDesc, DataArrayIdType *revDescIndx, DimM1DescNbrer nbrer) const;
    static void DistanceToPoint3DSurfAlg(const double *pt, const mcIdType *cellIdsBg, const mcIdType *cellIdsEnd, const double *coords, const mcIdType *nc, const mcIdType *ncI, double& ret0, mcIdType& cellId);
    static void DistanceToPoint2DCurveAlg(const double *pt, const mcIdType *cellIdsBg, const mcIdType *cellIdsEnd, const double *coords, const mcIdType *nc, const mcIdType *ncI, double& ret0, mcIdType& cellId);
    static DataArrayIdType *ComputeSpreadZoneGraduallyFromSeedAlg(std::vector<bool>& fetched, const mcIdType *seedBg, const mcIdType *seedEnd, const DataArrayIdType *arrIn, const DataArrayIdType *arrIndxIn, mcIdType nbOfDepthPeeling, mcIdType& nbOfDepthPeelingPerformed);
    static void FillInCompact3DMode(int spaceDim, mcIdType nbOfNodesInCell, const mcIdType *conn, const double *coo, double *zipFrmt);
    static void AppendExtrudedCell(const mcIdType *connBg, const mcIdType *connEnd, mcIdType nbOfNodesPerLev, bool isQuad, std::vector<mcIdType>& ret);
    static void Intersect1DMeshes(const MEDCouplingUMesh *m1Desc, const MEDCouplingUMesh *m2Desc, double eps, std::vector< std::vector<mcIdType> >& intersectEdge1, std::vector< std::vector<mcIdType> >& colinear2, std::vector< std::vector<mcIdType> >& subDiv2, std::vector<double>& addCoo, std::map<mcIdType,mcIdType>& mergedNodes);
    static void IntersectDescending2DMeshes(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2, double eps,
                                            std::vector< std::vector<mcIdType> >& intersectEdge1, std::vector< std::vector<mcIdType> >& colinear2, std::vector< std::vector<mcIdType> >& subDiv2,
                                            MEDCouplingUMesh *& m1Desc, DataArrayIdType *&desc1, DataArrayIdType *&descIndx1, DataArrayIdType *&revDesc1, DataArrayIdType *&revDescIndx1,
                                            std::vector<double>& addCoo,
                                            MEDCouplingUMesh *& m2Desc, DataArrayIdType *&desc2, DataArrayIdType *&descIndx2, DataArrayIdType *&revDesc2, DataArrayIdType *&revDescIndx2);
    static void BuildIntersectEdges(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2, const std::vector<double>& addCoo, const std::vector< std::vector<mcIdType> >& subDiv, std::vector< std::vector<mcIdType> >& intersectEdge);
    static void BuildIntersecting2DCellsFromEdges(double eps, const MEDCouplingUMesh *m1, const mcIdType *desc1, const mcIdType *descIndx1, const std::vector<std::vector<mcIdType> >& intesctEdges1, const std::vector< std::vector<mcIdType> >& colinear2,
                                                  const MEDCouplingUMesh *m2, const mcIdType *desc2, const mcIdType *descIndx2, const std::vector<std::vector<mcIdType> >& intesctEdges2,
                                                  const std::vector<double>& addCoords,
                                                  std::vector<double>& addCoordsQuadratic, std::vector<mcIdType>& cr, std::vector<mcIdType>& crI, std::vector<mcIdType>& cNb1, std::vector<mcIdType>& cNb2);
    static void AssemblyForSplitFrom3DCurve(const std::vector<mcIdType>& cut3DCurve, std::vector<mcIdType>& nodesOnPlane, const mcIdType *nodal3DSurf, const mcIdType *nodalIndx3DSurf,
                                              const mcIdType *nodal3DCurve, const mcIdType *nodalIndx3DCurve,
                                              const mcIdType *desc, const mcIdType *descIndx, std::vector< std::pair<mcIdType,mcIdType> >& cut3DSurf);
    void assemblyForSplitFrom3DSurf(const std::vector< std::pair<mcIdType,mcIdType> >& cut3DSurf,
                                    const mcIdType *desc, const mcIdType *descIndx, DataArrayIdType *nodalRes, DataArrayIdType *nodalResIndx, DataArrayIdType *cellIds) const;
    void buildSubCellsFromCut(const std::vector< std::pair<mcIdType,mcIdType> >& cut3DSurf, const mcIdType *desc, const mcIdType *descIndx, const double *coords, double eps, std::vector<std::vector<mcIdType> >& res) const;
    void split2DCellsLinear(const DataArrayIdType *desc, const DataArrayIdType *descI, const DataArrayIdType *subNodesInSeg, const DataArrayIdType *subNodesInSegI);
    mcIdType split2DCellsQuadratic(const DataArrayIdType *desc, const DataArrayIdType *descI, const DataArrayIdType *subNodesInSeg, const DataArrayIdType *subNodesInSegI, const DataArrayIdType *mid, const DataArrayIdType *midI);
    static bool Colinearize2DCell(const double *coords, const mcIdType *connBg, const mcIdType *connEnd, mcIdType offset, const std::map<mcIdType, bool>& forbiddenPoints, DataArrayIdType *newConnOfCell, DataArrayDouble *appendedCoords);
    static void ComputeAllTypesInternal(std::set<INTERP_KERNEL::NormalizedCellType>& types, const DataArrayIdType *nodalConnec, const DataArrayIdType *nodalConnecIndex);
    static bool OrderPointsAlongLine(const double * coo, mcIdType startNode, mcIdType endNode,
                                                const mcIdType * c, const mcIdType * cI, const mcIdType *idsBg, const mcIdType *endBg,
                                                std::vector<mcIdType> & pointIds, std::vector<mcIdType> & hitSegs);
    static void ReplaceEdgeInFace(const mcIdType * sIdxConn, const mcIdType * sIdxConnE, mcIdType startNode, mcIdType endNode,
                                      const std::vector<mcIdType>& insidePoints, std::vector<mcIdType>& modifiedFace);
    void attractSeg3MidPtsAroundNodesUnderground(double ratio, const mcIdType *nodeIdsBg, const mcIdType *nodeIdsEnd);
    DataArrayIdType *internalColinearize2D(double eps, bool stayConform);
    template<class MAPCLS>
    void renumberNodesInConnT(const MAPCLS& newNodeNumbersO2N);
  public:
    MEDCOUPLING_EXPORT static DataArrayIdType *ComputeRangesFromTypeDistribution(const std::vector<mcIdType>& code);
    MEDCOUPLING_EXPORT static const int N_MEDMEM_ORDER=25;
    MEDCOUPLING_EXPORT static const INTERP_KERNEL::NormalizedCellType MEDMEM_ORDER[N_MEDMEM_ORDER];
    /// @endcond
  private:
    int _mesh_dim;
    DataArrayIdType *_nodal_connec;
    DataArrayIdType *_nodal_connec_index;
    std::set<INTERP_KERNEL::NormalizedCellType> _types;
  public:
    static double EPS_FOR_POLYH_ORIENTATION;
  };

  class MEDCouplingUMeshCell;

  class MEDCouplingUMeshCellIterator
  {
  public:
    MEDCOUPLING_EXPORT MEDCouplingUMeshCellIterator(MEDCouplingUMesh *mesh);
    MEDCOUPLING_EXPORT MEDCouplingUMeshCellIterator(MEDCouplingUMesh *mesh, MEDCouplingUMeshCell *itc, mcIdType bg, mcIdType end);
    MEDCOUPLING_EXPORT ~MEDCouplingUMeshCellIterator();
    MEDCOUPLING_EXPORT MEDCouplingUMeshCell *nextt();
  private:
    MEDCouplingUMesh *_mesh;
    MEDCouplingUMeshCell *_cell;
    bool _own_cell;
    mcIdType _cell_id;
    mcIdType _nb_cell;
  };

  class MEDCouplingUMeshCellByTypeIterator;

  class MEDCouplingUMeshCellByTypeEntry
  {
  public:
    MEDCOUPLING_EXPORT MEDCouplingUMeshCellByTypeEntry(MEDCouplingUMesh *mesh);
    MEDCOUPLING_EXPORT MEDCouplingUMeshCellByTypeIterator *iterator();
    MEDCOUPLING_EXPORT ~MEDCouplingUMeshCellByTypeEntry();
  private:
    MEDCouplingUMesh *_mesh;
  };

  class MEDCouplingUMeshCellEntry
  {
  public:
    MEDCOUPLING_EXPORT MEDCouplingUMeshCellEntry(MEDCouplingUMesh *mesh,  INTERP_KERNEL::NormalizedCellType type, MEDCouplingUMeshCell *itc, mcIdType bg, mcIdType end);
    MEDCOUPLING_EXPORT ~MEDCouplingUMeshCellEntry();
    MEDCOUPLING_EXPORT INTERP_KERNEL::NormalizedCellType getType() const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfElems() const;
    MEDCOUPLING_EXPORT MEDCouplingUMeshCellIterator *iterator();
  private:
    MEDCouplingUMesh *_mesh;
    INTERP_KERNEL::NormalizedCellType _type;
    MEDCouplingUMeshCell *_itc;
    mcIdType _bg;
    mcIdType _end;
  };

  class MEDCouplingUMeshCellByTypeIterator
  {
  public:
    MEDCOUPLING_EXPORT MEDCouplingUMeshCellByTypeIterator(MEDCouplingUMesh *mesh);
    MEDCOUPLING_EXPORT ~MEDCouplingUMeshCellByTypeIterator();
    MEDCOUPLING_EXPORT MEDCouplingUMeshCellEntry *nextt();
  private:
    MEDCouplingUMesh *_mesh;
    MEDCouplingUMeshCell *_cell;
    mcIdType _cell_id;
    mcIdType _nb_cell;
  };

  class MEDCouplingUMeshCell
  {
  public:
    MEDCOUPLING_EXPORT MEDCouplingUMeshCell(MEDCouplingUMesh *mesh);
    MEDCOUPLING_EXPORT void next();
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT INTERP_KERNEL::NormalizedCellType getType() const;
    MEDCOUPLING_EXPORT const mcIdType *getAllConn(mcIdType& lgth) const;
  private:
    mcIdType *_conn;
    mcIdType *_conn_indx;
    mcIdType _conn_lgth;
    static const int NOTICABLE_FIRST_VAL=-7;
  };

  template<typename T, typename TOUT>
  class UMeshGenIterator
  {
  public:
    using iterator_category = std::input_iterator_tag;
    using value_type = const TOUT *;
    using difference_type = mcIdType;
    using pointer = const TOUT **;
    using reference = const TOUT *;

    explicit UMeshGenIterator(std::size_t num , std::vector<const MEDCouplingUMesh *> *data) : _num(num),_data(data) {}
    UMeshGenIterator<T,TOUT>& operator++() { ++_num; return *this; }
    bool operator==(const UMeshGenIterator& other) const { return _num == other._num; }
    bool operator!=(const UMeshGenIterator& other) const { return !(*this == other); }
    reference operator*() const { T tt; return tt((*_data)[_num]); }

  private:
    std::size_t _num = 0;
    std::vector<const MEDCouplingUMesh *> *_data = nullptr;
  };

  struct UMeshIndexConnectivityFunctor { const DataArrayIdType *operator()(const MEDCouplingUMesh *um) { return um->getNodalConnectivityIndex(); } };

  using UMeshConnectivityIndexIterator = UMeshGenIterator<UMeshIndexConnectivityFunctor,DataArrayIdType>;

  struct UMeshConnectivityFunctor { const DataArrayIdType *operator()(const MEDCouplingUMesh *um) { return um->getNodalConnectivity(); } };

  using UMeshConnectivityIterator = UMeshGenIterator<UMeshConnectivityFunctor,DataArrayIdType>;

  struct UMeshCoordsFunctor { const DataArrayDouble *operator()(const MEDCouplingUMesh *um) { return um->getCoords(); } };

  using UMeshCoordsIterator = UMeshGenIterator<UMeshCoordsFunctor,DataArrayDouble>;
}

