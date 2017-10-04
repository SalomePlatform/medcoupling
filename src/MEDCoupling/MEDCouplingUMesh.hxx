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

#ifndef __PARAMEDMEM_MEDCOUPLINGUMESH_HXX__
#define __PARAMEDMEM_MEDCOUPLINGUMESH_HXX__

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
    // Copy methods
    MEDCOUPLING_EXPORT MEDCouplingUMesh *deepCopy() const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *clone(bool recDeepCpy) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *deepCopyConnectivityOnly() const;

    MEDCOUPLING_EXPORT void shallowCopyConnectivityFrom(const MEDCouplingPointSet *other);
    MEDCOUPLING_EXPORT void updateTime() const;
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDCOUPLING_EXPORT MEDCouplingMeshType getType() const { return UNSTRUCTURED; }
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const;
    MEDCOUPLING_EXPORT void checkFastEquivalWith(const MEDCouplingMesh *other, double prec) const;
    MEDCOUPLING_EXPORT void checkConsistencyLight() const;
    MEDCOUPLING_EXPORT void checkConsistency(double eps=1e-12) const;
    MEDCOUPLING_EXPORT void setMeshDimension(int meshDim);
    MEDCOUPLING_EXPORT void allocateCells(int nbOfCells=0);
    MEDCOUPLING_EXPORT void insertNextCell(INTERP_KERNEL::NormalizedCellType type, int size, const int *nodalConnOfCell);
    MEDCOUPLING_EXPORT void finishInsertingCells();
    MEDCOUPLING_EXPORT MEDCouplingUMeshCellIterator *cellIterator();
    MEDCOUPLING_EXPORT MEDCouplingUMeshCellByTypeEntry *cellsByType();
    MEDCOUPLING_EXPORT std::set<INTERP_KERNEL::NormalizedCellType> getAllGeoTypes() const;
    MEDCOUPLING_EXPORT std::vector<INTERP_KERNEL::NormalizedCellType> getAllGeoTypesSorted() const;
    MEDCOUPLING_EXPORT std::set<INTERP_KERNEL::NormalizedCellType> getTypesOfPart(const int *begin, const int *end) const;
    MEDCOUPLING_EXPORT void setConnectivity(DataArrayInt *conn, DataArrayInt *connIndex, bool isComputingTypes=true);
    MEDCOUPLING_EXPORT const DataArrayInt *getNodalConnectivity() const { return _nodal_connec; }
    MEDCOUPLING_EXPORT const DataArrayInt *getNodalConnectivityIndex() const { return _nodal_connec_index; }
    MEDCOUPLING_EXPORT DataArrayInt *getNodalConnectivity() { return _nodal_connec; }
    MEDCOUPLING_EXPORT DataArrayInt *getNodalConnectivityIndex() { return _nodal_connec_index; }
    MEDCOUPLING_EXPORT INTERP_KERNEL::NormalizedCellType getTypeOfCell(std::size_t cellId) const;
    MEDCOUPLING_EXPORT DataArrayInt *giveCellsWithType(INTERP_KERNEL::NormalizedCellType type) const;
    MEDCOUPLING_EXPORT std::size_t getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const;
    MEDCOUPLING_EXPORT void getNodeIdsOfCell(std::size_t cellId, std::vector<int>& conn) const;
    MEDCOUPLING_EXPORT std::string simpleRepr() const;
    MEDCOUPLING_EXPORT std::string advancedRepr() const;
    MEDCOUPLING_EXPORT std::string cppRepr() const;
    MEDCOUPLING_EXPORT std::string reprConnectivityOfThis() const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildSetInstanceFromThis(int spaceDim) const;
    MEDCOUPLING_EXPORT int getNumberOfNodesInCell(int cellId) const;
    MEDCOUPLING_EXPORT std::size_t getNumberOfCells() const;
    MEDCOUPLING_EXPORT int getMeshDimension() const;
    MEDCOUPLING_EXPORT int getNodalConnectivityArrayLen() const;
    MEDCOUPLING_EXPORT void computeTypes();
    //! size of returned tinyInfo must be always the same.
    MEDCOUPLING_EXPORT void getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const;
    MEDCOUPLING_EXPORT bool isEmptyMesh(const std::vector<int>& tinyInfo) const;
    MEDCOUPLING_EXPORT void resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const;
    MEDCOUPLING_EXPORT void serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const;
    MEDCOUPLING_EXPORT void unserialization(const std::vector<double>& tinyInfoD, const std::vector<int>& tinyInfo, const DataArrayInt *a1, DataArrayDouble *a2, const std::vector<std::string>& littleStrings);
    MEDCOUPLING_EXPORT std::string getVTKDataSetType() const;
    MEDCOUPLING_EXPORT std::string getVTKFileExtension() const;
    MEDCOUPLING_EXPORT void writeVTKLL(std::ostream& ofs, const std::string& cellData, const std::string& pointData, DataArrayByte *byteData) const;
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
    //tools
    MEDCOUPLING_EXPORT static int AreCellsEqual(const int *conn, const int *connI, int cell1, int cell2, int compType);
    MEDCOUPLING_EXPORT static int AreCellsEqualPolicy0(const int *conn, const int *connI, int cell1, int cell2);
    MEDCOUPLING_EXPORT static int AreCellsEqualPolicy1(const int *conn, const int *connI, int cell1, int cell2);
    MEDCOUPLING_EXPORT static int AreCellsEqualPolicy2(const int *conn, const int *connI, int cell1, int cell2);
    MEDCOUPLING_EXPORT static int AreCellsEqualPolicy2NoType(const int *conn, const int *connI, int cell1, int cell2);
    MEDCOUPLING_EXPORT static int AreCellsEqualPolicy7(const int *conn, const int *connI, int cell1, int cell2);
    MEDCOUPLING_EXPORT void convertToPolyTypes(const int *cellIdsToConvertBg, const int *cellIdsToConvertEnd);
    MEDCOUPLING_EXPORT void convertAllToPoly();
    MEDCOUPLING_EXPORT void convertExtrudedPolyhedra();
    MEDCOUPLING_EXPORT bool unPolyze();
    MEDCOUPLING_EXPORT void simplifyPolyhedra(double eps);
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildSpreadZonesWithPoly() const;
    MEDCOUPLING_EXPORT std::vector<DataArrayInt *> partitionBySpreadZone() const;
    MEDCOUPLING_EXPORT DataArrayInt *computeFetchedNodeIds() const;
    MEDCOUPLING_EXPORT DataArrayInt *getNodeIdsInUse(int& nbrOfNodesInUse) const;
    MEDCOUPLING_EXPORT void computeNodeIdsAlg(std::vector<bool>& nodeIdsInUse) const;
    MEDCOUPLING_EXPORT DataArrayInt *computeNbOfNodesPerCell() const;
    MEDCOUPLING_EXPORT DataArrayInt *computeNbOfFacesPerCell() const;
    MEDCOUPLING_EXPORT DataArrayInt *computeEffectiveNbOfNodesPerCell() const;
    MEDCOUPLING_EXPORT DataArrayInt *zipCoordsTraducer();
    MEDCOUPLING_EXPORT void findCommonCells(int compType, int startCellId, DataArrayInt *& commonCellsArr, DataArrayInt *& commonCellsIArr) const;
    MEDCOUPLING_EXPORT bool areCellsIncludedIn(const MEDCouplingUMesh *other, int compType, DataArrayInt *& arr) const;
    MEDCOUPLING_EXPORT bool areCellsIncludedInPolicy7(const MEDCouplingUMesh *other, DataArrayInt *& arr) const;
    MEDCOUPLING_EXPORT void getReverseNodalConnectivity(DataArrayInt *revNodal, DataArrayInt *revNodalIndx) const;
    MEDCOUPLING_EXPORT MCAuto<MEDCouplingUMesh> explodeIntoEdges(MCAuto<DataArrayInt>& desc, MCAuto<DataArrayInt>& descIndex, MCAuto<DataArrayInt>& revDesc, MCAuto<DataArrayInt>& revDescIndx) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *explode3DMeshTo1D(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildDescendingConnectivity(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildDescendingConnectivity2(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *explodeMeshIntoMicroEdges(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx) const;
    MEDCOUPLING_EXPORT void computeNeighborsOfCells(DataArrayInt *&neighbors, DataArrayInt *&neighborsIdx) const;
    MEDCOUPLING_EXPORT void computeCellNeighborhoodFromNodesOne(const DataArrayInt *nodeNeigh, const DataArrayInt *nodeNeighI, MCAuto<DataArrayInt>& cellNeigh, MCAuto<DataArrayInt>& cellNeighIndex) const;
    MEDCOUPLING_EXPORT static void ComputeNeighborsOfCellsAdv(const DataArrayInt *desc, const DataArrayInt *descI, const DataArrayInt *revDesc, const DataArrayInt *revDescI,
                                                              DataArrayInt *&neighbors, DataArrayInt *&neighborsIdx);
    MEDCOUPLING_EXPORT void computeNeighborsOfNodes(DataArrayInt *&neighbors, DataArrayInt *&neighborsIdx) const;
    MEDCOUPLING_EXPORT void computeEnlargedNeighborsOfNodes(MCAuto<DataArrayInt> &neighbors, MCAuto<DataArrayInt>& neighborsIdx) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *mergeMyselfWithOnSameCoords(const MEDCouplingPointSet *other) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildPartOfMySelf(const int *begin, const int *end, bool keepCoords=true) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildPartOfMySelfSlice(int start, int end, int step, bool keepCoords=true) const;
    MEDCOUPLING_EXPORT void setPartOfMySelf(const int *cellIdsBg, const int *cellIdsEnd, const MEDCouplingUMesh& otherOnSameCoordsThanThis);
    MEDCOUPLING_EXPORT void setPartOfMySelfSlice(int start, int end, int step, const MEDCouplingUMesh& otherOnSameCoordsThanThis);
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildFacePartOfMySelfNode(const int *begin, const int *end, bool fullyIn) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildUnstructured() const;
    MEDCOUPLING_EXPORT DataArrayInt *findBoundaryNodes() const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildBoundaryMesh(bool keepCoords) const;
    MEDCOUPLING_EXPORT DataArrayInt *findCellIdsOnBoundary() const;
    MEDCOUPLING_EXPORT void findCellIdsLyingOn(const MEDCouplingUMesh& otherDimM1OnSameCoords, DataArrayInt *&cellIdsRk0, DataArrayInt *&cellIdsRk1) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *computeSkin() const;
    MEDCOUPLING_EXPORT void findNodesToDuplicate(const MEDCouplingUMesh& otherDimM1OnSameCoords, DataArrayInt *& nodeIdsToDuplicate,
                                                 DataArrayInt *& cellIdsNeededToBeRenum, DataArrayInt *& cellIdsNotModified) const;
    MEDCOUPLING_EXPORT void duplicateNodes(const int *nodeIdsToDuplicateBg, const int *nodeIdsToDuplicateEnd);
    MEDCOUPLING_EXPORT void renumberNodesWithOffsetInConn(int offset);
    MEDCOUPLING_EXPORT void renumberNodesInConn(const INTERP_KERNEL::HashMap<int,int>& newNodeNumbersO2N);
    MEDCOUPLING_EXPORT void renumberNodesInConn(const int *newNodeNumbersO2N);
    MEDCOUPLING_EXPORT void shiftNodeNumbersInConn(int delta);
    MEDCOUPLING_EXPORT void duplicateNodesInConn(const int *nodeIdsToDuplicateBg, const int *nodeIdsToDuplicateEnd, int offset);
    MEDCOUPLING_EXPORT void renumberCells(const int *old2NewBg, bool check=true);
    MEDCOUPLING_EXPORT DataArrayInt *getCellsInBoundingBox(const double *bbox, double eps) const;
    MEDCOUPLING_EXPORT DataArrayInt *getCellsInBoundingBox(const INTERP_KERNEL::DirectedBoundingBox& bbox, double eps);
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureField(bool isAbs) const;
    MEDCOUPLING_EXPORT DataArrayDouble *getPartMeasureField(bool isAbs, const int *begin, const int *end) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureFieldOnNode(bool isAbs) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildOrthogonalField() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildPartOrthogonalField(const int *begin, const int *end) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildDirectionVectorField() const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildSlice3D(const double *origin, const double *vec, double eps, DataArrayInt *&cellIds) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildSlice3DSurf(const double *origin, const double *vec, double eps, DataArrayInt *&cellIds) const;
    MEDCOUPLING_EXPORT MCAuto<MEDCouplingUMesh> clipSingle3DCellByPlane(const double origin[3], const double vec[3], double eps) const;
    MEDCOUPLING_EXPORT DataArrayInt *getCellIdsCrossingPlane(const double *origin, const double *vec, double eps) const;
    MEDCOUPLING_EXPORT bool isContiguous1D() const;
    MEDCOUPLING_EXPORT void project1D(const double *pt, const double *v, double eps, double *res) const;
    MEDCOUPLING_EXPORT double distanceToPoint(const double *ptBg, const double *ptEnd, int& cellId) const;
    MEDCOUPLING_EXPORT DataArrayDouble *distanceToPoints(const DataArrayDouble *pts, DataArrayInt *& cellIds) const;
    MEDCOUPLING_EXPORT int getCellContainingPoint(const double *pos, double eps) const;
    MEDCOUPLING_EXPORT void getCellsContainingPoint(const double *pos, double eps, std::vector<int>& elts) const;
    MEDCOUPLING_EXPORT void getCellsContainingPoints(const double *pos, int nbOfPoints, double eps, MCAuto<DataArrayInt>& elts, MCAuto<DataArrayInt>& eltsIndex) const;
    MEDCOUPLING_EXPORT void checkButterflyCells(std::vector<int>& cells, double eps=1e-12) const;
    MEDCOUPLING_EXPORT DataArrayInt *convexEnvelop2D();
    MEDCOUPLING_EXPORT DataArrayInt *findAndCorrectBadOriented3DExtrudedCells();
    MEDCOUPLING_EXPORT DataArrayInt *findAndCorrectBadOriented3DCells();
    MEDCOUPLING_EXPORT DataArrayDouble *getBoundingBoxForBBTree(double arcDetEps=1e-12) const;
    MEDCOUPLING_EXPORT DataArrayDouble *getBoundingBoxForBBTreeFast() const;
    MEDCOUPLING_EXPORT DataArrayDouble *getBoundingBoxForBBTree2DQuadratic(double arcDetEps=1e-12) const;
    MEDCOUPLING_EXPORT DataArrayDouble *getBoundingBoxForBBTree1DQuadratic(double arcDetEps=1e-12) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildExtrudedMesh(const MEDCouplingUMesh *mesh1D, int policy);
    MEDCOUPLING_EXPORT bool isFullyQuadratic() const;
    MEDCOUPLING_EXPORT bool isPresenceOfQuadratic() const;
    MEDCOUPLING_EXPORT void convertQuadraticCellsToLinear();
    MEDCOUPLING_EXPORT DataArrayInt *convertLinearCellsToQuadratic(int conversionType=0);
    MEDCOUPLING_EXPORT void tessellate2D(double eps);
    MEDCOUPLING_EXPORT MEDCoupling1SGTUMesh *tetrahedrize(int policy, DataArrayInt *& n2oCells, int& nbOfAdditionalPoints) const;
    MEDCOUPLING_EXPORT DataArrayInt *simplexize(int policy);
    MEDCOUPLING_EXPORT bool areOnlySimplexCells() const;
    MEDCOUPLING_EXPORT void convertDegeneratedCells();
    MEDCOUPLING_EXPORT void are2DCellsNotCorrectlyOriented(const double *vec, bool polyOnly, std::vector<int>& cells) const;
    MEDCOUPLING_EXPORT void orientCorrectly2DCells(const double *vec, bool polyOnly);
    MEDCOUPLING_EXPORT void changeOrientationOfCells();
    MEDCOUPLING_EXPORT void arePolyhedronsNotCorrectlyOriented(std::vector<int>& cells) const;
    MEDCOUPLING_EXPORT void orientCorrectlyPolyhedrons();
    MEDCOUPLING_EXPORT void invertOrientationOfAllCells();
    MEDCOUPLING_EXPORT void getFastAveragePlaneOfThis(double *vec, double *pos) const;
    //Mesh quality
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getEdgeRatioField() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getAspectRatioField() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getWarpField() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getSkewField() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *computeDiameterField() const;
    //utilities for MED File RW
    MEDCOUPLING_EXPORT std::vector<int> getDistributionOfTypes() const;
    MEDCOUPLING_EXPORT DataArrayInt *checkTypeConsistencyAndContig(const std::vector<int>& code, const std::vector<const DataArrayInt *>& idsPerType) const;
    MEDCOUPLING_EXPORT void splitProfilePerType(const DataArrayInt *profile, std::vector<int>& code, std::vector<DataArrayInt *>& idsInPflPerType, std::vector<DataArrayInt *>& idsPerType) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *emulateMEDMEMBDC(const MEDCouplingUMesh *nM1LevMesh, DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *&revDesc, DataArrayInt *&revDescIndx, DataArrayInt *& nM1LevMeshIds, DataArrayInt *&meshnM1Old2New) const;
    MEDCOUPLING_EXPORT DataArrayInt *sortCellsInMEDFileFrmt();
    MEDCOUPLING_EXPORT bool checkConsecutiveCellTypes() const;
    MEDCOUPLING_EXPORT bool checkConsecutiveCellTypesForMEDFileFrmt() const;
    MEDCOUPLING_EXPORT bool checkConsecutiveCellTypesAndOrder(const INTERP_KERNEL::NormalizedCellType *orderBg, const INTERP_KERNEL::NormalizedCellType *orderEnd) const;
    MEDCOUPLING_EXPORT DataArrayInt *getLevArrPerCellTypes(const INTERP_KERNEL::NormalizedCellType *orderBg, const INTERP_KERNEL::NormalizedCellType *orderEnd, DataArrayInt *&nbPerType) const;
    MEDCOUPLING_EXPORT DataArrayInt *getRenumArrForMEDFileFrmt() const;
    MEDCOUPLING_EXPORT DataArrayInt *getRenumArrForConsecutiveCellTypesSpec(const INTERP_KERNEL::NormalizedCellType *orderBg, const INTERP_KERNEL::NormalizedCellType *orderEnd) const;
    MEDCOUPLING_EXPORT DataArrayInt *rearrange2ConsecutiveCellTypes();
    MEDCOUPLING_EXPORT std::vector<MEDCouplingUMesh *> splitByType() const;
    MEDCOUPLING_EXPORT MEDCoupling1GTUMesh *convertIntoSingleGeoTypeMesh() const;
    MEDCOUPLING_EXPORT DataArrayInt *convertNodalConnectivityToStaticGeoTypeMesh() const;
    MEDCOUPLING_EXPORT void convertNodalConnectivityToDynamicGeoTypeMesh(DataArrayInt *&nodalConn, DataArrayInt *&nodalConnIndex) const;
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *AggregateSortedByTypeMeshesOnSameCoords(const std::vector<const MEDCouplingUMesh *>& ms,
                                                                                        DataArrayInt *&szOfCellGrpOfSameType,
                                                                                        DataArrayInt *&idInMsOfCellGrpOfSameType);
    MEDCOUPLING_EXPORT DataArrayInt *keepCellIdsByType(INTERP_KERNEL::NormalizedCellType type, const int *begin, const int *end) const;
    MEDCOUPLING_EXPORT DataArrayInt *convertCellArrayPerGeoType(const DataArrayInt *da) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *keepSpecifiedCells(INTERP_KERNEL::NormalizedCellType type, const int *idsPerGeoTypeBg, const int *idsPerGeoTypeEnd) const;
    MEDCOUPLING_EXPORT std::vector<bool> getQuadraticStatus() const;
    //
    MEDCOUPLING_EXPORT MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const;
    MEDCOUPLING_EXPORT DataArrayDouble *computeCellCenterOfMass() const;
    MEDCOUPLING_EXPORT DataArrayDouble *computeIsoBarycenterOfNodesPerCell() const;
    MEDCOUPLING_EXPORT DataArrayDouble *getPartBarycenterAndOwner(const int *begin, const int *end) const;
    MEDCOUPLING_EXPORT DataArrayDouble *computePlaneEquationOf3DFaces() const;
    MEDCOUPLING_EXPORT DataArrayInt *conformize2D(double eps);
    MEDCOUPLING_EXPORT DataArrayInt *colinearize2D(double eps);
    MEDCOUPLING_EXPORT DataArrayInt *conformize3D(double eps);
    MEDCOUPLING_EXPORT int split2DCells(const DataArrayInt *desc, const DataArrayInt *descI, const DataArrayInt *subNodesInSeg, const DataArrayInt *subNodesInSegI, const DataArrayInt *midOpt=0, const DataArrayInt *midOptI=0);
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *Build0DMeshFromCoords(DataArrayDouble *da);
    MEDCOUPLING_EXPORT static MCAuto<MEDCouplingUMesh> Build1DMeshFromCoords(DataArrayDouble *da);
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *MergeUMeshes(const MEDCouplingUMesh *mesh1, const MEDCouplingUMesh *mesh2);
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *MergeUMeshes(const std::vector<const MEDCouplingUMesh *>& a);
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *MergeUMeshesOnSameCoords(const MEDCouplingUMesh *mesh1, const MEDCouplingUMesh *mesh2);
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *MergeUMeshesOnSameCoords(const std::vector<const MEDCouplingUMesh *>& meshes);
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *FuseUMeshesOnSameCoords(const std::vector<const MEDCouplingUMesh *>& meshes, int compType, std::vector<DataArrayInt *>& corr);
    MEDCOUPLING_EXPORT static void PutUMeshesOnSameAggregatedCoords(const std::vector<MEDCouplingUMesh *>& meshes);
    MEDCOUPLING_EXPORT static void MergeNodesOnUMeshesSharingSameCoords(const std::vector<MEDCouplingUMesh *>& meshes, double eps);
    MEDCOUPLING_EXPORT static bool IsPolygonWellOriented(bool isQuadratic, const double *vec, const int *begin, const int *end, const double *coords);
    MEDCOUPLING_EXPORT static bool IsPolyhedronWellOriented(const int *begin, const int *end, const double *coords);
    MEDCOUPLING_EXPORT static bool Is3DExtrudedStaticCellWellOriented(const int *begin, const int *end, const double *coords);
    MEDCOUPLING_EXPORT static void CorrectExtrudedStaticCell(int *begin, int *end);
    MEDCOUPLING_EXPORT static bool IsTetra4WellOriented(const int *begin, const int *end, const double *coords);
    MEDCOUPLING_EXPORT static bool IsPyra5WellOriented(const int *begin, const int *end, const double *coords);
    MEDCOUPLING_EXPORT static void SimplifyPolyhedronCell(double eps, const DataArrayDouble *coords, int index, DataArrayInt *res, MEDCouplingUMesh *faces,
                                                          DataArrayInt *E_Fi, DataArrayInt *E_F, DataArrayInt *F_Ei, DataArrayInt *F_E);
    MEDCOUPLING_EXPORT static void ComputeVecAndPtOfFace(double eps, const double *coords, const int *begin, const int *end, double *v, double *p);
    MEDCOUPLING_EXPORT static void TryToCorrectPolyhedronOrientation(int *begin, int *end, const double *coords);
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *Intersect2DMeshes(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2, double eps, DataArrayInt *&cellNb1, DataArrayInt *&cellNb2);
    MEDCOUPLING_EXPORT static void Intersect2DMeshWith1DLine(const MEDCouplingUMesh *mesh2D, const MEDCouplingUMesh *mesh1D,
                                                             double eps, MEDCouplingUMesh *&splitMesh2D, MEDCouplingUMesh *&splitMesh1D, DataArrayInt *&cellIdInMesh2D, DataArrayInt *&cellIdInMesh1D);
    MEDCOUPLING_EXPORT static bool BuildConvexEnvelopOf2DCellJarvis(const double *coords, const int *nodalConnBg, const int *nodalConnEnd, DataArrayInt *nodalConnecOut);
    MEDCOUPLING_EXPORT static bool RemoveIdsFromIndexedArrays(const int *idsToRemoveBg, const int *idsToRemoveEnd, DataArrayInt *arr, DataArrayInt *arrIndx, int offsetForRemoval=0);
    MEDCOUPLING_EXPORT static void ExtractFromIndexedArrays(const int *idsOfSelectBg, const int *idsOfSelectEnd, const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn,
                                                            DataArrayInt* &arrOut, DataArrayInt* &arrIndexOut);
    MEDCOUPLING_EXPORT static void ExtractFromIndexedArraysSlice(int idsOfSelectStart, int idsOfSelectStop, int idsOfSelectStep, const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn,
                                                             DataArrayInt* &arrOut, DataArrayInt* &arrIndexOut);
    MEDCOUPLING_EXPORT static void SetPartOfIndexedArrays(const int *idsOfSelectBg, const int *idsOfSelectEnd, const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn,
                                                          const DataArrayInt *srcArr, const DataArrayInt *srcArrIndex,
                                                          DataArrayInt* &arrOut, DataArrayInt* &arrIndexOut);
    MEDCOUPLING_EXPORT static void SetPartOfIndexedArraysSameIdx(const int *idsOfSelectBg, const int *idsOfSelectEnd, DataArrayInt *arrInOut, const DataArrayInt *arrIndxIn,
                                                                 const DataArrayInt *srcArr, const DataArrayInt *srcArrIndex);
    MEDCOUPLING_EXPORT static void SetPartOfIndexedArraysSlice(int start, int end, int step, const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn,
                                                           const DataArrayInt *srcArr, const DataArrayInt *srcArrIndex,
                                                           DataArrayInt* &arrOut, DataArrayInt* &arrIndexOut);
    MEDCOUPLING_EXPORT static void SetPartOfIndexedArraysSameIdxSlice(int start, int end, int step, DataArrayInt *arrInOut, const DataArrayInt *arrIndxIn,
                                                                  const DataArrayInt *srcArr, const DataArrayInt *srcArrIndex);
    MEDCOUPLING_EXPORT static std::vector<DataArrayInt *> PartitionBySpreadZone(const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn);
    MEDCOUPLING_EXPORT static DataArrayInt *ComputeSpreadZoneGradually(const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn);
    MEDCOUPLING_EXPORT static DataArrayInt *ComputeSpreadZoneGraduallyFromSeed(const int *seedBg, const int *seedEnd, const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn, int nbOfDepthPeeling, int& nbOfDepthPeelingPerformed);
    MEDCOUPLING_EXPORT static void FindCommonCellsAlg(int compType, int startCellId, const DataArrayInt *nodal, const DataArrayInt *nodalI, const DataArrayInt *revNodal, const DataArrayInt *revNodalI,
                                                      DataArrayInt *& commonCellsArr, DataArrayInt *& commonCellsIArr);
    MEDCOUPLING_EXPORT DataArrayInt *buildUnionOf2DMesh() const;
    MEDCOUPLING_EXPORT DataArrayInt *buildUnionOf3DMesh() const;
    MEDCOUPLING_EXPORT DataArrayInt *orderConsecutiveCells1D() const;
    MEDCOUPLING_EXPORT MEDCouplingSkyLineArray* generateGraph() const;
  private: // all private methods are impl in MEDCouplingUMesh_internal.cxx

    MEDCouplingUMesh();
    MEDCouplingUMesh(const MEDCouplingUMesh& other, bool deepCpy);
    ~MEDCouplingUMesh();
    void checkFullyDefined() const;
    void checkConnectivityFullyDefined() const;
    void reprConnectivityOfThisLL(std::ostringstream& stream) const;
    //tools
    DataArrayInt *simplexizePol0();
    DataArrayInt *simplexizePol1();
    DataArrayInt *simplexizePlanarFace5();
    DataArrayInt *simplexizePlanarFace6();
    void tessellate2DInternal(double eps);
    void tessellate2DCurveInternal(double eps);
    void subDivide2DMesh(const int *nodeSubdived, const int *nodeIndxSubdived, const int *desc, const int *descIndex);
    void fillCellIdsToKeepFromNodeIds(const int *begin, const int *end, bool fullyIn, DataArrayInt *&cellIdsKeptArr) const;
    void split3DCurveWithPlane(const double *origin, const double *vec, double eps, std::vector<int>& cut3DCurve);
    MEDCouplingUMesh *buildExtrudedMeshFromThisLowLev(int nbOfNodesOf1Lev, bool isQuad) const;
    DataArrayDouble *fillExtCoordsUsingTranslation(const MEDCouplingUMesh *mesh1D, bool isQuad) const;
    DataArrayDouble *fillExtCoordsUsingTranslAndAutoRotation(const MEDCouplingUMesh *mesh1D, bool isQuad) const;
    DataArrayDouble *fillExtCoordsUsingTranslAndAutoRotation2D(const MEDCouplingUMesh *mesh1D, bool isQuad) const;
    DataArrayDouble *fillExtCoordsUsingTranslAndAutoRotation3D(const MEDCouplingUMesh *mesh1D, bool isQuad) const;
    static bool AreCellsEqualInPool(const std::vector<int>& candidates, int compType, const int *conn, const int *connI, DataArrayInt *result) ;
    MEDCouplingUMesh *buildPartOfMySelfKeepCoords(const int *begin, const int *end) const;
    MEDCouplingUMesh *buildPartOfMySelfKeepCoordsSlice(int start, int end, int step) const;
    DataArrayInt *convertLinearCellsToQuadratic1D0(DataArrayInt *&conn, DataArrayInt *&connI, DataArrayDouble *& coords, std::set<INTERP_KERNEL::NormalizedCellType>& types) const;
    DataArrayInt *convertLinearCellsToQuadratic2DAnd3D0(const MEDCouplingUMesh *m1D, const DataArrayInt *desc, const DataArrayInt *descI, DataArrayInt *&conn, DataArrayInt *&connI, DataArrayDouble *& coords, std::set<INTERP_KERNEL::NormalizedCellType>& types) const;
    DataArrayInt *convertLinearCellsToQuadratic2D0(DataArrayInt *&conn, DataArrayInt *&connI, DataArrayDouble *& coords, std::set<INTERP_KERNEL::NormalizedCellType>& types) const;
    DataArrayInt *convertLinearCellsToQuadratic2D1(DataArrayInt *&conn, DataArrayInt *&connI, DataArrayDouble *& coords, std::set<INTERP_KERNEL::NormalizedCellType>& types) const;
    DataArrayInt *convertLinearCellsToQuadratic3D0(DataArrayInt *&conn, DataArrayInt *&connI, DataArrayDouble *& coords, std::set<INTERP_KERNEL::NormalizedCellType>& types) const;
    DataArrayInt *convertLinearCellsToQuadratic3D1(DataArrayInt *&conn, DataArrayInt *&connI, DataArrayDouble *& coords, std::set<INTERP_KERNEL::NormalizedCellType>& types) const;
    DataArrayInt *buildUnionOf2DMeshLinear(const MEDCouplingUMesh *skin, const DataArrayInt *n2o) const;
    DataArrayInt *buildUnionOf2DMeshQuadratic(const MEDCouplingUMesh *skin, const DataArrayInt *n2o) const;
    template<int SPACEDIM>
    void getCellsContainingPointsAlg(const double *coords, const double *pos, int nbOfPoints,
                                     double eps, MCAuto<DataArrayInt>& elts, MCAuto<DataArrayInt>& eltsIndex) const;
/// @cond INTERNAL
    static MEDCouplingUMesh *MergeUMeshesLL(const std::vector<const MEDCouplingUMesh *>& a);
    typedef int (*DimM1DescNbrer)(int id, unsigned nb, const INTERP_KERNEL::CellModel& cm, bool compute, const int *conn1, const int *conn2);
    template<class SonsGenerator>
    MEDCouplingUMesh *buildDescendingConnectivityGen(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx, DimM1DescNbrer nbrer) const;
    static void DistanceToPoint3DSurfAlg(const double *pt, const int *cellIdsBg, const int *cellIdsEnd, const double *coords, const int *nc, const int *ncI, double& ret0, int& cellId);
    static void DistanceToPoint2DCurveAlg(const double *pt, const int *cellIdsBg, const int *cellIdsEnd, const double *coords, const int *nc, const int *ncI, double& ret0, int& cellId);
    static DataArrayInt *ComputeSpreadZoneGraduallyFromSeedAlg(std::vector<bool>& fetched, const int *seedBg, const int *seedEnd, const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn, int nbOfDepthPeeling, int& nbOfDepthPeelingPerformed);
    static void FillInCompact3DMode(int spaceDim, int nbOfNodesInCell, const int *conn, const double *coo, double *zipFrmt);
    static void AppendExtrudedCell(const int *connBg, const int *connEnd, int nbOfNodesPerLev, bool isQuad, std::vector<int>& ret);
    static void Intersect1DMeshes(const MEDCouplingUMesh *m1Desc, const MEDCouplingUMesh *m2Desc, double eps, std::vector< std::vector<int> >& intersectEdge1, std::vector< std::vector<int> >& colinear2, std::vector< std::vector<int> >& subDiv2, std::vector<double>& addCoo, std::map<int,int>& mergedNodes);
    static void IntersectDescending2DMeshes(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2, double eps,
                                            std::vector< std::vector<int> >& intersectEdge1, std::vector< std::vector<int> >& colinear2, std::vector< std::vector<int> >& subDiv2,
                                            MEDCouplingUMesh *& m1Desc, DataArrayInt *&desc1, DataArrayInt *&descIndx1, DataArrayInt *&revDesc1, DataArrayInt *&revDescIndx1,
                                            std::vector<double>& addCoo,
                                            MEDCouplingUMesh *& m2Desc, DataArrayInt *&desc2, DataArrayInt *&descIndx2, DataArrayInt *&revDesc2, DataArrayInt *&revDescIndx2);
    static void BuildIntersectEdges(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2, const std::vector<double>& addCoo, const std::vector< std::vector<int> >& subDiv, std::vector< std::vector<int> >& intersectEdge);
    static void BuildIntersecting2DCellsFromEdges(double eps, const MEDCouplingUMesh *m1, const int *desc1, const int *descIndx1, const std::vector<std::vector<int> >& intesctEdges1, const std::vector< std::vector<int> >& colinear2,
                                                  const MEDCouplingUMesh *m2, const int *desc2, const int *descIndx2, const std::vector<std::vector<int> >& intesctEdges2,
                                                  const std::vector<double>& addCoords,
                                                  std::vector<double>& addCoordsQuadratic, std::vector<int>& cr, std::vector<int>& crI, std::vector<int>& cNb1, std::vector<int>& cNb2);
    static void AssemblyForSplitFrom3DCurve(const std::vector<int>& cut3DCurve, std::vector<int>& nodesOnPlane, const int *nodal3DSurf, const int *nodalIndx3DSurf,
                                              const int *nodal3DCurve, const int *nodalIndx3DCurve,
                                              const int *desc, const int *descIndx, std::vector< std::pair<int,int> >& cut3DSurf);
    void assemblyForSplitFrom3DSurf(const std::vector< std::pair<int,int> >& cut3DSurf,
                                    const int *desc, const int *descIndx, DataArrayInt *nodalRes, DataArrayInt *nodalResIndx, DataArrayInt *cellIds) const;
    void buildSubCellsFromCut(const std::vector< std::pair<int,int> >& cut3DSurf, const int *desc, const int *descIndx, const double *coords, double eps, std::vector<std::vector<int> >& res) const;
    void split2DCellsLinear(const DataArrayInt *desc, const DataArrayInt *descI, const DataArrayInt *subNodesInSeg, const DataArrayInt *subNodesInSegI);
    int split2DCellsQuadratic(const DataArrayInt *desc, const DataArrayInt *descI, const DataArrayInt *subNodesInSeg, const DataArrayInt *subNodesInSegI, const DataArrayInt *mid, const DataArrayInt *midI);
    static bool Colinearize2DCell(const double *coords, const int *connBg, const int *connEnd, int offset, DataArrayInt *newConnOfCell, DataArrayDouble *appendedCoords);
    static void ComputeAllTypesInternal(std::set<INTERP_KERNEL::NormalizedCellType>& types, const DataArrayInt *nodalConnec, const DataArrayInt *nodalConnecIndex);
    static bool OrderPointsAlongLine(const double * coo, int startNode, int endNode,
                                                const int * c, const int * cI, const int *idsBg, const int *endBg,
                                                std::vector<int> & pointIds, std::vector<int> & hitSegs);
    static void ReplaceEdgeInFace(const int * sIdxConn, const int * sIdxConnE, int startNode, int endNode,
                                      const std::vector<int>& insidePoints, std::vector<int>& modifiedFace);
  public:
    MEDCOUPLING_EXPORT static DataArrayInt *ComputeRangesFromTypeDistribution(const std::vector<int>& code);
    MEDCOUPLING_EXPORT static const int N_MEDMEM_ORDER=25;
    MEDCOUPLING_EXPORT static const INTERP_KERNEL::NormalizedCellType MEDMEM_ORDER[N_MEDMEM_ORDER];
    MEDCOUPLING_EXPORT static const int MEDCOUPLING2VTKTYPETRADUCER[INTERP_KERNEL::NORM_MAXTYPE+1];
    /// @endcond
  private:
    int _mesh_dim;
    DataArrayInt *_nodal_connec;
    DataArrayInt *_nodal_connec_index;
    std::set<INTERP_KERNEL::NormalizedCellType> _types;
  public:
    static double EPS_FOR_POLYH_ORIENTATION;
  };

  class MEDCouplingUMeshCell;

  class MEDCouplingUMeshCellIterator
  {
  public:
    MEDCOUPLING_EXPORT MEDCouplingUMeshCellIterator(MEDCouplingUMesh *mesh);
    MEDCOUPLING_EXPORT MEDCouplingUMeshCellIterator(MEDCouplingUMesh *mesh, MEDCouplingUMeshCell *itc, int bg, int end);
    MEDCOUPLING_EXPORT ~MEDCouplingUMeshCellIterator();
    MEDCOUPLING_EXPORT MEDCouplingUMeshCell *nextt();
  private:
    MEDCouplingUMesh *_mesh;
    MEDCouplingUMeshCell *_cell;
    bool _own_cell;
    int _cell_id;
    int _nb_cell;
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
    MEDCOUPLING_EXPORT MEDCouplingUMeshCellEntry(MEDCouplingUMesh *mesh,  INTERP_KERNEL::NormalizedCellType type, MEDCouplingUMeshCell *itc, int bg, int end);
    MEDCOUPLING_EXPORT ~MEDCouplingUMeshCellEntry();
    MEDCOUPLING_EXPORT INTERP_KERNEL::NormalizedCellType getType() const;
    MEDCOUPLING_EXPORT int getNumberOfElems() const;
    MEDCOUPLING_EXPORT MEDCouplingUMeshCellIterator *iterator();
  private:
    MEDCouplingUMesh *_mesh;
    INTERP_KERNEL::NormalizedCellType _type;
    MEDCouplingUMeshCell *_itc;
    int _bg;
    int _end;
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
    int _cell_id;
    int _nb_cell;
  };

  class MEDCouplingUMeshCell
  {
  public:
    MEDCOUPLING_EXPORT MEDCouplingUMeshCell(MEDCouplingUMesh *mesh);
    MEDCOUPLING_EXPORT void next();
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT INTERP_KERNEL::NormalizedCellType getType() const;
    MEDCOUPLING_EXPORT const int *getAllConn(int& lgth) const;
  private:
    int *_conn;
    int *_conn_indx;
    int _conn_lgth;
    static const int NOTICABLE_FIRST_VAL=-7;
  };
}

#endif
