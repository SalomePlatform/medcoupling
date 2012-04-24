// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

#ifndef __PARAMEDMEM_MEDCOUPLINGUMESH_HXX__
#define __PARAMEDMEM_MEDCOUPLINGUMESH_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingPointSet.hxx"
#include "MEDCouplingMemArray.hxx"

#include "CellModel.hxx"

#include <set>

namespace ParaMEDMEM
{
  class MEDCouplingUMeshCellIterator;
  class MEDCouplingUMeshCellByTypeEntry;

  class MEDCouplingUMesh : public MEDCouplingPointSet
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *New();
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *New(const char *meshName, int meshDim);
    MEDCOUPLING_EXPORT MEDCouplingMesh *deepCpy() const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *clone(bool recDeepCpy) const;
    MEDCOUPLING_EXPORT void updateTime() const;
    MEDCOUPLING_EXPORT MEDCouplingMeshType getType() const { return UNSTRUCTURED; }
    MEDCOUPLING_EXPORT bool isEqual(const MEDCouplingMesh *other, double prec) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const;
    MEDCOUPLING_EXPORT void checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                 DataArrayInt *&cellCor, DataArrayInt *&nodeCor) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                            DataArrayInt *&cellCor) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void checkFastEquivalWith(const MEDCouplingMesh *other, double prec) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void checkCoherency() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void checkCoherency1(double eps=1e-12) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void checkCoherency2(double eps=1e-12) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setMeshDimension(int meshDim);
    MEDCOUPLING_EXPORT void allocateCells(int nbOfCells);
    MEDCOUPLING_EXPORT void insertNextCell(INTERP_KERNEL::NormalizedCellType type, int size, const int *nodalConnOfCell) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void finishInsertingCells();
    MEDCOUPLING_EXPORT MEDCouplingUMeshCellIterator *cellIterator();
    MEDCOUPLING_EXPORT MEDCouplingUMeshCellByTypeEntry *cellsByType() throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT const std::set<INTERP_KERNEL::NormalizedCellType>& getAllTypes() const { return _types; }
    MEDCOUPLING_EXPORT std::set<INTERP_KERNEL::NormalizedCellType> getAllGeoTypes() const;
    MEDCOUPLING_EXPORT std::set<INTERP_KERNEL::NormalizedCellType> getTypesOfPart(const int *begin, const int *end) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setConnectivity(DataArrayInt *conn, DataArrayInt *connIndex, bool isComputingTypes=true);
    MEDCOUPLING_EXPORT const DataArrayInt *getNodalConnectivity() const { return _nodal_connec; }
    MEDCOUPLING_EXPORT const DataArrayInt *getNodalConnectivityIndex() const { return _nodal_connec_index; }
    MEDCOUPLING_EXPORT DataArrayInt *getNodalConnectivity() { return _nodal_connec; }
    MEDCOUPLING_EXPORT DataArrayInt *getNodalConnectivityIndex() { return _nodal_connec_index; }
    MEDCOUPLING_EXPORT INTERP_KERNEL::NormalizedCellType getTypeOfCell(int cellId) const;
    MEDCOUPLING_EXPORT int getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const;
    MEDCOUPLING_EXPORT void getNodeIdsOfCell(int cellId, std::vector<int>& conn) const;
    MEDCOUPLING_EXPORT DataArrayInt *getCellIdsFullyIncludedInNodeIds(const int *partBg, const int *partEnd) const;
    MEDCOUPLING_EXPORT std::string simpleRepr() const;
    MEDCOUPLING_EXPORT std::string advancedRepr() const;
    MEDCOUPLING_EXPORT std::string reprConnectivityOfThis() const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildSetInstanceFromThis(int spaceDim) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT int getNumberOfNodesInCell(int cellId) const;
    MEDCOUPLING_EXPORT int getNumberOfCells() const;
    MEDCOUPLING_EXPORT int getMeshDimension() const;
    MEDCOUPLING_EXPORT int getMeshLength() const;
    MEDCOUPLING_EXPORT void computeTypes();
    //! size of returned tinyInfo must be always the same.
    MEDCOUPLING_EXPORT void getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const;
    MEDCOUPLING_EXPORT bool isEmptyMesh(const std::vector<int>& tinyInfo) const;
    MEDCOUPLING_EXPORT void resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const;
    MEDCOUPLING_EXPORT void serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const;
    MEDCOUPLING_EXPORT void unserialization(const std::vector<double>& tinyInfoD, const std::vector<int>& tinyInfo, const DataArrayInt *a1, DataArrayDouble *a2, const std::vector<std::string>& littleStrings);
    MEDCOUPLING_EXPORT std::string getVTKDataSetType() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void writeVTKLL(std::ostream& ofs, const std::string& cellData, const std::string& pointData) const throw(INTERP_KERNEL::Exception);
    //tools
    MEDCOUPLING_EXPORT bool areCellsEqual(int cell1, int cell2, int compType) const;
    MEDCOUPLING_EXPORT bool areCellsEqual0(int cell1, int cell2) const;
    MEDCOUPLING_EXPORT bool areCellsEqual1(int cell1, int cell2) const;
    MEDCOUPLING_EXPORT bool areCellsEqual2(int cell1, int cell2) const;
    MEDCOUPLING_EXPORT bool areCellsFrom2MeshEqual(const MEDCouplingUMesh *other, int cellId, double prec) const;
    MEDCOUPLING_EXPORT void convertToPolyTypes(const int *cellIdsToConvertBg, const int *cellIdsToConvertEnd);
    MEDCOUPLING_EXPORT void convertAllToPoly();
    MEDCOUPLING_EXPORT void convertExtrudedPolyhedra() throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void unPolyze();
    MEDCOUPLING_EXPORT DataArrayInt *getNodeIdsInUse(int& nbrOfNodesInUse) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *zipCoordsTraducer() throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *zipConnectivityTraducer(int compType) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool areCellsIncludedIn(const MEDCouplingUMesh *other, int compType, DataArrayInt *& arr) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getReverseNodalConnectivity(DataArrayInt *revNodal, DataArrayInt *revNodalIndx) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildDescendingConnectivity(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildDescendingConnectivity2(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *mergeNodes(double precision, bool& areNodesMerged, int& newNbOfNodes);
    MEDCOUPLING_EXPORT DataArrayInt *mergeNodes2(double precision, bool& areNodesMerged, int& newNbOfNodes);
    MEDCOUPLING_EXPORT void tryToShareSameCoordsPermute(const MEDCouplingPointSet& other, double epsilon) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingPointSet *buildPartOfMySelf(const int *begin, const int *end, bool keepCoords) const;
    MEDCOUPLING_EXPORT MEDCouplingPointSet *buildPartOfMySelfNode(const int *begin, const int *end, bool fullyIn) const;
    MEDCOUPLING_EXPORT DataArrayInt *getCellIdsLyingOnNodes(const int *begin, const int *end, bool fullyIn) const;
    MEDCOUPLING_EXPORT MEDCouplingPointSet *buildFacePartOfMySelfNode(const int *begin, const int *end, bool fullyIn) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildUnstructured() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void findBoundaryNodes(std::vector<int>& nodes) const;
    MEDCOUPLING_EXPORT MEDCouplingPointSet *buildBoundaryMesh(bool keepCoords) const;
    MEDCOUPLING_EXPORT DataArrayInt *findCellsIdsOnBoundary() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void renumberNodes(const int *newNodeNumbers, int newNbOfNodes);
    MEDCOUPLING_EXPORT void renumberNodes2(const int *newNodeNumbers, int newNbOfNodes);
    MEDCOUPLING_EXPORT void renumberNodesInConn(const int *newNodeNumbersO2N);
    MEDCOUPLING_EXPORT void shiftNodeNumbersInConn(int delta) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void renumberCells(const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getCellsInBoundingBox(const double *bbox, double eps, std::vector<int>& elems) const;
    MEDCOUPLING_EXPORT void getCellsInBoundingBox(const INTERP_KERNEL::DirectedBoundingBox& bbox, double eps, std::vector<int>& elems);
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureField(bool isAbs) const;
    MEDCOUPLING_EXPORT DataArrayDouble *getPartMeasureField(bool isAbs, const int *begin, const int *end) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureFieldOnNode(bool isAbs) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildOrthogonalField() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildPartOrthogonalField(const int *begin, const int *end) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildDirectionVectorField() const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildSlice3D(const double *origin, const double *vec, double eps, DataArrayInt *&cellIds) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildSlice3DSurf(const double *origin, const double *vec, double eps, DataArrayInt *&cellIds) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *getCellIdsCrossingPlane(const double *origin, const double *vec, double eps) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool isContiguous1D() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void project1D(const double *pt, const double *v, double eps, double *res) const;
    MEDCOUPLING_EXPORT int getCellContainingPoint(const double *pos, double eps) const;
    MEDCOUPLING_EXPORT void getCellsContainingPoint(const double *pos, double eps, std::vector<int>& elts) const;
    MEDCOUPLING_EXPORT void getCellsContainingPoints(const double *pos, int nbOfPoints, double eps, std::vector<int>& elts, std::vector<int>& eltsIndex) const;
    MEDCOUPLING_EXPORT void checkButterflyCells(std::vector<int>& cells) const;
    MEDCOUPLING_EXPORT void findAndCorrectBadOriented3DExtrudedCells(std::vector<int>& cells) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getBoundingBoxForBBTree(std::vector<double>& bbox) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildExtrudedMesh(const MEDCouplingUMesh *mesh1D, int policy);
    MEDCOUPLING_EXPORT bool isFullyQuadratic() const;
    MEDCOUPLING_EXPORT bool isPresenceOfQuadratic() const;
    MEDCOUPLING_EXPORT void convertQuadraticCellsToLinear() throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void tessellate2D(double eps) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void tessellate2DCurve(double eps) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *simplexize(int policy) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool areOnlySimplexCells() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void convertDegeneratedCells() throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void are2DCellsNotCorrectlyOriented(const double *vec, bool polyOnly, std::vector<int>& cells) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void orientCorrectly2DCells(const double *vec, bool polyOnly) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void arePolyhedronsNotCorrectlyOriented(std::vector<int>& cells) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void orientCorrectlyPolyhedrons() throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getFastAveragePlaneOfThis(double *vec, double *pos) const throw(INTERP_KERNEL::Exception);
    //Mesh quality
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getEdgeRatioField() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getAspectRatioField() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getWarpField() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getSkewField() const throw(INTERP_KERNEL::Exception);
    //utilities for MED File RW
    MEDCOUPLING_EXPORT std::vector<int> getDistributionOfTypes() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *checkTypeConsistencyAndContig(const std::vector<int>& code, const std::vector<const DataArrayInt *>& idsPerType) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void splitProfilePerType(const DataArrayInt *profile, std::vector<int>& code, std::vector<DataArrayInt *>& idsInPflPerType, std::vector<DataArrayInt *>& idsPerType) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingUMesh *emulateMEDMEMBDC(const MEDCouplingUMesh *nM1LevMesh, DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *&revDesc, DataArrayInt *&revDescIndx, DataArrayInt *& nM1LevMeshIds, DataArrayInt *&meshnM1Old2New) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *sortCellsInMEDFileFrmt() throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool checkConsecutiveCellTypes() const;
    MEDCOUPLING_EXPORT bool checkConsecutiveCellTypesAndOrder(const INTERP_KERNEL::NormalizedCellType *orderBg, const INTERP_KERNEL::NormalizedCellType *orderEnd) const;
    MEDCOUPLING_EXPORT DataArrayInt *getLevArrPerCellTypes(const INTERP_KERNEL::NormalizedCellType *orderBg, const INTERP_KERNEL::NormalizedCellType *orderEnd, DataArrayInt *&nbPerType) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *getRenumArrForConsecutiveCellTypesSpec(const INTERP_KERNEL::NormalizedCellType *orderBg, const INTERP_KERNEL::NormalizedCellType *orderEnd) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *rearrange2ConsecutiveCellTypes();
    MEDCOUPLING_EXPORT std::vector<MEDCouplingUMesh *> splitByType() const;
    MEDCOUPLING_EXPORT DataArrayInt *keepCellIdsByType(INTERP_KERNEL::NormalizedCellType type, const int *begin, const int *end) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *convertCellArrayPerGeoType(const DataArrayInt *da) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingUMesh *keepSpecifiedCells(INTERP_KERNEL::NormalizedCellType type, const int *idsPerGeoTypeBg, const int *idsPerGeoTypeEnd) const;
    MEDCOUPLING_EXPORT std::vector<bool> getQuadraticStatus() const throw(INTERP_KERNEL::Exception);
    //
    MEDCOUPLING_EXPORT MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const;
    MEDCOUPLING_EXPORT DataArrayDouble *getBarycenterAndOwner() const;
    MEDCOUPLING_EXPORT DataArrayDouble *getPartBarycenterAndOwner(const int *begin, const int *end) const;
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *Build0DMeshFromCoords(DataArrayDouble *da) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *MergeUMeshes(const MEDCouplingUMesh *mesh1, const MEDCouplingUMesh *mesh2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *MergeUMeshes(std::vector<const MEDCouplingUMesh *>& a) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *MergeUMeshesOnSameCoords(const MEDCouplingUMesh *mesh1, const MEDCouplingUMesh *mesh2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *MergeUMeshesOnSameCoords(const std::vector<const MEDCouplingUMesh *>& meshes);
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *FuseUMeshesOnSameCoords(const std::vector<const MEDCouplingUMesh *>& meshes, int compType, std::vector<DataArrayInt *>& corr);
    MEDCOUPLING_EXPORT static void PutUMeshesOnSameAggregatedCoords(const std::vector<MEDCouplingUMesh *>& meshes) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static bool IsPolygonWellOriented(bool isQuadratic, const double *vec, const int *begin, const int *end, const double *coords);
    MEDCOUPLING_EXPORT static bool IsPolyhedronWellOriented(const int *begin, const int *end, const double *coords);
    MEDCOUPLING_EXPORT static void TryToCorrectPolyhedronOrientation(int *begin, int *end, const double *coords) throw(INTERP_KERNEL::Exception);
/// @cond INTERNAL
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *Intersect2DMeshes(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2, double eps, DataArrayInt *&cellNb1, DataArrayInt *&cellNb2) throw(INTERP_KERNEL::Exception);
/// @endcond
  private:
    MEDCouplingUMesh();
    MEDCouplingUMesh(const MEDCouplingUMesh& other, bool deepCopy);
    ~MEDCouplingUMesh();
    void checkFullyDefined() const throw(INTERP_KERNEL::Exception);
    void checkConnectivityFullyDefined() const throw(INTERP_KERNEL::Exception);
    void reprConnectivityOfThisLL(std::ostringstream& stream) const;
    //tools
    DataArrayInt *simplexizePol0() throw(INTERP_KERNEL::Exception);
    DataArrayInt *simplexizePol1() throw(INTERP_KERNEL::Exception);
    void subDivide2DMesh(const int *nodeSubdived, const int *nodeIndxSubdived, const int *desc, const int *descIndex) throw(INTERP_KERNEL::Exception);
    void fillCellIdsToKeepFromNodeIds(const int *begin, const int *end, bool fullyIn, std::vector<int>& cellIdsKept) const;
    void split3DCurveWithPlane(const double *origin, const double *vec, double eps, std::vector<int>& cut3DCurve) throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *buildExtrudedMeshFromThisLowLev(int nbOfNodesOf1Lev, bool isQuad) const;
    DataArrayDouble *fillExtCoordsUsingTranslation(const MEDCouplingUMesh *mesh1D, bool isQuad) const;
    DataArrayDouble *fillExtCoordsUsingTranslAndAutoRotation(const MEDCouplingUMesh *mesh1D, bool isQuad) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *fillExtCoordsUsingTranslAndAutoRotation2D(const MEDCouplingUMesh *mesh1D, bool isQuad) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *fillExtCoordsUsingTranslAndAutoRotation3D(const MEDCouplingUMesh *mesh1D, bool isQuad) const throw(INTERP_KERNEL::Exception);
    template<int SPACEDIM>
    void findCommonCellsBase(int compType, std::vector<int>& res, std::vector<int>& resI) const;
    bool areCellsEqualInPool(const std::vector<int>& candidates, int compType, std::vector<int>& result) const;
    MEDCouplingUMesh *buildPartOfMySelfKeepCoords(const int *begin, const int *end) const;
    template<int SPACEDIM>
    void getCellsContainingPointsAlg(const double *coords, const double *pos, int nbOfPoints,
                                     double eps, std::vector<int>& elts, std::vector<int>& eltsIndex) const;
/// @cond INTERNAL
    static MEDCouplingUMesh *MergeUMeshesLL(std::vector<const MEDCouplingUMesh *>& a) throw(INTERP_KERNEL::Exception);
    typedef int (*DimM1DescNbrer)(int id, unsigned nb, const INTERP_KERNEL::CellModel& cm, bool compute, const int *conn1, const int *conn2);
    MEDCouplingUMesh *buildDescendingConnectivityGen(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx, DimM1DescNbrer nbrer) const throw(INTERP_KERNEL::Exception);
    static void FillInCompact3DMode(int spaceDim, int nbOfNodesInCell, const int *conn, const double *coo, double *zipFrmt) throw(INTERP_KERNEL::Exception);
    static void AppendExtrudedCell(const int *connBg, const int *connEnd, int nbOfNodesPerLev, bool isQuad, std::vector<int>& ret);
    static void IntersectDescending2DMeshes(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2, double eps,
                                            std::vector< std::vector<int> >& intersectEdge1, std::vector< std::vector<int> >& colinear2, std::vector< std::vector<int> >& subDiv2,
                                            MEDCouplingUMesh *& m1Desc, DataArrayInt *&desc1, DataArrayInt *&descIndx1, DataArrayInt *&revDesc1, DataArrayInt *&revDescIndx1,
                                            MEDCouplingUMesh *& m2Desc, DataArrayInt *&desc2, DataArrayInt *&descIndx2, DataArrayInt *&revDesc2, DataArrayInt *&revDescIndx2,
                                            std::vector<double>& addCoo) throw(INTERP_KERNEL::Exception);
    static void BuildIntersectEdges(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2, const std::vector<double>& addCoo, const std::vector< std::vector<int> >& subDiv, std::vector< std::vector<int> >& intersectEdge) throw(INTERP_KERNEL::Exception);
    static void BuildIntersecting2DCellsFromEdges(double eps, const MEDCouplingUMesh *m1, const int *desc1, const int *descIndx1, const std::vector<std::vector<int> >& intesctEdges1, const std::vector< std::vector<int> >& colinear2,
                                                  const MEDCouplingUMesh *m2, const int *desc2, const int *descIndx2, const std::vector<std::vector<int> >& intesctEdges2,
                                                  const std::vector<double>& addCoords,
                                                  std::vector<double>& addCoordsQuadratic, std::vector<int>& cr, std::vector<int>& crI, std::vector<int>& cNb1, std::vector<int>& cNb2);
    static void BuildUnionOf2DMesh(const std::vector<int>& conn2D, const std::vector<int>& connI2D, std::vector<int>& polyUnion);
    static void AssemblyForSplitFrom3DCurve(const std::vector<int>& cut3DCurve, std::vector<int>& nodesOnPlane, const int *nodal3DSurf, const int *nodalIndx3DSurf,
                                              const int *nodal3DCurve, const int *nodalIndx3DCurve,
                                              const int *desc, const int *descIndx, std::vector< std::pair<int,int> >& cut3DSurf) throw(INTERP_KERNEL::Exception);
    void assemblyForSplitFrom3DSurf(const std::vector< std::pair<int,int> >& cut3DSurf,
                                    const int *desc, const int *descIndx, std::vector<int>& nodalRes, std::vector<int>& nodalResIndx, std::vector<int>& cellIds) const throw(INTERP_KERNEL::Exception);
    /// @endcond
  private:
    //! this iterator stores current position in _nodal_connec array.
    mutable int _iterator;
    int _mesh_dim;
    DataArrayInt *_nodal_connec;
    DataArrayInt *_nodal_connec_index;
    std::set<INTERP_KERNEL::NormalizedCellType> _types;
  private:
    static const char PART_OF_NAME[];
  public:
    static double EPS_FOR_POLYH_ORIENTATION;
  };

  class MEDCouplingUMeshCell;

  class MEDCOUPLING_EXPORT MEDCouplingUMeshCellIterator
  {
  public:
    MEDCouplingUMeshCellIterator(MEDCouplingUMesh *mesh);
    MEDCouplingUMeshCellIterator(MEDCouplingUMesh *mesh, MEDCouplingUMeshCell *itc, int bg, int end);
    ~MEDCouplingUMeshCellIterator();
    MEDCouplingUMeshCell *nextt();
  private:
    MEDCouplingUMesh *_mesh;
    MEDCouplingUMeshCell *_cell;
    bool _own_cell;
    int _cell_id;
    int _nb_cell;
  };

  class MEDCouplingUMeshCellByTypeIterator;

  class MEDCOUPLING_EXPORT MEDCouplingUMeshCellByTypeEntry
  {
  public:
    MEDCouplingUMeshCellByTypeEntry(MEDCouplingUMesh *mesh);
    MEDCouplingUMeshCellByTypeIterator *iterator();
    ~MEDCouplingUMeshCellByTypeEntry();
  private:
    MEDCouplingUMesh *_mesh;
  };

  class MEDCOUPLING_EXPORT MEDCouplingUMeshCellEntry
  {
  public:
    MEDCouplingUMeshCellEntry(MEDCouplingUMesh *mesh,  INTERP_KERNEL::NormalizedCellType type, MEDCouplingUMeshCell *itc, int bg, int end);
    ~MEDCouplingUMeshCellEntry();
    INTERP_KERNEL::NormalizedCellType getType() const;
    int getNumberOfElems() const;
    MEDCouplingUMeshCellIterator *iterator();
  private:
    MEDCouplingUMesh *_mesh;
    INTERP_KERNEL::NormalizedCellType _type;
    MEDCouplingUMeshCell *_itc;
    int _bg;
    int _end;
  };

  class MEDCOUPLING_EXPORT MEDCouplingUMeshCellByTypeIterator
  {
  public:
    MEDCouplingUMeshCellByTypeIterator(MEDCouplingUMesh *mesh);
    ~MEDCouplingUMeshCellByTypeIterator();
    MEDCouplingUMeshCellEntry *nextt();
  private:
    MEDCouplingUMesh *_mesh;
    MEDCouplingUMeshCell *_cell;
    int _cell_id;
    int _nb_cell;
  };

  class MEDCOUPLING_EXPORT MEDCouplingUMeshCell
  {
  public:
    MEDCouplingUMeshCell(MEDCouplingUMesh *mesh);
    void next();
    std::string repr() const;
    INTERP_KERNEL::NormalizedCellType getType() const;
    const int *getAllConn(int& lgth) const;
  private:
    int *_conn;
    int *_conn_indx;
    int _conn_lgth;
    static const int NOTICABLE_FIRST_VAL=-7;
  };
}

#endif
