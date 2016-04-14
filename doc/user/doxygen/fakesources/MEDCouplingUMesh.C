// Copyright (C) 2013-2016  CEA/DEN, EDF R&D, OPEN CASCADE
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

// Copyright (C) 2013  CEA/DEN, EDF R&D, OPEN CASCADE
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

// This file contains some code used only for
// * generation of documentation for inline methods of MEDCouplingUMesh class,
// * groupping methods into "Basic API", "Advanced" and "Others..." sections

namespace MEDCoupling
{
  /*!
   * Returns the nodal connectivity array. For more info on how data in stored in
   * this array, see \ref MEDCouplingUMeshAdvBuild.
   *  \return const DataArrayInt * - a pointer to the nodal connectivity array
   *          referred by \a this mesh.
   */
  const DataArrayInt * MEDCouplingUMesh::getNodalConnectivity() const {}
  /*!
   * Returns the nodal connectivity index array. For more info on how data in stored in
   * this array, see \ref MEDCouplingUMeshAdvBuild.
   *  \return const DataArrayInt * - a pointer to the nodal connectivity index array
   *          referred by \a this mesh.
   */
  const DataArrayInt * MEDCouplingUMesh::getNodalConnectivityIndex() const {}
  /*!
   * Returns the nodal connectivity array. For more info on how data in stored in
   * this array, see \ref MEDCouplingUMeshAdvBuild.
   *  \return const DataArrayInt * - a pointer to the nodal connectivity array
   *          referred by \a this mesh.
   */
  DataArrayInt * MEDCouplingUMesh::getNodalConnectivity() {}
  /*!
   * Returns the nodal connectivity index array. For more info on how data in stored in
   * this array, see \ref MEDCouplingUMeshAdvBuild.
   *  \return const DataArrayInt * - a pointer to the nodal connectivity index array
   *          referred by \a this mesh.
   */
  DataArrayInt * MEDCouplingUMesh::getNodalConnectivityIndex() {}
}

namespace MEDCoupling
{
//================================================================================
/////////////////////// MEDCouplingUMesh GROUPPING ///////////////////////////////
//================================================================================

/*! \name Basic API   */
///@{
MEDCouplingUMesh::FuseUMeshesOnSameCoords(const std::vector<const MEDCouplingUMesh *>& meshes, int compType, std::vector<DataArrayInt *>& corr);
MEDCouplingUMesh::Intersect2DMeshes(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2, double eps, DataArrayInt *&cellNb1, DataArrayInt *&cellNb2);
MEDCouplingUMesh::MergeNodesOnUMeshesSharingSameCoords(const std::vector<MEDCouplingUMesh *>& meshes, double eps);
MEDCouplingUMesh::MergeUMeshes(const MEDCouplingUMesh *mesh1, const MEDCouplingUMesh *mesh2);
MEDCouplingUMesh::MergeUMeshes(std::vector<const MEDCouplingUMesh *>& a);
MEDCouplingUMesh::MergeUMeshesOnSameCoords(const MEDCouplingUMesh *mesh1, const MEDCouplingUMesh *mesh2);
MEDCouplingUMesh::MergeUMeshesOnSameCoords(const std::vector<const MEDCouplingUMesh *>& meshes);
MEDCouplingUMesh::PutUMeshesOnSameAggregatedCoords(const std::vector<MEDCouplingUMesh *>& meshes);
MEDCouplingUMesh::allocateCells(int nbOfCells);
MEDCouplingUMesh::are2DCellsNotCorrectlyOriented(const double *vec, bool polyOnly, std::vector<int>& cells) const;
MEDCouplingUMesh::areCellsIncludedIn(const MEDCouplingUMesh *other, int compType, DataArrayInt *& arr) const;
MEDCouplingUMesh::arePolyhedronsNotCorrectlyOriented(std::vector<int>& cells) const;
MEDCouplingUMesh::buildBoundaryMesh(bool keepCoords) const;
MEDCouplingUMesh::buildDescendingConnectivity(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx) const;
MEDCouplingUMesh::buildDescendingConnectivity2(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx) const;
MEDCouplingUMesh::buildDirectionVectorField() const;
MEDCouplingUMesh::buildFacePartOfMySelfNode(const int *begin, const int *end, bool fullyIn) const;
MEDCouplingUMesh::buildOrthogonalField() const;
MEDCouplingUMesh::buildPartOfMySelf(const int *begin, const int *end, bool keepCoords=true) const;
//MEDCouplingUMesh::buildPartOfMySelfNode(const int *begin, const int *end, bool fullyIn) const;
MEDCouplingUMesh::buildPartOrthogonalField(const int *begin, const int *end) const;
MEDCouplingUMesh::buildSlice3D(const double *origin, const double *vec, double eps, DataArrayInt *&cellIds) const;
MEDCouplingUMesh::buildSlice3DSurf(const double *origin, const double *vec, double eps, DataArrayInt *&cellIds) const;
MEDCouplingUMesh::checkCoherency() const;
MEDCouplingUMesh::checkCoherency1(double eps=1e-12) const;
//MEDCouplingUMesh::checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec,DataArrayInt *&cellCor) const;
//MEDCouplingUMesh::checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec,DataArrayInt *&cellCor, DataArrayInt *&nodeCor) const;
MEDCouplingUMesh::checkFastEquivalWith(const MEDCouplingMesh *other, double prec) const;
MEDCouplingUMesh::clone(bool recDeepCpy) const;
MEDCouplingUMesh::computeIsoBarycenterOfNodesPerCell() const;
MEDCouplingUMesh::convertAllToPoly();
MEDCouplingUMesh::convertQuadraticCellsToLinear();
MEDCouplingUMesh::convertToPolyTypes(const int *cellIdsToConvertBg, const int *cellIdsToConvertEnd);
MEDCouplingUMesh::deepCpy() const;
MEDCouplingUMesh::findAndCorrectBadOriented3DExtrudedCells();
MEDCouplingUMesh::findBoundaryNodes() const;
MEDCouplingUMesh::finishInsertingCells();
MEDCouplingUMesh::getAllGeoTypes() const;
MEDCouplingUMesh::getAspectRatioField() const;
MEDCouplingUMesh::getBarycenterAndOwner() const;
MEDCouplingUMesh::getCellContainingPoint(const double *pos, double eps) const;
MEDCouplingUMesh::getCellIdsCrossingPlane(const double *origin, const double *vec, double eps) const;
//MEDCouplingUMesh::getCellIdsFullyIncludedInNodeIds(const int *partBg, const int *partEnd) const;
//MEDCouplingUMesh::getCellIdsLyingOnNodes(const int *begin, const int *end, bool fullyIn) const;
MEDCouplingUMesh::getCellsContainingPoint(const double *pos, double eps, std::vector<int>& elts) const;
MEDCouplingUMesh::getCellsContainingPoints(const double *pos, int nbOfPoints, double eps, std::vector<int>& elts, std::vector<int>& eltsIndex) const;
MEDCouplingUMesh::getCellsInBoundingBox(const INTERP_KERNEL::DirectedBoundingBox& bbox, double eps);
MEDCouplingUMesh::getCellsInBoundingBox(const double *bbox, double eps) const;
MEDCouplingUMesh::getEdgeRatioField() const;
MEDCouplingUMesh::getMeasureField(bool isAbs) const;
MEDCouplingUMesh::getMeasureFieldOnNode(bool isAbs) const;
MEDCouplingUMesh::getMeshDimension() const;
MEDCouplingUMesh::getNodalConnectivity() const;
MEDCouplingUMesh::getNodalConnectivity();
MEDCouplingUMesh::getNodalConnectivityIndex() const;
MEDCouplingUMesh::getNodalConnectivityIndex();
MEDCouplingUMesh::getNodeIdsInUse(int& nbrOfNodesInUse) const;
MEDCouplingUMesh::getNodeIdsOfCell(int cellId, std::vector<int>& conn) const;
MEDCouplingUMesh::getNumberOfCells() const;
MEDCouplingUMesh::getPartBarycenterAndOwner(const int *begin, const int *end) const;
MEDCouplingUMesh::getPartMeasureField(bool isAbs, const int *begin, const int *end) const;
MEDCouplingUMesh::getReverseNodalConnectivity(DataArrayInt *revNodal, DataArrayInt *revNodalIndx) const;
MEDCouplingUMesh::getSkewField() const;
MEDCouplingUMesh::getTypeOfCell(int cellId) const;
MEDCouplingUMesh::getTypesOfPart(const int *begin, const int *end) const;
MEDCouplingUMesh::getWarpField() const;
MEDCouplingUMesh::insertNextCell(INTERP_KERNEL::NormalizedCellType type, int size, const int *nodalConnOfCell);
MEDCouplingUMesh::isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const;
//MEDCouplingUMesh::mergeNodes(double precision, bool& areNodesMerged, int& newNbOfNodes);
//MEDCouplingUMesh::mergeNodes2(double precision, bool& areNodesMerged, int& newNbOfNodes);
MEDCouplingUMesh::orientCorrectly2DCells(const double *vec, bool polyOnly);
MEDCouplingUMesh::orientCorrectlyPolyhedrons();
//MEDCouplingUMesh::renumberNodes(const int *newNodeNumbers, int newNbOfNodes);
//MEDCouplingUMesh::renumberNodes2(const int *newNodeNumbers, int newNbOfNodes);
MEDCouplingUMesh::renumberNodesInConn(const int *newNodeNumbersO2N);
MEDCouplingUMesh::reprQuickOverview(std::ostream& stream) const;
MEDCouplingUMesh::setConnectivity(DataArrayInt *conn, DataArrayInt *connIndex, bool isComputingTypes=true);
MEDCouplingUMesh::setMeshDimension(int meshDim);
MEDCouplingUMesh::sortCellsInMEDFileFrmt();
//MEDCouplingUMesh::tryToShareSameCoordsPermute(const MEDCouplingPointSet& other, double epsilon);
MEDCouplingUMesh::unPolyze();
//MEDCouplingUMesh::zipConnectivityTraducer(int compType, int startCellId=0);
MEDCouplingUMesh::zipCoordsTraducer();
  ///@} 

  /*! \name Advanced API  */
///@{
MEDCouplingUMesh::areOnlySimplexCells() const;
MEDCouplingUMesh::checkButterflyCells(std::vector<int>& cells, double eps=1e-12) const;
MEDCouplingUMesh::computeTypes();
MEDCouplingUMesh::convertDegeneratedCells();
MEDCouplingUMesh::convertExtrudedPolyhedra();
MEDCouplingUMesh::getMeshLength() const;
MEDCouplingUMesh::isFullyQuadratic() const;
MEDCouplingUMesh::isPresenceOfQuadratic() const;
MEDCouplingUMesh::simplexize(int policy);
MEDCouplingUMesh::tessellate2D(double eps);
  ///@

/*! \name Others... */
///@{
MEDCouplingUMesh::AggregateSortedByTypeMeshesOnSameCoords(const std::vector<const MEDCouplingUMesh *>& ms,DataArrayInt *&szOfCellGrpOfSameType,DataArrayInt *&idInMsOfCellGrpOfSameType);
//MEDCouplingUMesh::AppendExtrudedCell(const int *connBg, const int *connEnd, int nbOfNodesPerLev, bool isQuad, std::vector<int>& ret);
MEDCouplingUMesh::AreCellsEqual(const int *conn, const int *connI, int cell1, int cell2, int compType);
MEDCouplingUMesh::AreCellsEqual0(const int *conn, const int *connI, int cell1, int cell2);
MEDCouplingUMesh::AreCellsEqual1(const int *conn, const int *connI, int cell1, int cell2);
MEDCouplingUMesh::AreCellsEqual2(const int *conn, const int *connI, int cell1, int cell2);
MEDCouplingUMesh::AreCellsEqual3(const int *conn, const int *connI, int cell1, int cell2);
MEDCouplingUMesh::AreCellsEqual7(const int *conn, const int *connI, int cell1, int cell2);
MEDCouplingUMesh::AreCellsEqualInPool(const std::vector<int>& candidates, int compType, const int *conn, const int *connI, DataArrayInt *result) ;
//MEDCouplingUMesh::AssemblyForSplitFrom3DCurve(const std::vector<int>& cut3DCurve, std::vector<int>& nodesOnPlane, const int *nodal3DSurf, const int *nodalIndx3DSurf,const int *nodal3DCurve, const int *nodalIndx3DCurve,const int *desc, const int *descIndx, std::vector< std::pair<int,int> >& cut3DSurf);
MEDCouplingUMesh::Build0DMeshFromCoords(DataArrayDouble *da);
MEDCouplingUMesh::BuildConvexEnvelopOf2DCellJarvis(const double *coords, const int *nodalConnBg, const int *nodalConnEnd, DataArrayInt *nodalConnecOut);
//MEDCouplingUMesh::BuildIntersectEdges(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2, const std::vector<double>& addCoo, const std::vector< std::vector<int> >& subDiv, std::vector< std::vector<int> >& intersectEdge);
//MEDCouplingUMesh::BuildIntersecting2DCellsFromEdges(double eps, const MEDCouplingUMesh *m1, const int *desc1, const int *descIndx1, const std::vector<std::vector<int> >& intesctEdges1, const std::vector< std::vector<int> >& colinear2,const MEDCouplingUMesh *m2, const int *desc2, const int *descIndx2, const std::vector<std::vector<int> >& intesctEdges2,const std::vector<double>& addCoords,std::vector<double>& addCoordsQuadratic, std::vector<int>& cr, std::vector<int>& crI, std::vector<int>& cNb1, std::vector<int>& cNb2);
MEDCouplingUMesh::ComputeNeighborsOfCellsAdv(const DataArrayInt *desc, const DataArrayInt *descI, const DataArrayInt *revDesc, const DataArrayInt *revDescI,DataArrayInt *&neighbors, DataArrayInt *&neighborsIdx);
//MEDCouplingUMesh::ComputeRangesFromTypeDistribution(const std::vector<int>& code);
MEDCouplingUMesh::ComputeSpreadZoneGradually(const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn);
MEDCouplingUMesh::ComputeSpreadZoneGraduallyFromSeed(const int *seedBg, const int *seedEnd, const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn, int nbOfDepthPeeling, int& nbOfDepthPeelingPerformed);
//MEDCouplingUMesh::ComputeSpreadZoneGraduallyFromSeedAlg(std::vector<bool>& fetched, const int *seedBg, const int *seedEnd, const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn, int nbOfDepthPeeling, int& nbOfDepthPeelingPerformed);
MEDCouplingUMesh::ComputeVecAndPtOfFace(double eps, const double *coords, const int *begin, const int *end, double *v, double *p);
//MEDCouplingUMesh::CorrectExtrudedCell(int *begin, int *end);
MEDCouplingUMesh::CorrectExtrudedStaticCell(int *begin, int *end);
MEDCouplingUMesh::ExtractFromIndexedArrays(const int *idsOfSelectBg, const int *idsOfSelectEnd, const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn,DataArrayInt* &arrOut, DataArrayInt* &arrIndexOut);
//MEDCouplingUMesh::FillInCompact3DMode(int spaceDim, int nbOfNodesInCell, const int *conn, const double *coo, double *zipFrmt);
MEDCouplingUMesh::FindCommonCellsAlg(int compType, int startCellId, const DataArrayInt *nodal, const DataArrayInt *nodalI, const DataArrayInt *revNodal, const DataArrayInt *revNodalI,DataArrayInt *& commonCellsArr, DataArrayInt *& commonCellsIArr);
//MEDCouplingUMesh::IntersectDescending2DMeshes(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2, double eps,std::vector< std::vector<int> >& intersectEdge1, std::vector< std::vector<int> >& colinear2, std::vector< std::vector<int> >& subDiv2,MEDCouplingUMesh *& m1Desc, DataArrayInt *&desc1, DataArrayInt *&descIndx1, DataArrayInt *&revDesc1, DataArrayInt *&revDescIndx1,MEDCouplingUMesh *& m2Desc, DataArrayInt *&desc2, DataArrayInt *&descIndx2, DataArrayInt *&revDesc2, DataArrayInt *&revDescIndx2,std::vector<double>& addCoo);
//MEDCouplingUMesh::Is3DExtrudedCellWellOriented(const int *begin, const int *end, const double *coords);
MEDCouplingUMesh::Is3DExtrudedStaticCellWellOriented(const int *begin, const int *end, const double *coords);
MEDCouplingUMesh::IsPolygonWellOriented(bool isQuadratic, const double *vec, const int *begin, const int *end, const double *coords);
MEDCouplingUMesh::IsPolyhedronWellOriented(const int *begin, const int *end, const double *coords);
MEDCouplingUMesh::IsPyra5WellOriented(const int *begin, const int *end, const double *coords);
MEDCouplingUMesh::IsTetra4WellOriented(const int *begin, const int *end, const double *coords);
MEDCouplingUMesh::MEDCouplingUMesh();
MEDCouplingUMesh::MEDCouplingUMesh(const MEDCouplingUMesh& other, bool deepCopy);
//MEDCouplingUMesh::MergeUMeshesLL(std::vector<const MEDCouplingUMesh *>& a);
MEDCouplingUMesh::New();
MEDCouplingUMesh::New(const std::string& meshName, int meshDim);
MEDCouplingUMesh::RemoveIdsFromIndexedArrays(const int *idsToRemoveBg, const int *idsToRemoveEnd, DataArrayInt *arr, DataArrayInt *arrIndx, int offsetForRemoval=0);
MEDCouplingUMesh::SetPartOfIndexedArrays(const int *idsOfSelectBg, const int *idsOfSelectEnd, const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn,const DataArrayInt *srcArr, const DataArrayInt *srcArrIndex,DataArrayInt* &arrOut, DataArrayInt* &arrIndexOut);
MEDCouplingUMesh::SetPartOfIndexedArrays2(int start, int end, int step, const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn,const DataArrayInt *srcArr, const DataArrayInt *srcArrIndex,DataArrayInt* &arrOut, DataArrayInt* &arrIndexOut);
MEDCouplingUMesh::SetPartOfIndexedArraysSameIdx(const int *idsOfSelectBg, const int *idsOfSelectEnd, DataArrayInt *arrInOut, const DataArrayInt *arrIndxIn,const DataArrayInt *srcArr, const DataArrayInt *srcArrIndex);
MEDCouplingUMesh::SetPartOfIndexedArraysSameIdx2(int start, int end, int step, DataArrayInt *arrInOut, const DataArrayInt *arrIndxIn,const DataArrayInt *srcArr, const DataArrayInt *srcArrIndex);
MEDCouplingUMesh::SimplifyPolyhedronCell(double eps, const DataArrayDouble *coords, const int *begin, const int *end, DataArrayInt *res);
MEDCouplingUMesh::TryToCorrectPolyhedronOrientation(int *begin, int *end, const double *coords);
MEDCouplingUMesh::advancedRepr() const;
//MEDCouplingUMesh::areCellsFrom2MeshEqual(const MEDCouplingUMesh *other, int cellId, double prec) const;
MEDCouplingUMesh::areCellsIncludedIn2(const MEDCouplingUMesh *other, DataArrayInt *& arr) const;
//MEDCouplingUMesh::assemblyForSplitFrom3DSurf(const std::vector< std::pair<int,int> >& cut3DSurf,const int *desc, const int *descIndx, DataArrayInt *nodalRes, DataArrayInt *nodalResIndx, DataArrayInt *cellIds) const;
MEDCouplingUMesh::buildExtrudedMesh(const MEDCouplingUMesh *mesh1D, int policy);
MEDCouplingUMesh::buildExtrudedMeshFromThisLowLev(int nbOfNodesOf1Lev, bool isQuad) const;
MEDCouplingUMesh::buildPartOfMySelf2(int start, int end, int step, bool keepCoords=true) const;
MEDCouplingUMesh::buildPartOfMySelfKeepCoords(const int *begin, const int *end) const;
MEDCouplingUMesh::buildPartOfMySelfKeepCoords2(int start, int end, int step) const;
MEDCouplingUMesh::buildSetInstanceFromThis(int spaceDim) const;
MEDCouplingUMesh::buildSpreadZonesWithPoly() const;
MEDCouplingUMesh::buildUnionOf2DMesh() const;
MEDCouplingUMesh::buildUnionOf3DMesh() const;
MEDCouplingUMesh::buildUnstructured() const;
MEDCouplingUMesh::cellIterator();
MEDCouplingUMesh::cellsByType();
MEDCouplingUMesh::checkConnectivityFullyDefined() const;
MEDCouplingUMesh::checkConsecutiveCellTypes() const;
MEDCouplingUMesh::checkConsecutiveCellTypesAndOrder(const INTERP_KERNEL::NormalizedCellType *orderBg, const INTERP_KERNEL::NormalizedCellType *orderEnd) const;
MEDCouplingUMesh::checkConsecutiveCellTypesForMEDFileFrmt() const;
MEDCouplingUMesh::checkFullyDefined() const;
MEDCouplingUMesh::checkTypeConsistencyAndContig(const std::vector<int>& code, const std::vector<const DataArrayInt *>& idsPerType) const;
MEDCouplingUMesh::computeFetchedNodeIds() const;
MEDCouplingUMesh::computeNbOfNodesPerCell() const;
MEDCouplingUMesh::computeNeighborsOfCells(DataArrayInt *&neighbors, DataArrayInt *&neighborsIdx) const;
MEDCouplingUMesh::computeNodeIdsAlg(std::vector<bool>& nodeIdsInUse) const;
MEDCouplingUMesh::computeSkin() const;
MEDCouplingUMesh::convertCellArrayPerGeoType(const DataArrayInt *da) const;
MEDCouplingUMesh::convertLinearCellsToQuadratic(int conversionType=0);
MEDCouplingUMesh::convertLinearCellsToQuadratic1D0(DataArrayInt *&conn, DataArrayInt *&connI, DataArrayDouble *& coords, std::set<INTERP_KERNEL::NormalizedCellType>& types) const;
MEDCouplingUMesh::convertLinearCellsToQuadratic2D0(DataArrayInt *&conn, DataArrayInt *&connI, DataArrayDouble *& coords, std::set<INTERP_KERNEL::NormalizedCellType>& types) const;
MEDCouplingUMesh::convertLinearCellsToQuadratic2D1(DataArrayInt *&conn, DataArrayInt *&connI, DataArrayDouble *& coords, std::set<INTERP_KERNEL::NormalizedCellType>& types) const;
MEDCouplingUMesh::convertLinearCellsToQuadratic2DAnd3D0(const MEDCouplingUMesh *m1D, const DataArrayInt *desc, const DataArrayInt *descI, DataArrayInt *&conn, DataArrayInt *&connI, DataArrayDouble *& coords, std::set<INTERP_KERNEL::NormalizedCellType>& types) const;
MEDCouplingUMesh::convertLinearCellsToQuadratic3D0(DataArrayInt *&conn, DataArrayInt *&connI, DataArrayDouble *& coords, std::set<INTERP_KERNEL::NormalizedCellType>& types) const;
MEDCouplingUMesh::convertLinearCellsToQuadratic3D1(DataArrayInt *&conn, DataArrayInt *&connI, DataArrayDouble *& coords, std::set<INTERP_KERNEL::NormalizedCellType>& types) const;
MEDCouplingUMesh::convexEnvelop2D();
MEDCouplingUMesh::cppRepr() const;
MEDCouplingUMesh::distanceToPoint(const double *ptBg, const double *ptEnd, int& cellId, int& nodeId) const;
//MEDCouplingUMesh::distanceToPoint2DCurveAlg(const double *pt, const DataArrayInt *cellIds, double& ret0, int& cellId) const;
//MEDCouplingUMesh::distanceToPoint3DSurfAlg(const double *pt, const DataArrayInt *cellIds, double& ret0, int& cellId) const;
MEDCouplingUMesh::duplicateNodes(const int *nodeIdsToDuplicateBg, const int *nodeIdsToDuplicateEnd);
MEDCouplingUMesh::duplicateNodesInConn(const int *nodeIdsToDuplicateBg, const int *nodeIdsToDuplicateEnd, int offset);
MEDCouplingUMesh::emulateMEDMEMBDC(const MEDCouplingUMesh *nM1LevMesh, DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *&revDesc, DataArrayInt *&revDescIndx, DataArrayInt *& nM1LevMeshIds, DataArrayInt *&meshnM1Old2New) const;
MEDCouplingUMesh::explode3DMeshTo1D(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx) const;
MEDCouplingUMesh::fillCellIdsToKeepFromNodeIds(const int *begin, const int *end, bool fullyIn, DataArrayInt *&cellIdsKeptArr) const;
MEDCouplingUMesh::fillExtCoordsUsingTranslAndAutoRotation(const MEDCouplingUMesh *mesh1D, bool isQuad) const;
MEDCouplingUMesh::fillExtCoordsUsingTranslAndAutoRotation2D(const MEDCouplingUMesh *mesh1D, bool isQuad) const;
MEDCouplingUMesh::fillExtCoordsUsingTranslAndAutoRotation3D(const MEDCouplingUMesh *mesh1D, bool isQuad) const;
MEDCouplingUMesh::fillExtCoordsUsingTranslation(const MEDCouplingUMesh *mesh1D, bool isQuad) const;
MEDCouplingUMesh::findAndCorrectBadOriented3DCells();
MEDCouplingUMesh::findCellIdsLyingOn(const MEDCouplingUMesh& otherDimM1OnSameCoords, DataArrayInt *&cellIdsRk0, DataArrayInt *&cellIdsRk1) const;
MEDCouplingUMesh::findCellIdsOnBoundary() const;
MEDCouplingUMesh::findCommonCells(int compType, int startCellId, DataArrayInt *& commonCellsArr, DataArrayInt *& commonCellsIArr) const;
MEDCouplingUMesh::findNodesToDuplicate(const MEDCouplingUMesh& otherDimM1OnSameCoords, DataArrayInt *& nodeIdsToDuplicate,DataArrayInt *& cellIdsNeededToBeRenum, DataArrayInt *& cellIdsNotModified) const;
//MEDCouplingUMesh::getAllTypes() const;
MEDCouplingUMesh::getBoundingBoxForBBTree(std::vector<double>& bbox) const;
MEDCouplingUMesh::getDistributionOfTypes() const;
MEDCouplingUMesh::getFastAveragePlaneOfThis(double *vec, double *pos) const;
//MEDCouplingUMesh::getHeapMemorySize() const;
MEDCouplingUMesh::getLevArrPerCellTypes(const INTERP_KERNEL::NormalizedCellType *orderBg, const INTERP_KERNEL::NormalizedCellType *orderEnd, DataArrayInt *&nbPerType) const;
MEDCouplingUMesh::getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const;
MEDCouplingUMesh::getNumberOfNodesInCell(int cellId) const;
MEDCouplingUMesh::getQuadraticStatus() const;
MEDCouplingUMesh::getRenumArrForConsecutiveCellTypesSpec(const INTERP_KERNEL::NormalizedCellType *orderBg, const INTERP_KERNEL::NormalizedCellType *orderEnd) const;
MEDCouplingUMesh::getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const;
MEDCouplingUMesh::getType() const { return UNSTRUCTURED; }
MEDCouplingUMesh::getVTKDataSetType() const;
MEDCouplingUMesh::giveCellsWithType(INTERP_KERNEL::NormalizedCellType type) const;
MEDCouplingUMesh::isContiguous1D() const;
MEDCouplingUMesh::isEmptyMesh(const std::vector<int>& tinyInfo) const;
MEDCouplingUMesh::isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const;
MEDCouplingUMesh::keepCellIdsByType(INTERP_KERNEL::NormalizedCellType type, const int *begin, const int *end) const;
MEDCouplingUMesh::keepSpecifiedCells(INTERP_KERNEL::NormalizedCellType type, const int *idsPerGeoTypeBg, const int *idsPerGeoTypeEnd) const;
MEDCouplingUMesh::mergeMyselfWith(const MEDCouplingMesh *other) const;
MEDCouplingUMesh::partitionBySpreadZone() const;
MEDCouplingUMesh::project1D(const double *pt, const double *v, double eps, double *res) const;
MEDCouplingUMesh::rearrange2ConsecutiveCellTypes();
MEDCouplingUMesh::renumberCells(const int *old2NewBg, bool check=true);
MEDCouplingUMesh::reprConnectivityOfThis() const;
MEDCouplingUMesh::reprConnectivityOfThisLL(std::ostringstream& stream) const;
MEDCouplingUMesh::resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const;
MEDCouplingUMesh::serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const;
MEDCouplingUMesh::setPartOfMySelf(const int *cellIdsBg, const int *cellIdsEnd, const MEDCouplingUMesh& otherOnSameCoordsThanThis);
MEDCouplingUMesh::setPartOfMySelf2(int start, int end, int step, const MEDCouplingUMesh& otherOnSameCoordsThanThis);
MEDCouplingUMesh::shiftNodeNumbersInConn(int delta);
MEDCouplingUMesh::simpleRepr() const;
MEDCouplingUMesh::simplexizePlanarFace5();
MEDCouplingUMesh::simplexizePlanarFace6();
MEDCouplingUMesh::simplexizePol0();
MEDCouplingUMesh::simplexizePol1();
MEDCouplingUMesh::simplifyPolyhedra(double eps);
MEDCouplingUMesh::split3DCurveWithPlane(const double *origin, const double *vec, double eps, std::vector<int>& cut3DCurve);
MEDCouplingUMesh::splitByType() const;
MEDCouplingUMesh::splitProfilePerType(const DataArrayInt *profile, std::vector<int>& code, std::vector<DataArrayInt *>& idsInPflPerType, std::vector<DataArrayInt *>& idsPerType) const;
MEDCouplingUMesh::subDivide2DMesh(const int *nodeSubdived, const int *nodeIndxSubdived, const int *desc, const int *descIndex);
MEDCouplingUMesh::unserialization(const std::vector<double>& tinyInfoD, const std::vector<int>& tinyInfo, const DataArrayInt *a1, DataArrayDouble *a2, const std::vector<std::string>& littleStrings);
MEDCouplingUMesh::updateTime() const;
MEDCouplingUMesh::writeVTKLL(std::ostream& ofs, const std::string& cellData, const std::string& pointData) const;
MEDCouplingUMesh::~MEDCouplingUMesh();
//template<class SonsGenerator> MEDCouplingUMesh * MEDCouplingUMesh::buildDescendingConnectivityGen(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx, DimM1DescNbrer nbrer) const;
template<int SPACEDIM> void MEDCouplingUMesh::getCellsContainingPointsAlg
(const double *coords, const double *pos, int nbOfPoints,double eps, std::vector<int>& elts,
 std::vector<int>& eltsIndex) const;

//const INTERP_KERNEL::NormalizedCellType MEDCouplingUMesh::MEDMEM_ORDER[N_MEDMEM_ORDER];
//const int MEDCouplingUMesh::N_MEDMEM_ORDER=24;
double MEDCouplingUMesh::EPS_FOR_POLYH_ORIENTATION;
int MEDCouplingUMesh::_mesh_dim;
std::set<INTERP_KERNEL::NormalizedCellType> MEDCouplingUMesh::_types;
DataArrayInt * MEDCouplingUMesh::_nodal_connec;
DataArrayInt * MEDCouplingUMesh::_nodal_connec_index;
  ///@} 
}

