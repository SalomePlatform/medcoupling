// Copyright (C) 2007-2019  CEA/DEN, EDF R&D
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

#ifndef __PARAMEDMEM_MEDCOUPLINGPOINTSET_HXX__
#define __PARAMEDMEM_MEDCOUPLINGPOINTSET_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingMesh.hxx"

#include "InterpKernelHashMap.hxx"

#include <vector>

namespace INTERP_KERNEL
{
  class DirectedBoundingBox;
}

namespace MEDCoupling
{
  class DataArrayIdType;
  class DataArrayDouble;
  
  /*!
   * This class is abstract and not instanciable.
   * MEDCoupling::MEDCouplingUMesh class inherits from this class.
   * This class aggregates an array '_coords' containing nodes coordinates.
   * So all operations on coordinates are managed by this class.
   * This is the case for example for following methods :
   * rotation, translation, scaling, getNodeIdsNearPoint, boundingbox...
   */
  class MEDCouplingPointSet : public MEDCouplingMesh
  {
  protected:
    MEDCOUPLING_EXPORT MEDCouplingPointSet();
    MEDCOUPLING_EXPORT MEDCouplingPointSet(const MEDCouplingPointSet& other, bool deepCpy);
    MEDCOUPLING_EXPORT ~MEDCouplingPointSet();
  public:
    MEDCOUPLING_EXPORT void updateTime() const;
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfNodes() const;
    MEDCOUPLING_EXPORT int getSpaceDimension() const;
    MEDCOUPLING_EXPORT void setCoords(const DataArrayDouble *coords);
    MEDCOUPLING_EXPORT const DataArrayDouble *getCoords() const { return _coords; }
    MEDCOUPLING_EXPORT DataArrayDouble *getCoords() { return _coords; }
    MEDCOUPLING_EXPORT DataArrayDouble *getCoordinatesAndOwner() const;
    MEDCOUPLING_EXPORT const DataArrayDouble *getDirectAccessOfCoordsArrIfInStructure() const { return _coords; }
    MEDCOUPLING_EXPORT void copyTinyStringsFrom(const MEDCouplingMesh *other);
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const;
    MEDCOUPLING_EXPORT void checkFastEquivalWith(const MEDCouplingMesh *other, double prec) const;
    MEDCOUPLING_EXPORT void checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                 DataArrayIdType *&cellCor, DataArrayIdType *&nodeCor) const;
    MEDCOUPLING_EXPORT void checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                            DataArrayIdType *&cellCor) const;
    MEDCOUPLING_EXPORT bool areCoordsEqualIfNotWhy(const MEDCouplingPointSet& other, double prec, std::string& reason) const;
    MEDCOUPLING_EXPORT bool areCoordsEqual(const MEDCouplingPointSet& other, double prec) const;
    MEDCOUPLING_EXPORT bool areCoordsEqualWithoutConsideringStr(const MEDCouplingPointSet& other, double prec) const;
    MEDCOUPLING_EXPORT virtual MEDCouplingPointSet *deepCopyConnectivityOnly() const = 0;
    MEDCOUPLING_EXPORT virtual void shallowCopyConnectivityFrom(const MEDCouplingPointSet *other) = 0;
    MEDCOUPLING_EXPORT virtual DataArrayIdType *mergeNodes(double precision, bool& areNodesMerged, mcIdType& newNbOfNodes);
    MEDCOUPLING_EXPORT virtual DataArrayIdType *mergeNodesCenter(double precision, bool& areNodesMerged, mcIdType& newNbOfNodes);
    MEDCOUPLING_EXPORT virtual MEDCouplingPointSet *mergeMyselfWithOnSameCoords(const MEDCouplingPointSet *other) const = 0;
    MEDCOUPLING_EXPORT virtual void computeNodeIdsAlg(std::vector<bool>& nodeIdsInUse) const = 0;
    MEDCOUPLING_EXPORT void getCoordinatesOfNode(mcIdType nodeId, std::vector<double>& coo) const;
    MEDCOUPLING_EXPORT DataArrayIdType *buildPermArrayForMergeNode(double precision, mcIdType limitNodeId, bool& areNodesMerged, mcIdType& newNbOfNodes) const;
    MEDCOUPLING_EXPORT DataArrayIdType *getNodeIdsNearPoint(const double *pos, double eps) const;
    MEDCOUPLING_EXPORT void getNodeIdsNearPoints(const double *pos, mcIdType nbOfPoints, double eps, DataArrayIdType *& c, DataArrayIdType *& cI) const;
    MEDCOUPLING_EXPORT void findCommonNodes(double prec, mcIdType limitNodeId, DataArrayIdType *&comm, DataArrayIdType *&commIndex) const;
    MEDCOUPLING_EXPORT virtual void findCommonCells(int compType, mcIdType startCellId, DataArrayIdType *& commonCellsArr, DataArrayIdType *& commonCellsIArr) const = 0;
    MEDCOUPLING_EXPORT DataArrayIdType *buildNewNumberingFromCommonNodesFormat(const DataArrayIdType *comm, const DataArrayIdType *commIndex,
                                                                            mcIdType& newNbOfNodes) const;
    MEDCOUPLING_EXPORT void getBoundingBox(double *bbox) const;
    MEDCOUPLING_EXPORT void zipCoords();
    MEDCOUPLING_EXPORT double getCaracteristicDimension() const;
    MEDCOUPLING_EXPORT void recenterForMaxPrecision(double eps);
    MEDCOUPLING_EXPORT void rotate(const double *center, const double *vector, double angle);
    MEDCOUPLING_EXPORT void translate(const double *vector);
    MEDCOUPLING_EXPORT void scale(const double *point, double factor);
    MEDCOUPLING_EXPORT void changeSpaceDimension(int newSpaceDim, double dftVal=0.);
    MEDCOUPLING_EXPORT void tryToShareSameCoords(const MEDCouplingPointSet& other, double epsilon);
    MEDCOUPLING_EXPORT void duplicateNodesInCoords(const mcIdType *nodeIdsToDuplicateBg, const mcIdType *nodeIdsToDuplicateEnd);
    MEDCOUPLING_EXPORT virtual void tryToShareSameCoordsPermute(const MEDCouplingPointSet& other, double epsilon);
    MEDCOUPLING_EXPORT void findNodesOnPlane(const double *pt, const double *vec, double eps, std::vector<mcIdType>& nodes) const;
    MEDCOUPLING_EXPORT void findNodesOnLine(const double *pt, const double *vec, double eps, std::vector<mcIdType>& nodes) const;
    MEDCOUPLING_EXPORT static DataArrayDouble *MergeNodesArray(const MEDCouplingPointSet *m1, const MEDCouplingPointSet *m2);
    MEDCOUPLING_EXPORT static DataArrayDouble *MergeNodesArray(const std::vector<const MEDCouplingPointSet *>& ms);
    MEDCOUPLING_EXPORT static MEDCouplingPointSet *BuildInstanceFromMeshType(MEDCouplingMeshType type);
    MEDCOUPLING_EXPORT static DataArrayIdType *ComputeNbOfInteractionsWithSrcCells(const MEDCouplingPointSet *srcMesh, const MEDCouplingPointSet *trgMesh, double eps);
    MEDCOUPLING_EXPORT MEDCouplingMesh *buildPart(const mcIdType *start, const mcIdType *end) const;
    MEDCOUPLING_EXPORT MEDCouplingMesh *buildPartAndReduceNodes(const mcIdType *start, const mcIdType *end, DataArrayIdType*& arr) const;
    MEDCOUPLING_EXPORT MEDCouplingMesh *buildPartRange(mcIdType beginCellIds, mcIdType endCellIds, mcIdType stepCellIds) const;
    MEDCOUPLING_EXPORT MEDCouplingMesh *buildPartRangeAndReduceNodes(mcIdType beginCellIds, mcIdType endCellIds, mcIdType stepCellIds, mcIdType& beginOut, mcIdType& endOut, mcIdType& stepOut, DataArrayIdType*& arr) const;
    MEDCOUPLING_EXPORT DataArrayIdType *getCellIdsFullyIncludedInNodeIds(const mcIdType *partBg, const mcIdType *partEnd) const;
    MEDCOUPLING_EXPORT DataArrayIdType *getCellIdsLyingOnNodes(const mcIdType *begin, const mcIdType *end, bool fullyIn) const;
    MEDCOUPLING_EXPORT virtual MEDCouplingPointSet *buildPartOfMySelf(const mcIdType *start, const mcIdType *end, bool keepCoords=true) const;
    MEDCOUPLING_EXPORT virtual MEDCouplingPointSet *buildPartOfMySelfSlice(mcIdType start, mcIdType end, mcIdType step, bool keepCoords=true) const;
    MEDCOUPLING_EXPORT virtual MEDCouplingPointSet *buildPartOfMySelfKeepCoords(const mcIdType *begin, const mcIdType *end) const = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingPointSet *buildPartOfMySelfKeepCoordsSlice(mcIdType start, mcIdType end, mcIdType step) const = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingPointSet *buildPartOfMySelfNode(const mcIdType *start, const mcIdType *end, bool fullyIn) const;
    MEDCOUPLING_EXPORT virtual MEDCouplingPointSet *buildFacePartOfMySelfNode(const mcIdType *start, const mcIdType *end, bool fullyIn) const = 0;
    MEDCOUPLING_EXPORT virtual DataArrayIdType *findBoundaryNodes() const = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingPointSet *buildBoundaryMesh(bool keepCoords) const = 0;
    MEDCOUPLING_EXPORT virtual mcIdType getNumberOfNodesInCell(mcIdType cellId) const = 0;
    MEDCOUPLING_EXPORT virtual DataArrayIdType *computeFetchedNodeIds() const = 0;
    MEDCOUPLING_EXPORT virtual DataArrayIdType *getNodeIdsInUse(mcIdType& nbrOfNodesInUse) const = 0;
    MEDCOUPLING_EXPORT virtual void fillCellIdsToKeepFromNodeIds(const mcIdType *begin, const mcIdType *end, bool fullyIn, DataArrayIdType *&cellIdsKeptArr) const = 0;
    MEDCOUPLING_EXPORT virtual void renumberNodesInConn(const mcIdType *newNodeNumbersO2N) = 0;
    MEDCOUPLING_EXPORT virtual void renumberNodesInConn(const INTERP_KERNEL::HashMap<mcIdType,mcIdType>& newNodeNumbersO2N) = 0;
    MEDCOUPLING_EXPORT virtual void renumberNodesInConn(const std::map<mcIdType,mcIdType>& newNodeNumbersO2N) = 0;
    MEDCOUPLING_EXPORT virtual void renumberNodesWithOffsetInConn(mcIdType offset) = 0;
    MEDCOUPLING_EXPORT virtual void renumberNodes(const mcIdType *newNodeNumbers, mcIdType newNbOfNodes);
    MEDCOUPLING_EXPORT virtual void renumberNodesCenter(const mcIdType *newNodeNumbers, mcIdType newNbOfNodes);
    MEDCOUPLING_EXPORT virtual bool isEmptyMesh(const std::vector<mcIdType>& tinyInfo) const = 0;
    MEDCOUPLING_EXPORT virtual void invertOrientationOfAllCells() = 0;
    MEDCOUPLING_EXPORT virtual void checkFullyDefined() const = 0;
    MEDCOUPLING_EXPORT void getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<mcIdType>& tinyInfo, std::vector<std::string>& littleStrings) const;
    MEDCOUPLING_EXPORT void resizeForUnserialization(const std::vector<mcIdType>& tinyInfo, DataArrayIdType *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const;
    MEDCOUPLING_EXPORT void serialize(DataArrayIdType *&a1, DataArrayDouble *&a2) const;
    MEDCOUPLING_EXPORT void unserialization(const std::vector<double>& tinyInfoD, const std::vector<mcIdType>& tinyInfo, const DataArrayIdType *a1, DataArrayDouble *a2,
                                            const std::vector<std::string>& littleStrings);
    MEDCOUPLING_EXPORT virtual DataArrayDouble *getBoundingBoxForBBTree(double arcDetEps=1e-12) const = 0;
    MEDCOUPLING_EXPORT virtual DataArrayIdType *getCellsInBoundingBox(const double *bbox, double eps) const = 0;
    MEDCOUPLING_EXPORT virtual DataArrayIdType *getCellsInBoundingBox(const INTERP_KERNEL::DirectedBoundingBox& bbox, double eps) = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingFieldDouble *computeDiameterField() const = 0;
    MEDCOUPLING_EXPORT virtual DataArrayIdType *zipCoordsTraducer();
    MEDCOUPLING_EXPORT virtual DataArrayIdType *zipConnectivityTraducer(int compType, mcIdType startCellId=0);
    MEDCOUPLING_EXPORT virtual bool areAllNodesFetched() const;
    //tools
  public:
    MEDCOUPLING_EXPORT bool areCellsFrom2MeshEqual(const MEDCouplingPointSet *other, mcIdType cellId, double prec) const;
  protected:
    MEDCOUPLING_EXPORT void checkConsistencyLight() const;
    MEDCOUPLING_EXPORT static bool intersectsBoundingBox(const double* bb1, const double* bb2, int dim, double eps);
    MEDCOUPLING_EXPORT static bool intersectsBoundingBox(const INTERP_KERNEL::DirectedBoundingBox& bb1, const double* bb2, int dim, double eps);
    MEDCOUPLING_EXPORT void rotate2D(const double *center, double angle);
    MEDCOUPLING_EXPORT void rotate3D(const double *center, const double *vect, double angle);
    MEDCOUPLING_EXPORT void project2DCellOnXY(const mcIdType *startConn, const mcIdType *endConn, std::vector<double>& res) const;
    MEDCOUPLING_EXPORT static bool isButterfly2DCell(const std::vector<double>& res, bool isQuad, double eps);
  protected:
    DataArrayDouble *_coords;
  };
}

#endif
