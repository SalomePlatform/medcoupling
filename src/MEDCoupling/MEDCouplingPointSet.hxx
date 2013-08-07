// Copyright (C) 2007-2013  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (CEA/DEN)

#ifndef __PARAMEDMEM_MEDCOUPLINGPOINTSET_HXX__
#define __PARAMEDMEM_MEDCOUPLINGPOINTSET_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingMesh.hxx"

#include <vector>

namespace INTERP_KERNEL
{
  class DirectedBoundingBox;
}

namespace ParaMEDMEM
{
  class DataArrayInt;
  class DataArrayDouble;
  
  /*!
   * This class is abstract and not instanciable.
   * ParaMEDMEM::MEDCouplingUMesh class inherits from this class.
   * This class aggregates an array '_coords' containing nodes coordinates.
   * So all operations on coordinates are managed by this class.
   * This is the case for example for following methods :
   * rotation, translation, scaling, getNodeIdsNearPoint, boundingbox...
   */
  class MEDCOUPLING_EXPORT MEDCouplingPointSet : public MEDCouplingMesh
  {
  protected:
    MEDCouplingPointSet();
    MEDCouplingPointSet(const MEDCouplingPointSet& other, bool deepCopy);
    ~MEDCouplingPointSet();
  public:
    void updateTime() const;
    std::size_t getHeapMemorySize() const;
    int getNumberOfNodes() const;
    int getSpaceDimension() const;
    void setCoords(const DataArrayDouble *coords);
    const DataArrayDouble *getCoords() const { return _coords; }
    DataArrayDouble *getCoords() { return _coords; }
    DataArrayDouble *getCoordinatesAndOwner() const;
    void copyTinyStringsFrom(const MEDCouplingMesh *other) throw(INTERP_KERNEL::Exception);
    bool isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const throw(INTERP_KERNEL::Exception);
    bool isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const;
    void checkFastEquivalWith(const MEDCouplingMesh *other, double prec) const throw(INTERP_KERNEL::Exception);
    void checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                              DataArrayInt *&cellCor, DataArrayInt *&nodeCor) const throw(INTERP_KERNEL::Exception);
    void checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                         DataArrayInt *&cellCor) const throw(INTERP_KERNEL::Exception);
    bool areCoordsEqualIfNotWhy(const MEDCouplingPointSet& other, double prec, std::string& reason) const;
    bool areCoordsEqual(const MEDCouplingPointSet& other, double prec) const;
    bool areCoordsEqualWithoutConsideringStr(const MEDCouplingPointSet& other, double prec) const;
    virtual MEDCouplingPointSet *deepCpyConnectivityOnly() const throw(INTERP_KERNEL::Exception) = 0;
    virtual void shallowCopyConnectivityFrom(const MEDCouplingPointSet *other) throw(INTERP_KERNEL::Exception) = 0;
    virtual DataArrayInt *mergeNodes(double precision, bool& areNodesMerged, int& newNbOfNodes);
    virtual DataArrayInt *mergeNodes2(double precision, bool& areNodesMerged, int& newNbOfNodes);
    virtual MEDCouplingPointSet *mergeMyselfWithOnSameCoords(const MEDCouplingPointSet *other) const = 0;
    virtual void computeNodeIdsAlg(std::vector<bool>& nodeIdsInUse) const throw(INTERP_KERNEL::Exception) = 0;
    void getCoordinatesOfNode(int nodeId, std::vector<double>& coo) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *buildPermArrayForMergeNode(double precision, int limitNodeId, bool& areNodesMerged, int& newNbOfNodes) const;
    DataArrayInt *getNodeIdsNearPoint(const double *pos, double eps) const throw(INTERP_KERNEL::Exception);
    void getNodeIdsNearPoints(const double *pos, int nbOfPoints, double eps, DataArrayInt *& c, DataArrayInt *& cI) const throw(INTERP_KERNEL::Exception);
    void findCommonNodes(double prec, int limitNodeId, DataArrayInt *&comm, DataArrayInt *&commIndex) const;
    virtual void findCommonCells(int compType, int startCellId, DataArrayInt *& commonCellsArr, DataArrayInt *& commonCellsIArr) const throw(INTERP_KERNEL::Exception) = 0;
    DataArrayInt *buildNewNumberingFromCommonNodesFormat(const DataArrayInt *comm, const DataArrayInt *commIndex,
                                                         int& newNbOfNodes) const;
    void getBoundingBox(double *bbox) const throw(INTERP_KERNEL::Exception);
    void zipCoords();
    double getCaracteristicDimension() const;
    void recenterForMaxPrecision(double eps) throw(INTERP_KERNEL::Exception);
    void rotate(const double *center, const double *vector, double angle);
    void translate(const double *vector);
    void scale(const double *point, double factor);
    void changeSpaceDimension(int newSpaceDim, double dftVal=0.) throw(INTERP_KERNEL::Exception);
    void tryToShareSameCoords(const MEDCouplingPointSet& other, double epsilon) throw(INTERP_KERNEL::Exception);
    void duplicateNodesInCoords(const int *nodeIdsToDuplicateBg, const int *nodeIdsToDuplicateEnd) throw(INTERP_KERNEL::Exception);
    virtual void tryToShareSameCoordsPermute(const MEDCouplingPointSet& other, double epsilon) throw(INTERP_KERNEL::Exception);
    void findNodesOnPlane(const double *pt, const double *vec, double eps, std::vector<int>& nodes) const throw(INTERP_KERNEL::Exception);
    void findNodesOnLine(const double *pt, const double *vec, double eps, std::vector<int>& nodes) const throw(INTERP_KERNEL::Exception);
    static DataArrayDouble *MergeNodesArray(const MEDCouplingPointSet *m1, const MEDCouplingPointSet *m2) throw(INTERP_KERNEL::Exception);
    static DataArrayDouble *MergeNodesArray(const std::vector<const MEDCouplingPointSet *>& ms) throw(INTERP_KERNEL::Exception);
    static MEDCouplingPointSet *BuildInstanceFromMeshType(MEDCouplingMeshType type);
    static void Rotate2DAlg(const double *center, double angle, int nbNodes, double *coords);
    static void Rotate3DAlg(const double *center, const double *vect, double angle, int nbNodes, double *coords);
    MEDCouplingMesh *buildPart(const int *start, const int *end) const;
    MEDCouplingMesh *buildPartAndReduceNodes(const int *start, const int *end, DataArrayInt*& arr) const;
    MEDCouplingMesh *buildPartRange(int beginCellIds, int endCellIds, int stepCellIds) const throw(INTERP_KERNEL::Exception);
    MEDCouplingMesh *buildPartRangeAndReduceNodes(int beginCellIds, int endCellIds, int stepCellIds, int& beginOut, int& endOut, int& stepOut, DataArrayInt*& arr) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getCellIdsFullyIncludedInNodeIds(const int *partBg, const int *partEnd) const;
    DataArrayInt *getCellIdsLyingOnNodes(const int *begin, const int *end, bool fullyIn) const;
    virtual MEDCouplingPointSet *buildPartOfMySelf(const int *start, const int *end, bool keepCoords=true) const;
    virtual MEDCouplingPointSet *buildPartOfMySelf2(int start, int end, int step, bool keepCoords=true) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingPointSet *buildPartOfMySelfKeepCoords(const int *begin, const int *end) const = 0;
    virtual MEDCouplingPointSet *buildPartOfMySelfKeepCoords2(int start, int end, int step) const = 0;
    virtual MEDCouplingPointSet *buildPartOfMySelfNode(const int *start, const int *end, bool fullyIn) const;
    virtual MEDCouplingPointSet *buildFacePartOfMySelfNode(const int *start, const int *end, bool fullyIn) const = 0;
    virtual DataArrayInt *findBoundaryNodes() const = 0;
    virtual MEDCouplingPointSet *buildBoundaryMesh(bool keepCoords) const = 0;
    virtual int getNumberOfNodesInCell(int cellId) const = 0;
    virtual DataArrayInt *getNodeIdsInUse(int& nbrOfNodesInUse) const throw(INTERP_KERNEL::Exception) = 0;
    virtual void fillCellIdsToKeepFromNodeIds(const int *begin, const int *end, bool fullyIn, DataArrayInt *&cellIdsKeptArr) const = 0;
    virtual void renumberNodesInConn(const int *newNodeNumbersO2N) = 0;
    virtual void renumberNodes(const int *newNodeNumbers, int newNbOfNodes);
    virtual void renumberNodes2(const int *newNodeNumbers, int newNbOfNodes);
    virtual bool isEmptyMesh(const std::vector<int>& tinyInfo) const = 0;
    virtual void checkFullyDefined() const throw(INTERP_KERNEL::Exception) = 0;
    void getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const;
    void resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const;
    void serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const;
    void unserialization(const std::vector<double>& tinyInfoD, const std::vector<int>& tinyInfo, const DataArrayInt *a1, DataArrayDouble *a2,
                         const std::vector<std::string>& littleStrings);
    virtual DataArrayInt *getCellsInBoundingBox(const double *bbox, double eps) const = 0;
    virtual DataArrayInt *getCellsInBoundingBox(const INTERP_KERNEL::DirectedBoundingBox& bbox, double eps) = 0;
    virtual DataArrayInt *zipCoordsTraducer() throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *zipConnectivityTraducer(int compType, int startCellId=0) throw(INTERP_KERNEL::Exception);
    virtual void getReverseNodalConnectivity(DataArrayInt *revNodal, DataArrayInt *revNodalIndx) const throw(INTERP_KERNEL::Exception) = 0;
    //tools
  public:
    bool areCellsFrom2MeshEqual(const MEDCouplingPointSet *other, int cellId, double prec) const;
  protected:
    void checkCoherency() const throw(INTERP_KERNEL::Exception);
    static bool intersectsBoundingBox(const double* bb1, const double* bb2, int dim, double eps);
    static bool intersectsBoundingBox(const INTERP_KERNEL::DirectedBoundingBox& bb1, const double* bb2, int dim, double eps);
    void rotate2D(const double *center, double angle);
    void rotate3D(const double *center, const double *vect, double angle);
    void project2DCellOnXY(const int *startConn, const int *endConn, std::vector<double>& res) const;
    static bool isButterfly2DCell(const std::vector<double>& res, bool isQuad, double eps);
  protected:
    DataArrayDouble *_coords;
  };
}

#endif
