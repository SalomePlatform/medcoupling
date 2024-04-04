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

#ifndef __PARAMEDMEM_MEDCOUPLING1GTUMESH_HXX__
#define __PARAMEDMEM_MEDCOUPLING1GTUMESH_HXX__

#include "InterpKernelHashMap.hxx"
#include "MEDCoupling.hxx"
#include "MCType.hxx"
#include "MEDCouplingMesh.hxx"
#include "MEDCouplingPointSet.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MCAuto.hxx"

#include "CellModel.hxx"
#include <string>
#include "NormalizedGeometricTypes"
#include <set>
#include <vector>
#include <ostream>
#include <cstddef>
#include "MEDCouplingRefCountObject.hxx"
#include <map>

namespace MEDCoupling
{
  class MEDCoupling1GTUUMeshCellIterator;

  class MEDCoupling1GTUMesh : public MEDCouplingPointSet
  {
  public:
    MEDCOUPLING_EXPORT static MEDCoupling1GTUMesh *New(const std::string& name, INTERP_KERNEL::NormalizedCellType type);
    MEDCOUPLING_EXPORT static MEDCoupling1GTUMesh *New(const MEDCouplingUMesh *m);
    MEDCOUPLING_EXPORT const INTERP_KERNEL::CellModel& getCellModel() const;
    MEDCOUPLING_EXPORT INTERP_KERNEL::NormalizedCellType getCellModelEnum() const;
    MEDCOUPLING_EXPORT int getMeshDimension() const override;
    MEDCOUPLING_EXPORT DataArrayIdType *giveCellsWithType(INTERP_KERNEL::NormalizedCellType type) const override;
    MEDCOUPLING_EXPORT mcIdType getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const override;
    MEDCOUPLING_EXPORT INTERP_KERNEL::NormalizedCellType getTypeOfCell(mcIdType cellId) const override;
    MEDCOUPLING_EXPORT std::set<INTERP_KERNEL::NormalizedCellType> getAllGeoTypes() const override;
    MEDCOUPLING_EXPORT std::vector<mcIdType> getDistributionOfTypes() const override;
    MEDCOUPLING_EXPORT void splitProfilePerType(const DataArrayIdType *profile, std::vector<mcIdType>& code, std::vector<DataArrayIdType *>& idsInPflPerType, std::vector<DataArrayIdType *>& idsPerType, bool smartPflKiller=true) const override;
    MEDCOUPLING_EXPORT DataArrayIdType *checkTypeConsistencyAndContig(const std::vector<mcIdType>& code, const std::vector<const DataArrayIdType *>& idsPerType) const override;
    MEDCOUPLING_EXPORT void writeVTKLL(std::ostream& ofs, const std::string& cellData, const std::string& pointData, DataArrayByte *byteData) const override;
    MEDCOUPLING_EXPORT std::string getVTKDataSetType() const override;
    MEDCOUPLING_EXPORT std::string getVTKFileExtension() const override;
    //
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const override;
    MEDCOUPLING_EXPORT mcIdType getNodalConnectivityLength() const;
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const override;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const override;
    MEDCOUPLING_EXPORT void checkConsistencyLight() const override;
    MEDCOUPLING_EXPORT DataArrayDouble *computeCellCenterOfMass() const override;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureField(bool isAbs) const override;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureFieldOnNode(bool isAbs) const override;
    MEDCOUPLING_EXPORT mcIdType getCellContainingPoint(const double *pos, double eps) const override;
    MEDCOUPLING_EXPORT void getCellsContainingPoint(const double *pos, double eps, std::vector<mcIdType>& elts) const override;
    MEDCOUPLING_EXPORT void getCellsContainingPoints(const double *pos, mcIdType nbOfPoints, double eps, MCAuto<DataArrayIdType>& elts, MCAuto<DataArrayIdType>& eltsIndex) const override;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildOrthogonalField() const override;
    MEDCOUPLING_EXPORT DataArrayIdType *getCellsInBoundingBox(const double *bbox, double eps) const override;
    MEDCOUPLING_EXPORT DataArrayIdType *getCellsInBoundingBox(const INTERP_KERNEL::DirectedBoundingBox& bbox, double eps) override;
    MEDCOUPLING_EXPORT MEDCouplingPointSet *buildFacePartOfMySelfNode(const mcIdType *start, const mcIdType *end, bool fullyIn) const override;
    MEDCOUPLING_EXPORT DataArrayIdType *findBoundaryNodes() const override;
    MEDCOUPLING_EXPORT MEDCouplingPointSet *buildBoundaryMesh(bool keepCoords) const override;
    MEDCOUPLING_EXPORT void findCommonCells(int compType, mcIdType startCellId, DataArrayIdType *& commonCellsArr, DataArrayIdType *& commonCellsIArr) const override;
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *AggregateOnSameCoordsToUMesh(const std::vector< const MEDCoupling1GTUMesh *>& parts);
  public:
    MEDCOUPLING_EXPORT virtual void allocateCells(mcIdType nbOfCells=0) = 0;
    MEDCOUPLING_EXPORT virtual void insertNextCell(const mcIdType *nodalConnOfCellBg, const mcIdType *nodalConnOfCellEnd) = 0;
    MEDCOUPLING_EXPORT virtual DataArrayIdType *getNodalConnectivity() const = 0;
    MEDCOUPLING_EXPORT virtual void checkConsistencyOfConnectivity() const = 0;
  protected:
    MEDCoupling1GTUMesh(const std::string& name, const INTERP_KERNEL::CellModel& cm);
    MEDCoupling1GTUMesh(const MEDCoupling1GTUMesh& other, bool recDeepCpy);
    MEDCoupling1GTUMesh();
  protected:
    const INTERP_KERNEL::CellModel *_cm;
  };

  class MEDCoupling1DGTUMesh;
  class MEDCouplingCMesh;

  class MEDCoupling1SGTUMesh : public MEDCoupling1GTUMesh
  {
  public:
    MEDCOUPLING_EXPORT static MEDCoupling1SGTUMesh *New(const std::string& name, INTERP_KERNEL::NormalizedCellType type);
    MEDCOUPLING_EXPORT static MEDCoupling1SGTUMesh *New(const MEDCouplingUMesh *m);
    //! useless constructor only for CORBA -> not swigged
    MEDCOUPLING_EXPORT static MEDCoupling1SGTUMesh *New();
    MEDCOUPLING_EXPORT std::string getClassName() const override { return std::string("MEDCoupling1SGTUMesh"); }
    // Copy methods
    MEDCOUPLING_EXPORT MEDCoupling1SGTUMesh *clone(bool recDeepCpy) const override;
    MEDCOUPLING_EXPORT MEDCoupling1SGTUMesh *deepCopy() const override;
    MEDCOUPLING_EXPORT MEDCoupling1SGTUMesh *deepCopyConnectivityOnly() const override;
    // overload of TimeLabel and RefCountObject
    MEDCOUPLING_EXPORT void updateTime() const override;
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const override;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const override;
    // overload of MEDCouplingMesh
    MEDCOUPLING_EXPORT MEDCouplingMeshType getType() const override { return SINGLE_STATIC_GEO_TYPE_UNSTRUCTURED; }

    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const override;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const override;
    MEDCOUPLING_EXPORT void checkFastEquivalWith(const MEDCouplingMesh *other, double prec) const override;
    MEDCOUPLING_EXPORT void checkConsistencyLight() const override;
    MEDCOUPLING_EXPORT void checkConsistency(double eps=1e-12) const override;
    MEDCOUPLING_EXPORT mcIdType getNumberOfCells() const override;
    MEDCOUPLING_EXPORT DataArrayIdType *computeNbOfNodesPerCell() const override;
    MEDCOUPLING_EXPORT DataArrayIdType *computeNbOfFacesPerCell() const override;
    MEDCOUPLING_EXPORT DataArrayIdType *computeEffectiveNbOfNodesPerCell() const override;
    MEDCOUPLING_EXPORT void getNodeIdsOfCell(mcIdType cellId, std::vector<mcIdType>& conn) const override;
    MEDCOUPLING_EXPORT std::string simpleRepr() const override;
    MEDCOUPLING_EXPORT std::string advancedRepr() const override;
    MEDCOUPLING_EXPORT DataArrayDouble *computeIsoBarycenterOfNodesPerCell() const override;
    MEDCOUPLING_EXPORT void renumberCells(const mcIdType *old2NewBg, bool check=true) override;
    MEDCOUPLING_EXPORT MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const override;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildUnstructured() const override;
    MEDCOUPLING_EXPORT DataArrayIdType *simplexize(int policy) override;
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const override;
    // overload of MEDCouplingPointSet
    MEDCOUPLING_EXPORT void shallowCopyConnectivityFrom(const MEDCouplingPointSet *other) override;
    MEDCOUPLING_EXPORT MEDCouplingPointSet *mergeMyselfWithOnSameCoords(const MEDCouplingPointSet *other) const override;
    MEDCOUPLING_EXPORT MEDCouplingPointSet *buildPartOfMySelfKeepCoords(const mcIdType *begin, const mcIdType *end) const override;
    MEDCOUPLING_EXPORT MEDCouplingPointSet *buildPartOfMySelfKeepCoordsSlice(mcIdType start, mcIdType end, mcIdType step) const override;
    MEDCOUPLING_EXPORT void computeNodeIdsAlg(std::vector<bool>& nodeIdsInUse) const override;
    MEDCOUPLING_EXPORT void getReverseNodalConnectivity(DataArrayIdType *revNodal, DataArrayIdType *revNodalIndx) const override;
    MEDCOUPLING_EXPORT void checkFullyDefined() const override;
    MEDCOUPLING_EXPORT bool isEmptyMesh(const std::vector<mcIdType>& tinyInfo) const override;
    MEDCOUPLING_EXPORT DataArrayIdType *computeFetchedNodeIds() const override;
    MEDCOUPLING_EXPORT DataArrayIdType *getNodeIdsInUse(mcIdType& nbrOfNodesInUse) const override;
    MEDCOUPLING_EXPORT void renumberNodesWithOffsetInConn(mcIdType offset) override;
    MEDCOUPLING_EXPORT void renumberNodesInConn(const INTERP_KERNEL::HashMap<mcIdType,mcIdType>& newNodeNumbersO2N) override;
    MEDCOUPLING_EXPORT void renumberNodesInConn(const std::map<mcIdType,mcIdType>& newNodeNumbersO2N) override;
    MEDCOUPLING_EXPORT void renumberNodesInConn(const mcIdType *newNodeNumbersO2N) override;
    MEDCOUPLING_EXPORT void fillCellIdsToKeepFromNodeIds(const mcIdType *begin, const mcIdType *end, bool fullyIn, DataArrayIdType *&cellIdsKeptArr) const override;
    MEDCOUPLING_EXPORT mcIdType getNumberOfNodesInCell(mcIdType cellId) const override;
    MEDCOUPLING_EXPORT DataArrayDouble *getBoundingBoxForBBTree(double arcDetEps=1e-12) const override;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *computeDiameterField() const override;
    MEDCOUPLING_EXPORT void invertOrientationOfAllCells() override;
    // overload of MEDCoupling1GTUMesh
    MEDCOUPLING_EXPORT void checkConsistencyOfConnectivity() const override;
    MEDCOUPLING_EXPORT void allocateCells(mcIdType nbOfCells=0) override;
    MEDCOUPLING_EXPORT void insertNextCell(const mcIdType *nodalConnOfCellBg, const mcIdType *nodalConnOfCellEnd) override;
  public://specific
    MEDCOUPLING_EXPORT void setNodalConnectivity(DataArrayIdType *nodalConn);
    MEDCOUPLING_EXPORT DataArrayIdType *getNodalConnectivity() const override;
    MEDCOUPLING_EXPORT mcIdType getNumberOfNodesPerCell() const;
    MEDCOUPLING_EXPORT static MEDCoupling1SGTUMesh *Merge1SGTUMeshes(const MEDCoupling1SGTUMesh *mesh1, const MEDCoupling1SGTUMesh *mesh2);
    MEDCOUPLING_EXPORT static MEDCoupling1SGTUMesh *Merge1SGTUMeshes(std::vector<const MEDCoupling1SGTUMesh *>& a);
    MEDCOUPLING_EXPORT static MEDCoupling1SGTUMesh *Merge1SGTUMeshesOnSameCoords(std::vector<const MEDCoupling1SGTUMesh *>& a);
    MEDCOUPLING_EXPORT MEDCoupling1SGTUMesh *buildSetInstanceFromThis(std::size_t spaceDim) const;
    MEDCOUPLING_EXPORT MEDCoupling1GTUMesh *computeDualMesh() const;
    MEDCOUPLING_EXPORT DataArrayIdType *sortHexa8EachOther();
    MEDCOUPLING_EXPORT MEDCoupling1SGTUMesh *explodeEachHexa8To6Quad4() const;
    MEDCOUPLING_EXPORT MCAuto<DataArrayDouble> computeTriangleHeight() const;
    MEDCOUPLING_EXPORT MEDCouplingCMesh *structurizeMe(DataArrayIdType *& cellPerm, DataArrayIdType *& nodePerm, double eps=1e-12) const;
  public://serialization
    MEDCOUPLING_EXPORT void getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<mcIdType>& tinyInfo, std::vector<std::string>& littleStrings) const override;
    MEDCOUPLING_EXPORT void resizeForUnserialization(const std::vector<mcIdType>& tinyInfo, DataArrayIdType *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const override;
    MEDCOUPLING_EXPORT void serialize(DataArrayIdType *&a1, DataArrayDouble *&a2) const override;
    MEDCOUPLING_EXPORT void unserialization(const std::vector<double>& tinyInfoD, const std::vector<mcIdType>& tinyInfo, const DataArrayIdType *a1, DataArrayDouble *a2,
                                            const std::vector<std::string>& littleStrings) override;
  private:
    MEDCoupling1SGTUMesh(const std::string& name, const INTERP_KERNEL::CellModel& cm);
    MEDCoupling1SGTUMesh(const MEDCoupling1SGTUMesh& other, bool recDeepCpy);
    MEDCoupling1SGTUMesh();
  private:
    void checkNonDynamicGeoType() const;
    static MEDCoupling1SGTUMesh *Merge1SGTUMeshesLL(std::vector<const MEDCoupling1SGTUMesh *>& a);
    DataArrayIdType *simplexizePol0();
    DataArrayIdType *simplexizePol1();
    DataArrayIdType *simplexizePlanarFace5();
    DataArrayIdType *simplexizePlanarFace6();
    MEDCoupling1DGTUMesh *computeDualMesh3D() const;
    MEDCoupling1DGTUMesh *computeDualMesh2D() const;
    template<class MAPCLS>
    void renumberNodesInConnT(const MAPCLS& newNodeNumbersO2N);
  private:
    MCAuto<DataArrayIdType> _conn;
  public:
    static const int HEXA8_FACE_PAIRS[6];
  };

  class MEDCoupling1DGTUMesh : public MEDCoupling1GTUMesh
  {
  public:
    MEDCOUPLING_EXPORT static MEDCoupling1DGTUMesh *New(const std::string& name, INTERP_KERNEL::NormalizedCellType type);
    MEDCOUPLING_EXPORT static MEDCoupling1DGTUMesh *New(const MEDCouplingUMesh *m);
    //! useless constructor only for CORBA -> not swigged
    MEDCOUPLING_EXPORT static MEDCoupling1DGTUMesh *New();
    MEDCOUPLING_EXPORT std::string getClassName() const override { return std::string("MEDCoupling1DGTUMesh"); }
    // Copy methods
    MEDCOUPLING_EXPORT MEDCoupling1DGTUMesh *clone(bool recDeepCpy) const override;
    MEDCOUPLING_EXPORT MEDCoupling1DGTUMesh *deepCopy() const override;
    MEDCOUPLING_EXPORT MEDCoupling1DGTUMesh *deepCopyConnectivityOnly() const override;

    // overload of TimeLabel and RefCountObject
    MEDCOUPLING_EXPORT void updateTime() const override;
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const override;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const override;
    // overload of MEDCouplingMesh
    MEDCOUPLING_EXPORT MEDCouplingMeshType getType() const override { return SINGLE_DYNAMIC_GEO_TYPE_UNSTRUCTURED; }
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const override;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const override;
    MEDCOUPLING_EXPORT void checkFastEquivalWith(const MEDCouplingMesh *other, double prec) const override;
    MEDCOUPLING_EXPORT void checkConsistencyLight() const override;
    MEDCOUPLING_EXPORT void checkConsistency(double eps=1e-12) const override;
    MEDCOUPLING_EXPORT mcIdType getNumberOfCells() const override;
    MEDCOUPLING_EXPORT DataArrayIdType *computeNbOfNodesPerCell() const override;
    MEDCOUPLING_EXPORT DataArrayIdType *computeNbOfFacesPerCell() const override;
    MEDCOUPLING_EXPORT DataArrayIdType *computeEffectiveNbOfNodesPerCell() const override;
    MEDCOUPLING_EXPORT void getNodeIdsOfCell(mcIdType cellId, std::vector<mcIdType>& conn) const override;
    MEDCOUPLING_EXPORT std::string simpleRepr() const override;
    MEDCOUPLING_EXPORT std::string advancedRepr() const override;
    MEDCOUPLING_EXPORT DataArrayDouble *computeIsoBarycenterOfNodesPerCell() const override;
    MEDCOUPLING_EXPORT void renumberCells(const mcIdType *old2NewBg, bool check=true) override;
    MEDCOUPLING_EXPORT MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const override;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildUnstructured() const override;
    MEDCOUPLING_EXPORT DataArrayIdType *simplexize(int policy) override;
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const override;
    // overload of MEDCouplingPointSet
    MEDCOUPLING_EXPORT void shallowCopyConnectivityFrom(const MEDCouplingPointSet *other) override;
    MEDCOUPLING_EXPORT MEDCouplingPointSet *mergeMyselfWithOnSameCoords(const MEDCouplingPointSet *other) const override;
    MEDCOUPLING_EXPORT MEDCouplingPointSet *buildPartOfMySelfKeepCoords(const mcIdType *begin, const mcIdType *end) const override;
    MEDCOUPLING_EXPORT MEDCouplingPointSet *buildPartOfMySelfKeepCoordsSlice(mcIdType start, mcIdType end, mcIdType step) const override;
    MEDCOUPLING_EXPORT void computeNodeIdsAlg(std::vector<bool>& nodeIdsInUse) const override;
    MEDCOUPLING_EXPORT void getReverseNodalConnectivity(DataArrayIdType *revNodal, DataArrayIdType *revNodalIndx) const override;
    MEDCOUPLING_EXPORT void checkFullyDefined() const override;
    MEDCOUPLING_EXPORT bool isEmptyMesh(const std::vector<mcIdType>& tinyInfo) const override;
    MEDCOUPLING_EXPORT DataArrayIdType *computeFetchedNodeIds() const override;
    MEDCOUPLING_EXPORT DataArrayIdType *getNodeIdsInUse(mcIdType& nbrOfNodesInUse) const override;
    MEDCOUPLING_EXPORT void renumberNodesWithOffsetInConn(mcIdType offset) override;
    MEDCOUPLING_EXPORT void renumberNodesInConn(const INTERP_KERNEL::HashMap<mcIdType,mcIdType>& newNodeNumbersO2N) override;
    MEDCOUPLING_EXPORT void renumberNodesInConn(const std::map<mcIdType,mcIdType>& newNodeNumbersO2N) override;
    MEDCOUPLING_EXPORT void renumberNodesInConn(const mcIdType *newNodeNumbersO2N) override;
    MEDCOUPLING_EXPORT void fillCellIdsToKeepFromNodeIds(const mcIdType *begin, const mcIdType *end, bool fullyIn, DataArrayIdType *&cellIdsKeptArr) const override;
    MEDCOUPLING_EXPORT mcIdType getNumberOfNodesInCell(mcIdType cellId) const override;
    MEDCOUPLING_EXPORT DataArrayDouble *getBoundingBoxForBBTree(double arcDetEps=1e-12) const override;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *computeDiameterField() const override;
    MEDCOUPLING_EXPORT void invertOrientationOfAllCells() override;
    // overload of MEDCoupling1GTUMesh
    MEDCOUPLING_EXPORT void checkConsistencyOfConnectivity() const override;
    MEDCOUPLING_EXPORT void allocateCells(mcIdType nbOfCells=0) override;
    MEDCOUPLING_EXPORT void insertNextCell(const mcIdType *nodalConnOfCellBg, const mcIdType *nodalConnOfCellEnd) override;
  public://specific
    MEDCOUPLING_EXPORT void setNodalConnectivity(DataArrayIdType *nodalConn, DataArrayIdType *nodalConnIndex);
    MEDCOUPLING_EXPORT DataArrayIdType *getNodalConnectivity() const override;
    MEDCOUPLING_EXPORT DataArrayIdType *getNodalConnectivityIndex() const;
    MEDCOUPLING_EXPORT MEDCoupling1DGTUMesh *copyWithNodalConnectivityPacked(bool& isShallowCpyOfNodalConnn) const;
    MEDCOUPLING_EXPORT bool retrievePackedNodalConnectivity(DataArrayIdType *&nodalConn, DataArrayIdType *&nodalConnIndx) const;
    MEDCOUPLING_EXPORT bool isPacked() const;
    MEDCOUPLING_EXPORT static MEDCoupling1DGTUMesh *Merge1DGTUMeshes(const MEDCoupling1DGTUMesh *mesh1, const MEDCoupling1DGTUMesh *mesh2);
    MEDCOUPLING_EXPORT static MEDCoupling1DGTUMesh *Merge1DGTUMeshes(std::vector<const MEDCoupling1DGTUMesh *>& a);
    MEDCOUPLING_EXPORT static MEDCoupling1DGTUMesh *Merge1DGTUMeshesOnSameCoords(std::vector<const MEDCoupling1DGTUMesh *>& a);
    MEDCOUPLING_EXPORT static DataArrayIdType *AggregateNodalConnAndShiftNodeIds(const std::vector<const DataArrayIdType *>& nodalConns, const std::vector<mcIdType>& offsetInNodeIdsPerElt);
    MEDCOUPLING_EXPORT static std::vector<mcIdType> BuildAPolygonFromParts(const std::vector< std::vector<mcIdType> >& parts);
    MEDCOUPLING_EXPORT MEDCoupling1DGTUMesh *buildSetInstanceFromThis(std::size_t spaceDim) const;
  public://serialization
    MEDCOUPLING_EXPORT void getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<mcIdType>& tinyInfo, std::vector<std::string>& littleStrings) const override;
    MEDCOUPLING_EXPORT void resizeForUnserialization(const std::vector<mcIdType>& tinyInfo, DataArrayIdType *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const override;
    MEDCOUPLING_EXPORT void serialize(DataArrayIdType *&a1, DataArrayDouble *&a2) const override;
    MEDCOUPLING_EXPORT void unserialization(const std::vector<double>& tinyInfoD, const std::vector<mcIdType>& tinyInfo, const DataArrayIdType *a1, DataArrayDouble *a2,
                                            const std::vector<std::string>& littleStrings) override;
  private:
    MEDCoupling1DGTUMesh(const std::string& name, const INTERP_KERNEL::CellModel& cm);
    MEDCoupling1DGTUMesh(const MEDCoupling1DGTUMesh& other, bool recDeepCpy);
    MEDCoupling1DGTUMesh();
  private:
    void checkDynamicGeoT2ype() const;
    static MEDCoupling1DGTUMesh *Merge1DGTUMeshesLL(std::vector<const MEDCoupling1DGTUMesh *>& a);
    template<class MAPCLS>
    void renumberNodesInConnT(const MAPCLS& newNodeNumbersO2N);
  private:
    MCAuto<DataArrayIdType> _conn_indx;
    MCAuto<DataArrayIdType> _conn;
  };
}

#endif
