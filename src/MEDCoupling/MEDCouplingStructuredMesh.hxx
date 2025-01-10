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

#ifndef __PARAMEDMEM_MEDCOUPLINGSTRUCTUREDMESH_HXX__
#define __PARAMEDMEM_MEDCOUPLINGSTRUCTUREDMESH_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingMesh.hxx"

namespace MEDCoupling
{
  class MEDCoupling1SGTUMesh;

  class MEDCouplingStructuredMesh : public MEDCouplingMesh
  {
  public:
    MEDCOUPLING_EXPORT INTERP_KERNEL::NormalizedCellType getTypeOfCell(mcIdType cellId) const;
    MEDCOUPLING_EXPORT std::set<INTERP_KERNEL::NormalizedCellType> getAllGeoTypes() const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const;
    MEDCOUPLING_EXPORT DataArrayIdType *giveCellsWithType(INTERP_KERNEL::NormalizedCellType type) const;
    MEDCOUPLING_EXPORT DataArrayIdType *computeNbOfNodesPerCell() const;
    MEDCOUPLING_EXPORT DataArrayIdType *computeNbOfFacesPerCell() const;
    MEDCOUPLING_EXPORT DataArrayIdType *computeEffectiveNbOfNodesPerCell() const;
    MEDCOUPLING_EXPORT std::vector<mcIdType> getLocationFromCellId(mcIdType cellId) const;
    MEDCOUPLING_EXPORT std::vector<mcIdType> getLocationFromNodeId(mcIdType nodeId) const;
    MEDCOUPLING_EXPORT static void GetPosFromId(mcIdType eltId, int meshDim, const mcIdType *split, mcIdType *res);
    MEDCOUPLING_EXPORT static INTERP_KERNEL::NormalizedCellType GetGeoTypeGivenMeshDimension( int meshDim);
    MEDCOUPLING_EXPORT void getNodeIdsOfCell(mcIdType cellId, std::vector<mcIdType>& conn) const;
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT void copyTinyStringsFrom(const MEDCouplingMesh *other);
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const;
    //tools
    MEDCOUPLING_EXPORT std::vector<mcIdType> getDistributionOfTypes() const;
    MEDCOUPLING_EXPORT DataArrayIdType *checkTypeConsistencyAndContig(const std::vector<mcIdType>& code, const std::vector<const DataArrayIdType *>& idsPerType) const;
    MEDCOUPLING_EXPORT void splitProfilePerType(const DataArrayIdType *profile, std::vector<mcIdType>& code, std::vector<DataArrayIdType *>& idsInPflPerType, std::vector<DataArrayIdType *>& idsPerType, bool smartPflKiller=true) const;
    MEDCOUPLING_EXPORT MEDCoupling1SGTUMesh *build1SGTUnstructured() const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildUnstructured() const;
    MEDCOUPLING_EXPORT MEDCouplingMesh *buildPart(const mcIdType *start, const mcIdType *end) const;
    MEDCOUPLING_EXPORT MEDCouplingMesh *buildPartAndReduceNodes(const mcIdType *start, const mcIdType *end, DataArrayIdType*& arr) const;
    MEDCOUPLING_EXPORT DataArrayIdType *simplexize(int policy);
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildOrthogonalField() const;
    MEDCOUPLING_EXPORT void getReverseNodalConnectivity(DataArrayIdType *revNodal, DataArrayIdType *revNodalIndx) const;
    //some useful methods
    MEDCOUPLING_EXPORT MEDCoupling1SGTUMesh *build1SGTSubLevelMesh() const;
    MEDCOUPLING_EXPORT mcIdType getCellIdFromPos(mcIdType i, mcIdType j, mcIdType k) const;
    MEDCOUPLING_EXPORT mcIdType getNodeIdFromPos(mcIdType i, mcIdType j, mcIdType k) const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfCells() const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfNodes() const;
    MEDCOUPLING_EXPORT int getMeshDimension() const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfCellsOfSubLevelMesh() const;
    MEDCOUPLING_EXPORT int getSpaceDimensionOnNodeStruct() const;
    MEDCOUPLING_EXPORT virtual void getNodeGridStructure(mcIdType *res) const = 0;
    MEDCOUPLING_EXPORT virtual void getSplitCellValues(mcIdType *res) const;
    MEDCOUPLING_EXPORT virtual void getSplitNodeValues(mcIdType *res) const;
    MEDCOUPLING_EXPORT virtual std::vector<mcIdType> getNodeGridStructure() const = 0;
    MEDCOUPLING_EXPORT std::vector<mcIdType> getCellGridStructure() const;
    MEDCOUPLING_EXPORT double computeSquareness() const;
    MEDCOUPLING_EXPORT virtual MEDCouplingStructuredMesh *buildStructuredSubPart(const std::vector< std::pair<mcIdType,mcIdType> >& cellPart) const = 0;
    MEDCOUPLING_EXPORT static std::vector<mcIdType> GetSplitVectFromStruct(const std::vector<mcIdType>& strct);
    MEDCOUPLING_EXPORT static bool IsPartStructured(const mcIdType *startIds, const mcIdType *stopIds, const std::vector<mcIdType>& st, std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat);
    MEDCOUPLING_EXPORT static std::vector<mcIdType> GetDimensionsFromCompactFrmt(const std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat);
    MEDCOUPLING_EXPORT static std::vector< std::pair<mcIdType,mcIdType> > GetCompactFrmtFromDimensions(const std::vector<mcIdType>& dims);
    MEDCOUPLING_EXPORT static std::vector< std::pair<mcIdType,mcIdType> > IntersectRanges(const std::vector< std::pair<mcIdType,mcIdType> >& r1, const std::vector< std::pair<mcIdType,mcIdType> >& r2);
    MEDCOUPLING_EXPORT static bool AreRangesIntersect(const std::vector< std::pair<mcIdType,mcIdType> >& r1, const std::vector< std::pair<mcIdType,mcIdType> >& r2);
    MEDCOUPLING_EXPORT static void SwitchOnIdsFrom(const std::vector<mcIdType>& st, const std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat, std::vector<bool>& vectToSwitchOn);
    MEDCOUPLING_EXPORT static void ExtractFieldOfBoolFrom(const std::vector<mcIdType>& st, const std::vector<bool>& fieldOfBool, const std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat, std::vector<bool>& fieldOut);
    MEDCOUPLING_EXPORT static DataArrayDouble *ExtractFieldOfDoubleFrom(const std::vector<mcIdType>& st, const DataArrayDouble *fieldOfDbl, const std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat);
    MEDCOUPLING_EXPORT static void AssignPartOfFieldOfDoubleUsing(const std::vector<mcIdType>& st, DataArrayDouble *fieldOfDbl, const std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat, const DataArrayDouble *other);
    MEDCOUPLING_EXPORT static void ChangeReferenceFromGlobalOfCompactFrmt(const std::vector< std::pair<mcIdType,mcIdType> >& bigInAbs, const std::vector< std::pair<mcIdType,mcIdType> >& partOfBigInAbs, std::vector< std::pair<mcIdType,mcIdType> >& partOfBigRelativeToBig, bool check=true);
    MEDCOUPLING_EXPORT static void ChangeReferenceToGlobalOfCompactFrmt(const std::vector< std::pair<mcIdType,mcIdType> >& bigInAbs, const std::vector< std::pair<mcIdType,mcIdType> >& partOfBigRelativeToBig, std::vector< std::pair<mcIdType,mcIdType> >& partOfBigInAbs, bool check=true);
    MEDCOUPLING_EXPORT static std::vector< std::pair<mcIdType,mcIdType> > TranslateCompactFrmt(const std::vector< std::pair<mcIdType,mcIdType> >& part, const std::vector<mcIdType>& translation);
    MEDCOUPLING_EXPORT static std::vector<mcIdType> FindTranslationFrom(const std::vector< std::pair<mcIdType,mcIdType> >& startingFrom, const std::vector< std::pair<mcIdType,mcIdType> >& goingTo);
    MEDCOUPLING_EXPORT static DataArrayIdType *BuildExplicitIdsFrom(const std::vector<mcIdType>& st, const std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat);
    MEDCOUPLING_EXPORT static void MultiplyPartOf(const std::vector<mcIdType>& st, const std::vector< std::pair<mcIdType,mcIdType> >& part, double factor, DataArrayDouble *da);
    MEDCOUPLING_EXPORT static void MultiplyPartOfByGhost(const std::vector<mcIdType>& st, const std::vector< std::pair<mcIdType,mcIdType> >& part, mcIdType ghostSize, double factor, DataArrayDouble *da);
    MEDCOUPLING_EXPORT static void PutInGhostFormat(mcIdType ghostSize, const std::vector<mcIdType>& st, const std::vector< std::pair<mcIdType,mcIdType> >& part, std::vector<mcIdType>& stWithGhost, std::vector< std::pair<mcIdType,mcIdType> >&partWithGhost);
    MEDCOUPLING_EXPORT static void ApplyGhostOnCompactFrmt(std::vector< std::pair<mcIdType,mcIdType> >& partBeforeFact, mcIdType ghostSize);
    MEDCOUPLING_EXPORT static DataArrayIdType *Build1GTNodalConnectivity(const mcIdType *nodeStBg, const mcIdType *nodeStEnd);
    MEDCOUPLING_EXPORT static DataArrayIdType *Build1GTNodalConnectivityOfSubLevelMesh(const mcIdType *nodeStBg, const mcIdType *nodeStEnd);
    MEDCOUPLING_EXPORT static DataArrayIdType *ComputeCornersGhost(const std::vector<mcIdType>& st, mcIdType ghostLev);
    MEDCOUPLING_EXPORT static mcIdType DeduceNumberOfGivenRangeInCompactFrmt(const std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat);
    MEDCOUPLING_EXPORT static mcIdType DeduceNumberOfGivenStructure(const std::vector<mcIdType>& st);
    MEDCOUPLING_EXPORT static void FindTheWidestAxisOfGivenRangeInCompactFrmt(const std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat, int& axisId, mcIdType& sizeOfRange);
    MEDCOUPLING_EXPORT static mcIdType FindMinimalPartOf(mcIdType minPatchLgth, const std::vector<mcIdType>& st, const std::vector<bool>& crit, std::vector<bool>& reducedCrit, std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat);
    MEDCOUPLING_EXPORT static std::vector< std::vector<mcIdType> > ComputeSignaturePerAxisOf(const std::vector<mcIdType>& st, const std::vector<bool>& crit);
  private:
    static mcIdType GetNumberOfCellsOfSubLevelMesh(const std::vector<mcIdType>& cgs, int mdim);
    static void GetReverseNodalConnectivity1(const std::vector<mcIdType>& ngs, DataArrayIdType *revNodal, DataArrayIdType *revNodalIndx);
    static void GetReverseNodalConnectivity2(const std::vector<mcIdType>& ngs, DataArrayIdType *revNodal, DataArrayIdType *revNodalIndx);
    static void GetReverseNodalConnectivity3(const std::vector<mcIdType>& ngs, DataArrayIdType *revNodal, DataArrayIdType *revNodalIndx);
    static DataArrayIdType *Build1GTNodalConnectivity1D(const mcIdType *nodeStBg);
    static DataArrayIdType *Build1GTNodalConnectivity2D(const mcIdType *nodeStBg);
    static DataArrayIdType *Build1GTNodalConnectivity3D(const mcIdType *nodeStBg);
    static DataArrayIdType *Build1GTNodalConnectivityOfSubLevelMesh2D(const mcIdType *nodeStBg);
    static DataArrayIdType *Build1GTNodalConnectivityOfSubLevelMesh3D(const mcIdType *nodeStBg);
    static mcIdType FindMinimalPartOf1D(const std::vector<mcIdType>& st, const std::vector<bool>& crit, std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat);
    static mcIdType FindMinimalPartOf2D(const std::vector<mcIdType>& st, const std::vector<bool>& crit, std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat);
    static mcIdType FindMinimalPartOf3D(const std::vector<mcIdType>& st, const std::vector<bool>& crit, std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat);
  protected:
    static int ZipNodeStructure(const mcIdType *nodeStBg, const mcIdType *nodeStEnd, mcIdType zipNodeSt[3]);
  protected:
    MEDCOUPLING_EXPORT MEDCouplingStructuredMesh();
    MEDCOUPLING_EXPORT MEDCouplingStructuredMesh(const MEDCouplingStructuredMesh& other, bool deepCpy);
    MEDCOUPLING_EXPORT ~MEDCouplingStructuredMesh();
  };
}

#endif
