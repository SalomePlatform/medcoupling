// Copyright (C) 2007-2014  CEA/DEN, EDF R&D
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

namespace ParaMEDMEM
{
  class MEDCoupling1SGTUMesh;

  class MEDCouplingStructuredMesh : public MEDCouplingMesh
  {
  public:
    MEDCOUPLING_EXPORT INTERP_KERNEL::NormalizedCellType getTypeOfCell(int cellId) const;
    MEDCOUPLING_EXPORT std::set<INTERP_KERNEL::NormalizedCellType> getAllGeoTypes() const;
    MEDCOUPLING_EXPORT int getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const;
    MEDCOUPLING_EXPORT DataArrayInt *giveCellsWithType(INTERP_KERNEL::NormalizedCellType type) const;
    MEDCOUPLING_EXPORT DataArrayInt *computeNbOfNodesPerCell() const;
    MEDCOUPLING_EXPORT DataArrayInt *computeNbOfFacesPerCell() const;
    MEDCOUPLING_EXPORT DataArrayInt *computeEffectiveNbOfNodesPerCell() const;
    MEDCOUPLING_EXPORT static void GetPosFromId(int nodeId, int meshDim, const int *split, int *res);
    MEDCOUPLING_EXPORT static INTERP_KERNEL::NormalizedCellType GetGeoTypeGivenMeshDimension(int meshDim);
    MEDCOUPLING_EXPORT void getNodeIdsOfCell(int cellId, std::vector<int>& conn) const;
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT void copyTinyStringsFrom(const MEDCouplingMesh *other);
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const;
    //tools
    MEDCOUPLING_EXPORT std::vector<int> getDistributionOfTypes() const;
    MEDCOUPLING_EXPORT DataArrayInt *checkTypeConsistencyAndContig(const std::vector<int>& code, const std::vector<const DataArrayInt *>& idsPerType) const;
    MEDCOUPLING_EXPORT void splitProfilePerType(const DataArrayInt *profile, std::vector<int>& code, std::vector<DataArrayInt *>& idsInPflPerType, std::vector<DataArrayInt *>& idsPerType) const;
    MEDCOUPLING_EXPORT MEDCoupling1SGTUMesh *build1SGTUnstructured() const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildUnstructured() const;
    MEDCOUPLING_EXPORT MEDCouplingMesh *buildPart(const int *start, const int *end) const;
    MEDCOUPLING_EXPORT MEDCouplingMesh *buildPartAndReduceNodes(const int *start, const int *end, DataArrayInt*& arr) const;
    MEDCOUPLING_EXPORT DataArrayInt *simplexize(int policy);
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildOrthogonalField() const;
    MEDCOUPLING_EXPORT void getReverseNodalConnectivity(DataArrayInt *revNodal, DataArrayInt *revNodalIndx) const;
    //some useful methods
    MEDCOUPLING_EXPORT MEDCoupling1SGTUMesh *build1SGTSubLevelMesh() const;
    MEDCOUPLING_EXPORT int getCellIdFromPos(int i, int j, int k) const;
    MEDCOUPLING_EXPORT int getNodeIdFromPos(int i, int j, int k) const;
    MEDCOUPLING_EXPORT int getNumberOfCells() const;
    MEDCOUPLING_EXPORT int getNumberOfNodes() const;
    MEDCOUPLING_EXPORT int getMeshDimension() const;
    MEDCOUPLING_EXPORT int getNumberOfCellsOfSubLevelMesh() const;
    MEDCOUPLING_EXPORT int getSpaceDimensionOnNodeStruct() const;
    MEDCOUPLING_EXPORT virtual void getNodeGridStructure(int *res) const = 0;
    MEDCOUPLING_EXPORT virtual void getSplitCellValues(int *res) const;
    MEDCOUPLING_EXPORT virtual void getSplitNodeValues(int *res) const;
    MEDCOUPLING_EXPORT virtual std::vector<int> getNodeGridStructure() const = 0;
    MEDCOUPLING_EXPORT std::vector<int> getCellGridStructure() const;
    MEDCOUPLING_EXPORT virtual MEDCouplingStructuredMesh *buildStructuredSubPart(const std::vector< std::pair<int,int> >& cellPart) const = 0;
    MEDCOUPLING_EXPORT static std::vector<int> GetSplitVectFromStruct(const std::vector<int>& strct);
    MEDCOUPLING_EXPORT static bool IsPartStructured(const int *startIds, const int *stopIds, const std::vector<int>& st, std::vector< std::pair<int,int> >& partCompactFormat);
    MEDCOUPLING_EXPORT static std::vector<int> GetDimensionsFromCompactFrmt(const std::vector< std::pair<int,int> >& partCompactFormat);
    MEDCOUPLING_EXPORT static std::vector< std::pair<int,int> > GetCompactFrmtFromDimensions(const std::vector<int>& dims);
    MEDCOUPLING_EXPORT static std::vector< std::pair<int,int> > IntersectRanges(const std::vector< std::pair<int,int> >& r1, const std::vector< std::pair<int,int> >& r2);
    MEDCOUPLING_EXPORT static void SwitchOnIdsFrom(const std::vector<int>& st, const std::vector< std::pair<int,int> >& partCompactFormat, std::vector<bool>& vectToSwitchOn);
    MEDCOUPLING_EXPORT static void ExtractFieldOfBoolFrom(const std::vector<int>& st, const std::vector<bool>& fieldOfBool, const std::vector< std::pair<int,int> >& partCompactFormat, std::vector<bool>& fieldOut);
    MEDCOUPLING_EXPORT static DataArrayDouble *ExtractFieldOfDoubleFrom(const std::vector<int>& st, const DataArrayDouble *fieldOfDbl, const std::vector< std::pair<int,int> >& partCompactFormat);
    MEDCOUPLING_EXPORT static void ChangeReferenceFromGlobalOfCompactFrmt(const std::vector< std::pair<int,int> >& bigInAbs, const std::vector< std::pair<int,int> >& partOfBigInAbs, std::vector< std::pair<int,int> >& partOfBigRelativeToBig);
    MEDCOUPLING_EXPORT static void ChangeReferenceToGlobalOfCompactFrmt(const std::vector< std::pair<int,int> >& bigInAbs, const std::vector< std::pair<int,int> >& partOfBigRelativeToBig, std::vector< std::pair<int,int> >& partOfBigInAbs);
    MEDCOUPLING_EXPORT static DataArrayInt *BuildExplicitIdsFrom(const std::vector<int>& st, const std::vector< std::pair<int,int> >& partCompactFormat);
    MEDCOUPLING_EXPORT static DataArrayInt *Build1GTNodalConnectivity(const int *nodeStBg, const int *nodeStEnd);
    MEDCOUPLING_EXPORT static DataArrayInt *Build1GTNodalConnectivityOfSubLevelMesh(const int *nodeStBg, const int *nodeStEnd);
    MEDCOUPLING_EXPORT static int DeduceNumberOfGivenRangeInCompactFrmt(const std::vector< std::pair<int,int> >& partCompactFormat);
    MEDCOUPLING_EXPORT static int DeduceNumberOfGivenStructure(const std::vector<int>& st);
    MEDCOUPLING_EXPORT static void FindTheWidestAxisOfGivenRangeInCompactFrmt(const std::vector< std::pair<int,int> >& partCompactFormat, int& axisId, int& sizeOfRange);
    MEDCOUPLING_EXPORT static int FindMinimalPartOf(const std::vector<int>& st, const std::vector<bool>& crit, std::vector<bool>& reducedCrit, std::vector< std::pair<int,int> >& partCompactFormat);
    MEDCOUPLING_EXPORT static std::vector< std::vector<int> > ComputeSignaturePerAxisOf(const std::vector<int>& st, const std::vector<bool>& crit);
  private:
    static int GetNumberOfCellsOfSubLevelMesh(const std::vector<int>& cgs, int mdim);
    static void GetReverseNodalConnectivity1(const std::vector<int>& ngs, DataArrayInt *revNodal, DataArrayInt *revNodalIndx);
    static void GetReverseNodalConnectivity2(const std::vector<int>& ngs, DataArrayInt *revNodal, DataArrayInt *revNodalIndx);
    static void GetReverseNodalConnectivity3(const std::vector<int>& ngs, DataArrayInt *revNodal, DataArrayInt *revNodalIndx);
    static DataArrayInt *Build1GTNodalConnectivity1D(const int *nodeStBg);
    static DataArrayInt *Build1GTNodalConnectivity2D(const int *nodeStBg);
    static DataArrayInt *Build1GTNodalConnectivity3D(const int *nodeStBg);
    static DataArrayInt *Build1GTNodalConnectivityOfSubLevelMesh2D(const int *nodeStBg);
    static DataArrayInt *Build1GTNodalConnectivityOfSubLevelMesh3D(const int *nodeStBg);
    static int FindMinimalPartOf1D(const std::vector<int>& st, const std::vector<bool>& crit, std::vector< std::pair<int,int> >& partCompactFormat);
    static int FindMinimalPartOf2D(const std::vector<int>& st, const std::vector<bool>& crit, std::vector< std::pair<int,int> >& partCompactFormat);
    static int FindMinimalPartOf3D(const std::vector<int>& st, const std::vector<bool>& crit, std::vector< std::pair<int,int> >& partCompactFormat);
  protected:
    static int ZipNodeStructure(const int *nodeStBg, const int *nodeStEnd, int zipNodeSt[3]);
  protected:
    MEDCOUPLING_EXPORT MEDCouplingStructuredMesh();
    MEDCOUPLING_EXPORT MEDCouplingStructuredMesh(const MEDCouplingStructuredMesh& other, bool deepCpy);
    MEDCOUPLING_EXPORT ~MEDCouplingStructuredMesh();
  };
}

#endif
