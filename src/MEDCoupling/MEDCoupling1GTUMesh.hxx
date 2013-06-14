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

#ifndef __PARAMEDMEM_MEDCOUPLING1GTUMESH_HXX__
#define __PARAMEDMEM_MEDCOUPLING1GTUMESH_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingPointSet.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"

#include "CellModel.hxx"

namespace ParaMEDMEM
{
  class MEDCoupling1GTUUMeshCellIterator;

  class MEDCoupling1GTUMesh : public MEDCouplingPointSet
  {
  public:
    MEDCOUPLING_EXPORT static MEDCoupling1GTUMesh *New(const char *meshName, INTERP_KERNEL::NormalizedCellType type) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT const INTERP_KERNEL::CellModel& getCellModel() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT INTERP_KERNEL::NormalizedCellType getCellModelEnum() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT int getMeshDimension() const;
    MEDCOUPLING_EXPORT DataArrayInt *giveCellsWithType(INTERP_KERNEL::NormalizedCellType type) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT int getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const;
    MEDCOUPLING_EXPORT INTERP_KERNEL::NormalizedCellType getTypeOfCell(int cellId) const;
    MEDCOUPLING_EXPORT std::set<INTERP_KERNEL::NormalizedCellType> getAllGeoTypes() const;
    MEDCOUPLING_EXPORT std::vector<int> getDistributionOfTypes() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void splitProfilePerType(const DataArrayInt *profile, std::vector<int>& code, std::vector<DataArrayInt *>& idsInPflPerType, std::vector<DataArrayInt *>& idsPerType) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *checkTypeConsistencyAndContig(const std::vector<int>& code, const std::vector<const DataArrayInt *>& idsPerType) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void writeVTKLL(std::ostream& ofs, const std::string& cellData, const std::string& pointData) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT std::string getVTKDataSetType() const throw(INTERP_KERNEL::Exception);
  protected:
    const INTERP_KERNEL::CellModel& _cm;
  };

  class MEDCoupling1SGTUMesh : public MEDCoupling1GTUMesh
  {
  public:
    MEDCOUPLING_EXPORT static MEDCoupling1GTUMesh *New();
    MEDCOUPLING_EXPORT static MEDCoupling1GTUMesh *New(const char *meshName, INTERP_KERNEL::NormalizedCellType type) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCoupling1SGTUMesh *clone(bool recDeepCpy) const;
    MEDCOUPLING_EXPORT MEDCouplingMeshType getType() const { return SINGLE_STATIC_GEO_TYPE_UNSTRUCTURED; }
    MEDCOUPLING_EXPORT void updateTime() const;
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySize() const;
    MEDCOUPLING_EXPORT MEDCouplingMesh *deepCpy() const;
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const;
    MEDCOUPLING_EXPORT void checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                 DataArrayInt *&cellCor, DataArrayInt *&nodeCor) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                            DataArrayInt *&cellCor) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void checkFastEquivalWith(const MEDCouplingMesh *other, double prec) const throw(INTERP_KERNEL::Exception);//tony
    MEDCOUPLING_EXPORT void checkCoherency1(double eps=1e-12) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void checkCoherency2(double eps=1e-12) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT int getNumberOfCells() const;
    MEDCOUPLING_EXPORT DataArrayDouble *getBarycenterAndOwner() const;
    MEDCOUPLING_EXPORT DataArrayDouble *computeIsoBarycenterOfNodesPerCell() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *computeNbOfNodesPerCell() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *computeNbOfFacesPerCell() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getNodeIdsOfCell(int cellId, std::vector<int>& conn) const;
    MEDCOUPLING_EXPORT std::string simpleRepr() const;
    MEDCOUPLING_EXPORT std::string advancedRepr() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureField(bool isAbs) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureFieldOnNode(bool isAbs) const;
    MEDCOUPLING_EXPORT int getCellContainingPoint(const double *pos, double eps) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildOrthogonalField() const;
    MEDCOUPLING_EXPORT void renumberCells(const int *old2NewBg, bool check=true) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildUnstructured() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *simplexize(int policy) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const throw(INTERP_KERNEL::Exception);
    //
    MEDCOUPLING_EXPORT DataArrayInt *mergeNodes(double precision, bool& areNodesMerged, int& newNbOfNodes);
    MEDCOUPLING_EXPORT DataArrayInt *mergeNodes2(double precision, bool& areNodesMerged, int& newNbOfNodes);
    MEDCOUPLING_EXPORT void tryToShareSameCoordsPermute(const MEDCouplingPointSet& other, double epsilon) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingPointSet *buildPartOfMySelf(const int *start, const int *end, bool keepCoords=true) const;
    MEDCOUPLING_EXPORT MEDCouplingPointSet *buildPartOfMySelf2(int start, int end, int step, bool keepCoords=true) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingPointSet *buildPartOfMySelfNode(const int *start, const int *end, bool fullyIn) const;
    MEDCOUPLING_EXPORT MEDCouplingPointSet *buildFacePartOfMySelfNode(const int *start, const int *end, bool fullyIn) const;
    MEDCOUPLING_EXPORT DataArrayInt *findBoundaryNodes() const;
    MEDCOUPLING_EXPORT MEDCouplingPointSet *buildBoundaryMesh(bool keepCoords) const;
    MEDCOUPLING_EXPORT DataArrayInt *getCellsInBoundingBox(const double *bbox, double eps) const;
    MEDCOUPLING_EXPORT DataArrayInt *getCellsInBoundingBox(const INTERP_KERNEL::DirectedBoundingBox& bbox, double eps);
    MEDCOUPLING_EXPORT  DataArrayInt *zipCoordsTraducer();
    MEDCOUPLING_EXPORT void checkFullyDefined() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool isEmptyMesh(const std::vector<int>& tinyInfo) const;
  private:
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _conn;
  };
}

#endif
