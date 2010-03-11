//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
#ifndef __PARAMEDMEM_MEDCOUPLINGUMESH_HXX__
#define __PARAMEDMEM_MEDCOUPLINGUMESH_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingPointSet.hxx"
#include "MEDCouplingMemArray.hxx"

#include <set>

namespace ParaMEDMEM
{
  class MEDCouplingUMesh : public MEDCouplingPointSet
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *New();
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *New(const char *meshName, int meshDim);
    MEDCOUPLING_EXPORT MEDCouplingUMesh *clone(bool recDeepCpy) const;
    MEDCOUPLING_EXPORT void updateTime();
    MEDCOUPLING_EXPORT MEDCouplingMeshType getType() const { return UNSTRUCTURED; }
    MEDCOUPLING_EXPORT bool isEqual(const MEDCouplingMesh *other, double prec) const;
    MEDCOUPLING_EXPORT void checkCoherency() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setMeshDimension(int meshDim);
    MEDCOUPLING_EXPORT void allocateCells(int nbOfCells);
    MEDCOUPLING_EXPORT void insertNextCell(INTERP_KERNEL::NormalizedCellType type, int size, const int *nodalConnOfCell);
    MEDCOUPLING_EXPORT void finishInsertingCells();
    MEDCOUPLING_EXPORT const std::set<INTERP_KERNEL::NormalizedCellType>& getAllTypes() const { return _types; }
    MEDCOUPLING_EXPORT void setConnectivity(DataArrayInt *conn, DataArrayInt *connIndex, bool isComputingTypes=true);
    MEDCOUPLING_EXPORT DataArrayInt *getNodalConnectivity() const { return _nodal_connec; }
    MEDCOUPLING_EXPORT DataArrayInt *getNodalConnectivityIndex() const { return _nodal_connec_index; }
    MEDCOUPLING_EXPORT INTERP_KERNEL::NormalizedCellType getTypeOfCell(int cellId) const;
    MEDCOUPLING_EXPORT void getNodeIdsOfCell(int cellId, std::vector<int>& conn) const;
    MEDCOUPLING_EXPORT void getCoordinatesOfNode(int nodeId, std::vector<double>& coo) const;
    MEDCOUPLING_EXPORT int getNumberOfNodesInCell(int cellId) const;
    MEDCOUPLING_EXPORT int getNumberOfCells() const;
    MEDCOUPLING_EXPORT int getMeshDimension() const;
    MEDCOUPLING_EXPORT int getMeshLength() const;
    //! size of returned tinyInfo must be always the same.
    MEDCOUPLING_EXPORT void getTinySerializationInformation(std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const;
    MEDCOUPLING_EXPORT bool isEmptyMesh(const std::vector<int>& tinyInfo) const;
    MEDCOUPLING_EXPORT void resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings);
    MEDCOUPLING_EXPORT void serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const;
    MEDCOUPLING_EXPORT void unserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, const std::vector<std::string>& littleStrings);
    //tools
    MEDCOUPLING_EXPORT void convertToPolyTypes(const std::vector<int>& cellIdsToConvert);
    MEDCOUPLING_EXPORT DataArrayInt *zipCoordsTraducer();
    MEDCOUPLING_EXPORT void getReverseNodalConnectivity(DataArrayInt *revNodal, DataArrayInt *revNodalIndx) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildDescendingConnectivity(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx) const;
    MEDCOUPLING_EXPORT DataArrayInt *mergeNodes(double precision, bool& areNodesMerged);
    MEDCOUPLING_EXPORT MEDCouplingPointSet *buildPartOfMySelf(const int *start, const int *end, bool keepCoords) const;
    MEDCOUPLING_EXPORT MEDCouplingPointSet *buildPartOfMySelfNode(const int *start, const int *end, bool fullyIn) const;
    MEDCOUPLING_EXPORT void findBoundaryNodes(std::vector<int>& nodes) const;
    MEDCOUPLING_EXPORT MEDCouplingPointSet *buildBoundaryMesh(bool keepCoords) const;
    MEDCOUPLING_EXPORT void renumberNodes(const int *newNodeNumbers, int newNbOfNodes);
    MEDCOUPLING_EXPORT void giveElemsInBoundingBox(const double *bbox, double eps, std::vector<int>& elems);
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureField(bool isAbs) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureFieldOnNode(bool isAbs) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildOrthogonalField() const;
    MEDCOUPLING_EXPORT int getCellContainingPoint(const double *pos, double eps) const;
    MEDCOUPLING_EXPORT void getCellsContainingPoint(const double *pos, double eps, std::vector<int>& elts) const;
    MEDCOUPLING_EXPORT void getCellsContainingPoints(const double *pos, int nbOfPoints, double eps, std::vector<int>& elts, std::vector<int>& eltsIndex) const;
    MEDCOUPLING_EXPORT void checkButterflyCells(std::vector<int>& cells) const;
    MEDCOUPLING_EXPORT void getBoundingBoxForBBTree(std::vector<double>& bbox) const;
    MEDCOUPLING_EXPORT bool checkConsecutiveCellTypes() const;
    MEDCOUPLING_EXPORT MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const;
    MEDCOUPLING_EXPORT DataArrayDouble *getBarycenterAndOwner() const;
    MEDCOUPLING_EXPORT static MEDCouplingUMesh *mergeUMeshes(const MEDCouplingUMesh *mesh1, const MEDCouplingUMesh *mesh2);
  private:
    MEDCouplingUMesh();
    MEDCouplingUMesh(const MEDCouplingUMesh& other, bool deepCpy);
    ~MEDCouplingUMesh();
    void computeTypes();
    void checkFullyDefined() const throw(INTERP_KERNEL::Exception);
    //tools
    MEDCouplingUMesh *buildPartOfMySelfKeepCoords(const int *start, const int *end) const;
    template<int SPACEDIM>
    void getCellsContainingPointsAlg(const double *coords, const double *pos, int nbOfPoints,
                                     double eps, std::vector<int>& elts, std::vector<int>& eltsIndex) const;
  private:
    //! this iterator stores current position in _nodal_connec array.
    mutable int _iterator;
    int _mesh_dim;
    DataArrayInt *_nodal_connec;
    DataArrayInt *_nodal_connec_index;
    std::set<INTERP_KERNEL::NormalizedCellType> _types;
  private:
    static const char PART_OF_NAME[];
  };
}

#endif
