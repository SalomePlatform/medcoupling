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
  class MEDCOUPLING_EXPORT MEDCouplingUMesh : public MEDCouplingPointSet
  {
  public:
    static MEDCouplingUMesh *New();
    static MEDCouplingUMesh *New(const char *meshName, int meshDim);
    MEDCouplingUMesh *clone(bool recDeepCpy) const;
    void updateTime();
    MEDCouplingMeshType getType() const { return UNSTRUCTURED; }
    bool isEqual(const MEDCouplingMesh *other, double prec) const;
    void checkCoherency() const throw(INTERP_KERNEL::Exception);
    void setMeshDimension(int meshDim);
    void allocateCells(int nbOfCells);
    void insertNextCell(INTERP_KERNEL::NormalizedCellType type, int size, const int *nodalConnOfCell);
    void finishInsertingCells();
    const std::set<INTERP_KERNEL::NormalizedCellType>& getAllTypes() const { return _types; }
    void setConnectivity(DataArrayInt *conn, DataArrayInt *connIndex, bool isComputingTypes=true);
    DataArrayInt *getNodalConnectivity() const { return _nodal_connec; }
    DataArrayInt *getNodalConnectivityIndex() const { return _nodal_connec_index; }
    INTERP_KERNEL::NormalizedCellType getTypeOfCell(int cellId) const;
    void getNodeIdsOfCell(int cellId, std::vector<int>& conn) const;
    void getCoordinatesOfNode(int nodeId, std::vector<double>& coo) const;
    int getNumberOfNodesInCell(int cellId) const;
    int getNumberOfCells() const;
    int getMeshDimension() const;
    int getMeshLength() const;
    //! size of returned tinyInfo must be always the same.
    void getTinySerializationInformation(std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const;
    bool isEmptyMesh(const std::vector<int>& tinyInfo) const;
    void resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings);
    void serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const;
    void unserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2,
                         const std::vector<std::string>& littleStrings);
    //tools
    void convertToPolyTypes(const std::vector<int>& cellIdsToConvert);
    DataArrayInt *zipCoordsTraducer();
    void getReverseNodalConnectivity(DataArrayInt *revNodal, DataArrayInt *revNodalIndx) const;
    MEDCouplingUMesh *buildDescendingConnectivity(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx) const;
    DataArrayInt *mergeNodes(double precision, bool& areNodesMerged);
    MEDCouplingPointSet *buildPartOfMySelf(const int *start, const int *end, bool keepCoords) const;
    MEDCouplingPointSet *buildPartOfMySelfNode(const int *start, const int *end, bool fullyIn) const;
    void findBoundaryNodes(std::vector<int>& nodes) const;
    MEDCouplingPointSet *buildBoundaryMesh(bool keepCoords) const;
    void renumberNodes(const int *newNodeNumbers, int newNbOfNodes);
    void giveElemsInBoundingBox(const double *bbox, double eps, std::vector<int>& elems);
    MEDCouplingFieldDouble *getMeasureField(bool isAbs) const;
    MEDCouplingFieldDouble *getMeasureFieldOnNode(bool isAbs) const;
    MEDCouplingFieldDouble *buildOrthogonalField() const;
    int getElementContainingPoint(const double *pos, double eps) const;
    void getElementsContainingPoint(const double *pos, double eps, std::vector<int>& elts) const;
    void getElementsContainingPoints(const double *pos, int nbOfPoints, double eps,
                                     std::vector<int>& elts, std::vector<int>& eltsIndex) const;
    void checkButterflyCells(std::vector<int>& cells) const;
    void getBoundingBoxForBBTree(std::vector<double>& bbox) const;
    bool checkConsecutiveCellTypes() const;
    MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const;
    DataArrayDouble *getBarycenterAndOwner() const;
    static MEDCouplingUMesh *mergeUMeshes(const MEDCouplingUMesh *mesh1, const MEDCouplingUMesh *mesh2);
  private:
    MEDCouplingUMesh();
    MEDCouplingUMesh(const MEDCouplingUMesh& other, bool deepCpy);
    ~MEDCouplingUMesh();
    void computeTypes();
    void checkFullyDefined() const throw(INTERP_KERNEL::Exception);
    //tools
    MEDCouplingUMesh *buildPartOfMySelfKeepCoords(const int *start, const int *end) const;
    template<int SPACEDIM>
    void getElementsContainingPointsAlg(const double *coords, const double *pos, int nbOfPoints,
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
