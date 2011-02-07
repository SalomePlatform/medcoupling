//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
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

#ifndef __PARAMEDMEM_MEDCOUPLINGUMESHDESC_HXX__
#define __PARAMEDMEM_MEDCOUPLINGUMESHDESC_HXX__

#include "MEDCouplingPointSet.hxx"
#include "MEDCoupling.hxx"
#include "NormalizedUnstructuredMesh.hxx"

#include <set>

namespace ParaMEDMEM
{
  class MEDCouplingUMeshDesc : public MEDCouplingPointSet
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingUMeshDesc *New();
    MEDCOUPLING_EXPORT static MEDCouplingUMeshDesc *New(const char *meshName, int meshDim);
    MEDCOUPLING_EXPORT MEDCouplingMesh *deepCpy() const;
    MEDCOUPLING_EXPORT void checkCoherency() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                 DataArrayInt *&cellCor, DataArrayInt *&nodeCor) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                            DataArrayInt *&cellCor) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setMeshDimension(unsigned meshDim);
    MEDCOUPLING_EXPORT int getNumberOfCells() const;
    MEDCOUPLING_EXPORT int getNumberOfFaces() const;
    MEDCOUPLING_EXPORT int getCellMeshLength() const;
    MEDCOUPLING_EXPORT int getFaceMeshLength() const;
    MEDCOUPLING_EXPORT int getMeshDimension() const { return _mesh_dim; }
    MEDCOUPLING_EXPORT INTERP_KERNEL::NormalizedCellType getTypeOfCell(int cellId) const;
    MEDCOUPLING_EXPORT int getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const;
    MEDCOUPLING_EXPORT void getNodeIdsOfCell(int cellId, std::vector<int>& conn) const;
    MEDCOUPLING_EXPORT void getCoordinatesOfNode(int nodeId, std::vector<double>& coo) const;
    MEDCOUPLING_EXPORT std::string simpleRepr() const;
    MEDCOUPLING_EXPORT std::string advancedRepr() const;
    MEDCOUPLING_EXPORT MEDCouplingMeshType getType() const { return UNSTRUCTURED_DESC; }
    MEDCOUPLING_EXPORT void setConnectivity(DataArrayInt *descConn, DataArrayInt *descConnIndex, DataArrayInt *nodalFaceConn, DataArrayInt *nodalFaceConnIndx);
    //tools to overload
    MEDCOUPLING_EXPORT void getTinySerializationInformation(std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const;
    MEDCOUPLING_EXPORT bool isEmptyMesh(const std::vector<int>& tinyInfo) const;
    MEDCOUPLING_EXPORT void resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings);
    MEDCOUPLING_EXPORT void serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const;
    MEDCOUPLING_EXPORT void unserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, const std::vector<std::string>& littleStrings);
    MEDCOUPLING_EXPORT void getCellsInBoundingBox(const double *bbox, double eps, std::vector<int>& elems);
    MEDCOUPLING_EXPORT void getCellsInBoundingBox(const INTERP_KERNEL::DirectedBoundingBox &bbox, double eps, std::vector<int>& elems);
    MEDCOUPLING_EXPORT DataArrayInt *mergeNodes(double precision, bool& areNodesMerged, int& newNbOfNodes);
    MEDCOUPLING_EXPORT DataArrayInt *mergeNodes2(double precision, bool& areNodesMerged, int& newNbOfNodes);
    MEDCOUPLING_EXPORT void tryToShareSameCoordsPermute(const MEDCouplingPointSet& other, double epsilon) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingPointSet *buildPartOfMySelf(const int *start, const int *end, bool keepCoords) const;
    MEDCOUPLING_EXPORT MEDCouplingPointSet *buildPartOfMySelfNode(const int *start, const int *end, bool fullyIn) const;
    MEDCOUPLING_EXPORT MEDCouplingPointSet *buildFacePartOfMySelfNode(const int *start, const int *end, bool fullyIn) const;
    MEDCOUPLING_EXPORT DataArrayInt *simplexize(int policy) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void findBoundaryNodes(std::vector<int>& nodes) const;
    MEDCOUPLING_EXPORT MEDCouplingPointSet *buildBoundaryMesh(bool keepCoords) const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildUnstructured() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void renumberCells(const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void renumberNodes(const int *newNodeNumbers, int newNbOfNodes);
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureField(bool isAbs) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureFieldOnNode(bool isAbs) const;
    MEDCOUPLING_EXPORT int getCellContainingPoint(const double *pos, double eps) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildOrthogonalField() const;
    MEDCOUPLING_EXPORT DataArrayInt *zipCoordsTraducer();
    MEDCOUPLING_EXPORT MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const;
    MEDCOUPLING_EXPORT DataArrayDouble *getBarycenterAndOwner() const;
  private:
    MEDCouplingUMeshDesc();
    ~MEDCouplingUMeshDesc();
    void computeTypes();
    void checkFullyDefined() const throw(INTERP_KERNEL::Exception);
  private:
    unsigned _mesh_dim;
    DataArrayInt *_desc_connec;
    DataArrayInt *_desc_connec_index;
    DataArrayInt *_nodal_connec_face;
    DataArrayInt *_nodal_connec_face_index;
    std::set<INTERP_KERNEL::NormalizedCellType> _types;
  };
}

#endif
