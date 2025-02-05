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

#ifndef __PARAMEDMEM_MEDCOUPLINGEXTRUDEDMESH_HXX__
#define __PARAMEDMEM_MEDCOUPLINGEXTRUDEDMESH_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingMesh.hxx"

#include <vector>

namespace MEDCoupling
{
  class DataArrayIdType;
  class DataArrayDouble;
  class MEDCouplingUMesh;
  class MEDCouplingCMesh;
  class MEDCouplingFieldDouble;

  class MEDCouplingMappedExtrudedMesh : public MEDCouplingMesh
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingMappedExtrudedMesh *New(const MEDCouplingUMesh *mesh3D, const MEDCouplingUMesh *mesh2D, mcIdType cell2DId);
    MEDCOUPLING_EXPORT static MEDCouplingMappedExtrudedMesh *New(const MEDCouplingCMesh *mesh3D);
    MEDCOUPLING_EXPORT static MEDCouplingMappedExtrudedMesh *New();
    std::string getClassName() const override { return std::string("MEDCouplingMappedExtrudedMesh"); }
    MEDCOUPLING_EXPORT MEDCouplingMeshType getType() const;
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDCOUPLING_EXPORT void copyTinyStringsFrom(const MEDCouplingMesh *other);
    MEDCOUPLING_EXPORT mcIdType getNumberOfCells() const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfNodes() const;
    MEDCOUPLING_EXPORT int getSpaceDimension() const;
    MEDCOUPLING_EXPORT int getMeshDimension() const;
    MEDCOUPLING_EXPORT MEDCouplingMappedExtrudedMesh *deepCopy() const;
    MEDCOUPLING_EXPORT MEDCouplingMappedExtrudedMesh *clone(bool recDeepCpy) const;
    MEDCOUPLING_EXPORT const DataArrayDouble *getDirectAccessOfCoordsArrIfInStructure() const;
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const;
    MEDCOUPLING_EXPORT void checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                 DataArrayIdType *&cellCor, DataArrayIdType *&nodeCor) const;
    MEDCOUPLING_EXPORT void checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                            DataArrayIdType *&cellCor) const;
    MEDCOUPLING_EXPORT INTERP_KERNEL::NormalizedCellType getTypeOfCell(mcIdType cellId) const;
    MEDCOUPLING_EXPORT std::set<INTERP_KERNEL::NormalizedCellType> getAllGeoTypes() const;
    MEDCOUPLING_EXPORT DataArrayIdType *giveCellsWithType(INTERP_KERNEL::NormalizedCellType type) const;
    MEDCOUPLING_EXPORT DataArrayIdType *computeNbOfNodesPerCell() const;
    MEDCOUPLING_EXPORT DataArrayIdType *computeNbOfFacesPerCell() const;
    MEDCOUPLING_EXPORT DataArrayIdType *computeEffectiveNbOfNodesPerCell() const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const;
    MEDCOUPLING_EXPORT void getNodeIdsOfCell(mcIdType cellId, std::vector<mcIdType>& conn) const;
    MEDCOUPLING_EXPORT void getCoordinatesOfNode(mcIdType nodeId, std::vector<double>& coo) const;
    MEDCOUPLING_EXPORT std::string simpleRepr() const;
    MEDCOUPLING_EXPORT std::string advancedRepr() const;
    MEDCOUPLING_EXPORT void checkConsistencyLight() const;
    MEDCOUPLING_EXPORT void checkConsistency(double eps=1e-12) const;
    MEDCOUPLING_EXPORT void getBoundingBox(double *bbox) const;
    MEDCOUPLING_EXPORT void updateTime() const;
    MEDCOUPLING_EXPORT void renumberCells(const mcIdType *old2NewBg, bool check=true);
    MEDCouplingUMesh *getMesh2D() const { return _mesh2D.iAmATrollConstCast(); }
    MEDCouplingUMesh *getMesh1D() const { return _mesh1D.iAmATrollConstCast(); }
    DataArrayIdType *getMesh3DIds() const { return _mesh3D_ids.iAmATrollConstCast(); }
    MEDCOUPLING_EXPORT MEDCouplingUMesh *build3DUnstructuredMesh() const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildUnstructured() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureField(bool) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureFieldOnNode(bool) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildOrthogonalField() const;
    MEDCOUPLING_EXPORT mcIdType getCellContainingPoint(const double *pos, double eps) const;
    MEDCOUPLING_EXPORT void getCellsContainingPoint(const double *pos, double eps, std::vector<mcIdType>& elts) const;
    MEDCOUPLING_EXPORT static mcIdType FindCorrespCellByNodalConn(const std::vector<mcIdType>& nodalConnec,
                                                             const mcIdType *revNodalPtr, const mcIdType *revNodalIndxPtr);
    MEDCOUPLING_EXPORT static void Project1DMeshes(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2, double eps,
                                                   MEDCouplingUMesh *&m1r, MEDCouplingUMesh *&m2r, double *v);
    MEDCOUPLING_EXPORT void rotate(const double *center, const double *vector, double angle);
    MEDCOUPLING_EXPORT void translate(const double *vector);
    MEDCOUPLING_EXPORT void scale(const double *point, double factor);
    MEDCOUPLING_EXPORT std::vector<mcIdType> getDistributionOfTypes() const;
    MEDCOUPLING_EXPORT DataArrayIdType *checkTypeConsistencyAndContig(const std::vector<mcIdType>& code, const std::vector<const DataArrayIdType *>& idsPerType) const;
    MEDCOUPLING_EXPORT void splitProfilePerType(const DataArrayIdType *profile, std::vector<mcIdType>& code, std::vector<DataArrayIdType *>& idsInPflPerType, std::vector<DataArrayIdType *>& idsPerType, bool smartPflKiller=true) const;
    MEDCOUPLING_EXPORT MEDCouplingMesh *buildPart(const mcIdType *start, const mcIdType *end) const;
    MEDCOUPLING_EXPORT MEDCouplingMesh *buildPartAndReduceNodes(const mcIdType *start, const mcIdType *end, DataArrayIdType*& arr) const;
    MEDCOUPLING_EXPORT DataArrayIdType *simplexize(int policy);
    MEDCOUPLING_EXPORT MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const;
    MEDCOUPLING_EXPORT DataArrayDouble *getCoordinatesAndOwner() const;
    MEDCOUPLING_EXPORT DataArrayDouble *computeCellCenterOfMass() const;
    MEDCOUPLING_EXPORT DataArrayDouble *computeIsoBarycenterOfNodesPerCell() const;
    MEDCOUPLING_EXPORT void getReverseNodalConnectivity(DataArrayIdType *revNodal, DataArrayIdType *revNodalIndx) const;
    //Serialization unserialisation
    MEDCOUPLING_EXPORT void getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<mcIdType>& tinyInfo, std::vector<std::string>& littleStrings) const;
    MEDCOUPLING_EXPORT void resizeForUnserialization(const std::vector<mcIdType>& tinyInfo, DataArrayIdType *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const;
    MEDCOUPLING_EXPORT void serialize(DataArrayIdType *&a1, DataArrayDouble *&a2) const;
    MEDCOUPLING_EXPORT void unserialization(const std::vector<double>& tinyInfoD, const std::vector<mcIdType>& tinyInfo, const DataArrayIdType *a1, DataArrayDouble *a2,
                                            const std::vector<std::string>& littleStrings);
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
    MEDCOUPLING_EXPORT std::string getVTKFileExtension() const;
    mcIdType get2DCellIdForExtrusion() const { return _cell_2D_id; }
  private:
    MEDCouplingMappedExtrudedMesh(const MEDCouplingUMesh *mesh3D, const MEDCouplingUMesh *mesh2D, mcIdType cell2DId);
    MEDCouplingMappedExtrudedMesh(const MEDCouplingCMesh *mesh3D);
    MEDCouplingMappedExtrudedMesh(const MEDCouplingMappedExtrudedMesh& other, bool deepCpy);
    MEDCouplingMappedExtrudedMesh();
    void computeExtrusion(const MEDCouplingUMesh *mesh3D);
    void computeExtrusionAlg(const MEDCouplingUMesh *mesh3D);
    void build1DExtrusion(mcIdType idIn3DDesc, mcIdType newId, mcIdType nbOf1DLev, MEDCouplingUMesh *subMesh,
                          const mcIdType *desc3D, const mcIdType *descIndx3D,
                          const mcIdType *revDesc3D, const mcIdType *revDescIndx3D,
                          bool computeMesh1D);
    mcIdType findOppositeFaceOf(mcIdType current2DCell, mcIdType current3DCell, const std::vector<mcIdType>& connSorted,
                           const mcIdType *desc3D, const mcIdType *descIndx3D,
                           const mcIdType *conn2D, const mcIdType *conn2DIndx);
    void computeBaryCenterOfFace(const std::vector<mcIdType>& nodalConnec, mcIdType lev1DId);
    ~MEDCouplingMappedExtrudedMesh();
    void writeVTKLL(std::ostream& ofs, const std::string& cellData, const std::string& pointData, DataArrayByte *byteData) const;
    std::string getVTKDataSetType() const;
  private:
    MCAuto<MEDCouplingUMesh> _mesh2D;
    MCAuto<MEDCouplingUMesh> _mesh1D;
    //! New to old 3D cell Ids Array
    MCAuto<DataArrayIdType> _mesh3D_ids;
    mcIdType _cell_2D_id;
  };
}

#endif
