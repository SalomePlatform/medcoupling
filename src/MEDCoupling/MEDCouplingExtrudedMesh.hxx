// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
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

namespace ParaMEDMEM
{
  class DataArrayInt;
  class DataArrayDouble;
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;

  class MEDCouplingExtrudedMesh : public MEDCouplingMesh
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingExtrudedMesh *New(const MEDCouplingUMesh *mesh3D, const MEDCouplingUMesh *mesh2D, int cell2DId);
    MEDCOUPLING_EXPORT static MEDCouplingExtrudedMesh *New();
    MEDCOUPLING_EXPORT MEDCouplingMeshType getType() const;
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDCOUPLING_EXPORT void copyTinyStringsFrom(const MEDCouplingMesh *other);
    MEDCOUPLING_EXPORT int getNumberOfCells() const;
    MEDCOUPLING_EXPORT int getNumberOfNodes() const;
    MEDCOUPLING_EXPORT int getSpaceDimension() const;
    MEDCOUPLING_EXPORT int getMeshDimension() const;
    MEDCOUPLING_EXPORT MEDCouplingExtrudedMesh *deepCpy() const;
    MEDCouplingExtrudedMesh *clone(bool recDeepCpy) const;
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const;
    MEDCOUPLING_EXPORT void checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                 DataArrayInt *&cellCor, DataArrayInt *&nodeCor) const;
    MEDCOUPLING_EXPORT void checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                            DataArrayInt *&cellCor) const;
    MEDCOUPLING_EXPORT INTERP_KERNEL::NormalizedCellType getTypeOfCell(int cellId) const;
    MEDCOUPLING_EXPORT std::set<INTERP_KERNEL::NormalizedCellType> getAllGeoTypes() const;
    MEDCOUPLING_EXPORT DataArrayInt *giveCellsWithType(INTERP_KERNEL::NormalizedCellType type) const;
    MEDCOUPLING_EXPORT DataArrayInt *computeNbOfNodesPerCell() const;
    MEDCOUPLING_EXPORT DataArrayInt *computeNbOfFacesPerCell() const;
    MEDCOUPLING_EXPORT DataArrayInt *computeEffectiveNbOfNodesPerCell() const;
    MEDCOUPLING_EXPORT int getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const;
    MEDCOUPLING_EXPORT void getNodeIdsOfCell(int cellId, std::vector<int>& conn) const;
    MEDCOUPLING_EXPORT void getCoordinatesOfNode(int nodeId, std::vector<double>& coo) const;
    MEDCOUPLING_EXPORT std::string simpleRepr() const;
    MEDCOUPLING_EXPORT std::string advancedRepr() const;
    MEDCOUPLING_EXPORT void checkCoherency() const;
    MEDCOUPLING_EXPORT void checkCoherency1(double eps=1e-12) const;
    MEDCOUPLING_EXPORT void getBoundingBox(double *bbox) const;
    MEDCOUPLING_EXPORT void updateTime() const;
    MEDCOUPLING_EXPORT void renumberCells(const int *old2NewBg, bool check=true);
    MEDCOUPLING_EXPORT MEDCouplingUMesh *getMesh2D() const { return _mesh2D; }
    MEDCOUPLING_EXPORT MEDCouplingUMesh *getMesh1D() const { return _mesh1D; }
    MEDCOUPLING_EXPORT DataArrayInt *getMesh3DIds() const { return _mesh3D_ids; }
    MEDCOUPLING_EXPORT MEDCouplingUMesh *build3DUnstructuredMesh() const;
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildUnstructured() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureField(bool) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureFieldOnNode(bool) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildOrthogonalField() const;
    MEDCOUPLING_EXPORT int getCellContainingPoint(const double *pos, double eps) const;
    MEDCOUPLING_EXPORT static int FindCorrespCellByNodalConn(const std::vector<int>& nodalConnec,
                                                             const int *revNodalPtr, const int *revNodalIndxPtr);
    MEDCOUPLING_EXPORT static void Project1DMeshes(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2, double eps,
                                                   MEDCouplingUMesh *&m1r, MEDCouplingUMesh *&m2r, double *v);
    MEDCOUPLING_EXPORT void rotate(const double *center, const double *vector, double angle);
    MEDCOUPLING_EXPORT void translate(const double *vector);
    MEDCOUPLING_EXPORT void scale(const double *point, double factor);
    MEDCOUPLING_EXPORT std::vector<int> getDistributionOfTypes() const;
    MEDCOUPLING_EXPORT DataArrayInt *checkTypeConsistencyAndContig(const std::vector<int>& code, const std::vector<const DataArrayInt *>& idsPerType) const;
    MEDCOUPLING_EXPORT void splitProfilePerType(const DataArrayInt *profile, std::vector<int>& code, std::vector<DataArrayInt *>& idsInPflPerType, std::vector<DataArrayInt *>& idsPerType) const;
    MEDCOUPLING_EXPORT MEDCouplingMesh *buildPart(const int *start, const int *end) const;
    MEDCOUPLING_EXPORT MEDCouplingMesh *buildPartAndReduceNodes(const int *start, const int *end, DataArrayInt*& arr) const;
    MEDCOUPLING_EXPORT DataArrayInt *simplexize(int policy);
    MEDCOUPLING_EXPORT MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const;
    MEDCOUPLING_EXPORT DataArrayDouble *getCoordinatesAndOwner() const;
    MEDCOUPLING_EXPORT DataArrayDouble *getBarycenterAndOwner() const;
    MEDCOUPLING_EXPORT DataArrayDouble *computeIsoBarycenterOfNodesPerCell() const;
    MEDCOUPLING_EXPORT void getReverseNodalConnectivity(DataArrayInt *revNodal, DataArrayInt *revNodalIndx) const;
    //Serialization unserialisation
    MEDCOUPLING_EXPORT void getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const;
    MEDCOUPLING_EXPORT void resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const;
    MEDCOUPLING_EXPORT void serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const;
    MEDCOUPLING_EXPORT void unserialization(const std::vector<double>& tinyInfoD, const std::vector<int>& tinyInfo, const DataArrayInt *a1, DataArrayDouble *a2,
                                            const std::vector<std::string>& littleStrings);
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
    MEDCOUPLING_EXPORT std::string getVTKFileExtension() const;
  private:
    MEDCouplingExtrudedMesh(const MEDCouplingUMesh *mesh3D, const MEDCouplingUMesh *mesh2D, int cell2DId);
    MEDCouplingExtrudedMesh(const MEDCouplingExtrudedMesh& other, bool deepCopy);
    MEDCouplingExtrudedMesh();
    void computeExtrusion(const MEDCouplingUMesh *mesh3D);
    void computeExtrusionAlg(const MEDCouplingUMesh *mesh3D);
    void build1DExtrusion(int idIn3DDesc, int newId, int nbOf1DLev, MEDCouplingUMesh *subMesh,
                          const int *desc3D, const int *descIndx3D,
                          const int *revDesc3D, const int *revDescIndx3D,
                          bool computeMesh1D);
    int findOppositeFaceOf(int current2DCell, int current3DCell, const std::vector<int>& connSorted,
                           const int *desc3D, const int *descIndx3D,
                           const int *conn2D, const int *conn2DIndx);
    void computeBaryCenterOfFace(const std::vector<int>& nodalConnec, int lev1DId);
    ~MEDCouplingExtrudedMesh();
    void writeVTKLL(std::ostream& ofs, const std::string& cellData, const std::string& pointData, DataArrayByte *byteData) const;
    std::string getVTKDataSetType() const;
  private:
    MEDCouplingUMesh *_mesh2D;
    MEDCouplingUMesh *_mesh1D;
    //! New to old 3D cell Ids Array
    DataArrayInt *_mesh3D_ids;
    int _cell_2D_id;
  };
}

#endif
