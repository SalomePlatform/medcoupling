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

  class MEDCOUPLING_EXPORT MEDCouplingExtrudedMesh : public MEDCouplingMesh
  {
  public:
    static MEDCouplingExtrudedMesh *New(const MEDCouplingUMesh *mesh3D, MEDCouplingUMesh *mesh2D, int cell2DId) throw(INTERP_KERNEL::Exception);
    static MEDCouplingExtrudedMesh *New();
    MEDCouplingMeshType getType() const;
    void copyTinyStringsFrom(const MEDCouplingMesh *other) throw(INTERP_KERNEL::Exception);
    int getNumberOfCells() const;
    int getNumberOfNodes() const;
    int getSpaceDimension() const;
    int getMeshDimension() const;
    MEDCouplingMesh *deepCpy() const;
    MEDCouplingExtrudedMesh *clone(bool recDeepCpy) const;
    bool isEqual(const MEDCouplingMesh *other, double prec) const;
    bool isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const;
    void checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                              DataArrayInt *&cellCor, DataArrayInt *&nodeCor) const throw(INTERP_KERNEL::Exception);
    void checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                         DataArrayInt *&cellCor) const throw(INTERP_KERNEL::Exception);
    INTERP_KERNEL::NormalizedCellType getTypeOfCell(int cellId) const;
    int getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const;
    void getNodeIdsOfCell(int cellId, std::vector<int>& conn) const;
    void getCoordinatesOfNode(int nodeId, std::vector<double>& coo) const;
    std::string simpleRepr() const;
    std::string advancedRepr() const;
    void checkCoherency() const throw (INTERP_KERNEL::Exception);
    void getBoundingBox(double *bbox) const;
    void updateTime();
    void renumberCells(const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getMesh2D() const { return _mesh2D; }
    MEDCouplingUMesh *getMesh1D() const { return _mesh1D; }
    DataArrayInt *getMesh3DIds() const { return _mesh3D_ids; }
    MEDCouplingUMesh *build3DUnstructuredMesh() const;
    MEDCouplingFieldDouble *getMeasureField(bool) const;
    MEDCouplingFieldDouble *getMeasureFieldOnNode(bool) const;
    MEDCouplingFieldDouble *buildOrthogonalField() const;
    int getCellContainingPoint(const double *pos, double eps) const;
    static int findCorrespCellByNodalConn(const std::vector<int>& nodalConnec,
                                          const int *revNodalPtr, const int *revNodalIndxPtr) throw(INTERP_KERNEL::Exception);
    static void project1DMeshes(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2, double eps,
                                MEDCouplingUMesh *&m1r, MEDCouplingUMesh *&m2r, double *v) throw(INTERP_KERNEL::Exception);
    void rotate(const double *center, const double *vector, double angle);
    void translate(const double *vector);
    void scale(const double *point, double factor);
    MEDCouplingMesh *buildPart(const int *start, const int *end) const;
    MEDCouplingMesh *buildPartAndReduceNodes(const int *start, const int *end, DataArrayInt*& arr) const;
    MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const;
    DataArrayDouble *getCoordinatesAndOwner() const;
    DataArrayDouble *getBarycenterAndOwner() const;
    //Serialization unserialisation
    void getTinySerializationInformation(std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const;
    void resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const;
    void serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const;
    void unserialization(const std::vector<int>& tinyInfo, const DataArrayInt *a1, DataArrayDouble *a2,
                         const std::vector<std::string>& littleStrings);
  private:
    MEDCouplingExtrudedMesh(const MEDCouplingUMesh *mesh3D, MEDCouplingUMesh *mesh2D, int cell2DId) throw(INTERP_KERNEL::Exception);
    MEDCouplingExtrudedMesh(const MEDCouplingExtrudedMesh& other, bool deepCpy);
    MEDCouplingExtrudedMesh();
    void computeExtrusion(const MEDCouplingUMesh *mesh3D) throw(INTERP_KERNEL::Exception);
    void computeExtrusionAlg(const MEDCouplingUMesh *mesh3D) throw(INTERP_KERNEL::Exception);
    void build1DExtrusion(int idIn3DDesc, int newId, int nbOf1DLev, MEDCouplingUMesh *subMesh,
                          const int *desc3D, const int *descIndx3D,
                          const int *revDesc3D, const int *revDescIndx3D,
                          bool computeMesh1D) throw(INTERP_KERNEL::Exception);
    int findOppositeFaceOf(int current2DCell, int current3DCell, const std::vector<int>& connSorted,
                           const int *desc3D, const int *descIndx3D,
                           const int *conn2D, const int *conn2DIndx) throw(INTERP_KERNEL::Exception);
    void computeBaryCenterOfFace(const std::vector<int>& nodalConnec, int lev1DId);
    ~MEDCouplingExtrudedMesh();
  private:
    MEDCouplingUMesh *_mesh2D;
    MEDCouplingUMesh *_mesh1D;
    //! New to old 3D cell Ids Array
    DataArrayInt *_mesh3D_ids;
    int _cell_2D_id;
  };
}

#endif
