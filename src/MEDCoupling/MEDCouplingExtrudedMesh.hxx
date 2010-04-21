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
    MEDCouplingMeshType getType() const;
    bool isStructured() const;
    int getNumberOfCells() const;
    int getNumberOfNodes() const;
    int getSpaceDimension() const;
    int getMeshDimension() const;
    INTERP_KERNEL::NormalizedCellType getTypeOfCell(int cellId) const;
    void getNodeIdsOfCell(int cellId, std::vector<int>& conn) const;
    void getCoordinatesOfNode(int nodeId, std::vector<double>& coo) const;
    void checkCoherency() const throw (INTERP_KERNEL::Exception);
    void getBoundingBox(double *bbox) const;
    void updateTime();
    MEDCouplingUMesh *getMesh2D() const { return _mesh2D; }
    MEDCouplingUMesh *getMesh1D() const { return _mesh1D; }
    DataArrayInt *getMesh3DIds() const { return _mesh3D_ids; }
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
    MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const;
    DataArrayDouble *getCoordinatesAndOwner() const;
    DataArrayDouble *getBarycenterAndOwner() const;
  private:
    MEDCouplingExtrudedMesh(const MEDCouplingUMesh *mesh3D, MEDCouplingUMesh *mesh2D, int cell2DId) throw(INTERP_KERNEL::Exception);
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
    DataArrayInt *_mesh3D_ids;
    int _cell_2D_id;
  };
}

#endif
