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

#ifndef __PARAMEDMEM_MEDCOUPLINGCMESH_HXX__
#define __PARAMEDMEM_MEDCOUPLINGCMESH_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingMesh.hxx"

namespace ParaMEDMEM
{
  class DataArrayDouble;

  class MEDCOUPLING_EXPORT MEDCouplingCMesh : public MEDCouplingMesh
  {
  public:
    static MEDCouplingCMesh *New();
    void updateTime();
    MEDCouplingMeshType getType() const { return CARTESIAN; }
    bool isEqual(const MEDCouplingMesh *other, double prec) const;
    void checkCoherency() const throw(INTERP_KERNEL::Exception);
    bool isStructured() const;
    int getNumberOfCells() const;
    int getNumberOfNodes() const;
    int getSpaceDimension() const;
    int getMeshDimension() const;
    int getCellIdFromPos(int i, int j, int k) const;
    int getNodeIdFromPos(int i, int j, int k) const;
    static void getPosFromId(int nodeId, int spaceDim, const int *split, int *res);
    INTERP_KERNEL::NormalizedCellType getTypeOfCell(int cellId) const;
    void getNodeIdsOfCell(int cellId, std::vector<int>& conn) const;
    void getCoordinatesOfNode(int nodeId, std::vector<double>& coo) const;
    DataArrayDouble *getCoordsAt(int i) const throw(INTERP_KERNEL::Exception);
    void setCoords(DataArrayDouble *coordsX,
                   DataArrayDouble *coordsY=0,
                   DataArrayDouble *coordsZ=0);
    // tools
    void getBoundingBox(double *bbox) const;
    MEDCouplingFieldDouble *getMeasureField(bool isAbs) const;
    MEDCouplingFieldDouble *getMeasureFieldOnNode(bool isAbs) const;
    MEDCouplingFieldDouble *buildOrthogonalField() const;
    int getCellContainingPoint(const double *pos, double eps) const;
    void rotate(const double *center, const double *vector, double angle);
    void translate(const double *vector);
    MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const;
    DataArrayDouble *getCoordinatesAndOwner() const;
    DataArrayDouble *getBarycenterAndOwner() const;
    //some useful methods
    void getSplitCellValues(int *res) const;
    void getSplitNodeValues(int *res) const;
  private:
    MEDCouplingCMesh();
    ~MEDCouplingCMesh();
  private:
    DataArrayDouble *_x_array;
    DataArrayDouble *_y_array;
    DataArrayDouble *_z_array;
  };
}

#endif
