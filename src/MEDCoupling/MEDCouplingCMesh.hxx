// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
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

#ifndef __PARAMEDMEM_MEDCOUPLINGCMESH_HXX__
#define __PARAMEDMEM_MEDCOUPLINGCMESH_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingMesh.hxx"

namespace ParaMEDMEM
{
  class DataArrayDouble;
  class MEDCouplingUMesh;

  class MEDCOUPLING_EXPORT MEDCouplingCMesh : public MEDCouplingMesh
  {
  public:
    static MEDCouplingCMesh *New();
    static MEDCouplingCMesh *New(const char *meshName);
    MEDCouplingMesh *deepCpy() const;
    MEDCouplingCMesh *clone(bool recDeepCpy) const;
    void updateTime() const;
    MEDCouplingMeshType getType() const { return CARTESIAN; }
    void copyTinyStringsFrom(const MEDCouplingMesh *other) throw(INTERP_KERNEL::Exception);
    bool isEqual(const MEDCouplingMesh *other, double prec) const;
    bool isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const;
    void checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                              DataArrayInt *&cellCor, DataArrayInt *&nodeCor) const throw(INTERP_KERNEL::Exception);
    void checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                         DataArrayInt *&cellCor) const throw(INTERP_KERNEL::Exception);
    void checkCoherency() const throw(INTERP_KERNEL::Exception);
    void checkCoherency1(double eps=1e-12) const throw(INTERP_KERNEL::Exception);
    void checkCoherency2(double eps=1e-12) const throw(INTERP_KERNEL::Exception);
    int getNumberOfCells() const;
    int getNumberOfNodes() const;
    int getSpaceDimension() const;
    int getMeshDimension() const;
    int getCellIdFromPos(int i, int j, int k) const;
    int getNodeIdFromPos(int i, int j, int k) const;
    static void GetPosFromId(int nodeId, int spaceDim, const int *split, int *res);
    INTERP_KERNEL::NormalizedCellType getTypeOfCell(int cellId) const;
    std::set<INTERP_KERNEL::NormalizedCellType> getAllGeoTypes() const;
    int getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const;
    void getNodeIdsOfCell(int cellId, std::vector<int>& conn) const;
    void getCoordinatesOfNode(int nodeId, std::vector<double>& coo) const throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const;
    std::string advancedRepr() const;
    const DataArrayDouble *getCoordsAt(int i) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getCoordsAt(int i) throw(INTERP_KERNEL::Exception);
    void setCoordsAt(int i, const DataArrayDouble *arr) throw(INTERP_KERNEL::Exception);
    void setCoords(const DataArrayDouble *coordsX,
                   const DataArrayDouble *coordsY=0,
                   const DataArrayDouble *coordsZ=0);
    // tools
    std::vector<int> getDistributionOfTypes() const throw(INTERP_KERNEL::Exception);
    DataArrayInt *checkTypeConsistencyAndContig(const std::vector<int>& code, const std::vector<const DataArrayInt *>& idsPerType) const throw(INTERP_KERNEL::Exception);
    void splitProfilePerType(const DataArrayInt *profile, std::vector<int>& code, std::vector<DataArrayInt *>& idsInPflPerType, std::vector<DataArrayInt *>& idsPerType) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *buildUnstructured() const throw(INTERP_KERNEL::Exception);
    MEDCouplingMesh *buildPart(const int *start, const int *end) const;
    MEDCouplingMesh *buildPartAndReduceNodes(const int *start, const int *end, DataArrayInt*& arr) const;
    DataArrayInt *simplexize(int policy) throw(INTERP_KERNEL::Exception);
    void getBoundingBox(double *bbox) const;
    MEDCouplingFieldDouble *getMeasureField(bool isAbs) const;
    MEDCouplingFieldDouble *getMeasureFieldOnNode(bool isAbs) const;
    MEDCouplingFieldDouble *buildOrthogonalField() const;
    int getCellContainingPoint(const double *pos, double eps) const;
    void rotate(const double *center, const double *vector, double angle);
    void translate(const double *vector);
    void scale(const double *point, double factor);
    MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const;
    DataArrayDouble *getCoordinatesAndOwner() const;
    DataArrayDouble *getBarycenterAndOwner() const;
    void renumberCells(const int *old2NewBg, bool check=true) throw(INTERP_KERNEL::Exception);
    void fill1DUnstructuredMesh(MEDCouplingUMesh *m) const;
    void fill2DUnstructuredMesh(MEDCouplingUMesh *m) const;
    void fill3DUnstructuredMesh(MEDCouplingUMesh *m) const;
    //some useful methods
    void getSplitCellValues(int *res) const;
    void getSplitNodeValues(int *res) const;
    //serialisation-unserialization
    void getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const;
    void resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const;
    void serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const;
    void unserialization(const std::vector<double>& tinyInfoD, const std::vector<int>& tinyInfo, const DataArrayInt *a1, DataArrayDouble *a2,
                         const std::vector<std::string>& littleStrings);
  private:
    MEDCouplingCMesh();
    MEDCouplingCMesh(const MEDCouplingCMesh& other, bool deepCpy);
    ~MEDCouplingCMesh();
    void writeVTKLL(std::ostream& ofs, const std::string& cellData, const std::string& pointData) const throw(INTERP_KERNEL::Exception);
    std::string getVTKDataSetType() const throw(INTERP_KERNEL::Exception);
  private:
    DataArrayDouble *_x_array;
    DataArrayDouble *_y_array;
    DataArrayDouble *_z_array;
  };
}

#endif
