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

#ifndef __PARAMEDMEM_MEDCOUPLINGCURVELINEARMESH_HXX__
#define __PARAMEDMEM_MEDCOUPLINGCURVELINEARMESH_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingStructuredMesh.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"

namespace ParaMEDMEM
{
  class DataArrayDouble;
  class MEDCouplingUMesh;

  class MEDCOUPLING_EXPORT MEDCouplingCurveLinearMesh : public MEDCouplingStructuredMesh
  {
  public:
    static MEDCouplingCurveLinearMesh *New();
    static MEDCouplingCurveLinearMesh *New(const char *meshName);
    MEDCouplingMesh *deepCpy() const;
    MEDCouplingCurveLinearMesh *clone(bool recDeepCpy) const;
    void updateTime() const;
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<RefCountObject *> getDirectChildren() const;
    MEDCouplingMeshType getType() const { return CURVE_LINEAR; }
    void copyTinyStringsFrom(const MEDCouplingMesh *other) throw(INTERP_KERNEL::Exception);
    bool isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const throw(INTERP_KERNEL::Exception);
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
    void getCoordinatesOfNode(int nodeId, std::vector<double>& coo) const throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const;
    std::string advancedRepr() const;
    DataArrayDouble *getCoords() throw(INTERP_KERNEL::Exception);
    const DataArrayDouble *getCoords() const throw(INTERP_KERNEL::Exception);
    void setCoords(const DataArrayDouble *coords) throw(INTERP_KERNEL::Exception);
    void setNodeGridStructure(const int *gridStructBg, const int *gridStructEnd) throw(INTERP_KERNEL::Exception);
    std::vector<int> getNodeGridStructure() const throw(INTERP_KERNEL::Exception);
    MEDCouplingStructuredMesh *buildStructuredSubPart(const std::vector< std::pair<int,int> >& cellPart) const throw(INTERP_KERNEL::Exception);
    // tools
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
    DataArrayDouble *computeIsoBarycenterOfNodesPerCell() const throw(INTERP_KERNEL::Exception);
    void renumberCells(const int *old2NewBg, bool check=true) throw(INTERP_KERNEL::Exception);
    //some useful methods
    void getSplitCellValues(int *res) const;
    void getSplitNodeValues(int *res) const;
    void getNodeGridStructure(int *res) const;
    //serialisation-unserialization
    void getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const;
    void resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const;
    void serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const;
    void unserialization(const std::vector<double>& tinyInfoD, const std::vector<int>& tinyInfo, const DataArrayInt *a1, DataArrayDouble *a2,
                         const std::vector<std::string>& littleStrings);
    void reprQuickOverview(std::ostream& stream) const throw(INTERP_KERNEL::Exception);
  private:
    void getMeasureFieldMeshDim1(bool isAbs, MEDCouplingFieldDouble *field) const throw(INTERP_KERNEL::Exception);
    void getMeasureFieldMeshDim2(bool isAbs, MEDCouplingFieldDouble *field) const throw(INTERP_KERNEL::Exception);
    void getMeasureFieldMeshDim3(bool isAbs, MEDCouplingFieldDouble *field) const throw(INTERP_KERNEL::Exception);
    void getBarycenterAndOwnerMeshDim3(DataArrayDouble *bary) const;
    void getBarycenterAndOwnerMeshDim2(DataArrayDouble *bary) const;
    void getBarycenterAndOwnerMeshDim1(DataArrayDouble *bary) const;
  private:
    MEDCouplingCurveLinearMesh();
    MEDCouplingCurveLinearMesh(const MEDCouplingCurveLinearMesh& other, bool deepCpy);
    ~MEDCouplingCurveLinearMesh();
    void writeVTKLL(std::ostream& ofs, const std::string& cellData, const std::string& pointData, DataArrayByte *byteData) const throw(INTERP_KERNEL::Exception);
    std::string getVTKDataSetType() const throw(INTERP_KERNEL::Exception);
  private:
    MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> _coords;
    std::vector<int> _structure;
  };
}

#endif
