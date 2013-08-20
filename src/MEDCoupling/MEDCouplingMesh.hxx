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

#ifndef __PARAMEDMEM_MEDCOUPLINGMESH_HXX__
#define __PARAMEDMEM_MEDCOUPLINGMESH_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingTimeLabel.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"

#include "InterpKernelException.hxx"

#include <set>
#include <vector>

namespace ParaMEDMEM
{
  typedef enum
    {
      UNSTRUCTURED = 5,
      CARTESIAN = 7,
      EXTRUDED = 8,
      CURVE_LINEAR = 9,
      SINGLE_STATIC_GEO_TYPE_UNSTRUCTURED = 10,
      SINGLE_DYNAMIC_GEO_TYPE_UNSTRUCTURED = 11
    } MEDCouplingMeshType;
  // -- WARNING this enum must be synchronized with MEDCouplingCommon.i file ! --

  class DataArrayInt;
  class DataArrayByte;
  class DataArrayDouble;
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;

  class MEDCOUPLING_EXPORT MEDCouplingMesh : public RefCountObject, public TimeLabel
  {
  public:
    std::size_t getHeapMemorySize() const;
    void setName(const char *name) { _name=name; }
    std::string getName() const { return _name; }
    void setDescription(const char *descr) { _description=descr; }
    std::string getDescription() const { return _description; }
    double getTime(int& iteration, int& order) const { iteration=_iteration; order=_order; return _time; }
    void setTime(double val, int iteration, int order) { _time=val; _iteration=iteration; _order=order; }
    void setTimeUnit(const char *unit) { _time_unit=unit; }
    const char *getTimeUnit() const { return _time_unit.c_str(); }
    virtual MEDCouplingMesh *deepCpy() const = 0;
    virtual MEDCouplingMeshType getType() const = 0;
    bool isStructured() const;
    virtual void copyTinyStringsFrom(const MEDCouplingMesh *other) throw(INTERP_KERNEL::Exception);
    virtual void copyTinyInfoFrom(const MEDCouplingMesh *other) throw(INTERP_KERNEL::Exception);
    // comparison methods
    virtual bool isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const throw(INTERP_KERNEL::Exception);
    virtual bool isEqual(const MEDCouplingMesh *other, double prec) const throw(INTERP_KERNEL::Exception);
    virtual bool isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const = 0;
    virtual void checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                      DataArrayInt *&cellCor, DataArrayInt *&nodeCor) const throw(INTERP_KERNEL::Exception) = 0;
    virtual void checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                 DataArrayInt *&cellCor) const throw(INTERP_KERNEL::Exception) = 0;
    virtual void checkFastEquivalWith(const MEDCouplingMesh *other, double prec) const throw(INTERP_KERNEL::Exception);
    void checkGeoEquivalWith(const MEDCouplingMesh *other, int levOfCheck, double prec,
                             DataArrayInt *&cellCor, DataArrayInt *&nodeCor) const throw(INTERP_KERNEL::Exception);
    //
    virtual void checkCoherency() const throw(INTERP_KERNEL::Exception) = 0;
    virtual void checkCoherency1(double eps=1e-12) const throw(INTERP_KERNEL::Exception) = 0;
    virtual void checkCoherency2(double eps=1e-12) const throw(INTERP_KERNEL::Exception) = 0;
    virtual int getNumberOfCells() const = 0;
    virtual int getNumberOfNodes() const = 0;
    virtual int getSpaceDimension() const = 0;
    virtual int getMeshDimension() const = 0;
    virtual DataArrayDouble *getCoordinatesAndOwner() const = 0;
    virtual DataArrayDouble *getBarycenterAndOwner() const = 0;
    virtual DataArrayDouble *computeIsoBarycenterOfNodesPerCell() const throw(INTERP_KERNEL::Exception) = 0;
    virtual DataArrayInt *giveCellsWithType(INTERP_KERNEL::NormalizedCellType type) const throw(INTERP_KERNEL::Exception) = 0;
    virtual DataArrayInt *computeNbOfNodesPerCell() const throw(INTERP_KERNEL::Exception) = 0;
    virtual DataArrayInt *computeEffectiveNbOfNodesPerCell() const throw(INTERP_KERNEL::Exception) = 0;
    virtual DataArrayInt *computeNbOfFacesPerCell() const throw(INTERP_KERNEL::Exception) = 0;
    virtual int getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const = 0;
    virtual INTERP_KERNEL::NormalizedCellType getTypeOfCell(int cellId) const = 0;
    virtual std::set<INTERP_KERNEL::NormalizedCellType> getAllGeoTypes() const = 0;
    virtual void getNodeIdsOfCell(int cellId, std::vector<int>& conn) const = 0;
    virtual DataArrayInt *getCellIdsFullyIncludedInNodeIds(const int *partBg, const int *partEnd) const;
    virtual void getCoordinatesOfNode(int nodeId, std::vector<double>& coo) const throw(INTERP_KERNEL::Exception) = 0;
    virtual std::string simpleRepr() const = 0;
    virtual std::string advancedRepr() const = 0;
    // tools
    virtual std::vector<int> getDistributionOfTypes() const throw(INTERP_KERNEL::Exception) = 0;
    virtual DataArrayInt *checkTypeConsistencyAndContig(const std::vector<int>& code, const std::vector<const DataArrayInt *>& idsPerType) const throw(INTERP_KERNEL::Exception) = 0;
    virtual void splitProfilePerType(const DataArrayInt *profile, std::vector<int>& code, std::vector<DataArrayInt *>& idsInPflPerType, std::vector<DataArrayInt *>& idsPerType) const throw(INTERP_KERNEL::Exception) = 0;
    virtual void getBoundingBox(double *bbox) const = 0;
    virtual MEDCouplingFieldDouble *getMeasureField(bool isAbs) const = 0;
    virtual MEDCouplingFieldDouble *getMeasureFieldOnNode(bool isAbs) const = 0;
    virtual int getCellContainingPoint(const double *pos, double eps) const = 0;
    virtual void getCellsContainingPoint(const double *pos, double eps, std::vector<int>& elts) const;
    virtual void getCellsContainingPoints(const double *pos, int nbOfPoints, double eps, MEDCouplingAutoRefCountObjectPtr<DataArrayInt>& elts, MEDCouplingAutoRefCountObjectPtr<DataArrayInt>& eltsIndex) const;
    virtual MEDCouplingFieldDouble *fillFromAnalytic(TypeOfField t, int nbOfComp, FunctionToEvaluate func) const;
    virtual MEDCouplingFieldDouble *fillFromAnalytic(TypeOfField t, int nbOfComp, const char *func) const;
    virtual MEDCouplingFieldDouble *fillFromAnalytic2(TypeOfField t, int nbOfComp, const char *func) const;
    virtual MEDCouplingFieldDouble *fillFromAnalytic3(TypeOfField t, int nbOfComp, const std::vector<std::string>& varsOrder, const char *func) const;
    virtual MEDCouplingFieldDouble *buildOrthogonalField() const = 0;
    virtual void rotate(const double *center, const double *vector, double angle) = 0;
    virtual void translate(const double *vector) = 0;
    virtual void scale(const double *point, double factor) = 0;
    virtual void renumberCells(const int *old2NewBg, bool check=true) throw(INTERP_KERNEL::Exception) = 0;
    virtual MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const = 0;
    virtual MEDCouplingMesh *buildPart(const int *start, const int *end) const = 0;
    virtual MEDCouplingMesh *buildPartAndReduceNodes(const int *start, const int *end, DataArrayInt*& arr) const = 0;
    virtual MEDCouplingMesh *buildPartRange(int beginCellIds, int endCellIds, int stepCellIds) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingMesh *buildPartRangeAndReduceNodes(int beginCellIds, int endCellIds, int stepCellIds, int& beginOut, int& endOut, int& stepOut, DataArrayInt*& arr) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingUMesh *buildUnstructured() const throw(INTERP_KERNEL::Exception) = 0;
    virtual DataArrayInt *simplexize(int policy) throw(INTERP_KERNEL::Exception) = 0;
    virtual bool areCompatibleForMerge(const MEDCouplingMesh *other) const;
    static MEDCouplingMesh *MergeMeshes(const MEDCouplingMesh *mesh1, const MEDCouplingMesh *mesh2) throw(INTERP_KERNEL::Exception);
    static MEDCouplingMesh *MergeMeshes(std::vector<const MEDCouplingMesh *>& meshes) throw(INTERP_KERNEL::Exception);
    static bool IsStaticGeometricType(INTERP_KERNEL::NormalizedCellType type) throw(INTERP_KERNEL::Exception);
    static bool IsLinearGeometricType(INTERP_KERNEL::NormalizedCellType type) throw(INTERP_KERNEL::Exception);
    static INTERP_KERNEL::NormalizedCellType GetCorrespondingPolyType(INTERP_KERNEL::NormalizedCellType type) throw(INTERP_KERNEL::Exception);
    static int GetNumberOfNodesOfGeometricType(INTERP_KERNEL::NormalizedCellType type) throw(INTERP_KERNEL::Exception);
    static int GetDimensionOfGeometricType(INTERP_KERNEL::NormalizedCellType type) throw(INTERP_KERNEL::Exception);
    static const char *GetReprOfGeometricType(INTERP_KERNEL::NormalizedCellType type) throw(INTERP_KERNEL::Exception);
    //serialisation-unserialization
    virtual void getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const = 0;
    virtual void resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const = 0;
    virtual void serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const = 0;
    virtual void unserialization(const std::vector<double>& tinyInfoD, const std::vector<int>& tinyInfo, const DataArrayInt *a1, DataArrayDouble *a2,
                                 const std::vector<std::string>& littleStrings) = 0;
    void writeVTK(const char *fileName, bool isBinary=true) const throw(INTERP_KERNEL::Exception);
    /// @cond INTERNAL
    void writeVTKAdvanced(const char *fileName, const std::string& cda, const std::string& pda, DataArrayByte *byteData) const throw(INTERP_KERNEL::Exception);
    /// @endcond
    virtual void writeVTKLL(std::ostream& ofs, const std::string& cellData, const std::string& pointData, DataArrayByte *byteData) const throw(INTERP_KERNEL::Exception) = 0;
    virtual void reprQuickOverview(std::ostream& stream) const throw(INTERP_KERNEL::Exception) = 0;
  protected:
    MEDCouplingMesh();
    MEDCouplingMesh(const MEDCouplingMesh& other);
    virtual std::string getVTKDataSetType() const throw(INTERP_KERNEL::Exception) = 0;
    virtual ~MEDCouplingMesh() { }
  private:
    std::string _name;
    std::string _description;
    double _time;
    int _iteration;
    int _order;
    std::string _time_unit;
  };
}

#endif
