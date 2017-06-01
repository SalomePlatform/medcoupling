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

#ifndef __PARAMEDMEM_MEDCOUPLINGMESH_HXX__
#define __PARAMEDMEM_MEDCOUPLINGMESH_HXX__

#include "MEDCoupling.hxx"
#include "MCType.hxx"
#include "MEDCouplingTimeLabel.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "MCAuto.hxx"

#include "InterpKernelException.hxx"

#include <set>
#include <vector>

namespace MEDCoupling
{
  typedef enum
    {
      UNSTRUCTURED = 5,
      CARTESIAN = 7,
      EXTRUDED = 8,
      CURVE_LINEAR = 9,
      SINGLE_STATIC_GEO_TYPE_UNSTRUCTURED = 10,
      SINGLE_DYNAMIC_GEO_TYPE_UNSTRUCTURED = 11,
      IMAGE_GRID = 12
    } MEDCouplingMeshType;
  // -- WARNING this enum must be synchronized with MEDCouplingCommon.i file ! --

  class DataArrayInt;
  class DataArrayByte;
  class DataArrayDouble;
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;

  class MEDCouplingMesh : public RefCountObject, public TimeLabel
  {
  public:
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT void setName(const std::string& name) { _name=name; }
    MEDCOUPLING_EXPORT std::string getName() const { return _name; }
    MEDCOUPLING_EXPORT void setDescription(const std::string& descr) { _description=descr; }
    MEDCOUPLING_EXPORT std::string getDescription() const { return _description; }
    MEDCOUPLING_EXPORT double getTime(int& iteration, int& order) const { iteration=_iteration; order=_order; return _time; }
    MEDCOUPLING_EXPORT void setTime(double val, int iteration, int order) { _time=val; _iteration=iteration; _order=order; }
    MEDCOUPLING_EXPORT void setTimeUnit(const std::string& unit) { _time_unit=unit; }
    MEDCOUPLING_EXPORT std::string getTimeUnit() const { return _time_unit; }
    MEDCOUPLING_EXPORT virtual MEDCouplingMeshType getType() const = 0;
    MEDCOUPLING_EXPORT bool isStructured() const;
    // Copy methods
    MEDCOUPLING_EXPORT virtual MEDCouplingMesh *deepCopy() const = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingMesh *clone(bool recDeepCpy) const = 0;
    MEDCOUPLING_EXPORT virtual void copyTinyStringsFrom(const MEDCouplingMesh *other);
    MEDCOUPLING_EXPORT virtual void copyTinyInfoFrom(const MEDCouplingMesh *other);
    // comparison methods
    MEDCOUPLING_EXPORT virtual bool isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const;
    MEDCOUPLING_EXPORT virtual bool isEqual(const MEDCouplingMesh *other, double prec) const;
    MEDCOUPLING_EXPORT virtual bool isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const = 0;
    MEDCOUPLING_EXPORT virtual void checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                         DataArrayInt *&cellCor, DataArrayInt *&nodeCor) const = 0;
    MEDCOUPLING_EXPORT virtual void checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                                    DataArrayInt *&cellCor) const = 0;
    MEDCOUPLING_EXPORT virtual void checkFastEquivalWith(const MEDCouplingMesh *other, double prec) const;
    MEDCOUPLING_EXPORT void checkGeoEquivalWith(const MEDCouplingMesh *other, int levOfCheck, double prec,
                                                DataArrayInt *&cellCor, DataArrayInt *&nodeCor) const;
    //
    MEDCOUPLING_EXPORT virtual void checkConsistencyLight() const = 0;
    MEDCOUPLING_EXPORT virtual void checkConsistency(double eps=1e-12) const = 0;
    MEDCOUPLING_EXPORT virtual std::size_t getNumberOfCells() const = 0;
    MEDCOUPLING_EXPORT virtual int getNumberOfNodes() const = 0;
    MEDCOUPLING_EXPORT virtual int getSpaceDimension() const = 0;
    MEDCOUPLING_EXPORT virtual int getMeshDimension() const = 0;
    MEDCOUPLING_EXPORT virtual DataArrayDouble *getCoordinatesAndOwner() const = 0;
    MEDCOUPLING_EXPORT virtual DataArrayDouble *computeCellCenterOfMass() const = 0;
    MEDCOUPLING_EXPORT virtual DataArrayDouble *computeIsoBarycenterOfNodesPerCell() const = 0;
    MEDCOUPLING_EXPORT virtual DataArrayInt *giveCellsWithType(INTERP_KERNEL::NormalizedCellType type) const = 0;
    MEDCOUPLING_EXPORT virtual DataArrayInt *computeNbOfNodesPerCell() const = 0;
    MEDCOUPLING_EXPORT virtual DataArrayInt *computeEffectiveNbOfNodesPerCell() const = 0;
    MEDCOUPLING_EXPORT virtual DataArrayInt *computeNbOfFacesPerCell() const = 0;
    MEDCOUPLING_EXPORT virtual std::size_t getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const = 0;
    MEDCOUPLING_EXPORT virtual INTERP_KERNEL::NormalizedCellType getTypeOfCell(std::size_t cellId) const = 0;
    MEDCOUPLING_EXPORT virtual std::set<INTERP_KERNEL::NormalizedCellType> getAllGeoTypes() const = 0;
    MEDCOUPLING_EXPORT virtual void getNodeIdsOfCell(std::size_t cellId, std::vector<int>& conn) const = 0;
    MEDCOUPLING_EXPORT virtual DataArrayInt *getCellIdsFullyIncludedInNodeIds(const int *partBg, const int *partEnd) const;
    MEDCOUPLING_EXPORT virtual void getCoordinatesOfNode(int nodeId, std::vector<double>& coo) const = 0;
    MEDCOUPLING_EXPORT virtual std::string simpleRepr() const = 0;
    MEDCOUPLING_EXPORT virtual std::string advancedRepr() const = 0;
    // tools
    MEDCOUPLING_EXPORT virtual const DataArrayDouble *getDirectAccessOfCoordsArrIfInStructure() const = 0;
    MEDCOUPLING_EXPORT virtual std::vector<int> getDistributionOfTypes() const = 0;
    MEDCOUPLING_EXPORT virtual DataArrayInt *checkTypeConsistencyAndContig(const std::vector<int>& code, const std::vector<const DataArrayInt *>& idsPerType) const = 0;
    MEDCOUPLING_EXPORT virtual void splitProfilePerType(const DataArrayInt *profile, std::vector<int>& code, std::vector<DataArrayInt *>& idsInPflPerType, std::vector<DataArrayInt *>& idsPerType) const = 0;
    MEDCOUPLING_EXPORT virtual void getBoundingBox(double *bbox) const = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingFieldDouble *getMeasureField(bool isAbs) const = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingFieldDouble *getMeasureFieldOnNode(bool isAbs) const = 0;
    MEDCOUPLING_EXPORT virtual int getCellContainingPoint(const double *pos, double eps) const = 0;
    MEDCOUPLING_EXPORT virtual void getCellsContainingPoint(const double *pos, double eps, std::vector<int>& elts) const = 0;
    MEDCOUPLING_EXPORT virtual void getCellsContainingPoints(const double *pos, int nbOfPoints, double eps, MCAuto<DataArrayInt>& elts, MCAuto<DataArrayInt>& eltsIndex) const;
    MEDCOUPLING_EXPORT virtual MEDCouplingFieldDouble *fillFromAnalytic(TypeOfField t, int nbOfComp, FunctionToEvaluate func) const;
    MEDCOUPLING_EXPORT virtual MEDCouplingFieldDouble *fillFromAnalytic(TypeOfField t, int nbOfComp, const std::string& func) const;
    MEDCOUPLING_EXPORT virtual MEDCouplingFieldDouble *fillFromAnalyticCompo(TypeOfField t, int nbOfComp, const std::string& func) const;
    MEDCOUPLING_EXPORT virtual MEDCouplingFieldDouble *fillFromAnalyticNamedCompo(TypeOfField t, int nbOfComp, const std::vector<std::string>& varsOrder, const std::string& func) const;
    MEDCOUPLING_EXPORT virtual MEDCouplingFieldDouble *buildOrthogonalField() const = 0;
    MEDCOUPLING_EXPORT virtual void rotate(const double *center, const double *vector, double angle) = 0;
    MEDCOUPLING_EXPORT virtual void translate(const double *vector) = 0;
    MEDCOUPLING_EXPORT virtual void scale(const double *point, double factor) = 0;
    MEDCOUPLING_EXPORT virtual void renumberCells(const int *old2NewBg, bool check=true) = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingMesh *buildPart(const int *start, const int *end) const = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingMesh *buildPartAndReduceNodes(const int *start, const int *end, DataArrayInt*& arr) const = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingMesh *buildPartRange(int beginCellIds, int endCellIds, int stepCellIds) const;
    MEDCOUPLING_EXPORT virtual MEDCouplingMesh *buildPartRangeAndReduceNodes(int beginCellIds, int endCellIds, int stepCellIds, int& beginOut, int& endOut, int& stepOut, DataArrayInt*& arr) const;
    MEDCOUPLING_EXPORT virtual MEDCouplingUMesh *buildUnstructured() const = 0;
    MEDCOUPLING_EXPORT virtual DataArrayInt *simplexize(int policy) = 0;
    MEDCOUPLING_EXPORT virtual void getReverseNodalConnectivity(DataArrayInt *revNodal, DataArrayInt *revNodalIndx) const = 0;
    MEDCOUPLING_EXPORT virtual bool areCompatibleForMerge(const MEDCouplingMesh *other) const;
    MEDCOUPLING_EXPORT static MEDCouplingMesh *MergeMeshes(const MEDCouplingMesh *mesh1, const MEDCouplingMesh *mesh2);
    MEDCOUPLING_EXPORT static MEDCouplingMesh *MergeMeshes(std::vector<const MEDCouplingMesh *>& meshes);
    MEDCOUPLING_EXPORT static bool IsStaticGeometricType(INTERP_KERNEL::NormalizedCellType type);
    MEDCOUPLING_EXPORT static bool IsLinearGeometricType(INTERP_KERNEL::NormalizedCellType type);
    MEDCOUPLING_EXPORT static INTERP_KERNEL::NormalizedCellType GetCorrespondingPolyType(INTERP_KERNEL::NormalizedCellType type);
    MEDCOUPLING_EXPORT static int GetNumberOfNodesOfGeometricType(INTERP_KERNEL::NormalizedCellType type);
    MEDCOUPLING_EXPORT static int GetDimensionOfGeometricType(INTERP_KERNEL::NormalizedCellType type);
    MEDCOUPLING_EXPORT static const char *GetReprOfGeometricType(INTERP_KERNEL::NormalizedCellType type);
    //serialisation-unserialization
    MEDCOUPLING_EXPORT virtual void getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const = 0;
    MEDCOUPLING_EXPORT virtual void resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const = 0;
    MEDCOUPLING_EXPORT virtual void serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const = 0;
    MEDCOUPLING_EXPORT virtual void unserialization(const std::vector<double>& tinyInfoD, const std::vector<int>& tinyInfo, const DataArrayInt *a1, DataArrayDouble *a2,
                                                    const std::vector<std::string>& littleStrings) = 0;
    MEDCOUPLING_EXPORT std::string writeVTK(const std::string& fileName, bool isBinary=true) const;
    MEDCOUPLING_EXPORT std::string getVTKFileNameOf(const std::string& fileName) const;
    MEDCOUPLING_EXPORT virtual std::string getVTKFileExtension() const = 0;
    /// @cond INTERNAL
    MEDCOUPLING_EXPORT void writeVTKAdvanced(const std::string& fileName, const std::string& cda, const std::string& pda, DataArrayByte *byteData) const;
    MEDCOUPLING_EXPORT static void SplitExtension(const std::string& fileName, std::string& baseName, std::string& extension);
    /// @endcond
    MEDCOUPLING_EXPORT virtual void writeVTKLL(std::ostream& ofs, const std::string& cellData, const std::string& pointData, DataArrayByte *byteData) const = 0;
    MEDCOUPLING_EXPORT virtual void reprQuickOverview(std::ostream& stream) const = 0;
  protected:
    MEDCOUPLING_EXPORT MEDCouplingMesh();
    MEDCOUPLING_EXPORT MEDCouplingMesh(const MEDCouplingMesh& other);
    MEDCOUPLING_EXPORT virtual std::string getVTKDataSetType() const = 0;
    MEDCOUPLING_EXPORT virtual ~MEDCouplingMesh() { }
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
