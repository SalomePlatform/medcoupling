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

#ifndef __PARAMEDMEM_MEDCOUPLINGFIELDDISCRETIZATION_HXX__
#define __PARAMEDMEM_MEDCOUPLINGFIELDDISCRETIZATION_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "InterpKernelException.hxx"
#include "MEDCouplingTimeLabel.hxx"
#include "MEDCouplingNatureOfField.hxx"
#include "MEDCouplingGaussLocalization.hxx"
#include "MCAuto.hxx"

#include <set>
#include <vector>

namespace MEDCoupling
{
  class DataArray;
  class DataArrayIdType;
  class MEDCouplingMesh;
  class DataArrayDouble;
  class MEDCouplingFieldDouble;

  class MEDCouplingFieldDiscretization : public RefCountObject, public TimeLabel
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingFieldDiscretization *New(TypeOfField type);
    double getPrecision() const { return _precision; }
    void setPrecision(double val) { _precision=val; }
    MEDCOUPLING_EXPORT void updateTime() const;
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDCOUPLING_EXPORT static TypeOfField GetTypeOfFieldFromStringRepr(const std::string& repr);
    MEDCOUPLING_EXPORT static std::string GetTypeOfFieldRepr(TypeOfField type);
    MEDCOUPLING_EXPORT virtual TypeOfField getEnum() const = 0;
    MEDCOUPLING_EXPORT virtual bool isEqual(const MEDCouplingFieldDiscretization *other, double eps) const;
    MEDCOUPLING_EXPORT virtual bool isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const = 0;
    MEDCOUPLING_EXPORT virtual bool isEqualWithoutConsideringStr(const MEDCouplingFieldDiscretization *other, double eps) const;
    MEDCOUPLING_EXPORT virtual MEDCouplingFieldDiscretization *deepCopy() const;
    MEDCOUPLING_EXPORT virtual MEDCouplingFieldDiscretization *clone() const = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingFieldDiscretization *clonePart(const mcIdType *startCellIds, const mcIdType *endCellIds) const;
    MEDCOUPLING_EXPORT virtual MEDCouplingFieldDiscretization *clonePartRange(mcIdType beginCellIds, mcIdType endCellIds, mcIdType stepCellIds) const;
    MEDCOUPLING_EXPORT virtual std::string getStringRepr() const = 0;
    MEDCOUPLING_EXPORT virtual const char *getRepr() const = 0;
    MEDCOUPLING_EXPORT virtual mcIdType getNumberOfTuplesExpectedRegardingCode(const std::vector<mcIdType>& code, const std::vector<const DataArrayIdType *>& idsPerType) const = 0;
    MEDCOUPLING_EXPORT virtual mcIdType getNumberOfTuples(const MEDCouplingMesh *mesh) const = 0;
    MEDCOUPLING_EXPORT virtual mcIdType getNumberOfMeshPlaces(const MEDCouplingMesh *mesh) const = 0;
    MEDCOUPLING_EXPORT virtual DataArrayIdType *getOffsetArr(const MEDCouplingMesh *mesh) const = 0;
    MEDCOUPLING_EXPORT virtual void normL1(const MEDCouplingMesh *mesh, const DataArrayDouble *arr, double *res) const;
    MEDCOUPLING_EXPORT virtual void normL2(const MEDCouplingMesh *mesh, const DataArrayDouble *arr, double *res) const;
    MEDCOUPLING_EXPORT virtual void integral(const MEDCouplingMesh *mesh, const DataArrayDouble *arr, bool isWAbs, double *res) const;
    MEDCOUPLING_EXPORT virtual DataArrayDouble *getLocalizationOfDiscValues(const MEDCouplingMesh *mesh) const = 0;
    MEDCOUPLING_EXPORT virtual void computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, const mcIdType *tupleIdsBg, const mcIdType *tupleIdsEnd,
                                                                       DataArrayIdType *&cellRestriction, DataArrayIdType *&trueTupleRestriction) const = 0;
    MEDCOUPLING_EXPORT virtual void checkCompatibilityWithNature(NatureOfField nat) const = 0;
    MEDCOUPLING_EXPORT virtual void renumberCells(const mcIdType *old2NewBg, bool check=true);
    MEDCOUPLING_EXPORT virtual void renumberArraysForCell(const MEDCouplingMesh *mesh, const std::vector<DataArray *>& arrays,
                                                          const mcIdType *old2NewBg, bool check) = 0;
    MEDCOUPLING_EXPORT virtual double getIJK(const MEDCouplingMesh *mesh, const DataArrayDouble *da, mcIdType cellId, mcIdType nodeIdInCell, int compoId) const;
    MEDCOUPLING_EXPORT virtual void checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArray *da) const = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingFieldDouble *getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const = 0;
    MEDCOUPLING_EXPORT virtual void getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, double *res) const = 0;
    MEDCOUPLING_EXPORT virtual void getValueOnPos(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, mcIdType i, mcIdType j, mcIdType k, double *res) const = 0;
    MEDCOUPLING_EXPORT virtual DataArrayDouble *getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, mcIdType nbOfPoints) const = 0;
    MEDCOUPLING_EXPORT virtual DataArrayIdType *computeTupleIdsToSelectFromCellIds(const MEDCouplingMesh *mesh, const mcIdType *startCellIds, const mcIdType *endCellIds) const = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingMesh *buildSubMeshData(const MEDCouplingMesh *mesh, const mcIdType *start, const mcIdType *end, DataArrayIdType *&di) const = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingMesh *buildSubMeshDataRange(const MEDCouplingMesh *mesh, mcIdType beginCellIds, mcIdType endCellIds, mcIdType stepCellIds, mcIdType& beginOut, mcIdType& endOut, mcIdType& stepOut, DataArrayIdType *&di) const;
    MEDCOUPLING_EXPORT virtual void renumberValuesOnNodes(double epsOnVals, const mcIdType *old2New, mcIdType newNbOfNodes, DataArrayDouble *arr) const = 0;
    MEDCOUPLING_EXPORT virtual void renumberValuesOnCells(double epsOnVals, const MEDCouplingMesh *mesh, const mcIdType *old2New, mcIdType newSz, DataArrayDouble *arr) const = 0;
    MEDCOUPLING_EXPORT virtual void renumberValuesOnCellsR(const MEDCouplingMesh *mesh, const mcIdType *new2old, mcIdType newSz, DataArrayDouble *arr) const = 0;
    MEDCOUPLING_EXPORT virtual MCAuto<MEDCouplingFieldDiscretization> aggregate(std::vector<const MEDCouplingFieldDiscretization *>& fds) const = 0;
    MEDCOUPLING_EXPORT virtual void getSerializationIntArray(DataArrayIdType *& arr) const;
    MEDCOUPLING_EXPORT virtual void getTinySerializationIntInformation(std::vector<mcIdType>& tinyInfo) const;
    MEDCOUPLING_EXPORT virtual void getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const;
    MEDCOUPLING_EXPORT virtual void finishUnserialization(const std::vector<double>& tinyInfo);
    MEDCOUPLING_EXPORT virtual void resizeForUnserialization(const std::vector<mcIdType>& tinyInfo, DataArrayIdType *& arr);
    MEDCOUPLING_EXPORT virtual void checkForUnserialization(const std::vector<mcIdType>& tinyInfo, const DataArrayIdType *arr);
    MEDCOUPLING_EXPORT virtual void setGaussLocalizationOnType(const MEDCouplingMesh *m, INTERP_KERNEL::NormalizedCellType type, const std::vector<double>& refCoo,
                                                               const std::vector<double>& gsCoo, const std::vector<double>& wg);
    MEDCOUPLING_EXPORT virtual void setGaussLocalizationOnCells(const MEDCouplingMesh *m, const mcIdType *begin, const mcIdType *end, const std::vector<double>& refCoo,
                                                                const std::vector<double>& gsCoo, const std::vector<double>& wg);
    MEDCOUPLING_EXPORT virtual void clearGaussLocalizations();
    MEDCOUPLING_EXPORT virtual MEDCouplingGaussLocalization& getGaussLocalization(mcIdType locId);
    MEDCOUPLING_EXPORT virtual mcIdType getNbOfGaussLocalization() const;
    MEDCOUPLING_EXPORT virtual mcIdType getGaussLocalizationIdOfOneCell(mcIdType cellId) const;
    MEDCOUPLING_EXPORT virtual mcIdType getGaussLocalizationIdOfOneType(INTERP_KERNEL::NormalizedCellType type) const;
    MEDCOUPLING_EXPORT virtual std::set<mcIdType> getGaussLocalizationIdsOfOneType(INTERP_KERNEL::NormalizedCellType type) const;
    MEDCOUPLING_EXPORT virtual void getCellIdsHavingGaussLocalization(mcIdType locId, std::vector<mcIdType>& cellIds) const;
    MEDCOUPLING_EXPORT virtual const MEDCouplingGaussLocalization& getGaussLocalization(mcIdType locId) const;
    MEDCOUPLING_EXPORT virtual void reprQuickOverview(std::ostream& stream) const = 0;
    MEDCOUPLING_EXPORT virtual ~MEDCouplingFieldDiscretization();
  protected:
    MEDCOUPLING_EXPORT MEDCouplingFieldDiscretization();
    MEDCOUPLING_EXPORT static void RenumberEntitiesFromO2NArr(double epsOnVals, const mcIdType *old2NewPtr, mcIdType newNbOfEntity, DataArrayDouble *arr, const std::string& msg);
    MEDCOUPLING_EXPORT static void RenumberEntitiesFromN2OArr(const mcIdType *new2OldPtr, mcIdType new2OldSz, DataArrayDouble *arr, const std::string& msg);
    template<class FIELD_DISC>
    static MCAuto<MEDCouplingFieldDiscretization> EasyAggregate(std::vector<const MEDCouplingFieldDiscretization *>& fds);
  protected:
    double _precision;
    static const double DFLT_PRECISION;
  };

  class MEDCouplingFieldDiscretizationP0 : public MEDCouplingFieldDiscretization
  {
  public:
    MEDCOUPLING_EXPORT TypeOfField getEnum() const;
    std::string getClassName() const override { return std::string("MEDCouplingFieldDiscretizationP0"); }
    MEDCOUPLING_EXPORT MEDCouplingFieldDiscretization *clone() const override;
    MEDCOUPLING_EXPORT std::string getStringRepr() const;
    MEDCOUPLING_EXPORT const char *getRepr() const;
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const override;
    MEDCOUPLING_EXPORT mcIdType getNumberOfTuplesExpectedRegardingCode(const std::vector<mcIdType>& code, const std::vector<const DataArrayIdType *>& idsPerType) const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfTuples(const MEDCouplingMesh *mesh) const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfMeshPlaces(const MEDCouplingMesh *mesh) const;
    MEDCOUPLING_EXPORT DataArrayIdType *getOffsetArr(const MEDCouplingMesh *mesh) const;
    MEDCOUPLING_EXPORT void renumberArraysForCell(const MEDCouplingMesh *mesh, const std::vector<DataArray *>& arrays,
                                                  const mcIdType *old2NewBg, bool check);
    MEDCOUPLING_EXPORT DataArrayDouble *getLocalizationOfDiscValues(const MEDCouplingMesh *mesh) const;
    MEDCOUPLING_EXPORT void checkCompatibilityWithNature(NatureOfField nat) const override;
    MEDCOUPLING_EXPORT void computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, const mcIdType *tupleIdsBg, const mcIdType *tupleIdsEnd,
                                                               DataArrayIdType *&cellRestriction, DataArrayIdType *&trueTupleRestriction) const;
    MEDCOUPLING_EXPORT void checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArray *da) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const override;
    MEDCOUPLING_EXPORT void getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, double *res) const override;
    MEDCOUPLING_EXPORT void getValueOnPos(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, mcIdType i, mcIdType j, mcIdType k, double *res) const;
    MEDCOUPLING_EXPORT DataArrayDouble *getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, mcIdType nbOfPoints) const override;
    MEDCOUPLING_EXPORT void renumberValuesOnNodes(double epsOnVals, const mcIdType *old2New, mcIdType newNbOfNodes, DataArrayDouble *arr) const;
    MEDCOUPLING_EXPORT void renumberValuesOnCells(double epsOnVals, const MEDCouplingMesh *mesh, const mcIdType *old2New, mcIdType newSz, DataArrayDouble *arr) const;
    MEDCOUPLING_EXPORT void renumberValuesOnCellsR(const MEDCouplingMesh *mesh, const mcIdType *new2old, mcIdType newSz, DataArrayDouble *arr) const;
    MEDCOUPLING_EXPORT MCAuto<MEDCouplingFieldDiscretization> aggregate(std::vector<const MEDCouplingFieldDiscretization *>& fds) const override;
    MEDCOUPLING_EXPORT MEDCouplingMesh *buildSubMeshData(const MEDCouplingMesh *mesh, const mcIdType *start, const mcIdType *end, DataArrayIdType *&di) const;
    MEDCOUPLING_EXPORT MEDCouplingMesh *buildSubMeshDataRange(const MEDCouplingMesh *mesh, mcIdType beginCellIds, mcIdType endCellIds, mcIdType stepCellIds, mcIdType& beginOut, mcIdType& endOut, mcIdType& stepOut, DataArrayIdType *&di) const;
    MEDCOUPLING_EXPORT DataArrayIdType *computeTupleIdsToSelectFromCellIds(const MEDCouplingMesh *mesh, const mcIdType *startCellIds, const mcIdType *endCellIds) const;
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
  public:
    static const char REPR[];
    static constexpr TypeOfField TYPE = ON_CELLS;
  };

  class MEDCouplingFieldDiscretizationOnNodes : public MEDCouplingFieldDiscretization
  {
  public:
    MEDCOUPLING_EXPORT mcIdType getNumberOfTuples(const MEDCouplingMesh *mesh) const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfTuplesExpectedRegardingCode(const std::vector<mcIdType>& code, const std::vector<const DataArrayIdType *>& idsPerType) const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfMeshPlaces(const MEDCouplingMesh *mesh) const;
    MEDCOUPLING_EXPORT DataArrayIdType *getOffsetArr(const MEDCouplingMesh *mesh) const;
    MEDCOUPLING_EXPORT void renumberArraysForCell(const MEDCouplingMesh *mesh, const std::vector<DataArray *>& arrays,
                                                  const mcIdType *old2NewBg, bool check);
    MEDCOUPLING_EXPORT DataArrayDouble *getLocalizationOfDiscValues(const MEDCouplingMesh *mesh) const;
    MEDCOUPLING_EXPORT void computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, const mcIdType *tupleIdsBg, const mcIdType *tupleIdsEnd,
                                                               DataArrayIdType *&cellRestriction, DataArrayIdType *&trueTupleRestriction) const;
    MEDCOUPLING_EXPORT void checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArray *da) const;
    MEDCOUPLING_EXPORT MEDCouplingMesh *buildSubMeshData(const MEDCouplingMesh *mesh, const mcIdType *start, const mcIdType *end, DataArrayIdType *&di) const;
    MEDCOUPLING_EXPORT MEDCouplingMesh *buildSubMeshDataRange(const MEDCouplingMesh *mesh, mcIdType beginCellIds, mcIdType endCellIds, mcIdType stepCellIds, mcIdType& beginOut, mcIdType& endOut, mcIdType& stepOut, DataArrayIdType *&di) const;
    MEDCOUPLING_EXPORT DataArrayIdType *computeTupleIdsToSelectFromCellIds(const MEDCouplingMesh *mesh, const mcIdType *startCellIds, const mcIdType *endCellIds) const;
    MEDCOUPLING_EXPORT void renumberValuesOnNodes(double epsOnVals, const mcIdType *old2New, mcIdType newNbOfNodes, DataArrayDouble *arr) const;
    MEDCOUPLING_EXPORT void renumberValuesOnCells(double epsOnVals, const MEDCouplingMesh *mesh, const mcIdType *old2New, mcIdType newSz, DataArrayDouble *arr) const;
    MEDCOUPLING_EXPORT void renumberValuesOnCellsR(const MEDCouplingMesh *mesh, const mcIdType *new2old, mcIdType newSz, DataArrayDouble *arr) const;
  public:
    MEDCOUPLING_EXPORT void getValueOnPos(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, mcIdType i, mcIdType j, mcIdType k, double *res) const;
  };

  class MEDCouplingFieldDiscretizationP1 : public MEDCouplingFieldDiscretizationOnNodes
  {
  public:
    MEDCOUPLING_EXPORT TypeOfField getEnum() const;
    std::string getClassName() const override { return std::string("MEDCouplingFieldDiscretizationP1"); }
    MEDCOUPLING_EXPORT MEDCouplingFieldDiscretization *clone() const override;
    MEDCOUPLING_EXPORT std::string getStringRepr() const;
    MEDCOUPLING_EXPORT const char *getRepr() const;
    MEDCOUPLING_EXPORT void checkCompatibilityWithNature(NatureOfField nat) const override;
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const override;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const override;
    MEDCOUPLING_EXPORT void getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, double *res) const override;
    MEDCOUPLING_EXPORT DataArrayDouble *getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, mcIdType nbOfPoints) const override;
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
    MEDCOUPLING_EXPORT MCAuto<MEDCouplingFieldDiscretization> aggregate(std::vector<const MEDCouplingFieldDiscretization *>& fds) const override;
  public:
    static const char REPR[];
    static constexpr TypeOfField TYPE = ON_NODES;
  protected:
    MEDCOUPLING_EXPORT void getValueInCell(const MEDCouplingMesh *mesh, mcIdType cellId, const DataArrayDouble *arr, const double *loc, double *res) const;
  };

  /*!
   * This class abstracts MEDCouplingFieldDiscretization that needs an information on each cell to perform their job.
   * All classes that inherits from this are more linked to mesh.
   */
  class MEDCouplingFieldDiscretizationPerCell : public MEDCouplingFieldDiscretization
  {
  public:
    MEDCOUPLING_EXPORT const DataArrayIdType *getArrayOfDiscIds() const;
    MEDCOUPLING_EXPORT void setArrayOfDiscIds(const DataArrayIdType *adids);
    MEDCOUPLING_EXPORT void checkNoOrphanCells() const;
    MEDCOUPLING_EXPORT std::vector<DataArrayIdType *> splitIntoSingleGaussDicrPerCellType(std::vector< mcIdType >& locIds) const;
  protected:
    MEDCouplingFieldDiscretizationPerCell();
    MEDCouplingFieldDiscretizationPerCell(const MEDCouplingFieldDiscretizationPerCell& other, const mcIdType *startCellIds, const mcIdType *endCellIds);
    MEDCouplingFieldDiscretizationPerCell(const MEDCouplingFieldDiscretizationPerCell& other, mcIdType beginCellIds, mcIdType endCellIds, mcIdType stepCellIds);
    MEDCouplingFieldDiscretizationPerCell(DataArrayIdType *dpc);
    ~MEDCouplingFieldDiscretizationPerCell();
    void updateTime() const;
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    void checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArray *da) const;
    bool isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const override;
    bool isEqualWithoutConsideringStr(const MEDCouplingFieldDiscretization *other, double eps) const;
    void renumberCells(const mcIdType *old2NewBg, bool check);
  protected:
    void buildDiscrPerCellIfNecessary(const MEDCouplingMesh *mesh);
  protected:
    DataArrayIdType *_discr_per_cell;
    static const mcIdType DFT_INVALID_LOCID_VALUE;
  };

  class MEDCouplingFieldDiscretizationGauss : public MEDCouplingFieldDiscretizationPerCell
  {
  public:
    MEDCOUPLING_EXPORT MEDCouplingFieldDiscretizationGauss();
    MEDCOUPLING_EXPORT TypeOfField getEnum() const;
    std::string getClassName() const override { return std::string("MEDCouplingFieldDiscretizationGauss"); }
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const override;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const MEDCouplingFieldDiscretization *other, double eps) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDiscretization *clone() const override;
    MEDCOUPLING_EXPORT MEDCouplingFieldDiscretization *clonePart(const mcIdType *startCellIds, const mcIdType *endCellIds) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDiscretization *clonePartRange(mcIdType beginCellIds, mcIdType endCellIds, mcIdType stepCellIds) const;
    MEDCOUPLING_EXPORT std::string getStringRepr() const;
    MEDCOUPLING_EXPORT const char *getRepr() const;
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfTuplesExpectedRegardingCode(const std::vector<mcIdType>& code, const std::vector<const DataArrayIdType *>& idsPerType) const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfTuples(const MEDCouplingMesh *mesh) const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfMeshPlaces(const MEDCouplingMesh *mesh) const;
    MEDCOUPLING_EXPORT DataArrayIdType *getOffsetArr(const MEDCouplingMesh *mesh) const;
    MEDCOUPLING_EXPORT void renumberArraysForCell(const MEDCouplingMesh *mesh, const std::vector<DataArray *>& arrays,
                                                  const mcIdType *old2NewBg, bool check);
    MEDCOUPLING_EXPORT DataArrayDouble *getLocalizationOfDiscValues(const MEDCouplingMesh *mesh) const;
    MEDCOUPLING_EXPORT void computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, const mcIdType *tupleIdsBg, const mcIdType *tupleIdsEnd,
                                                               DataArrayIdType *&cellRestriction, DataArrayIdType *&trueTupleRestriction) const;
    MEDCOUPLING_EXPORT void checkCompatibilityWithNature(NatureOfField nat) const override;
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation(std::vector<mcIdType>& tinyInfo) const;
    MEDCOUPLING_EXPORT void getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const;
    MEDCOUPLING_EXPORT void finishUnserialization(const std::vector<double>& tinyInfo);
    MEDCOUPLING_EXPORT void getSerializationIntArray(DataArrayIdType *& arr) const;
    MEDCOUPLING_EXPORT void resizeForUnserialization(const std::vector<mcIdType>& tinyInfo, DataArrayIdType *& arr);
    MEDCOUPLING_EXPORT void checkForUnserialization(const std::vector<mcIdType>& tinyInfo, const DataArrayIdType *arr);
    MEDCOUPLING_EXPORT double getIJK(const MEDCouplingMesh *mesh, const DataArrayDouble *da, mcIdType cellId, mcIdType nodeIdInCell, int compoId) const;
    MEDCOUPLING_EXPORT void checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArray *da) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const override;
    MEDCOUPLING_EXPORT void getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, double *res) const override;
    MEDCOUPLING_EXPORT void getValueOnPos(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, mcIdType i, mcIdType j, mcIdType k, double *res) const;
    MEDCOUPLING_EXPORT DataArrayDouble *getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, mcIdType nbOfPoints) const override;
    MEDCOUPLING_EXPORT MEDCouplingMesh *buildSubMeshData(const MEDCouplingMesh *mesh, const mcIdType *start, const mcIdType *end, DataArrayIdType *&di) const;
    MEDCOUPLING_EXPORT MEDCouplingMesh *buildSubMeshDataRange(const MEDCouplingMesh *mesh, mcIdType beginCellIds, mcIdType endCellIds, mcIdType stepCellIds, mcIdType& beginOut, mcIdType& endOut, mcIdType& stepOut, DataArrayIdType *&di) const;
    MEDCOUPLING_EXPORT DataArrayIdType *computeTupleIdsToSelectFromCellIds(const MEDCouplingMesh *mesh, const mcIdType *startCellIds, const mcIdType *endCellIds) const;
    MEDCOUPLING_EXPORT void renumberValuesOnNodes(double epsOnVals, const mcIdType *old2New, mcIdType newNbOfNodes, DataArrayDouble *arr) const;
    MEDCOUPLING_EXPORT void renumberValuesOnCells(double epsOnVals, const MEDCouplingMesh *mesh, const mcIdType *old2New, mcIdType newSz, DataArrayDouble *arr) const;
    MEDCOUPLING_EXPORT void renumberValuesOnCellsR(const MEDCouplingMesh *mesh, const mcIdType *new2old, mcIdType newSz, DataArrayDouble *arr) const;
    MEDCOUPLING_EXPORT MCAuto<MEDCouplingFieldDiscretization> aggregate(std::vector<const MEDCouplingFieldDiscretization *>& fds) const override;
    MEDCOUPLING_EXPORT void setGaussLocalizationOnType(const MEDCouplingMesh *mesh, INTERP_KERNEL::NormalizedCellType type, const std::vector<double>& refCoo,
                                                       const std::vector<double>& gsCoo, const std::vector<double>& wg);
    MEDCOUPLING_EXPORT void setGaussLocalizationOnCells(const MEDCouplingMesh *mesh, const mcIdType *begin, const mcIdType *end, const std::vector<double>& refCoo,
                                                        const std::vector<double>& gsCoo, const std::vector<double>& wg);
    MEDCOUPLING_EXPORT void clearGaussLocalizations();
    MEDCOUPLING_EXPORT void setGaussLocalization(mcIdType locId, const MEDCouplingGaussLocalization& loc);
    MEDCOUPLING_EXPORT void resizeLocalizationVector(mcIdType newSz);
    MEDCOUPLING_EXPORT MEDCouplingGaussLocalization& getGaussLocalization(mcIdType locId);
    MEDCOUPLING_EXPORT mcIdType getNbOfGaussLocalization() const;
    MEDCOUPLING_EXPORT mcIdType getGaussLocalizationIdOfOneCell(mcIdType cellId) const;
    MEDCOUPLING_EXPORT mcIdType getGaussLocalizationIdOfOneType(INTERP_KERNEL::NormalizedCellType type) const;
    MEDCOUPLING_EXPORT std::set<mcIdType> getGaussLocalizationIdsOfOneType(INTERP_KERNEL::NormalizedCellType type) const;
    MEDCOUPLING_EXPORT void getCellIdsHavingGaussLocalization(mcIdType locId, std::vector<mcIdType>& cellIds) const;
    MEDCOUPLING_EXPORT const MEDCouplingGaussLocalization& getGaussLocalization(mcIdType locId) const;
    MEDCOUPLING_EXPORT DataArrayIdType *buildNbOfGaussPointPerCellField() const;
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
  protected:
    MEDCouplingFieldDiscretizationGauss(const MEDCouplingFieldDiscretizationGauss& other, const mcIdType *startCellIds=0, const mcIdType *endCellIds=0);
    MEDCouplingFieldDiscretizationGauss(const MEDCouplingFieldDiscretizationGauss& other, mcIdType beginCellIds, mcIdType endCellIds, mcIdType stepCellIds);
    MEDCouplingFieldDiscretizationGauss(DataArrayIdType *dpc, const std::vector<MEDCouplingGaussLocalization>& loc):MEDCouplingFieldDiscretizationPerCell(dpc),_loc(loc) { }
    void zipGaussLocalizations();
    mcIdType getOffsetOfCell(mcIdType cellId) const;
    void checkLocalizationId(mcIdType locId) const;
    void commonUnserialization(const std::vector<mcIdType>& tinyInfo);
  public:
    static const char REPR[];
    static constexpr TypeOfField TYPE = ON_GAUSS_PT;
  private:
    std::vector<MEDCouplingGaussLocalization> _loc;
  };

  /*!
   * Gauss with points of values located on nodes of element. This is a specialization of MEDCouplingFieldDiscretizationGauss.
   */
  class MEDCouplingFieldDiscretizationGaussNE : public MEDCouplingFieldDiscretization
  {
  public:
    MEDCOUPLING_EXPORT MEDCouplingFieldDiscretizationGaussNE();
    MEDCOUPLING_EXPORT TypeOfField getEnum() const;
    std::string getClassName() const override { return std::string("MEDCouplingFieldDiscretizationGaussNE"); }
    MEDCOUPLING_EXPORT MEDCouplingFieldDiscretization *clone() const override;
    MEDCOUPLING_EXPORT std::string getStringRepr() const;
    MEDCOUPLING_EXPORT const char *getRepr() const;
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const override;
    MEDCOUPLING_EXPORT mcIdType getNumberOfTuplesExpectedRegardingCode(const std::vector<mcIdType>& code, const std::vector<const DataArrayIdType *>& idsPerType) const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfTuples(const MEDCouplingMesh *mesh) const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfMeshPlaces(const MEDCouplingMesh *mesh) const;
    MEDCOUPLING_EXPORT DataArrayIdType *getOffsetArr(const MEDCouplingMesh *mesh) const;
    MEDCOUPLING_EXPORT void renumberArraysForCell(const MEDCouplingMesh *mesh, const std::vector<DataArray *>& arrays,
                                                  const mcIdType *old2NewBg, bool check);
    MEDCOUPLING_EXPORT DataArrayDouble *getLocalizationOfDiscValues(const MEDCouplingMesh *mesh) const;
    MEDCOUPLING_EXPORT void integral(const MEDCouplingMesh *mesh, const DataArrayDouble *arr, bool isWAbs, double *res) const;
    MEDCOUPLING_EXPORT void computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, const mcIdType *tupleIdsBg, const mcIdType *tupleIdsEnd,
                                                               DataArrayIdType *&cellRestriction, DataArrayIdType *&trueTupleRestriction) const;
    MEDCOUPLING_EXPORT void checkCompatibilityWithNature(NatureOfField nat) const override;
    MEDCOUPLING_EXPORT double getIJK(const MEDCouplingMesh *mesh, const DataArrayDouble *da, mcIdType cellId, mcIdType nodeIdInCell, int compoId) const;
    MEDCOUPLING_EXPORT void checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArray *da) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const override;
    MEDCOUPLING_EXPORT void getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, double *res) const override;
    MEDCOUPLING_EXPORT void getValueOnPos(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, mcIdType i, mcIdType j, mcIdType k, double *res) const;
    MEDCOUPLING_EXPORT DataArrayDouble *getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, mcIdType nbOfPoints) const override;
    MEDCOUPLING_EXPORT MEDCouplingMesh *buildSubMeshData(const MEDCouplingMesh *mesh, const mcIdType *start, const mcIdType *end, DataArrayIdType *&di) const;
    MEDCOUPLING_EXPORT MEDCouplingMesh *buildSubMeshDataRange(const MEDCouplingMesh *mesh, mcIdType beginCellIds, mcIdType endCellIds, mcIdType stepCellIds, mcIdType& beginOut, mcIdType& endOut, mcIdType& stepOut, DataArrayIdType *&di) const;
    MEDCOUPLING_EXPORT DataArrayIdType *computeTupleIdsToSelectFromCellIds(const MEDCouplingMesh *mesh, const mcIdType *startCellIds, const mcIdType *endCellIds) const;
    MEDCOUPLING_EXPORT void renumberValuesOnNodes(double epsOnVals, const mcIdType *old2New, mcIdType newNbOfNodes, DataArrayDouble *arr) const;
    MEDCOUPLING_EXPORT void renumberValuesOnCells(double epsOnVals, const MEDCouplingMesh *mesh, const mcIdType *old2New, mcIdType newSz, DataArrayDouble *arr) const;
    MEDCOUPLING_EXPORT void renumberValuesOnCellsR(const MEDCouplingMesh *mesh, const mcIdType *new2old, mcIdType newSz, DataArrayDouble *arr) const;
    MEDCOUPLING_EXPORT MCAuto<MEDCouplingFieldDiscretization> aggregate(std::vector<const MEDCouplingFieldDiscretization *>& fds) const override;
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
    MEDCOUPLING_EXPORT static const double *GetWeightArrayFromGeometricType(INTERP_KERNEL::NormalizedCellType geoType, std::size_t& lgth);
    MEDCOUPLING_EXPORT static const double *GetRefCoordsFromGeometricType(INTERP_KERNEL::NormalizedCellType geoType, std::size_t& lgth);
    MEDCOUPLING_EXPORT static const double *GetLocsFromGeometricType(INTERP_KERNEL::NormalizedCellType geoType, std::size_t& lgth);
  protected:
    MEDCOUPLING_EXPORT MEDCouplingFieldDiscretizationGaussNE(const MEDCouplingFieldDiscretizationGaussNE& other);
  public:
    static const char REPR[];
    static constexpr TypeOfField TYPE = ON_GAUSS_NE;
    static const double FGP_POINT1[1];
    static const double FGP_SEG2[2];
    static const double FGP_SEG3[3];
    static const double FGP_SEG4[4];
    static const double FGP_TRI3[3];
    static const double FGP_TRI6[6];
    static const double FGP_TRI7[7];
    static const double FGP_QUAD4[4];
    static const double FGP_QUAD8[8];
    static const double FGP_QUAD9[9];
    static const double FGP_TETRA4[4];
    static const double FGP_TETRA10[10];//to check
    static const double FGP_PENTA6[6];
    static const double FGP_PENTA15[15];//to check
    static const double FGP_PENTA18[18];//to check
    static const double FGP_HEXA8[8];
    static const double FGP_HEXA20[20];//to check
    static const double FGP_HEXA27[27];
    static const double FGP_PYRA5[5];
    static const double FGP_PYRA13[13];//to check
    static const double REF_SEG2[2];
    static const double REF_SEG3[3];
    static const double REF_SEG4[4];
    static const double REF_TRI3[6];
    static const double REF_TRI6[12];
    static const double REF_TRI7[14];
    static const double REF_QUAD4[8];
    static const double REF_QUAD8[16];
    static const double REF_QUAD9[18];
    static const double REF_TETRA4[12];
    static const double REF_TETRA10[30];
    static const double REF_PENTA6[18];
    static const double REF_PENTA15[45];
    static const double REF_PENTA18[54];
    static const double REF_HEXA8[24];
    static const double REF_HEXA20[60];
    static const double REF_HEXA27[81];
    static const double REF_PYRA5[15];
    static const double REF_PYRA13[39];
    static const double LOC_SEG2[2];
    static const double LOC_SEG3[3];
    static const double LOC_SEG4[4];
    static const double LOC_TRI3[6];
    static const double LOC_TRI6[12];
    static const double LOC_TRI7[14];
    static const double LOC_QUAD4[8];
    static const double LOC_QUAD8[16];
    static const double LOC_QUAD9[18];
    static const double LOC_TETRA4[12];
    static const double LOC_TETRA10[30];//to check
    static const double LOC_PENTA6[18];
    static const double LOC_PENTA15[45];//to check
    static const double LOC_PENTA18[54];//to check
    static const double LOC_HEXA8[24];
    static const double LOC_HEXA20[60];//to check
    static const double LOC_HEXA27[81];
    static const double LOC_PYRA5[15];
    static const double LOC_PYRA13[39];//to check
  };

  class MEDCouplingFieldDiscretizationKriging : public MEDCouplingFieldDiscretizationOnNodes
  {
  public:
    MEDCOUPLING_EXPORT TypeOfField getEnum() const;
    std::string getClassName() const override { return std::string("MEDCouplingFieldDiscretizationKriging"); }
    MEDCOUPLING_EXPORT const char *getRepr() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDiscretization *clone() const override;
    MEDCOUPLING_EXPORT std::string getStringRepr() const;
    MEDCOUPLING_EXPORT void checkCompatibilityWithNature(NatureOfField nat) const override;
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const override;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const override;
    MEDCOUPLING_EXPORT void getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, double *res) const override;
    MEDCOUPLING_EXPORT DataArrayDouble *getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, mcIdType nbOfPoints) const override;
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
    MEDCOUPLING_EXPORT MCAuto<MEDCouplingFieldDiscretization> aggregate(std::vector<const MEDCouplingFieldDiscretization *>& fds) const override;
  public://specific part
    MEDCOUPLING_EXPORT DataArrayDouble *computeEvaluationMatrixOnGivenPts(const MEDCouplingMesh *mesh, const double *loc, mcIdType nbOfTargetPoints, mcIdType& nbCols) const;
    MEDCOUPLING_EXPORT DataArrayDouble *computeInverseMatrix(const MEDCouplingMesh *mesh, mcIdType& isDrift, mcIdType& matSz) const;
    MEDCOUPLING_EXPORT DataArrayDouble *computeMatrix(const MEDCouplingMesh *mesh, mcIdType& isDrift, mcIdType& matSz) const;
    MEDCOUPLING_EXPORT DataArrayDouble *computeVectorOfCoefficients(const MEDCouplingMesh *mesh, const DataArrayDouble *arr, mcIdType& isDrift) const;
    MEDCOUPLING_EXPORT void operateOnDenseMatrix(int spaceDimension, mcIdType nbOfElems, double *matrixPtr) const;
    MEDCOUPLING_EXPORT DataArrayDouble *performDrift(const DataArrayDouble *matr, const DataArrayDouble *arr, mcIdType& delta) const;
    MEDCOUPLING_EXPORT static void OperateOnDenseMatrixH3(mcIdType nbOfElems, double *matrixPtr);
    MEDCOUPLING_EXPORT static void OperateOnDenseMatrixH2Ln(mcIdType nbOfElems, double *matrixPtr);
    MEDCOUPLING_EXPORT static DataArrayDouble *PerformDriftRect(const DataArrayDouble *matr, const DataArrayDouble *arr, mcIdType& delta);
    MEDCOUPLING_EXPORT static DataArrayDouble *PerformDriftOfVec(const DataArrayDouble *arr, mcIdType isDrift);
  public:
    static const char REPR[];
    static constexpr TypeOfField TYPE = ON_NODES_KR;
  };
  
  template<class FIELD_DISC>
  MCAuto<MEDCouplingFieldDiscretization> MEDCouplingFieldDiscretization::EasyAggregate(std::vector<const MEDCouplingFieldDiscretization *>& fds)
  {
    if(fds.empty())
      throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretization::aggregate : input array is empty");
    for(const MEDCouplingFieldDiscretization * it : fds)
      {
        const FIELD_DISC *itc(dynamic_cast<const FIELD_DISC *>(it));
        if(!itc)
          throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretization::aggregate : same field discretization expected for all input discretizations !");
      }
    return fds[0]->clone();
  }
}

#endif
