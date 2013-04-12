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

#ifndef __PARAMEDMEM_MEDCOUPLINGFIELDDISCRETIZATION_HXX__
#define __PARAMEDMEM_MEDCOUPLINGFIELDDISCRETIZATION_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "InterpKernelException.hxx"
#include "MEDCouplingTimeLabel.hxx"
#include "MEDCouplingNatureOfField.hxx"
#include "MEDCouplingGaussLocalization.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"

#include <set>
#include <vector>

namespace ParaMEDMEM
{
  class DataArrayInt;
  class MEDCouplingMesh;
  class DataArrayDouble;
  class MEDCouplingFieldDouble;

  class MEDCOUPLING_EXPORT MEDCouplingFieldDiscretization : public RefCountObject, public TimeLabel
  {
  public:
    static MEDCouplingFieldDiscretization *New(TypeOfField type);
    double getPrecision() const { return _precision; }
    void setPrecision(double val) { _precision=val; }
    void updateTime() const;
    std::size_t getHeapMemorySize() const;
    static TypeOfField GetTypeOfFieldFromStringRepr(const char *repr) throw(INTERP_KERNEL::Exception);
    virtual TypeOfField getEnum() const = 0;
    virtual bool isEqual(const MEDCouplingFieldDiscretization *other, double eps) const;
    virtual bool isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const = 0;
    virtual bool isEqualWithoutConsideringStr(const MEDCouplingFieldDiscretization *other, double eps) const;
    virtual MEDCouplingFieldDiscretization *deepCpy() const;
    virtual MEDCouplingFieldDiscretization *clone() const = 0;
    virtual MEDCouplingFieldDiscretization *clonePart(const int *startCellIds, const int *endCellIds) const;
    virtual MEDCouplingFieldDiscretization *clonePartRange(int beginCellIds, int endCellIds, int stepCellIds) const;
    virtual std::string getStringRepr() const = 0;
    virtual const char *getRepr() const = 0;
    virtual int getNumberOfTuples(const MEDCouplingMesh *mesh) const throw(INTERP_KERNEL::Exception) = 0;
    virtual int getNumberOfMeshPlaces(const MEDCouplingMesh *mesh) const = 0;
    virtual DataArrayInt *getOffsetArr(const MEDCouplingMesh *mesh) const = 0;
    virtual void normL1(const MEDCouplingMesh *mesh, const DataArrayDouble *arr, double *res) const throw(INTERP_KERNEL::Exception);
    virtual void normL2(const MEDCouplingMesh *mesh, const DataArrayDouble *arr, double *res) const throw(INTERP_KERNEL::Exception);
    virtual void integral(const MEDCouplingMesh *mesh, const DataArrayDouble *arr, bool isWAbs, double *res) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayDouble *getLocalizationOfDiscValues(const MEDCouplingMesh *mesh) const = 0;
    virtual void computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, const int *tupleIdsBg, const int *tupleIdsEnd,
                                                    DataArrayInt *&cellRestriction, DataArrayInt *&trueTupleRestriction) const throw(INTERP_KERNEL::Exception) = 0;
    virtual void checkCompatibilityWithNature(NatureOfField nat) const throw(INTERP_KERNEL::Exception) = 0;
    virtual void renumberCells(const int *old2NewBg, bool check=true) throw(INTERP_KERNEL::Exception);
    virtual void renumberArraysForCell(const MEDCouplingMesh *mesh, const std::vector<DataArrayDouble *>& arrays,
                                       const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception) = 0;
    virtual double getIJK(const MEDCouplingMesh *mesh, const DataArrayDouble *da, int cellId, int nodeIdInCell, int compoId) const throw(INTERP_KERNEL::Exception);
    virtual void checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArrayDouble *da) const throw(INTERP_KERNEL::Exception) = 0;
    virtual MEDCouplingFieldDouble *getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const = 0;
    virtual void getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, double *res) const = 0;
    virtual void getValueOnPos(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, int i, int j, int k, double *res) const = 0;
    virtual DataArrayDouble *getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, int nbOfPoints) const = 0;
    virtual DataArrayInt *computeTupleIdsToSelectFromCellIds(const MEDCouplingMesh *mesh, const int *startCellIds, const int *endCellIds) const = 0;
    virtual MEDCouplingMesh *buildSubMeshData(const MEDCouplingMesh *mesh, const int *start, const int *end, DataArrayInt *&di) const = 0;
    virtual MEDCouplingMesh *buildSubMeshDataRange(const MEDCouplingMesh *mesh, int beginCellIds, int endCellIds, int stepCellIds, int& beginOut, int& endOut, int& stepOut, DataArrayInt *&di) const;
    virtual void renumberValuesOnNodes(double epsOnVals, const int *old2New, int newNbOfNodes, DataArrayDouble *arr) const = 0;
    virtual void renumberValuesOnCells(double epsOnVals, const MEDCouplingMesh *mesh, const int *old2New, int newSz, DataArrayDouble *arr) const = 0;
    virtual void renumberValuesOnCellsR(const MEDCouplingMesh *mesh, const int *new2old, int newSz, DataArrayDouble *arr) const = 0;
    virtual void getSerializationIntArray(DataArrayInt *& arr) const;
    virtual void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    virtual void getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const;
    virtual void finishUnserialization(const std::vector<double>& tinyInfo);
    virtual void resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *& arr);
    virtual void setGaussLocalizationOnType(const MEDCouplingMesh *m, INTERP_KERNEL::NormalizedCellType type, const std::vector<double>& refCoo,
                                            const std::vector<double>& gsCoo, const std::vector<double>& wg) throw(INTERP_KERNEL::Exception);
    virtual void setGaussLocalizationOnCells(const MEDCouplingMesh *m, const int *begin, const int *end, const std::vector<double>& refCoo,
                                             const std::vector<double>& gsCoo, const std::vector<double>& wg) throw(INTERP_KERNEL::Exception);
    virtual void clearGaussLocalizations() throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingGaussLocalization& getGaussLocalization(int locId) throw(INTERP_KERNEL::Exception);
    virtual int getNbOfGaussLocalization() const throw(INTERP_KERNEL::Exception);
    virtual int getGaussLocalizationIdOfOneCell(int cellId) const throw(INTERP_KERNEL::Exception);
    virtual int getGaussLocalizationIdOfOneType(INTERP_KERNEL::NormalizedCellType type) const throw(INTERP_KERNEL::Exception);
    virtual std::set<int> getGaussLocalizationIdsOfOneType(INTERP_KERNEL::NormalizedCellType type) const throw(INTERP_KERNEL::Exception);
    virtual void getCellIdsHavingGaussLocalization(int locId, std::vector<int>& cellIds) const throw(INTERP_KERNEL::Exception);
    virtual const MEDCouplingGaussLocalization& getGaussLocalization(int locId) const throw(INTERP_KERNEL::Exception);
    virtual void reprQuickOverview(std::ostream& stream) const throw(INTERP_KERNEL::Exception) = 0;
    virtual ~MEDCouplingFieldDiscretization();
  protected:
    MEDCouplingFieldDiscretization();
    static void RenumberEntitiesFromO2NArr(double epsOnVals, const int *old2NewPtr, int newNbOfEntity, DataArrayDouble *arr, const char *msg);
    static void RenumberEntitiesFromN2OArr(const int *new2OldPtr, int new2OldSz, DataArrayDouble *arr, const char *msg);
  protected:
    double _precision;
    static const double DFLT_PRECISION;
  };

  class MEDCOUPLING_EXPORT MEDCouplingFieldDiscretizationP0 : public MEDCouplingFieldDiscretization
  {
  public:
    TypeOfField getEnum() const;
    MEDCouplingFieldDiscretization *clone() const;
    std::string getStringRepr() const;
    const char *getRepr() const;
    bool isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const;
    int getNumberOfTuples(const MEDCouplingMesh *mesh) const throw(INTERP_KERNEL::Exception);
    int getNumberOfMeshPlaces(const MEDCouplingMesh *mesh) const;
    DataArrayInt *getOffsetArr(const MEDCouplingMesh *mesh) const;
    void renumberArraysForCell(const MEDCouplingMesh *mesh, const std::vector<DataArrayDouble *>& arrays,
                               const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getLocalizationOfDiscValues(const MEDCouplingMesh *mesh) const;
    void checkCompatibilityWithNature(NatureOfField nat) const throw(INTERP_KERNEL::Exception);
    void computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, const int *tupleIdsBg, const int *tupleIdsEnd,
                                            DataArrayInt *&cellRestriction, DataArrayInt *&trueTupleRestriction) const throw(INTERP_KERNEL::Exception);
    void checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArrayDouble *da) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const;
    void getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, double *res) const;
    void getValueOnPos(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, int i, int j, int k, double *res) const;
    DataArrayDouble *getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, int nbOfPoints) const;
    void renumberValuesOnNodes(double epsOnVals, const int *old2New, int newNbOfNodes, DataArrayDouble *arr) const;
    void renumberValuesOnCells(double epsOnVals, const MEDCouplingMesh *mesh, const int *old2New, int newSz, DataArrayDouble *arr) const;
    void renumberValuesOnCellsR(const MEDCouplingMesh *mesh, const int *new2old, int newSz, DataArrayDouble *arr) const;
    MEDCouplingMesh *buildSubMeshData(const MEDCouplingMesh *mesh, const int *start, const int *end, DataArrayInt *&di) const;
    MEDCouplingMesh *buildSubMeshDataRange(const MEDCouplingMesh *mesh, int beginCellIds, int endCellIds, int stepCellIds, int& beginOut, int& endOut, int& stepOut, DataArrayInt *&di) const;
    DataArrayInt *computeTupleIdsToSelectFromCellIds(const MEDCouplingMesh *mesh, const int *startCellIds, const int *endCellIds) const;
    void reprQuickOverview(std::ostream& stream) const throw(INTERP_KERNEL::Exception);
  public:
    static const char REPR[];
    static const TypeOfField TYPE;
  };

  class MEDCOUPLING_EXPORT MEDCouplingFieldDiscretizationOnNodes : public MEDCouplingFieldDiscretization
  {
  public:
    int getNumberOfTuples(const MEDCouplingMesh *mesh) const throw(INTERP_KERNEL::Exception);
    int getNumberOfMeshPlaces(const MEDCouplingMesh *mesh) const;
    DataArrayInt *getOffsetArr(const MEDCouplingMesh *mesh) const;
    void renumberArraysForCell(const MEDCouplingMesh *mesh, const std::vector<DataArrayDouble *>& arrays,
                               const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getLocalizationOfDiscValues(const MEDCouplingMesh *mesh) const;
    void computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, const int *tupleIdsBg, const int *tupleIdsEnd,
                                            DataArrayInt *&cellRestriction, DataArrayInt *&trueTupleRestriction) const throw(INTERP_KERNEL::Exception);
    void checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArrayDouble *da) const throw(INTERP_KERNEL::Exception);
    MEDCouplingMesh *buildSubMeshData(const MEDCouplingMesh *mesh, const int *start, const int *end, DataArrayInt *&di) const;
    MEDCouplingMesh *buildSubMeshDataRange(const MEDCouplingMesh *mesh, int beginCellIds, int endCellIds, int stepCellIds, int& beginOut, int& endOut, int& stepOut, DataArrayInt *&di) const;
    DataArrayInt *computeTupleIdsToSelectFromCellIds(const MEDCouplingMesh *mesh, const int *startCellIds, const int *endCellIds) const;
    void renumberValuesOnNodes(double epsOnVals, const int *old2New, int newNbOfNodes, DataArrayDouble *arr) const;
    void renumberValuesOnCells(double epsOnVals, const MEDCouplingMesh *mesh, const int *old2New, int newSz, DataArrayDouble *arr) const;
    void renumberValuesOnCellsR(const MEDCouplingMesh *mesh, const int *new2old, int newSz, DataArrayDouble *arr) const;
  public:
    void getValueOnPos(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, int i, int j, int k, double *res) const;
  };

  class MEDCOUPLING_EXPORT MEDCouplingFieldDiscretizationP1 : public MEDCouplingFieldDiscretizationOnNodes
  {
  public:
    TypeOfField getEnum() const;
    MEDCouplingFieldDiscretization *clone() const;
    std::string getStringRepr() const;
    const char *getRepr() const;
    void checkCompatibilityWithNature(NatureOfField nat) const throw(INTERP_KERNEL::Exception);
    bool isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const;
    MEDCouplingFieldDouble *getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const;
    void getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, double *res) const;
    DataArrayDouble *getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, int nbOfPoints) const;
    void reprQuickOverview(std::ostream& stream) const throw(INTERP_KERNEL::Exception);
  public:
    static const char REPR[];
    static const TypeOfField TYPE;
  protected:
    void getValueInCell(const MEDCouplingMesh *mesh, int cellId, const DataArrayDouble *arr, const double *loc, double *res) const;
  };

  /*!
   * This class abstracts MEDCouplingFieldDiscretization that needs an information on each cell to perform their job.
   * All classes that inherits from this are more linked to mesh.
   */
  class MEDCOUPLING_EXPORT MEDCouplingFieldDiscretizationPerCell : public MEDCouplingFieldDiscretization
  {
  public:
    const DataArrayInt *getArrayOfDiscIds() const;
    void setArrayOfDiscIds(const DataArrayInt *adids) throw(INTERP_KERNEL::Exception);
    void checkNoOrphanCells() const throw(INTERP_KERNEL::Exception);
    std::vector<DataArrayInt *> splitIntoSingleGaussDicrPerCellType(std::vector< int >& locIds) const throw(INTERP_KERNEL::Exception);
  protected:
    MEDCouplingFieldDiscretizationPerCell();
    MEDCouplingFieldDiscretizationPerCell(const MEDCouplingFieldDiscretizationPerCell& other, const int *startCellIds, const int *endCellIds);
    MEDCouplingFieldDiscretizationPerCell(const MEDCouplingFieldDiscretizationPerCell& other, int beginCellIds, int endCellIds, int stepCellIds);
    ~MEDCouplingFieldDiscretizationPerCell();
    void updateTime() const;
    std::size_t getHeapMemorySize() const;
    void checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArrayDouble *da) const throw(INTERP_KERNEL::Exception);
    bool isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const;
    bool isEqualWithoutConsideringStr(const MEDCouplingFieldDiscretization *other, double eps) const;
    void renumberCells(const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception);
  protected:
    void buildDiscrPerCellIfNecessary(const MEDCouplingMesh *mesh);
  protected:
    DataArrayInt *_discr_per_cell;
    static const int DFT_INVALID_LOCID_VALUE;
  };

  class MEDCOUPLING_EXPORT MEDCouplingFieldDiscretizationGauss : public MEDCouplingFieldDiscretizationPerCell
  {
  public:
    MEDCouplingFieldDiscretizationGauss();
    TypeOfField getEnum() const;
    bool isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const;
    bool isEqualWithoutConsideringStr(const MEDCouplingFieldDiscretization *other, double eps) const;
    MEDCouplingFieldDiscretization *clone() const;
    MEDCouplingFieldDiscretization *clonePart(const int *startCellIds, const int *endCellIds) const;
    MEDCouplingFieldDiscretization *clonePartRange(int beginCellIds, int endCellIds, int stepCellIds) const;
    std::string getStringRepr() const;
    const char *getRepr() const;
    std::size_t getHeapMemorySize() const;
    int getNumberOfTuples(const MEDCouplingMesh *mesh) const throw(INTERP_KERNEL::Exception);
    int getNumberOfMeshPlaces(const MEDCouplingMesh *mesh) const;
    DataArrayInt *getOffsetArr(const MEDCouplingMesh *mesh) const;
    void renumberArraysForCell(const MEDCouplingMesh *mesh, const std::vector<DataArrayDouble *>& arrays,
                               const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getLocalizationOfDiscValues(const MEDCouplingMesh *mesh) const;
    void computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, const int *tupleIdsBg, const int *tupleIdsEnd,
                                            DataArrayInt *&cellRestriction, DataArrayInt *&trueTupleRestriction) const throw(INTERP_KERNEL::Exception);
    void checkCompatibilityWithNature(NatureOfField nat) const throw(INTERP_KERNEL::Exception);
    void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    void getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const;
    void finishUnserialization(const std::vector<double>& tinyInfo);
    void getSerializationIntArray(DataArrayInt *& arr) const;
    void resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *& arr);
    double getIJK(const MEDCouplingMesh *mesh, const DataArrayDouble *da, int cellId, int nodeIdInCell, int compoId) const throw(INTERP_KERNEL::Exception);
    void checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArrayDouble *da) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const;
    void getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, double *res) const;
    void getValueOnPos(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, int i, int j, int k, double *res) const;
    DataArrayDouble *getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, int nbOfPoints) const;
    MEDCouplingMesh *buildSubMeshData(const MEDCouplingMesh *mesh, const int *start, const int *end, DataArrayInt *&di) const;
    MEDCouplingMesh *buildSubMeshDataRange(const MEDCouplingMesh *mesh, int beginCellIds, int endCellIds, int stepCellIds, int& beginOut, int& endOut, int& stepOut, DataArrayInt *&di) const;
    DataArrayInt *computeTupleIdsToSelectFromCellIds(const MEDCouplingMesh *mesh, const int *startCellIds, const int *endCellIds) const;
    void renumberValuesOnNodes(double epsOnVals, const int *old2New, int newNbOfNodes, DataArrayDouble *arr) const;
    void renumberValuesOnCells(double epsOnVals, const MEDCouplingMesh *mesh, const int *old2New, int newSz, DataArrayDouble *arr) const;
    void renumberValuesOnCellsR(const MEDCouplingMesh *mesh, const int *new2old, int newSz, DataArrayDouble *arr) const;
    void setGaussLocalizationOnType(const MEDCouplingMesh *mesh, INTERP_KERNEL::NormalizedCellType type, const std::vector<double>& refCoo,
                                    const std::vector<double>& gsCoo, const std::vector<double>& wg) throw(INTERP_KERNEL::Exception);
    void setGaussLocalizationOnCells(const MEDCouplingMesh *mesh, const int *begin, const int *end, const std::vector<double>& refCoo,
                                     const std::vector<double>& gsCoo, const std::vector<double>& wg) throw(INTERP_KERNEL::Exception);
    void clearGaussLocalizations() throw(INTERP_KERNEL::Exception);
    void setGaussLocalization(int locId, const MEDCouplingGaussLocalization& loc) throw(INTERP_KERNEL::Exception);
    void resizeLocalizationVector(int newSz) throw(INTERP_KERNEL::Exception);
    MEDCouplingGaussLocalization& getGaussLocalization(int locId) throw(INTERP_KERNEL::Exception);
    int getNbOfGaussLocalization() const throw(INTERP_KERNEL::Exception);
    int getGaussLocalizationIdOfOneCell(int cellId) const throw(INTERP_KERNEL::Exception);
    int getGaussLocalizationIdOfOneType(INTERP_KERNEL::NormalizedCellType type) const throw(INTERP_KERNEL::Exception);
    std::set<int> getGaussLocalizationIdsOfOneType(INTERP_KERNEL::NormalizedCellType type) const throw(INTERP_KERNEL::Exception);
    void getCellIdsHavingGaussLocalization(int locId, std::vector<int>& cellIds) const throw(INTERP_KERNEL::Exception);
    const MEDCouplingGaussLocalization& getGaussLocalization(int locId) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *buildNbOfGaussPointPerCellField() const throw(INTERP_KERNEL::Exception);
    void reprQuickOverview(std::ostream& stream) const throw(INTERP_KERNEL::Exception);
  protected:
    MEDCouplingFieldDiscretizationGauss(const MEDCouplingFieldDiscretizationGauss& other, const int *startCellIds=0, const int *endCellIds=0);
    MEDCouplingFieldDiscretizationGauss(const MEDCouplingFieldDiscretizationGauss& other, int beginCellIds, int endCellIds, int stepCellIds);
    void zipGaussLocalizations();
    int getOffsetOfCell(int cellId) const throw(INTERP_KERNEL::Exception);
    void checkLocalizationId(int locId) const throw(INTERP_KERNEL::Exception);
  public:
    static const char REPR[];
    static const TypeOfField TYPE;
  private:
    std::vector<MEDCouplingGaussLocalization> _loc;
  };

  /*!
   * Gauss with points of values located on nodes of element. This is a specialization of MEDCouplingFieldDiscretizationGauss.
   */
  class MEDCOUPLING_EXPORT MEDCouplingFieldDiscretizationGaussNE : public MEDCouplingFieldDiscretization
  {
  public:
    MEDCouplingFieldDiscretizationGaussNE();
    TypeOfField getEnum() const;
    MEDCouplingFieldDiscretization *clone() const;
    std::string getStringRepr() const;
    const char *getRepr() const;
    bool isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const;
    int getNumberOfTuples(const MEDCouplingMesh *mesh) const throw(INTERP_KERNEL::Exception);
    int getNumberOfMeshPlaces(const MEDCouplingMesh *mesh) const;
    DataArrayInt *getOffsetArr(const MEDCouplingMesh *mesh) const;
    void renumberArraysForCell(const MEDCouplingMesh *mesh, const std::vector<DataArrayDouble *>& arrays,
                               const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getLocalizationOfDiscValues(const MEDCouplingMesh *mesh) const;
    void integral(const MEDCouplingMesh *mesh, const DataArrayDouble *arr, bool isWAbs, double *res) const throw(INTERP_KERNEL::Exception);
    void computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, const int *tupleIdsBg, const int *tupleIdsEnd,
                                            DataArrayInt *&cellRestriction, DataArrayInt *&trueTupleRestriction) const throw(INTERP_KERNEL::Exception);
    void checkCompatibilityWithNature(NatureOfField nat) const throw(INTERP_KERNEL::Exception);
    double getIJK(const MEDCouplingMesh *mesh, const DataArrayDouble *da, int cellId, int nodeIdInCell, int compoId) const throw(INTERP_KERNEL::Exception);
    void checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArrayDouble *da) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const;
    void getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, double *res) const;
    void getValueOnPos(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, int i, int j, int k, double *res) const;
    DataArrayDouble *getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, int nbOfPoints) const;
    MEDCouplingMesh *buildSubMeshData(const MEDCouplingMesh *mesh, const int *start, const int *end, DataArrayInt *&di) const;
    MEDCouplingMesh *buildSubMeshDataRange(const MEDCouplingMesh *mesh, int beginCellIds, int endCellIds, int stepCellIds, int& beginOut, int& endOut, int& stepOut, DataArrayInt *&di) const;
    DataArrayInt *computeTupleIdsToSelectFromCellIds(const MEDCouplingMesh *mesh, const int *startCellIds, const int *endCellIds) const;
    void renumberValuesOnNodes(double epsOnVals, const int *old2New, int newNbOfNodes, DataArrayDouble *arr) const;
    void renumberValuesOnCells(double epsOnVals, const MEDCouplingMesh *mesh, const int *old2New, int newSz, DataArrayDouble *arr) const;
    void renumberValuesOnCellsR(const MEDCouplingMesh *mesh, const int *new2old, int newSz, DataArrayDouble *arr) const;
    void reprQuickOverview(std::ostream& stream) const throw(INTERP_KERNEL::Exception);
    static const double *GetWeightArrayFromGeometricType(INTERP_KERNEL::NormalizedCellType geoType, std::size_t& lgth) throw(INTERP_KERNEL::Exception);
    static const double *GetRefCoordsFromGeometricType(INTERP_KERNEL::NormalizedCellType geoType, std::size_t& lgth) throw(INTERP_KERNEL::Exception);
  protected:
    MEDCouplingFieldDiscretizationGaussNE(const MEDCouplingFieldDiscretizationGaussNE& other);
  public:
    static const char REPR[];
    static const TypeOfField TYPE;
    static const double FGP_SEG2[2];
    static const double FGP_SEG3[3];
    static const double FGP_SEG4[4];
    static const double FGP_TRI3[3];
    static const double FGP_TRI6[6];
    static const double FGP_TRI7[7];
    static const double FGP_QUAD4[4];
    //static const double FGP_QUAD8[8];
    static const double FGP_QUAD9[9];
    static const double FGP_TETRA4[4];
    //static const double FGP_TETRA10[10];
    static const double FGP_PENTA6[6];
    //static const double FGP_PENTA15[15];
    static const double FGP_HEXA8[8];
    static const double FGP_HEXA27[27];
    static const double FGP_PYRA5[5];
    //static const double FGP_PYRA13[13];
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
    static const double REF_HEXA8[24];
    static const double REF_HEXA20[60];
    static const double REF_HEXA27[81];
    static const double REF_PYRA5[15];
    static const double REF_PYRA13[39];
  };

  class MEDCOUPLING_EXPORT MEDCouplingFieldDiscretizationKriging : public MEDCouplingFieldDiscretizationOnNodes
  {
  public:
    TypeOfField getEnum() const;
    const char *getRepr() const;
    MEDCouplingFieldDiscretization *clone() const;
    std::string getStringRepr() const;
    void checkCompatibilityWithNature(NatureOfField nat) const throw(INTERP_KERNEL::Exception);
    bool isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const;
    MEDCouplingFieldDouble *getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const;
    void getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, double *res) const;
    DataArrayDouble *getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, int nbOfPoints) const;
    void reprQuickOverview(std::ostream& stream) const throw(INTERP_KERNEL::Exception);
  public://specific part
    DataArrayDouble *computeVectorOfCoefficients(const MEDCouplingMesh *mesh, const DataArrayDouble *arr, int& isDrift) const;
  protected:
    void operateOnDenseMatrix(int spaceDimension, int nbOfElems, double *matrixPtr) const;
    DataArrayDouble *performDrift(const DataArrayDouble *matr, const DataArrayDouble *arr, int& delta) const;
  public:
    static const char REPR[];
    static const TypeOfField TYPE;
  };
}

#endif
