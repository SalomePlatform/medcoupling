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

#ifndef __PARAMEDMEM_MEDCOUPLINGTIMEDISCRETIZATION_HXX__
#define __PARAMEDMEM_MEDCOUPLINGTIMEDISCRETIZATION_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingTimeLabel.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "InterpKernelException.hxx"

#include <vector>

namespace ParaMEDMEM
{
  class MEDCouplingMesh;
  class DataArrayDouble;
  class TimeLabel;

  class MEDCouplingTimeDiscretization : public TimeLabel, public BigMemoryObject
  {
  protected:
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization();
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization(const MEDCouplingTimeDiscretization& other, bool deepCpy);
  public:
    MEDCOUPLING_EXPORT void updateTime() const;
    MEDCOUPLING_EXPORT virtual std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT virtual std::vector<const BigMemoryObject *> getDirectChildren() const;
    MEDCOUPLING_EXPORT static MEDCouplingTimeDiscretization *New(TypeOfTimeDiscretization type) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setTimeUnit(const char *unit) { _time_unit=unit; }
    MEDCOUPLING_EXPORT const char *getTimeUnit() const { return _time_unit.c_str(); }
    MEDCOUPLING_EXPORT virtual void copyTinyAttrFrom(const MEDCouplingTimeDiscretization& other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void copyTinyStringsFrom(const MEDCouplingTimeDiscretization& other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void checkCoherency() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual bool areCompatible(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual bool areStrictlyCompatible(const MEDCouplingTimeDiscretization *other, std::string& reason) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual bool areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual bool areStrictlyCompatibleForDiv(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual bool areCompatibleForMeld(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual bool isEqualIfNotWhy(const MEDCouplingTimeDiscretization *other, double prec, std::string& reason) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual bool isEqual(const MEDCouplingTimeDiscretization *other, double prec) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual bool isEqualWithoutConsideringStr(const MEDCouplingTimeDiscretization *other, double prec) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *buildNewTimeReprFromThis(TypeOfTimeDiscretization type, bool deepCpy) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual std::string getStringRepr() const throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual TypeOfTimeDiscretization getEnum() const throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual void synchronizeTimeWith(const MEDCouplingMesh *mesh) throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *aggregate(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *aggregate(const std::vector<const MEDCouplingTimeDiscretization *>& other) const throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *meld(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *dot(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *crossProduct(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *max(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *min(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *add(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual void addEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *substract(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual void substractEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *multiply(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual void multiplyEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *divide(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual void divideEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *pow(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual void powEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void resizeForUnserialization(const std::vector<int>& tinyInfoI, std::vector<DataArrayDouble *>& arrays) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void getTinySerializationIntInformation2(std::vector<int>& tinyInfo) const throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual void getTinySerializationDbleInformation2(std::vector<double>& tinyInfo) const throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual void finishUnserialization2(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD) throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *performCpy(bool deepCpy) const throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT void setTimeTolerance(double val) { _time_tolerance=val; }
    MEDCOUPLING_EXPORT double getTimeTolerance() const { return _time_tolerance; }
    MEDCOUPLING_EXPORT virtual void checkNoTimePresence() const throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual void checkTimePresence(double time) const throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual void setArray(DataArrayDouble *array, TimeLabel *owner) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void setEndArray(DataArrayDouble *array, TimeLabel *owner) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void setArrays(const std::vector<DataArrayDouble *>& arrays, TimeLabel *owner) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *getArray() throw(INTERP_KERNEL::Exception) { return _array; }
    MEDCOUPLING_EXPORT const DataArrayDouble *getArray() const throw(INTERP_KERNEL::Exception) { return _array; }
    MEDCOUPLING_EXPORT virtual const DataArrayDouble *getEndArray() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual DataArrayDouble *getEndArray() throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual std::vector< const DataArrayDouble *> getArraysForTime(double time) const throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual void getValueForTime(double time, const std::vector<double>& vals, double *res) const throw(INTERP_KERNEL::Exception) = 0; 
    MEDCOUPLING_EXPORT virtual void getArrays(std::vector<DataArrayDouble *>& arrays) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual bool isBefore(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual bool isStrictlyBefore(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT double getTime(int& iteration, int& order) const throw(INTERP_KERNEL::Exception) { return getStartTime(iteration,order); }
    MEDCOUPLING_EXPORT virtual double getStartTime(int& iteration, int& order) const throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual double getEndTime(int& iteration, int& order) const throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT void setTime(double time, int iteration, int order) throw(INTERP_KERNEL::Exception) { setStartTime(time,iteration,order); }
    MEDCOUPLING_EXPORT void setIteration(int it) throw(INTERP_KERNEL::Exception) { setStartIteration(it); }
    MEDCOUPLING_EXPORT void setOrder(int order) throw(INTERP_KERNEL::Exception) { setStartOrder(order); }
    MEDCOUPLING_EXPORT void setTimeValue(double val) throw(INTERP_KERNEL::Exception) { setStartTimeValue(val); }
    MEDCOUPLING_EXPORT virtual void setStartIteration(int it) throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual void setEndIteration(int it) throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual void setStartOrder(int order) throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual void setEndOrder(int order) throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual void setStartTimeValue(double time) throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual void setEndTimeValue(double time) throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual void setStartTime(double time, int iteration, int order) throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual void setEndTime(double time, int iteration, int order) throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual void getValueOnTime(int eltId, double time, double *value) const throw(INTERP_KERNEL::Exception) = 0;
    MEDCOUPLING_EXPORT virtual void getValueOnDiscTime(int eltId, int iteration, int order, double *value) const throw(INTERP_KERNEL::Exception) = 0;
    //
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *doublyContractedProduct() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *determinant() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *eigenValues() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *eigenVectors() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *inverse() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *trace() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *deviator() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *magnitude() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *negate() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *maxPerTuple() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *keepSelectedComponents(const std::vector<int>& compoIds) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void setSelectedComponents(const MEDCouplingTimeDiscretization *other, const std::vector<int>& compoIds) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void changeNbOfComponents(int newNbOfComp, double dftValue) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void sortPerTuple(bool asc) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void setUniformValue(int nbOfTuple, int nbOfCompo, double value) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void setOrCreateUniformValueOnAllComponents(int nbOfTuple, double value) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void applyLin(double a, double b, int compoId) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void applyFunc(int nbOfComp, FunctionToEvaluate func) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void applyFunc(int nbOfComp, const char *func) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void applyFunc2(int nbOfComp, const char *func) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void applyFunc3(int nbOfComp, const std::vector<std::string>& varsOrder, const char *func) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void applyFunc(const char *func) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void applyFuncFast32(const char *func) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void applyFuncFast64(const char *func) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void fillFromAnalytic(const DataArrayDouble *loc, int nbOfComp, FunctionToEvaluate func) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void fillFromAnalytic(const DataArrayDouble *loc, int nbOfComp, const char *func) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void fillFromAnalytic2(const DataArrayDouble *loc, int nbOfComp, const char *func) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT virtual void fillFromAnalytic3(const DataArrayDouble *loc, int nbOfComp, const std::vector<std::string>& varsOrder, const char *func) throw(INTERP_KERNEL::Exception);
    //
    MEDCOUPLING_EXPORT virtual ~MEDCouplingTimeDiscretization();
  protected:
    std::string _time_unit;
    double _time_tolerance;
    DataArrayDouble *_array;
  protected:
    static const double TIME_TOLERANCE_DFT;
  };

  class MEDCouplingNoTimeLabel : public MEDCouplingTimeDiscretization
  {
  public:
    MEDCOUPLING_EXPORT MEDCouplingNoTimeLabel();
    MEDCOUPLING_EXPORT MEDCouplingNoTimeLabel(const MEDCouplingTimeDiscretization& other, bool deepCpy);
    MEDCOUPLING_EXPORT std::string getStringRepr() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT TypeOfTimeDiscretization getEnum() const throw(INTERP_KERNEL::Exception) { return DISCRETIZATION; }
    MEDCOUPLING_EXPORT void synchronizeTimeWith(const MEDCouplingMesh *mesh) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *aggregate(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *aggregate(const std::vector<const MEDCouplingTimeDiscretization *>& other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *meld(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *dot(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *crossProduct(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *max(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *min(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *add(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void addEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *substract(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void substractEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *multiply(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void multiplyEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *divide(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void divideEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *pow(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void powEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingTimeDiscretization *other, double prec, std::string& reason) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const MEDCouplingTimeDiscretization *other, double prec) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool areCompatible(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool areStrictlyCompatible(const MEDCouplingTimeDiscretization *other, std::string& reason) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool areStrictlyCompatibleForDiv(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool areCompatibleForMeld(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *performCpy(bool deepCpy) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void checkNoTimePresence() const throw(INTERP_KERNEL::Exception) { }
    MEDCOUPLING_EXPORT void checkTimePresence(double time) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT std::vector< const DataArrayDouble *> getArraysForTime(double time) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getValueForTime(double time, const std::vector<double>& vals, double *res) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool isBefore(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool isStrictlyBefore(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT double getStartTime(int& iteration, int& order) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT double getEndTime(int& iteration, int& order) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setStartIteration(int it) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setEndIteration(int it) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setStartOrder(int order) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setEndOrder(int order) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setStartTimeValue(double time) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setEndTimeValue(double time) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setStartTime(double time, int iteration, int order) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setEndTime(double time, int iteration, int order) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getValueOnTime(int eltId, double time, double *value) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getValueOnDiscTime(int eltId, int iteration, int order, double *value) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation2(std::vector<int>& tinyInfo) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getTinySerializationDbleInformation2(std::vector<double>& tinyInfo) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void finishUnserialization2(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD) throw(INTERP_KERNEL::Exception);
  public:
    static const TypeOfTimeDiscretization DISCRETIZATION=NO_TIME;
    static const char REPR[];
  private:
    static const char EXCEPTION_MSG[];
  };

  class MEDCouplingWithTimeStep : public MEDCouplingTimeDiscretization
  {
  protected:
    MEDCOUPLING_EXPORT MEDCouplingWithTimeStep(const MEDCouplingWithTimeStep& other, bool deepCpy);
  public:
    MEDCOUPLING_EXPORT MEDCouplingWithTimeStep();
    MEDCOUPLING_EXPORT std::string getStringRepr() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void copyTinyAttrFrom(const MEDCouplingTimeDiscretization& other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT TypeOfTimeDiscretization getEnum() const throw(INTERP_KERNEL::Exception) { return DISCRETIZATION; }
    MEDCOUPLING_EXPORT void synchronizeTimeWith(const MEDCouplingMesh *mesh) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *aggregate(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *aggregate(const std::vector<const MEDCouplingTimeDiscretization *>& other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *meld(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *dot(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *crossProduct(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *max(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *min(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *add(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void addEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *substract(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void substractEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *multiply(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void multiplyEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *divide(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void divideEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *pow(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void powEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingTimeDiscretization *other, double prec, std::string& reason) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const MEDCouplingTimeDiscretization *other, double prec) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool areCompatible(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool areStrictlyCompatible(const MEDCouplingTimeDiscretization *other, std::string& reason) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool areStrictlyCompatibleForDiv(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool areCompatibleForMeld(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation2(std::vector<int>& tinyInfo) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getTinySerializationDbleInformation2(std::vector<double>& tinyInfo) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void finishUnserialization2(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *performCpy(bool deepCpy) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void checkNoTimePresence() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void checkTimePresence(double time) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setStartTime(double time, int iteration, int order) throw(INTERP_KERNEL::Exception) { _time=time; _iteration=iteration; _order=order; }
    MEDCOUPLING_EXPORT void setEndTime(double time, int iteration, int order) throw(INTERP_KERNEL::Exception) { _time=time; _iteration=iteration; _order=order; }
    MEDCOUPLING_EXPORT double getStartTime(int& iteration, int& order) const throw(INTERP_KERNEL::Exception) { iteration=_iteration; order=_order; return _time; }
    MEDCOUPLING_EXPORT double getEndTime(int& iteration, int& order) const throw(INTERP_KERNEL::Exception) { iteration=_iteration; order=_order; return _time; }
    MEDCOUPLING_EXPORT void setStartIteration(int it) throw(INTERP_KERNEL::Exception) { _iteration=it; }
    MEDCOUPLING_EXPORT void setEndIteration(int it) throw(INTERP_KERNEL::Exception) { _iteration=it; }
    MEDCOUPLING_EXPORT void setStartOrder(int order) throw(INTERP_KERNEL::Exception) { _order=order; }
    MEDCOUPLING_EXPORT void setEndOrder(int order) throw(INTERP_KERNEL::Exception) { _order=order; }
    MEDCOUPLING_EXPORT void setStartTimeValue(double time) throw(INTERP_KERNEL::Exception) { _time=time; }
    MEDCOUPLING_EXPORT void setEndTimeValue(double time) throw(INTERP_KERNEL::Exception) { _time=time; }
    MEDCOUPLING_EXPORT std::vector< const DataArrayDouble *> getArraysForTime(double time) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getValueForTime(double time, const std::vector<double>& vals, double *res) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getValueOnTime(int eltId, double time, double *value) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getValueOnDiscTime(int eltId, int iteration, int order, double *value) const throw(INTERP_KERNEL::Exception);
  public:
    static const TypeOfTimeDiscretization DISCRETIZATION=ONE_TIME;
    static const char REPR[];
  private:
    static const char EXCEPTION_MSG[];
  protected:
    double _time;
    int _iteration;
    int _order;
  };

  class MEDCouplingConstOnTimeInterval : public MEDCouplingTimeDiscretization
  {
  protected:
    MEDCOUPLING_EXPORT MEDCouplingConstOnTimeInterval(const MEDCouplingConstOnTimeInterval& other, bool deepCpy);
  public:
    MEDCOUPLING_EXPORT MEDCouplingConstOnTimeInterval();
    MEDCOUPLING_EXPORT void copyTinyAttrFrom(const MEDCouplingTimeDiscretization& other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation2(std::vector<int>& tinyInfo) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getTinySerializationDbleInformation2(std::vector<double>& tinyInfo) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void finishUnserialization2(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *performCpy(bool deepCpy) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool areCompatible(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool areStrictlyCompatible(const MEDCouplingTimeDiscretization *other, std::string& reason) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool areStrictlyCompatibleForDiv(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool areCompatibleForMeld(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingTimeDiscretization *other, double prec, std::string& reason) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const MEDCouplingTimeDiscretization *other, double prec) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT std::vector< const DataArrayDouble *> getArraysForTime(double time) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getValueForTime(double time, const std::vector<double>& vals, double *res) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getValueOnTime(int eltId, double time, double *value) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getValueOnDiscTime(int eltId, int iteration, int order, double *value) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT TypeOfTimeDiscretization getEnum() const throw(INTERP_KERNEL::Exception) { return DISCRETIZATION; }
    MEDCOUPLING_EXPORT void synchronizeTimeWith(const MEDCouplingMesh *mesh) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT std::string getStringRepr() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *aggregate(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *aggregate(const std::vector<const MEDCouplingTimeDiscretization *>& other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *meld(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *dot(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *crossProduct(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *max(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *min(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *add(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void addEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *substract(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void substractEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *multiply(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void multiplyEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *divide(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void divideEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception);
    MEDCouplingTimeDiscretization *pow(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void powEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setStartTime(double time, int iteration, int order) throw(INTERP_KERNEL::Exception) { _start_time=time; _start_iteration=iteration; _start_order=order; }
    MEDCOUPLING_EXPORT void setEndTime(double time, int iteration, int order) throw(INTERP_KERNEL::Exception) { _end_time=time; _end_iteration=iteration; _end_order=order; }
    MEDCOUPLING_EXPORT double getStartTime(int& iteration, int& order) const throw(INTERP_KERNEL::Exception) { iteration=_start_iteration; order=_start_order; return _start_time; }
    MEDCOUPLING_EXPORT double getEndTime(int& iteration, int& order) const throw(INTERP_KERNEL::Exception) { iteration=_end_iteration; order=_end_order; return _end_time; }
    MEDCOUPLING_EXPORT void setStartIteration(int it) throw(INTERP_KERNEL::Exception) { _start_iteration=it; }
    MEDCOUPLING_EXPORT void setEndIteration(int it) throw(INTERP_KERNEL::Exception) { _end_iteration=it; }
    MEDCOUPLING_EXPORT void setStartOrder(int order) throw(INTERP_KERNEL::Exception) { _start_order=order; }
    MEDCOUPLING_EXPORT void setEndOrder(int order) throw(INTERP_KERNEL::Exception) { _end_order=order; }
    MEDCOUPLING_EXPORT void setStartTimeValue(double time) throw(INTERP_KERNEL::Exception) { _start_time=time; }
    MEDCOUPLING_EXPORT void setEndTimeValue(double time) throw(INTERP_KERNEL::Exception) { _end_time=time; }
    MEDCOUPLING_EXPORT void checkNoTimePresence() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void checkTimePresence(double time) const throw(INTERP_KERNEL::Exception);
  public:
    static const TypeOfTimeDiscretization DISCRETIZATION=CONST_ON_TIME_INTERVAL;
    static const char REPR[];
  private:
    static const char EXCEPTION_MSG[];
  protected:
    double _start_time;
    double _end_time;
    int _start_iteration;
    int _end_iteration;
    int _start_order;
    int _end_order;
  };

  class MEDCouplingTwoTimeSteps : public MEDCouplingTimeDiscretization
  {
  protected:
    MEDCOUPLING_EXPORT MEDCouplingTwoTimeSteps(const MEDCouplingTwoTimeSteps& other, bool deepCpy);
    MEDCOUPLING_EXPORT MEDCouplingTwoTimeSteps();
    MEDCOUPLING_EXPORT ~MEDCouplingTwoTimeSteps();
  public:
    MEDCOUPLING_EXPORT void updateTime() const;
    MEDCOUPLING_EXPORT void synchronizeTimeWith(const MEDCouplingMesh *mesh) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildren() const;
    MEDCOUPLING_EXPORT void copyTinyAttrFrom(const MEDCouplingTimeDiscretization& other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void copyTinyStringsFrom(const MEDCouplingTimeDiscretization& other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT const DataArrayDouble *getEndArray() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *getEndArray() throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void checkCoherency() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingTimeDiscretization *other, double prec, std::string& reason) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const MEDCouplingTimeDiscretization *other, double prec) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void checkNoTimePresence() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void checkTimePresence(double time) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getArrays(std::vector<DataArrayDouble *>& arrays) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setEndArray(DataArrayDouble *array, TimeLabel *owner) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setStartTime(double time, int iteration, int order) throw(INTERP_KERNEL::Exception) { _start_time=time; _start_iteration=iteration; _start_order=order; }
    MEDCOUPLING_EXPORT void setEndTime(double time, int iteration, int order) throw(INTERP_KERNEL::Exception) { _end_time=time; _end_iteration=iteration; _end_order=order; }
    MEDCOUPLING_EXPORT double getStartTime(int& iteration, int& order) const throw(INTERP_KERNEL::Exception) { iteration=_start_iteration; order=_start_order; return _start_time; }
    MEDCOUPLING_EXPORT double getEndTime(int& iteration, int& order) const throw(INTERP_KERNEL::Exception) { iteration=_end_iteration; order=_end_order; return _end_time; }
    MEDCOUPLING_EXPORT void setStartIteration(int it) throw(INTERP_KERNEL::Exception) { _start_iteration=it; }
    MEDCOUPLING_EXPORT void setEndIteration(int it) throw(INTERP_KERNEL::Exception) { _end_iteration=it; }
    MEDCOUPLING_EXPORT void setStartOrder(int order) throw(INTERP_KERNEL::Exception) { _start_order=order; }
    MEDCOUPLING_EXPORT void setEndOrder(int order) throw(INTERP_KERNEL::Exception) { _end_order=order; }
    MEDCOUPLING_EXPORT void setStartTimeValue(double time) throw(INTERP_KERNEL::Exception) { _start_time=time; }
    MEDCOUPLING_EXPORT void setEndTimeValue(double time) throw(INTERP_KERNEL::Exception) { _end_time=time; }
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void resizeForUnserialization(const std::vector<int>& tinyInfoI, std::vector<DataArrayDouble *>& arrays) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation2(std::vector<int>& tinyInfo) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getTinySerializationDbleInformation2(std::vector<double>& tinyInfo) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void finishUnserialization2(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT std::vector< const DataArrayDouble *> getArraysForTime(double time) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setArrays(const std::vector<DataArrayDouble *>& arrays, TimeLabel *owner) throw(INTERP_KERNEL::Exception);
  protected:
    static const char EXCEPTION_MSG[];
  protected:
    double _start_time;
    double _end_time;
    int _start_iteration;
    int _end_iteration;
    int _start_order;
    int _end_order;
    DataArrayDouble *_end_array;
  };

  class MEDCouplingLinearTime : public MEDCouplingTwoTimeSteps
  {
  protected:
    MEDCOUPLING_EXPORT MEDCouplingLinearTime(const MEDCouplingLinearTime& other, bool deepCpy);
  public:
    MEDCOUPLING_EXPORT MEDCouplingLinearTime();
    MEDCOUPLING_EXPORT std::string getStringRepr() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT TypeOfTimeDiscretization getEnum() const throw(INTERP_KERNEL::Exception) { return DISCRETIZATION; }
    MEDCOUPLING_EXPORT void checkCoherency() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *performCpy(bool deepCpy) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool areCompatible(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool areStrictlyCompatible(const MEDCouplingTimeDiscretization *other, std::string& reason) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool areStrictlyCompatibleForDiv(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool areCompatibleForMeld(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getValueForTime(double time, const std::vector<double>& vals, double *res) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getValueOnTime(int eltId, double time, double *value) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getValueOnDiscTime(int eltId, int iteration, int order, double *value) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *aggregate(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *aggregate(const std::vector<const MEDCouplingTimeDiscretization *>& other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *meld(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *dot(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *crossProduct(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *max(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *min(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *add(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void addEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *substract(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void substractEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *multiply(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void multiplyEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *divide(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void divideEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *pow(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void powEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception);
  public:
    static const TypeOfTimeDiscretization DISCRETIZATION=LINEAR_TIME;
    static const char REPR[];
  };
}

#endif
