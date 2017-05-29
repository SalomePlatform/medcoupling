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

#ifndef __MEDCOUPLINGTIMEDISCRETIZATION_HXX__
#define __MEDCOUPLINGTIMEDISCRETIZATION_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingTimeLabel.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "InterpKernelException.hxx"
#include "MEDCouplingTraits.hxx"

#include <vector>

namespace MEDCoupling
{
  class MEDCouplingMesh;
  class DataArrayDouble;
  class TimeLabel;

  class TimeHolder
  {
  public:
    MEDCOUPLING_EXPORT std::string getTimeUnit() const { return _time_unit; }
    MEDCOUPLING_EXPORT void setTimeUnit(const std::string& unit) { _time_unit=unit; }
    MEDCOUPLING_EXPORT double getTime(int& iteration, int& order) const { return getStartTime(iteration,order); }
    MEDCOUPLING_EXPORT virtual double getStartTime(int& iteration, int& order) const = 0;
    void copyTinyAttrFrom(const TimeHolder& other) { _time_unit=other._time_unit; }
  protected:
    TimeHolder() { }
    TimeHolder(const TimeHolder& other):_time_unit(other._time_unit) { }
    virtual ~TimeHolder() { }
  private:
    std::string _time_unit;
  };

  template<class T>
  class MEDCouplingTimeDiscretizationTemplate : public TimeLabel, public BigMemoryObject, public TimeHolder
  {
  public:
    MEDCOUPLING_EXPORT void updateTime() const;
    MEDCOUPLING_EXPORT virtual void setArray(typename Traits<T>::ArrayType *array, TimeLabel *owner);
    MEDCOUPLING_EXPORT typename Traits<T>::ArrayType *getArray() { return _array; }
    MEDCOUPLING_EXPORT const typename Traits<T>::ArrayType *getArray() const { return _array; }
    MEDCOUPLING_EXPORT void setTimeTolerance(double val) { _time_tolerance=val; }
    MEDCOUPLING_EXPORT double getTimeTolerance() const { return _time_tolerance; }
    MEDCOUPLING_EXPORT void setTime(double time, int iteration, int order) { setStartTime(time,iteration,order); }
    MEDCOUPLING_EXPORT void setIteration(int it) { setStartIteration(it); }
    MEDCOUPLING_EXPORT void setOrder(int order) { setStartOrder(order); }
    MEDCOUPLING_EXPORT void setTimeValue(double val) { setStartTimeValue(val); }
    MEDCOUPLING_EXPORT virtual void setStartIteration(int it) = 0;
    MEDCOUPLING_EXPORT virtual void setStartOrder(int order) = 0;
    MEDCOUPLING_EXPORT virtual void setStartTime(double time, int iteration, int order) = 0;
    MEDCOUPLING_EXPORT virtual void setStartTimeValue(double time) = 0;
    MEDCOUPLING_EXPORT virtual double getEndTime(int& iteration, int& order) const = 0;
    MEDCOUPLING_EXPORT virtual void setEndIteration(int it) = 0;
    MEDCOUPLING_EXPORT virtual void setEndOrder(int order) = 0;
    MEDCOUPLING_EXPORT virtual void setEndTimeValue(double time) = 0;
    MEDCOUPLING_EXPORT virtual void setEndTime(double time, int iteration, int order) = 0;
    
    MEDCOUPLING_EXPORT virtual void setEndArray(typename Traits<T>::ArrayType *array, TimeLabel *owner);
    MEDCOUPLING_EXPORT virtual void setArrays(const std::vector< typename Traits<T>::ArrayType *>& arrays, TimeLabel *owner);
    MEDCOUPLING_EXPORT virtual const typename Traits<T>::ArrayType *getEndArray() const;
    MEDCOUPLING_EXPORT virtual typename Traits<T>::ArrayType *getEndArray();
    MEDCOUPLING_EXPORT virtual void getArrays(std::vector<typename Traits<T>::ArrayType *>& arrays) const;
    
    MEDCOUPLING_EXPORT virtual std::string getStringRepr() const = 0;
    MEDCOUPLING_EXPORT virtual TypeOfTimeDiscretization getEnum() const = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretizationTemplate<T> *performCopyOrIncrRef(bool deepCopy) const = 0;
    //
    MEDCOUPLING_EXPORT virtual void copyTinyAttrFrom(const MEDCouplingTimeDiscretizationTemplate<T>& other);
    MEDCOUPLING_EXPORT virtual void copyTinyStringsFrom(const MEDCouplingTimeDiscretizationTemplate<T>& other);
    MEDCOUPLING_EXPORT virtual void checkConsistencyLight() const;
    MEDCOUPLING_EXPORT virtual std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT virtual std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    //
    MEDCOUPLING_EXPORT virtual bool isEqualIfNotWhy(const MEDCouplingTimeDiscretizationTemplate<T> *other, T prec, std::string& reason) const = 0;
    MEDCOUPLING_EXPORT virtual bool isEqualWithoutConsideringStr(const MEDCouplingTimeDiscretizationTemplate<T> *other, T prec) const = 0;
    //
    MEDCOUPLING_EXPORT virtual bool areStrictlyCompatible(const MEDCouplingTimeDiscretizationTemplate<T> *other, std::string& reason) const;
    MEDCOUPLING_EXPORT virtual bool areCompatible(const MEDCouplingTimeDiscretizationTemplate<T> *other) const;
    MEDCOUPLING_EXPORT virtual bool areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretizationTemplate<T> *other) const;
    MEDCOUPLING_EXPORT virtual bool areStrictlyCompatibleForDiv(const MEDCouplingTimeDiscretizationTemplate<T> *other) const;
    MEDCOUPLING_EXPORT virtual ~MEDCouplingTimeDiscretizationTemplate();
    MEDCOUPLING_EXPORT virtual void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    MEDCOUPLING_EXPORT virtual void getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const;
    MEDCOUPLING_EXPORT virtual void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    MEDCOUPLING_EXPORT virtual void resizeForUnserialization(const std::vector<int>& tinyInfoI, std::vector<typename Traits<T>::ArrayType *>& arrays);
    MEDCOUPLING_EXPORT virtual void checkForUnserialization(const std::vector<int>& tinyInfoI, const std::vector<typename Traits<T>::ArrayType *>& arrays);
    MEDCOUPLING_EXPORT virtual void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS);
  protected:
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretizationTemplate();
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretizationTemplate(const MEDCouplingTimeDiscretizationTemplate<T>& other, bool deepCopy);
  protected:
    double _time_tolerance;
    typename Traits<T>::ArrayType *_array;
  protected:
    static const double TIME_TOLERANCE_DFT;
  };

  class MEDCouplingTimeKeeper
  {
  public:
    MEDCOUPLING_EXPORT MEDCouplingTimeKeeper():_time(0.),_iteration(-1),_order(-1) { }
    MEDCOUPLING_EXPORT double getAllInfo(int& iteration, int& order) const { iteration=_iteration; order=_order; return _time; }
    MEDCOUPLING_EXPORT void setAllInfo(double time, int iteration, int order) { _time=time; _iteration=iteration; _order=order; }
    MEDCOUPLING_EXPORT int getIteration() const { return _iteration; }
    MEDCOUPLING_EXPORT void setIteration(int it) { _iteration=it; }
    MEDCOUPLING_EXPORT int getOrder() const { return _order; }
    MEDCOUPLING_EXPORT void setOrder(int order) { _order=order; }
    MEDCOUPLING_EXPORT double getTimeValue() const { return _time; }
    MEDCOUPLING_EXPORT void setTimeValue(double time) { _time=time; }
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingTimeKeeper& other, double prec, std::string& reason) const;
    MEDCOUPLING_EXPORT bool isEqual(const MEDCouplingTimeKeeper& other, double prec) const;
    MEDCOUPLING_EXPORT void copyFrom(const MEDCouplingTimeKeeper& other);
    MEDCOUPLING_EXPORT void checkTimePresence(double time, double eps) const;
  private:
    double _time;
    int _iteration;
    int _order;
  };

  class MEDCouplingTimeDiscretization : public MEDCouplingTimeDiscretizationTemplate<double>
  {
  protected:
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization() { }
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization(const MEDCouplingTimeDiscretization& other, bool deepCopy):MEDCouplingTimeDiscretizationTemplate<double>(other,deepCopy) { }
  public:
    MEDCOUPLING_EXPORT static MEDCouplingTimeDiscretization *New(TypeOfTimeDiscretization type);
    MEDCOUPLING_EXPORT virtual bool areCompatibleForMeld(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT virtual bool isEqualWithoutConsideringStr(const MEDCouplingTimeDiscretizationTemplate<double> *other, double prec) const;
    MEDCOUPLING_EXPORT virtual bool isEqualIfNotWhy(const MEDCouplingTimeDiscretizationTemplate<double> *other, double prec, std::string& reason) const;
    MEDCOUPLING_EXPORT virtual bool isEqual(const MEDCouplingTimeDiscretization *other, double prec) const;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *buildNewTimeReprFromThis(TypeOfTimeDiscretization type, bool deepCopy) const;
    MEDCOUPLING_EXPORT virtual void synchronizeTimeWith(const MEDCouplingMesh *mesh) = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *aggregate(const MEDCouplingTimeDiscretization *other) const = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *aggregate(const std::vector<const MEDCouplingTimeDiscretization *>& other) const = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *meld(const MEDCouplingTimeDiscretization *other) const = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *dot(const MEDCouplingTimeDiscretization *other) const = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *crossProduct(const MEDCouplingTimeDiscretization *other) const = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *max(const MEDCouplingTimeDiscretization *other) const = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *min(const MEDCouplingTimeDiscretization *other) const = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *add(const MEDCouplingTimeDiscretization *other) const = 0;
    MEDCOUPLING_EXPORT virtual void addEqual(const MEDCouplingTimeDiscretization *other) = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *substract(const MEDCouplingTimeDiscretization *other) const = 0;
    MEDCOUPLING_EXPORT virtual void substractEqual(const MEDCouplingTimeDiscretization *other) = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *multiply(const MEDCouplingTimeDiscretization *other) const = 0;
    MEDCOUPLING_EXPORT virtual void multiplyEqual(const MEDCouplingTimeDiscretization *other) = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *divide(const MEDCouplingTimeDiscretization *other) const = 0;
    MEDCOUPLING_EXPORT virtual void divideEqual(const MEDCouplingTimeDiscretization *other) = 0;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *pow(const MEDCouplingTimeDiscretization *other) const = 0;
    MEDCOUPLING_EXPORT virtual void powEqual(const MEDCouplingTimeDiscretization *other) = 0;
    MEDCOUPLING_EXPORT virtual void getTinySerializationIntInformation2(std::vector<int>& tinyInfo) const = 0;
    MEDCOUPLING_EXPORT virtual void getTinySerializationDbleInformation2(std::vector<double>& tinyInfo) const = 0;
    MEDCOUPLING_EXPORT virtual void finishUnserialization2(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD) = 0;
    MEDCOUPLING_EXPORT virtual void checkNoTimePresence() const = 0;
    MEDCOUPLING_EXPORT virtual void checkTimePresence(double time) const = 0;
    MEDCOUPLING_EXPORT virtual std::vector< const DataArrayDouble *> getArraysForTime(double time) const = 0;
    MEDCOUPLING_EXPORT virtual void getValueForTime(double time, const std::vector<double>& vals, double *res) const = 0;
    MEDCOUPLING_EXPORT virtual bool isBefore(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT virtual bool isStrictlyBefore(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT virtual void getValueOnTime(int eltId, double time, double *value) const = 0;
    MEDCOUPLING_EXPORT virtual void getValueOnDiscTime(int eltId, int iteration, int order, double *value) const = 0;
    //
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *computeVectorFieldCyl(const DataArrayDouble *coords, const double center[3], const double vect[3]) const;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *doublyContractedProduct() const;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *determinant() const;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *eigenValues() const;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *eigenVectors() const;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *inverse() const;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *trace() const;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *deviator() const;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *magnitude() const;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *negate() const;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *maxPerTuple() const;
    MEDCOUPLING_EXPORT virtual MEDCouplingTimeDiscretization *keepSelectedComponents(const std::vector<int>& compoIds) const;
    MEDCOUPLING_EXPORT virtual void setSelectedComponents(const MEDCouplingTimeDiscretization *other, const std::vector<int>& compoIds);
    MEDCOUPLING_EXPORT virtual void changeNbOfComponents(int newNbOfComp, double dftValue);
    MEDCOUPLING_EXPORT virtual void sortPerTuple(bool asc);
    MEDCOUPLING_EXPORT virtual void setUniformValue(int nbOfTuple, int nbOfCompo, double value);
    MEDCOUPLING_EXPORT virtual void setOrCreateUniformValueOnAllComponents(int nbOfTuple, double value);
    MEDCOUPLING_EXPORT virtual void applyLin(double a, double b, int compoId);
    MEDCOUPLING_EXPORT virtual void applyLin(double a, double b);
    MEDCOUPLING_EXPORT virtual void applyFunc(int nbOfComp, FunctionToEvaluate func);
    MEDCOUPLING_EXPORT virtual void applyFunc(int nbOfComp, const std::string& func);
    MEDCOUPLING_EXPORT virtual void applyFuncCompo(int nbOfComp, const std::string& func);
    MEDCOUPLING_EXPORT virtual void applyFuncNamedCompo(int nbOfComp, const std::vector<std::string>& varsOrder, const std::string& func);
    MEDCOUPLING_EXPORT virtual void applyFunc(const std::string& func);
    MEDCOUPLING_EXPORT virtual void applyFuncFast32(const std::string& func);
    MEDCOUPLING_EXPORT virtual void applyFuncFast64(const std::string& func);
    MEDCOUPLING_EXPORT virtual void fillFromAnalytic(const DataArrayDouble *loc, int nbOfComp, FunctionToEvaluate func);
    MEDCOUPLING_EXPORT virtual void fillFromAnalytic(const DataArrayDouble *loc, int nbOfComp, const std::string& func);
    MEDCOUPLING_EXPORT virtual void fillFromAnalyticCompo(const DataArrayDouble *loc, int nbOfComp, const std::string& func);
    MEDCOUPLING_EXPORT virtual void fillFromAnalyticNamedCompo(const DataArrayDouble *loc, int nbOfComp, const std::vector<std::string>& varsOrder, const std::string& func);
  };

  template<class T>
  class MEDCouplingTimeDiscretizationSimple : public MEDCouplingTimeDiscretizationTemplate<T>
  {
  public:
    std::string getStringRepr() const;
    double getStartTime(int& iteration, int& order) const { return _tk.getAllInfo(iteration,order); }
    void setStartIteration(int it) { _tk.setIteration(it); }
    void setStartOrder(int order) { _tk.setOrder(order); }
    void setStartTimeValue(double time) { _tk.setTimeValue(time); } 
    void setStartTime(double time, int iteration, int order) { _tk.setAllInfo(time,iteration,order); }
    double getEndTime(int& iteration, int& order) const;
    void setEndIteration(int it);
    void setEndOrder(int order);
    void setEndTimeValue(double time);
    void setEndTime(double time, int iteration, int order);
    TypeOfTimeDiscretization getEnum() const { return DISCRETIZATION; }
  protected:
    MEDCouplingTimeDiscretizationSimple(const MEDCouplingTimeDiscretizationSimple<T>& other, bool deepCopy);
    MEDCouplingTimeDiscretizationSimple() { }
  protected:
    MEDCouplingTimeKeeper _tk;
  private:
    MEDCOUPLING_EXPORT static const char REPR[];
  protected:
    static const TypeOfTimeDiscretization DISCRETIZATION=ONE_TIME;
  };
  
  class MEDCouplingTimeDiscretizationInt : public MEDCouplingTimeDiscretizationSimple<int>
  {
  public:
    MEDCouplingTimeDiscretizationInt() { }
    MEDCouplingTimeDiscretizationInt(const MEDCouplingTimeDiscretizationInt& other, bool deepCopy);
    static MEDCouplingTimeDiscretizationInt *New(TypeOfTimeDiscretization type);
    MEDCouplingTimeDiscretizationInt *performCopyOrIncrRef(bool deepCopy) const;
    bool isEqualIfNotWhy(const MEDCouplingTimeDiscretizationTemplate<int> *other, int prec, std::string& reason) const;
    bool isEqualWithoutConsideringStr(const MEDCouplingTimeDiscretizationTemplate<int> *other, int prec) const;
  };

  class MEDCouplingTimeDiscretizationFloat : public MEDCouplingTimeDiscretizationSimple<float>
  {
  public:
    MEDCouplingTimeDiscretizationFloat() { }
    MEDCouplingTimeDiscretizationFloat(const MEDCouplingTimeDiscretizationFloat& other, bool deepCopy);
    static MEDCouplingTimeDiscretizationFloat *New(TypeOfTimeDiscretization type);
    MEDCouplingTimeDiscretizationFloat *performCopyOrIncrRef(bool deepCopy) const;
    bool isEqualIfNotWhy(const MEDCouplingTimeDiscretizationTemplate<float> *other, float prec, std::string& reason) const;
    bool isEqualWithoutConsideringStr(const MEDCouplingTimeDiscretizationTemplate<float> *other, float prec) const;
  private:
    static const TypeOfTimeDiscretization DISCRETIZATION=ONE_TIME;
  };
  
  class MEDCouplingNoTimeLabel : public MEDCouplingTimeDiscretization
  {
  public:
    MEDCOUPLING_EXPORT MEDCouplingNoTimeLabel();
    MEDCOUPLING_EXPORT MEDCouplingNoTimeLabel(const MEDCouplingTimeDiscretization& other, bool deepCopy);
    MEDCOUPLING_EXPORT std::string getStringRepr() const;
    MEDCOUPLING_EXPORT TypeOfTimeDiscretization getEnum() const { return DISCRETIZATION; }
    MEDCOUPLING_EXPORT void synchronizeTimeWith(const MEDCouplingMesh *mesh);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *aggregate(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *aggregate(const std::vector<const MEDCouplingTimeDiscretization *>& other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *meld(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *dot(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *crossProduct(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *max(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *min(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *add(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT void addEqual(const MEDCouplingTimeDiscretization *other);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *substract(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT void substractEqual(const MEDCouplingTimeDiscretization *other);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *multiply(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT void multiplyEqual(const MEDCouplingTimeDiscretization *other);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *divide(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT void divideEqual(const MEDCouplingTimeDiscretization *other);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *pow(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT void powEqual(const MEDCouplingTimeDiscretization *other);
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingTimeDiscretizationTemplate<double> *other, double prec, std::string& reason) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const MEDCouplingTimeDiscretizationTemplate<double> *other, double prec) const;
    MEDCOUPLING_EXPORT bool areCompatible(const MEDCouplingTimeDiscretizationTemplate<double> *other) const;
    MEDCOUPLING_EXPORT bool areStrictlyCompatible(const MEDCouplingTimeDiscretizationTemplate<double> *other, std::string& reason) const;
    MEDCOUPLING_EXPORT bool areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretizationTemplate<double> *other) const;
    MEDCOUPLING_EXPORT bool areStrictlyCompatibleForDiv(const MEDCouplingTimeDiscretizationTemplate<double> *other) const;
    MEDCOUPLING_EXPORT bool areCompatibleForMeld(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *performCopyOrIncrRef(bool deepCopy) const;
    MEDCOUPLING_EXPORT void checkNoTimePresence() const { }
    MEDCOUPLING_EXPORT void checkTimePresence(double time) const;
    MEDCOUPLING_EXPORT std::vector< const DataArrayDouble *> getArraysForTime(double time) const;
    MEDCOUPLING_EXPORT void getValueForTime(double time, const std::vector<double>& vals, double *res) const;
    MEDCOUPLING_EXPORT bool isBefore(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT bool isStrictlyBefore(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT double getStartTime(int& iteration, int& order) const;
    MEDCOUPLING_EXPORT double getEndTime(int& iteration, int& order) const;
    MEDCOUPLING_EXPORT void setStartIteration(int it);
    MEDCOUPLING_EXPORT void setEndIteration(int it);
    MEDCOUPLING_EXPORT void setStartOrder(int order);
    MEDCOUPLING_EXPORT void setEndOrder(int order);
    MEDCOUPLING_EXPORT void setStartTimeValue(double time);
    MEDCOUPLING_EXPORT void setEndTimeValue(double time);
    MEDCOUPLING_EXPORT void setStartTime(double time, int iteration, int order);
    MEDCOUPLING_EXPORT void setEndTime(double time, int iteration, int order);
    MEDCOUPLING_EXPORT void getValueOnTime(int eltId, double time, double *value) const;
    MEDCOUPLING_EXPORT void getValueOnDiscTime(int eltId, int iteration, int order, double *value) const;
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation2(std::vector<int>& tinyInfo) const;
    MEDCOUPLING_EXPORT void getTinySerializationDbleInformation2(std::vector<double>& tinyInfo) const;
    MEDCOUPLING_EXPORT void finishUnserialization2(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD);
  public:
    static const TypeOfTimeDiscretization DISCRETIZATION=NO_TIME;
    MEDCOUPLING_EXPORT static const char REPR[];
  private:
    static const char EXCEPTION_MSG[];
  };

  class MEDCouplingWithTimeStep : public MEDCouplingTimeDiscretization
  {
  protected:
    MEDCOUPLING_EXPORT MEDCouplingWithTimeStep(const MEDCouplingWithTimeStep& other, bool deepCopy);
  public:
    MEDCOUPLING_EXPORT MEDCouplingWithTimeStep();
    MEDCOUPLING_EXPORT std::string getStringRepr() const;
    MEDCOUPLING_EXPORT void copyTinyAttrFrom(const MEDCouplingTimeDiscretizationTemplate<double>& other);
    MEDCOUPLING_EXPORT TypeOfTimeDiscretization getEnum() const { return DISCRETIZATION; }
    MEDCOUPLING_EXPORT void synchronizeTimeWith(const MEDCouplingMesh *mesh);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *aggregate(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *aggregate(const std::vector<const MEDCouplingTimeDiscretization *>& other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *meld(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *dot(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *crossProduct(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *max(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *min(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *add(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT void addEqual(const MEDCouplingTimeDiscretization *other);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *substract(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT void substractEqual(const MEDCouplingTimeDiscretization *other);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *multiply(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT void multiplyEqual(const MEDCouplingTimeDiscretization *other);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *divide(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT void divideEqual(const MEDCouplingTimeDiscretization *other);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *pow(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT void powEqual(const MEDCouplingTimeDiscretization *other);
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingTimeDiscretizationTemplate<double> *other, double prec, std::string& reason) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const MEDCouplingTimeDiscretizationTemplate<double> *other, double prec) const;
    MEDCOUPLING_EXPORT bool areCompatible(const MEDCouplingTimeDiscretizationTemplate<double> *other) const;
    MEDCOUPLING_EXPORT bool areStrictlyCompatible(const MEDCouplingTimeDiscretizationTemplate<double> *other, std::string& reason) const;
    MEDCOUPLING_EXPORT bool areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretizationTemplate<double> *other) const;
    MEDCOUPLING_EXPORT bool areStrictlyCompatibleForDiv(const MEDCouplingTimeDiscretizationTemplate<double> *other) const;
    MEDCOUPLING_EXPORT bool areCompatibleForMeld(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    MEDCOUPLING_EXPORT void getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const;
    MEDCOUPLING_EXPORT void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS);
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation2(std::vector<int>& tinyInfo) const;
    MEDCOUPLING_EXPORT void getTinySerializationDbleInformation2(std::vector<double>& tinyInfo) const;
    MEDCOUPLING_EXPORT void finishUnserialization2(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *performCopyOrIncrRef(bool deepCopy) const;
    MEDCOUPLING_EXPORT void checkNoTimePresence() const;
    MEDCOUPLING_EXPORT void checkTimePresence(double time) const;
    MEDCOUPLING_EXPORT void setStartTime(double time, int iteration, int order) { _tk.setAllInfo(time,iteration,order); }
    MEDCOUPLING_EXPORT void setEndTime(double time, int iteration, int order) { _tk.setAllInfo(time,iteration,order); }
    MEDCOUPLING_EXPORT double getStartTime(int& iteration, int& order) const { return _tk.getAllInfo(iteration,order); }
    MEDCOUPLING_EXPORT double getEndTime(int& iteration, int& order) const { return _tk.getAllInfo(iteration,order); }
    MEDCOUPLING_EXPORT void setStartIteration(int it) { _tk.setIteration(it); }
    MEDCOUPLING_EXPORT void setEndIteration(int it) { _tk.setIteration(it); }
    MEDCOUPLING_EXPORT void setStartOrder(int order) { _tk.setOrder(order); }
    MEDCOUPLING_EXPORT void setEndOrder(int order) { _tk.setOrder(order); }
    MEDCOUPLING_EXPORT void setStartTimeValue(double time) { _tk.setTimeValue(time); }
    MEDCOUPLING_EXPORT void setEndTimeValue(double time) { _tk.setTimeValue(time); }
    MEDCOUPLING_EXPORT std::vector< const DataArrayDouble *> getArraysForTime(double time) const;
    MEDCOUPLING_EXPORT void getValueForTime(double time, const std::vector<double>& vals, double *res) const;
    MEDCOUPLING_EXPORT void getValueOnTime(int eltId, double time, double *value) const;
    MEDCOUPLING_EXPORT void getValueOnDiscTime(int eltId, int iteration, int order, double *value) const;
  public:
    static const TypeOfTimeDiscretization DISCRETIZATION=ONE_TIME;
    MEDCOUPLING_EXPORT static const char REPR[];
  private:
    static const char EXCEPTION_MSG[];
  protected:
    MEDCouplingTimeKeeper _tk;
  };

  class MEDCouplingConstOnTimeInterval : public MEDCouplingTimeDiscretization
  {
  protected:
    MEDCOUPLING_EXPORT MEDCouplingConstOnTimeInterval(const MEDCouplingConstOnTimeInterval& other, bool deepCopy);
  public:
    MEDCOUPLING_EXPORT MEDCouplingConstOnTimeInterval();
    MEDCOUPLING_EXPORT void copyTinyAttrFrom(const MEDCouplingTimeDiscretizationTemplate<double>& other);
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    MEDCOUPLING_EXPORT void getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const;
    MEDCOUPLING_EXPORT void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS);
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation2(std::vector<int>& tinyInfo) const;
    MEDCOUPLING_EXPORT void getTinySerializationDbleInformation2(std::vector<double>& tinyInfo) const;
    MEDCOUPLING_EXPORT void finishUnserialization2(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *performCopyOrIncrRef(bool deepCopy) const;
    MEDCOUPLING_EXPORT bool areCompatible(const MEDCouplingTimeDiscretizationTemplate<double> *other) const;
    MEDCOUPLING_EXPORT bool areStrictlyCompatible(const MEDCouplingTimeDiscretizationTemplate<double> *other, std::string& reason) const;
    MEDCOUPLING_EXPORT bool areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretizationTemplate<double> *other) const;
    MEDCOUPLING_EXPORT bool areStrictlyCompatibleForDiv(const MEDCouplingTimeDiscretizationTemplate<double> *other) const;
    MEDCOUPLING_EXPORT bool areCompatibleForMeld(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingTimeDiscretizationTemplate<double> *other, double prec, std::string& reason) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const MEDCouplingTimeDiscretizationTemplate<double> *other, double prec) const;
    MEDCOUPLING_EXPORT std::vector< const DataArrayDouble *> getArraysForTime(double time) const;
    MEDCOUPLING_EXPORT void getValueForTime(double time, const std::vector<double>& vals, double *res) const;
    MEDCOUPLING_EXPORT void getValueOnTime(int eltId, double time, double *value) const;
    MEDCOUPLING_EXPORT void getValueOnDiscTime(int eltId, int iteration, int order, double *value) const;
    MEDCOUPLING_EXPORT TypeOfTimeDiscretization getEnum() const { return DISCRETIZATION; }
    MEDCOUPLING_EXPORT void synchronizeTimeWith(const MEDCouplingMesh *mesh);
    MEDCOUPLING_EXPORT std::string getStringRepr() const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *aggregate(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *aggregate(const std::vector<const MEDCouplingTimeDiscretization *>& other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *meld(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *dot(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *crossProduct(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *max(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *min(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *add(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT void addEqual(const MEDCouplingTimeDiscretization *other);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *substract(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT void substractEqual(const MEDCouplingTimeDiscretization *other);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *multiply(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT void multiplyEqual(const MEDCouplingTimeDiscretization *other);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *divide(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT void divideEqual(const MEDCouplingTimeDiscretization *other);
    MEDCouplingTimeDiscretization *pow(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT void powEqual(const MEDCouplingTimeDiscretization *other);
    MEDCOUPLING_EXPORT void setStartTime(double time, int iteration, int order) { _start.setAllInfo(time,iteration,order); }
    MEDCOUPLING_EXPORT void setEndTime(double time, int iteration, int order) { _end.setAllInfo(time,iteration,order); }
    MEDCOUPLING_EXPORT double getStartTime(int& iteration, int& order) const { return _start.getAllInfo(iteration,order); }
    MEDCOUPLING_EXPORT double getEndTime(int& iteration, int& order) const { return _end.getAllInfo(iteration,order); }
    MEDCOUPLING_EXPORT void setStartIteration(int it) { _start.setIteration(it); }
    MEDCOUPLING_EXPORT void setEndIteration(int it) { _end.setIteration(it); }
    MEDCOUPLING_EXPORT void setStartOrder(int order) { _start.setOrder(order); }
    MEDCOUPLING_EXPORT void setEndOrder(int order) { _end.setOrder(order); }
    MEDCOUPLING_EXPORT void setStartTimeValue(double time) { _start.setTimeValue(time); }
    MEDCOUPLING_EXPORT void setEndTimeValue(double time) { _end.setTimeValue(time); }
    MEDCOUPLING_EXPORT void checkNoTimePresence() const;
    MEDCOUPLING_EXPORT void checkTimePresence(double time) const;
  public:
    static const TypeOfTimeDiscretization DISCRETIZATION=CONST_ON_TIME_INTERVAL;
    MEDCOUPLING_EXPORT static const char REPR[];
  private:
    static const char EXCEPTION_MSG[];
  protected:
    MEDCouplingTimeKeeper _start;
    MEDCouplingTimeKeeper _end;
  };

  class MEDCouplingTwoTimeSteps : public MEDCouplingTimeDiscretization
  {
  protected:
    MEDCOUPLING_EXPORT MEDCouplingTwoTimeSteps(const MEDCouplingTwoTimeSteps& other, bool deepCopy);
    MEDCOUPLING_EXPORT MEDCouplingTwoTimeSteps();
    MEDCOUPLING_EXPORT ~MEDCouplingTwoTimeSteps();
  public:
    MEDCOUPLING_EXPORT void updateTime() const;
    MEDCOUPLING_EXPORT void synchronizeTimeWith(const MEDCouplingMesh *mesh);
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDCOUPLING_EXPORT void copyTinyAttrFrom(const MEDCouplingTimeDiscretizationTemplate<double>& other);
    MEDCOUPLING_EXPORT void copyTinyStringsFrom(const MEDCouplingTimeDiscretizationTemplate<double>& other);
    MEDCOUPLING_EXPORT const DataArrayDouble *getEndArray() const;
    MEDCOUPLING_EXPORT DataArrayDouble *getEndArray();
    MEDCOUPLING_EXPORT void checkConsistencyLight() const;
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingTimeDiscretizationTemplate<double> *other, double prec, std::string& reason) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const MEDCouplingTimeDiscretizationTemplate<double> *other, double prec) const;
    MEDCOUPLING_EXPORT void checkNoTimePresence() const;
    MEDCOUPLING_EXPORT void checkTimePresence(double time) const;
    MEDCOUPLING_EXPORT void getArrays(std::vector<DataArrayDouble *>& arrays) const;
    MEDCOUPLING_EXPORT void setEndArray(DataArrayDouble *array, TimeLabel *owner);
    MEDCOUPLING_EXPORT void setStartTime(double time, int iteration, int order) { _start.setAllInfo(time,iteration,order); }
    MEDCOUPLING_EXPORT void setEndTime(double time, int iteration, int order) { _end.setAllInfo(time,iteration,order); }
    MEDCOUPLING_EXPORT double getStartTime(int& iteration, int& order) const { return _start.getAllInfo(iteration,order); }
    MEDCOUPLING_EXPORT double getEndTime(int& iteration, int& order) const { return _end.getAllInfo(iteration,order); }
    MEDCOUPLING_EXPORT void setStartIteration(int it) { _start.setIteration(it); }
    MEDCOUPLING_EXPORT void setEndIteration(int it) { _end.setIteration(it); }
    MEDCOUPLING_EXPORT void setStartOrder(int order) { _start.setOrder(order); }
    MEDCOUPLING_EXPORT void setEndOrder(int order) { _end.setOrder(order); }
    MEDCOUPLING_EXPORT void setStartTimeValue(double time) { _start.setTimeValue(time); }
    MEDCOUPLING_EXPORT void setEndTimeValue(double time) { _end.setTimeValue(time); }
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    MEDCOUPLING_EXPORT void getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const;
    MEDCOUPLING_EXPORT void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    MEDCOUPLING_EXPORT void resizeForUnserialization(const std::vector<int>& tinyInfoI, std::vector<DataArrayDouble *>& arrays);
    MEDCOUPLING_EXPORT void checkForUnserialization(const std::vector<int>& tinyInfoI, const std::vector<DataArrayDouble *>& arrays);
    MEDCOUPLING_EXPORT void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS);
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation2(std::vector<int>& tinyInfo) const;
    MEDCOUPLING_EXPORT void getTinySerializationDbleInformation2(std::vector<double>& tinyInfo) const;
    MEDCOUPLING_EXPORT void finishUnserialization2(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD);
    MEDCOUPLING_EXPORT std::vector< const DataArrayDouble *> getArraysForTime(double time) const;
    MEDCOUPLING_EXPORT void setArrays(const std::vector<DataArrayDouble *>& arrays, TimeLabel *owner);
  protected:
    static const char EXCEPTION_MSG[];
  protected:
    MEDCouplingTimeKeeper _start;
    MEDCouplingTimeKeeper _end;
    DataArrayDouble *_end_array;
  };

  class MEDCouplingLinearTime : public MEDCouplingTwoTimeSteps
  {
  protected:
    MEDCOUPLING_EXPORT MEDCouplingLinearTime(const MEDCouplingLinearTime& other, bool deepCopy);
  public:
    MEDCOUPLING_EXPORT MEDCouplingLinearTime();
    MEDCOUPLING_EXPORT std::string getStringRepr() const;
    MEDCOUPLING_EXPORT TypeOfTimeDiscretization getEnum() const { return DISCRETIZATION; }
    MEDCOUPLING_EXPORT void checkConsistencyLight() const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *performCopyOrIncrRef(bool deepCopy) const;
    MEDCOUPLING_EXPORT bool areCompatible(const MEDCouplingTimeDiscretizationTemplate<double> *other) const;
    MEDCOUPLING_EXPORT bool areStrictlyCompatible(const MEDCouplingTimeDiscretizationTemplate<double> *other, std::string& reason) const;
    MEDCOUPLING_EXPORT bool areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretizationTemplate<double> *other) const;
    MEDCOUPLING_EXPORT bool areStrictlyCompatibleForDiv(const MEDCouplingTimeDiscretizationTemplate<double> *other) const;
    MEDCOUPLING_EXPORT bool areCompatibleForMeld(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT void getValueForTime(double time, const std::vector<double>& vals, double *res) const;
    MEDCOUPLING_EXPORT void getValueOnTime(int eltId, double time, double *value) const;
    MEDCOUPLING_EXPORT void getValueOnDiscTime(int eltId, int iteration, int order, double *value) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *aggregate(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *aggregate(const std::vector<const MEDCouplingTimeDiscretization *>& other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *meld(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *dot(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *crossProduct(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *max(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *min(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *add(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT void addEqual(const MEDCouplingTimeDiscretization *other);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *substract(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT void substractEqual(const MEDCouplingTimeDiscretization *other);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *multiply(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT void multiplyEqual(const MEDCouplingTimeDiscretization *other);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *divide(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT void divideEqual(const MEDCouplingTimeDiscretization *other);
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *pow(const MEDCouplingTimeDiscretization *other) const;
    MEDCOUPLING_EXPORT void powEqual(const MEDCouplingTimeDiscretization *other);
  public:
    static const TypeOfTimeDiscretization DISCRETIZATION=LINEAR_TIME;
    MEDCOUPLING_EXPORT static const char REPR[];
  };
}

#endif
