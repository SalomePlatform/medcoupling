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

#ifndef __PARAMEDMEM_MEDCOUPLINGTIMEDISCRETIZATION_HXX__
#define __PARAMEDMEM_MEDCOUPLINGTIMEDISCRETIZATION_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingTimeLabel.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "InterpKernelException.hxx"

#include <vector>

namespace ParaMEDMEM
{
  class DataArrayDouble;
  class TimeLabel;

  class MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization : public TimeLabel
  {
  protected:
    MEDCouplingTimeDiscretization();
    MEDCouplingTimeDiscretization(const MEDCouplingTimeDiscretization& other, bool deepCpy);
  public:
    void updateTime();
    static MEDCouplingTimeDiscretization *New(TypeOfTimeDiscretization type);
    virtual void copyTinyAttrFrom(const MEDCouplingTimeDiscretization& other);
    virtual void copyTinyStringsFrom(const MEDCouplingTimeDiscretization& other);
    virtual void checkCoherency() const throw(INTERP_KERNEL::Exception);
    virtual bool areCompatible(const MEDCouplingTimeDiscretization *other) const;
    virtual bool areStrictlyCompatible(const MEDCouplingTimeDiscretization *other) const;
    virtual bool areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretization *other) const;
    virtual bool isEqual(const MEDCouplingTimeDiscretization *other, double prec) const;
    virtual bool isEqualWithoutConsideringStr(const MEDCouplingTimeDiscretization *other, double prec) const;
    virtual MEDCouplingTimeDiscretization *buildNewTimeReprFromThis(const MEDCouplingTimeDiscretization *other,
                                                                    TypeOfTimeDiscretization type, bool deepCpy) const;
    virtual std::string getStringRepr() const = 0;
    virtual TypeOfTimeDiscretization getEnum() const = 0;
    virtual MEDCouplingTimeDiscretization *aggregate(const MEDCouplingTimeDiscretization *other) const = 0;
    virtual MEDCouplingTimeDiscretization *dot(const MEDCouplingTimeDiscretization *other) const = 0;
    virtual MEDCouplingTimeDiscretization *crossProduct(const MEDCouplingTimeDiscretization *other) const = 0;
    virtual MEDCouplingTimeDiscretization *max(const MEDCouplingTimeDiscretization *other) const = 0;
    virtual MEDCouplingTimeDiscretization *min(const MEDCouplingTimeDiscretization *other) const = 0;
    virtual MEDCouplingTimeDiscretization *add(const MEDCouplingTimeDiscretization *other) const = 0;
    virtual void addEqual(const MEDCouplingTimeDiscretization *other) = 0;
    virtual MEDCouplingTimeDiscretization *substract(const MEDCouplingTimeDiscretization *other) const = 0;
    virtual void substractEqual(const MEDCouplingTimeDiscretization *other) = 0;
    virtual MEDCouplingTimeDiscretization *multiply(const MEDCouplingTimeDiscretization *other) const = 0;
    virtual void multiplyEqual(const MEDCouplingTimeDiscretization *other) = 0;
    virtual MEDCouplingTimeDiscretization *divide(const MEDCouplingTimeDiscretization *other) const = 0;
    virtual void divideEqual(const MEDCouplingTimeDiscretization *other) = 0;
    virtual void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    virtual void getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const;
    virtual void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    virtual void resizeForUnserialization(const std::vector<int>& tinyInfoI, std::vector<DataArrayDouble *>& arrays);
    virtual void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS);
    virtual MEDCouplingTimeDiscretization *performCpy(bool deepCpy) const = 0;
    void setTimeTolerance(double val) { _time_tolerance=val; }
    double getTimeTolerance() const { return _time_tolerance; }
    virtual void checkNoTimePresence() const throw(INTERP_KERNEL::Exception) = 0;
    virtual void checkTimePresence(double time) const throw(INTERP_KERNEL::Exception) = 0;
    virtual void setArray(DataArrayDouble *array, TimeLabel *owner);
    virtual void setEndArray(DataArrayDouble *array, TimeLabel *owner);
    virtual void setArrays(const std::vector<DataArrayDouble *>& arrays, TimeLabel *owner) throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getArray() const { return _array; }
    virtual DataArrayDouble *getEndArray() const;
    virtual std::vector< const DataArrayDouble *> getArraysForTime(double time) const throw(INTERP_KERNEL::Exception) = 0;
    virtual void getValueForTime(double time, const std::vector<double>& vals, double *res) const = 0; 
    virtual void getArrays(std::vector<DataArrayDouble *>& arrays) const;
    virtual bool isBefore(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    virtual bool isStrictlyBefore(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    double getTime(int& iteration, int& order) const throw(INTERP_KERNEL::Exception) { return getStartTime(iteration,order); }
    virtual double getStartTime(int& iteration, int& order) const throw(INTERP_KERNEL::Exception) = 0;
    virtual double getEndTime(int& iteration, int& order) const throw(INTERP_KERNEL::Exception) = 0;
    void setTime(double time, int iteration, int order) throw(INTERP_KERNEL::Exception) { setStartTime(time,iteration,order); }
    virtual void setStartTime(double time, int iteration, int order) throw(INTERP_KERNEL::Exception) = 0;
    virtual void setEndTime(double time, int iteration, int order) throw(INTERP_KERNEL::Exception) = 0;
    virtual void getValueOnTime(int eltId, double time, double *value) const throw(INTERP_KERNEL::Exception) = 0;
    virtual void getValueOnDiscTime(int eltId, int iteration, int order, double *value) const throw(INTERP_KERNEL::Exception) = 0;
    //
    virtual MEDCouplingTimeDiscretization *doublyContractedProduct() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *determinant() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *eigenValues() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *eigenVectors() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *inverse() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *trace() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *deviator() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *magnitude() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *maxPerTuple() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *keepSelectedComponents(const std::vector<int>& compoIds) const throw(INTERP_KERNEL::Exception);
    virtual void setSelectedComponents(const MEDCouplingTimeDiscretization *other, const std::vector<int>& compoIds) throw(INTERP_KERNEL::Exception);
    virtual void changeNbOfComponents(int newNbOfComp, double dftValue) throw(INTERP_KERNEL::Exception);
    virtual void sortPerTuple(bool asc) throw(INTERP_KERNEL::Exception);
    virtual void setUniformValue(int nbOfTuple, int nbOfCompo, double value);
    virtual void applyLin(double a, double b, int compoId);
    virtual void applyFunc(int nbOfComp, FunctionToEvaluate func);
    virtual void applyFunc(int nbOfComp, const char *func);
    virtual void applyFunc(const char *func);
    virtual void applyFuncFast32(const char *func);
    virtual void applyFuncFast64(const char *func);
    virtual void fillFromAnalytic(const DataArrayDouble *loc, int nbOfComp, FunctionToEvaluate func) throw(INTERP_KERNEL::Exception);
    virtual void fillFromAnalytic(const DataArrayDouble *loc, int nbOfComp, const char *func) throw(INTERP_KERNEL::Exception);
    //
    virtual ~MEDCouplingTimeDiscretization();
  protected:
    double _time_tolerance;
    DataArrayDouble *_array;
  protected:
    static const double TIME_TOLERANCE_DFT;
  };

  class MEDCOUPLING_EXPORT MEDCouplingNoTimeLabel : public MEDCouplingTimeDiscretization
  {
  public:
    MEDCouplingNoTimeLabel();
    MEDCouplingNoTimeLabel(const MEDCouplingTimeDiscretization& other, bool deepCpy);
    std::string getStringRepr() const;
    TypeOfTimeDiscretization getEnum() const { return DISCRETIZATION; }
    MEDCouplingTimeDiscretization *aggregate(const MEDCouplingTimeDiscretization *other) const;
    MEDCouplingTimeDiscretization *dot(const MEDCouplingTimeDiscretization *other) const;
    MEDCouplingTimeDiscretization *crossProduct(const MEDCouplingTimeDiscretization *other) const;
    MEDCouplingTimeDiscretization *max(const MEDCouplingTimeDiscretization *other) const;
    MEDCouplingTimeDiscretization *min(const MEDCouplingTimeDiscretization *other) const;
    MEDCouplingTimeDiscretization *add(const MEDCouplingTimeDiscretization *other) const;
    void addEqual(const MEDCouplingTimeDiscretization *other);
    MEDCouplingTimeDiscretization *substract(const MEDCouplingTimeDiscretization *other) const;
    void substractEqual(const MEDCouplingTimeDiscretization *other);
    MEDCouplingTimeDiscretization *multiply(const MEDCouplingTimeDiscretization *other) const;
    void multiplyEqual(const MEDCouplingTimeDiscretization *other);
    MEDCouplingTimeDiscretization *divide(const MEDCouplingTimeDiscretization *other) const;
    void divideEqual(const MEDCouplingTimeDiscretization *other);
    bool isEqual(const MEDCouplingTimeDiscretization *other, double prec) const;
    bool isEqualWithoutConsideringStr(const MEDCouplingTimeDiscretization *other, double prec) const;
    bool areCompatible(const MEDCouplingTimeDiscretization *other) const;
    bool areStrictlyCompatible(const MEDCouplingTimeDiscretization *other) const;
    bool areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretization *other) const;
    MEDCouplingTimeDiscretization *performCpy(bool deepCpy) const;
    void checkNoTimePresence() const throw(INTERP_KERNEL::Exception) { }
    void checkTimePresence(double time) const throw(INTERP_KERNEL::Exception);
    std::vector< const DataArrayDouble *> getArraysForTime(double time) const throw(INTERP_KERNEL::Exception);
    void getValueForTime(double time, const std::vector<double>& vals, double *res) const;
    bool isBefore(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    bool isStrictlyBefore(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    double getStartTime(int& iteration, int& order) const throw(INTERP_KERNEL::Exception);
    double getEndTime(int& iteration, int& order) const throw(INTERP_KERNEL::Exception);
    void setStartTime(double time, int iteration, int order) throw(INTERP_KERNEL::Exception);
    void setEndTime(double time, int iteration, int order) throw(INTERP_KERNEL::Exception);
    void getValueOnTime(int eltId, double time, double *value) const throw(INTERP_KERNEL::Exception);
    void getValueOnDiscTime(int eltId, int iteration, int order, double *value) const throw(INTERP_KERNEL::Exception);
  public:
    static const TypeOfTimeDiscretization DISCRETIZATION=NO_TIME;
    static const char REPR[];
  private:
    static const char EXCEPTION_MSG[];
  };

  class MEDCOUPLING_EXPORT MEDCouplingWithTimeStep : public MEDCouplingTimeDiscretization
  {
  protected:
    MEDCouplingWithTimeStep(const MEDCouplingWithTimeStep& other, bool deepCpy);
  public:
    MEDCouplingWithTimeStep();
    std::string getStringRepr() const;
    void copyTinyAttrFrom(const MEDCouplingTimeDiscretization& other);
    TypeOfTimeDiscretization getEnum() const { return DISCRETIZATION; }
    MEDCouplingTimeDiscretization *aggregate(const MEDCouplingTimeDiscretization *other) const;
    MEDCouplingTimeDiscretization *dot(const MEDCouplingTimeDiscretization *other) const;
    MEDCouplingTimeDiscretization *crossProduct(const MEDCouplingTimeDiscretization *other) const;
    MEDCouplingTimeDiscretization *max(const MEDCouplingTimeDiscretization *other) const;
    MEDCouplingTimeDiscretization *min(const MEDCouplingTimeDiscretization *other) const;
    MEDCouplingTimeDiscretization *add(const MEDCouplingTimeDiscretization *other) const;
    void addEqual(const MEDCouplingTimeDiscretization *other);
    MEDCouplingTimeDiscretization *substract(const MEDCouplingTimeDiscretization *other) const;
    void substractEqual(const MEDCouplingTimeDiscretization *other);
    MEDCouplingTimeDiscretization *multiply(const MEDCouplingTimeDiscretization *other) const;
    void multiplyEqual(const MEDCouplingTimeDiscretization *other);
    MEDCouplingTimeDiscretization *divide(const MEDCouplingTimeDiscretization *other) const;
    void divideEqual(const MEDCouplingTimeDiscretization *other);
    bool isEqual(const MEDCouplingTimeDiscretization *other, double prec) const;
    bool isEqualWithoutConsideringStr(const MEDCouplingTimeDiscretization *other, double prec) const;
    bool areCompatible(const MEDCouplingTimeDiscretization *other) const;
    bool areStrictlyCompatible(const MEDCouplingTimeDiscretization *other) const;
    bool areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretization *other) const;
    void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    void getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const;
    void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS);
    MEDCouplingTimeDiscretization *performCpy(bool deepCpy) const;
    void checkNoTimePresence() const throw(INTERP_KERNEL::Exception);
    void checkTimePresence(double time) const throw(INTERP_KERNEL::Exception);
    void setStartTime(double time, int iteration, int order) throw(INTERP_KERNEL::Exception) { _time=time; _iteration=iteration; _order=order; }
    void setEndTime(double time, int iteration, int order) throw(INTERP_KERNEL::Exception) { _time=time; _iteration=iteration; _order=order; }
    double getStartTime(int& iteration, int& order) const throw(INTERP_KERNEL::Exception) { iteration=_iteration; order=_order; return _time; }
    double getEndTime(int& iteration, int& order) const throw(INTERP_KERNEL::Exception) { iteration=_iteration; order=_order; return _time; }
    std::vector< const DataArrayDouble *> getArraysForTime(double time) const throw(INTERP_KERNEL::Exception);
    void getValueForTime(double time, const std::vector<double>& vals, double *res) const;
    void getValueOnTime(int eltId, double time, double *value) const throw(INTERP_KERNEL::Exception);
    void getValueOnDiscTime(int eltId, int iteration, int order, double *value) const throw(INTERP_KERNEL::Exception);
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

  class MEDCOUPLING_EXPORT MEDCouplingConstOnTimeInterval : public MEDCouplingTimeDiscretization
  {
  protected:
    MEDCouplingConstOnTimeInterval(const MEDCouplingConstOnTimeInterval& other, bool deepCpy);
  public:
    MEDCouplingConstOnTimeInterval();
    void copyTinyAttrFrom(const MEDCouplingTimeDiscretization& other);
    void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    void getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const;
    void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS);
    MEDCouplingTimeDiscretization *performCpy(bool deepCpy) const;
    bool areCompatible(const MEDCouplingTimeDiscretization *other) const;
    bool areStrictlyCompatible(const MEDCouplingTimeDiscretization *other) const;
    bool areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretization *other) const;
    bool isEqual(const MEDCouplingTimeDiscretization *other, double prec) const;
    bool isEqualWithoutConsideringStr(const MEDCouplingTimeDiscretization *other, double prec) const;
    std::vector< const DataArrayDouble *> getArraysForTime(double time) const throw(INTERP_KERNEL::Exception);
    void getValueForTime(double time, const std::vector<double>& vals, double *res) const;
    void getValueOnTime(int eltId, double time, double *value) const throw(INTERP_KERNEL::Exception);
    void getValueOnDiscTime(int eltId, int iteration, int order, double *value) const throw(INTERP_KERNEL::Exception);
    TypeOfTimeDiscretization getEnum() const { return DISCRETIZATION; }
    std::string getStringRepr() const;
    MEDCouplingTimeDiscretization *aggregate(const MEDCouplingTimeDiscretization *other) const;
    MEDCouplingTimeDiscretization *dot(const MEDCouplingTimeDiscretization *other) const;
    MEDCouplingTimeDiscretization *crossProduct(const MEDCouplingTimeDiscretization *other) const;
    MEDCouplingTimeDiscretization *max(const MEDCouplingTimeDiscretization *other) const;
    MEDCouplingTimeDiscretization *min(const MEDCouplingTimeDiscretization *other) const;
    MEDCouplingTimeDiscretization *add(const MEDCouplingTimeDiscretization *other) const;
    void addEqual(const MEDCouplingTimeDiscretization *other);
    MEDCouplingTimeDiscretization *substract(const MEDCouplingTimeDiscretization *other) const;
    void substractEqual(const MEDCouplingTimeDiscretization *other);
    MEDCouplingTimeDiscretization *multiply(const MEDCouplingTimeDiscretization *other) const;
    void multiplyEqual(const MEDCouplingTimeDiscretization *other);
    MEDCouplingTimeDiscretization *divide(const MEDCouplingTimeDiscretization *other) const;
    void divideEqual(const MEDCouplingTimeDiscretization *other);
    void setStartTime(double time, int iteration, int order) throw(INTERP_KERNEL::Exception) { _start_time=time; _start_iteration=iteration; _start_order=order; }
    void setEndTime(double time, int iteration, int order) throw(INTERP_KERNEL::Exception) { _end_time=time; _end_iteration=iteration; _end_order=order; }
    double getStartTime(int& iteration, int& order) const throw(INTERP_KERNEL::Exception) { iteration=_start_iteration; order=_start_order; return _start_time; }
    double getEndTime(int& iteration, int& order) const throw(INTERP_KERNEL::Exception) { iteration=_end_iteration; order=_end_order; return _end_time; }
    void checkNoTimePresence() const throw(INTERP_KERNEL::Exception);
    void checkTimePresence(double time) const throw(INTERP_KERNEL::Exception);
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

  class MEDCOUPLING_EXPORT MEDCouplingTwoTimeSteps : public MEDCouplingTimeDiscretization
  {
  protected:
    MEDCouplingTwoTimeSteps(const MEDCouplingTwoTimeSteps& other, bool deepCpy);
    MEDCouplingTwoTimeSteps();
    ~MEDCouplingTwoTimeSteps();
  public:
    void updateTime();
    void copyTinyAttrFrom(const MEDCouplingTimeDiscretization& other);
    void copyTinyStringsFrom(const MEDCouplingTimeDiscretization& other);
    DataArrayDouble *getEndArray() const;
    void checkCoherency() const throw(INTERP_KERNEL::Exception);
    bool isEqual(const MEDCouplingTimeDiscretization *other, double prec) const;
    bool isEqualWithoutConsideringStr(const MEDCouplingTimeDiscretization *other, double prec) const;
    void checkNoTimePresence() const throw(INTERP_KERNEL::Exception);
    void checkTimePresence(double time) const throw(INTERP_KERNEL::Exception);
    void getArrays(std::vector<DataArrayDouble *>& arrays) const;
    void setEndArray(DataArrayDouble *array, TimeLabel *owner);
    void setStartTime(double time, int iteration, int order) throw(INTERP_KERNEL::Exception) { _start_time=time; _start_iteration=iteration; _start_order=order; }
    void setEndTime(double time, int iteration, int order) throw(INTERP_KERNEL::Exception) { _end_time=time; _end_iteration=iteration; _end_order=order; }
    double getStartTime(int& iteration, int& order) const throw(INTERP_KERNEL::Exception) { iteration=_start_iteration; order=_start_order; return _start_time; }
    double getEndTime(int& iteration, int& order) const throw(INTERP_KERNEL::Exception) { iteration=_end_iteration; order=_end_order; return _end_time; }
    void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    void getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const;
    void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    void resizeForUnserialization(const std::vector<int>& tinyInfoI, std::vector<DataArrayDouble *>& arrays);
    void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS);
    std::vector< const DataArrayDouble *> getArraysForTime(double time) const throw(INTERP_KERNEL::Exception);
    void setArrays(const std::vector<DataArrayDouble *>& arrays, TimeLabel *owner) throw(INTERP_KERNEL::Exception);
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

  class MEDCOUPLING_EXPORT MEDCouplingLinearTime : public MEDCouplingTwoTimeSteps
  {
  protected:
    MEDCouplingLinearTime(const MEDCouplingLinearTime& other, bool deepCpy);
  public:
    MEDCouplingLinearTime();
    std::string getStringRepr() const;
    TypeOfTimeDiscretization getEnum() const { return DISCRETIZATION; }
    void checkCoherency() const throw(INTERP_KERNEL::Exception);
    MEDCouplingTimeDiscretization *performCpy(bool deepCpy) const;
    bool areCompatible(const MEDCouplingTimeDiscretization *other) const;
    bool areStrictlyCompatible(const MEDCouplingTimeDiscretization *other) const;
    bool areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretization *other) const;
    void getValueForTime(double time, const std::vector<double>& vals, double *res) const;
    void getValueOnTime(int eltId, double time, double *value) const throw(INTERP_KERNEL::Exception);
    void getValueOnDiscTime(int eltId, int iteration, int order, double *value) const throw(INTERP_KERNEL::Exception);
    MEDCouplingTimeDiscretization *aggregate(const MEDCouplingTimeDiscretization *other) const;
    MEDCouplingTimeDiscretization *dot(const MEDCouplingTimeDiscretization *other) const;
    MEDCouplingTimeDiscretization *crossProduct(const MEDCouplingTimeDiscretization *other) const;
    MEDCouplingTimeDiscretization *max(const MEDCouplingTimeDiscretization *other) const;
    MEDCouplingTimeDiscretization *min(const MEDCouplingTimeDiscretization *other) const;
    MEDCouplingTimeDiscretization *add(const MEDCouplingTimeDiscretization *other) const;
    void addEqual(const MEDCouplingTimeDiscretization *other);
    MEDCouplingTimeDiscretization *substract(const MEDCouplingTimeDiscretization *other) const;
    void substractEqual(const MEDCouplingTimeDiscretization *other);
    MEDCouplingTimeDiscretization *multiply(const MEDCouplingTimeDiscretization *other) const;
    void multiplyEqual(const MEDCouplingTimeDiscretization *other);
    MEDCouplingTimeDiscretization *divide(const MEDCouplingTimeDiscretization *other) const;
    void divideEqual(const MEDCouplingTimeDiscretization *other);
  public:
    static const TypeOfTimeDiscretization DISCRETIZATION=LINEAR_TIME;
    static const char REPR[];
  };
}

#endif
