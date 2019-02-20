// Copyright (C) 2007-2019  CEA/DEN, EDF R&D
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

namespace MEDCoupling
{
  class MEDCouplingTimeDiscretization : public TimeLabel, public BigMemoryObject
  {
  public:
    static MEDCouplingTimeDiscretization *New(TypeOfTimeDiscretization type);
    void setTimeUnit(const char *unit);
    const char *getTimeUnit() const;
    virtual void copyTinyAttrFrom(const MEDCouplingTimeDiscretization& other);
    virtual void copyTinyStringsFrom(const MEDCouplingTimeDiscretization& other);
    virtual void checkConsistencyLight() const;
    virtual bool areCompatible(const MEDCouplingTimeDiscretization *other) const;
    virtual bool areStrictlyCompatible(const MEDCouplingTimeDiscretization *other, std::string& reason) const;
    virtual bool areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretization *other) const;
    virtual bool areStrictlyCompatibleForDiv(const MEDCouplingTimeDiscretization *other) const;
    virtual bool areCompatibleForMeld(const MEDCouplingTimeDiscretization *other) const;
    virtual bool isEqualIfNotWhy(const MEDCouplingTimeDiscretization *other, double prec, std::string& reason) const;
    virtual bool isEqual(const MEDCouplingTimeDiscretization *other, double prec) const;
    virtual bool isEqualWithoutConsideringStr(const MEDCouplingTimeDiscretization *other, double prec) const;
    virtual MEDCouplingTimeDiscretization *buildNewTimeReprFromThis(TypeOfTimeDiscretization type, bool deepCopy) const;
    virtual std::string getStringRepr() const;
    virtual TypeOfTimeDiscretization getEnum() const;
    virtual void synchronizeTimeWith(const MEDCouplingMesh *mesh);
    virtual MEDCouplingTimeDiscretization *aggregate(const MEDCouplingTimeDiscretization *other) const;
    virtual MEDCouplingTimeDiscretization *aggregate(const std::vector<const MEDCouplingTimeDiscretization *>& other) const;
    virtual MEDCouplingTimeDiscretization *meld(const MEDCouplingTimeDiscretization *other) const;
    virtual MEDCouplingTimeDiscretization *dot(const MEDCouplingTimeDiscretization *other) const;
    virtual MEDCouplingTimeDiscretization *crossProduct(const MEDCouplingTimeDiscretization *other) const;
    virtual MEDCouplingTimeDiscretization *max(const MEDCouplingTimeDiscretization *other) const;
    virtual MEDCouplingTimeDiscretization *min(const MEDCouplingTimeDiscretization *other) const;
    virtual MEDCouplingTimeDiscretization *add(const MEDCouplingTimeDiscretization *other) const;
    virtual void addEqual(const MEDCouplingTimeDiscretization *other);
    virtual MEDCouplingTimeDiscretization *substract(const MEDCouplingTimeDiscretization *other) const;
    virtual void substractEqual(const MEDCouplingTimeDiscretization *other);
    virtual MEDCouplingTimeDiscretization *multiply(const MEDCouplingTimeDiscretization *other) const;
    virtual void multiplyEqual(const MEDCouplingTimeDiscretization *other);
    virtual MEDCouplingTimeDiscretization *divide(const MEDCouplingTimeDiscretization *other) const;
    virtual void divideEqual(const MEDCouplingTimeDiscretization *other);
    virtual MEDCouplingTimeDiscretization *pow(const MEDCouplingTimeDiscretization *other) const;
    virtual void powEqual(const MEDCouplingTimeDiscretization *other);
    virtual MEDCouplingTimeDiscretization *performCopyOrIncrRef(bool deepCopy) const;
    void setTimeTolerance(double val);
    double getTimeTolerance() const;
    virtual void checkNoTimePresence() const;
    virtual void checkTimePresence(double time) const;
    virtual void setArray(DataArrayDouble *array, TimeLabel *owner);
    virtual void setEndArray(DataArrayDouble *array, TimeLabel *owner);
    virtual void setArrays(const std::vector<DataArrayDouble *>& arrays, TimeLabel *owner);
    DataArrayDouble *getArray();
    const DataArrayDouble *getArray() const;
    virtual const DataArrayDouble *getEndArray() const;
    virtual DataArrayDouble *getEndArray();
    virtual std::vector< const DataArrayDouble *> getArraysForTime(double time) const;
    virtual void getValueForTime(double time, const std::vector<double>& vals, double *res) const; 
    virtual void getArrays(std::vector<DataArrayDouble *>& arrays) const;
    virtual bool isBefore(const MEDCouplingTimeDiscretization *other) const;
    virtual bool isStrictlyBefore(const MEDCouplingTimeDiscretization *other) const;
    double getTime(int& iteration, int& order) const;
    virtual double getStartTime(int& iteration, int& order) const;
    virtual double getEndTime(int& iteration, int& order) const;
    void setTime(double time, int iteration, int order);
    void setIteration(int it);
    void setOrder(int order);
    void setTimeValue(double val);
    virtual void setStartIteration(int it);
    virtual void setEndIteration(int it);
    virtual void setStartOrder(int order);
    virtual void setEndOrder(int order);
    virtual void setStartTimeValue(double time);
    virtual void setEndTimeValue(double time);
    virtual void setStartTime(double time, int iteration, int order);
    virtual void setEndTime(double time, int iteration, int order);
    virtual void getValueOnTime(int eltId, double time, double *value) const;
    virtual void getValueOnDiscTime(int eltId, int iteration, int order, double *value) const;
    //
    virtual MEDCouplingTimeDiscretization *doublyContractedProduct() const;
    virtual MEDCouplingTimeDiscretization *determinant() const;
    virtual MEDCouplingTimeDiscretization *eigenValues() const;
    virtual MEDCouplingTimeDiscretization *eigenVectors() const;
    virtual MEDCouplingTimeDiscretization *inverse() const;
    virtual MEDCouplingTimeDiscretization *trace() const;
    virtual MEDCouplingTimeDiscretization *deviator() const;
    virtual MEDCouplingTimeDiscretization *magnitude() const;
    virtual MEDCouplingTimeDiscretization *negate() const;
    virtual MEDCouplingTimeDiscretization *maxPerTuple() const;
    virtual MEDCouplingTimeDiscretization *keepSelectedComponents(const std::vector<int>& compoIds) const;
    virtual void setSelectedComponents(const MEDCouplingTimeDiscretization *other, const std::vector<int>& compoIds);
    virtual void changeNbOfComponents(int newNbOfComp, double dftValue);
    virtual void sortPerTuple(bool asc);
    virtual void setUniformValue(int nbOfTuple, int nbOfCompo, double value);
    virtual void setOrCreateUniformValueOnAllComponents(int nbOfTuple, double value);
    virtual void applyLin(double a, double b, int compoId);
    virtual void applyFunc(int nbOfComp, FunctionToEvaluate func);
    virtual void applyFunc(int nbOfComp, const char *func);
    virtual void applyFuncCompo(int nbOfComp, const char *func);
    virtual void applyFuncNamedCompo(int nbOfComp, const std::vector<std::string>& varsOrder, const char *func);
    virtual void applyFunc(const char *func);
    virtual void applyFuncFast32(const char *func);
    virtual void applyFuncFast64(const char *func);
    virtual void fillFromAnalytic(const DataArrayDouble *loc, int nbOfComp, FunctionToEvaluate func);
    virtual void fillFromAnalytic(const DataArrayDouble *loc, int nbOfComp, const char *func);
    virtual void fillFromAnalyticCompo(const DataArrayDouble *loc, int nbOfComp, const char *func);
    virtual void fillFromAnalyticNamedCompo(const DataArrayDouble *loc, int nbOfComp, const std::vector<std::string>& varsOrder, const char *func);
    //
    virtual ~MEDCouplingTimeDiscretization();
  };

  class MEDCouplingNoTimeLabel : public MEDCouplingTimeDiscretization
  {
  public:
    MEDCouplingNoTimeLabel();
    MEDCouplingNoTimeLabel(const MEDCouplingTimeDiscretization& other, bool deepCopy);
  public:
    static const TypeOfTimeDiscretization DISCRETIZATION=NO_TIME;
    static const char REPR[];
  };

  class MEDCouplingWithTimeStep : public MEDCouplingTimeDiscretization
  {
  public:
    MEDCouplingWithTimeStep();
  public:
    static const TypeOfTimeDiscretization DISCRETIZATION=ONE_TIME;
    static const char REPR[];
  };

  class MEDCouplingConstOnTimeInterval : public MEDCouplingTimeDiscretization
  {
  protected:
    MEDCouplingConstOnTimeInterval();
    MEDCouplingConstOnTimeInterval(const MEDCouplingConstOnTimeInterval& other, bool deepCopy);
  public:
    static const TypeOfTimeDiscretization DISCRETIZATION=CONST_ON_TIME_INTERVAL;
    static const char REPR[];
  };

  class MEDCouplingTwoTimeSteps : public MEDCouplingTimeDiscretization
  {
  };

  class MEDCouplingLinearTime : public MEDCouplingTwoTimeSteps
  {
  public:
    MEDCouplingLinearTime();
  public:
    static const TypeOfTimeDiscretization DISCRETIZATION=LINEAR_TIME;
    static const char REPR[];
  };
}

#endif
