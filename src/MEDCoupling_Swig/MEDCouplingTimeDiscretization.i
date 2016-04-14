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

namespace MEDCoupling
{
  class MEDCouplingTimeDiscretization : public TimeLabel, public BigMemoryObject
  {
  public:
    static MEDCouplingTimeDiscretization *New(TypeOfTimeDiscretization type) throw(INTERP_KERNEL::Exception);
    void setTimeUnit(const char *unit);
    const char *getTimeUnit() const;
    virtual void copyTinyAttrFrom(const MEDCouplingTimeDiscretization& other) throw(INTERP_KERNEL::Exception);
    virtual void copyTinyStringsFrom(const MEDCouplingTimeDiscretization& other) throw(INTERP_KERNEL::Exception);
    virtual void checkConsistencyLight() const throw(INTERP_KERNEL::Exception);
    virtual bool areCompatible(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    virtual bool areStrictlyCompatible(const MEDCouplingTimeDiscretization *other, std::string& reason) const throw(INTERP_KERNEL::Exception);
    virtual bool areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    virtual bool areStrictlyCompatibleForDiv(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    virtual bool areCompatibleForMeld(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    virtual bool isEqualIfNotWhy(const MEDCouplingTimeDiscretization *other, double prec, std::string& reason) const throw(INTERP_KERNEL::Exception);
    virtual bool isEqual(const MEDCouplingTimeDiscretization *other, double prec) const throw(INTERP_KERNEL::Exception);
    virtual bool isEqualWithoutConsideringStr(const MEDCouplingTimeDiscretization *other, double prec) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *buildNewTimeReprFromThis(TypeOfTimeDiscretization type, bool deepCopy) const throw(INTERP_KERNEL::Exception);
    virtual std::string getStringRepr() const throw(INTERP_KERNEL::Exception);
    virtual TypeOfTimeDiscretization getEnum() const throw(INTERP_KERNEL::Exception);
    virtual void synchronizeTimeWith(const MEDCouplingMesh *mesh) throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *aggregate(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *aggregate(const std::vector<const MEDCouplingTimeDiscretization *>& other) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *meld(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *dot(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *crossProduct(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *max(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *min(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *add(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    virtual void addEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *substract(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    virtual void substractEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *multiply(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    virtual void multiplyEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *divide(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    virtual void divideEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *pow(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    virtual void powEqual(const MEDCouplingTimeDiscretization *other) throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *performCopyOrIncrRef(bool deepCopy) const throw(INTERP_KERNEL::Exception);
    void setTimeTolerance(double val);
    double getTimeTolerance() const;
    virtual void checkNoTimePresence() const throw(INTERP_KERNEL::Exception);
    virtual void checkTimePresence(double time) const throw(INTERP_KERNEL::Exception);
    virtual void setArray(DataArrayDouble *array, TimeLabel *owner) throw(INTERP_KERNEL::Exception);
    virtual void setEndArray(DataArrayDouble *array, TimeLabel *owner) throw(INTERP_KERNEL::Exception);
    virtual void setArrays(const std::vector<DataArrayDouble *>& arrays, TimeLabel *owner) throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getArray() throw(INTERP_KERNEL::Exception);
    const DataArrayDouble *getArray() const throw(INTERP_KERNEL::Exception);
    virtual const DataArrayDouble *getEndArray() const throw(INTERP_KERNEL::Exception);
    virtual DataArrayDouble *getEndArray() throw(INTERP_KERNEL::Exception);
    virtual std::vector< const DataArrayDouble *> getArraysForTime(double time) const throw(INTERP_KERNEL::Exception);
    virtual void getValueForTime(double time, const std::vector<double>& vals, double *res) const throw(INTERP_KERNEL::Exception); 
    virtual void getArrays(std::vector<DataArrayDouble *>& arrays) const throw(INTERP_KERNEL::Exception);
    virtual bool isBefore(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    virtual bool isStrictlyBefore(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    double getTime(int& iteration, int& order) const throw(INTERP_KERNEL::Exception);
    virtual double getStartTime(int& iteration, int& order) const throw(INTERP_KERNEL::Exception);
    virtual double getEndTime(int& iteration, int& order) const throw(INTERP_KERNEL::Exception);
    void setTime(double time, int iteration, int order) throw(INTERP_KERNEL::Exception);
    void setIteration(int it) throw(INTERP_KERNEL::Exception);
    void setOrder(int order) throw(INTERP_KERNEL::Exception);
    void setTimeValue(double val) throw(INTERP_KERNEL::Exception);
    virtual void setStartIteration(int it) throw(INTERP_KERNEL::Exception);
    virtual void setEndIteration(int it) throw(INTERP_KERNEL::Exception);
    virtual void setStartOrder(int order) throw(INTERP_KERNEL::Exception);
    virtual void setEndOrder(int order) throw(INTERP_KERNEL::Exception);
    virtual void setStartTimeValue(double time) throw(INTERP_KERNEL::Exception);
    virtual void setEndTimeValue(double time) throw(INTERP_KERNEL::Exception);
    virtual void setStartTime(double time, int iteration, int order) throw(INTERP_KERNEL::Exception);
    virtual void setEndTime(double time, int iteration, int order) throw(INTERP_KERNEL::Exception);
    virtual void getValueOnTime(int eltId, double time, double *value) const throw(INTERP_KERNEL::Exception);
    virtual void getValueOnDiscTime(int eltId, int iteration, int order, double *value) const throw(INTERP_KERNEL::Exception);
    //
    virtual MEDCouplingTimeDiscretization *doublyContractedProduct() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *determinant() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *eigenValues() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *eigenVectors() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *inverse() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *trace() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *deviator() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *magnitude() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *negate() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *maxPerTuple() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingTimeDiscretization *keepSelectedComponents(const std::vector<int>& compoIds) const throw(INTERP_KERNEL::Exception);
    virtual void setSelectedComponents(const MEDCouplingTimeDiscretization *other, const std::vector<int>& compoIds) throw(INTERP_KERNEL::Exception);
    virtual void changeNbOfComponents(int newNbOfComp, double dftValue) throw(INTERP_KERNEL::Exception);
    virtual void sortPerTuple(bool asc) throw(INTERP_KERNEL::Exception);
    virtual void setUniformValue(int nbOfTuple, int nbOfCompo, double value) throw(INTERP_KERNEL::Exception);
    virtual void setOrCreateUniformValueOnAllComponents(int nbOfTuple, double value) throw(INTERP_KERNEL::Exception);
    virtual void applyLin(double a, double b, int compoId) throw(INTERP_KERNEL::Exception);
    virtual void applyFunc(int nbOfComp, FunctionToEvaluate func) throw(INTERP_KERNEL::Exception);
    virtual void applyFunc(int nbOfComp, const char *func) throw(INTERP_KERNEL::Exception);
    virtual void applyFuncCompo(int nbOfComp, const char *func) throw(INTERP_KERNEL::Exception);
    virtual void applyFuncNamedCompo(int nbOfComp, const std::vector<std::string>& varsOrder, const char *func) throw(INTERP_KERNEL::Exception);
    virtual void applyFunc(const char *func) throw(INTERP_KERNEL::Exception);
    virtual void applyFuncFast32(const char *func) throw(INTERP_KERNEL::Exception);
    virtual void applyFuncFast64(const char *func) throw(INTERP_KERNEL::Exception);
    virtual void fillFromAnalytic(const DataArrayDouble *loc, int nbOfComp, FunctionToEvaluate func) throw(INTERP_KERNEL::Exception);
    virtual void fillFromAnalytic(const DataArrayDouble *loc, int nbOfComp, const char *func) throw(INTERP_KERNEL::Exception);
    virtual void fillFromAnalyticCompo(const DataArrayDouble *loc, int nbOfComp, const char *func) throw(INTERP_KERNEL::Exception);
    virtual void fillFromAnalyticNamedCompo(const DataArrayDouble *loc, int nbOfComp, const std::vector<std::string>& varsOrder, const char *func) throw(INTERP_KERNEL::Exception);
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
