// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
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

#ifndef __PARAMEDMEM_MEDCOUPLINGFIELDDOUBLE_HXX__
#define __PARAMEDMEM_MEDCOUPLINGFIELDDOUBLE_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingField.hxx"
#include "MEDCouplingTimeDiscretization.hxx"
#include "MEDCouplingMemArray.hxx"

namespace ParaMEDMEM
{
  class MEDCouplingFieldTemplate;

  class MEDCOUPLING_EXPORT MEDCouplingFieldDouble : public MEDCouplingField
  {
  public:
    static MEDCouplingFieldDouble *New(TypeOfField type, TypeOfTimeDiscretization td=NO_TIME);
    static MEDCouplingFieldDouble *New(const MEDCouplingFieldTemplate *ft, TypeOfTimeDiscretization td=NO_TIME);
    void setTimeUnit(const char *unit);
    const char *getTimeUnit() const;
    void copyTinyStringsFrom(const MEDCouplingFieldDouble *other) throw(INTERP_KERNEL::Exception);
    void copyTinyAttrFrom(const MEDCouplingFieldDouble *other) throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const;
    std::string advancedRepr() const;
    bool isEqualIfNotWhy(const MEDCouplingField *other, double meshPrec, double valsPrec, std::string& reason) const throw(INTERP_KERNEL::Exception);
    bool isEqualWithoutConsideringStr(const MEDCouplingField *other, double meshPrec, double valsPrec) const;
    bool areCompatibleForMerge(const MEDCouplingField *other) const;
    bool areStrictlyCompatible(const MEDCouplingField *other) const;
    bool areCompatibleForMul(const MEDCouplingField *other) const;
    bool areCompatibleForDiv(const MEDCouplingField *other) const;
    bool areCompatibleForMeld(const MEDCouplingFieldDouble *other) const;
    void renumberCells(const int *old2NewBg, bool check=true) throw(INTERP_KERNEL::Exception);
    void renumberCellsWithoutMesh(const int *old2NewBg, bool check=true) throw(INTERP_KERNEL::Exception);
    void renumberNodes(const int *old2NewBg) throw(INTERP_KERNEL::Exception);
    void renumberNodesWithoutMesh(const int *old2NewBg, double eps=1e-15) throw(INTERP_KERNEL::Exception);
    DataArrayInt *getIdsInRange(double vmin, double vmax) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *buildSubPart(const DataArrayInt *part) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *buildSubPart(const int *partBg, const int *partEnd) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *deepCpy() const;
    MEDCouplingFieldDouble *clone(bool recDeepCpy) const;
    MEDCouplingFieldDouble *cloneWithMesh(bool recDeepCpy) const;
    MEDCouplingFieldDouble *buildNewTimeReprFromThis(TypeOfTimeDiscretization td, bool deepCopy) const;
    TypeOfTimeDiscretization getTimeDiscretization() const;
    void checkCoherency() const throw(INTERP_KERNEL::Exception);
    void setNature(NatureOfField nat) throw(INTERP_KERNEL::Exception);
    void setTimeTolerance(double val) { _time_discr->setTimeTolerance(val); }
    double getTimeTolerance() const { return _time_discr->getTimeTolerance(); }
    void setIteration(int it) throw(INTERP_KERNEL::Exception) { _time_discr->setIteration(it); }
    void setEndIteration(int it) throw(INTERP_KERNEL::Exception) { _time_discr->setEndIteration(it); }
    void setOrder(int order) throw(INTERP_KERNEL::Exception) { _time_discr->setOrder(order); }
    void setEndOrder(int order) throw(INTERP_KERNEL::Exception) { _time_discr->setEndOrder(order); }
    void setTimeValue(double val) throw(INTERP_KERNEL::Exception) { _time_discr->setTimeValue(val); }
    void setEndTimeValue(double val) throw(INTERP_KERNEL::Exception) { _time_discr->setEndTimeValue(val); }
    void setTime(double val, int iteration, int order) { _time_discr->setTime(val,iteration,order); }
    void setStartTime(double val, int iteration, int order) { _time_discr->setStartTime(val,iteration,order); }
    void setEndTime(double val, int iteration, int order) { _time_discr->setEndTime(val,iteration,order); }
    double getTime(int& iteration, int& order) const { return _time_discr->getTime(iteration,order); }
    double getStartTime(int& iteration, int& order) const { return _time_discr->getStartTime(iteration,order); }
    double getEndTime(int& iteration, int& order) const { return _time_discr->getEndTime(iteration,order); }
    double getIJ(int tupleId, int compoId) const { return getArray()->getIJ(tupleId,compoId); }
    double getIJK(int cellId, int nodeIdInCell, int compoId) const;
    void setArray(DataArrayDouble *array);
    void setEndArray(DataArrayDouble *array);
    void setArrays(const std::vector<DataArrayDouble *>& arrs) throw(INTERP_KERNEL::Exception);
    const DataArrayDouble *getArray() const { return _time_discr->getArray(); }
    DataArrayDouble *getArray() { return _time_discr->getArray(); }
    const DataArrayDouble *getEndArray() const { return _time_discr->getEndArray(); }
    DataArrayDouble *getEndArray() { return _time_discr->getEndArray(); }
    std::vector<DataArrayDouble *> getArrays() const { std::vector<DataArrayDouble *> ret; _time_discr->getArrays(ret); return ret; }
    double accumulate(int compId) const;
    void accumulate(double *res) const;
    double getMaxValue() const throw(INTERP_KERNEL::Exception);
    double getMaxValue2(DataArrayInt*& tupleIds) const throw(INTERP_KERNEL::Exception);
    double getMinValue() const throw(INTERP_KERNEL::Exception);
    double getMinValue2(DataArrayInt*& tupleIds) const throw(INTERP_KERNEL::Exception);
    double getAverageValue() const throw(INTERP_KERNEL::Exception);
    double norm2() const throw(INTERP_KERNEL::Exception);
    double normMax() const throw(INTERP_KERNEL::Exception);
    double getWeightedAverageValue() const throw(INTERP_KERNEL::Exception);
    double normL1(int compId) const throw(INTERP_KERNEL::Exception);
    void normL1(double *res) const throw(INTERP_KERNEL::Exception);
    double normL2(int compId) const throw(INTERP_KERNEL::Exception);
    void normL2(double *res) const throw(INTERP_KERNEL::Exception);
    double integral(int compId, bool isWAbs) const throw(INTERP_KERNEL::Exception);
    void integral(bool isWAbs, double *res) const throw(INTERP_KERNEL::Exception);
    void getValueOnPos(int i, int j, int k, double *res) const throw(INTERP_KERNEL::Exception);
    void getValueOn(const double *spaceLoc, double *res) const throw(INTERP_KERNEL::Exception);
    void getValueOn(const double *spaceLoc, double time, double *res) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getValueOnMulti(const double *spaceLoc, int nbOfPoints) const throw(INTERP_KERNEL::Exception);
    //! \b temporary
    void applyLin(double a, double b, int compoId);
    MEDCouplingFieldDouble &operator=(double value) throw(INTERP_KERNEL::Exception);
    void fillFromAnalytic(int nbOfComp, FunctionToEvaluate func) throw(INTERP_KERNEL::Exception);
    void fillFromAnalytic(int nbOfComp, const char *func) throw(INTERP_KERNEL::Exception);
    void fillFromAnalytic2(int nbOfComp, const char *func) throw(INTERP_KERNEL::Exception);
    void fillFromAnalytic3(int nbOfComp, const std::vector<std::string>& varsOrder, const char *func) throw(INTERP_KERNEL::Exception);
    void applyFunc(int nbOfComp, FunctionToEvaluate func);
    void applyFunc(int nbOfComp, double val);
    void applyFunc(int nbOfComp, const char *func) throw(INTERP_KERNEL::Exception);
    void applyFunc2(int nbOfComp, const char *func) throw(INTERP_KERNEL::Exception);
    void applyFunc3(int nbOfComp, const std::vector<std::string>& varsOrder, const char *func) throw(INTERP_KERNEL::Exception);
    void applyFunc(const char *func) throw(INTERP_KERNEL::Exception);
    void applyFuncFast32(const char *func) throw(INTERP_KERNEL::Exception);
    void applyFuncFast64(const char *func) throw(INTERP_KERNEL::Exception);
    int getNumberOfComponents() const throw(INTERP_KERNEL::Exception);
    int getNumberOfTuples() const throw(INTERP_KERNEL::Exception);
    int getNumberOfValues() const throw(INTERP_KERNEL::Exception);
    void updateTime() const;
    //
    void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    void getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const;
    void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    void resizeForUnserialization(const std::vector<int>& tinyInfoI, DataArrayInt *&dataInt, std::vector<DataArrayDouble *>& arrays);
    void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS);
    void serialize(DataArrayInt *&dataInt, std::vector<DataArrayDouble *>& arrays) const;
    //
    void changeUnderlyingMesh(const MEDCouplingMesh *other, int levOfCheck, double prec) throw(INTERP_KERNEL::Exception);
    void substractInPlaceDM(const MEDCouplingFieldDouble *f, int levOfCheck, double prec) throw(INTERP_KERNEL::Exception);
    bool mergeNodes(double eps, double epsOnVals=1e-15) throw(INTERP_KERNEL::Exception);
    bool mergeNodes2(double eps, double epsOnVals=1e-15) throw(INTERP_KERNEL::Exception);
    bool zipCoords(double epsOnVals=1e-15) throw(INTERP_KERNEL::Exception);
    bool zipConnectivity(int compType, double epsOnVals=1e-15) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *extractSlice3D(const double *origin, const double *vec, double eps) const throw(INTERP_KERNEL::Exception);
    bool simplexize(int policy) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *doublyContractedProduct() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *determinant() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *eigenValues() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *eigenVectors() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *inverse() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *trace() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *deviator() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *magnitude() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *maxPerTuple() const throw(INTERP_KERNEL::Exception);
    void changeNbOfComponents(int newNbOfComp, double dftValue=0.) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *keepSelectedComponents(const std::vector<int>& compoIds) const throw(INTERP_KERNEL::Exception);
    void setSelectedComponents(const MEDCouplingFieldDouble *f, const std::vector<int>& compoIds) throw(INTERP_KERNEL::Exception);
    void sortPerTuple(bool asc) throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldDouble *MergeFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldDouble *MergeFields(const std::vector<const MEDCouplingFieldDouble *>& a) throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldDouble *MeldFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldDouble *DotFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *dot(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception) { return DotFields(this,&other); }
    static MEDCouplingFieldDouble *CrossProductFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *crossProduct(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception) { return CrossProductFields(this,&other); }
    static MEDCouplingFieldDouble *MaxFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *max(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception) { return MaxFields(this,&other); }
    static MEDCouplingFieldDouble *MinFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *min(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception) { return MinFields(this,&other); }
    MEDCouplingFieldDouble *operator+(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception) { return AddFields(this,&other); }
    const MEDCouplingFieldDouble &operator+=(const MEDCouplingFieldDouble& other) throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldDouble *AddFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *operator-(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception) { return SubstractFields(this,&other); }
    const MEDCouplingFieldDouble &operator-=(const MEDCouplingFieldDouble& other) throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldDouble *SubstractFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *operator*(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception) { return MultiplyFields(this,&other); }
    const MEDCouplingFieldDouble &operator*=(const MEDCouplingFieldDouble& other) throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldDouble *MultiplyFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *operator/(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception) { return DivideFields(this,&other); }
    const MEDCouplingFieldDouble &operator/=(const MEDCouplingFieldDouble& other) throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldDouble *DivideFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
    static void WriteVTK(const char *fileName, const std::vector<const MEDCouplingFieldDouble *>& fs) throw(INTERP_KERNEL::Exception);
  public:
    const MEDCouplingTimeDiscretization *getTimeDiscretizationUnderGround() const { return _time_discr; }
    MEDCouplingTimeDiscretization *getTimeDiscretizationUnderGround() { return _time_discr; }
  private:
    MEDCouplingFieldDouble(TypeOfField type, TypeOfTimeDiscretization td);
    MEDCouplingFieldDouble(const MEDCouplingFieldTemplate *ft, TypeOfTimeDiscretization td);
    MEDCouplingFieldDouble(const MEDCouplingFieldDouble& other, bool deepCopy);
    MEDCouplingFieldDouble(NatureOfField n, MEDCouplingTimeDiscretization *td, MEDCouplingFieldDiscretization *type);
    ~MEDCouplingFieldDouble();
  private:
    MEDCouplingTimeDiscretization *_time_discr;
  };
}

#endif
