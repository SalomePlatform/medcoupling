// Copyright (C) 2007-2015  CEA/DEN, EDF R&D
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

#ifndef __PARAMEDMEM_MEDCOUPLINGFIELDDOUBLE_HXX__
#define __PARAMEDMEM_MEDCOUPLINGFIELDDOUBLE_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingField.hxx"
#include "MEDCouplingTimeDiscretization.hxx"
#include "MEDCouplingMemArray.hxx"

namespace ParaMEDMEM
{
  class MEDCouplingFieldTemplate;

  class MEDCouplingFieldDouble : public MEDCouplingField
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingFieldDouble *New(TypeOfField type, TypeOfTimeDiscretization td=ONE_TIME);
    MEDCOUPLING_EXPORT static MEDCouplingFieldDouble *New(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td=ONE_TIME);
    MEDCOUPLING_EXPORT void setTimeUnit(const std::string& unit);
    MEDCOUPLING_EXPORT std::string getTimeUnit() const;
    MEDCOUPLING_EXPORT void synchronizeTimeWithSupport();
    MEDCOUPLING_EXPORT void copyTinyStringsFrom(const MEDCouplingField *other);
    MEDCOUPLING_EXPORT void copyTinyAttrFrom(const MEDCouplingFieldDouble *other);
    MEDCOUPLING_EXPORT void copyAllTinyAttrFrom(const MEDCouplingFieldDouble *other);
    MEDCOUPLING_EXPORT std::string simpleRepr() const;
    MEDCOUPLING_EXPORT std::string advancedRepr() const;
    MEDCOUPLING_EXPORT std::string writeVTK(const std::string& fileName, bool isBinary=true) const;
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingField *other, double meshPrec, double valsPrec, std::string& reason) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const MEDCouplingField *other, double meshPrec, double valsPrec) const;
    MEDCOUPLING_EXPORT bool areCompatibleForMerge(const MEDCouplingField *other) const;
    MEDCOUPLING_EXPORT bool areStrictlyCompatible(const MEDCouplingField *other) const;
    MEDCOUPLING_EXPORT bool areCompatibleForMul(const MEDCouplingField *other) const;
    MEDCOUPLING_EXPORT bool areCompatibleForDiv(const MEDCouplingField *other) const;
    MEDCOUPLING_EXPORT bool areCompatibleForMeld(const MEDCouplingFieldDouble *other) const;
    MEDCOUPLING_EXPORT void renumberCells(const int *old2NewBg, bool check=true);
    MEDCOUPLING_EXPORT void renumberCellsWithoutMesh(const int *old2NewBg, bool check=true);
    MEDCOUPLING_EXPORT void renumberNodes(const int *old2NewBg, double eps=1e-15);
    MEDCOUPLING_EXPORT void renumberNodesWithoutMesh(const int *old2NewBg, int newNbOfNodes, double eps=1e-15);
    MEDCOUPLING_EXPORT DataArrayInt *getIdsInRange(double vmin, double vmax) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildSubPart(const DataArrayInt *part) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildSubPart(const int *partBg, const int *partEnd) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildSubPartRange(int begin, int end, int step) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *deepCpy() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *clone(bool recDeepCpy) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *cloneWithMesh(bool recDeepCpy) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildNewTimeReprFromThis(TypeOfTimeDiscretization td, bool deepCopy) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *nodeToCellDiscretization() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *cellToNodeDiscretization() const;
    MEDCOUPLING_EXPORT TypeOfTimeDiscretization getTimeDiscretization() const;
    MEDCOUPLING_EXPORT void checkCoherency() const;
    MEDCOUPLING_EXPORT void setNature(NatureOfField nat);
    MEDCOUPLING_EXPORT void setTimeTolerance(double val) { _time_discr->setTimeTolerance(val); }
    MEDCOUPLING_EXPORT double getTimeTolerance() const { return _time_discr->getTimeTolerance(); }
    MEDCOUPLING_EXPORT void setIteration(int it) { _time_discr->setIteration(it); }
    MEDCOUPLING_EXPORT void setEndIteration(int it) { _time_discr->setEndIteration(it); }
    MEDCOUPLING_EXPORT void setOrder(int order) { _time_discr->setOrder(order); }
    MEDCOUPLING_EXPORT void setEndOrder(int order) { _time_discr->setEndOrder(order); }
    MEDCOUPLING_EXPORT void setTimeValue(double val) { _time_discr->setTimeValue(val); }
    MEDCOUPLING_EXPORT void setEndTimeValue(double val) { _time_discr->setEndTimeValue(val); }
    MEDCOUPLING_EXPORT void setTime(double val, int iteration, int order) { _time_discr->setTime(val,iteration,order); }
    MEDCOUPLING_EXPORT void synchronizeTimeWithMesh();
    MEDCOUPLING_EXPORT void setStartTime(double val, int iteration, int order) { _time_discr->setStartTime(val,iteration,order); }
    MEDCOUPLING_EXPORT void setEndTime(double val, int iteration, int order) { _time_discr->setEndTime(val,iteration,order); }
    MEDCOUPLING_EXPORT double getTime(int& iteration, int& order) const { return _time_discr->getTime(iteration,order); }
    MEDCOUPLING_EXPORT double getStartTime(int& iteration, int& order) const { return _time_discr->getStartTime(iteration,order); }
    MEDCOUPLING_EXPORT double getEndTime(int& iteration, int& order) const { return _time_discr->getEndTime(iteration,order); }
    MEDCOUPLING_EXPORT double getIJ(int tupleId, int compoId) const { return getArray()->getIJ(tupleId,compoId); }
    MEDCOUPLING_EXPORT double getIJK(int cellId, int nodeIdInCell, int compoId) const;
    MEDCOUPLING_EXPORT void setArray(DataArrayDouble *array);
    MEDCOUPLING_EXPORT void setEndArray(DataArrayDouble *array);
    MEDCOUPLING_EXPORT void setArrays(const std::vector<DataArrayDouble *>& arrs);
    MEDCOUPLING_EXPORT const DataArrayDouble *getArray() const { return _time_discr->getArray(); }
    MEDCOUPLING_EXPORT DataArrayDouble *getArray() { return _time_discr->getArray(); }
    MEDCOUPLING_EXPORT const DataArrayDouble *getEndArray() const { return _time_discr->getEndArray(); }
    MEDCOUPLING_EXPORT DataArrayDouble *getEndArray() { return _time_discr->getEndArray(); }
    MEDCOUPLING_EXPORT std::vector<DataArrayDouble *> getArrays() const { std::vector<DataArrayDouble *> ret; _time_discr->getArrays(ret); return ret; }
    MEDCOUPLING_EXPORT double accumulate(int compId) const;
    MEDCOUPLING_EXPORT void accumulate(double *res) const;
    MEDCOUPLING_EXPORT double getMaxValue() const;
    MEDCOUPLING_EXPORT double getMaxValue2(DataArrayInt*& tupleIds) const;
    MEDCOUPLING_EXPORT double getMinValue() const;
    MEDCOUPLING_EXPORT double getMinValue2(DataArrayInt*& tupleIds) const;
    MEDCOUPLING_EXPORT double getAverageValue() const;
    MEDCOUPLING_EXPORT double norm2() const;
    MEDCOUPLING_EXPORT double normMax() const;
    MEDCOUPLING_EXPORT void getWeightedAverageValue(double *res, bool isWAbs=true) const;
    MEDCOUPLING_EXPORT double getWeightedAverageValue(int compId, bool isWAbs=true) const;
    MEDCOUPLING_EXPORT double normL1(int compId) const;
    MEDCOUPLING_EXPORT void normL1(double *res) const;
    MEDCOUPLING_EXPORT double normL2(int compId) const;
    MEDCOUPLING_EXPORT void normL2(double *res) const;
    MEDCOUPLING_EXPORT double integral(int compId, bool isWAbs) const;
    MEDCOUPLING_EXPORT void integral(bool isWAbs, double *res) const;
    MEDCOUPLING_EXPORT void getValueOnPos(int i, int j, int k, double *res) const;
    MEDCOUPLING_EXPORT void getValueOn(const double *spaceLoc, double *res) const;
    MEDCOUPLING_EXPORT void getValueOn(const double *spaceLoc, double time, double *res) const;
    MEDCOUPLING_EXPORT DataArrayDouble *getValueOnMulti(const double *spaceLoc, int nbOfPoints) const;
    MEDCOUPLING_EXPORT void applyLin(double a, double b, int compoId);
    MEDCOUPLING_EXPORT void applyLin(double a, double b);
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble &operator=(double value);
    MEDCOUPLING_EXPORT void fillFromAnalytic(int nbOfComp, FunctionToEvaluate func);
    MEDCOUPLING_EXPORT void fillFromAnalytic(int nbOfComp, const std::string& func);
    MEDCOUPLING_EXPORT void fillFromAnalytic2(int nbOfComp, const std::string& func);
    MEDCOUPLING_EXPORT void fillFromAnalytic3(int nbOfComp, const std::vector<std::string>& varsOrder, const std::string& func);
    MEDCOUPLING_EXPORT void applyFunc(int nbOfComp, FunctionToEvaluate func);
    MEDCOUPLING_EXPORT void applyFunc(int nbOfComp, double val);
    MEDCOUPLING_EXPORT void applyFunc(int nbOfComp, const std::string& func);
    MEDCOUPLING_EXPORT void applyFunc2(int nbOfComp, const std::string& func);
    MEDCOUPLING_EXPORT void applyFunc3(int nbOfComp, const std::vector<std::string>& varsOrder, const std::string& func);
    MEDCOUPLING_EXPORT void applyFunc(const std::string& func);
    MEDCOUPLING_EXPORT void applyFuncFast32(const std::string& func);
    MEDCOUPLING_EXPORT void applyFuncFast64(const std::string& func);
    MEDCOUPLING_EXPORT int getNumberOfComponents() const;
    MEDCOUPLING_EXPORT int getNumberOfTuples() const;
    MEDCOUPLING_EXPORT int getNumberOfValues() const;
    MEDCOUPLING_EXPORT void updateTime() const;
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    //
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    MEDCOUPLING_EXPORT void getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const;
    MEDCOUPLING_EXPORT void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    MEDCOUPLING_EXPORT void resizeForUnserialization(const std::vector<int>& tinyInfoI, DataArrayInt *&dataInt, std::vector<DataArrayDouble *>& arrays);
    MEDCOUPLING_EXPORT void checkForUnserialization(const std::vector<int>& tinyInfoI, const DataArrayInt *dataInt, const std::vector<DataArrayDouble *>& arrays);
    MEDCOUPLING_EXPORT void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS);
    MEDCOUPLING_EXPORT void serialize(DataArrayInt *&dataInt, std::vector<DataArrayDouble *>& arrays) const;
    //
    MEDCOUPLING_EXPORT void changeUnderlyingMesh(const MEDCouplingMesh *other, int levOfCheck, double precOnMesh, double eps=1e-15);
    MEDCOUPLING_EXPORT void substractInPlaceDM(const MEDCouplingFieldDouble *f, int levOfCheck, double precOnMesh, double eps=1e-15);
    MEDCOUPLING_EXPORT bool mergeNodes(double eps, double epsOnVals=1e-15);
    MEDCOUPLING_EXPORT bool mergeNodes2(double eps, double epsOnVals=1e-15);
    MEDCOUPLING_EXPORT bool zipCoords(double epsOnVals=1e-15);
    MEDCOUPLING_EXPORT bool zipConnectivity(int compType, double epsOnVals=1e-15);
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *extractSlice3D(const double *origin, const double *vec, double eps) const;
    MEDCOUPLING_EXPORT bool simplexize(int policy);
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *doublyContractedProduct() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *determinant() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *eigenValues() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *eigenVectors() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *inverse() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *trace() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *deviator() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *magnitude() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *maxPerTuple() const;
    MEDCOUPLING_EXPORT void changeNbOfComponents(int newNbOfComp, double dftValue=0.);
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *keepSelectedComponents(const std::vector<int>& compoIds) const;
    MEDCOUPLING_EXPORT void setSelectedComponents(const MEDCouplingFieldDouble *f, const std::vector<int>& compoIds);
    MEDCOUPLING_EXPORT void sortPerTuple(bool asc);
    MEDCOUPLING_EXPORT static MEDCouplingFieldDouble *MergeFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCOUPLING_EXPORT static MEDCouplingFieldDouble *MergeFields(const std::vector<const MEDCouplingFieldDouble *>& a);
    MEDCOUPLING_EXPORT static MEDCouplingFieldDouble *MeldFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCOUPLING_EXPORT static MEDCouplingFieldDouble *DotFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *dot(const MEDCouplingFieldDouble& other) const { return DotFields(this,&other); }
    MEDCOUPLING_EXPORT static MEDCouplingFieldDouble *CrossProductFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *crossProduct(const MEDCouplingFieldDouble& other) const { return CrossProductFields(this,&other); }
    MEDCOUPLING_EXPORT static MEDCouplingFieldDouble *MaxFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *max(const MEDCouplingFieldDouble& other) const { return MaxFields(this,&other); }
    MEDCOUPLING_EXPORT static MEDCouplingFieldDouble *MinFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *min(const MEDCouplingFieldDouble& other) const { return MinFields(this,&other); }
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *negate() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *operator+(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception) { return AddFields(this,&other); }
    MEDCOUPLING_EXPORT const MEDCouplingFieldDouble &operator+=(const MEDCouplingFieldDouble& other);
    MEDCOUPLING_EXPORT static MEDCouplingFieldDouble *AddFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *operator-(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception) { return SubstractFields(this,&other); }
    MEDCOUPLING_EXPORT const MEDCouplingFieldDouble &operator-=(const MEDCouplingFieldDouble& other);
    MEDCOUPLING_EXPORT static MEDCouplingFieldDouble *SubstractFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCouplingFieldDouble *operator*(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception) { return MultiplyFields(this,&other); }
    MEDCOUPLING_EXPORT const MEDCouplingFieldDouble &operator*=(const MEDCouplingFieldDouble& other);
    MEDCOUPLING_EXPORT static MEDCouplingFieldDouble *MultiplyFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *operator/(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception) { return DivideFields(this,&other); }
    MEDCOUPLING_EXPORT const MEDCouplingFieldDouble &operator/=(const MEDCouplingFieldDouble& other);
    MEDCOUPLING_EXPORT static MEDCouplingFieldDouble *DivideFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *operator^(const MEDCouplingFieldDouble& other) const;
    MEDCOUPLING_EXPORT const MEDCouplingFieldDouble &operator^=(const MEDCouplingFieldDouble& other);
    MEDCOUPLING_EXPORT static MEDCouplingFieldDouble *PowFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCOUPLING_EXPORT static std::string WriteVTK(const std::string& fileName, const std::vector<const MEDCouplingFieldDouble *>& fs, bool isBinary=true);
  public:
    MEDCOUPLING_EXPORT const MEDCouplingTimeDiscretization *getTimeDiscretizationUnderGround() const { return _time_discr; }
    MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization *getTimeDiscretizationUnderGround() { return _time_discr; }
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
  private:
    MEDCouplingFieldDouble(TypeOfField type, TypeOfTimeDiscretization td);
    MEDCouplingFieldDouble(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td);
    MEDCouplingFieldDouble(const MEDCouplingFieldDouble& other, bool deepCopy);
    MEDCouplingFieldDouble(NatureOfField n, MEDCouplingTimeDiscretization *td, MEDCouplingFieldDiscretization *type);
    ~MEDCouplingFieldDouble();
  private:
    MEDCouplingTimeDiscretization *_time_discr;
  };
}

#endif
