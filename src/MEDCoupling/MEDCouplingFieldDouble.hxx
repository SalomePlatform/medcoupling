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

#ifndef __PARAMEDMEM_MEDCOUPLINGFIELDDOUBLE_HXX__
#define __PARAMEDMEM_MEDCOUPLINGFIELDDOUBLE_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingField.hxx"
#include "MEDCouplingTimeDiscretization.hxx"
#include "MEDCouplingNatureOfField.hxx"
#include "MEDCouplingMemArray.hxx"

namespace ParaMEDMEM
{
  class MEDCOUPLING_EXPORT MEDCouplingFieldDouble : public MEDCouplingField
  {
  public:
    static MEDCouplingFieldDouble *New(TypeOfField type, TypeOfTimeDiscretization td=NO_TIME);
    void copyTinyStringsFrom(const MEDCouplingFieldDouble *other) throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const;
    std::string advancedRepr() const;
    bool isEqual(const MEDCouplingField *other, double meshPrec, double valsPrec) const;
    bool isEqualWithoutConsideringStr(const MEDCouplingField *other, double meshPrec, double valsPrec) const;
    bool areCompatibleForMerge(const MEDCouplingField *other) const;
    bool areStrictlyCompatible(const MEDCouplingField *other) const;
    bool areCompatibleForMul(const MEDCouplingField *other) const;
    void renumberCells(const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception);
    void renumberCellsWithoutMesh(const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception);
    void renumberNodes(const int *old2NewBg) throw(INTERP_KERNEL::Exception);
    void renumberNodesWithoutMesh(const int *old2NewBg) throw(INTERP_KERNEL::Exception);
    DataArrayInt *getIdsInRange(double vmin, double vmax) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *buildSubPart(const DataArrayInt *part) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *buildSubPart(const int *partBg, const int *partEnd) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *clone(bool recDeepCpy) const;
    MEDCouplingFieldDouble *cloneWithMesh(bool recDeepCpy) const;
    MEDCouplingFieldDouble *buildNewTimeReprFromThis(TypeOfTimeDiscretization td, bool deepCpy) const;
    TypeOfTimeDiscretization getTimeDiscretization() const;
    void checkCoherency() const throw(INTERP_KERNEL::Exception);
    NatureOfField getNature() const { return _nature; }
    void setNature(NatureOfField nat) throw(INTERP_KERNEL::Exception);
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
    DataArrayDouble *getArray() const { return _time_discr->getArray(); }
    DataArrayDouble *getEndArray() const { return _time_discr->getEndArray(); }
    double accumulate(int compId) const;
    void accumulate(double *res) const;
    double getMaxValue() const throw(INTERP_KERNEL::Exception);
    double getMaxValue2(DataArrayInt*& tupleIds) const throw(INTERP_KERNEL::Exception);
    double getMinValue() const throw(INTERP_KERNEL::Exception);
    double getMinValue2(DataArrayInt*& tupleIds) const throw(INTERP_KERNEL::Exception);
    double getAverageValue() const throw(INTERP_KERNEL::Exception);
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
    //! \b temporary
    void applyLin(double a, double b, int compoId);
    MEDCouplingFieldDouble &operator=(double value) throw(INTERP_KERNEL::Exception);
    void fillFromAnalytic(int nbOfComp, FunctionToEvaluate func) throw(INTERP_KERNEL::Exception);
    void fillFromAnalytic(int nbOfComp, const char *func) throw(INTERP_KERNEL::Exception);
    void applyFunc(int nbOfComp, FunctionToEvaluate func);
    void applyFunc(int nbOfComp, double val);
    void applyFunc(int nbOfComp, const char *func);
    void applyFunc(const char *func);
    void applyFuncFast32(const char *func) throw(INTERP_KERNEL::Exception);
    void applyFuncFast64(const char *func) throw(INTERP_KERNEL::Exception);
    int getNumberOfComponents() const throw(INTERP_KERNEL::Exception);
    int getNumberOfTuples() const throw(INTERP_KERNEL::Exception);
    int getNumberOfValues() const throw(INTERP_KERNEL::Exception);
    void updateTime();
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
    bool mergeNodes(double eps) throw(INTERP_KERNEL::Exception);
    bool zipCoords() throw(INTERP_KERNEL::Exception);
    bool zipConnectivity(int compType) throw(INTERP_KERNEL::Exception);
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
    static MEDCouplingFieldDouble *mergeFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    static MEDCouplingFieldDouble *dotFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCouplingFieldDouble *dot(const MEDCouplingFieldDouble& other) const { return dotFields(this,&other); }
    static MEDCouplingFieldDouble *crossProductFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCouplingFieldDouble *crossProduct(const MEDCouplingFieldDouble& other) const { return crossProductFields(this,&other); }
    static MEDCouplingFieldDouble *maxFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCouplingFieldDouble *max(const MEDCouplingFieldDouble& other) const { return maxFields(this,&other); }
    static MEDCouplingFieldDouble *minFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCouplingFieldDouble *min(const MEDCouplingFieldDouble& other) const { return minFields(this,&other); }
    MEDCouplingFieldDouble *operator+(const MEDCouplingFieldDouble& other) const { return addFields(this,&other); }
    const MEDCouplingFieldDouble &operator+=(const MEDCouplingFieldDouble& other);
    static MEDCouplingFieldDouble *addFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCouplingFieldDouble *operator-(const MEDCouplingFieldDouble& other) const { return substractFields(this,&other); }
    const MEDCouplingFieldDouble &operator-=(const MEDCouplingFieldDouble& other);
    static MEDCouplingFieldDouble *substractFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCouplingFieldDouble *operator*(const MEDCouplingFieldDouble& other) const { return multiplyFields(this,&other); }
    const MEDCouplingFieldDouble &operator*=(const MEDCouplingFieldDouble& other);
    static MEDCouplingFieldDouble *multiplyFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCouplingFieldDouble *operator/(const MEDCouplingFieldDouble& other) const { return divideFields(this,&other); }
    const MEDCouplingFieldDouble &operator/=(const MEDCouplingFieldDouble& other);
    static MEDCouplingFieldDouble *divideFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
  private:
    MEDCouplingFieldDouble(TypeOfField type, TypeOfTimeDiscretization td);
    MEDCouplingFieldDouble(const MEDCouplingFieldDouble& other, bool deepCpy);
    MEDCouplingFieldDouble(NatureOfField n, MEDCouplingTimeDiscretization *td, MEDCouplingFieldDiscretization *type);
    ~MEDCouplingFieldDouble();
  private:
    NatureOfField _nature;
    MEDCouplingTimeDiscretization *_time_discr;
  };
}

#endif
