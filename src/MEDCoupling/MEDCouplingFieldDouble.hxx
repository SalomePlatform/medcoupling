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
    bool isEqual(const MEDCouplingField *other, double meshPrec, double valsPrec) const;
    bool areCompatibleForMerge(const MEDCouplingField *other) const;
    bool areStrictlyCompatible(const MEDCouplingField *other) const;
    bool areCompatibleForMul(const MEDCouplingField *other) const;
    void renumberCells(const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception);
    void renumberCellsWithoutMesh(const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception);
    void renumberNodes(const int *old2NewBg) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *clone(bool recDeepCpy) const;
    MEDCouplingFieldDouble *cloneWithMesh(bool recDeepCpy) const;
    MEDCouplingFieldDouble *buildNewTimeReprFromThis(TypeOfTimeDiscretization td, bool deepCpy) const;
    TypeOfTimeDiscretization getTimeDiscretization() const;
    void checkCoherency() const throw(INTERP_KERNEL::Exception);
    NatureOfField getNature() const { return _nature; }
    void setNature(NatureOfField nat) throw(INTERP_KERNEL::Exception);
    void setTime(double val, int dt, int it) { _time_discr->setTime(val,dt,it); }
    void setStartTime(double val, int dt, int it) { _time_discr->setStartTime(val,dt,it); }
    void setEndTime(double val, int dt, int it) { _time_discr->setEndTime(val,dt,it); }
    double getTime(int& dt, int& it) const { return _time_discr->getTime(dt,it); }
    double getStartTime(int& dt, int& it) const { return _time_discr->getStartTime(dt,it); }
    double getEndTime(int& dt, int& it) const { return _time_discr->getEndTime(dt,it); }
    double getIJ(int tupleId, int compoId) const { return getArray()->getIJ(tupleId,compoId); }
    double getIJK(int cellId, int nodeIdInCell, int compoId) const;
    void setArray(DataArrayDouble *array);
    void setEndArray(DataArrayDouble *array);
    DataArrayDouble *getArray() const { return _time_discr->getArray(); }
    DataArrayDouble *getEndArray() const { return _time_discr->getEndArray(); }
    double accumulate(int compId) const;
    void accumulate(double *res) const;
    double normL1(int compId, bool isWAbs) const throw(INTERP_KERNEL::Exception);
    void normL1(bool isWAbs, double *res) const throw(INTERP_KERNEL::Exception);
    double normL2(int compId, bool isWAbs) const throw(INTERP_KERNEL::Exception);
    void normL2(bool isWAbs, double *res) const throw(INTERP_KERNEL::Exception);
    double integral(int compId, bool isWAbs) const throw(INTERP_KERNEL::Exception);
    void integral(bool isWAbs, double *res) const throw(INTERP_KERNEL::Exception);
    void getValueOnPos(int i, int j, int k, double *res) const throw(INTERP_KERNEL::Exception);
    void getValueOn(const double *spaceLoc, double *res) const throw(INTERP_KERNEL::Exception);
    void getValueOn(const double *spaceLoc, double time, double *res) const throw(INTERP_KERNEL::Exception);
    //! \b temporary
    void applyLin(double a, double b, int compoId);
    void applyFunc(int nbOfComp, FunctionToEvaluate func);
    void applyFunc(int nbOfComp, const char *func);
    void applyFunc(const char *func);
    int getNumberOfComponents() const;
    int getNumberOfTuples() const throw(INTERP_KERNEL::Exception);
    void updateTime();
    //
    void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    void getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const;
    void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    void resizeForUnserialization(const std::vector<int>& tinyInfoI, DataArrayInt *&dataInt, std::vector<DataArrayDouble *>& arrays);
    void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS);
    void serialize(DataArrayInt *&dataInt, std::vector<DataArrayDouble *>& arrays) const;
    bool mergeNodes(double eps) throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldDouble *mergeFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
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
    MEDCouplingFieldDouble(NatureOfField n, MEDCouplingTimeDiscretization *td, TypeOfField type);
    ~MEDCouplingFieldDouble();
  private:
    NatureOfField _nature;
    MEDCouplingTimeDiscretization *_time_discr;
  };
}

#endif
