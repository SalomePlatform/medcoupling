//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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
    bool isEqual(const MEDCouplingField *other, double meshPrec, double valsPrec) const;
    bool areCompatible(const MEDCouplingField *other) const;
    MEDCouplingFieldDouble *clone(bool recDeepCpy) const;
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
    void setArray(DataArrayDouble *array);
    DataArrayDouble *getArray() const { return _time_discr->getArray(); }
    double accumulate(int compId) const;
    double measureAccumulate(int compId, bool isWAbs) const;
    void getValueOn(const double *spaceLoc, double *res) const throw(INTERP_KERNEL::Exception);
    void getValueOn(const double *spaceLoc, double time, double *res) const throw(INTERP_KERNEL::Exception);
    //! \b temporary
    void applyLin(double a, double b, int compoId);
    void applyFunc(int nbOfComp, FunctionToEvaluate func);
    int getNumberOfComponents() const;
    int getNumberOfTuples() const throw(INTERP_KERNEL::Exception);
    void updateTime();
    //
    void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    void getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const;
    void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    void resizeForUnserialization(const std::vector<int>& tinyInfoI, std::vector<DataArrayDouble *>& arrays);
    void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS);
    void serialize(std::vector<DataArrayDouble *>& arrays) const;
    static MEDCouplingFieldDouble *mergeFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCouplingFieldDouble *operator+(const MEDCouplingFieldDouble& other) const { return addFields(this,&other); }
    static MEDCouplingFieldDouble *addFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCouplingFieldDouble *operator-(const MEDCouplingFieldDouble& other) const { return substractFields(this,&other); }
    static MEDCouplingFieldDouble *substractFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCouplingFieldDouble *operator*(const MEDCouplingFieldDouble& other) const { return multiplyFields(this,&other); }
    static MEDCouplingFieldDouble *multiplyFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCouplingFieldDouble *operator/(const MEDCouplingFieldDouble& other) const { return divideFields(this,&other); }
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
