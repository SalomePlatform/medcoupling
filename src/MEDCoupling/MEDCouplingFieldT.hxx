// Copyright (C) 2016  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (EDF R&D)

#ifndef __MEDCOUPLINGFIELDT_HXX__
#define __MEDCOUPLINGFIELDT_HXX__

#include "MEDCouplingField.hxx"
#include "MEDCouplingTraits.hxx"
#include "MEDCouplingTimeDiscretization.hxx"

#include <sstream>

namespace MEDCoupling
{
  template<class T>
  class MEDCouplingTimeDiscretizationTemplate;
  
  template<class T>
  class MEDCouplingFieldT : public MEDCouplingField
  {
  protected:
    MEDCouplingFieldT(const MEDCouplingFieldT<T>& other, bool deepCopy);
    MEDCouplingFieldT(const MEDCouplingField& other, MEDCouplingTimeDiscretizationTemplate<T> *timeDiscr, bool deepCopy=true);
    MEDCouplingFieldT(TypeOfField type, MEDCouplingTimeDiscretizationTemplate<T> *timeDiscr);
    MEDCouplingFieldT(MEDCouplingFieldDiscretization *type, NatureOfField n, MEDCouplingTimeDiscretizationTemplate<T> *timeDiscr);
    ~MEDCouplingFieldT();
  public:
    MEDCOUPLING_EXPORT TypeOfTimeDiscretization getTimeDiscretization() const;
    MEDCOUPLING_EXPORT virtual typename Traits<T>::FieldType *clone(bool recDeepCpy) const = 0;
    MEDCOUPLING_EXPORT void checkConsistencyLight() const;
    MEDCOUPLING_EXPORT typename Traits<T>::FieldType *cloneWithMesh(bool recDeepCpy) const;
    MEDCOUPLING_EXPORT typename Traits<T>::FieldType *buildSubPart(const DataArrayInt *part) const;
    MEDCOUPLING_EXPORT typename Traits<T>::FieldType *buildSubPart(const int *partBg, const int *partEnd) const;
    MEDCOUPLING_EXPORT typename Traits<T>::FieldType *buildSubPartRange(int begin, int end, int step) const;
    MEDCOUPLING_EXPORT void setArray(typename Traits<T>::ArrayType *array) { _time_discr->setArray(array,this); }
    MEDCOUPLING_EXPORT void setEndArray(typename Traits<T>::ArrayType *array) { _time_discr->setEndArray(array,this); }
    MEDCOUPLING_EXPORT const typename Traits<T>::ArrayType *getArray() const { return _time_discr->getArray(); }
    MEDCOUPLING_EXPORT typename Traits<T>::ArrayType *getArray() { return _time_discr->getArray(); }
    MEDCOUPLING_EXPORT const typename Traits<T>::ArrayType *getEndArray() const { return _time_discr->getEndArray(); }
    MEDCOUPLING_EXPORT typename Traits<T>::ArrayType *getEndArray() { return _time_discr->getEndArray(); }
    MEDCOUPLING_EXPORT void setArrays(const std::vector<typename Traits<T>::ArrayType *>& arrs) { _time_discr->setArrays(arrs,this); }
    MEDCOUPLING_EXPORT std::vector<typename Traits<T>::ArrayType *> getArrays() const { std::vector<typename Traits<T>::ArrayType *> ret; _time_discr->getArrays(ret); return ret; }
    MEDCOUPLING_EXPORT void setTimeUnit(const std::string& unit) { _time_discr->setTimeUnit(unit); }
    MEDCOUPLING_EXPORT std::string getTimeUnit() const { return _time_discr->getTimeUnit(); }
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
    MEDCOUPLING_EXPORT T getIJ(int tupleId, int compoId) const { return getArray()->getIJ(tupleId,compoId); }
    MEDCOUPLING_EXPORT virtual bool isEqual(const MEDCouplingFieldT<T> *other, double meshPrec, T valsPrec) const;
    MEDCOUPLING_EXPORT virtual bool isEqualIfNotWhy(const MEDCouplingFieldT<T> *other, double meshPrec, T valsPrec, std::string& reason) const;
    MEDCOUPLING_EXPORT virtual bool isEqualWithoutConsideringStr(const MEDCouplingFieldT<T> *other, double meshPrec, T valsPrec) const;
    MEDCOUPLING_EXPORT void copyTinyStringsFrom(const MEDCouplingField *other);
    MEDCOUPLING_EXPORT bool areStrictlyCompatible(const MEDCouplingField *other) const;
    MEDCOUPLING_EXPORT bool areStrictlyCompatibleForMulDiv(const MEDCouplingField *other) const;
    MEDCOUPLING_EXPORT std::string simpleRepr() const;
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
    // specific
    MEDCOUPLING_EXPORT bool areCompatibleForMul(const MEDCouplingField *other) const;
    MEDCOUPLING_EXPORT bool areCompatibleForDiv(const MEDCouplingField *other) const;
    MEDCOUPLING_EXPORT void copyTinyAttrFrom(const MEDCouplingFieldT<T> *other);
    MEDCOUPLING_EXPORT void copyAllTinyAttrFrom(const MEDCouplingFieldT<T> *other);
    MEDCOUPLING_EXPORT void renumberCells(const int *old2NewBg, bool check=true);
    MEDCOUPLING_EXPORT void renumberCellsWithoutMesh(const int *old2NewBg, bool check=true);
    //
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    MEDCOUPLING_EXPORT void getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const;
    MEDCOUPLING_EXPORT void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    MEDCOUPLING_EXPORT void resizeForUnserialization(const std::vector<int>& tinyInfoI, DataArrayInt *&dataInt, std::vector<typename Traits<T>::ArrayType *>& arrays);
    MEDCOUPLING_EXPORT void checkForUnserialization(const std::vector<int>& tinyInfoI, const DataArrayInt *dataInt, const std::vector<typename Traits<T>::ArrayType *>& arrays);
    MEDCOUPLING_EXPORT void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS);
    MEDCOUPLING_EXPORT void serialize(DataArrayInt *&dataInt, std::vector<typename Traits<T>::ArrayType *>& arrays) const;
    MEDCOUPLING_EXPORT const MEDCouplingTimeDiscretizationTemplate<T> *timeDiscrSafe() const;
  protected:
    MEDCouplingTimeDiscretizationTemplate<T> *timeDiscrSafe();
  protected:
    MEDCouplingTimeDiscretizationTemplate<T> *_time_discr;
  };
}

#endif
