// Copyright (C) 2016-2025  CEA, EDF
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

#pragma once

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
    MEDCOUPLING_EXPORT typename Traits<T>::FieldType *buildSubPart(const DataArrayIdType *part) const;
    MEDCOUPLING_EXPORT typename Traits<T>::FieldType *buildSubPart(const mcIdType *partBg, const mcIdType *partEnd) const;
    MEDCOUPLING_EXPORT typename Traits<T>::FieldType *buildSubPartRange(mcIdType begin, mcIdType end, mcIdType step) const;
    void setArray(typename Traits<T>::ArrayType *array) { _time_discr->setArray(array,this); }
    void setEndArray(typename Traits<T>::ArrayType *array) { _time_discr->setEndArray(array,this); }
    const typename Traits<T>::ArrayType *getArray() const { return _time_discr->getArray(); }
    typename Traits<T>::ArrayType *getArray() { return _time_discr->getArray(); }
    const typename Traits<T>::ArrayType *getEndArray() const { return _time_discr->getEndArray(); }
    typename Traits<T>::ArrayType *getEndArray() { return _time_discr->getEndArray(); }
    void setArrays(const std::vector<typename Traits<T>::ArrayType *>& arrs) { _time_discr->setArrays(arrs,this); }
    std::vector<typename Traits<T>::ArrayType *> getArrays() const { std::vector<typename Traits<T>::ArrayType *> ret; _time_discr->getArrays(ret); return ret; }
    void setTimeUnit(const std::string& unit) { _time_discr->setTimeUnit(unit); }
    std::string getTimeUnit() const { return _time_discr->getTimeUnit(); }
    void setTimeTolerance(double val) { _time_discr->setTimeTolerance(val); }
    double getTimeTolerance() const { return _time_discr->getTimeTolerance(); }
    void setIteration(int it) { _time_discr->setIteration(it); }
    void setEndIteration(int it) { _time_discr->setEndIteration(it); }
    void setOrder(int order) { _time_discr->setOrder(order); }
    void setEndOrder(int order) { _time_discr->setEndOrder(order); }
    void setTimeValue(double val) { _time_discr->setTimeValue(val); }
    void setEndTimeValue(double val) { _time_discr->setEndTimeValue(val); }
    void setTime(double val, int iteration, int order) { _time_discr->setTime(val,iteration,order); }
    MEDCOUPLING_EXPORT void synchronizeTimeWithMesh();
    void setStartTime(double val, int iteration, int order) { _time_discr->setStartTime(val,iteration,order); }
    void setEndTime(double val, int iteration, int order) { _time_discr->setEndTime(val,iteration,order); }
    double getTime(int& iteration, int& order) const { return _time_discr->getTime(iteration,order); }
    double getStartTime(int& iteration, int& order) const { return _time_discr->getStartTime(iteration,order); }
    double getEndTime(int& iteration, int& order) const { return _time_discr->getEndTime(iteration,order); }
    T getIJ(mcIdType tupleId, std::size_t compoId) const { return getArray()->getIJ(tupleId,compoId); }
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
    MEDCOUPLING_EXPORT void renumberCells(const mcIdType *old2NewBg, bool check=true);
    MEDCOUPLING_EXPORT void renumberCellsWithoutMesh(const mcIdType *old2NewBg, bool check=true);
    //
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation(std::vector<mcIdType>& tinyInfo) const;
    MEDCOUPLING_EXPORT void getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const;
    MEDCOUPLING_EXPORT void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    MEDCOUPLING_EXPORT void resizeForUnserialization(const std::vector<mcIdType>& tinyInfoI, DataArrayIdType *&dataInt, std::vector<typename Traits<T>::ArrayType *>& arrays);
    MEDCOUPLING_EXPORT void checkForUnserialization(const std::vector<mcIdType>& tinyInfoI, const DataArrayIdType *dataInt, const std::vector<typename Traits<T>::ArrayType *>& arrays);
    MEDCOUPLING_EXPORT void finishUnserialization(const std::vector<mcIdType>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS);
    MEDCOUPLING_EXPORT void serialize(DataArrayIdType *&dataInt, std::vector<typename Traits<T>::ArrayType *>& arrays) const;
    MEDCOUPLING_EXPORT const MEDCouplingTimeDiscretizationTemplate<T> *timeDiscrSafe() const;
  protected:
    MEDCouplingTimeDiscretizationTemplate<T> *timeDiscrSafe();
  protected:
    MEDCouplingTimeDiscretizationTemplate<T> *_time_discr;
  };
}
