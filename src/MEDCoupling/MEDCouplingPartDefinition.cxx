// Copyright (C) 2007-2024  CEA, EDF
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

#include "MEDCouplingPartDefinition.hxx"
#include "MCType.hxx"
#include "MCAuto.hxx"
#include "MCIdType.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "MEDCouplingMemArray.hxx"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <sstream>
#include <vector>
#include <string>

using namespace MEDCoupling;

PartDefinition *PartDefinition::New(mcIdType start, mcIdType stop, mcIdType step)
{
  return SlicePartDefinition::New(start,stop,step);
}

PartDefinition *PartDefinition::New(DataArrayIdType *listOfIds)
{
  return DataArrayPartDefinition::New(listOfIds);
}

PartDefinition *PartDefinition::Unserialize(std::vector<mcIdType>& tinyInt, std::vector< MCAuto<DataArrayIdType> >& bigArraysI)
{
  if(tinyInt.empty())
    {
      MCAuto<PartDefinition> ret(DataArrayPartDefinition::New(bigArraysI.back()));
      bigArraysI.pop_back();
      return ret.retn();
    }
  else if(tinyInt.size()==3)
    {
      MCAuto<PartDefinition> ret(SlicePartDefinition::New(tinyInt[0],tinyInt[1],tinyInt[2]));
      tinyInt.erase(tinyInt.begin(),tinyInt.begin()+3);
      return ret.retn();
    }
  else
    throw INTERP_KERNEL::Exception("PartDefinition::Unserialize");
}

PartDefinition::~PartDefinition()
= default;

DataArrayPartDefinition *DataArrayPartDefinition::New(DataArrayIdType *listOfIds)
{
  return new DataArrayPartDefinition(listOfIds);
}

bool DataArrayPartDefinition::isEqual(const PartDefinition *other, std::string& what) const
{
  if(!other)
    {
      what="DataArrayPartDefinition::isEqual : other is null, this is not null !";
      return false;
    }
  const auto *otherC(dynamic_cast<const DataArrayPartDefinition *>(other));
  if(!otherC)
    {
      what="DataArrayPartDefinition::isEqual : other is not DataArrayPartDefinition !";
      return false;
    }
  const DataArrayIdType *arr0(_arr),*arr1(otherC->_arr);
  if(!arr0 && !arr1)
    return true;
  if((arr0 && !arr1) || (!arr0 && arr1))
    {
      what="DataArrayPartDefinition::isEqual : array is not defined both in other and this !";
      return false;
    }
  std::string what1;
  bool const ret(arr0->isEqualIfNotWhy(*arr1,what1));
  if(!ret)
    {
      what=std::string("DataArrayPartDefinition::isEqual : arrays are not equal :\n")+what1;
      return false;
    }
  return true;
}

DataArrayPartDefinition *DataArrayPartDefinition::deepCopy() const
{
  const DataArrayIdType *arr(_arr);
  if(!arr)
    throw INTERP_KERNEL::Exception("DataArrayPartDefinition::deepCopy : array is null !");
  return DataArrayPartDefinition::New(const_cast<DataArrayIdType *>(arr));
}

mcIdType DataArrayPartDefinition::getNumberOfElems() const
{
  checkInternalArrayOK();
  return _arr->getNumberOfTuples();
}

PartDefinition *DataArrayPartDefinition::operator+(const PartDefinition& other) const
{
  const PartDefinition *otherPt(&other);
  if(!otherPt)
    throw INTERP_KERNEL::Exception("DataArrayPartDefinition::operator+ : NULL input !");
  const auto *other1(dynamic_cast<const DataArrayPartDefinition *>(otherPt));
  if(other1)
    return add1(other1);
  const auto *other2(dynamic_cast<const SlicePartDefinition *>(otherPt));
  if(other2)
    return add2(other2);
  throw INTERP_KERNEL::Exception("DataArrayPartDefinition::operator+ : unrecognized type in input !");
}

std::string DataArrayPartDefinition::getRepr() const
{
  std::ostringstream oss; oss << "DataArray Part : ";
  const DataArrayIdType *arr(_arr);
  if(arr)
    arr->reprQuickOverview(oss);
  else
    oss << "No Data !";
  return oss.str();
}

/*!
 * This method operates FoG where F is \a this and G is \a other.
 * Example : if \a other is SlicePart(4,14,1) and if \a this is DataArrayPartDefinition([0,1,2,3,6,7,8,9]) -> DataArrayPartDefinition([4,5,6,7,11,12,13]) will be returned
 */
PartDefinition *DataArrayPartDefinition::composeWith(const PartDefinition *other) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayPartDefinition::composeWith : input PartDef must be not NULL !");
  checkConsistencyLight();
  other->checkConsistencyLight();
  const auto *spd(dynamic_cast<const SlicePartDefinition *>(other));
  if(spd)
    {//special case for optim
      mcIdType a(0),b(0),c(0);
      spd->getSlice(a,b,c);
      if(c==1)
        {
          MCAuto<DataArrayIdType> arr(DataArrayIdType::New());
          arr->alloc(_arr->getNumberOfTuples(),1);
          std::transform(_arr->begin(),_arr->end(),arr->getPointer(),std::bind(std::plus<mcIdType>(),std::placeholders::_1,a));
          return DataArrayPartDefinition::New(arr);
        }
    }
  //
  MCAuto<DataArrayIdType> arr1(other->toDAI());
  MCAuto<DataArrayIdType> arr2(arr1->selectByTupleIdSafe(_arr->begin(),_arr->end()));
  return DataArrayPartDefinition::New(arr2);
}

void DataArrayPartDefinition::checkConsistencyLight() const
{
  CheckInternalArrayOK(_arr);
}

/*!
 * This method tries to simplify \a this if possible.
 * 
 * \return a new reference (equal to this) to be decrRefed.
 */
PartDefinition *DataArrayPartDefinition::tryToSimplify() const
{
  checkConsistencyLight();
  mcIdType a(0),b(0),c(0);
  if(_arr->isRange(a,b,c))
    {
      return SlicePartDefinition::New(a,b,c);
    }
  else
    {
      PartDefinition *ret(const_cast<DataArrayPartDefinition *>(this));
      ret->incrRef();
      return ret;
    }
}

void DataArrayPartDefinition::serialize(std::vector<mcIdType>&  /*tinyInt*/, std::vector< MCAuto<DataArrayIdType> >& bigArraysI) const
{
  bigArraysI.push_back(_arr);
}

DataArrayIdType *DataArrayPartDefinition::toDAI() const
{
  checkInternalArrayOK();
  const DataArrayIdType *arr(_arr);
  auto *arr2(const_cast<DataArrayIdType *>(arr));
  arr2->incrRef();
  return arr2;
}

DataArrayPartDefinition::DataArrayPartDefinition(DataArrayIdType *listOfIds)
{
  CheckInternalArrayOK(listOfIds);
  _arr=listOfIds;
  _arr->incrRef();
}

void DataArrayPartDefinition::checkInternalArrayOK() const
{
  CheckInternalArrayOK(_arr);
}

void DataArrayPartDefinition::CheckInternalArrayOK(const DataArrayIdType *listOfIds)
{
  if(!listOfIds || !listOfIds->isAllocated() || listOfIds->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayPartDefinition::CheckInternalArrayOK : Input list must be not null allocated and with one components !");
}

void DataArrayPartDefinition::updateTime() const
{
  if((const DataArrayIdType *)_arr)
    updateTimeWith(*_arr);
}

std::size_t DataArrayPartDefinition::getHeapMemorySizeWithoutChildren() const
{
  return sizeof(DataArrayPartDefinition);
}

std::vector<const BigMemoryObject *> DataArrayPartDefinition::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret(1,(const DataArrayIdType *)_arr);
  return ret;
}

DataArrayPartDefinition *DataArrayPartDefinition::add1(const DataArrayPartDefinition *other) const
{
  MCAuto<DataArrayIdType> a1(toDAI()),a2(other->toDAI());
  MCAuto<DataArrayIdType> a3(DataArrayIdType::Aggregate(a1,a2,0));
  a3->sort();
  return DataArrayPartDefinition::New(a3);
}

DataArrayPartDefinition *DataArrayPartDefinition::add2(const SlicePartDefinition *other) const
{
  MCAuto<DataArrayIdType> a1(toDAI()),a2(other->toDAI());
  MCAuto<DataArrayIdType> a3(DataArrayIdType::Aggregate(a1,a2,0));
  a3->sort();
  return DataArrayPartDefinition::New(a3);
}

DataArrayPartDefinition::~DataArrayPartDefinition()
= default;

SlicePartDefinition *SlicePartDefinition::New(mcIdType start, mcIdType stop, mcIdType step)
{
  return new SlicePartDefinition(start,stop,step);
}

bool SlicePartDefinition::isEqual(const PartDefinition *other, std::string& what) const
{
  if(!other)
    {
      what="SlicePartDefinition::isEqual : other is null, this is not null !";
      return false;
    }
  const auto *otherC(dynamic_cast<const SlicePartDefinition *>(other));
  if(!otherC)
    {
      what="SlicePartDefinition::isEqual : other is not SlicePartDefinition !";
      return false;
    }
  bool const ret((_start==otherC->_start) && (_stop==otherC->_stop) && (_step==otherC->_step));
  if(!ret)
    {
      what="SlicePartDefinition::isEqual : values are not the same !";
      return false;
    }
  return true;
}

SlicePartDefinition *SlicePartDefinition::deepCopy() const
{
  return SlicePartDefinition::New(_start,_stop,_step);
}

DataArrayIdType *SlicePartDefinition::toDAI() const
{
  return DataArrayIdType::Range(_start,_stop,_step);
}

mcIdType SlicePartDefinition::getNumberOfElems() const
{
  return DataArray::GetNumberOfItemGivenBES(_start,_stop,_step,"SlicePartDefinition::getNumberOfElems");
}

PartDefinition *SlicePartDefinition::operator+(const PartDefinition& other) const
{
  const PartDefinition *otherPt(&other);
  if(!otherPt)
    throw INTERP_KERNEL::Exception("DataArrayPartDefinition::operator+ : NULL input !");
  const auto *other1(dynamic_cast<const DataArrayPartDefinition *>(otherPt));
  if(other1)
    return add1(other1);
  const auto *other2(dynamic_cast<const SlicePartDefinition *>(otherPt));
  if(other2)
    return add2(other2);
  throw INTERP_KERNEL::Exception("SlicePartDefinition::operator+ : unrecognized type in input !");
}

/*!
 * This method operates FoG where F is \a this and G is \a other.
 * Example : if \a this is SlicePart(4,6,1) and if \a other is DataArrayPartDefinition([12,13,17,18,22,28,34,44]) -> DataArrayPartDefinition([22,28]) will be returned
 */
PartDefinition *SlicePartDefinition::composeWith(const PartDefinition *other) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("SlicePartDefinition::composeWith : input PartDef must be not NULL !");
  checkConsistencyLight();
  other->checkConsistencyLight();
  MCAuto<DataArrayIdType> arr(other->toDAI());
  MCAuto<DataArrayIdType> arr1(arr->selectByTupleIdSafeSlice(_start,_stop,_step));
  return DataArrayPartDefinition::New(arr1);
}

/*!
 * Do nothing it is not a bug.
 */
void SlicePartDefinition::checkConsistencyLight() const
{
}

/*!
 * Return \a this (because it cannot be simplified)
 * 
 * \return a new reference (equal to this) to be decrRefed.
 */
PartDefinition *SlicePartDefinition::tryToSimplify() const
{
  PartDefinition *ret(const_cast<SlicePartDefinition *>(this));
  ret->incrRef();
  return ret;
}

void SlicePartDefinition::serialize(std::vector<mcIdType>& tinyInt, std::vector< MCAuto<DataArrayIdType> >&  /*bigArraysI*/) const
{
  tinyInt.push_back(_start);
  tinyInt.push_back(_stop);
  tinyInt.push_back(_step);
}

std::string SlicePartDefinition::getRepr() const
{
  std::ostringstream oss;
  oss << "Slice is defined with : start=" << _start << " stop=" << _stop << " step=" << _step;
  return oss.str();
}

mcIdType SlicePartDefinition::getEffectiveStop() const
{
  mcIdType const nbElems(DataArray::GetNumberOfItemGivenBES(_start,_stop,_step,"SlicePartDefinition::getEffectiveStop"));
  return _start+nbElems*_step;
}

void SlicePartDefinition::getSlice(mcIdType& start, mcIdType& stop, mcIdType& step) const
{
  start=_start;
  stop=_stop;
  step=_step;
}

SlicePartDefinition::SlicePartDefinition(mcIdType start, mcIdType stop, mcIdType step):_start(start),_stop(stop),_step(step)
{
}

/*!
 * No child ! It is the leaf ! So no implementation.
 */
void SlicePartDefinition::updateTime() const
{
}

std::size_t SlicePartDefinition::getHeapMemorySizeWithoutChildren() const
{
  return sizeof(SlicePartDefinition);
}

std::vector<const BigMemoryObject *> SlicePartDefinition::getDirectChildrenWithNull() const
{
  return std::vector<const BigMemoryObject *>();
}

DataArrayPartDefinition *SlicePartDefinition::add1(const DataArrayPartDefinition *other) const
{
  MCAuto<DataArrayIdType> a1(toDAI()),a2(other->toDAI());
  MCAuto<DataArrayIdType> a3(DataArrayIdType::Aggregate(a1,a2,0));
  a3->sort();
  return DataArrayPartDefinition::New(a3);
}

PartDefinition *SlicePartDefinition::add2(const SlicePartDefinition *other) const
{
  if(_step==other->_step && getEffectiveStop()==other->_start)
    {
      return SlicePartDefinition::New(_start,other->_stop,_step);
    }
  else
    {
      MCAuto<DataArrayIdType> a1(toDAI()),a2(other->toDAI());
      MCAuto<DataArrayIdType> a3(DataArrayIdType::Aggregate(a1,a2,0));
      a3->sort();
      return DataArrayPartDefinition::New(a3);
    }
}

SlicePartDefinition::~SlicePartDefinition()
= default;
