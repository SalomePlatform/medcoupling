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
// Author : Anthony Geay (EDF R&D)

#include "MEDCouplingPartDefinition.hxx"

using namespace ParaMEDMEM;

PartDefinition *PartDefinition::New(int start, int stop, int step)
{
  return SlicePartDefinition::New(start,stop,step);
}

PartDefinition *PartDefinition::New(DataArrayInt *listOfIds)
{
  return DataArrayPartDefinition::New(listOfIds);
}

PartDefinition::~PartDefinition()
{
}

DataArrayPartDefinition *DataArrayPartDefinition::New(DataArrayInt *listOfIds)
{
  return new DataArrayPartDefinition(listOfIds);
}

int DataArrayPartDefinition::getNumberOfElems() const
{
  checkInternalArrayOK();
  return _arr->getNumberOfTuples();
}

PartDefinition *DataArrayPartDefinition::operator+(const PartDefinition& other) const
{
  const PartDefinition *otherPt(&other);
  if(!otherPt)
    throw INTERP_KERNEL::Exception("DataArrayPartDefinition::operator+ : NULL input !");
  const DataArrayPartDefinition *other1(dynamic_cast<const DataArrayPartDefinition *>(otherPt));
  if(other1)
    return add1(other1);
  const SlicePartDefinition *other2(dynamic_cast<const SlicePartDefinition *>(otherPt));
  if(other2)
    return add2(other2);
  throw INTERP_KERNEL::Exception("DataArrayPartDefinition::operator+ : unrecognized type in input !");
}

std::string DataArrayPartDefinition::getRepr() const
{
  std::ostringstream oss; oss << "DataArray Part : ";
  const DataArrayInt *arr(_arr);
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
  checkCoherency();
  other->checkCoherency();
  const SlicePartDefinition *spd(dynamic_cast<const SlicePartDefinition *>(other));
  if(spd)
    {//special case for optim
      int a(0),b(0),c(0);
      spd->getSlice(a,b,c);
      if(c==1)
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr(DataArrayInt::New());
          arr->alloc(_arr->getNumberOfTuples(),1);
          std::transform(_arr->begin(),_arr->end(),arr->getPointer(),std::bind2nd(std::plus<int>(),a));
          return DataArrayPartDefinition::New(arr);
        }
    }
  //
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr1(other->toDAI());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr2(arr1->selectByTupleIdSafe(_arr->begin(),_arr->end()));
  return DataArrayPartDefinition::New(arr2);
}

void DataArrayPartDefinition::checkCoherency() const
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
  checkCoherency();
  int a(0),b(0),c(0);
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

DataArrayInt *DataArrayPartDefinition::toDAI() const
{
  checkInternalArrayOK();
  const DataArrayInt *arr(_arr);
  DataArrayInt *arr2(const_cast<DataArrayInt *>(arr));
  arr2->incrRef();
  return arr2;
}

DataArrayPartDefinition::DataArrayPartDefinition(DataArrayInt *listOfIds)
{
  CheckInternalArrayOK(listOfIds);
  _arr=listOfIds;
  _arr->incrRef();
}

void DataArrayPartDefinition::checkInternalArrayOK() const
{
  CheckInternalArrayOK(_arr);
}

void DataArrayPartDefinition::CheckInternalArrayOK(const DataArrayInt *listOfIds)
{
  if(!listOfIds || !listOfIds->isAllocated() || listOfIds->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayPartDefinition::CheckInternalArrayOK : Input list must be not null allocated and with one components !");
}

void DataArrayPartDefinition::updateTime() const
{
  if((const DataArrayInt *)_arr)
    updateTimeWith(*_arr);
}

std::size_t DataArrayPartDefinition::getHeapMemorySizeWithoutChildren() const
{
  return sizeof(DataArrayPartDefinition);
}

std::vector<const BigMemoryObject *> DataArrayPartDefinition::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret(1,(const DataArrayInt *)_arr);
  return ret;
}

DataArrayPartDefinition *DataArrayPartDefinition::add1(const DataArrayPartDefinition *other) const
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> a1(toDAI()),a2(other->toDAI());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> a3(DataArrayInt::Aggregate(a1,a2,0));
  a3->sort();
  return DataArrayPartDefinition::New(a3);
}

DataArrayPartDefinition *DataArrayPartDefinition::add2(const SlicePartDefinition *other) const
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> a1(toDAI()),a2(other->toDAI());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> a3(DataArrayInt::Aggregate(a1,a2,0));
  a3->sort();
  return DataArrayPartDefinition::New(a3);
}

DataArrayPartDefinition::~DataArrayPartDefinition()
{
}

SlicePartDefinition *SlicePartDefinition::New(int start, int stop, int step)
{
  return new SlicePartDefinition(start,stop,step);
}

DataArrayInt *SlicePartDefinition::toDAI() const
{
  return DataArrayInt::Range(_start,_stop,_step);
}

int SlicePartDefinition::getNumberOfElems() const
{
  return DataArray::GetNumberOfItemGivenBES(_start,_stop,_step,"SlicePartDefinition::getNumberOfElems");
}

PartDefinition *SlicePartDefinition::operator+(const PartDefinition& other) const
{
  const PartDefinition *otherPt(&other);
  if(!otherPt)
    throw INTERP_KERNEL::Exception("DataArrayPartDefinition::operator+ : NULL input !");
  const DataArrayPartDefinition *other1(dynamic_cast<const DataArrayPartDefinition *>(otherPt));
  if(other1)
    return add1(other1);
  const SlicePartDefinition *other2(dynamic_cast<const SlicePartDefinition *>(otherPt));
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
  checkCoherency();
  other->checkCoherency();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr(other->toDAI());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr1(arr->selectByTupleId2(_start,_stop,_step));
  return DataArrayPartDefinition::New(arr1);
}

/*!
 * Do nothing it is not a bug.
 */
void SlicePartDefinition::checkCoherency() const
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

std::string SlicePartDefinition::getRepr() const
{
  std::ostringstream oss;
  oss << "Slice is defined with : start=" << _start << " stop=" << _stop << " step=" << _step;
  return oss.str();
}

int SlicePartDefinition::getEffectiveStop() const
{
  int nbElems(DataArray::GetNumberOfItemGivenBES(_start,_stop,_step,"SlicePartDefinition::getEffectiveStop"));
  return _start+nbElems*_step;
}

void SlicePartDefinition::getSlice(int& start, int& stop, int& step) const
{
  start=_start;
  stop=_stop;
  step=_step;
}

SlicePartDefinition::SlicePartDefinition(int start, int stop, int step):_start(start),_stop(stop),_step(step)
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
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> a1(toDAI()),a2(other->toDAI());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> a3(DataArrayInt::Aggregate(a1,a2,0));
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
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> a1(toDAI()),a2(other->toDAI());
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> a3(DataArrayInt::Aggregate(a1,a2,0));
      a3->sort();
      return DataArrayPartDefinition::New(a3);
    }
}

SlicePartDefinition::~SlicePartDefinition()
{
}
