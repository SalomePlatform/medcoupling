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

#include "MEDCouplingTimeDiscretization.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingMesh.hxx"

#include <cmath>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <functional>

using namespace ParaMEDMEM;

const double MEDCouplingTimeDiscretization::TIME_TOLERANCE_DFT=1.e-12;

const char MEDCouplingNoTimeLabel::EXCEPTION_MSG[]="MEDCouplingNoTimeLabel::setTime : no time info attached.";

const char MEDCouplingNoTimeLabel::REPR[]="No time label defined.";

const char MEDCouplingWithTimeStep::EXCEPTION_MSG[]="No data on this time.";

const char MEDCouplingWithTimeStep::REPR[]="One time label.";

const char MEDCouplingConstOnTimeInterval::EXCEPTION_MSG[]="No data on this time.";

const char MEDCouplingConstOnTimeInterval::REPR[]="Constant on a time interval.";

const char MEDCouplingTwoTimeSteps::EXCEPTION_MSG[]="No data on this time.";

const char MEDCouplingLinearTime::REPR[]="Linear time between 2 time steps.";

MEDCouplingTimeDiscretization *MEDCouplingTimeDiscretization::New(TypeOfTimeDiscretization type)
{
  switch(type)
  {
    case MEDCouplingNoTimeLabel::DISCRETIZATION:
      return new MEDCouplingNoTimeLabel;
    case MEDCouplingWithTimeStep::DISCRETIZATION:
      return new MEDCouplingWithTimeStep;
    case MEDCouplingConstOnTimeInterval::DISCRETIZATION:
      return new MEDCouplingConstOnTimeInterval;
    case MEDCouplingLinearTime::DISCRETIZATION:
      return new MEDCouplingLinearTime;
    default:
      throw INTERP_KERNEL::Exception("Time discretization not implemented yet");
  }
}

void MEDCouplingTimeDiscretization::copyTinyAttrFrom(const MEDCouplingTimeDiscretization& other)
{
  _time_tolerance=other._time_tolerance;
  _time_unit=other._time_unit;
}

void MEDCouplingTimeDiscretization::copyTinyStringsFrom(const MEDCouplingTimeDiscretization& other)
{
  _time_unit=other._time_unit;
  if(_array && other._array)
    _array->copyStringInfoFrom(*other._array);
}

void MEDCouplingTimeDiscretization::checkCoherency() const
{
  if(!_array)
    throw INTERP_KERNEL::Exception("Field invalid because no values set !");
  if(_time_tolerance<0.)
    throw INTERP_KERNEL::Exception("time tolerance is expected to be greater than 0. !");
}

void MEDCouplingTimeDiscretization::updateTime() const
{
  if(_array)
    updateTimeWith(*_array);
}

std::size_t MEDCouplingTimeDiscretization::getHeapMemorySizeWithoutChildren() const
{
  std::size_t ret(_time_unit.capacity());
  return ret;
}

std::vector<const BigMemoryObject *> MEDCouplingTimeDiscretization::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  ret.push_back(_array);
  return ret;
}

bool MEDCouplingTimeDiscretization::areCompatible(const MEDCouplingTimeDiscretization *other) const
{
  if(std::fabs(_time_tolerance-other->_time_tolerance)>1.e-16)
    return false;
  if(_array==0 && other->_array==0)
    return true;
  if(_array==0 || other->_array==0)
    return false;
  if(_array->getNumberOfComponents()!=other->_array->getNumberOfComponents())
    return false;
  return true;
}

bool MEDCouplingTimeDiscretization::areStrictlyCompatible(const MEDCouplingTimeDiscretization *other, std::string& reason) const
{
  std::ostringstream oss; oss.precision(15);
  if(_time_unit!=other->_time_unit)
    {
      oss << "Field discretizations differ : this time unit = \"" << _time_unit << "\" and other time unit = \"" << other->_time_unit << "\" !";
      reason=oss.str();
      return false;
    }
  if(std::fabs(_time_tolerance-other->_time_tolerance)>1.e-16)
    {
      oss << "Field discretizations differ : this time tolerance = \"" << _time_tolerance << "\" and other time tolerance = \"" << other->_time_tolerance << "\" !";
      reason=oss.str();
      return false;
    }
  if(_array==0 && other->_array==0)
    return true;
  if(_array==0 || other->_array==0)
    {
      reason="Field discretizations differ : Only one timediscretization between the two this and other has a DataArrayDouble for values defined";
      return false;
    }
  if(_array->getNumberOfComponents()!=other->_array->getNumberOfComponents())
    return false;
  if(_array->getNumberOfTuples()!=other->_array->getNumberOfTuples())
    return false;
  return true;
}

bool MEDCouplingTimeDiscretization::areCompatibleForMeld(const MEDCouplingTimeDiscretization *other) const
{
  if(std::fabs(_time_tolerance-other->_time_tolerance)>1.e-16)
    return false;
  if(_array==0 && other->_array==0)
    return true;
  if(_array==0 || other->_array==0)
    return false;
  if(_array->getNumberOfTuples()!=other->_array->getNumberOfTuples())
    return false;
  return true;
}

bool MEDCouplingTimeDiscretization::areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretization *other) const
{
  if(std::fabs(_time_tolerance-other->_time_tolerance)>1.e-16)
    return false;
  if(_array==0 && other->_array==0)
    return true;
  if(_array==0 || other->_array==0)
    return false;
  int nbC1=_array->getNumberOfComponents();
  int nbC2=other->_array->getNumberOfComponents();
  int nbMin=std::min(nbC1,nbC2);
  if(nbC1!=nbC2 && nbMin!=1)
    return false;
  return true;
}

bool MEDCouplingTimeDiscretization::areStrictlyCompatibleForDiv(const MEDCouplingTimeDiscretization *other) const
{
  if(std::fabs(_time_tolerance-other->_time_tolerance)>1.e-16)
    return false;
  if(_array==0 && other->_array==0)
    return true;
  if(_array==0 || other->_array==0)
    return false;
  int nbC1=_array->getNumberOfComponents();
  int nbC2=other->_array->getNumberOfComponents();
  if(nbC1!=nbC2 && nbC2!=1)
    return false;
  return true;
}

bool MEDCouplingTimeDiscretization::isEqualIfNotWhy(const MEDCouplingTimeDiscretization *other, double prec, std::string& reason) const
{
  if(!areStrictlyCompatible(other,reason))
    return false;
  if(_array==other->_array)
    return true;
  return _array->isEqualIfNotWhy(*other->_array,prec,reason);
}

bool MEDCouplingTimeDiscretization::isEqual(const MEDCouplingTimeDiscretization *other, double prec) const
{
  std::string reason;
  return isEqualIfNotWhy(other,prec,reason);
}

bool MEDCouplingTimeDiscretization::isEqualWithoutConsideringStr(const MEDCouplingTimeDiscretization *other, double prec) const
{
  std::string tmp;
  if(!areStrictlyCompatible(other,tmp))
    return false;
  if(_array==other->_array)
    return true;
  return _array->isEqualWithoutConsideringStr(*other->_array,prec);
}

MEDCouplingTimeDiscretization *MEDCouplingTimeDiscretization::buildNewTimeReprFromThis(TypeOfTimeDiscretization type, bool deepCpy) const
{
  MEDCouplingTimeDiscretization *ret=MEDCouplingTimeDiscretization::New(type);
  ret->setTimeUnit(getTimeUnit());
  const DataArrayDouble *arrSrc=getArray();
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr;
  if(arrSrc)
    arr=arrSrc->performCpy(deepCpy);
  ret->setArray(arr,0);
  return ret;
}

void MEDCouplingTimeDiscretization::getTinySerializationIntInformation(std::vector<int>& tinyInfo) const
{
  if(_array)
    {
      tinyInfo.push_back(_array->getNumberOfTuples());
      tinyInfo.push_back(_array->getNumberOfComponents());
    }
  else
    {
      tinyInfo.push_back(-1);
      tinyInfo.push_back(-1);
    }
}

void MEDCouplingTimeDiscretization::resizeForUnserialization(const std::vector<int>& tinyInfoI, std::vector<DataArrayDouble *>& arrays)
{
  arrays.resize(1);
  if(_array!=0)
    _array->decrRef();
  DataArrayDouble *arr=0;
  if(tinyInfoI[0]!=-1 && tinyInfoI[1]!=-1)
    {
      arr=DataArrayDouble::New();
      arr->alloc(tinyInfoI[0],tinyInfoI[1]);
    }
  _array=arr;
  arrays[0]=arr;
}

void MEDCouplingTimeDiscretization::checkForUnserialization(const std::vector<int>& tinyInfoI, const std::vector<DataArrayDouble *>& arrays)
{
  static const char MSG[]="MEDCouplingTimeDiscretization::checkForUnserialization : arrays in input is expected to have size one !";
  if(arrays.size()!=1)
    throw INTERP_KERNEL::Exception(MSG);
  if(_array!=0)
    _array->decrRef();
  _array=0;
  if(tinyInfoI[0]!=-1 && tinyInfoI[1]!=-1)
    {
      if(!arrays[0])
        throw INTERP_KERNEL::Exception(MSG);
      arrays[0]->checkNbOfTuplesAndComp(tinyInfoI[0],tinyInfoI[1],MSG);
      _array=arrays[0];
      _array->incrRef();
    }
}

void MEDCouplingTimeDiscretization::finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS)
{
  _time_tolerance=tinyInfoD[0];
  int nbOfCompo=_array->getNumberOfComponents();
  for(int i=0;i<nbOfCompo;i++)
    _array->setInfoOnComponent(i,tinyInfoS[i]);
}

void MEDCouplingTimeDiscretization::getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const
{
  tinyInfo.push_back(_time_tolerance);
}

void MEDCouplingTimeDiscretization::getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const
{
  int nbOfCompo=_array->getNumberOfComponents();
  for(int i=0;i<nbOfCompo;i++)
    tinyInfo.push_back(_array->getInfoOnComponent(i));
}

MEDCouplingTimeDiscretization::MEDCouplingTimeDiscretization():_time_tolerance(TIME_TOLERANCE_DFT),_array(0)
{
}

MEDCouplingTimeDiscretization::MEDCouplingTimeDiscretization(const MEDCouplingTimeDiscretization& other, bool deepCpy):_time_unit(other._time_unit),_time_tolerance(other._time_tolerance)
{
  if(other._array)
    _array=other._array->performCpy(deepCpy);
  else
    _array=0;
}

MEDCouplingTimeDiscretization::~MEDCouplingTimeDiscretization()
{
  if(_array)
    _array->decrRef();
}

void MEDCouplingTimeDiscretization::setArray(DataArrayDouble *array, TimeLabel *owner)
{
  if(array!=_array)
    {
      if(_array)
        _array->decrRef();
      _array=array;
      if(_array)
        _array->incrRef();
      if(owner)
        owner->declareAsNew();
    }
}

const DataArrayDouble *MEDCouplingTimeDiscretization::getEndArray() const
{
  throw INTERP_KERNEL::Exception("getEndArray not available for this type of time discretization !");
}

DataArrayDouble *MEDCouplingTimeDiscretization::getEndArray()
{
  throw INTERP_KERNEL::Exception("getEndArray not available for this type of time discretization !");
}

void MEDCouplingTimeDiscretization::setEndArray(DataArrayDouble *array, TimeLabel *owner)
{
  throw INTERP_KERNEL::Exception("setEndArray not available for this type of time discretization !");
}

void MEDCouplingTimeDiscretization::setArrays(const std::vector<DataArrayDouble *>& arrays, TimeLabel *owner)
{
  if(arrays.size()!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingTimeDiscretization::setArrays : number of arrays must be one.");
  setArray(arrays.back(),owner);
}

void MEDCouplingTimeDiscretization::getArrays(std::vector<DataArrayDouble *>& arrays) const
{
  arrays.resize(1);
  arrays[0]=_array;
}

bool MEDCouplingTimeDiscretization::isBefore(const MEDCouplingTimeDiscretization *other) const
{
  int iteration,order;
  double time1=getEndTime(iteration,order)-_time_tolerance;
  double time2=other->getStartTime(iteration,order)+other->getTimeTolerance();
  return time1<=time2;
}

bool MEDCouplingTimeDiscretization::isStrictlyBefore(const MEDCouplingTimeDiscretization *other) const
{
  int iteration,order;
  double time1=getEndTime(iteration,order)+_time_tolerance;
  double time2=other->getStartTime(iteration,order)-other->getTimeTolerance();
  return time1<time2;
}

MEDCouplingTimeDiscretization *MEDCouplingTimeDiscretization::doublyContractedProduct() const
{
  MEDCouplingTimeDiscretization *ret=MEDCouplingTimeDiscretization::New(getEnum());
  ret->setTimeUnit(getTimeUnit());
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > arrays2(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    {
      if(arrays[j])
        arrays2[j]=arrays[j]->doublyContractedProduct();
      else
        arrays2[j]=0;
    }
  std::vector<DataArrayDouble *> arrays3(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays3[j]=arrays2[j];
  ret->setArrays(arrays3,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingTimeDiscretization::determinant() const
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > arrays2(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    {
      if(arrays[j])
        arrays2[j]=arrays[j]->determinant();
      else
        arrays2[j]=0;
    }
  std::vector<DataArrayDouble *> arrays3(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays3[j]=arrays2[j];
  MEDCouplingTimeDiscretization *ret=MEDCouplingTimeDiscretization::New(getEnum());
  ret->setTimeUnit(getTimeUnit());
  ret->setArrays(arrays3,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingTimeDiscretization::eigenValues() const
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > arrays2(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    {
      if(arrays[j])
        arrays2[j]=arrays[j]->eigenValues();
      else
        arrays2[j]=0;
    }
  std::vector<DataArrayDouble *> arrays3(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays3[j]=arrays2[j];
  MEDCouplingTimeDiscretization *ret=MEDCouplingTimeDiscretization::New(getEnum());
  ret->setTimeUnit(getTimeUnit());
  ret->setArrays(arrays3,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingTimeDiscretization::eigenVectors() const
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > arrays2(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    {
      if(arrays[j])
        arrays2[j]=arrays[j]->eigenVectors();
      else
        arrays2[j]=0;
    }
  std::vector<DataArrayDouble *> arrays3(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays3[j]=arrays2[j];
  MEDCouplingTimeDiscretization *ret=MEDCouplingTimeDiscretization::New(getEnum());
  ret->setTimeUnit(getTimeUnit());
  ret->setArrays(arrays3,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingTimeDiscretization::inverse() const
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > arrays2(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    {
      if(arrays[j])
        arrays2[j]=arrays[j]->inverse();
      else
        arrays2[j]=0;
    }
  std::vector<DataArrayDouble *> arrays3(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays3[j]=arrays2[j];
  MEDCouplingTimeDiscretization *ret=MEDCouplingTimeDiscretization::New(getEnum());
  ret->setTimeUnit(getTimeUnit());
  ret->setArrays(arrays3,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingTimeDiscretization::trace() const
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > arrays2(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    {
      if(arrays[j])
        arrays2[j]=arrays[j]->trace();
      else
        arrays2[j]=0;
    }
  std::vector<DataArrayDouble *> arrays3(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays3[j]=arrays2[j];
  MEDCouplingTimeDiscretization *ret=MEDCouplingTimeDiscretization::New(getEnum());
  ret->setTimeUnit(getTimeUnit());
  ret->setArrays(arrays3,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingTimeDiscretization::deviator() const
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > arrays2(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    {
      if(arrays[j])
        arrays2[j]=arrays[j]->deviator();
      else
        arrays2[j]=0;
    }
  std::vector<DataArrayDouble *> arrays3(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays3[j]=arrays2[j];
  MEDCouplingTimeDiscretization *ret=MEDCouplingTimeDiscretization::New(getEnum());
  ret->setTimeUnit(getTimeUnit());
  ret->setArrays(arrays3,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingTimeDiscretization::magnitude() const
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > arrays2(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    {
      if(arrays[j])
        arrays2[j]=arrays[j]->magnitude();
      else
        arrays2[j]=0;
    }
  std::vector<DataArrayDouble *> arrays3(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays3[j]=arrays2[j];
  MEDCouplingTimeDiscretization *ret=MEDCouplingTimeDiscretization::New(getEnum());
  ret->setTimeUnit(getTimeUnit());
  ret->setArrays(arrays3,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingTimeDiscretization::negate() const
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > arrays2(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    {
      if(arrays[j])
        arrays2[j]=arrays[j]->negate();
      else
        arrays2[j]=0;
    }
  std::vector<DataArrayDouble *> arrays3(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays3[j]=arrays2[j];
  MEDCouplingTimeDiscretization *ret=MEDCouplingTimeDiscretization::New(getEnum());
  ret->setTimeUnit(getTimeUnit());
  ret->setArrays(arrays3,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingTimeDiscretization::maxPerTuple() const
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > arrays2(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    {
      if(arrays[j])
        arrays2[j]=arrays[j]->maxPerTuple();
      else
        arrays2[j]=0;
    }
  std::vector<DataArrayDouble *> arrays3(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays3[j]=arrays2[j];
  MEDCouplingTimeDiscretization *ret=MEDCouplingTimeDiscretization::New(getEnum());
  ret->setTimeUnit(getTimeUnit());
  ret->setArrays(arrays3,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingTimeDiscretization::keepSelectedComponents(const std::vector<int>& compoIds) const
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > arrays2(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    {
      if(arrays[j])
        arrays2[j]=static_cast<DataArrayDouble *>(arrays[j]->keepSelectedComponents(compoIds));
      else
        arrays2[j]=0;
    }
  std::vector<DataArrayDouble *> arrays3(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays3[j]=arrays2[j];
  MEDCouplingTimeDiscretization *ret=MEDCouplingTimeDiscretization::New(getEnum());
  ret->setTimeUnit(getTimeUnit());
  ret->setArrays(arrays3,0);
  return ret;
}

void MEDCouplingTimeDiscretization::setSelectedComponents(const MEDCouplingTimeDiscretization *other, const std::vector<int>& compoIds)
{
  std::vector<DataArrayDouble *> arrays1,arrays2;
  getArrays(arrays1);
  other->getArrays(arrays2);
  if(arrays1.size()!=arrays2.size())
    throw INTERP_KERNEL::Exception("TimeDiscretization::setSelectedComponents : number of arrays mismatch !");
  for(std::size_t i=0;i<arrays1.size();i++)
    {
      if(arrays1[i]!=0 && arrays2[i]!=0)
        arrays1[i]->setSelectedComponents(arrays2[i],compoIds);
      else if(arrays1[i]!=0 || arrays2[i]!=0)
        throw INTERP_KERNEL::Exception("TimeDiscretization::setSelectedComponents : some time array in correspondance are not defined symetrically !");
    }
}

void MEDCouplingTimeDiscretization::changeNbOfComponents(int newNbOfComp, double dftValue)
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > arrays2(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    {
      if(arrays[j])
        arrays2[j]=arrays[j]->changeNbOfComponents(newNbOfComp,dftValue);
      else
        arrays2[j]=0;
    }
  std::vector<DataArrayDouble *> arrays3(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays3[j]=arrays2[j];
  setArrays(arrays3,0);
}

void MEDCouplingTimeDiscretization::sortPerTuple(bool asc)
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  for(std::size_t j=0;j<arrays.size();j++)
    {
      if(arrays[j])
        arrays[j]->sortPerTuple(asc);
    }
}

void MEDCouplingTimeDiscretization::setUniformValue(int nbOfTuple, int nbOfCompo, double value)
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > arrays2(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    {
      if(arrays[j])
        {
          arrays2[j]=arrays[j]->changeNbOfComponents(nbOfCompo,value);
          arrays2[j]->fillWithValue(value);
        }
      else
        {
          arrays2[j]=DataArrayDouble::New();
          arrays2[j]->alloc(nbOfTuple,nbOfCompo);
          arrays2[j]->fillWithValue(value);
        }
    }
  std::vector<DataArrayDouble *> arrays3(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays3[j]=arrays2[j];
  setArrays(arrays3,0);
}

void MEDCouplingTimeDiscretization::setOrCreateUniformValueOnAllComponents(int nbOfTuple, double value)
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > arrays2(arrays.size());
  bool newArr=false;
  for(std::size_t j=0;j<arrays.size();j++)
    {
      if(arrays[j])
        {
          arrays2[j]=arrays[j]; arrays2[j]->incrRef();
          arrays2[j]->fillWithValue(value);
        }
      else
        {
          newArr=true;
          arrays2[j]=DataArrayDouble::New();
          arrays2[j]->alloc(nbOfTuple,1);
          arrays2[j]->fillWithValue(value);
        }
    }
  if(newArr)
    {
      std::vector<DataArrayDouble *> arrays3(arrays.size());
      for(std::size_t j=0;j<arrays.size();j++)
        arrays3[j]=arrays2[j];
      setArrays(arrays3,0);
    }
}

void MEDCouplingTimeDiscretization::applyLin(double a, double b, int compoId)
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  for(std::size_t j=0;j<arrays.size();j++)
    {
      if(arrays[j])
        arrays[j]->applyLin(a,b,compoId);
    }
}

void MEDCouplingTimeDiscretization::applyLin(double a, double b)
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  for(std::size_t j=0;j<arrays.size();j++)
    {
      if(arrays[j])
        arrays[j]->applyLin(a,b);
    }
}

void MEDCouplingTimeDiscretization::applyFunc(int nbOfComp, FunctionToEvaluate func)
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > arrays2(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    {
      if(arrays[j])
        arrays2[j]=arrays[j]->applyFunc(nbOfComp,func);
      else
        arrays2[j]=0;
    }
  std::vector<DataArrayDouble *> arrays3(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays3[j]=arrays2[j];
  setArrays(arrays3,0);
}

void MEDCouplingTimeDiscretization::applyFunc(int nbOfComp, const std::string& func)
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > arrays2(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    {
      if(arrays[j])
        arrays2[j]=arrays[j]->applyFunc(nbOfComp,func);
      else
        arrays2[j]=0;
    }
  std::vector<DataArrayDouble *> arrays3(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays3[j]=arrays2[j];
  setArrays(arrays3,0);
}

void MEDCouplingTimeDiscretization::applyFunc2(int nbOfComp, const std::string& func)
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > arrays2(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    {
      if(arrays[j])
        arrays2[j]=arrays[j]->applyFunc2(nbOfComp,func);
      else
        arrays2[j]=0;
    }
  std::vector<DataArrayDouble *> arrays3(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays3[j]=arrays2[j];
  setArrays(arrays3,0);
}

void MEDCouplingTimeDiscretization::applyFunc3(int nbOfComp, const std::vector<std::string>& varsOrder, const std::string& func)
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > arrays2(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    {
      if(arrays[j])
        arrays2[j]=arrays[j]->applyFunc3(nbOfComp,varsOrder,func);
      else
        arrays2[j]=0;
    }
  std::vector<DataArrayDouble *> arrays3(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays3[j]=arrays2[j];
  setArrays(arrays3,0);
}

void MEDCouplingTimeDiscretization::applyFunc(const std::string& func)
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > arrays2(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    {
      if(arrays[j])
        arrays2[j]=arrays[j]->applyFunc(func);
      else
        arrays2[j]=0;
    }
  std::vector<DataArrayDouble *> arrays3(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays3[j]=arrays2[j];
  setArrays(arrays3,0);
}

void MEDCouplingTimeDiscretization::applyFuncFast32(const std::string& func)
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  for(std::size_t j=0;j<arrays.size();j++)
    {
      if(arrays[j])
        arrays[j]->applyFuncFast32(func);
    }
}

void MEDCouplingTimeDiscretization::applyFuncFast64(const std::string& func)
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  for(std::size_t j=0;j<arrays.size();j++)
    {
      if(arrays[j])
        arrays[j]->applyFuncFast64(func);
    }
}

void MEDCouplingTimeDiscretization::fillFromAnalytic(const DataArrayDouble *loc, int nbOfComp, FunctionToEvaluate func)
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > arrays2(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays2[j]=loc->applyFunc(nbOfComp,func);
  std::vector<DataArrayDouble *> arrays3(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays3[j]=arrays2[j];
  setArrays(arrays3,0);
}

void MEDCouplingTimeDiscretization::fillFromAnalytic(const DataArrayDouble *loc, int nbOfComp, const std::string& func)
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > arrays2(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays2[j]=loc->applyFunc(nbOfComp,func);
  std::vector<DataArrayDouble *> arrays3(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays3[j]=arrays2[j];
  setArrays(arrays3,0);
}

void MEDCouplingTimeDiscretization::fillFromAnalytic2(const DataArrayDouble *loc, int nbOfComp, const std::string& func)
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > arrays2(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays2[j]=loc->applyFunc2(nbOfComp,func);
  std::vector<DataArrayDouble *> arrays3(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays3[j]=arrays2[j];
  setArrays(arrays3,0);
}

void MEDCouplingTimeDiscretization::fillFromAnalytic3(const DataArrayDouble *loc, int nbOfComp, const std::vector<std::string>& varsOrder, const std::string& func)
{
  std::vector<DataArrayDouble *> arrays;
  getArrays(arrays);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > arrays2(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays2[j]=loc->applyFunc3(nbOfComp,varsOrder,func);
  std::vector<DataArrayDouble *> arrays3(arrays.size());
  for(std::size_t j=0;j<arrays.size();j++)
    arrays3[j]=arrays2[j];
  setArrays(arrays3,0);
}

MEDCouplingNoTimeLabel::MEDCouplingNoTimeLabel()
{
}

MEDCouplingNoTimeLabel::MEDCouplingNoTimeLabel(const MEDCouplingTimeDiscretization& other, bool deepCpy):MEDCouplingTimeDiscretization(other,deepCpy)
{
}

std::string MEDCouplingNoTimeLabel::getStringRepr() const
{
  std::ostringstream stream;
  stream << REPR;
  stream << "\nTime unit is : \"" << _time_unit << "\"";
  return stream.str();
}

void MEDCouplingNoTimeLabel::synchronizeTimeWith(const MEDCouplingMesh *mesh)
{
  throw INTERP_KERNEL::Exception("MEDCouplingNoTimeLabel::synchronizeTimeWith : impossible to synchronize time with a MEDCouplingMesh because the time discretization is incompatible with it !");
}

bool MEDCouplingNoTimeLabel::areCompatible(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areCompatible(other))
    return false;
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  return otherC!=0;
}

bool MEDCouplingNoTimeLabel::areStrictlyCompatible(const MEDCouplingTimeDiscretization *other, std::string& reason) const
{
  if(!MEDCouplingTimeDiscretization::areStrictlyCompatible(other,reason))
    return false;
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  bool ret=otherC!=0;
  if(!ret)
    reason.insert(0,"time discretization of this is NO_TIME, other has a different time discretization.");
  return ret;
}

bool MEDCouplingNoTimeLabel::areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areStrictlyCompatibleForMul(other))
    return false;
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  return otherC!=0;
}

bool MEDCouplingNoTimeLabel::areStrictlyCompatibleForDiv(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areStrictlyCompatibleForDiv(other))
    return false;
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  return otherC!=0;
}

bool MEDCouplingNoTimeLabel::areCompatibleForMeld(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areCompatibleForMeld(other))
    return false;
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  return otherC!=0;
}

bool MEDCouplingNoTimeLabel::isEqualIfNotWhy(const MEDCouplingTimeDiscretization *other, double prec, std::string& reason) const
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    {
      reason="This has time discretization NO_TIME, other not.";
      return false;
    }
  return MEDCouplingTimeDiscretization::isEqualIfNotWhy(other,prec,reason);
}

bool MEDCouplingNoTimeLabel::isEqualWithoutConsideringStr(const MEDCouplingTimeDiscretization *other, double prec) const
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    return false;
  return MEDCouplingTimeDiscretization::isEqualWithoutConsideringStr(other,prec);
}

MEDCouplingTimeDiscretization *MEDCouplingNoTimeLabel::aggregate(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("NoTimeLabel::aggregation on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Aggregate(getArray(),other->getArray());
  MEDCouplingNoTimeLabel *ret=new MEDCouplingNoTimeLabel;
  ret->setArray(arr,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingNoTimeLabel::aggregate(const std::vector<const MEDCouplingTimeDiscretization *>& other) const
{
  std::vector<const DataArrayDouble *> a(other.size());
  int i=0;
  for(std::vector<const MEDCouplingTimeDiscretization *>::const_iterator it=other.begin();it!=other.end();it++,i++)
    {
      const MEDCouplingNoTimeLabel *itC=dynamic_cast<const MEDCouplingNoTimeLabel *>(*it);
      if(!itC)
        throw INTERP_KERNEL::Exception("NoTimeLabel::aggregate on mismatched time discretization !");
      a[i]=itC->getArray();
    }
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Aggregate(a);
  MEDCouplingNoTimeLabel *ret=new MEDCouplingNoTimeLabel;
  ret->setArray(arr,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingNoTimeLabel::meld(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("NoTimeLabel::meld on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Meld(getArray(),other->getArray());
  MEDCouplingNoTimeLabel *ret=new MEDCouplingNoTimeLabel;
  ret->setTimeTolerance(getTimeTolerance());
  ret->setArray(arr,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingNoTimeLabel::dot(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("NoTimeLabel::dot on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Dot(getArray(),other->getArray());
  MEDCouplingNoTimeLabel *ret=new MEDCouplingNoTimeLabel;
  ret->setArray(arr,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingNoTimeLabel::crossProduct(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("NoTimeLabel::crossProduct on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::CrossProduct(getArray(),other->getArray());
  MEDCouplingNoTimeLabel *ret=new MEDCouplingNoTimeLabel;
  ret->setArray(arr,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingNoTimeLabel::max(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("NoTimeLabel::max on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Max(getArray(),other->getArray());
  MEDCouplingNoTimeLabel *ret=new MEDCouplingNoTimeLabel;
  ret->setArray(arr,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingNoTimeLabel::min(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("NoTimeLabel::max on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Min(getArray(),other->getArray());
  MEDCouplingNoTimeLabel *ret=new MEDCouplingNoTimeLabel;
  ret->setArray(arr,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingNoTimeLabel::add(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("NoTimeLabel::add on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Add(getArray(),other->getArray());
  MEDCouplingNoTimeLabel *ret=new MEDCouplingNoTimeLabel;
  ret->setArray(arr,0);
  return ret;
}

void MEDCouplingNoTimeLabel::addEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("NoTimeLabel::addEqual on mismatched time discretization !");
  if(!getArray())
    throw INTERP_KERNEL::Exception("MEDCouplingNoTimeLabel::addEqual : Data Array is NULL !");
  getArray()->addEqual(other->getArray());
}

MEDCouplingTimeDiscretization *MEDCouplingNoTimeLabel::substract(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("NoTimeLabel::substract on mismatched time discretization !");
  if(!getArray())
    throw INTERP_KERNEL::Exception("MEDCouplingNoTimeLabel::substract : Data Array is NULL !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Substract(getArray(),other->getArray());
  MEDCouplingNoTimeLabel *ret=new MEDCouplingNoTimeLabel;
  ret->setArray(arr,0);
  return ret;
}

void MEDCouplingNoTimeLabel::substractEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("NoTimeLabel::substractEqual on mismatched time discretization !");
  if(!getArray())
    throw INTERP_KERNEL::Exception("MEDCouplingNoTimeLabel::substractEqual : Data Array is NULL !");
  getArray()->substractEqual(other->getArray());
}

MEDCouplingTimeDiscretization *MEDCouplingNoTimeLabel::multiply(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("NoTimeLabel::multiply on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Multiply(getArray(),other->getArray());
  MEDCouplingNoTimeLabel *ret=new MEDCouplingNoTimeLabel;
  ret->setArray(arr,0);
  return ret;
}

void MEDCouplingNoTimeLabel::multiplyEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("NoTimeLabel::multiplyEqual on mismatched time discretization !");
  if(!getArray())
    throw INTERP_KERNEL::Exception("MEDCouplingNoTimeLabel::multiplyEqual : Data Array is NULL !");
  getArray()->multiplyEqual(other->getArray());
}

MEDCouplingTimeDiscretization *MEDCouplingNoTimeLabel::divide(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("divide on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Divide(getArray(),other->getArray());
  MEDCouplingNoTimeLabel *ret=new MEDCouplingNoTimeLabel;
  ret->setArray(arr,0);
  return ret;
}

void MEDCouplingNoTimeLabel::divideEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("NoTimeLabel::divideEqual on mismatched time discretization !");
  if(!getArray())
    throw INTERP_KERNEL::Exception("MEDCouplingNoTimeLabel::divideEqual : Data Array is NULL !");
  getArray()->divideEqual(other->getArray());
}

MEDCouplingTimeDiscretization *MEDCouplingNoTimeLabel::pow(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("pow on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Pow(getArray(),other->getArray());
  MEDCouplingNoTimeLabel *ret=new MEDCouplingNoTimeLabel;
  ret->setArray(arr,0);
  return ret;
}

void MEDCouplingNoTimeLabel::powEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("NoTimeLabel::powEqual on mismatched time discretization !");
  if(!getArray())
    throw INTERP_KERNEL::Exception("MEDCouplingNoTimeLabel::powEqual : Data Array is NULL !");
  getArray()->powEqual(other->getArray());
}

MEDCouplingTimeDiscretization *MEDCouplingNoTimeLabel::performCpy(bool deepCpy) const
{
  return new MEDCouplingNoTimeLabel(*this,deepCpy);
}

void MEDCouplingNoTimeLabel::checkTimePresence(double time) const
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

std::vector< const DataArrayDouble *> MEDCouplingNoTimeLabel::getArraysForTime(double time) const
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingNoTimeLabel::getValueForTime(double time, const std::vector<double>& vals, double *res) const
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

bool MEDCouplingNoTimeLabel::isBefore(const MEDCouplingTimeDiscretization *other) const
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

bool MEDCouplingNoTimeLabel::isStrictlyBefore(const MEDCouplingTimeDiscretization *other) const
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

double MEDCouplingNoTimeLabel::getStartTime(int& iteration, int& order) const
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

double MEDCouplingNoTimeLabel::getEndTime(int& iteration, int& order) const
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingNoTimeLabel::setStartIteration(int it)
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingNoTimeLabel::setEndIteration(int it)
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingNoTimeLabel::setStartOrder(int order)
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingNoTimeLabel::setEndOrder(int order)
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingNoTimeLabel::setStartTimeValue(double time)
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingNoTimeLabel::setEndTimeValue(double time)
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingNoTimeLabel::setStartTime(double time, int iteration, int order)
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingNoTimeLabel::setEndTime(double time, int iteration, int order)
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingNoTimeLabel::getValueOnTime(int eltId, double time, double *value) const
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingNoTimeLabel::getValueOnDiscTime(int eltId, int iteration, int order, double *value) const
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

/*!
 * idem getTinySerializationIntInformation except that it is for multi field fetch
 */
void MEDCouplingNoTimeLabel::getTinySerializationIntInformation2(std::vector<int>& tinyInfo) const
{
  tinyInfo.clear();
}

/*!
 * idem getTinySerializationDbleInformation except that it is for multi field fetch
 */
void MEDCouplingNoTimeLabel::getTinySerializationDbleInformation2(std::vector<double>& tinyInfo) const
{
  tinyInfo.resize(1);
  tinyInfo[0]=_time_tolerance;
}

/*!
 * idem finishUnserialization except that it is for multi field fetch
 */
void MEDCouplingNoTimeLabel::finishUnserialization2(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD)
{
  _time_tolerance=tinyInfoD[0];
}

MEDCouplingWithTimeStep::MEDCouplingWithTimeStep(const MEDCouplingWithTimeStep& other, bool deepCpy):MEDCouplingTimeDiscretization(other,deepCpy),
    _time(other._time),_iteration(other._iteration),_order(other._order)
{
}

MEDCouplingWithTimeStep::MEDCouplingWithTimeStep():_time(0.),_iteration(-1),_order(-1)
{
}

std::string MEDCouplingWithTimeStep::getStringRepr() const
{
  std::ostringstream stream;
  stream << REPR << " Time is defined by iteration=" << _iteration << " order=" << _order << " and time=" << _time << ".";
  stream << "\nTime unit is : \"" << _time_unit << "\"";
  return stream.str();
}

void MEDCouplingWithTimeStep::synchronizeTimeWith(const MEDCouplingMesh *mesh)
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingWithTimeStep::synchronizeTimeWith : mesh instance is NULL ! Impossible to synchronize time !");
  int it=-1,order=-1;
  double val=mesh->getTime(it,order);
  _time=val; _iteration=it; _order=order;
  std::string tUnit=mesh->getTimeUnit();
  _time_unit=tUnit;
}

void MEDCouplingWithTimeStep::getTinySerializationIntInformation(std::vector<int>& tinyInfo) const
{
  MEDCouplingTimeDiscretization::getTinySerializationIntInformation(tinyInfo);
  tinyInfo.push_back(_iteration);
  tinyInfo.push_back(_order);
}

void MEDCouplingWithTimeStep::getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const
{
  MEDCouplingTimeDiscretization::getTinySerializationDbleInformation(tinyInfo);
  tinyInfo.push_back(_time);
}

void MEDCouplingWithTimeStep::finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS)
{
  MEDCouplingTimeDiscretization::finishUnserialization(tinyInfoI,tinyInfoD,tinyInfoS);
  _time=tinyInfoD[1];
  _iteration=tinyInfoI[2];
  _order=tinyInfoI[3];
}

/*!
 * idem getTinySerializationIntInformation except that it is for multi field fetch
 */
void MEDCouplingWithTimeStep::getTinySerializationIntInformation2(std::vector<int>& tinyInfo) const
{
  tinyInfo.resize(2);
  tinyInfo[0]=_iteration;
  tinyInfo[1]=_order;
}

/*!
 * idem getTinySerializationDbleInformation except that it is for multi field fetch
 */
void MEDCouplingWithTimeStep::getTinySerializationDbleInformation2(std::vector<double>& tinyInfo) const
{
  tinyInfo.resize(2);
  tinyInfo[0]=_time_tolerance;
  tinyInfo[1]=_time;
}

/*!
 * idem finishUnserialization except that it is for multi field fetch
 */
void MEDCouplingWithTimeStep::finishUnserialization2(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD)
{
  _iteration=tinyInfoI[0];
  _order=tinyInfoI[1];
  _time_tolerance=tinyInfoD[0];
  _time=tinyInfoD[1];
}

bool MEDCouplingWithTimeStep::areCompatible(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areCompatible(other))
    return false;
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  return otherC!=0;
}

bool MEDCouplingWithTimeStep::areStrictlyCompatible(const MEDCouplingTimeDiscretization *other, std::string& reason) const
{
  if(!MEDCouplingTimeDiscretization::areStrictlyCompatible(other,reason))
    return false;
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  bool ret=otherC!=0;
  if(!ret)
    reason.insert(0,"time discretization of this is ONE_TIME, other has a different time discretization.");
  return ret;
}

bool MEDCouplingWithTimeStep::areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areStrictlyCompatibleForMul(other))
    return false;
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  return otherC!=0;
}

bool MEDCouplingWithTimeStep::areStrictlyCompatibleForDiv(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areStrictlyCompatibleForDiv(other))
    return false;
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  return otherC!=0;
}

bool MEDCouplingWithTimeStep::areCompatibleForMeld(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areCompatibleForMeld(other))
    return false;
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  return otherC!=0;
}

bool MEDCouplingWithTimeStep::isEqualIfNotWhy(const MEDCouplingTimeDiscretization *other, double prec, std::string& reason) const
{
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  std::ostringstream oss; oss.precision(15);
  if(!otherC)
    {
      reason="This has time discretization ONE_TIME, other not.";
      return false;
    }
  if(_iteration!=otherC->_iteration)
    {
      oss << "iterations differ. this iteration=" << _iteration << " other iteration=" << otherC->_iteration;
      reason=oss.str();
      return false;
    }
  if(_order!=otherC->_order)
    {
      oss << "orders differ. this order=" << _order << " other order=" << otherC->_order;
      reason=oss.str();
      return false;
    }
  if(std::fabs(_time-otherC->_time)>_time_tolerance)
    {
      oss << "times differ. this time=" << _time << " other time=" << otherC->_time;
      reason=oss.str();
      return false;
    }
  return MEDCouplingTimeDiscretization::isEqualIfNotWhy(other,prec,reason);
}

bool MEDCouplingWithTimeStep::isEqualWithoutConsideringStr(const MEDCouplingTimeDiscretization *other, double prec) const
{
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  if(!otherC)
    return false;
  if(_iteration!=otherC->_iteration)
    return false;
  if(_order!=otherC->_order)
    return false;
  if(std::fabs(_time-otherC->_time)>_time_tolerance)
    return false;
  return MEDCouplingTimeDiscretization::isEqualWithoutConsideringStr(other,prec);
}

void MEDCouplingWithTimeStep::copyTinyAttrFrom(const MEDCouplingTimeDiscretization& other)
{
  MEDCouplingTimeDiscretization::copyTinyAttrFrom(other);
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(&other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("MEDCouplingWithTimeStep::copyTinyAttrFrom : mismatch of time discretization !");
  _time=otherC->_time;
  _iteration=otherC->_iteration;
  _order=otherC->_order;
}

MEDCouplingTimeDiscretization *MEDCouplingWithTimeStep::aggregate(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("WithTimeStep::aggregation on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Aggregate(getArray(),other->getArray());
  MEDCouplingWithTimeStep *ret=new MEDCouplingWithTimeStep;
  ret->setArray(arr,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingWithTimeStep::aggregate(const std::vector<const MEDCouplingTimeDiscretization *>& other) const
{
  std::vector<const DataArrayDouble *> a(other.size());
  int i=0;
  for(std::vector<const MEDCouplingTimeDiscretization *>::const_iterator it=other.begin();it!=other.end();it++,i++)
    {
      const MEDCouplingWithTimeStep *itC=dynamic_cast<const MEDCouplingWithTimeStep *>(*it);
      if(!itC)
        throw INTERP_KERNEL::Exception("WithTimeStep::aggregate on mismatched time discretization !");
      a[i]=itC->getArray();
    }
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Aggregate(a);
  MEDCouplingWithTimeStep *ret=new MEDCouplingWithTimeStep;
  ret->setArray(arr,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingWithTimeStep::meld(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("WithTimeStep::meld on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Meld(getArray(),other->getArray());
  MEDCouplingWithTimeStep *ret=new MEDCouplingWithTimeStep;
  ret->setArray(arr,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingWithTimeStep::dot(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("WithTimeStep::dot on mismatched time discretization !");
  MEDCouplingWithTimeStep *ret=new MEDCouplingWithTimeStep;
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Dot(getArray(),other->getArray());
  ret->setArray(arr,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingWithTimeStep::crossProduct(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("WithTimeStep::crossProduct on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::CrossProduct(getArray(),other->getArray());
  MEDCouplingWithTimeStep *ret=new MEDCouplingWithTimeStep;
  ret->setArray(arr,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingWithTimeStep::max(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("WithTimeStep::max on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Max(getArray(),other->getArray());
  MEDCouplingWithTimeStep *ret=new MEDCouplingWithTimeStep;
  ret->setArray(arr,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingWithTimeStep::min(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("WithTimeStep::min on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Min(getArray(),other->getArray());
  MEDCouplingWithTimeStep *ret=new MEDCouplingWithTimeStep;
  ret->setArray(arr,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingWithTimeStep::add(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("WithTimeStep::add on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Add(getArray(),other->getArray());
  MEDCouplingWithTimeStep *ret=new MEDCouplingWithTimeStep;
  ret->setArray(arr,0);
  int tmp1,tmp2;
  double tmp3=getStartTime(tmp1,tmp2);
  ret->setStartTime(tmp3,tmp1,tmp2);
  return ret;
}

void MEDCouplingWithTimeStep::addEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("WithTimeStep::addEqual on mismatched time discretization !");
  if(!getArray())
    throw INTERP_KERNEL::Exception("MEDCouplingWithTimeLabel::addEqual : Data Array is NULL !");
  getArray()->addEqual(other->getArray());
}

MEDCouplingTimeDiscretization *MEDCouplingWithTimeStep::substract(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("WithTimeStep::substract on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Substract(getArray(),other->getArray());
  MEDCouplingWithTimeStep *ret=new MEDCouplingWithTimeStep;
  ret->setArray(arr,0);
  int tmp1,tmp2;
  double tmp3=getStartTime(tmp1,tmp2);
  ret->setStartTime(tmp3,tmp1,tmp2);
  return ret;
}

void MEDCouplingWithTimeStep::substractEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("WithTimeStep::substractEqual on mismatched time discretization !");
  if(!getArray())
    throw INTERP_KERNEL::Exception("MEDCouplingWithTimeLabel::substractEqual : Data Array is NULL !");
  getArray()->substractEqual(other->getArray());
}

MEDCouplingTimeDiscretization *MEDCouplingWithTimeStep::multiply(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("WithTimeStep::multiply on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Multiply(getArray(),other->getArray());
  MEDCouplingWithTimeStep *ret=new MEDCouplingWithTimeStep;
  ret->setArray(arr,0);
  int tmp1,tmp2;
  double tmp3=getStartTime(tmp1,tmp2);
  ret->setStartTime(tmp3,tmp1,tmp2);
  return ret;
}

void MEDCouplingWithTimeStep::multiplyEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("WithTimeStep::multiplyEqual on mismatched time discretization !");
  if(!getArray())
    throw INTERP_KERNEL::Exception("MEDCouplingWithTimeLabel::multiplyEqual : Data Array is NULL !");
  getArray()->multiplyEqual(other->getArray());
}

MEDCouplingTimeDiscretization *MEDCouplingWithTimeStep::divide(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("WithTimeStep::divide on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Divide(getArray(),other->getArray());
  MEDCouplingWithTimeStep *ret=new MEDCouplingWithTimeStep;
  ret->setArray(arr,0);
  int tmp1,tmp2;
  double tmp3=getStartTime(tmp1,tmp2);
  ret->setStartTime(tmp3,tmp1,tmp2);
  return ret;
}

void MEDCouplingWithTimeStep::divideEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("WithTimeStep::divideEqual on mismatched time discretization !");
  if(!getArray())
    throw INTERP_KERNEL::Exception("MEDCouplingWithTimeLabel::divideEqual : Data Array is NULL !");
  getArray()->divideEqual(other->getArray());
}

MEDCouplingTimeDiscretization *MEDCouplingWithTimeStep::pow(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("WithTimeStep::pow on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Pow(getArray(),other->getArray());
  MEDCouplingWithTimeStep *ret=new MEDCouplingWithTimeStep;
  ret->setArray(arr,0);
  int tmp1,tmp2;
  double tmp3=getStartTime(tmp1,tmp2);
  ret->setStartTime(tmp3,tmp1,tmp2);
  return ret;
}

void MEDCouplingWithTimeStep::powEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("WithTimeStep::powEqual on mismatched time discretization !");
  if(!getArray())
    throw INTERP_KERNEL::Exception("MEDCouplingWithTimeLabel::powEqual : Data Array is NULL !");
  getArray()->powEqual(other->getArray());
}

MEDCouplingTimeDiscretization *MEDCouplingWithTimeStep::performCpy(bool deepCpy) const
{
  return new MEDCouplingWithTimeStep(*this,deepCpy);
}

void MEDCouplingWithTimeStep::checkNoTimePresence() const
{
  throw INTERP_KERNEL::Exception("No time specified on a field defined on one time");
}

void MEDCouplingWithTimeStep::checkTimePresence(double time) const
{
  if(std::fabs(time-_time)>_time_tolerance)
    {
      std::ostringstream stream;
      stream << "The field is defined on time " << _time << " with eps=" << _time_tolerance << " and asking time = " << time << " !";
      throw INTERP_KERNEL::Exception(stream.str().c_str());
    }
}

std::vector< const DataArrayDouble *> MEDCouplingWithTimeStep::getArraysForTime(double time) const
{
  if(std::fabs(time-_time)<=_time_tolerance)
    {
      std::vector< const DataArrayDouble *> ret(1);
      ret[0]=_array;
      return ret;
    }
  else
    throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingWithTimeStep::getValueForTime(double time, const std::vector<double>& vals, double *res) const
{
  std::copy(vals.begin(),vals.end(),res);
}

void MEDCouplingWithTimeStep::getValueOnTime(int eltId, double time, double *value) const
{
  if(std::fabs(time-_time)<=_time_tolerance)
    if(_array)
      _array->getTuple(eltId,value);
    else
      throw INTERP_KERNEL::Exception("No array existing.");
  else
    throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingWithTimeStep::getValueOnDiscTime(int eltId, int iteration, int order, double *value) const
{
  if(_iteration==iteration && _order==order)
    if(_array)
      _array->getTuple(eltId,value);
    else
      throw INTERP_KERNEL::Exception("No array existing.");
  else
    throw INTERP_KERNEL::Exception("No data on this discrete time.");
}

MEDCouplingConstOnTimeInterval::MEDCouplingConstOnTimeInterval():_start_time(0.),_end_time(0.),_start_iteration(-1),_end_iteration(-1),_start_order(-1),_end_order(-1)
{
}

void MEDCouplingConstOnTimeInterval::copyTinyAttrFrom(const MEDCouplingTimeDiscretization& other)
{
  MEDCouplingTimeDiscretization::copyTinyAttrFrom(other);
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(&other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("MEDCouplingConstOnTimeInterval::copyTinyAttrFrom : mismatch of time discretization !");
  _start_time=otherC->_start_time;
  _end_time=otherC->_end_time;
  _start_iteration=otherC->_start_iteration;
  _end_iteration=otherC->_end_iteration;
  _start_order=otherC->_start_order;
  _end_order=otherC->_end_order;
}

void MEDCouplingConstOnTimeInterval::getTinySerializationIntInformation(std::vector<int>& tinyInfo) const
{
  MEDCouplingTimeDiscretization::getTinySerializationIntInformation(tinyInfo);
  tinyInfo.push_back(_start_iteration);
  tinyInfo.push_back(_start_order);
  tinyInfo.push_back(_end_iteration);
  tinyInfo.push_back(_end_order);
}

void MEDCouplingConstOnTimeInterval::getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const
{
  MEDCouplingTimeDiscretization::getTinySerializationDbleInformation(tinyInfo);
  tinyInfo.push_back(_start_time);
  tinyInfo.push_back(_end_time);
}

void MEDCouplingConstOnTimeInterval::finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS)
{
  MEDCouplingTimeDiscretization::finishUnserialization(tinyInfoI,tinyInfoD,tinyInfoS);
  _start_time=tinyInfoD[1];
  _end_time=tinyInfoD[2];
  _start_iteration=tinyInfoI[2];
  _start_order=tinyInfoI[3];
  _end_iteration=tinyInfoI[4];
  _end_order=tinyInfoI[5];
}

/*!
 * idem getTinySerializationIntInformation except that it is for multi field fetch
 */
void MEDCouplingConstOnTimeInterval::getTinySerializationIntInformation2(std::vector<int>& tinyInfo) const
{
  tinyInfo.resize(4);
  tinyInfo[0]=_start_iteration;
  tinyInfo[1]=_start_order;
  tinyInfo[2]=_end_iteration;
  tinyInfo[3]=_end_order;
}

/*!
 * idem getTinySerializationDbleInformation except that it is for multi field fetch
 */
void MEDCouplingConstOnTimeInterval::getTinySerializationDbleInformation2(std::vector<double>& tinyInfo) const
{
  tinyInfo.resize(3);
  tinyInfo[0]=_time_tolerance;
  tinyInfo[1]=_start_time;
  tinyInfo[2]=_end_time;
}

/*!
 * idem finishUnserialization except that it is for multi field fetch
 */
void MEDCouplingConstOnTimeInterval::finishUnserialization2(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD)
{
  _start_iteration=tinyInfoI[0];
  _start_order=tinyInfoI[1];
  _end_iteration=tinyInfoI[2];
  _end_order=tinyInfoI[3];
  _time_tolerance=tinyInfoD[0];
  _start_time=tinyInfoD[1];
  _end_time=tinyInfoD[2];
}

MEDCouplingConstOnTimeInterval::MEDCouplingConstOnTimeInterval(const MEDCouplingConstOnTimeInterval& other, bool deepCpy):
      MEDCouplingTimeDiscretization(other,deepCpy),_start_time(other._start_time),_end_time(other._end_time),_start_iteration(other._start_iteration),
      _end_iteration(other._end_iteration),_start_order(other._start_order),_end_order(other._end_order)
{
}

std::string MEDCouplingConstOnTimeInterval::getStringRepr() const
{
  std::ostringstream stream;
  stream << REPR << " Time interval is defined by :\niteration_start=" << _start_iteration << " order_start=" << _start_order << " and time_start=" << _start_time << "\n";
  stream << "iteration_end=" << _end_iteration << " order_end=" << _end_order << " and end_time=" << _end_time << "\n";
  stream << "\nTime unit is : \"" << _time_unit << "\"";
  return stream.str();
}

void MEDCouplingConstOnTimeInterval::synchronizeTimeWith(const MEDCouplingMesh *mesh)
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingWithTimeStep::synchronizeTimeWith : mesh instance is NULL ! Impossible to synchronize time !");
  int it=-1,order=-1;
  double val=mesh->getTime(it,order);
  _start_time=val; _start_iteration=it; _start_order=order;
  _end_time=val; _end_iteration=it; _end_order=order;
  std::string tUnit=mesh->getTimeUnit();
  _time_unit=tUnit;
}

MEDCouplingTimeDiscretization *MEDCouplingConstOnTimeInterval::performCpy(bool deepCpy) const
{
  return new MEDCouplingConstOnTimeInterval(*this,deepCpy);
}

std::vector< const DataArrayDouble *> MEDCouplingConstOnTimeInterval::getArraysForTime(double time) const
{
  if(time>_start_time-_time_tolerance && time<_end_time+_time_tolerance)
    {
      std::vector< const DataArrayDouble *> ret(1);
      ret[0]=_array;
      return ret;
    }
  else
    throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingConstOnTimeInterval::getValueForTime(double time, const std::vector<double>& vals, double *res) const
{
  std::copy(vals.begin(),vals.end(),res);
}

bool MEDCouplingConstOnTimeInterval::areCompatible(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areCompatible(other))
    return false;
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  return otherC!=0;
}

bool MEDCouplingConstOnTimeInterval::areStrictlyCompatible(const MEDCouplingTimeDiscretization *other, std::string& reason) const
{
  if(!MEDCouplingTimeDiscretization::areStrictlyCompatible(other,reason))
    return false;
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  bool ret=otherC!=0;
  if(!ret)
    reason.insert(0,"time discretization of this is CONST_ON_TIME_INTERVAL, other has a different time discretization.");
  return ret;
}

bool MEDCouplingConstOnTimeInterval::areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areStrictlyCompatibleForMul(other))
    return false;
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  return otherC!=0;
}

bool MEDCouplingConstOnTimeInterval::areStrictlyCompatibleForDiv(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areStrictlyCompatibleForDiv(other))
    return false;
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  return otherC!=0;
}

bool MEDCouplingConstOnTimeInterval::areCompatibleForMeld(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areCompatibleForMeld(other))
    return false;
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  return otherC!=0;
}

bool MEDCouplingConstOnTimeInterval::isEqualIfNotWhy(const MEDCouplingTimeDiscretization *other, double prec, std::string& reason) const
{
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  std::ostringstream oss; oss.precision(15);
  if(!otherC)
    {
      reason="This has time discretization CONST_ON_TIME_INTERVAL, other not.";
      return false;
    }
  if(_start_iteration!=otherC->_start_iteration)
    {
      oss << "start iterations differ. this start iteration=" << _start_iteration << " other start iteration=" << otherC->_start_iteration;
      reason=oss.str();
      return false;
    }
  if(_start_order!=otherC->_start_order)
    {
      oss << "start orders differ. this start order=" << _start_order << " other start order=" << otherC->_start_order;
      reason=oss.str();
      return false;
    }
  if(std::fabs(_start_time-otherC->_start_time)>_time_tolerance)
    {
      oss << "start times differ. this start time=" << _start_time << " other start time=" << otherC->_start_time;
      reason=oss.str();
      return false;
    }
  if(_end_iteration!=otherC->_end_iteration)
    {
      oss << "end iterations differ. this end iteration=" << _end_iteration << " other end iteration=" << otherC->_end_iteration;
      reason=oss.str();
      return false;
    }
  if(_end_order!=otherC->_end_order)
    {
      oss << "end orders differ. this end order=" << _end_order << " other end order=" << otherC->_end_order;
      reason=oss.str();
      return false;
    }
  if(std::fabs(_end_time-otherC->_end_time)>_time_tolerance)
    {
      oss << "end times differ. this end time=" << _end_time << " other end time=" << otherC->_end_time;
      reason=oss.str();
      return false;
    }
  return MEDCouplingTimeDiscretization::isEqualIfNotWhy(other,prec,reason);
}

bool MEDCouplingConstOnTimeInterval::isEqualWithoutConsideringStr(const MEDCouplingTimeDiscretization *other, double prec) const
{
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  if(!otherC)
    return false;
  if(_start_iteration!=otherC->_start_iteration)
    return false;
  if(_start_order!=otherC->_start_order)
    return false;
  if(std::fabs(_start_time-otherC->_start_time)>_time_tolerance)
    return false;
  if(_end_iteration!=otherC->_end_iteration)
    return false;
  if(_end_order!=otherC->_end_order)
    return false;
  if(std::fabs(_end_time-otherC->_end_time)>_time_tolerance)
    return false;
  return MEDCouplingTimeDiscretization::isEqualWithoutConsideringStr(other,prec);
}

void MEDCouplingConstOnTimeInterval::getValueOnTime(int eltId, double time, double *value) const
{
  if(time>_start_time-_time_tolerance && time<_end_time+_time_tolerance)
    if(_array)
      _array->getTuple(eltId,value);
    else
      throw INTERP_KERNEL::Exception("No array existing.");
  else
    throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingConstOnTimeInterval::getValueOnDiscTime(int eltId, int iteration, int order, double *value) const
{
  if(iteration>=_start_iteration && iteration<=_end_iteration)
    if(_array)
      _array->getTuple(eltId,value);
    else
      throw INTERP_KERNEL::Exception("No array existing.");
  else
    throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingConstOnTimeInterval::checkNoTimePresence() const
{
  throw INTERP_KERNEL::Exception("No time specified on a field defined as constant on one time interval");
}

void MEDCouplingConstOnTimeInterval::checkTimePresence(double time) const
{
  if(time<_start_time-_time_tolerance || time>_end_time+_time_tolerance)
    {
      std::ostringstream stream;
      stream << "The field is defined between times " << _start_time << " and " << _end_time << " worderh tolerance ";
      stream << _time_tolerance << " and trying to access on time = " << time;
      throw INTERP_KERNEL::Exception(stream.str().c_str());
    }
}

MEDCouplingTimeDiscretization *MEDCouplingConstOnTimeInterval::aggregate(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("ConstOnTimeInterval::aggregation on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Aggregate(getArray(),other->getArray());
  MEDCouplingConstOnTimeInterval *ret=new MEDCouplingConstOnTimeInterval;
  ret->setArray(arr,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingConstOnTimeInterval::aggregate(const std::vector<const MEDCouplingTimeDiscretization *>& other) const
{
  std::vector<const DataArrayDouble *> a(other.size());
  int i=0;
  for(std::vector<const MEDCouplingTimeDiscretization *>::const_iterator it=other.begin();it!=other.end();it++,i++)
    {
      const MEDCouplingConstOnTimeInterval *itC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(*it);
      if(!itC)
        throw INTERP_KERNEL::Exception("ConstOnTimeInterval::aggregate on mismatched time discretization !");
      a[i]=itC->getArray();
    }
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Aggregate(a);
  MEDCouplingConstOnTimeInterval *ret=new MEDCouplingConstOnTimeInterval;
  ret->setArray(arr,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingConstOnTimeInterval::meld(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("ConstOnTimeInterval::meld on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Meld(getArray(),other->getArray());
  MEDCouplingConstOnTimeInterval *ret=new MEDCouplingConstOnTimeInterval;
  ret->setTimeTolerance(getTimeTolerance());
  ret->setArray(arr,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingConstOnTimeInterval::dot(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("ConstOnTimeInterval::dot on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Dot(getArray(),other->getArray());
  MEDCouplingConstOnTimeInterval *ret=new MEDCouplingConstOnTimeInterval;
  ret->setArray(arr,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingConstOnTimeInterval::crossProduct(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("ConstOnTimeInterval::crossProduct on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::CrossProduct(getArray(),other->getArray());
  MEDCouplingConstOnTimeInterval *ret=new MEDCouplingConstOnTimeInterval;
  ret->setArray(arr,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingConstOnTimeInterval::max(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("ConstOnTimeInterval::max on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Max(getArray(),other->getArray());
  MEDCouplingConstOnTimeInterval *ret=new MEDCouplingConstOnTimeInterval;
  ret->setArray(arr,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingConstOnTimeInterval::min(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("ConstOnTimeInterval::min on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Min(getArray(),other->getArray());
  MEDCouplingConstOnTimeInterval *ret=new MEDCouplingConstOnTimeInterval;
  ret->setArray(arr,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingConstOnTimeInterval::add(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("ConstOnTimeInterval::add on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Add(getArray(),other->getArray());
  MEDCouplingConstOnTimeInterval *ret=new MEDCouplingConstOnTimeInterval;
  ret->setArray(arr,0);
  int tmp1,tmp2;
  double tmp3=getStartTime(tmp1,tmp2);
  ret->setStartTime(tmp3,tmp1,tmp2);
  tmp3=getEndTime(tmp1,tmp2);
  ret->setEndTime(tmp3,tmp1,tmp2);
  return ret;
}

void MEDCouplingConstOnTimeInterval::addEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("ConstOnTimeInterval::addEqual on mismatched time discretization !");
  if(!getArray())
    throw INTERP_KERNEL::Exception("MEDCouplingConstOnTimeInterval::substractaddEqual : Data Array is NULL !");
  getArray()->addEqual(other->getArray());
}

MEDCouplingTimeDiscretization *MEDCouplingConstOnTimeInterval::substract(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("ConstOnTimeInterval::substract on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Substract(getArray(),other->getArray());
  MEDCouplingConstOnTimeInterval *ret=new MEDCouplingConstOnTimeInterval;
  ret->setArray(arr,0);
  int tmp1,tmp2;
  double tmp3=getStartTime(tmp1,tmp2);
  ret->setStartTime(tmp3,tmp1,tmp2);
  tmp3=getEndTime(tmp1,tmp2);
  ret->setEndTime(tmp3,tmp1,tmp2);
  return ret;
}

void MEDCouplingConstOnTimeInterval::substractEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("ConstOnTimeInterval::substractEqual on mismatched time discretization !");
  if(!getArray())
    throw INTERP_KERNEL::Exception("MEDCouplingConstOnTimeInterval::substractEqual : Data Array is NULL !");
  getArray()->substractEqual(other->getArray());
}

MEDCouplingTimeDiscretization *MEDCouplingConstOnTimeInterval::multiply(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("multiply on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Multiply(getArray(),other->getArray());
  MEDCouplingConstOnTimeInterval *ret=new MEDCouplingConstOnTimeInterval;
  ret->setArray(arr,0);
  int tmp1,tmp2;
  double tmp3=getStartTime(tmp1,tmp2);
  ret->setStartTime(tmp3,tmp1,tmp2);
  tmp3=getEndTime(tmp1,tmp2);
  ret->setEndTime(tmp3,tmp1,tmp2);
  return ret;
}

void MEDCouplingConstOnTimeInterval::multiplyEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("ConstOnTimeInterval::multiplyEqual on mismatched time discretization !");
  if(!getArray())
    throw INTERP_KERNEL::Exception("MEDCouplingConstOnTimeInterval::multiplyEqual : Data Array is NULL !");
  getArray()->multiplyEqual(other->getArray());
}

MEDCouplingTimeDiscretization *MEDCouplingConstOnTimeInterval::divide(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("divide on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Divide(getArray(),other->getArray());
  MEDCouplingConstOnTimeInterval *ret=new MEDCouplingConstOnTimeInterval;
  ret->setArray(arr,0);
  int tmp1,tmp2;
  double tmp3=getStartTime(tmp1,tmp2);
  ret->setStartTime(tmp3,tmp1,tmp2);
  tmp3=getEndTime(tmp1,tmp2);
  ret->setEndTime(tmp3,tmp1,tmp2);
  return ret;
}

void MEDCouplingConstOnTimeInterval::divideEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("ConstOnTimeInterval::divideEqual on mismatched time discretization !");
  if(!getArray())
    throw INTERP_KERNEL::Exception("MEDCouplingConstOnTimeInterval::divideEqual : Data Array is NULL !");
  getArray()->divideEqual(other->getArray());
}

MEDCouplingTimeDiscretization *MEDCouplingConstOnTimeInterval::pow(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("pow on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Pow(getArray(),other->getArray());
  MEDCouplingConstOnTimeInterval *ret=new MEDCouplingConstOnTimeInterval;
  ret->setArray(arr,0);
  int tmp1,tmp2;
  double tmp3=getStartTime(tmp1,tmp2);
  ret->setStartTime(tmp3,tmp1,tmp2);
  tmp3=getEndTime(tmp1,tmp2);
  ret->setEndTime(tmp3,tmp1,tmp2);
  return ret;
}

void MEDCouplingConstOnTimeInterval::powEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("ConstOnTimeInterval::powEqual on mismatched time discretization !");
  if(!getArray())
    throw INTERP_KERNEL::Exception("MEDCouplingConstOnTimeInterval::powEqual : Data Array is NULL !");
  getArray()->powEqual(other->getArray());
}

MEDCouplingTwoTimeSteps::MEDCouplingTwoTimeSteps(const MEDCouplingTwoTimeSteps& other, bool deepCpy):MEDCouplingTimeDiscretization(other,deepCpy),
    _start_time(other._start_time),_end_time(other._end_time),
    _start_iteration(other._start_iteration),_end_iteration(other._end_iteration),
    _start_order(other._start_order),_end_order(other._end_order)
{
  if(other._end_array)
    _end_array=other._end_array->performCpy(deepCpy);
  else
    _end_array=0;
}

void MEDCouplingTwoTimeSteps::updateTime() const
{
  MEDCouplingTimeDiscretization::updateTime();
  if(_end_array)
    updateTimeWith(*_end_array);
}

void MEDCouplingTwoTimeSteps::synchronizeTimeWith(const MEDCouplingMesh *mesh)
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingTwoTimeSteps::synchronizeTimeWith : mesh instance is NULL ! Impossible to synchronize time !");
  int it=-1,order=-1;
  double val=mesh->getTime(it,order);
  _start_time=val; _start_iteration=it; _start_order=order;
  _end_time=val; _end_iteration=it; _end_order=order;
  std::string tUnit=mesh->getTimeUnit();
  _time_unit=tUnit;
}

std::size_t MEDCouplingTwoTimeSteps::getHeapMemorySizeWithoutChildren() const
{
  return MEDCouplingTimeDiscretization::getHeapMemorySizeWithoutChildren();
}

std::vector<const BigMemoryObject *> MEDCouplingTwoTimeSteps::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret(MEDCouplingTimeDiscretization::getDirectChildrenWithNull());
  ret.push_back(_end_array);
  return ret;
}

void MEDCouplingTwoTimeSteps::copyTinyAttrFrom(const MEDCouplingTimeDiscretization& other)
{
  MEDCouplingTimeDiscretization::copyTinyAttrFrom(other);
  const MEDCouplingTwoTimeSteps *otherC=dynamic_cast<const MEDCouplingTwoTimeSteps *>(&other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("MEDCouplingTwoTimeSteps::copyTinyAttrFrom : mismatch of time discretization !");
  _start_time=otherC->_start_time;
  _end_time=otherC->_end_time;
  _start_iteration=otherC->_start_iteration;
  _end_iteration=otherC->_end_iteration;
  _start_order=otherC->_start_order;
  _end_order=otherC->_end_order;
}

void MEDCouplingTwoTimeSteps::copyTinyStringsFrom(const MEDCouplingTimeDiscretization& other)
{
  MEDCouplingTimeDiscretization::copyTinyStringsFrom(other);
  const MEDCouplingTwoTimeSteps *otherC=dynamic_cast<const MEDCouplingTwoTimeSteps *>(&other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("Trying to operate copyTinyStringsFrom on different field type (two times//one time) !");
  if(_end_array && otherC->_end_array)
    _end_array->copyStringInfoFrom(*otherC->_end_array);
}

const DataArrayDouble *MEDCouplingTwoTimeSteps::getEndArray() const
{
  return _end_array;
}

DataArrayDouble *MEDCouplingTwoTimeSteps::getEndArray()
{
  return _end_array;
}

void MEDCouplingTwoTimeSteps::checkCoherency() const
{
  MEDCouplingTimeDiscretization::checkCoherency();
  if(!_end_array)
    throw INTERP_KERNEL::Exception("No end array specified !");
  if(_array->getNumberOfComponents()!=_end_array->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("The number of components mismatch between the start and the end arrays !");
  if(_array->getNumberOfTuples()!=_end_array->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("The number of tuples mismatch between the start and the end arrays !");
}

bool MEDCouplingTwoTimeSteps::isEqualIfNotWhy(const MEDCouplingTimeDiscretization *other, double prec, std::string& reason) const
{
  std::ostringstream oss;
  const MEDCouplingTwoTimeSteps *otherC=dynamic_cast<const MEDCouplingTwoTimeSteps *>(other);
  if(!otherC)
    {
      reason="This has time discretization LINEAR_TIME, other not.";
      return false;
    }
  if(_start_iteration!=otherC->_start_iteration)
    {
      oss << "start iterations differ. this start iteration=" << _start_iteration << " other start iteration=" << otherC->_start_iteration;
      reason=oss.str();
      return false;
    }
  if(_start_order!=otherC->_start_order)
    {
      oss << "start orders differ. this start order=" << _start_order << " other start order=" << otherC->_start_order;
      reason=oss.str();
      return false;
    }
  if(std::fabs(_start_time-otherC->_start_time)>_time_tolerance)
    {
      oss << "start times differ. this start time=" << _start_time << " other start time=" << otherC->_start_time;
      reason=oss.str();
      return false;
    }
  if(_end_iteration!=otherC->_end_iteration)
    {
      oss << "end iterations differ. this end iteration=" << _end_iteration << " other end iteration=" << otherC->_end_iteration;
      reason=oss.str();
      return false;
    }
  if(_end_order!=otherC->_end_order)
    {
      oss << "end orders differ. this end order=" << _end_order << " other end order=" << otherC->_end_order;
      reason=oss.str();
      return false;
    }
  if(std::fabs(_end_time-otherC->_end_time)>_time_tolerance)
    {
      oss << "end times differ. this end time=" << _end_time << " other end time=" << otherC->_end_time;
      reason=oss.str();
      return false;
    }
  if(_end_array!=otherC->_end_array)
    if(!_end_array->isEqualIfNotWhy(*otherC->_end_array,prec,reason))
      {
        reason.insert(0,"end arrays differ for linear time.");
        return false;
      }
  return MEDCouplingTimeDiscretization::isEqualIfNotWhy(other,prec,reason);
}

bool MEDCouplingTwoTimeSteps::isEqualWithoutConsideringStr(const MEDCouplingTimeDiscretization *other, double prec) const
{
  const MEDCouplingTwoTimeSteps *otherC=dynamic_cast<const MEDCouplingTwoTimeSteps *>(other);
  if(!otherC)
    return false;
  if(_start_iteration!=otherC->_start_iteration)
    return false;
  if(_end_iteration!=otherC->_end_iteration)
    return false;
  if(_start_order!=otherC->_start_order)
    return false;
  if(_end_order!=otherC->_end_order)
    return false;
  if(std::fabs(_start_time-otherC->_start_time)>_time_tolerance)
    return false;
  if(std::fabs(_end_time-otherC->_end_time)>_time_tolerance)
    return false;
  if(_end_array!=otherC->_end_array)
    if(!_end_array->isEqualWithoutConsideringStr(*otherC->_end_array,prec))
      return false;
  return MEDCouplingTimeDiscretization::isEqualWithoutConsideringStr(other,prec);
}

MEDCouplingTwoTimeSteps::MEDCouplingTwoTimeSteps():_start_time(0.),_end_time(0.),_start_iteration(-1),_end_iteration(-1),_start_order(-1),_end_order(-1),_end_array(0)
{
}

MEDCouplingTwoTimeSteps::~MEDCouplingTwoTimeSteps()
{
  if(_end_array)
    _end_array->decrRef();
}

void MEDCouplingTwoTimeSteps::checkNoTimePresence() const
{
  throw INTERP_KERNEL::Exception("The field presents a time to be specified in every access !");
}

void MEDCouplingTwoTimeSteps::checkTimePresence(double time) const
{
  if(time<_start_time-_time_tolerance || time>_end_time+_time_tolerance)
    {
      std::ostringstream stream;
      stream << "The field is defined between times " << _start_time << " and " << _end_time << " worderh tolerance ";
      stream << _time_tolerance << " and trying to access on time = " << time;
      throw INTERP_KERNEL::Exception(stream.str().c_str());
    }
}

void MEDCouplingTwoTimeSteps::getArrays(std::vector<DataArrayDouble *>& arrays) const
{
  arrays.resize(2);
  arrays[0]=_array;
  arrays[1]=_end_array;
}

void MEDCouplingTwoTimeSteps::setEndArray(DataArrayDouble *array, TimeLabel *owner)
{
  if(array!=_end_array)
    {
      if(_end_array)
        _end_array->decrRef();
      _end_array=array;
      if(_end_array)
        _end_array->incrRef();
      if(owner)
        owner->declareAsNew();
    }
}

void MEDCouplingTwoTimeSteps::getTinySerializationIntInformation(std::vector<int>& tinyInfo) const
{
  MEDCouplingTimeDiscretization::getTinySerializationIntInformation(tinyInfo);
  tinyInfo.push_back(_start_iteration);
  tinyInfo.push_back(_start_order);
  tinyInfo.push_back(_end_iteration);
  tinyInfo.push_back(_end_order);
  if(_end_array)
    {
      tinyInfo.push_back(_end_array->getNumberOfTuples());
      tinyInfo.push_back(_end_array->getNumberOfComponents());
    }
  else
    {
      tinyInfo.push_back(-1);
      tinyInfo.push_back(-1);
    }
}

void MEDCouplingTwoTimeSteps::getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const
{
  MEDCouplingTimeDiscretization::getTinySerializationDbleInformation(tinyInfo);
  tinyInfo.push_back(_start_time);
  tinyInfo.push_back(_end_time);
}

void MEDCouplingTwoTimeSteps::getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const
{
  int nbOfCompo=_array->getNumberOfComponents();
  for(int i=0;i<nbOfCompo;i++)
    tinyInfo.push_back(_array->getInfoOnComponent(i));
  for(int i=0;i<nbOfCompo;i++)
    tinyInfo.push_back(_end_array->getInfoOnComponent(i));
}

void MEDCouplingTwoTimeSteps::resizeForUnserialization(const std::vector<int>& tinyInfoI, std::vector<DataArrayDouble *>& arrays)
{
  arrays.resize(2);
  if(_array!=0)
    _array->decrRef();
  if(_end_array!=0)
    _end_array->decrRef();
  DataArrayDouble *arr=0;
  if(tinyInfoI[0]!=-1 && tinyInfoI[1]!=-1)
    {
      arr=DataArrayDouble::New();
      arr->alloc(tinyInfoI[0],tinyInfoI[1]);
    }
  _array=arr;
  arrays[0]=arr;
  arr=0;
  if(tinyInfoI[6]!=-1 && tinyInfoI[7]!=-1)
    {
      arr=DataArrayDouble::New();
      arr->alloc(tinyInfoI[6],tinyInfoI[7]);
    }
  _end_array=arr;
  arrays[1]=arr;
}

void MEDCouplingTwoTimeSteps::checkForUnserialization(const std::vector<int>& tinyInfoI, const std::vector<DataArrayDouble *>& arrays)
{
  static const char MSG[]="MEDCouplingTimeDiscretization::checkForUnserialization : arrays in input is expected to have size two !";
  if(arrays.size()!=2)
    throw INTERP_KERNEL::Exception(MSG);
  if(_array!=0)
    _array->decrRef();
  if(_end_array!=0)
    _end_array->decrRef();
  _array=0; _end_array=0;
  if(tinyInfoI[0]!=-1 && tinyInfoI[1]!=-1)
    {
      if(!arrays[0])
        throw INTERP_KERNEL::Exception(MSG);
      arrays[0]->checkNbOfTuplesAndComp(tinyInfoI[0],tinyInfoI[1],MSG);
      _array=arrays[0]; _array->incrRef();
    }
  if(tinyInfoI[6]!=-1 && tinyInfoI[7]!=-1)
    {
      if(!arrays[1])
        throw INTERP_KERNEL::Exception(MSG);
      arrays[1]->checkNbOfTuplesAndComp(tinyInfoI[0],tinyInfoI[1],MSG);
      _end_array=arrays[1]; _end_array->incrRef();
    }
}

void MEDCouplingTwoTimeSteps::finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS)
{
  MEDCouplingTimeDiscretization::finishUnserialization(tinyInfoI,tinyInfoD,tinyInfoS);
  _start_time=tinyInfoD[1];
  _end_time=tinyInfoD[2];
  _start_iteration=tinyInfoI[2];
  _start_order=tinyInfoI[3];
  _end_iteration=tinyInfoI[4];
  _end_order=tinyInfoI[5];
}

/*!
 * idem getTinySerializationIntInformation except that it is for multi field fetch
 */
void MEDCouplingTwoTimeSteps::getTinySerializationIntInformation2(std::vector<int>& tinyInfo) const
{
  tinyInfo.resize(4);
  tinyInfo[0]=_start_iteration;
  tinyInfo[1]=_start_order;
  tinyInfo[2]=_end_iteration;
  tinyInfo[3]=_end_order;
}

/*!
 * idem getTinySerializationDbleInformation except that it is for multi field fetch
 */
void MEDCouplingTwoTimeSteps::getTinySerializationDbleInformation2(std::vector<double>& tinyInfo) const
{
  tinyInfo.resize(3);
  tinyInfo[0]=_time_tolerance;
  tinyInfo[1]=_start_time;
  tinyInfo[2]=_end_time;
}

/*!
 * idem finishUnserialization except that it is for multi field fetch
 */
void MEDCouplingTwoTimeSteps::finishUnserialization2(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD)
{
  _start_iteration=tinyInfoI[0];
  _start_order=tinyInfoI[1];
  _end_iteration=tinyInfoI[2];
  _end_order=tinyInfoI[3];
  _time_tolerance=tinyInfoD[0];
  _start_time=tinyInfoD[1];
  _end_time=tinyInfoD[2];
}

std::vector< const DataArrayDouble *> MEDCouplingTwoTimeSteps::getArraysForTime(double time) const
{
  if(time>_start_time-_time_tolerance && time<_end_time+_time_tolerance)
    {
      std::vector< const DataArrayDouble *> ret(2);
      ret[0]=_array;
      ret[1]=_end_array;
      return ret;
    }
  else
    throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingTwoTimeSteps::setArrays(const std::vector<DataArrayDouble *>& arrays, TimeLabel *owner)
{
  if(arrays.size()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingTwoTimeSteps::setArrays : number of arrays must be two.");
  setArray(arrays.front(),owner);
  setEndArray(arrays.back(),owner);
}

MEDCouplingLinearTime::MEDCouplingLinearTime(const MEDCouplingLinearTime& other, bool deepCpy):MEDCouplingTwoTimeSteps(other,deepCpy)
{
}

MEDCouplingLinearTime::MEDCouplingLinearTime()
{
}

std::string MEDCouplingLinearTime::getStringRepr() const
{
  std::ostringstream stream;
  stream << REPR << " Time interval is defined by :\niteration_start=" << _start_iteration << " order_start=" << _start_order << " and time_start=" << _start_time << "\n";
  stream << "iteration_end=" << _end_iteration << " order_end=" << _end_order << " and end_time=" << _end_time << "\n";
  stream << "Time unit is : \"" << _time_unit << "\"";
  return stream.str();
}

void MEDCouplingLinearTime::checkCoherency() const
{
  MEDCouplingTwoTimeSteps::checkCoherency();
  if(std::fabs(_start_time-_end_time)<_time_tolerance)
    throw INTERP_KERNEL::Exception("Start time and end time are equals regarding time tolerance.");
}

MEDCouplingTimeDiscretization *MEDCouplingLinearTime::performCpy(bool deepCpy) const
{
  return new MEDCouplingLinearTime(*this,deepCpy);
}

bool MEDCouplingLinearTime::areCompatible(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areCompatible(other))
    return false;
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(otherC==0)
    return false;
  if(_end_array==0 && otherC->_end_array==0)
    return true;
  if(_end_array==0 || otherC->_end_array==0)
    return false;
  if(_end_array->getNumberOfComponents()!=otherC->_end_array->getNumberOfComponents())
    return false;
  return true;
}

bool MEDCouplingLinearTime::areStrictlyCompatible(const MEDCouplingTimeDiscretization *other, std::string& reason) const
{
  if(!MEDCouplingTimeDiscretization::areStrictlyCompatible(other,reason))
    return false;
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  bool ret=otherC!=0;
  if(!ret)
    reason.insert(0,"time discretization of this is LINEAR_TIME, other has a different time discretization.");
  return ret;
}

bool MEDCouplingLinearTime::areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areStrictlyCompatibleForMul(other))
    return false;
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  return otherC!=0;
}

bool MEDCouplingLinearTime::areStrictlyCompatibleForDiv(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areStrictlyCompatibleForDiv(other))
    return false;
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(otherC==0)
    return false;
  if(_end_array==0 && otherC->_end_array==0)
    return true;
  if(_end_array==0 || otherC->_end_array==0)
    return false;
  int nbC1=_end_array->getNumberOfComponents();
  int nbC2=otherC->_end_array->getNumberOfComponents();
  if(nbC1!=nbC2 && nbC2!=1)
    return false;
  return true;
}

bool MEDCouplingLinearTime::areCompatibleForMeld(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areCompatibleForMeld(other))
    return false;
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  return otherC!=0;
}

/*!
 * vals is expected to be of size 2*_array->getNumberOfTuples()==_array->getNumberOfTuples()+_end_array->getNumberOfTuples()
 */
void MEDCouplingLinearTime::getValueForTime(double time, const std::vector<double>& vals, double *res) const
{
  double alpha=(_end_time-time)/(_end_time-_start_time);
  std::size_t nbComp=vals.size()/2;
  std::transform(vals.begin(),vals.begin()+nbComp,res,std::bind2nd(std::multiplies<double>(),alpha));
  std::vector<double> tmp(nbComp);
  std::transform(vals.begin()+nbComp,vals.end(),tmp.begin(),std::bind2nd(std::multiplies<double>(),1-alpha));
  std::transform(tmp.begin(),tmp.end(),res,res,std::plus<double>());
}

void MEDCouplingLinearTime::getValueOnTime(int eltId, double time, double *value) const
{
  double alpha=(_end_time-time)/(_end_time-_start_time);
  int nbComp;
  if(_array)
    _array->getTuple(eltId,value);
  else
    throw INTERP_KERNEL::Exception("No start array existing.");
  nbComp=_array->getNumberOfComponents();
  std::transform(value,value+nbComp,value,std::bind2nd(std::multiplies<double>(),alpha));
  std::vector<double> tmp(nbComp);
  if(_end_array)
    _end_array->getTuple(eltId,&tmp[0]);
  else
    throw INTERP_KERNEL::Exception("No end array existing.");
  std::transform(tmp.begin(),tmp.end(),tmp.begin(),std::bind2nd(std::multiplies<double>(),1-alpha));
  std::transform(tmp.begin(),tmp.end(),value,value,std::plus<double>());
}

void MEDCouplingLinearTime::getValueOnDiscTime(int eltId, int iteration, int order, double *value) const
{
  if(iteration==_start_iteration && order==_start_order)
    {
      if(_array)
        _array->getTuple(eltId,value);
      else
        throw INTERP_KERNEL::Exception("iteration order match with start time but no start array existing.");
    }
  if(iteration==_end_iteration && order==_end_order)
    {
      if(_end_array)
        _end_array->getTuple(eltId,value);
      else
        throw INTERP_KERNEL::Exception("iteration order match with end time but no end array existing.");
    }
  else
    throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

MEDCouplingTimeDiscretization *MEDCouplingLinearTime::aggregate(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("LinearTime::aggregation on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr1=DataArrayDouble::Aggregate(getArray(),other->getArray());
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr2=DataArrayDouble::Aggregate(getEndArray(),other->getEndArray());
  MEDCouplingLinearTime *ret=new MEDCouplingLinearTime;
  ret->setArray(arr1,0);
  ret->setEndArray(arr2,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingLinearTime::aggregate(const std::vector<const MEDCouplingTimeDiscretization *>& other) const
{
  std::vector<const DataArrayDouble *> a(other.size());
  std::vector<const DataArrayDouble *> b(other.size());
  int i=0;
  for(std::vector<const MEDCouplingTimeDiscretization *>::const_iterator it=other.begin();it!=other.end();it++,i++)
    {
      const MEDCouplingLinearTime *itC=dynamic_cast<const MEDCouplingLinearTime *>(*it);
      if(!itC)
        throw INTERP_KERNEL::Exception("MEDCouplingLinearTime::aggregate on mismatched time discretization !");
      a[i]=itC->getArray();
      b[i]=itC->getEndArray();
    }
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::Aggregate(a);
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr2=DataArrayDouble::Aggregate(b);
  MEDCouplingLinearTime *ret=new MEDCouplingLinearTime;
  ret->setArray(arr,0);
  ret->setEndArray(arr2,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingLinearTime::meld(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("LinearTime::meld on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr1=DataArrayDouble::Meld(getArray(),other->getArray());
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr2=DataArrayDouble::Meld(getEndArray(),other->getEndArray());
  MEDCouplingLinearTime *ret=new MEDCouplingLinearTime;
  ret->setTimeTolerance(getTimeTolerance());
  ret->setArray(arr1,0);
  ret->setEndArray(arr2,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingLinearTime::dot(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("LinearTime::dot on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr1=DataArrayDouble::Dot(getArray(),other->getArray());
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr2=DataArrayDouble::Dot(getEndArray(),other->getEndArray());
  MEDCouplingLinearTime *ret=new MEDCouplingLinearTime;
  ret->setArray(arr1,0);
  ret->setEndArray(arr2,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingLinearTime::crossProduct(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("LinearTime::crossProduct on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr1=DataArrayDouble::CrossProduct(getArray(),other->getArray());
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr2=DataArrayDouble::CrossProduct(getEndArray(),other->getEndArray());
  MEDCouplingLinearTime *ret=new MEDCouplingLinearTime;
  ret->setArray(arr1,0);
  ret->setEndArray(arr2,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingLinearTime::max(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("LinearTime::max on mismatched time discretization !");
  MEDCouplingLinearTime *ret=new MEDCouplingLinearTime;
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr1=DataArrayDouble::Max(getArray(),other->getArray());
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr2=DataArrayDouble::Max(getEndArray(),other->getEndArray());
  ret->setArray(arr1,0);
  ret->setEndArray(arr2,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingLinearTime::min(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("LinearTime::min on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr1=DataArrayDouble::Min(getArray(),other->getArray());
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr2=DataArrayDouble::Min(getEndArray(),other->getEndArray());
  MEDCouplingLinearTime *ret=new MEDCouplingLinearTime;
  ret->setArray(arr1,0);
  ret->setEndArray(arr2,0);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingLinearTime::add(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("LinearTime::add on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr1=DataArrayDouble::Add(getArray(),other->getArray());
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr2=DataArrayDouble::Add(getEndArray(),other->getEndArray());
  MEDCouplingLinearTime *ret=new MEDCouplingLinearTime;
  ret->setArray(arr1,0);
  ret->setEndArray(arr2,0);
  return ret;
}

void MEDCouplingLinearTime::addEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("LinearTime::addEqual on mismatched time discretization !");
  if(!getArray())
    throw INTERP_KERNEL::Exception("MEDCouplingLinearTime::addEqual : Data Array is NULL !");
  if(!getEndArray())
    throw INTERP_KERNEL::Exception("MEDCouplingLinearTime::addEqual : Data Array (end) is NULL !");
  getArray()->addEqual(other->getArray());
  getEndArray()->addEqual(other->getEndArray());
}

MEDCouplingTimeDiscretization *MEDCouplingLinearTime::substract(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("LinearTime::substract on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr1=DataArrayDouble::Substract(getArray(),other->getArray());
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr2=DataArrayDouble::Substract(getEndArray(),other->getEndArray());
  MEDCouplingLinearTime *ret=new MEDCouplingLinearTime;
  ret->setArray(arr1,0);
  ret->setEndArray(arr2,0);
  return ret;
}

void MEDCouplingLinearTime::substractEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("LinearTime::addEqual on mismatched time discretization !");
  if(!getArray())
    throw INTERP_KERNEL::Exception("MEDCouplingLinearTime::substractEqual : Data Array is NULL !");
  if(!getEndArray())
    throw INTERP_KERNEL::Exception("MEDCouplingLinearTime::substractEqual : Data Array (end) is NULL !");
  getArray()->substractEqual(other->getArray());
  getEndArray()->substractEqual(other->getEndArray());
}

MEDCouplingTimeDiscretization *MEDCouplingLinearTime::multiply(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("LinearTime::multiply on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr1=DataArrayDouble::Multiply(getArray(),other->getArray());
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr2=DataArrayDouble::Multiply(getEndArray(),other->getEndArray());
  MEDCouplingLinearTime *ret=new MEDCouplingLinearTime;
  ret->setArray(arr1,0);
  ret->setEndArray(arr2,0);
  return ret;
}

void MEDCouplingLinearTime::multiplyEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("LinearTime::addEqual on mismatched time discretization !");
  if(!getArray())
    throw INTERP_KERNEL::Exception("MEDCouplingLinearTime::multiplyEqual : Data Array is NULL !");
  if(!getEndArray())
    throw INTERP_KERNEL::Exception("MEDCouplingLinearTime::multiplyEqual : Data Array (end) is NULL !");
  getArray()->multiplyEqual(other->getArray());
  getEndArray()->multiplyEqual(other->getEndArray());
}

MEDCouplingTimeDiscretization *MEDCouplingLinearTime::divide(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("LinearTime::divide on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr1=DataArrayDouble::Divide(getArray(),other->getArray());
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr2=DataArrayDouble::Divide(getEndArray(),other->getEndArray());
  MEDCouplingLinearTime *ret=new MEDCouplingLinearTime;
  ret->setArray(arr1,0);
  ret->setEndArray(arr2,0);
  return ret;
}

void MEDCouplingLinearTime::divideEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("LinearTime::addEqual on mismatched time discretization !");
  if(!getArray())
    throw INTERP_KERNEL::Exception("MEDCouplingLinearTime::divideEqual : Data Array is NULL !");
  if(!getEndArray())
    throw INTERP_KERNEL::Exception("MEDCouplingLinearTime::divideEqual : Data Array (end) is NULL !");
  getArray()->divideEqual(other->getArray());
  getEndArray()->divideEqual(other->getEndArray());
}

MEDCouplingTimeDiscretization *MEDCouplingLinearTime::pow(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("LinearTime::pow on mismatched time discretization !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr1=DataArrayDouble::Pow(getArray(),other->getArray());
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr2=DataArrayDouble::Pow(getEndArray(),other->getEndArray());
  MEDCouplingLinearTime *ret=new MEDCouplingLinearTime;
  ret->setArray(arr1,0);
  ret->setEndArray(arr2,0);
  return ret;
}

void MEDCouplingLinearTime::powEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("LinearTime::addEqual on mismatched time discretization !");
  if(!getArray())
    throw INTERP_KERNEL::Exception("MEDCouplingLinearTime::powEqual : Data Array is NULL !");
  if(!getEndArray())
    throw INTERP_KERNEL::Exception("MEDCouplingLinearTime::powEqual : Data Array (end) is NULL !");
  getArray()->powEqual(other->getArray());
  getEndArray()->powEqual(other->getEndArray());
}
