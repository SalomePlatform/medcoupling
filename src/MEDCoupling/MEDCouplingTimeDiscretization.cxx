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
#include "MEDCouplingTimeDiscretization.hxx"
#include "MEDCouplingMemArray.hxx"

#include <cmath>

using namespace ParaMEDMEM;

const char MEDCouplingNoTimeLabel::EXCEPTION_MSG[]="MEDCouplingNoTimeLabel::setTime : no time info attached.";

const char MEDCouplingWithTimeStep::EXCEPTION_MSG[]="No data on this time.";

const double MEDCouplingTimeDiscretization::TIME_TOLERANCE_DFT=1.e-12;

MEDCouplingTimeDiscretization *MEDCouplingTimeDiscretization::New(TypeOfTimeDiscretization type)
{
  switch(type)
    {
    case MEDCouplingNoTimeLabel::DISCRETIZATION:
      return new MEDCouplingNoTimeLabel;
    case MEDCouplingWithTimeStep::DISCRETIZATION:
      return new MEDCouplingWithTimeStep;
    default:
      throw INTERP_KERNEL::Exception("Time discretization not implemented yet");
    }
}

bool MEDCouplingTimeDiscretization::isEqual(const MEDCouplingTimeDiscretization *other, double prec) const
{
  if(std::fabs(_time_tolerance-other->_time_tolerance)>1.e-16)
    return false;
  if(_array==0 && other->_array==0)
    return true;
  if(_array==0 || other->_array==0)
    return false;
  if(_array==other->_array)
    return true;
  return _array->isEqual(*other->_array,prec);
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

void MEDCouplingTimeDiscretization::finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS)
{
  _time_tolerance=tinyInfoD[0];
  int nbOfCompo=_array->getNumberOfComponents();
  for(int i=0;i<nbOfCompo;i++)
    _array->setInfoOnComponent(i,tinyInfoS[i].c_str());
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

MEDCouplingTimeDiscretization::MEDCouplingTimeDiscretization(const MEDCouplingTimeDiscretization& other, bool deepCpy):_time_tolerance(other._time_tolerance)
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
      owner->declareAsNew();
    }
}

void MEDCouplingTimeDiscretization::setArrays(const std::vector<DataArrayDouble *>& arrays, TimeLabel *owner) throw(INTERP_KERNEL::Exception)
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

bool MEDCouplingTimeDiscretization::isBefore(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception)
{
  int dt,it;
  double time1=getEndTime(dt,it)-_time_tolerance;
  double time2=other->getStartTime(dt,it)+other->getTimeTolerance();
  return time1<=time2;
}

bool MEDCouplingTimeDiscretization::isStrictlyBefore(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception)
{
  int dt,it;
  double time1=getEndTime(dt,it)+_time_tolerance;
  double time2=other->getStartTime(dt,it)-other->getTimeTolerance();
  return time1<time2;
}

MEDCouplingNoTimeLabel::MEDCouplingNoTimeLabel()
{
}

MEDCouplingNoTimeLabel::MEDCouplingNoTimeLabel(const MEDCouplingTimeDiscretization& other, bool deepCpy):MEDCouplingTimeDiscretization(other,deepCpy)
{
}

bool MEDCouplingNoTimeLabel::isEqual(const MEDCouplingTimeDiscretization *other, double prec) const
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    return false;
  return MEDCouplingTimeDiscretization::isEqual(other,prec);
}

MEDCouplingTimeDiscretization *MEDCouplingNoTimeLabel::performCpy(bool deepCpy) const
{
  return new MEDCouplingNoTimeLabel(*this,deepCpy);
}

void MEDCouplingNoTimeLabel::checkTimePresence(double time) const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

DataArrayDouble *MEDCouplingNoTimeLabel::getArrayOnTime(double time) const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

bool MEDCouplingNoTimeLabel::isBefore(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

bool MEDCouplingNoTimeLabel::isStrictlyBefore(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

double MEDCouplingNoTimeLabel::getStartTime(int& dt, int& it) const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

double MEDCouplingNoTimeLabel::getEndTime(int& dt, int& it) const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingNoTimeLabel::setStartTime(double time, int dt, int it) throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingNoTimeLabel::setEndTime(double time, int dt, int it) throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingNoTimeLabel::getValueOnTime(int eltId, double time, double *value) const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingNoTimeLabel::getValueOnDiscTime(int eltId, int dt, int it, double *value) const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

MEDCouplingWithTimeStep::MEDCouplingWithTimeStep(const MEDCouplingWithTimeStep& other, bool deepCpy):MEDCouplingTimeDiscretization(other,deepCpy),
                                                                                                     _time(other._time),_dt(other._dt),_it(other._it)
{
}

MEDCouplingWithTimeStep::MEDCouplingWithTimeStep():_time(0.),_dt(-1),_it(-1)
{
}

void MEDCouplingWithTimeStep::getTinySerializationIntInformation(std::vector<int>& tinyInfo) const
{
  MEDCouplingTimeDiscretization::getTinySerializationIntInformation(tinyInfo);
  tinyInfo.push_back(_dt);
  tinyInfo.push_back(_it);
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
  _dt=tinyInfoI[2];
  _it=tinyInfoI[3];
}

bool MEDCouplingWithTimeStep::isEqual(const MEDCouplingTimeDiscretization *other, double prec) const
{
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  if(!otherC)
    return false;
  if(_dt!=otherC->_dt)
    return false;
  if(_it!=otherC->_it)
    return false;
  if(std::fabs(_time-otherC->_time)>_time_tolerance)
    return false;
  return MEDCouplingTimeDiscretization::isEqual(other,prec);
}

MEDCouplingTimeDiscretization *MEDCouplingWithTimeStep::performCpy(bool deepCpy) const
{
  return new MEDCouplingWithTimeStep(*this,deepCpy);
}

void MEDCouplingWithTimeStep::checkNoTimePresence() const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("No time specified on a field defined on one time");
}

void MEDCouplingWithTimeStep::checkTimePresence(double time) const throw(INTERP_KERNEL::Exception)
{
  if(std::fabs(time-_time)>_time_tolerance)
    {
      std::ostringstream stream;
      stream << "The field is defined on time " << _time << " with eps=" << _time_tolerance << " and asking time = " << time << " !";
      throw INTERP_KERNEL::Exception(stream.str().c_str());
    }
}

DataArrayDouble *MEDCouplingWithTimeStep::getArrayOnTime(double time) const throw(INTERP_KERNEL::Exception)
{
  if(std::fabs(time-_time)<=_time_tolerance)
    {
      if(_array)
        _array->incrRef();
      return _array;
    }
  else
    throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingWithTimeStep::getValueOnTime(int eltId, double time, double *value) const throw(INTERP_KERNEL::Exception)
{
  if(std::fabs(time-_time)<=_time_tolerance)
    if(_array)
      _array->getTuple(eltId,value);
    else
      throw INTERP_KERNEL::Exception("No array existing.");
  else
    throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingWithTimeStep::getValueOnDiscTime(int eltId, int dt, int it, double *value) const throw(INTERP_KERNEL::Exception)
{
  if(_dt==dt && _it==it)
    if(_array)
      _array->getTuple(eltId,value);
    else
      throw INTERP_KERNEL::Exception("No array existing.");
  else
    throw INTERP_KERNEL::Exception("No data on this discrete time.");
}

MEDCouplingTwoTimeSteps::MEDCouplingTwoTimeSteps():_start_time(0.),_end_time(0.),_start_dt(-1),_end_dt(-1),_start_it(-1),_end_it(-1),_end_array(0)
{
}

MEDCouplingTwoTimeSteps::~MEDCouplingTwoTimeSteps()
{
  if(_end_array)
    _end_array->decrRef();
}

void MEDCouplingTwoTimeSteps::checkNoTimePresence() const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("The field presents a time to be specified in every access !");
}

void MEDCouplingTwoTimeSteps::checkTimePresence(double time) const throw(INTERP_KERNEL::Exception)
{
  if(time<_start_time-_time_tolerance || time>_end_time+_time_tolerance)
    {
      std::ostringstream stream;
      stream << "The field is defined between times " << _start_time << " and " << _end_time << " with tolerance ";
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
  if(tinyInfoI[2]!=-1 && tinyInfoI[3]!=-1)
    {
      arr=DataArrayDouble::New();
      arr->alloc(tinyInfoI[2],tinyInfoI[3]);
    }
  _end_array=arr;
  arrays[1]=arr;
}

void MEDCouplingTwoTimeSteps::finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS)
{
  MEDCouplingTimeDiscretization::finishUnserialization(tinyInfoI,tinyInfoD,tinyInfoS);
  _start_time=tinyInfoD[1];
  _end_time=tinyInfoD[2];
  _start_dt=tinyInfoI[2];
  _end_dt=tinyInfoI[3];
  _start_it=tinyInfoI[4];
  _end_it=tinyInfoI[5];
}
