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

#include "MEDCouplingTimeDiscretization.hxx"
#include "MEDCouplingMemArray.hxx"

#include "InterpKernelExprParser.hxx"

#include <cmath>
#include <iterator>

using namespace ParaMEDMEM;

const char MEDCouplingNoTimeLabel::EXCEPTION_MSG[]="MEDCouplingNoTimeLabel::setTime : no time info attached.";

const char MEDCouplingWithTimeStep::EXCEPTION_MSG[]="No data on this time.";

const char MEDCouplingConstOnTimeInterval::EXCEPTION_MSG[]="No data on this time.";

const char MEDCouplingTwoTimeSteps::EXCEPTION_MSG[]="No data on this time.";

const double MEDCouplingTimeDiscretization::TIME_TOLERANCE_DFT=1.e-12;

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
}

void MEDCouplingTimeDiscretization::copyTinyStringsFrom(const MEDCouplingTimeDiscretization& other)
{
  if(_array && other._array)
    _array->copyStringInfoFrom(*other._array);
}

void MEDCouplingTimeDiscretization::checkCoherency() const throw(INTERP_KERNEL::Exception)
{
  if(!_array)
    throw INTERP_KERNEL::Exception("Field invalid because no values set !");
  if(_time_tolerance<0.)
    throw INTERP_KERNEL::Exception("time tolerance is expected to be greater than 0. !");
}

void MEDCouplingTimeDiscretization::updateTime()
{
  if(_array)
    updateTimeWith(*_array);
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

bool MEDCouplingTimeDiscretization::areStrictlyCompatible(const MEDCouplingTimeDiscretization *other) const
{
  if(std::fabs(_time_tolerance-other->_time_tolerance)>1.e-16)
    return false;
  if(_array==0 && other->_array==0)
    return true;
  if(_array==0 || other->_array==0)
    return false;
  if(_array->getNumberOfComponents()!=other->_array->getNumberOfComponents())
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

bool MEDCouplingTimeDiscretization::isEqual(const MEDCouplingTimeDiscretization *other, double prec) const
{
  if(!areStrictlyCompatible(other))
    return false;
  if(_array==other->_array)
    return true;
  return _array->isEqual(*other->_array,prec);
}

MEDCouplingTimeDiscretization *MEDCouplingTimeDiscretization::buildNewTimeReprFromThis(const MEDCouplingTimeDiscretization *other,
                                                                                       TypeOfTimeDiscretization type, bool deepCpy) const
{
  MEDCouplingTimeDiscretization *ret=MEDCouplingTimeDiscretization::New(type);
  DataArrayDouble *arrSrc=getArray();
  DataArrayDouble *arr=0;
  if(arrSrc)
    arr=arrSrc->performCpy(deepCpy);
  else
    arr=0;
  ret->setArray(arr,0);
  arr->decrRef();
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
      if(owner)
        owner->declareAsNew();
    }
}

DataArrayDouble *MEDCouplingTimeDiscretization::getEndArray() const
{
  throw INTERP_KERNEL::Exception("getEndArray not available for this type of time discretization !");
}

void MEDCouplingTimeDiscretization::setEndArray(DataArrayDouble *array, TimeLabel *owner)
{
  throw INTERP_KERNEL::Exception("setEndArray not available for this type of time discretization !");
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

void MEDCouplingTimeDiscretization::applyLin(double a, double b, int compoId)
{
  double *ptr=_array->getPointer()+compoId;
  int nbOfComp=_array->getNumberOfComponents();
  int nbOfTuple=_array->getNumberOfTuples();
  for(int i=0;i<nbOfTuple;i++,ptr+=nbOfComp)
    *ptr=a*(*ptr)+b;
}

void MEDCouplingTimeDiscretization::applyFunc(int nbOfComp, FunctionToEvaluate func)
{
  DataArrayDouble *newArr=DataArrayDouble::New();
  int nbOfTuples=_array->getNumberOfTuples();
  int oldNbOfComp=_array->getNumberOfComponents();
  newArr->alloc(nbOfTuples,nbOfComp);
  const double *ptr=_array->getConstPointer();
  double *ptrToFill=newArr->getPointer();
  for(int i=0;i<nbOfTuples;i++)
    {
      if(!func(ptr+i*oldNbOfComp,ptrToFill+i*nbOfComp))
        {
          std::ostringstream oss; oss << "For tuple # " << i << " with value (";
          std::copy(ptr+oldNbOfComp*i,ptr+oldNbOfComp*(i+1),std::ostream_iterator<double>(oss,", "));
          oss << ") : Evaluation of function failed !";
          newArr->decrRef();
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  _array->decrRef();
  _array=newArr;
}

void MEDCouplingTimeDiscretization::applyFunc(int nbOfComp, const char *func)
{
  INTERP_KERNEL::ExprParser expr(func);
  expr.parse();
  std::set<std::string> vars;
  expr.getTrueSetOfVars(vars);
  int oldNbOfComp=_array->getNumberOfComponents();
  if((int)vars.size()>oldNbOfComp)
    {
      std::ostringstream oss; oss << "The field has a " << oldNbOfComp << " components and there are ";
      oss << vars.size() << " variables : ";
      std::copy(vars.begin(),vars.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  std::vector<std::string> varsV(vars.begin(),vars.end());
  expr.prepareExprEvaluation(varsV);
  //
  DataArrayDouble *newArr=DataArrayDouble::New();
  int nbOfTuples=_array->getNumberOfTuples();
  newArr->alloc(nbOfTuples,nbOfComp);
  const double *ptr=_array->getConstPointer();
  double *ptrToFill=newArr->getPointer();
  for(int i=0;i<nbOfTuples;i++)
    {
      try
        {
          expr.evaluateExpr(nbOfComp,ptr+i*oldNbOfComp,ptrToFill+i*nbOfComp);
        }
      catch(INTERP_KERNEL::Exception& e)
        {
          std::ostringstream oss; oss << "For tuple # " << i << " with value (";
          std::copy(ptr+oldNbOfComp*i,ptr+oldNbOfComp*(i+1),std::ostream_iterator<double>(oss,", "));
          oss << ") : Evaluation of function failed !" << e.what();
          newArr->decrRef();
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  _array->decrRef();
  _array=newArr;
}

void MEDCouplingTimeDiscretization::applyFunc(const char *func)
{
  INTERP_KERNEL::ExprParser expr(func);
  expr.parse();
  expr.prepareExprEvaluationVec();
  //
  DataArrayDouble *newArr=DataArrayDouble::New();
  int nbOfTuples=_array->getNumberOfTuples();
  int nbOfComp=_array->getNumberOfComponents();
  newArr->alloc(nbOfTuples,nbOfComp);
  const double *ptr=_array->getConstPointer();
  double *ptrToFill=newArr->getPointer();
  for(int i=0;i<nbOfTuples;i++)
    {
      try
        {
          expr.evaluateExpr(nbOfComp,ptr+i*nbOfComp,ptrToFill+i*nbOfComp);
        }
      catch(INTERP_KERNEL::Exception& e)
        {
          std::ostringstream oss; oss << "For tuple # " << i << " with value (";
          std::copy(ptr+nbOfComp*i,ptr+nbOfComp*(i+1),std::ostream_iterator<double>(oss,", "));
          oss << ") : Evaluation of function failed ! " << e.what();
          newArr->decrRef();
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  _array->decrRef();
  _array=newArr;
}

MEDCouplingNoTimeLabel::MEDCouplingNoTimeLabel()
{
}

MEDCouplingNoTimeLabel::MEDCouplingNoTimeLabel(const MEDCouplingTimeDiscretization& other, bool deepCpy):MEDCouplingTimeDiscretization(other,deepCpy)
{
}

bool MEDCouplingNoTimeLabel::areCompatible(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areCompatible(other))
    return false;
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  return otherC!=0;
}

bool MEDCouplingNoTimeLabel::areStrictlyCompatible(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areStrictlyCompatible(other))
    return false;
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  return otherC!=0;
}

bool MEDCouplingNoTimeLabel::areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areStrictlyCompatibleForMul(other))
    return false;
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  return otherC!=0;
}

bool MEDCouplingNoTimeLabel::isEqual(const MEDCouplingTimeDiscretization *other, double prec) const
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    return false;
  return MEDCouplingTimeDiscretization::isEqual(other,prec);
}

MEDCouplingTimeDiscretization *MEDCouplingNoTimeLabel::aggregate(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("NoTimeLabel::aggregation on mismatched time discretization !");
  MEDCouplingNoTimeLabel *ret=new MEDCouplingNoTimeLabel;
  ret->setTimeTolerance(getTimeTolerance());
  DataArrayDouble *arr=DataArrayDouble::aggregate(getArray(),other->getArray());
  ret->setArray(arr,0);
  arr->decrRef();
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingNoTimeLabel::add(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("NoTimeLabel::add on mismatched time discretization !");
  MEDCouplingNoTimeLabel *ret=new MEDCouplingNoTimeLabel;
  DataArrayDouble *arr=DataArrayDouble::add(getArray(),other->getArray());
  ret->setArray(arr,0);
  arr->decrRef();
  return ret;
}

void MEDCouplingNoTimeLabel::addEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("NoTimeLabel::addEqual on mismatched time discretization !");
  getArray()->addEqual(other->getArray());
}

MEDCouplingTimeDiscretization *MEDCouplingNoTimeLabel::substract(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("NoTimeLabel::substract on mismatched time discretization !");
  MEDCouplingNoTimeLabel *ret=new MEDCouplingNoTimeLabel;
  DataArrayDouble *arr=DataArrayDouble::substract(getArray(),other->getArray());
  ret->setArray(arr,0);
  arr->decrRef();
  return ret;
}

void MEDCouplingNoTimeLabel::substractEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("NoTimeLabel::substractEqual on mismatched time discretization !");
  getArray()->substractEqual(other->getArray());
}

MEDCouplingTimeDiscretization *MEDCouplingNoTimeLabel::multiply(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("NoTimeLabel::multiply on mismatched time discretization !");
  MEDCouplingNoTimeLabel *ret=new MEDCouplingNoTimeLabel;
  DataArrayDouble *arr=DataArrayDouble::multiply(getArray(),other->getArray());
  ret->setArray(arr,0);
  arr->decrRef();
  return ret;
}

void MEDCouplingNoTimeLabel::multiplyEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("NoTimeLabel::multiplyEqual on mismatched time discretization !");
  getArray()->multiplyEqual(other->getArray());
}

MEDCouplingTimeDiscretization *MEDCouplingNoTimeLabel::divide(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("divide on mismatched time discretization !");
  MEDCouplingNoTimeLabel *ret=new MEDCouplingNoTimeLabel;
  DataArrayDouble *arr=DataArrayDouble::divide(getArray(),other->getArray());
  ret->setArray(arr,0);
  arr->decrRef();
  return ret;
}

void MEDCouplingNoTimeLabel::divideEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingNoTimeLabel *otherC=dynamic_cast<const MEDCouplingNoTimeLabel *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("NoTimeLabel::divideEqual on mismatched time discretization !");
  getArray()->divideEqual(other->getArray());
}

MEDCouplingTimeDiscretization *MEDCouplingNoTimeLabel::performCpy(bool deepCpy) const
{
  return new MEDCouplingNoTimeLabel(*this,deepCpy);
}

void MEDCouplingNoTimeLabel::checkTimePresence(double time) const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

std::vector< const DataArrayDouble *> MEDCouplingNoTimeLabel::getArraysForTime(double time) const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingNoTimeLabel::getValueForTime(double time, const std::vector<double>& vals, double *res) const
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

bool MEDCouplingWithTimeStep::areCompatible(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areCompatible(other))
    return false;
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  return otherC!=0;
}

bool MEDCouplingWithTimeStep::areStrictlyCompatible(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areStrictlyCompatible(other))
    return false;
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  return otherC!=0;
}

bool MEDCouplingWithTimeStep::areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areStrictlyCompatibleForMul(other))
    return false;
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  return otherC!=0;
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

void MEDCouplingWithTimeStep::copyTinyAttrFrom(const MEDCouplingTimeDiscretization& other)
{
  MEDCouplingTimeDiscretization::copyTinyAttrFrom(other);
  const MEDCouplingWithTimeStep& otherC=dynamic_cast<const MEDCouplingWithTimeStep& >(other);
  _time=otherC._time;
  _dt=otherC._dt;
  _it=otherC._it;
}

MEDCouplingTimeDiscretization *MEDCouplingWithTimeStep::aggregate(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("WithTimeStep::aggregation on mismatched time discretization !");
  MEDCouplingWithTimeStep *ret=new MEDCouplingWithTimeStep;
  ret->setTimeTolerance(getTimeTolerance());
  DataArrayDouble *arr=DataArrayDouble::aggregate(getArray(),other->getArray());
  ret->setArray(arr,0);
  arr->decrRef();
  int tmp1,tmp2;
  double tmp3=getStartTime(tmp1,tmp2);
  ret->setStartTime(tmp3,tmp1,tmp2);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingWithTimeStep::add(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("WithTimeStep::add on mismatched time discretization !");
  MEDCouplingWithTimeStep *ret=new MEDCouplingWithTimeStep;
  DataArrayDouble *arr=DataArrayDouble::add(getArray(),other->getArray());
  ret->setArray(arr,0);
  arr->decrRef();
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
  getArray()->addEqual(other->getArray());
}

MEDCouplingTimeDiscretization *MEDCouplingWithTimeStep::substract(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("WithTimeStep::substract on mismatched time discretization !");
  MEDCouplingWithTimeStep *ret=new MEDCouplingWithTimeStep;
  DataArrayDouble *arr=DataArrayDouble::substract(getArray(),other->getArray());
  ret->setArray(arr,0);
  arr->decrRef();
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
  getArray()->substractEqual(other->getArray());
}

MEDCouplingTimeDiscretization *MEDCouplingWithTimeStep::multiply(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("WithTimeStep::multiply on mismatched time discretization !");
  MEDCouplingWithTimeStep *ret=new MEDCouplingWithTimeStep;
  DataArrayDouble *arr=DataArrayDouble::multiply(getArray(),other->getArray());
  ret->setArray(arr,0);
  arr->decrRef();
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
  getArray()->multiplyEqual(other->getArray());
}

MEDCouplingTimeDiscretization *MEDCouplingWithTimeStep::divide(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingWithTimeStep *otherC=dynamic_cast<const MEDCouplingWithTimeStep *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("WithTimeStep::divide on mismatched time discretization !");
  MEDCouplingWithTimeStep *ret=new MEDCouplingWithTimeStep;
  DataArrayDouble *arr=DataArrayDouble::divide(getArray(),other->getArray());
  ret->setArray(arr,0);
  arr->decrRef();
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
  getArray()->divideEqual(other->getArray());
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

std::vector< const DataArrayDouble *> MEDCouplingWithTimeStep::getArraysForTime(double time) const throw(INTERP_KERNEL::Exception)
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

MEDCouplingConstOnTimeInterval::MEDCouplingConstOnTimeInterval():_start_time(0.),_end_time(0.),_start_dt(-1),_end_dt(-1),_start_it(-1),_end_it(-1)
{
}

void MEDCouplingConstOnTimeInterval::getTinySerializationIntInformation(std::vector<int>& tinyInfo) const
{
  MEDCouplingTimeDiscretization::getTinySerializationIntInformation(tinyInfo);
  tinyInfo.push_back(_start_dt);
  tinyInfo.push_back(_start_it);
  tinyInfo.push_back(_end_dt);
  tinyInfo.push_back(_end_it);
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
  _start_dt=tinyInfoI[2];
  _start_it=tinyInfoI[3];
  _end_dt=tinyInfoI[4];
  _end_it=tinyInfoI[5];
}

MEDCouplingConstOnTimeInterval::MEDCouplingConstOnTimeInterval(const MEDCouplingConstOnTimeInterval& other, bool deepCpy):
  MEDCouplingTimeDiscretization(other,deepCpy),_start_time(other._start_time),_end_time(other._end_time),_start_dt(other._start_dt),
  _end_dt(other._end_dt),_start_it(other._start_it),_end_it(other._end_it)
{
}

MEDCouplingTimeDiscretization *MEDCouplingConstOnTimeInterval::performCpy(bool deepCpy) const
{
  return new MEDCouplingConstOnTimeInterval(*this,deepCpy);
}

std::vector< const DataArrayDouble *> MEDCouplingConstOnTimeInterval::getArraysForTime(double time) const throw(INTERP_KERNEL::Exception)
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

bool MEDCouplingConstOnTimeInterval::areStrictlyCompatible(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areStrictlyCompatible(other))
    return false;
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  return otherC!=0;
}

bool MEDCouplingConstOnTimeInterval::areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areStrictlyCompatible(other))
    return false;
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  return otherC!=0;
}

bool MEDCouplingConstOnTimeInterval::isEqual(const MEDCouplingTimeDiscretization *other, double prec) const
{
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  if(!otherC)
    return false;
  if(_start_dt!=otherC->_start_dt)
    return false;
  if(_start_it!=otherC->_start_it)
    return false;
  if(std::fabs(_start_time-otherC->_start_time)>_time_tolerance)
    return false;
  if(_end_dt!=otherC->_end_dt)
    return false;
  if(_end_it!=otherC->_end_it)
    return false;
  if(std::fabs(_end_time-otherC->_end_time)>_time_tolerance)
    return false;
  return MEDCouplingTimeDiscretization::isEqual(other,prec);
}

void MEDCouplingConstOnTimeInterval::getValueOnTime(int eltId, double time, double *value) const throw(INTERP_KERNEL::Exception)
{
  if(time>_start_time-_time_tolerance && time<_end_time+_time_tolerance)
    if(_array)
      _array->getTuple(eltId,value);
    else
      throw INTERP_KERNEL::Exception("No array existing.");
  else
    throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingConstOnTimeInterval::getValueOnDiscTime(int eltId, int dt, int it, double *value) const throw(INTERP_KERNEL::Exception)
{
  if(dt>=_start_dt && dt<=_end_dt)
    if(_array)
      _array->getTuple(eltId,value);
    else
      throw INTERP_KERNEL::Exception("No array existing.");
  else
    throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

void MEDCouplingConstOnTimeInterval::checkNoTimePresence() const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("No time specified on a field defined as constant on one time interval");
}

void MEDCouplingConstOnTimeInterval::checkTimePresence(double time) const throw(INTERP_KERNEL::Exception)
{
  if(time<_start_time-_time_tolerance || time>_end_time+_time_tolerance)
    {
      std::ostringstream stream;
      stream << "The field is defined between times " << _start_time << " and " << _end_time << " with tolerance ";
      stream << _time_tolerance << " and trying to access on time = " << time;
      throw INTERP_KERNEL::Exception(stream.str().c_str());
    }
}

MEDCouplingTimeDiscretization *MEDCouplingConstOnTimeInterval::aggregate(const MEDCouplingTimeDiscretization *other) const
{
   const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("ConstOnTimeInterval::aggregation on mismatched time discretization !");
  MEDCouplingConstOnTimeInterval *ret=new MEDCouplingConstOnTimeInterval;
  ret->setTimeTolerance(getTimeTolerance());
  DataArrayDouble *arr=DataArrayDouble::aggregate(getArray(),other->getArray());
  ret->setArray(arr,0);
  arr->decrRef();
  int tmp1,tmp2;
  double tmp3=getStartTime(tmp1,tmp2);
  ret->setStartTime(tmp3,tmp1,tmp2);
  tmp3=getEndTime(tmp1,tmp2);
  ret->setEndTime(tmp3,tmp1,tmp2);
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingConstOnTimeInterval::add(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("ConstOnTimeInterval::add on mismatched time discretization !");
  MEDCouplingConstOnTimeInterval *ret=new MEDCouplingConstOnTimeInterval;
  DataArrayDouble *arr=DataArrayDouble::add(getArray(),other->getArray());
  ret->setArray(arr,0);
  arr->decrRef();
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
  getArray()->addEqual(other->getArray());
}
 
MEDCouplingTimeDiscretization *MEDCouplingConstOnTimeInterval::substract(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("ConstOnTimeInterval::substract on mismatched time discretization !");
  MEDCouplingConstOnTimeInterval *ret=new MEDCouplingConstOnTimeInterval;
  DataArrayDouble *arr=DataArrayDouble::substract(getArray(),other->getArray());
  ret->setArray(arr,0);
  arr->decrRef();
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
  getArray()->substractEqual(other->getArray());
}

MEDCouplingTimeDiscretization *MEDCouplingConstOnTimeInterval::multiply(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("multiply on mismatched time discretization !");
  MEDCouplingConstOnTimeInterval *ret=new MEDCouplingConstOnTimeInterval;
  DataArrayDouble *arr=DataArrayDouble::multiply(getArray(),other->getArray());
  ret->setArray(arr,0);
  arr->decrRef();
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
  getArray()->multiplyEqual(other->getArray());
}

MEDCouplingTimeDiscretization *MEDCouplingConstOnTimeInterval::divide(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingConstOnTimeInterval *otherC=dynamic_cast<const MEDCouplingConstOnTimeInterval *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("divide on mismatched time discretization !");
  MEDCouplingConstOnTimeInterval *ret=new MEDCouplingConstOnTimeInterval;
  DataArrayDouble *arr=DataArrayDouble::divide(getArray(),other->getArray());
  ret->setArray(arr,0);
  arr->decrRef();
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
  getArray()->divideEqual(other->getArray());
}

MEDCouplingTwoTimeSteps::MEDCouplingTwoTimeSteps(const MEDCouplingTwoTimeSteps& other, bool deepCpy):MEDCouplingTimeDiscretization(other,deepCpy),
                                                                                                     _start_time(other._start_time),_end_time(other._end_time),
                                                                                                     _start_dt(other._start_dt),_end_dt(other._end_dt),
                                                                                                     _start_it(other._start_it),_end_it(other._end_it)
{
  if(other._end_array)
    _end_array=other._end_array->performCpy(deepCpy);
  else
    _end_array=0;
}

void MEDCouplingTwoTimeSteps::updateTime()
{
  MEDCouplingTimeDiscretization::updateTime();
  if(_end_array)
    updateTimeWith(*_end_array);
}

void MEDCouplingTwoTimeSteps::copyTinyAttrFrom(const MEDCouplingTimeDiscretization& other)
{
  MEDCouplingTimeDiscretization::copyTinyAttrFrom(other);
  const MEDCouplingTwoTimeSteps& otherC=dynamic_cast<const MEDCouplingTwoTimeSteps& >(other);
  _start_time=otherC._start_time;
  _end_time=otherC._end_time;
  _start_dt=otherC._start_dt;
  _end_dt=otherC._end_dt;
  _start_it=otherC._start_it;
  _end_it=otherC._end_it;
}

void MEDCouplingTwoTimeSteps::copyTinyStringsFrom(const MEDCouplingTimeDiscretization& other)
{
  MEDCouplingTimeDiscretization::copyTinyStringsFrom(other);
  const MEDCouplingTwoTimeSteps* otherC=dynamic_cast<const MEDCouplingTwoTimeSteps* >(&other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("Trying to operate copyTinyStringsFrom on different field type (two times//one time) !");
  if(_end_array && otherC->_end_array)
    _end_array->copyStringInfoFrom(*otherC->_end_array);
}

DataArrayDouble *MEDCouplingTwoTimeSteps::getEndArray() const
{
  return _end_array;
}

void MEDCouplingTwoTimeSteps::checkCoherency() const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingTimeDiscretization::checkCoherency();
  if(!_end_array)
    throw INTERP_KERNEL::Exception("No end array specified !");
  if(_array->getNumberOfComponents()!=_end_array->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("The number of components mismatch between the start and the end arrays !");
  if(_array->getNumberOfTuples()!=_end_array->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("The number of tuples mismatch between the start and the end arrays !");
}

bool MEDCouplingTwoTimeSteps::isEqual(const MEDCouplingTimeDiscretization *other, double prec) const
{
  const MEDCouplingTwoTimeSteps *otherC=dynamic_cast<const MEDCouplingTwoTimeSteps *>(other);
  if(!otherC)
    return false;
  if(_start_dt!=otherC->_start_dt)
    return false;
  if(_end_dt!=otherC->_end_dt)
    return false;
  if(_start_it!=otherC->_start_it)
    return false;
  if(_end_it!=otherC->_end_it)
    return false;
  if(std::fabs(_start_time-otherC->_start_time)>_time_tolerance)
    return false;
  if(std::fabs(_end_time-otherC->_end_time)>_time_tolerance)
    return false;
  if(_end_array!=otherC->_end_array)
    if(!_end_array->isEqual(*otherC->_end_array,prec))
      return false;
  return MEDCouplingTimeDiscretization::isEqual(other,prec);
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
  tinyInfo.push_back(_start_dt);
  tinyInfo.push_back(_start_it);
  tinyInfo.push_back(_end_dt);
  tinyInfo.push_back(_end_it);
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

void MEDCouplingTwoTimeSteps::finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS)
{
  MEDCouplingTimeDiscretization::finishUnserialization(tinyInfoI,tinyInfoD,tinyInfoS);
  _start_time=tinyInfoD[1];
  _end_time=tinyInfoD[2];
  _start_dt=tinyInfoI[2];
  _start_it=tinyInfoI[3];
  _end_dt=tinyInfoI[4];
  _end_it=tinyInfoI[5];
}

std::vector< const DataArrayDouble *> MEDCouplingTwoTimeSteps::getArraysForTime(double time) const throw(INTERP_KERNEL::Exception)
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

void MEDCouplingTwoTimeSteps::setArrays(const std::vector<DataArrayDouble *>& arrays, TimeLabel *owner) throw(INTERP_KERNEL::Exception)
{
  if(arrays.size()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingTwoTimeSteps::setArrays : number of arrays must be two.");
  setArray(arrays.front(),owner);
  setArray(arrays.back(),owner);
}

MEDCouplingLinearTime::MEDCouplingLinearTime(const MEDCouplingLinearTime& other, bool deepCpy):MEDCouplingTwoTimeSteps(other,deepCpy)
{
}

MEDCouplingLinearTime::MEDCouplingLinearTime()
{
}

void MEDCouplingLinearTime::checkCoherency() const throw(INTERP_KERNEL::Exception)
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
  return otherC!=0;
}

bool MEDCouplingLinearTime::areStrictlyCompatible(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areStrictlyCompatible(other))
    return false;
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  return otherC!=0;
}

bool MEDCouplingLinearTime::areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretization *other) const
{
  if(!MEDCouplingTimeDiscretization::areStrictlyCompatibleForMul(other))
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
  int nbComp=vals.size()/2;
  std::transform(vals.begin(),vals.begin()+nbComp,res,std::bind2nd(std::multiplies<double>(),alpha));
  std::vector<double> tmp(nbComp);
  std::transform(vals.begin()+nbComp,vals.end(),tmp.begin(),std::bind2nd(std::multiplies<double>(),1-alpha));
  std::transform(tmp.begin(),tmp.end(),res,res,std::plus<double>());
}

void MEDCouplingLinearTime::getValueOnTime(int eltId, double time, double *value) const throw(INTERP_KERNEL::Exception)
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

void MEDCouplingLinearTime::getValueOnDiscTime(int eltId, int dt, int it, double *value) const throw(INTERP_KERNEL::Exception)
{
  if(dt==_start_dt && it==_start_it)
    {
      if(_array)
        _array->getTuple(eltId,value);
      else
        throw INTERP_KERNEL::Exception("dt it match with start time but no start array existing.");
    }
  if(dt==_end_dt && it==_end_it)
    {
      if(_end_array)
        _end_array->getTuple(eltId,value);
      else
        throw INTERP_KERNEL::Exception("dt it match with end time but no end array existing.");
    }
  else
    throw INTERP_KERNEL::Exception(EXCEPTION_MSG);
}

MEDCouplingTimeDiscretization *MEDCouplingLinearTime::aggregate(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("LinearTime::aggregation on mismatched time discretization !");
  MEDCouplingLinearTime *ret=new MEDCouplingLinearTime;
  ret->setTimeTolerance(getTimeTolerance());
  DataArrayDouble *arr1=DataArrayDouble::aggregate(getArray(),other->getArray());
  ret->setArray(arr1,0);
  arr1->decrRef();
  DataArrayDouble *arr2=DataArrayDouble::aggregate(getEndArray(),other->getEndArray());
  ret->setEndArray(arr2,0);
  arr2->decrRef();
  return ret;
}

MEDCouplingTimeDiscretization *MEDCouplingLinearTime::add(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("LinearTime::add on mismatched time discretization !");
  MEDCouplingLinearTime *ret=new MEDCouplingLinearTime;
  DataArrayDouble *arr1=DataArrayDouble::add(getArray(),other->getArray());
  ret->setArray(arr1,0);
  arr1->decrRef();
  DataArrayDouble *arr2=DataArrayDouble::add(getEndArray(),other->getEndArray());
  ret->setEndArray(arr2,0);
  arr2->decrRef();
  return ret;
}

void MEDCouplingLinearTime::addEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("LinearTime::addEqual on mismatched time discretization !");
  getArray()->addEqual(other->getArray());
  getEndArray()->addEqual(other->getEndArray());
}

MEDCouplingTimeDiscretization *MEDCouplingLinearTime::substract(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("LinearTime::substract on mismatched time discretization !");
  MEDCouplingLinearTime *ret=new MEDCouplingLinearTime;
  DataArrayDouble *arr1=DataArrayDouble::substract(getArray(),other->getArray());
  ret->setArray(arr1,0);
  arr1->decrRef();
  DataArrayDouble *arr2=DataArrayDouble::substract(getEndArray(),other->getEndArray());
  ret->setEndArray(arr2,0);
  arr2->decrRef();
  return ret;
}

void MEDCouplingLinearTime::substractEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("LinearTime::addEqual on mismatched time discretization !");
  getArray()->substractEqual(other->getArray());
  getEndArray()->substractEqual(other->getEndArray());
}

MEDCouplingTimeDiscretization *MEDCouplingLinearTime::multiply(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("LinearTime::multiply on mismatched time discretization !");
  MEDCouplingLinearTime *ret=new MEDCouplingLinearTime;
  DataArrayDouble *arr1=DataArrayDouble::multiply(getArray(),other->getArray());
  ret->setArray(arr1,0);
  arr1->decrRef();
  DataArrayDouble *arr2=DataArrayDouble::multiply(getEndArray(),other->getEndArray());
  ret->setEndArray(arr2,0);
  arr2->decrRef();
  return ret;
}

void MEDCouplingLinearTime::multiplyEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("LinearTime::addEqual on mismatched time discretization !");
  getArray()->multiplyEqual(other->getArray());
  getEndArray()->multiplyEqual(other->getEndArray());
}

MEDCouplingTimeDiscretization *MEDCouplingLinearTime::divide(const MEDCouplingTimeDiscretization *other) const
{
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("LinearTime::divide on mismatched time discretization !");
  MEDCouplingLinearTime *ret=new MEDCouplingLinearTime;
  DataArrayDouble *arr1=DataArrayDouble::divide(getArray(),other->getArray());
  ret->setArray(arr1,0);
  arr1->decrRef();
  DataArrayDouble *arr2=DataArrayDouble::divide(getEndArray(),other->getEndArray());
  ret->setEndArray(arr2,0);
  arr2->decrRef();
  return ret;
}

void MEDCouplingLinearTime::divideEqual(const MEDCouplingTimeDiscretization *other)
{
  const MEDCouplingLinearTime *otherC=dynamic_cast<const MEDCouplingLinearTime *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("LinearTime::addEqual on mismatched time discretization !");
  getArray()->divideEqual(other->getArray());
  getEndArray()->divideEqual(other->getEndArray());
}
