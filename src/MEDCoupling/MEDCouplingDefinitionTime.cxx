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

#include "MEDCouplingDefinitionTime.hxx"
#include "MEDCouplingFieldDouble.hxx"

#include <cmath>

using namespace MEDCoupling;

const double MEDCouplingDefinitionTime::EPS_DFT=1e-15;

MEDCouplingDefinitionTimeSlice *MEDCouplingDefinitionTimeSlice::New(const MEDCouplingFieldDouble *f, int meshId, const std::vector<int>& arrId, int fieldId)
{
  static const char msg[]="TimeSlice::New : mismatch of arrays number of a fieldDouble and its policy !!! Internal error !!!";
  if(!f)
    throw INTERP_KERNEL::Exception("MEDCouplingDefinitionTimeSlice::New : empty field !");
  switch(f->getTimeDiscretization())
  {
    case ONE_TIME:
      {
        if(arrId.size()!=1)
          throw INTERP_KERNEL::Exception(msg);
        return new MEDCouplingDefinitionTimeSliceInst(f,meshId,arrId[0],fieldId);
      }
    case CONST_ON_TIME_INTERVAL:
      {
        if(arrId.size()!=1)
          throw INTERP_KERNEL::Exception(msg);
        return new MEDCouplingDefinitionTimeSliceCstOnTI(f,meshId,arrId[0],fieldId);
      }
    case LINEAR_TIME:
      {
        if(arrId.size()!=2)
          throw INTERP_KERNEL::Exception(msg);
        return new MEDCouplingDefinitionTimeSliceLT(f,meshId,arrId[0],arrId[1],fieldId);
      }
    case NO_TIME:
      throw INTERP_KERNEL::Exception("Invalide time discretization ! NO_TIME ! Impossible to build a definition time slice !");
    default:
      throw INTERP_KERNEL::Exception("Invalide time discretization : Not recognized !");
  }
}

MEDCouplingDefinitionTimeSlice *MEDCouplingDefinitionTimeSlice::New(TypeOfTimeDiscretization type, const std::vector<int>& tiI, const std::vector<double>& tiD)
{
  switch(type)
  {
    case ONE_TIME:
      return MEDCouplingDefinitionTimeSliceInst::New(tiI,tiD);
    case CONST_ON_TIME_INTERVAL:
      return MEDCouplingDefinitionTimeSliceCstOnTI::New(tiI,tiD);
    case LINEAR_TIME:
      return MEDCouplingDefinitionTimeSliceLT::New(tiI,tiD);
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingDefinitionTimeSlice::New : unrecognized time discretization type !");
  }
}

bool MEDCouplingDefinitionTimeSlice::isEqual(const MEDCouplingDefinitionTimeSlice& other, double eps) const
{
  if(_mesh_id!=other._mesh_id)
    return false;
  if(_array_id!=other._array_id)
    return false;
  if(_field_id!=other._field_id)
    return false;
  return true;
}

int MEDCouplingDefinitionTimeSlice::getStartId() const
{
  return _array_id;
}

int MEDCouplingDefinitionTimeSlice::getEndId() const
{
  return _array_id;
}

void MEDCouplingDefinitionTimeSlice::appendRepr(std::ostream& stream) const
{
  stream << " *** MeshId : " << _mesh_id << " ArrayId : " << _array_id;
}

MEDCouplingDefinitionTimeSlice::MEDCouplingDefinitionTimeSlice(const MEDCouplingFieldDouble *f, int meshId, int arrId, int fieldId):_mesh_id(meshId),_array_id(arrId),_field_id(fieldId)
{
  int tmp1,tmp2;
  double t1=f->getStartTime(tmp1,tmp2);
  double t2=f->getEndTime(tmp1,tmp2);
  if(t2<t1)
    throw INTERP_KERNEL::Exception("MEDCouplingDefinitionTimeSlice : End time strictly before Start time ...");
}

std::size_t MEDCouplingDefinitionTimeSlice::getHeapMemorySizeWithoutChildren() const
{
  return 0;
}

std::vector<const BigMemoryObject *> MEDCouplingDefinitionTimeSlice::getDirectChildrenWithNull() const
{
  return std::vector<const BigMemoryObject *>();
}

bool MEDCouplingDefinitionTimeSlice::isFullyIncludedInMe(const MEDCouplingDefinitionTimeSlice *other, double eps) const
{
  double t1=getStartTime();
  double t2=getEndTime();
  double o1=other->getStartTime();
  double o2=other->getEndTime();
  return o1>t1-eps && o2<t2+eps;
}

bool MEDCouplingDefinitionTimeSlice::isOverllapingWithMe(const MEDCouplingDefinitionTimeSlice *other, double eps) const
{
  double t1=getStartTime();
  double t2=getEndTime();
  double o1=other->getStartTime();
  double o2=other->getEndTime();
  return (o1<t1+eps && o2<t1+eps) || (o1>t2-eps && o2>t2-eps);
}

bool MEDCouplingDefinitionTimeSlice::isAfterMe(const MEDCouplingDefinitionTimeSlice *other, double eps) const
{
  double t2=getEndTime();
  double o1=other->getStartTime();
  double o2=other->getEndTime();
  return (o1>t2-eps && o2>t2-eps);
}

bool MEDCouplingDefinitionTimeSlice::isBeforeMe(const MEDCouplingDefinitionTimeSlice *other, double eps) const
{
  double t1=getStartTime();
  double o1=other->getStartTime();
  double o2=other->getEndTime();
  return (o1<t1+eps && o2<t1+eps);
}

MEDCouplingDefinitionTimeSliceInst *MEDCouplingDefinitionTimeSliceInst::New(const std::vector<int>& tiI, const std::vector<double>& tiD)
{
  MEDCouplingDefinitionTimeSliceInst *ret=new MEDCouplingDefinitionTimeSliceInst;
  ret->unserialize(tiI,tiD);
  return ret;
}

void MEDCouplingDefinitionTimeSliceInst::getTinySerializationInformation(std::vector<int>& tiI, std::vector<double>& tiD) const
{
  tiI.resize(3);
  tiI[0]=_mesh_id; tiI[1]=_array_id; tiI[2]=_field_id;
  tiD.resize(1);
  tiD[0]=_instant;
}

void MEDCouplingDefinitionTimeSliceInst::unserialize(const std::vector<int>& tiI, const std::vector<double>& tiD)
{
  _mesh_id=tiI[0]; _array_id=tiI[1]; _field_id=tiI[2];
  _instant=tiD[0];
}

TypeOfTimeDiscretization MEDCouplingDefinitionTimeSliceInst::getTimeType() const
{
  return ONE_TIME;
}

MEDCouplingDefinitionTimeSlice *MEDCouplingDefinitionTimeSliceInst::copy() const
{
  return new MEDCouplingDefinitionTimeSliceInst(*this);
}

bool MEDCouplingDefinitionTimeSliceInst::isEqual(const MEDCouplingDefinitionTimeSlice& other, double eps) const
{
  if(!MEDCouplingDefinitionTimeSlice::isEqual(other,eps))
    return false;
  const MEDCouplingDefinitionTimeSliceInst *otherC=dynamic_cast<const MEDCouplingDefinitionTimeSliceInst *>(&other);
  if(!otherC)
    return false;
  return fabs(otherC->_instant-_instant)<eps;
}

void MEDCouplingDefinitionTimeSliceInst::getHotSpotsTime(std::vector<double>& ret) const
{
  ret.resize(1);
  ret[0]=_instant;
}

void MEDCouplingDefinitionTimeSliceInst::getIdsOnTime(double tm, double eps, int& meshId, int& arrId, int& arrIdInField, int& fieldId) const
{
  meshId=_mesh_id;
  arrId=_array_id;
  arrIdInField=0;
  fieldId=_field_id;
}

bool MEDCouplingDefinitionTimeSliceInst::isContaining(double tmp, double eps) const
{
  return fabs(tmp-_instant)<eps;
}

void MEDCouplingDefinitionTimeSliceInst::appendRepr(std::ostream& stream) const
{
  stream << "single point " << _instant;
  MEDCouplingDefinitionTimeSlice::appendRepr(stream);
}

double MEDCouplingDefinitionTimeSliceInst::getStartTime() const
{
  return _instant;
}

double MEDCouplingDefinitionTimeSliceInst::getEndTime() const
{
  return _instant;
}

MEDCouplingDefinitionTimeSliceInst::MEDCouplingDefinitionTimeSliceInst(const MEDCouplingFieldDouble *f, int meshId, int arrId, int fieldId):MEDCouplingDefinitionTimeSlice(f,meshId,arrId,fieldId)
{
  int tmp1,tmp2;
  double t1=f->getStartTime(tmp1,tmp2);
  double t2=f->getEndTime(tmp1,tmp2);
  double eps=f->getTimeTolerance();
  if(fabs(t1-t2)>eps)
    throw INTERP_KERNEL::Exception("MEDCouplingDefinitionTimeSliceInst : times differs in this");
  _instant=t1;
}

MEDCouplingDefinitionTimeSliceCstOnTI *MEDCouplingDefinitionTimeSliceCstOnTI::New(const std::vector<int>& tiI, const std::vector<double>& tiD)
{
  MEDCouplingDefinitionTimeSliceCstOnTI *ret=new MEDCouplingDefinitionTimeSliceCstOnTI;
  ret->unserialize(tiI,tiD);
  return ret;
}

MEDCouplingDefinitionTimeSlice *MEDCouplingDefinitionTimeSliceCstOnTI::copy() const
{
  return new MEDCouplingDefinitionTimeSliceCstOnTI(*this);
}

bool MEDCouplingDefinitionTimeSliceCstOnTI::isEqual(const MEDCouplingDefinitionTimeSlice& other, double eps) const
{
  if(!MEDCouplingDefinitionTimeSlice::isEqual(other,eps))
    return false;
  const MEDCouplingDefinitionTimeSliceCstOnTI *otherC=dynamic_cast<const MEDCouplingDefinitionTimeSliceCstOnTI *>(&other);
  if(!otherC)
    return false;
  if(fabs(otherC->_start-_start)>eps)
    return false;
  return fabs(otherC->_end-_end)<eps;
}

void MEDCouplingDefinitionTimeSliceCstOnTI::getHotSpotsTime(std::vector<double>& ret) const
{
  ret.resize(1);
  ret[0]=(_start+_end)/2.;
}

void MEDCouplingDefinitionTimeSliceCstOnTI::getIdsOnTime(double tm, double eps, int& meshId, int& arrId, int& arrIdInField, int& fieldId) const
{
  meshId=_mesh_id;
  arrId=_array_id;
  arrIdInField=0;
  fieldId=_field_id;
}

bool MEDCouplingDefinitionTimeSliceCstOnTI::isContaining(double tmp, double eps) const
{
  return _start-eps<tmp && _end+eps>tmp;
}

void MEDCouplingDefinitionTimeSliceCstOnTI::appendRepr(std::ostream& stream) const
{
  stream << "Constant on time interval [" << _start << "," << _end << "]";
  MEDCouplingDefinitionTimeSlice::appendRepr(stream);
}

double MEDCouplingDefinitionTimeSliceCstOnTI::getStartTime() const
{
  return _start;
}

double MEDCouplingDefinitionTimeSliceCstOnTI::getEndTime() const
{
  return _end;
}

void MEDCouplingDefinitionTimeSliceCstOnTI::getTinySerializationInformation(std::vector<int>& tiI, std::vector<double>& tiD) const
{
  tiI.resize(3);
  tiI[0]=_mesh_id; tiI[1]=_array_id; tiI[2]=_field_id;
  tiD.resize(2);
  tiD[0]=_start; tiD[1]=_end;
}

void MEDCouplingDefinitionTimeSliceCstOnTI::unserialize(const std::vector<int>& tiI, const std::vector<double>& tiD)
{
  _mesh_id=tiI[0]; _array_id=tiI[1]; _field_id=tiI[2];
  _start=tiD[0]; _end=tiD[1];
}

TypeOfTimeDiscretization MEDCouplingDefinitionTimeSliceCstOnTI::getTimeType() const
{
  return CONST_ON_TIME_INTERVAL;
}

MEDCouplingDefinitionTimeSliceCstOnTI::MEDCouplingDefinitionTimeSliceCstOnTI(const MEDCouplingFieldDouble *f, int meshId, int arrId, int fieldId):MEDCouplingDefinitionTimeSlice(f,meshId,arrId,fieldId)
{
  int tmp1,tmp2;
  double t1=f->getStartTime(tmp1,tmp2);
  double t2=f->getEndTime(tmp1,tmp2);
  _start=t1;
  _end=t2;
}

MEDCouplingDefinitionTimeSliceLT *MEDCouplingDefinitionTimeSliceLT::New(const std::vector<int>& tiI, const std::vector<double>& tiD)
{
  MEDCouplingDefinitionTimeSliceLT *ret=new MEDCouplingDefinitionTimeSliceLT;
  ret->unserialize(tiI,tiD);
  return ret;
}

MEDCouplingDefinitionTimeSlice *MEDCouplingDefinitionTimeSliceLT::copy() const
{
  return new MEDCouplingDefinitionTimeSliceLT(*this);
}

bool MEDCouplingDefinitionTimeSliceLT::isEqual(const MEDCouplingDefinitionTimeSlice& other, double eps) const
{
  if(!MEDCouplingDefinitionTimeSlice::isEqual(other,eps))
    return false;
  const MEDCouplingDefinitionTimeSliceLT *otherC=dynamic_cast<const MEDCouplingDefinitionTimeSliceLT *>(&other);
  if(!otherC)
    return false;
  if(_array_id_end!=otherC->_array_id_end)
    return false;
  if(fabs(otherC->_start-_start)>eps)
    return false;
  return fabs(otherC->_end-_end)<eps;
}

void MEDCouplingDefinitionTimeSliceLT::getHotSpotsTime(std::vector<double>& ret) const
{
  ret.resize(2);
  ret[0]=_start;
  ret[1]=_end;
}

void MEDCouplingDefinitionTimeSliceLT::getIdsOnTime(double tm, double eps, int& meshId, int& arrId, int& arrIdInField, int& fieldId) const
{
  if(fabs(tm-_start)<eps)
    {
      meshId=_mesh_id;
      arrId=_array_id;
      arrIdInField=0;
      fieldId=_field_id;
      return ;
    }
  if(fabs(tm-_end)<eps)
    {
      meshId=_mesh_id;
      arrId=_array_id_end;
      arrIdInField=1;
      fieldId=_field_id;
      return ;
    }
  throw INTERP_KERNEL::Exception("LinearTime request not in boundary of this ! use hot spots !");
}

bool MEDCouplingDefinitionTimeSliceLT::isContaining(double tmp, double eps) const
{
  return _start-eps<tmp && _end+eps>tmp;
}

void MEDCouplingDefinitionTimeSliceLT::appendRepr(std::ostream& stream) const
{
  stream << "Linear on time interval [" << _start << "," << _end << "]";
  MEDCouplingDefinitionTimeSlice::appendRepr(stream);
  stream << " EndArrayId : " << _array_id_end;
}

double MEDCouplingDefinitionTimeSliceLT::getStartTime() const
{
  return _start;
}

double MEDCouplingDefinitionTimeSliceLT::getEndTime() const
{
  return _end;
}

int MEDCouplingDefinitionTimeSliceLT::getEndId() const
{
  return _array_id_end;
}

void MEDCouplingDefinitionTimeSliceLT::getTinySerializationInformation(std::vector<int>& tiI, std::vector<double>& tiD) const
{
  tiI.resize(4);
  tiI[0]=_mesh_id; tiI[1]=_array_id; tiI[2]=_field_id; tiI[3]=_array_id_end;
  tiD.resize(2);
  tiD[0]=_start; tiD[1]=_end;
}

void MEDCouplingDefinitionTimeSliceLT::unserialize(const std::vector<int>& tiI, const std::vector<double>& tiD)
{
  _mesh_id=tiI[0]; _array_id=tiI[1]; _field_id=tiI[2]; _array_id_end=tiI[3];
  _start=tiD[0]; _end=tiD[1];
}

TypeOfTimeDiscretization MEDCouplingDefinitionTimeSliceLT::getTimeType() const
{
  return LINEAR_TIME;
}

MEDCouplingDefinitionTimeSliceLT::MEDCouplingDefinitionTimeSliceLT(const MEDCouplingFieldDouble *f, int meshId, int arrId, int arr2Id, int fieldId):MEDCouplingDefinitionTimeSlice(f,meshId,arrId,fieldId),_array_id_end(arr2Id)
{
  int tmp1,tmp2;
  double t1=f->getStartTime(tmp1,tmp2);
  double t2=f->getEndTime(tmp1,tmp2);
  _start=t1;
  _end=t2;
}

MEDCouplingDefinitionTime::MEDCouplingDefinitionTime():_eps(EPS_DFT)
{
}

MEDCouplingDefinitionTime::MEDCouplingDefinitionTime(const std::vector<const MEDCouplingFieldDouble *>& fs, const std::vector<int>& meshRefs, const std::vector<std::vector<int> >& arrRefs)
{
  std::size_t sz=fs.size();
  if(sz!=arrRefs.size())
    throw INTERP_KERNEL::Exception("MEDCouplingDefinitionTime constructor : internal error ! should never happen !");
  _slices.resize(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      if(arrRefs.empty())
        throw INTERP_KERNEL::Exception("MEDCouplingDefinitionTime constructor : A field is null in list impossible to build a time definition !");
      _slices[i]=MEDCouplingDefinitionTimeSlice::New(fs[i],meshRefs[i],arrRefs[i],(int)i);
    }
  if(sz<=1)
    return ;
  const MEDCouplingDefinitionTimeSlice *ref=_slices[0];
  _eps=fs[0]->getTimeTolerance();
  for(std::size_t i=1;i<sz;i++)
    {
      if(!ref->isAfterMe(_slices[i],_eps))
        throw INTERP_KERNEL::Exception("MEDCouplingDefinitionTime constructors : the sequences of fields does NOT defines a stricly ascendant monotonic time sequence !");
      // double t1=ref->getEndTime();
      // double t2=_slices[i]->getStartTime();
      // if(fabs(t1-t2)<_eps)
      //   if(ref->getEndId() != _slices[i]->getStartId())
      //     throw INTERP_KERNEL::Exception("MEDCouplingDefinitionTime constructor : 2 slices refers to the same time and underlying arrays differs !");
      ref=_slices[i];
    }
}

std::size_t MEDCouplingDefinitionTime::getHeapMemorySizeWithoutChildren() const
{
  return _slices.capacity()*(sizeof(MEDCouplingDefinitionTimeSlice)+sizeof(int));
}

std::vector<const BigMemoryObject *> MEDCouplingDefinitionTime::getDirectChildrenWithNull() const
{
  return std::vector<const BigMemoryObject *>();
}

void MEDCouplingDefinitionTime::assign(const MEDCouplingDefinitionTime& other)
{
  std::size_t sz=other._slices.size();
  _slices.resize(sz);
  for(std::size_t i=0;i<sz;i++)
    _slices[i]=other._slices[i]->copy();
}

bool MEDCouplingDefinitionTime::isEqual(const MEDCouplingDefinitionTime& other) const
{
  std::size_t sz=_slices.size();
  if(sz!=other._slices.size())
    return false;
  for(std::size_t i=0;i<sz;i++)
    if(!_slices[i]->isEqual(*other._slices[i],_eps))
      return false;
  return true;
}

void MEDCouplingDefinitionTime::getIdsOnTimeRight(double tm, int& meshId, int& arrId, int& arrIdInField, int& fieldId) const
{
  std::vector<int> meshIds;
  std::vector<int> arrIds;
  std::vector<int> arrIdsInField;
  std::vector<int> fieldIds;
  getIdsOnTime(tm,meshIds,arrIds,arrIdsInField,fieldIds);
  meshId=meshIds.back();
  arrId=arrIds.back();
  arrIdInField=arrIdsInField.back();
  fieldId=fieldIds.back();
}

void MEDCouplingDefinitionTime::getIdsOnTimeLeft(double tm, int& meshId, int& arrId, int& arrIdInField, int& fieldId) const
{
  std::vector<int> meshIds;
  std::vector<int> arrIds;
  std::vector<int> arrIdsInField;
  std::vector<int> fieldIds;
  getIdsOnTime(tm,meshIds,arrIds,arrIdsInField,fieldIds);
  meshId=meshIds.front();
  arrId=arrIds.front();
  arrIdInField=arrIdsInField.front();
  fieldId=fieldIds.front();
}

void MEDCouplingDefinitionTime::getIdsOnTime(double tm, std::vector<int>& meshIds, std::vector<int>& arrIds, std::vector<int>& arrIdsInField, std::vector<int>& fieldIds) const
{
  std::vector<int> ids;
  int id=0;
  for(std::vector< MCAuto<MEDCouplingDefinitionTimeSlice> >::const_iterator it=_slices.begin();it!=_slices.end();it++,id++)
    if((*it)->isContaining(tm,_eps))
      ids.push_back(id);
  if(ids.empty())
    throw INTERP_KERNEL::Exception("MEDCouplingDefinitionTime::getIdsOnTime : No matching slice for such time !");
  std::size_t sz=ids.size();
  if(sz>2)
    throw INTERP_KERNEL::Exception("MEDCouplingDefinitionTime::getIdsOnTime : Too many slices match this time !");
  //
  meshIds.resize(sz);
  arrIds.resize(sz);
  arrIdsInField.resize(sz);
  fieldIds.resize(sz);
  for(std::size_t i=0;i<sz;i++)
    _slices[ids[i]]->getIdsOnTime(tm,_eps,meshIds[i],arrIds[i],arrIdsInField[i],fieldIds[i]);
}

std::vector<double> MEDCouplingDefinitionTime::getHotSpotsTime() const
{
  std::vector<double> ret;
  for(std::vector< MCAuto<MEDCouplingDefinitionTimeSlice> >::const_iterator it=_slices.begin();it!=_slices.end();it++)
    {
      std::vector<double> tmp;
      (*it)->getHotSpotsTime(tmp);
      if(!ret.empty())
        {
          if(fabs(ret.back()-tmp.front())>_eps)
            ret.insert(ret.end(),tmp.begin(),tmp.end());
          else
            ret.insert(ret.end(),tmp.begin()+1,tmp.end());
        }
      else
        ret.insert(ret.end(),tmp.begin(),tmp.end());
    }
  return ret;
}

void MEDCouplingDefinitionTime::appendRepr(std::ostream& stream) const
{
  stream << "Time definition :\n";
  for(std::vector< MCAuto<MEDCouplingDefinitionTimeSlice> >::const_iterator it=_slices.begin();it!=_slices.end();it++)
    {
      stream << " - ";
      (*it)->appendRepr(stream);
      stream << std::endl;
    }
}

void MEDCouplingDefinitionTime::getTinySerializationInformation(std::vector<int>& tinyInfoI, std::vector<double>& tinyInfoD) const
{
  int sz=(int)_slices.size();
  tinyInfoD.resize(1);
  tinyInfoD[0]=_eps;
  tinyInfoI.resize(3*sz+2);
  tinyInfoI[0]=sz;
  std::vector<int> coreData;
  for(int i=0;i<sz;i++)
    {
      std::vector<int> tmp1;
      std::vector<double> tmp2;
      tinyInfoI[i+2]=(int)_slices[i]->getTimeType();
      _slices[i]->getTinySerializationInformation(tmp1,tmp2);
      tinyInfoI[i+sz+2]=(int)tmp1.size();
      tinyInfoI[i+2*sz+2]=(int)tmp2.size();
      coreData.insert(coreData.end(),tmp1.begin(),tmp1.end());
      tinyInfoD.insert(tinyInfoD.end(),tmp2.begin(),tmp2.end());
    }
  tinyInfoI[1]=(int)coreData.size();
  tinyInfoI.insert(tinyInfoI.end(),coreData.begin(),coreData.end());
}

void MEDCouplingDefinitionTime::unserialize(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD)
{
  int sz=tinyInfoI[0];
  _slices.resize(sz);
  _eps=tinyInfoD[0];
  int offset1=0;
  int offset2=1;
  for(int i=0;i<sz;i++)
    {
      TypeOfTimeDiscretization ty=(TypeOfTimeDiscretization) tinyInfoI[i+2];  
      int sz1=tinyInfoI[i+sz+2];
      int sz2=tinyInfoI[i+2*sz+2];
      std::vector<int> tmp1(tinyInfoI.begin()+3*sz+2+offset1,tinyInfoI.begin()+3*sz+2+offset1+sz1);
      std::vector<double> tmp2(tinyInfoD.begin()+offset2,tinyInfoD.begin()+offset2+sz2);
      MEDCouplingDefinitionTimeSlice *pt=MEDCouplingDefinitionTimeSlice::New(ty,tmp1,tmp2);
      _slices[i]=pt;
      offset1+=sz1;
      offset2+=sz2;
    }
}

