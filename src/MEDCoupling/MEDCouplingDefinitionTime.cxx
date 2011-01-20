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

#include "MEDCouplingDefinitionTime.hxx"
#include "MEDCouplingFieldDouble.hxx"

#include <cmath>

using namespace ParaMEDMEM;

MEDCouplingDefinitionTimeSlice *MEDCouplingDefinitionTimeSlice::New(const MEDCouplingFieldDouble *f, int meshId, const std::vector<int>& arrId, int fieldId) throw(INTERP_KERNEL::Exception)
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

MEDCouplingDefinitionTimeSlice::MEDCouplingDefinitionTimeSlice(const MEDCouplingFieldDouble *f, int meshId, int arrId, int fieldId) throw(INTERP_KERNEL::Exception):_mesh_id(meshId),_array_id(arrId),_field_id(fieldId)
{
  int tmp1,tmp2;
  double t1=f->getStartTime(tmp1,tmp2);
  double t2=f->getEndTime(tmp1,tmp2);
  if(t2<t1)
    throw INTERP_KERNEL::Exception("MEDCouplingDefinitionTimeSlice : End time strictly before Start time ...");
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

MEDCouplingDefinitionTimeSliceInst::MEDCouplingDefinitionTimeSliceInst(const MEDCouplingFieldDouble *f, int meshId, int arrId, int fieldId) throw(INTERP_KERNEL::Exception):MEDCouplingDefinitionTimeSlice(f,meshId,arrId,fieldId)
{
  int tmp1,tmp2;
  double t1=f->getStartTime(tmp1,tmp2);
  double t2=f->getEndTime(tmp1,tmp2);
  double eps=f->getTimeTolerance();
  if(fabs(t1-t2)>eps)
    throw INTERP_KERNEL::Exception("MEDCouplingDefinitionTimeSliceInst : times differs in this");
  _instant=t1;
}

bool MEDCouplingDefinitionTimeSliceCstOnTI::isContaining(double tmp, double eps) const
{
  return _start-eps>tmp && _end+eps<tmp;
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

MEDCouplingDefinitionTimeSliceCstOnTI::MEDCouplingDefinitionTimeSliceCstOnTI(const MEDCouplingFieldDouble *f, int meshId, int arrId, int fieldId) throw(INTERP_KERNEL::Exception):MEDCouplingDefinitionTimeSlice(f,meshId,arrId,fieldId)
{
  int tmp1,tmp2;
  double t1=f->getStartTime(tmp1,tmp2);
  double t2=f->getEndTime(tmp1,tmp2);
  _start=t1;
  _end=t2;
}

bool MEDCouplingDefinitionTimeSliceLT::isContaining(double tmp, double eps) const
{
  return _start-eps>tmp && _end+eps<tmp;
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

MEDCouplingDefinitionTimeSliceLT::MEDCouplingDefinitionTimeSliceLT(const MEDCouplingFieldDouble *f, int meshId, int arrId, int arr2Id, int fieldId) throw(INTERP_KERNEL::Exception):MEDCouplingDefinitionTimeSlice(f,meshId,arrId,fieldId),_array_id_end(arr2Id)
{
  int tmp1,tmp2;
  double t1=f->getStartTime(tmp1,tmp2);
  double t2=f->getEndTime(tmp1,tmp2);
  _start=t1;
  _end=t2;
}

MEDCouplingDefinitionTime::MEDCouplingDefinitionTime()
{
}

MEDCouplingDefinitionTime::MEDCouplingDefinitionTime(const std::vector<const MEDCouplingFieldDouble *>& fs, const std::vector<int>& meshRefs, const std::vector<std::vector<int> >& arrRefs) throw(INTERP_KERNEL::Exception)
{
  std::size_t sz=fs.size();
  if(sz!=arrRefs.size())
    throw INTERP_KERNEL::Exception("MEDCouplingDefinitionTime constructor : internal error ! should never happen !");
  _slices.resize(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      if(arrRefs.empty())
        throw INTERP_KERNEL::Exception("MEDCouplingDefinitionTime constructor : A field is null in list impossible to build a time definition !");
      _slices[i]=MEDCouplingDefinitionTimeSlice::New(fs[i],meshRefs[i],arrRefs[i],i);
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

void MEDCouplingDefinitionTime::getIdsOnTimeRight(double tm, int& meshId, int& arrId, int& arrIdInField, int& fieldId) const throw(INTERP_KERNEL::Exception)
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

void MEDCouplingDefinitionTime::getIdsOnTimeLeft(double tm, int& meshId, int& arrId, int& arrIdInField, int& fieldId) const throw(INTERP_KERNEL::Exception)
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

void MEDCouplingDefinitionTime::getIdsOnTime(double tm, std::vector<int>& meshIds, std::vector<int>& arrIds, std::vector<int>& arrIdsInField, std::vector<int>& fieldIds) const throw(INTERP_KERNEL::Exception)
{
  std::vector<int> ids;
  int id=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingDefinitionTimeSlice> >::const_iterator it=_slices.begin();it!=_slices.end();it++,id++)
    if((*it)->isContaining(tm,_eps))
      ids.push_back(id);
  if(ids.empty())
    throw INTERP_KERNEL::Exception("MEDCouplingDefinitionTime::getIdsOnTime : No matching slice for such time !");
  int sz=ids.size();
  if(sz>2)
    throw INTERP_KERNEL::Exception("MEDCouplingDefinitionTime::getIdsOnTime : Too many slices match this time !");
  //tony
}

void MEDCouplingDefinitionTime::appendRepr(std::ostream& stream) const
{
  stream << "Time definition :\n";
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingDefinitionTimeSlice> >::const_iterator it=_slices.begin();it!=_slices.end();it++)
    {
      stream << " - ";
      (*it)->appendRepr(stream);
      stream << std::endl;
    }
}
