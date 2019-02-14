// Copyright (C) 2017-2019  CEA/DEN, EDF R&D
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

#include "MEDFileFieldMultiTS.hxx"
#include "MEDFileFieldVisitor.hxx"
#include "MEDFileSafeCaller.txx"
#include "MEDLoaderBase.hxx"
#include "MEDFileField.txx"

#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldTemplate.hxx"

#include <sstream>

using namespace MEDCoupling;

template class MEDCoupling::MEDFileTemplateFieldMultiTSWithoutSDA<int>;
template class MEDCoupling::MEDFileTemplateFieldMultiTSWithoutSDA<float>;
template class MEDCoupling::MEDFileTemplateFieldMultiTSWithoutSDA<double>;
template class MEDCoupling::MEDFileNDTemplateFieldMultiTSWithoutSDA<int>;
template class MEDCoupling::MEDFileNDTemplateFieldMultiTSWithoutSDA<float>;
template class MEDCoupling::MEDFileTemplateFieldMultiTS<int>;
template class MEDCoupling::MEDFileTemplateFieldMultiTS<float>;
template class MEDCoupling::MEDFileTemplateFieldMultiTS<double>;
template class MEDCoupling::MEDFileNDTemplateFieldMultiTS<int>;
template class MEDCoupling::MEDFileNDTemplateFieldMultiTS<float>;

//= MEDFileAnyTypeFieldMultiTSWithoutSDA

MEDFileAnyTypeFieldMultiTSWithoutSDA::MEDFileAnyTypeFieldMultiTSWithoutSDA()
{
}

MEDFileAnyTypeFieldMultiTSWithoutSDA::MEDFileAnyTypeFieldMultiTSWithoutSDA(const std::string& fieldName, const std::string& meshName):MEDFileFieldNameScope(fieldName,meshName)
{
}

/*!
 * \param [in] fieldId field id in C mode
 */
MEDFileAnyTypeFieldMultiTSWithoutSDA::MEDFileAnyTypeFieldMultiTSWithoutSDA(med_idt fid, int fieldId, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities)
{
  med_field_type typcha;
  std::string dtunitOut,meshName;
  int nbOfStep(MEDFileAnyTypeField1TS::LocateField2(fid,fieldId,false,_name,typcha,_infos,dtunitOut,meshName));
  setMeshName(meshName);
  setDtUnit(dtunitOut.c_str());
  loadStructureOrStructureAndBigArraysRecursively(fid,nbOfStep,typcha,loadAll,ms,entities);
}

MEDFileAnyTypeFieldMultiTSWithoutSDA::MEDFileAnyTypeFieldMultiTSWithoutSDA(med_idt fid, const std::string& fieldName, const std::string& meshName, med_field_type fieldTyp, const std::vector<std::string>& infos, int nbOfStep, const std::string& dtunit, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities)
try:MEDFileFieldNameScope(fieldName,meshName),_infos(infos)
{
  setDtUnit(dtunit.c_str());
  loadStructureOrStructureAndBigArraysRecursively(fid,nbOfStep,fieldTyp,loadAll,ms,entities);
}
catch(INTERP_KERNEL::Exception& e)
{
    throw e;
}

std::size_t MEDFileAnyTypeFieldMultiTSWithoutSDA::getHeapMemorySizeWithoutChildren() const
{
  std::size_t ret(_mesh_name.capacity()+_name.capacity()+_infos.capacity()*sizeof(std::string)+_time_steps.capacity()*sizeof(MCAuto<MEDFileField1TSWithoutSDA>));
  for(std::vector<std::string>::const_iterator it=_infos.begin();it!=_infos.end();it++)
    ret+=(*it).capacity();
  return ret;
}

std::vector<const BigMemoryObject *> MEDFileAnyTypeFieldMultiTSWithoutSDA::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    ret.push_back((const MEDFileAnyTypeField1TSWithoutSDA *)*it);
  return ret;
}

/*!
 * If one of the id in [ \a startIds , \a endIds ) points to a null element, there is not throw. Simply, this empty element is added as if it were not
 * NULL.
 */
MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTSWithoutSDA::buildFromTimeStepIds(const int *startIds, const int *endIds) const
{
  MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> ret=createNew();
  ret->setInfo(_infos);
  int sz=(int)_time_steps.size();
  for(const int *id=startIds;id!=endIds;id++)
    {
      if(*id>=0 && *id<sz)
        {
          const MEDFileAnyTypeField1TSWithoutSDA *tse=_time_steps[*id];
          MCAuto<MEDFileAnyTypeField1TSWithoutSDA> tse2;
          if(tse)
            {
              tse->incrRef();
              tse2=(const_cast<MEDFileAnyTypeField1TSWithoutSDA *>(tse));
            }
          ret->pushBackTimeStep(tse2);
        }
      else
        {
          std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::buildFromTimeStepIds : At pos #" << std::distance(startIds,id) << " value is " << *id;
          oss << " ! Should be in [0," << sz << ") !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
    }
  if(ret->getNumberOfTS()>0)
    ret->synchronizeNameScope();
  ret->copyNameScope(*this);
  return ret.retn();
}

/*!
 * If one of the id in the input range points to a null element, there is not throw. Simply, this empty element is added as if it were not
 * NULL.
 */
MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTSWithoutSDA::buildFromTimeStepIds2(int bg, int end, int step) const
{
  static const char msg[]="MEDFileAnyTypeFieldMultiTSWithoutSDA::buildFromTimeStepIds2";
  int nbOfEntriesToKeep=DataArrayInt::GetNumberOfItemGivenBESRelative(bg,end,step,msg);
  MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> ret=createNew();
  ret->setInfo(_infos);
  int sz=(int)_time_steps.size();
  int j=bg;
  for(int i=0;i<nbOfEntriesToKeep;i++,j+=step)
    {
      if(j>=0 && j<sz)
        {
          const MEDFileAnyTypeField1TSWithoutSDA *tse=_time_steps[j];
          MCAuto<MEDFileAnyTypeField1TSWithoutSDA> tse2;
          if(tse)
            {
              tse->incrRef();
              tse2=(const_cast<MEDFileAnyTypeField1TSWithoutSDA *>(tse));
            }
          ret->pushBackTimeStep(tse2);
        }
      else
        {
          std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::buildFromTimeStepIds : At pos #" << i << " value is " << j;
          oss << " ! Should be in [0," << sz << ") !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
    }
  if(ret->getNumberOfTS()>0)
    ret->synchronizeNameScope();
  ret->copyNameScope(*this);
  return ret.retn();
}

MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTSWithoutSDA::partOfThisLyingOnSpecifiedTimeSteps(const std::vector< std::pair<int,int> >& timeSteps) const
{
  int id=0;
  MCAuto<DataArrayInt> ids=DataArrayInt::New(); ids->alloc(0,1);
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,id++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *cur(*it);
      if(!cur)
        continue;
      std::pair<int,int> p(cur->getIteration(),cur->getOrder());
      if(std::find(timeSteps.begin(),timeSteps.end(),p)!=timeSteps.end())
        ids->pushBackSilent(id);
    }
  return buildFromTimeStepIds(ids->begin(),ids->end());
}

MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTSWithoutSDA::partOfThisNotLyingOnSpecifiedTimeSteps(const std::vector< std::pair<int,int> >& timeSteps) const
{
  int id=0;
  MCAuto<DataArrayInt> ids=DataArrayInt::New(); ids->alloc(0,1);
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,id++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *cur(*it);
      if(!cur)
        continue;
      std::pair<int,int> p(cur->getIteration(),cur->getOrder());
      if(std::find(timeSteps.begin(),timeSteps.end(),p)==timeSteps.end())
        ids->pushBackSilent(id);
    }
  return buildFromTimeStepIds(ids->begin(),ids->end());
}

bool MEDFileAnyTypeFieldMultiTSWithoutSDA::presenceOfStructureElements() const
{
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    if((*it).isNotNull())
      if((*it)->presenceOfStructureElements())
        return true;
  return false;
}

bool MEDFileAnyTypeFieldMultiTSWithoutSDA::onlyStructureElements() const
{
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    if((*it).isNotNull())
      if(!(*it)->onlyStructureElements())
        return false;
  return true;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::killStructureElements()
{
  std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > ret;
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    if((*it).isNotNull())
      {
        if((*it)->presenceOfStructureElements())
          {
            if(!(*it)->onlyStructureElements())
              {
                (*it)->killStructureElements();
                ret.push_back(*it);
              }
          }
        else
          {
            ret.push_back(*it);
          }
      }
  _time_steps=ret;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::keepOnlyStructureElements()
{
  std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > ret;
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    if((*it).isNotNull())
      {
        if((*it)->presenceOfStructureElements())
          {
            if(!(*it)->onlyStructureElements())
              (*it)->keepOnlyStructureElements();
            ret.push_back(*it);
          }
      }
  _time_steps=ret;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::keepOnlyOnSE(const std::string& seName)
{
  std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > ret;
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    if((*it).isNotNull())
      (*it)->keepOnlyOnSE(seName);
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::getMeshSENames(std::vector< std::pair<std::string,std::string> >& ps) const
{
  std::vector< std::pair<std::string,std::string> > ps2;
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    if((*it).isNotNull())
      {
        (*it)->getMeshSENames(ps2);
        break;
      }
  if(ps2.empty())
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTSWithoutSDA::getMeshSENames : this appears to not contain SE only !");
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    if((*it).isNotNull())
      {
        std::vector< std::pair<std::string,std::string> > ps3;
        (*it)->getMeshSENames(ps3);
        if(ps2!=ps3)
          throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTSWithoutSDA::getMeshSENames : For the moment only homogeneous SE def through time managed !");
      }
  for(std::vector< std::pair<std::string,std::string> >::const_iterator it=ps2.begin();it!=ps2.end();it++)
    {
      std::vector< std::pair<std::string,std::string> >::iterator it2(std::find(ps.begin(),ps.end(),*it));
      if(it2==ps.end())
        ps.push_back(*it);
    }
}

bool MEDFileAnyTypeFieldMultiTSWithoutSDA::presenceOfMultiDiscPerGeoType() const
{
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *cur(*it);
      if(!cur)
        continue;
      if(cur->presenceOfMultiDiscPerGeoType())
        return true;
    }
  return false;
}

const std::vector<std::string>& MEDFileAnyTypeFieldMultiTSWithoutSDA::getInfo() const
{
  return _infos;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::setInfo(const std::vector<std::string>& info)
{
  _infos=info;
}

int MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepPos(int iteration, int order) const
{
  int ret=0;
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA>  >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,ret++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *pt(*it);
      if(pt->isDealingTS(iteration,order))
        return ret;
    }
  std::ostringstream oss; oss << "MEDFileFieldMultiTS::getTimeStepPos : Muli timestep field on time (" << iteration << "," << order << ") does not exist ! Available (iteration,order) are :\n";
  std::vector< std::pair<int,int> > vp=getIterations();
  for(std::vector< std::pair<int,int> >::const_iterator it2=vp.begin();it2!=vp.end();it2++)
    oss << "(" << (*it2).first << "," << (*it2).second << ") ";
  throw INTERP_KERNEL::Exception(oss.str());
}

const MEDFileAnyTypeField1TSWithoutSDA& MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepEntry(int iteration, int order) const
{
  return *_time_steps[getTimeStepPos(iteration,order)];
}

MEDFileAnyTypeField1TSWithoutSDA& MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepEntry(int iteration, int order)
{
  return *_time_steps[getTimeStepPos(iteration,order)];
}

bool MEDFileAnyTypeFieldMultiTSWithoutSDA::changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab)
{
  bool ret(false);
  for(std::vector< std::pair<std::string,std::string> >::const_iterator it=modifTab.begin();it!=modifTab.end();it++)
    {
      if((*it).first==getMeshName())
        {
          setMeshName((*it).second);
          ret=true;
        }
    }
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      MEDFileAnyTypeField1TSWithoutSDA *cur(*it);
      if(cur)
        ret=cur->changeMeshNames(modifTab) || ret;
    }
  return ret;
}

/*!
 * See doc at MEDFileField1TSWithoutSDA::getUndergroundDataArray
 */
DataArray *MEDFileAnyTypeFieldMultiTSWithoutSDA::getUndergroundDataArray(int iteration, int order) const
{
  return getTimeStepEntry(iteration,order).getUndergroundDataArray();
}

/*!
 * See doc at MEDFileField1TSWithoutSDA::getUndergroundDataArrayExt
 */
DataArray *MEDFileAnyTypeFieldMultiTSWithoutSDA::getUndergroundDataArrayExt(int iteration, int order, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const
{
  return getTimeStepEntry(iteration,order).getUndergroundDataArrayExt(entries);
}

bool MEDFileAnyTypeFieldMultiTSWithoutSDA::renumberEntitiesLyingOnMesh(const std::string& meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N,
                                                                       MEDFileFieldGlobsReal& glob)
{
  bool ret=false;
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      MEDFileAnyTypeField1TSWithoutSDA *f1ts(*it);
      if(f1ts)
        ret=f1ts->renumberEntitiesLyingOnMesh(meshName,oldCode,newCode,renumO2N,glob) || ret;
    }
  return ret;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::accept(MEDFileFieldVisitor& visitor) const
{
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    if((*it).isNotNull())
      {
        visitor.newTimeStepEntry(*it);
        (*it)->accept(visitor);
        visitor.endTimeStepEntry(*it);
      }
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::simpleRepr(int bkOffset, std::ostream& oss, int fmtsId) const
{
  std::string startLine(bkOffset,' ');
  oss << startLine << "Field multi time steps [Type=" << getTypeStr() << "]";
  if(fmtsId>=0)
    oss << " (" << fmtsId << ")";
  oss << " has the following name: \"" << _name << "\"." << std::endl;
  oss << startLine << "Field multi time steps has " << _infos.size() << " components with the following infos :" << std::endl;
  for(std::vector<std::string>::const_iterator it=_infos.begin();it!=_infos.end();it++)
    {
      oss << startLine << "  -  \"" << *it << "\"" << std::endl;
    }
  int i=0;
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,i++)
    {
      std::string chapter(17,'0'+i);
      oss << startLine << chapter << std::endl;
      const MEDFileAnyTypeField1TSWithoutSDA *cur=(*it);
      if(cur)
        cur->simpleRepr(bkOffset+2,oss,i);
      else
        oss << startLine << "  Field on one time step #" << i << " is not defined !" << std::endl;
      oss << startLine << chapter << std::endl;
    }
}

std::vector< std::pair<int,int> > MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeSteps(std::vector<double>& ret1) const
{
  std::size_t sz=_time_steps.size();
  std::vector< std::pair<int,int> > ret(sz);
  ret1.resize(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *f1ts=_time_steps[i];
      if(f1ts)
        {
          ret1[i]=f1ts->getTime(ret[i].first,ret[i].second);
        }
      else
        {
          std::ostringstream oss; oss << "MEDFileFieldMultiTSWithoutSDA::getTimeSteps : At rank #" << i << " time step is not defined. Invoke eraseEmptyTS method !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
    }
  return ret;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::pushBackTimeStep(MCAuto<MEDFileAnyTypeField1TSWithoutSDA>& tse)
{
  MEDFileAnyTypeField1TSWithoutSDA *tse2(tse);
  if(!tse2)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTSWithoutSDA::pushBackTimeStep : input content object is null !");
  checkCoherencyOfType(tse2);
  if(_time_steps.empty())
    {
      setName(tse2->getName());
      setMeshName(tse2->getMeshName());
      setInfo(tse2->getInfo());
    }
  checkThatComponentsMatch(tse2->getInfo());
  if(getDtUnit().empty() && !tse->getDtUnit().empty())
    setDtUnit(tse->getDtUnit());
  _time_steps.push_back(tse);
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::synchronizeNameScope()
{
  std::size_t nbOfCompo=_infos.size();
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      MEDFileAnyTypeField1TSWithoutSDA *cur=(*it);
      if(cur)
        {
          if((cur->getInfo()).size()!=nbOfCompo)
            {
              std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::synchronizeNameScope : Mismatch in the number of components of parts ! Should be " << nbOfCompo;
              oss << " ! but the field at iteration=" << cur->getIteration() << " order=" << cur->getOrder() << " has " << (cur->getInfo()).size() << " components !";
              throw INTERP_KERNEL::Exception(oss.str());
            }
          cur->copyNameScope(*this);
        }
    }
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::loadStructureOrStructureAndBigArraysRecursively(med_idt fid, int nbPdt, med_field_type fieldTyp, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities)
{
  _time_steps.resize(nbPdt);
  for(int i=0;i<nbPdt;i++)
    {
      std::vector< std::pair<int,int> > ts;
      med_int numdt=0,numo=0;
      med_float dt=0.0;
      MEDFILESAFECALLERRD0(MEDfieldComputingStepInfo,(fid,_name.c_str(),i+1,&numdt,&numo,&dt));
      switch(fieldTyp)
      {
        case MED_FLOAT64:
          {
            _time_steps[i]=MEDFileField1TSWithoutSDA::New(getName(),getMeshName(),i+1,numdt,numo,_infos);
            break;
          }
        case MED_INT32:
          {
            _time_steps[i]=MEDFileIntField1TSWithoutSDA::New(getName(),getMeshName(),i+1,numdt,numo,_infos);
            break;
          }
        case MED_FLOAT32:
          {
            _time_steps[i]=MEDFileFloatField1TSWithoutSDA::New(getName(),getMeshName(),i+1,numdt,numo,_infos);
            break;
          }
        case MED_INT:
          {
            if(sizeof(med_int)==sizeof(int))
              {
                _time_steps[i]=MEDFileIntField1TSWithoutSDA::New(getName(),getMeshName(),i+1,numdt,numo,_infos);
                break;
              }
          }
        default:
          throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTSWithoutSDA::loadStructureOrStructureAndBigArraysRecursively : managed field type are : FLOAT64, INT32, FLOAT32 !");
      }
      if(loadAll)
        _time_steps[i]->loadStructureAndBigArraysRecursively(fid,*this,ms,entities);
      else
        _time_steps[i]->loadOnlyStructureOfDataRecursively(fid,*this,ms,entities);
      synchronizeNameScope();
    }
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::writeLL(med_idt fid, const MEDFileWritable& opts) const
{
  if(_time_steps.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::writeLL : no time steps set !");
  checkThatNbOfCompoOfTSMatchThis();
  std::vector<std::string> infos(getInfo());
  int nbComp=infos.size();
  INTERP_KERNEL::AutoPtr<char> comp=MEDLoaderBase::buildEmptyString(nbComp*MED_SNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> unit=MEDLoaderBase::buildEmptyString(nbComp*MED_SNAME_SIZE);
  for(int i=0;i<nbComp;i++)
    {
      std::string info=infos[i];
      std::string c,u;
      MEDLoaderBase::splitIntoNameAndUnit(info,c,u);
      MEDLoaderBase::safeStrCpy2(c.c_str(),MED_SNAME_SIZE,comp+i*MED_SNAME_SIZE,opts.getTooLongStrPolicy());
      MEDLoaderBase::safeStrCpy2(u.c_str(),MED_SNAME_SIZE,unit+i*MED_SNAME_SIZE,opts.getTooLongStrPolicy());
    }
  if(_name.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::write : MED file does not accept field with empty name !");
  MEDFILESAFECALLERWR0(MEDfieldCr,(fid,_name.c_str(),getMEDFileFieldType(),nbComp,comp,unit,getDtUnit().c_str(),getMeshName().c_str()));
  int nbOfTS=_time_steps.size();
  for(int i=0;i<nbOfTS;i++)
    _time_steps[i]->writeLL(fid,opts,*this);
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::loadBigArraysRecursively(med_idt fid, const MEDFileFieldNameScope& nasc)
{
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      MEDFileAnyTypeField1TSWithoutSDA *elt(*it);
      if(elt)
        elt->loadBigArraysRecursively(fid,nasc);
    }
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::loadBigArraysRecursivelyIfNecessary(med_idt fid, const MEDFileFieldNameScope& nasc)
{
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      MEDFileAnyTypeField1TSWithoutSDA *elt(*it);
      if(elt)
        elt->loadBigArraysRecursivelyIfNecessary(fid,nasc);
    }
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::unloadArrays()
{
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      MEDFileAnyTypeField1TSWithoutSDA *elt(*it);
      if(elt)
        elt->unloadArrays();
    }
}

int MEDFileAnyTypeFieldMultiTSWithoutSDA::getNumberOfTS() const
{
  return _time_steps.size();
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::eraseEmptyTS()
{
  std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA>  > newTS;
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA>  >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *tmp=(*it);
      if(tmp)
        newTS.push_back(*it);
    }
  _time_steps=newTS;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::eraseTimeStepIds(const int *startIds, const int *endIds)
{
  std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > newTS;
  int maxId=(int)_time_steps.size();
  int ii=0;
  std::set<int> idsToDel;
  for(const int *id=startIds;id!=endIds;id++,ii++)
    {
      if(*id>=0 && *id<maxId)
        {
          idsToDel.insert(*id);
        }
      else
        {
          std::ostringstream oss; oss << "MEDFileFieldMultiTSWithoutSDA::eraseTimeStepIds : At pos #" << ii << " request for id=" << *id << " not in [0," << maxId << ") !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
    }
  for(int iii=0;iii<maxId;iii++)
    if(idsToDel.find(iii)==idsToDel.end())
      newTS.push_back(_time_steps[iii]);
  _time_steps=newTS;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::eraseTimeStepIds2(int bg, int end, int step)
{
  static const char msg[]="MEDFileAnyTypeFieldMultiTSWithoutSDA::eraseTimeStepIds2";
  int nbOfEntriesToKill=DataArrayInt::GetNumberOfItemGivenBESRelative(bg,end,step,msg);
  if(nbOfEntriesToKill==0)
    return ;
  std::size_t sz=_time_steps.size();
  std::vector<bool> b(sz,true);
  int j=bg;
  for(int i=0;i<nbOfEntriesToKill;i++,j+=step)
    b[j]=false;
  std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > newTS;
  for(std::size_t i=0;i<sz;i++)
    if(b[i])
      newTS.push_back(_time_steps[i]);
  _time_steps=newTS;
}

int MEDFileAnyTypeFieldMultiTSWithoutSDA::getPosOfTimeStep(int iteration, int order) const
{
  int ret=0;
  std::ostringstream oss; oss << "MEDFileFieldMultiTSWithoutSDA::getPosOfTimeStep : No such time step (" << iteration << "," << order << ") !\nPossibilities are : "; 
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA>  >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,ret++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *tmp(*it);
      if(tmp)
        {
          int it2,ord;
          tmp->getTime(it2,ord);
          if(it2==iteration && order==ord)
            return ret;
          else
            oss << "(" << it2 << ","  << ord << "), ";
        }
    }
  throw INTERP_KERNEL::Exception(oss.str());
}

int MEDFileAnyTypeFieldMultiTSWithoutSDA::getPosGivenTime(double time, double eps) const
{
  int ret=0;
  std::ostringstream oss; oss << "MEDFileFieldMultiTSWithoutSDA::getPosGivenTime : No such time step " << time << "! \nPossibilities are : ";
  oss.precision(15);
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA>  >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,ret++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *tmp(*it);
      if(tmp)
        {
          int it2,ord;
          double ti=tmp->getTime(it2,ord);
          if(fabs(time-ti)<eps)
            return ret;
          else
            oss << ti << ", ";
        }
    }
  throw INTERP_KERNEL::Exception(oss.str());
}

std::vector< std::pair<int,int> > MEDFileAnyTypeFieldMultiTSWithoutSDA::getIterations() const
{
  int lgth=_time_steps.size();
  std::vector< std::pair<int,int> > ret(lgth);
  for(int i=0;i<lgth;i++)
    _time_steps[i]->fillIteration(ret[i]);
  return ret;
}

/*!
 * This method has 3 inputs 'iteration' 'order' 'mname'. 'mname' can be null if the user is the general case where there is only one meshName lying on 'this'
 * This method returns two things.
 * - The absolute dimension of 'this' in first parameter. 
 * - The available ext levels relative to the absolute dimension returned in first parameter. These relative levels are relative
 *   to the first output parameter. The values in 'levs' will be returned in decreasing order.
 *
 * This method is designed for MEDFileFieldMultiTS instances that have a discritization ON_CELLS, ON_GAUSS_NE and ON_GAUSS.
 * Only these 3 discretizations will be taken into account here.
 *
 * If 'this' is empty this method will throw an INTERP_KERNEL::Exception.
 * If there is \b only node fields defined in 'this' -1 is returned and 'levs' output parameter will be empty. In this
 * case the caller has to know the underlying mesh it refers to. By default it is the level 0 of the corresponding mesh.
 *
 * This method is useful to make the link between meshDimension of the underlying mesh in 'this' and the levels on 'this'.
 * It is possible (even if it is not common) that the highest level in 'this' were not equal to the meshDimension of the underlying mesh in 'this'.
 * 
 * Let's consider the typical following case :
 * - a mesh 'm1' has a meshDimension 3 and has the following non empty levels
 * [0,-1,-2] for example 'm1' lies on TETRA4, HEXA8 TRI3 and SEG2
 * - 'f1' lies on 'm1' and is defined on 3D and 1D cells for example
 *   TETRA4 and SEG2
 * - 'f2' lies on 'm1' too and is defined on 2D and 1D cells for example TRI3 and SEG2
 *
 * In this case f1->getNonEmptyLevelsExt will return (3,[0,-2]) and f2->getNonEmptyLevelsExt will return (2,[0,-1])
 * 
 * To retrieve the highest level of f1 it should be done, f1->getFieldAtLevel(ON_CELLS,3-3+0);//absDim-meshDim+relativeLev
 * To retrieve the lowest level of f1 it should be done, f1->getFieldAtLevel(ON_CELLS,3-3+(-2));//absDim-meshDim+relativeLev
 * To retrieve the highest level of f2 it should be done, f1->getFieldAtLevel(ON_CELLS,2-3+0);//absDim-meshDim+relativeLev
 * To retrieve the lowest level of f2 it should be done, f1->getFieldAtLevel(ON_CELLS,2-3+(-1));//absDim-meshDim+relativeLev
 */
int MEDFileAnyTypeFieldMultiTSWithoutSDA::getNonEmptyLevels(int iteration, int order, const std::string& mname, std::vector<int>& levs) const
{
  return getTimeStepEntry(iteration,order).getNonEmptyLevels(mname,levs);
}

const MEDFileAnyTypeField1TSWithoutSDA *MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepAtPos2(int pos) const
{
  if(pos<0 || pos>=(int)_time_steps.size())
    {
      std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepAtPos2 : request for pos #" << pos << " whereas should be in [0," << _time_steps.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  const MEDFileAnyTypeField1TSWithoutSDA *item=_time_steps[pos];
  if(item==0)
    {
      std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepAtPos2 : request for pos #" << pos << ", this pos id exists but the underlying Field1TS is null !";
      oss << "\nTry to use following method eraseEmptyTS !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  return item;
}

MEDFileAnyTypeField1TSWithoutSDA *MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepAtPos2(int pos)
{
  if(pos<0 || pos>=(int)_time_steps.size())
    {
      std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepAtPos2 : request for pos #" << pos << " whereas should be in [0," << _time_steps.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  MEDFileAnyTypeField1TSWithoutSDA *item=_time_steps[pos];
  if(item==0)
    {
      std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepAtPos2 : request for pos #" << pos << ", this pos id exists but the underlying Field1TS is null !";
      oss << "\nTry to use following method eraseEmptyTS !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  return item;
}

std::vector<std::string> MEDFileAnyTypeFieldMultiTSWithoutSDA::getPflsReallyUsed2() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MCAuto< MEDFileAnyTypeField1TSWithoutSDA > >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getPflsReallyUsed2();
      for(std::vector<std::string>::const_iterator it2=tmp.begin();it2!=tmp.end();it2++)
        if(ret2.find(*it2)==ret2.end())
          {
            ret.push_back(*it2);
            ret2.insert(*it2);
          }
    }
  return ret;
}

std::vector<std::string> MEDFileAnyTypeFieldMultiTSWithoutSDA::getLocsReallyUsed2() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MCAuto< MEDFileAnyTypeField1TSWithoutSDA > >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getLocsReallyUsed2();
      for(std::vector<std::string>::const_iterator it2=tmp.begin();it2!=tmp.end();it2++)
        if(ret2.find(*it2)==ret2.end())
          {
            ret.push_back(*it2);
            ret2.insert(*it2);
          }
    }
  return ret;
}

std::vector<std::string> MEDFileAnyTypeFieldMultiTSWithoutSDA::getPflsReallyUsedMulti2() const
{
  std::vector<std::string> ret;
  for(std::vector< MCAuto< MEDFileAnyTypeField1TSWithoutSDA > >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getPflsReallyUsedMulti2();
      ret.insert(ret.end(),tmp.begin(),tmp.end());
    }
  return ret;
}

std::vector<std::string> MEDFileAnyTypeFieldMultiTSWithoutSDA::getLocsReallyUsedMulti2() const
{
  std::vector<std::string> ret;
  for(std::vector< MCAuto< MEDFileAnyTypeField1TSWithoutSDA > >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getLocsReallyUsedMulti2();
      ret.insert(ret.end(),tmp.begin(),tmp.end());
    }
  return ret;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::changePflsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
{
  for(std::vector< MCAuto< MEDFileAnyTypeField1TSWithoutSDA > >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    (*it)->changePflsRefsNamesGen2(mapOfModif);
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::changeLocsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
{
  for(std::vector< MCAuto< MEDFileAnyTypeField1TSWithoutSDA > >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    (*it)->changeLocsRefsNamesGen2(mapOfModif);
}

std::vector< std::vector<TypeOfField> > MEDFileAnyTypeFieldMultiTSWithoutSDA::getTypesOfFieldAvailable() const
{
  int lgth=_time_steps.size();
  std::vector< std::vector<TypeOfField> > ret(lgth);
  for(int i=0;i<lgth;i++)
    _time_steps[i]->fillTypesOfFieldAvailable(ret[i]);
  return ret;
}

/*!
 * entry point for users that want to iterate into MEDFile DataStructure without any overhead.
 */
std::vector< std::vector< std::pair<int,int> > > MEDFileAnyTypeFieldMultiTSWithoutSDA::getFieldSplitedByType(int iteration, int order, const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const
{
  return getTimeStepEntry(iteration,order).getFieldSplitedByType(mname,types,typesF,pfls,locs);
}

MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTSWithoutSDA::deepCopy() const
{
  MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> ret=shallowCpy();
  std::size_t i=0;
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,i++)
    {
      if((const MEDFileAnyTypeField1TSWithoutSDA *)*it)
        ret->_time_steps[i]=(*it)->deepCopy();
    }
  return ret.retn();
}

std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > MEDFileAnyTypeFieldMultiTSWithoutSDA::splitComponents() const
{
  std::size_t sz(_infos.size()),sz2(_time_steps.size());
  std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > ret(sz);
  std::vector< std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > > ts(sz2);
  for(std::size_t i=0;i<sz;i++)
    {
      ret[i]=shallowCpy();
      ret[i]->_infos.resize(1); ret[i]->_infos[0]=_infos[i];
    }
  for(std::size_t i=0;i<sz2;i++)
    {
      std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > ret1=_time_steps[i]->splitComponents();
      if(ret1.size()!=sz)
        {
          std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::splitComponents : At rank #" << i << " number of components is " << ret1.size() << " whereas it should be for all time steps " << sz << " !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
      ts[i]=ret1;
    }
  for(std::size_t i=0;i<sz;i++)
    for(std::size_t j=0;j<sz2;j++)
      ret[i]->_time_steps[j]=ts[j][i];
  return ret;
}

/*!
 * This method splits into discretization each time steps in \a this.
 * ** WARNING ** the returned instances are not compulsory defined on the same time steps series !
 */
std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > MEDFileAnyTypeFieldMultiTSWithoutSDA::splitDiscretizations() const
{
  std::size_t sz(_time_steps.size());
  std::vector< std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > > items(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *timeStep(_time_steps[i]);
      if(!timeStep)
        {
          std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::splitDiscretizations : time step #" << i << " is null !"; 
          throw INTERP_KERNEL::Exception(oss.str());
        }
      items[i]=timeStep->splitDiscretizations();  
    }
  //
  std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > ret;
  std::vector< std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > > ret2;
  std::vector< TypeOfField > types;
  for(std::vector< std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > >::const_iterator it0=items.begin();it0!=items.end();it0++)
    for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it1=(*it0).begin();it1!=(*it0).end();it1++)
      {
        std::vector<TypeOfField> ts=(*it1)->getTypesOfFieldAvailable();
        if(ts.size()!=1)
          throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTSWithoutSDA::splitDiscretizations : it appears that the splitting of MEDFileAnyTypeField1TSWithoutSDA::splitDiscretizations has returned invalid result !");
        std::vector< TypeOfField >::iterator it2=std::find(types.begin(),types.end(),ts[0]);
        if(it2==types.end())
          types.push_back(ts[0]);
      }
  ret.resize(types.size()); ret2.resize(types.size());
  for(std::vector< std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > >::const_iterator it0=items.begin();it0!=items.end();it0++)
    for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it1=(*it0).begin();it1!=(*it0).end();it1++)
      {
        TypeOfField typ=(*it1)->getTypesOfFieldAvailable()[0];
        std::size_t pos=std::distance(types.begin(),std::find(types.begin(),types.end(),typ));
        ret2[pos].push_back(*it1);
      }
  for(std::size_t i=0;i<types.size();i++)
    {
      MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> elt(createNew());
      for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::iterator it1=ret2[i].begin();it1!=ret2[i].end();it1++)
        elt->pushBackTimeStep(*it1);//also updates infos in elt
      ret[i]=elt;
      elt->MEDFileFieldNameScope::operator=(*this);
    }
  return ret;
}

/*!
 * Contrary to splitDiscretizations method this method makes the hypothesis that the times series are **NOT** impacted by the splitting of multi discretization.
 */
std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > MEDFileAnyTypeFieldMultiTSWithoutSDA::splitMultiDiscrPerGeoTypes() const
{
  std::size_t sz(_time_steps.size());
  std::vector< std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > > items(sz);
  std::size_t szOut(std::numeric_limits<std::size_t>::max());
  for(std::size_t i=0;i<sz;i++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *timeStep(_time_steps[i]);
      if(!timeStep)
        {
          std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::splitMultiDiscrPerGeoTypes : time step #" << i << " is null !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
      items[i]=timeStep->splitMultiDiscrPerGeoTypes();
      if(szOut==std::numeric_limits<std::size_t>::max())
        szOut=items[i].size();
      else
        if(items[i].size()!=szOut)
          throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTSWithoutSDA::splitMultiDiscrPerGeoTypes : The splitting per discretization is expected to be same among time steps !");
    }
  if(szOut==std::numeric_limits<std::size_t>::max())
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTSWithoutSDA::splitMultiDiscrPerGeoTypes : empty field !");
  std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > ret(szOut);
  for(std::size_t i=0;i<szOut;i++)
    {
      MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> elt(createNew());
      for(std::size_t j=0;j<sz;j++)
        elt->pushBackTimeStep(items[j][i]);
      ret[i]=elt;
      elt->MEDFileFieldNameScope::operator=(*this);
    }
  return ret;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::copyTinyInfoFrom(const MEDCouplingFieldDouble *field, const DataArray *arr)
{
  setName(field->getName());
  if(field->getMesh())
    setMeshName(field->getMesh()->getName());
  if(_name.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::copyTinyInfoFrom : unsupported fields with no name in MED file !");
  if(!arr)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::copyTinyInfoFrom : no array set !");
  _infos=arr->getInfoOnComponents();
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::checkCoherencyOfTinyInfo(const MEDCouplingFieldDouble *field, const DataArray *arr) const
{
  static const char MSG[]="MEDFileFieldMultiTSWithoutSDA::checkCoherencyOfTinyInfo : invalid ";
  if(_name!=field->getName())
    {
      std::ostringstream oss; oss << MSG << "name ! should be \"" << _name;
      oss << "\" and it is set in input field to \"" << field->getName() << "\" !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  if(!arr)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::checkCoherencyOfTinyInfo : no array set !");
  checkThatComponentsMatch(arr->getInfoOnComponents());
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::checkThatComponentsMatch(const std::vector<std::string>& compos) const
{
  static const char MSG[]="MEDFileFieldMultiTSWithoutSDA::checkThatComponentsMatch : ";
  if(getInfo().size()!=compos.size())
    {
      std::ostringstream oss; oss << MSG << "mismatch of number of components between this (" << getInfo().size() << ") and ";
      oss << " number of components of element to append (" << compos.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  if(_infos!=compos)
    {
      std::ostringstream oss; oss << MSG << "components have same size but are different ! should be \"";
      std::copy(_infos.begin(),_infos.end(),std::ostream_iterator<std::string>(oss,", "));
      oss << " But compo in input fields are : ";
      std::copy(compos.begin(),compos.end(),std::ostream_iterator<std::string>(oss,", "));
      oss << " !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::checkThatNbOfCompoOfTSMatchThis() const
{
  std::size_t sz=_infos.size();
  int j=0;
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,j++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *elt(*it);
      if(elt)
        if(elt->getInfo().size()!=sz)
          {
            std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::checkThatNbOfCompoOfTSMatchThis : At pos #" << j << " the number of components is equal to ";
            oss << elt->getInfo().size() << " whereas it is expected to be equal to " << sz << " !";
            throw INTERP_KERNEL::Exception(oss.str());
          }
    }
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::appendFieldNoProfileSBT(const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob)
{
  if(!field)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTSWithoutSDA::appendFieldNoProfileSBT : input field is NULL !");
  if(!_time_steps.empty())
    checkCoherencyOfTinyInfo(field,arr);
  MCAuto<MEDFileAnyTypeField1TSWithoutSDA> obj(createNew1TSWithoutSDAEmptyInstance());
  {
    MCAuto<MEDCouplingFieldTemplate> ft(MEDCouplingFieldTemplate::New(*field));
    obj->setFieldNoProfileSBT(field->timeDiscrSafe(),ft,arr,glob,*this);
  }
  copyTinyInfoFrom(field,arr);
  _time_steps.push_back(obj);
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::appendFieldProfile(const MEDCouplingFieldDouble *field, const DataArray *arr, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile, MEDFileFieldGlobsReal& glob, bool smartPflKiller)
{
  if(!field)
    throw INTERP_KERNEL::Exception("MEDFileIntFieldMultiTSWithoutSDA::appendFieldNoProfileSBT : input field is NULL !");
  if(!_time_steps.empty())
    checkCoherencyOfTinyInfo(field,arr);
  MCAuto<MEDFileAnyTypeField1TSWithoutSDA> obj(createNew1TSWithoutSDAEmptyInstance());
  {
    MCAuto<MEDCouplingFieldTemplate> ft(MEDCouplingFieldTemplate::NewWithoutCheck(*field));
    obj->setFieldProfile(field->timeDiscrSafe(),ft,arr,mesh,meshDimRelToMax,profile,glob,*this,smartPflKiller);
  }
  copyTinyInfoFrom(field,arr);
  setMeshName(obj->getMeshName());
  _time_steps.push_back(obj);
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::setIteration(int i, MCAuto<MEDFileAnyTypeField1TSWithoutSDA> ts)
{
  int sz=(int)_time_steps.size();
  if(i<0 || i>=sz)
    {
      std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::setIteration : trying to set element at place #" << i << " should be in [0," << sz << ") !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  const MEDFileAnyTypeField1TSWithoutSDA *tsPtr(ts);
  if(tsPtr)
    {
      if(tsPtr->getNumberOfComponents()!=(int)_infos.size())
        {
          std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::setIteration : trying to set element with " << tsPtr->getNumberOfComponents() << " components ! Should be " << _infos.size() <<  " !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
    }
  _time_steps[i]=ts;
}

//= MEDFileFieldMultiTSWithoutSDA

/*!
 * entry point for users that want to iterate into MEDFile DataStructure with a reduced overhead because output arrays are extracted (created) specially
 * for the call of this method. That's why the DataArrayDouble instance in returned vector of vector should be dealed by the caller.
 */
std::vector< std::vector<DataArrayDouble *> > MEDFileFieldMultiTSWithoutSDA::getFieldSplitedByType2(int iteration, int order, const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const
{
  const MEDFileAnyTypeField1TSWithoutSDA& myF1TS=getTimeStepEntry(iteration,order);
  const MEDFileField1TSWithoutSDA *myF1TSC=dynamic_cast<const MEDFileField1TSWithoutSDA *>(&myF1TS);
  if(!myF1TSC)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::getFieldSplitedByType2 : mismatch of type of field expecting FLOAT64 !");
  return myF1TSC->getFieldSplitedByType2(mname,types,typesF,pfls,locs);
}

MEDFileIntFieldMultiTSWithoutSDA *MEDFileFieldMultiTSWithoutSDA::convertToInt() const
{
  MCAuto<MEDFileIntFieldMultiTSWithoutSDA> ret(new MEDFileIntFieldMultiTSWithoutSDA);
  ret->MEDFileAnyTypeFieldMultiTSWithoutSDA::operator =(*this);
  int i=0;
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,i++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *eltToConv(*it);
      if(eltToConv)
        {
          const MEDFileField1TSWithoutSDA *eltToConvC=dynamic_cast<const MEDFileField1TSWithoutSDA *>(eltToConv);
          if(!eltToConvC)
            throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::convertToInt : presence of an invalid 1TS type ! Should be of type FLOAT64 !");
          MCAuto<MEDFileAnyTypeField1TSWithoutSDA> elt=eltToConvC->convertToInt();
          ret->setIteration(i,elt);
        }
    }
  return ret.retn();
}

//= MEDFileAnyTypeFieldMultiTS

MEDFileAnyTypeFieldMultiTS::MEDFileAnyTypeFieldMultiTS()
{
}

MEDFileAnyTypeFieldMultiTS::MEDFileAnyTypeFieldMultiTS(med_idt fid, bool loadAll, const MEDFileMeshes *ms)
try:MEDFileFieldGlobsReal(fid)
{
  _content=BuildContentFrom(fid,loadAll,ms);
  loadGlobals(fid);
}
catch(INTERP_KERNEL::Exception& e)
{
    throw e;
}

MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTS::BuildContentFrom(med_idt fid, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities)
{
  med_field_type typcha;
  std::vector<std::string> infos;
  std::string dtunit;
  std::string meshName;
  int i(-1);
  MEDFileAnyTypeField1TS::LocateField(fid,fieldName,i,typcha,infos,dtunit,meshName);
  MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> ret;
  switch(typcha)
  {
    case MED_FLOAT64:
      {
        ret=new MEDFileFieldMultiTSWithoutSDA(fid,i,loadAll,ms,entities);
        break;
      }
    case MED_INT32:
      {
        ret=new MEDFileIntFieldMultiTSWithoutSDA(fid,i,loadAll,ms,entities);
        break;
      }
    case MED_FLOAT32:
      {
        ret=new MEDFileFloatFieldMultiTSWithoutSDA(fid,i,loadAll,ms,entities);
        break;
      }
    case MED_INT:
      {
        if(sizeof(med_int)==sizeof(int))
          {
            ret=new MEDFileIntFieldMultiTSWithoutSDA(fid,i,loadAll,ms,entities);
            break;
          }
      }
    default:
      {
        std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTS::BuildContentFrom(fid,fieldName) : file \'" << FileNameFromFID(fid) << "\' contains field with name \'" << fieldName << "\' but the type of field is not in [MED_FLOAT64, MED_INT32, MED_FLOAT32] !";
        throw INTERP_KERNEL::Exception(oss.str());
      }
  }
  ret->setMeshName(meshName);
  ret->setDtUnit(dtunit.c_str());
  return ret.retn();
}

MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTS::BuildContentFrom(med_idt fid, bool loadAll, const MEDFileMeshes *ms)
{
  med_field_type typcha;
  //
  std::vector<std::string> infos;
  std::string dtunit,fieldName,meshName;
  MEDFileAnyTypeField1TS::LocateField2(fid,0,true,fieldName,typcha,infos,dtunit,meshName);
  MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> ret;
  switch(typcha)
  {
    case MED_FLOAT64:
      {
        ret=new MEDFileFieldMultiTSWithoutSDA(fid,0,loadAll,ms,0);
        break;
      }
    case MED_INT32:
      {
        ret=new MEDFileIntFieldMultiTSWithoutSDA(fid,0,loadAll,ms,0);
        break;
      }
    case MED_FLOAT32:
      {
        ret=new MEDFileFloatFieldMultiTSWithoutSDA(fid,0,loadAll,ms,0);
        break;
      }
    case MED_INT:
      {
        if(sizeof(med_int)==sizeof(int))
          {
            ret=new MEDFileIntFieldMultiTSWithoutSDA(fid,0,loadAll,ms,0);
            break;
          }
      }
    default:
      {
        std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTS::BuildContentFrom(fid) : file \'" << FileNameFromFID(fid) << "\' contains field with name \'" << fieldName << "\' but the type of the first field is not in [MED_FLOAT64, MED_INT32, MED_FLOAT32] !";
        throw INTERP_KERNEL::Exception(oss.str());
      }
  }
  ret->setMeshName(meshName);
  ret->setDtUnit(dtunit.c_str());
  return ret.retn();
}

MEDFileAnyTypeFieldMultiTS *MEDFileAnyTypeFieldMultiTS::BuildNewInstanceFromContent(MEDFileAnyTypeFieldMultiTSWithoutSDA *c)
{
  if(!c)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::BuildNewInstanceFromContent : empty content in input : unable to build a new instance !");
  if(dynamic_cast<const MEDFileFieldMultiTSWithoutSDA *>(c))
    {
      MCAuto<MEDFileFieldMultiTS> ret(MEDFileFieldMultiTS::New());
      ret->_content=c;  c->incrRef();
      return ret.retn();
    }
  if(dynamic_cast<const MEDFileIntFieldMultiTSWithoutSDA *>(c))
    {
      MCAuto<MEDFileIntFieldMultiTS> ret(MEDFileIntFieldMultiTS::New());
      ret->_content=c;  c->incrRef();
      return ret.retn();
    }
  if(dynamic_cast<const MEDFileFloatFieldMultiTSWithoutSDA *>(c))
    {
      MCAuto<MEDFileFloatFieldMultiTS> ret(MEDFileFloatFieldMultiTS::New());
      ret->_content=c;  c->incrRef();
      return ret.retn();
    }
  throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::BuildNewInstanceFromContent : internal error ! a content of type different from FLOAT64 FLOAT32 and INT32 has been built but not intercepted !");
}

MEDFileAnyTypeFieldMultiTS *MEDFileAnyTypeFieldMultiTS::BuildNewInstanceFromContent(MEDFileAnyTypeFieldMultiTSWithoutSDA *c, med_idt fid)
{
  MEDFileAnyTypeFieldMultiTS *ret(BuildNewInstanceFromContent(c));
  std::string fileName(FileNameFromFID(fid));
  ret->setFileName(fileName);
  return ret;
}

MEDFileAnyTypeFieldMultiTS::MEDFileAnyTypeFieldMultiTS(med_idt fid, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities)
try:MEDFileFieldGlobsReal(fid)
{
  _content=BuildContentFrom(fid,fieldName,loadAll,ms,entities);
  loadGlobals(fid);
}
catch(INTERP_KERNEL::Exception& e)
{
    throw e;
}

//= MEDFileAnyTypeFieldMultiTS

/*!
 * Returns a new instance of MEDFileFieldMultiTS or MEDFileIntFieldMultiTS holding data of the first field
 * that has been read from a specified MED file.
 *  \param [in] fileName - the name of the MED file to read.
 *  \return MEDFileFieldMultiTS * - a new instance of MEDFileFieldMultiTS or MEDFileIntFieldMultiTS. The caller
 *          is to delete this field using decrRef() as it is no more needed.
 *  \throw If reading the file fails.
 */
MEDFileAnyTypeFieldMultiTS *MEDFileAnyTypeFieldMultiTS::New(const std::string& fileName, bool loadAll)
{
  MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(fileName));
  return New(fid,loadAll);
}

MEDFileAnyTypeFieldMultiTS *MEDFileAnyTypeFieldMultiTS::New(med_idt fid, bool loadAll)
{
  MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> c(BuildContentFrom(fid,loadAll,0));
  MCAuto<MEDFileAnyTypeFieldMultiTS> ret(BuildNewInstanceFromContent(c,fid));
  ret->loadGlobals(fid);
  return ret.retn();
}

/*!
 * Returns a new instance of MEDFileFieldMultiTS or MEDFileIntFieldMultiTS holding data of a given field
 * that has been read from a specified MED file.
 *  \param [in] fileName - the name of the MED file to read.
 *  \param [in] fieldName - the name of the field to read.
 *  \return MEDFileFieldMultiTS * - a new instance of MEDFileFieldMultiTS or MEDFileIntFieldMultiTS. The caller
 *          is to delete this field using decrRef() as it is no more needed.
 *  \throw If reading the file fails.
 *  \throw If there is no field named \a fieldName in the file.
 */
MEDFileAnyTypeFieldMultiTS *MEDFileAnyTypeFieldMultiTS::New(const std::string& fileName, const std::string& fieldName, bool loadAll)
{
  MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(fileName));
  return New(fid,fieldName,loadAll);
}

MEDFileAnyTypeFieldMultiTS *MEDFileAnyTypeFieldMultiTS::New(med_idt fid, const std::string& fieldName, bool loadAll)
{
  MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> c(BuildContentFrom(fid,fieldName,loadAll,0,0));
  MCAuto<MEDFileAnyTypeFieldMultiTS> ret(BuildNewInstanceFromContent(c,fid));
  ret->loadGlobals(fid);
  return ret.retn();
}

/*!
 * This constructor is a shallow copy constructor. If \a shallowCopyOfContent is true the content of \a other is shallow copied.
 * If \a shallowCopyOfContent is false, \a other is taken to be the content of \a this.
 *
 * \warning this is a shallow copy constructor
 */
MEDFileAnyTypeFieldMultiTS::MEDFileAnyTypeFieldMultiTS(const MEDFileAnyTypeFieldMultiTSWithoutSDA& other, bool shallowCopyOfContent)
{
  if(!shallowCopyOfContent)
    {
      const MEDFileAnyTypeFieldMultiTSWithoutSDA *otherPtr(&other);
      otherPtr->incrRef();
      _content=const_cast<MEDFileAnyTypeFieldMultiTSWithoutSDA *>(otherPtr);
    }
  else
    {
      _content=other.shallowCpy();
    }
}

MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTS::contentNotNullBase()
{
  MEDFileAnyTypeFieldMultiTSWithoutSDA *ret=_content;
  if(!ret)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS : content is expected to be not null !");
  return ret;
}

const MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTS::contentNotNullBase() const
{
  const MEDFileAnyTypeFieldMultiTSWithoutSDA *ret=_content;
  if(!ret)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS : const content is expected to be not null !");
  return ret;
}

std::vector<std::string> MEDFileAnyTypeFieldMultiTS::getPflsReallyUsed() const
{
  return contentNotNullBase()->getPflsReallyUsed2();
}

std::vector<std::string> MEDFileAnyTypeFieldMultiTS::getLocsReallyUsed() const
{
  return contentNotNullBase()->getLocsReallyUsed2();
}

std::vector<std::string> MEDFileAnyTypeFieldMultiTS::getPflsReallyUsedMulti() const
{
  return contentNotNullBase()->getPflsReallyUsedMulti2();
}

std::vector<std::string> MEDFileAnyTypeFieldMultiTS::getLocsReallyUsedMulti() const
{
  return contentNotNullBase()->getLocsReallyUsedMulti2();
}

void MEDFileAnyTypeFieldMultiTS::changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
{
  contentNotNullBase()->changePflsRefsNamesGen2(mapOfModif);
}

void MEDFileAnyTypeFieldMultiTS::changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
{
  contentNotNullBase()->changeLocsRefsNamesGen2(mapOfModif);
}

int MEDFileAnyTypeFieldMultiTS::getNumberOfTS() const
{
  return contentNotNullBase()->getNumberOfTS();
}

void MEDFileAnyTypeFieldMultiTS::eraseEmptyTS()
{
  contentNotNullBase()->eraseEmptyTS();
}

void MEDFileAnyTypeFieldMultiTS::eraseTimeStepIds(const int *startIds, const int *endIds)
{
  contentNotNullBase()->eraseTimeStepIds(startIds,endIds);
}

void MEDFileAnyTypeFieldMultiTS::eraseTimeStepIds2(int bg, int end, int step)
{
  contentNotNullBase()->eraseTimeStepIds2(bg,end,step);
}

MEDFileAnyTypeFieldMultiTS *MEDFileAnyTypeFieldMultiTS::buildSubPart(const int *startIds, const int *endIds) const
{
  MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> c=contentNotNullBase()->buildFromTimeStepIds(startIds,endIds);
  MCAuto<MEDFileAnyTypeFieldMultiTS> ret=shallowCpy();
  ret->_content=c;
  return ret.retn();
}

MEDFileAnyTypeFieldMultiTS *MEDFileAnyTypeFieldMultiTS::buildSubPartSlice(int bg, int end, int step) const
{
  MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> c=contentNotNullBase()->buildFromTimeStepIds2(bg,end,step);
  MCAuto<MEDFileAnyTypeFieldMultiTS> ret=shallowCpy();
  ret->_content=c;
  return ret.retn();
}

std::vector< std::pair<int,int> > MEDFileAnyTypeFieldMultiTS::getIterations() const
{
  return contentNotNullBase()->getIterations();
}

void MEDFileAnyTypeFieldMultiTS::pushBackTimeSteps(const std::vector<MEDFileAnyTypeField1TS *>& f1ts)
{
  for(std::vector<MEDFileAnyTypeField1TS *>::const_iterator it=f1ts.begin();it!=f1ts.end();it++)
    pushBackTimeStep(*it);
}

void MEDFileAnyTypeFieldMultiTS::pushBackTimeSteps(MEDFileAnyTypeFieldMultiTS *fmts)
{
  if(!fmts)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::pushBackTimeSteps : Input fmts is NULL !");
  int nbOfTS(fmts->getNumberOfTS());
  for(int i=0;i<nbOfTS;i++)
    {
      MCAuto<MEDFileAnyTypeField1TS> elt(fmts->getTimeStepAtPos(i));
      pushBackTimeStep(elt);
    }
}

void MEDFileAnyTypeFieldMultiTS::pushBackTimeStep(MEDFileAnyTypeField1TS *f1ts)
{
  if(!f1ts)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTSWithoutSDA::pushBackTimeStep : input pointer is NULL !");
  checkCoherencyOfType(f1ts);
  f1ts->incrRef();
  MCAuto<MEDFileAnyTypeField1TS> f1tsSafe(f1ts);
  MEDFileAnyTypeField1TSWithoutSDA *c=f1ts->contentNotNullBase();
  c->incrRef();
  MCAuto<MEDFileAnyTypeField1TSWithoutSDA> cSafe(c);
  if(!((MEDFileAnyTypeFieldMultiTSWithoutSDA *)_content))
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTSWithoutSDA::pushBackTimeStep : no content in this !");
  _content->pushBackTimeStep(cSafe);
  appendGlobs(*f1ts,1e-12);
}

void MEDFileAnyTypeFieldMultiTS::synchronizeNameScope()
{
  contentNotNullBase()->synchronizeNameScope();
}

int MEDFileAnyTypeFieldMultiTS::getPosOfTimeStep(int iteration, int order) const
{
  return contentNotNullBase()->getPosOfTimeStep(iteration,order);
}

int MEDFileAnyTypeFieldMultiTS::getPosGivenTime(double time, double eps) const
{
  return contentNotNullBase()->getPosGivenTime(time,eps);
}

int MEDFileAnyTypeFieldMultiTS::getNonEmptyLevels(int iteration, int order, const std::string& mname, std::vector<int>& levs) const
{
  return contentNotNullBase()->getNonEmptyLevels(iteration,order,mname,levs);
}

std::vector< std::vector<TypeOfField> > MEDFileAnyTypeFieldMultiTS::getTypesOfFieldAvailable() const
{
  return contentNotNullBase()->getTypesOfFieldAvailable();
}

std::vector< std::vector< std::pair<int,int> > > MEDFileAnyTypeFieldMultiTS::getFieldSplitedByType(int iteration, int order, const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const
{
  return contentNotNullBase()->getFieldSplitedByType(iteration,order,mname,types,typesF,pfls,locs);
}

std::string MEDFileAnyTypeFieldMultiTS::getName() const
{
  return contentNotNullBase()->getName();
}

void MEDFileAnyTypeFieldMultiTS::setName(const std::string& name)
{
  contentNotNullBase()->setName(name);
}

std::string MEDFileAnyTypeFieldMultiTS::getDtUnit() const
{
  return contentNotNullBase()->getDtUnit();
}

void MEDFileAnyTypeFieldMultiTS::setDtUnit(const std::string& dtUnit)
{
  contentNotNullBase()->setDtUnit(dtUnit);
}

void MEDFileAnyTypeFieldMultiTS::simpleRepr(int bkOffset, std::ostream& oss, int fmtsId) const
{
  contentNotNullBase()->simpleRepr(bkOffset,oss,fmtsId);
}

std::vector< std::pair<int,int> > MEDFileAnyTypeFieldMultiTS::getTimeSteps(std::vector<double>& ret1) const
{
  return contentNotNullBase()->getTimeSteps(ret1);
}

std::string MEDFileAnyTypeFieldMultiTS::getMeshName() const
{
  return contentNotNullBase()->getMeshName();
}

void MEDFileAnyTypeFieldMultiTS::setMeshName(const std::string& newMeshName)
{
  contentNotNullBase()->setMeshName(newMeshName);
}

bool MEDFileAnyTypeFieldMultiTS::changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab)
{
  return contentNotNullBase()->changeMeshNames(modifTab);
}

const std::vector<std::string>& MEDFileAnyTypeFieldMultiTS::getInfo() const
{
  return contentNotNullBase()->getInfo();
}

bool MEDFileAnyTypeFieldMultiTS::presenceOfMultiDiscPerGeoType() const
{
  return contentNotNullBase()->presenceOfMultiDiscPerGeoType();
}

void MEDFileAnyTypeFieldMultiTS::setInfo(const std::vector<std::string>& info)
{
  return contentNotNullBase()->setInfo(info);
}

int MEDFileAnyTypeFieldMultiTS::getNumberOfComponents() const
{
  const std::vector<std::string> ret=getInfo();
  return (int)ret.size();
}

void MEDFileAnyTypeFieldMultiTS::writeLL(med_idt fid) const
{
  writeGlobals(fid,*this);
  contentNotNullBase()->writeLL(fid,*this);
}

/*!
 * This method alloc the arrays and load potentially huge arrays contained in this field.
 * This method should be called when a MEDFileAnyTypeFieldMultiTS::New constructor has been with false as the last parameter.
 * This method can be also called to refresh or reinit values from a file.
 * 
 * \throw If the fileName is not set or points to a non readable MED file.
 */
void MEDFileAnyTypeFieldMultiTS::loadArrays()
{
  if(getFileName().empty())
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::loadArrays : the structure does not come from a file !");
  MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(getFileName()));
  contentNotNullBase()->loadBigArraysRecursively(fid,*contentNotNullBase());
}

/*!
 * This method behaves as MEDFileAnyTypeFieldMultiTS::loadArrays does, the first call, if \a this was built using a file without loading big arrays.
 * But once data loaded once, this method does nothing.
 * 
 * \throw If the fileName is not set or points to a non readable MED file.
 * \sa MEDFileAnyTypeFieldMultiTS::loadArrays, MEDFileAnyTypeFieldMultiTS::unloadArrays
 */
void MEDFileAnyTypeFieldMultiTS::loadArraysIfNecessary()
{
  if(!getFileName().empty())
    {
      MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(getFileName()));
      contentNotNullBase()->loadBigArraysRecursivelyIfNecessary(fid,*contentNotNullBase());
    }
}

/*!
 * This method releases potentially big data arrays and so returns to the same heap memory than status loaded with 'loadAll' parameter set to false.
 * \b WARNING, this method does release arrays even if \a this does not come from a load of a MED file.
 * So this method can lead to a loss of data. If you want to unload arrays safely call MEDFileAnyTypeFieldMultiTS::unloadArraysWithoutDataLoss instead.
 * 
 * \sa MEDFileAnyTypeFieldMultiTS::loadArrays, MEDFileAnyTypeFieldMultiTS::loadArraysIfNecessary, MEDFileAnyTypeFieldMultiTS::unloadArraysWithoutDataLoss
 */
void MEDFileAnyTypeFieldMultiTS::unloadArrays()
{
  contentNotNullBase()->unloadArrays();
}

/*!
 * This method potentially releases big data arrays if \a this is coming from a file. If \a this has been built from scratch this method will have no effect.
 * This method is the symmetrical method of MEDFileAnyTypeFieldMultiTS::loadArraysIfNecessary.
 * This method is useful to reduce \b safely amount of heap memory necessary for \a this by using MED file as database.
 * 
 * \sa MEDFileAnyTypeFieldMultiTS::loadArraysIfNecessary
 */
void MEDFileAnyTypeFieldMultiTS::unloadArraysWithoutDataLoss()
{
  if(!getFileName().empty())
    contentNotNullBase()->unloadArrays();
}

std::string MEDFileAnyTypeFieldMultiTS::simpleRepr() const
{
  std::ostringstream oss;
  contentNotNullBase()->simpleRepr(0,oss,-1);
  simpleReprGlobs(oss);
  return oss.str();
}

std::size_t MEDFileAnyTypeFieldMultiTS::getHeapMemorySizeWithoutChildren() const
{
  return MEDFileFieldGlobsReal::getHeapMemorySizeWithoutChildren();
}

std::vector<const BigMemoryObject *> MEDFileAnyTypeFieldMultiTS::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret(MEDFileFieldGlobsReal::getDirectChildrenWithNull());
  ret.push_back((const MEDFileAnyTypeFieldMultiTSWithoutSDA *)_content);
  return ret;
}

/*!
 * This method returns as MEDFileAnyTypeFieldMultiTS new instances as number of components in \a this.
 * The returned instances are deep copy of \a this except that for globals that are shared with those contained in \a this.
 * ** WARNING ** do no forget to rename the output instances to avoid to write n-times in the same MED file field !
 */
std::vector< MCAuto< MEDFileAnyTypeFieldMultiTS > > MEDFileAnyTypeFieldMultiTS::splitComponents() const
{
  const MEDFileAnyTypeFieldMultiTSWithoutSDA *content(_content);
  if(!content)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::splitComponents : no content in this ! Unable to split components !");
  std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > contentsSplit=content->splitComponents();
  std::size_t sz(contentsSplit.size());
  std::vector< MCAuto< MEDFileAnyTypeFieldMultiTS > > ret(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      ret[i]=shallowCpy();
      ret[i]->_content=contentsSplit[i];
    }
  return ret;
}

/*!
 * This method returns as MEDFileAnyTypeFieldMultiTS new instances as number of discretizations over time steps in \a this.
 * The returned instances are shallow copied of \a this included globals that are shared with those contained in \a this.
 */
std::vector< MCAuto< MEDFileAnyTypeFieldMultiTS > > MEDFileAnyTypeFieldMultiTS::splitDiscretizations() const
{
  const MEDFileAnyTypeFieldMultiTSWithoutSDA *content(_content);
  if(!content)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::splitDiscretizations : no content in this ! Unable to split discretizations !");
  std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > contentsSplit(content->splitDiscretizations());
  std::size_t sz(contentsSplit.size());
  std::vector< MCAuto< MEDFileAnyTypeFieldMultiTS > > ret(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      ret[i]=shallowCpy();
      ret[i]->_content=contentsSplit[i];
    }
  return ret;
}

/*!
 * This method returns as MEDFileAnyTypeFieldMultiTS new instances as number of sub-discretizations over time steps in \a this.
 * The returned instances are shallow copied of \a this included globals that are shared with those contained in \a this.
 */
std::vector< MCAuto< MEDFileAnyTypeFieldMultiTS > > MEDFileAnyTypeFieldMultiTS::splitMultiDiscrPerGeoTypes() const
{
  const MEDFileAnyTypeFieldMultiTSWithoutSDA *content(_content);
  if(!content)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::splitMultiDiscrPerGeoTypes : no content in this ! Unable to split discretizations !");
  std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > contentsSplit(content->splitMultiDiscrPerGeoTypes());
  std::size_t sz(contentsSplit.size());
  std::vector< MCAuto< MEDFileAnyTypeFieldMultiTS > > ret(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      ret[i]=shallowCpy();
      ret[i]->_content=contentsSplit[i];
    }
  return ret;
}

MEDFileAnyTypeFieldMultiTS *MEDFileAnyTypeFieldMultiTS::deepCopy() const
{
  MCAuto<MEDFileAnyTypeFieldMultiTS> ret=shallowCpy();
  if((const MEDFileAnyTypeFieldMultiTSWithoutSDA *)_content)
    ret->_content=_content->deepCopy();
  ret->deepCpyGlobs(*this);
  return ret.retn();
}

MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> MEDFileAnyTypeFieldMultiTS::getContent()
{
  return _content;
}

/*!
 * Returns a new MEDFileField1TS or MEDFileIntField1TS holding data of a given time step of \a this field.
 *  \param [in] iteration - the iteration number of a required time step.
 *  \param [in] order - the iteration order number of required time step.
 *  \return MEDFileField1TS * or MEDFileIntField1TS *- a new instance of MEDFileField1TS or MEDFileIntField1TS. The caller is to
 *          delete this field using decrRef() as it is no more needed.
 *  \throw If there is no required time step in \a this field.
 */
MEDFileAnyTypeField1TS *MEDFileAnyTypeFieldMultiTS::getTimeStep(int iteration, int order) const
{
  int pos=getPosOfTimeStep(iteration,order);
  return getTimeStepAtPos(pos);
}

/*!
 * Returns a new MEDFileField1TS or MEDFileIntField1TS holding data of a given time step of \a this field.
 *  \param [in] time - the time of the time step of interest.
 *  \param [in] eps - a precision used to compare time values.
 *  \return MEDFileField1TS * - a new instance of MEDFileField1TS. The caller is to
 *          delete this field using decrRef() as it is no more needed.
 *  \throw If there is no required time step in \a this field.
 */
MEDFileAnyTypeField1TS *MEDFileAnyTypeFieldMultiTS::getTimeStepGivenTime(double time, double eps) const
{
  int pos=getPosGivenTime(time,eps);
  return getTimeStepAtPos(pos);
}

/*!
 * This method groups not null items in \a vectFMTS per time step series. Two time series are considered equal if the list of their pair of integers iteration,order are equal.
 * The float64 value of time attached to the pair of integers are not considered here.
 * WARNING the returned pointers are not incremented. The caller is \b not responsible to deallocate them ! This method only reorganizes entries in \a vectFMTS.
 *
 * \param [in] vectFMTS - vector of not null fields defined on a same global data pointer.
 * \throw If there is a null pointer in \a vectFMTS.
 */
std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > MEDFileAnyTypeFieldMultiTS::SplitIntoCommonTimeSeries(const std::vector<MEDFileAnyTypeFieldMultiTS *>& vectFMTS)
{
  static const char msg[]="MEDFileAnyTypeFieldMultiTS::SplitIntoCommonTimeSeries : presence of null instance in input vector !";
  std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > ret;
  std::list<MEDFileAnyTypeFieldMultiTS *> lstFMTS(vectFMTS.begin(),vectFMTS.end());
  while(!lstFMTS.empty())
    {
      std::list<MEDFileAnyTypeFieldMultiTS *>::iterator it(lstFMTS.begin());
      MEDFileAnyTypeFieldMultiTS *curIt(*it);
      if(!curIt)
        throw INTERP_KERNEL::Exception(msg);
      std::vector< std::pair<int,int> > refIts=curIt->getIterations();
      std::vector<MEDFileAnyTypeFieldMultiTS *> elt;
      elt.push_back(curIt); it=lstFMTS.erase(it);
      while(it!=lstFMTS.end())
        {
          curIt=*it;
          if(!curIt)
            throw INTERP_KERNEL::Exception(msg);
          std::vector< std::pair<int,int> > curIts=curIt->getIterations();
          if(refIts==curIts)
            { elt.push_back(curIt); it=lstFMTS.erase(it); }
          else
            it++;
        }
      ret.push_back(elt);
    }
  return ret;
}

/*!
 * This method splits the input list \a vectFMTS considering the aspect of the geometrical support over time.
 * All returned instances in a subvector can be safely loaded, rendered along time
 * All items must be defined on the same time step ids ( see MEDFileAnyTypeFieldMultiTS::SplitIntoCommonTimeSeries method ).
 * Each item in \a vectFMTS is expected to have one and exactly one spatial discretization along time.
 * All items in \a vectFMTS must lie on the mesh (located by meshname and time step) and compatible with the input mesh \a mesh (having the same name than those in items).
 * All items in \a vectFMTS whose spatial discretization is not ON_NODES will appear once.
 * For items in \a vectFMTS that are ON_NODES it is possible to appear several times (more than once or once) in the returned vector.
 *
 * \param [in] vectFMTS - list of multi times step part all defined each on a same spatial discretization along time and pointing to a mesh whose name is equal to \c mesh->getName().
 * \param [in] mesh - the mesh shared by all items in \a vectFMTS across time.
 * \param [out] fsc - A vector having same size than returned vector. It specifies the support comporator of the corresponding vector of MEDFileAnyTypeFieldMultiTS in returned vector of vector.
 * \return - A vector of vector of objects that contains the same pointers (objects) than thoose in \a vectFMTS except that there are organized differently. So pointers included in returned vector of vector should \b not been dealt by the caller.
 *
 * \throw If an element in \a vectFMTS has not only one spatial discretization set.
 * \throw If an element in \a vectFMTS change of spatial discretization along time.
 * \throw If an element in \a vectFMTS lies on a mesh with meshname different from those in \a mesh.
 * \thorw If some elements in \a vectFMTS do not have the same times steps.
 * \throw If mesh is null.
 * \throw If an element in \a vectFMTS is null.
 * \sa MEDFileAnyTypeFieldMultiTS::AreOnSameSupportAcrossTime
 */
std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > MEDFileAnyTypeFieldMultiTS::SplitPerCommonSupport(const std::vector<MEDFileAnyTypeFieldMultiTS *>& vectFMTS, const MEDFileMesh *mesh, std::vector< MCAuto<MEDFileFastCellSupportComparator> >& fsc)
{
  static const char msg[]="MEDFileAnyTypeFieldMultiTS::SplitPerCommonSupport : presence of a null instance in the input vector !";
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::SplitPerCommonSupport : input mesh is null !");
  std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > ret;
  if(vectFMTS.empty())
    return ret;
  std::vector<MEDFileAnyTypeFieldMultiTS *>::const_iterator it(vectFMTS.begin());
  MEDFileAnyTypeFieldMultiTS *frstElt(*it);
  if(!frstElt)
    throw INTERP_KERNEL::Exception(msg);
  std::size_t i=0;
  std::vector<MEDFileAnyTypeFieldMultiTS *> vectFMTSNotNodes;
  std::vector<MEDFileAnyTypeFieldMultiTS *> vectFMTSNodes;
  for(;it!=vectFMTS.end();it++,i++)
    {
      if(!(*it))
        throw INTERP_KERNEL::Exception(msg);
      TypeOfField tof0,tof1;
      if(CheckSupportAcrossTime(frstElt,*it,mesh,tof0,tof1)>0)
        {
          if(tof1!=ON_NODES)
            vectFMTSNotNodes.push_back(*it);
          else
            vectFMTSNodes.push_back(*it);
        }
      else
        vectFMTSNotNodes.push_back(*it);
    }
  std::vector< MCAuto<MEDFileFastCellSupportComparator> > cmps;
  std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > retCell=SplitPerCommonSupportNotNodesAlg(vectFMTSNotNodes,mesh,cmps);
  ret=retCell;
  for(std::vector<MEDFileAnyTypeFieldMultiTS *>::const_iterator it2=vectFMTSNodes.begin();it2!=vectFMTSNodes.end();it2++)
    {
      i=0;
      bool isFetched(false);
      for(std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> >::const_iterator it0=retCell.begin();it0!=retCell.end();it0++,i++)
        {
          if((*it0).empty())
            throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::SplitPerCommonSupport : internal error !");
          if(cmps[i]->isCompatibleWithNodesDiscr(*it2))
            { ret[i].push_back(*it2); isFetched=true; }
        }
      if(!isFetched)
        {
          std::vector<MEDFileAnyTypeFieldMultiTS *> tmp(1,*it2);
          MCAuto<MEDFileMeshStruct> tmp2(MEDFileMeshStruct::New(mesh));
          ret.push_back(tmp); retCell.push_back(tmp); cmps.push_back(MEDFileFastCellSupportComparator::New(tmp2,*it2));
        }
    }
  fsc=cmps;
  return ret;
}

/*!
 * WARNING no check here. The caller must be sure that all items in vectFMTS are coherent each other in time steps, only one same spatial discretization and not ON_NODES.
 * \param [out] cmps - same size than the returned vector.
 */
std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > MEDFileAnyTypeFieldMultiTS::SplitPerCommonSupportNotNodesAlg(const std::vector<MEDFileAnyTypeFieldMultiTS *>& vectFMTS, const MEDFileMesh *mesh, std::vector< MCAuto<MEDFileFastCellSupportComparator> >& cmps)
{
  std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > ret;
  std::list<MEDFileAnyTypeFieldMultiTS *> lstFMTS(vectFMTS.begin(),vectFMTS.end());
  while(!lstFMTS.empty())
    {
      std::list<MEDFileAnyTypeFieldMultiTS *>::iterator it(lstFMTS.begin());
      MEDFileAnyTypeFieldMultiTS *ref(*it);
      std::vector<MEDFileAnyTypeFieldMultiTS *> elt;
      elt.push_back(ref); it=lstFMTS.erase(it);
      MCAuto<MEDFileMeshStruct> mst(MEDFileMeshStruct::New(mesh));
      MCAuto<MEDFileFastCellSupportComparator> cmp(MEDFileFastCellSupportComparator::New(mst,ref));
      while(it!=lstFMTS.end())
        {
          MEDFileAnyTypeFieldMultiTS *curIt(*it);
          if(cmp->isEqual(curIt))
            { elt.push_back(curIt); it=lstFMTS.erase(it); }
          else
            it++;
        }
      ret.push_back(elt); cmps.push_back(cmp);
    }
  return ret;
}

/*!
 * This method scan the two main structs along time of \a f0 and \a f1 to see if there are all lying on the same mesh along time than those in \a mesh.
 * \a f0 and \a f1 must be defined each only on a same spatial discretization even if this can be different each other.
 *
 * \throw If \a f0 or \a f1 has not only one spatial discretization set.
 * \throw If \a f0 or \a f1 change of spatial discretization along time.
 * \throw If \a f0 or \a f1 on a mesh with meshname different from those in \a mesh.
 * \thorw If \a f0 and \a f1 do not have the same times steps.
 * \throw If mesh is null.
 * \throw If \a f0 or \a f1 is null.
 * \sa MEDFileAnyTypeFieldMultiTS::SplitPerCommonSupport
 */
int MEDFileAnyTypeFieldMultiTS::CheckSupportAcrossTime(MEDFileAnyTypeFieldMultiTS *f0, MEDFileAnyTypeFieldMultiTS *f1, const MEDFileMesh *mesh, TypeOfField& tof0, TypeOfField& tof1)
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::CheckSupportAcrossTime : input mesh is null !");
  if(!f0 || !f1)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::CheckSupportAcrossTime : presence of null instance in fields over time !");
  if(f0->getMeshName()!=mesh->getName())
    {
      std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTS::CheckSupportAcrossTime : first field points to mesh \""<< f0->getMeshName() << "\" and input mesh to compare has name \"" << mesh->getName() << "\" !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  if(f1->getMeshName()!=mesh->getName())
    {
      std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTS::CheckSupportAcrossTime : second field points to mesh \""<< f1->getMeshName() << "\" and input mesh to compare has name \"" << mesh->getName() << "\" !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  int nts=f0->getNumberOfTS();
  if(nts!=f1->getNumberOfTS())
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::CheckSupportAcrossTime : number of time steps are not the same !");
  if(nts==0)
    return nts;
  for(int i=0;i<nts;i++)
    {
      MCAuto<MEDFileAnyTypeField1TS> f0cur=f0->getTimeStepAtPos(i);
      MCAuto<MEDFileAnyTypeField1TS> f1cur=f1->getTimeStepAtPos(i);
      std::vector<TypeOfField> tofs0(f0cur->getTypesOfFieldAvailable()),tofs1(f1cur->getTypesOfFieldAvailable());
      if(tofs0.size()!=1 || tofs1.size()!=1)
        throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::CheckSupportAcrossTime : All time steps must be defined on only one spatial discretization !");
      if(i!=0)
        {
          if(tof0!=tofs0[0] || tof1!=tofs1[0])
            throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::CheckSupportAcrossTime : Across times steps MEDFileAnyTypeFieldMultiTS instances have to keep the same unique spatial discretization !");
        }
      else
        { tof0=tofs0[0]; tof1=tofs1[0]; }
      if(f0cur->getMeshIteration()!=mesh->getIteration() || f0cur->getMeshOrder()!=mesh->getOrder())
        {
          std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTS::CheckSupportAcrossTime : first field points to mesh time step (" << f0cur->getMeshIteration() << ","<< f0cur->getMeshOrder() << ") whereas input mesh points to time step (" << mesh->getIteration() << "," << mesh->getOrder() << ") !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
      if(f1cur->getMeshIteration()!=mesh->getIteration() || f1cur->getMeshOrder()!=mesh->getOrder())
        {
          std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTS::CheckSupportAcrossTime : second field points to mesh time step (" << f1cur->getMeshIteration() << ","<< f1cur->getMeshOrder() << ") whereas input mesh points to time step (" << mesh->getIteration() << "," << mesh->getOrder() << ") !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
      if(f0cur->getIteration()!=f1cur->getIteration() || f0cur->getOrder()!=f1cur->getOrder())
        {
          std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTS::CheckSupportAcrossTime : all the time steps must be the same ! it is not the case (" << f0cur->getIteration() << "," << f0cur->getOrder() << ")!=(" << f1cur->getIteration() << "," << f1cur->getOrder() << ") !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
    }
  return nts;
}

template<class T>
MCAuto<MEDFileAnyTypeField1TS> AggregateHelperF1TS(const std::vector< typename MLFieldTraits<T>::F1TSType const *>& f1tss, const std::vector< std::vector< std::pair<int,int> > >& dts)
{
  MCAuto< typename MLFieldTraits<T>::F1TSType > ret(MLFieldTraits<T>::F1TSType::New());
  if(f1tss.empty())
    throw INTERP_KERNEL::Exception("AggregateHelperF1TS : empty vector !");
  std::size_t sz(f1tss.size()),i(0);
  std::vector< typename MLFieldTraits<T>::F1TSWSDAType const *> f1tsw(sz);
  for(typename std::vector< typename MLFieldTraits<T>::F1TSType const *>::const_iterator it=f1tss.begin();it!=f1tss.end();it++,i++)
    {
      typename MLFieldTraits<T>::F1TSType const *elt(*it);
      if(!elt)
        throw INTERP_KERNEL::Exception("AggregateHelperF1TS : presence of a null pointer !");
      f1tsw[i]=dynamic_cast<typename MLFieldTraits<T>::F1TSWSDAType const *>(elt->contentNotNullBase());
    }
  typename MLFieldTraits<T>::F1TSWSDAType *retc(dynamic_cast<typename MLFieldTraits<T>::F1TSWSDAType *>(ret->contentNotNullBase()));
  if(!retc)
    throw INTERP_KERNEL::Exception("AggregateHelperF1TS : internal error 1 !");
  retc->aggregate(f1tsw,dts);
  ret->setDtUnit(f1tss[0]->getDtUnit());
  return DynamicCast<typename MLFieldTraits<T>::F1TSType , MEDFileAnyTypeField1TS>(ret);
}

template<class T>
MCAuto< MEDFileAnyTypeFieldMultiTS > AggregateHelperFMTS(const std::vector< typename MLFieldTraits<T>::FMTSType const *>& fmtss, const std::vector< std::vector< std::pair<int,int> > >& dts)
{
  MCAuto< typename MLFieldTraits<T>::FMTSType > ret(MLFieldTraits<T>::FMTSType::New());
  if(fmtss.empty())
    throw INTERP_KERNEL::Exception("AggregateHelperFMTS : empty vector !");
  std::size_t sz(fmtss.size());
  for(typename std::vector< typename MLFieldTraits<T>::FMTSType const *>::const_iterator it=fmtss.begin();it!=fmtss.end();it++)
    {
      typename MLFieldTraits<T>::FMTSType const *elt(*it);
      if(!elt)
        throw INTERP_KERNEL::Exception("AggregateHelperFMTS : presence of null pointer !");
    }
  int nbTS(fmtss[0]->getNumberOfTS());
  for(typename std::vector< typename MLFieldTraits<T>::FMTSType const *>::const_iterator it=fmtss.begin();it!=fmtss.end();it++)
    if((*it)->getNumberOfTS()!=nbTS)
      throw INTERP_KERNEL::Exception("AggregateHelperFMTS : all fields must have the same number of TS !");
  for(int iterTS=0;iterTS<nbTS;iterTS++)
    {
      std::size_t i(0);
      std::vector< typename MLFieldTraits<T>::F1TSType const *> f1tss(sz);
      std::vector< MCAuto<typename MLFieldTraits<T>::F1TSType> > f1tss2(sz);
      for(typename std::vector< typename MLFieldTraits<T>::FMTSType const *>::const_iterator it=fmtss.begin();it!=fmtss.end();it++,i++)
        { f1tss2[i]=(*it)->getTimeStepAtPos(iterTS); f1tss[i]=f1tss2[i]; }
      MCAuto<MEDFileAnyTypeField1TS> f1ts(AggregateHelperF1TS<T>(f1tss,dts));
      ret->pushBackTimeStep(f1ts);
      ret->setDtUnit(f1ts->getDtUnit());
    }
  return DynamicCast<typename MLFieldTraits<T>::FMTSType , MEDFileAnyTypeFieldMultiTS>(ret);
}

/*!
 * \a dts and \a ftmss are expected to have same size.
 */
MCAuto<MEDFileAnyTypeFieldMultiTS> MEDFileAnyTypeFieldMultiTS::Aggregate(const std::vector<const MEDFileAnyTypeFieldMultiTS *>& fmtss, const std::vector< std::vector< std::pair<int,int> > >& dts)
{
  if(fmtss.empty())
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::Aggregate : input vector is empty !");
  std::size_t sz(fmtss.size());
  std::vector<const MEDFileFieldMultiTS *> fmtss1;
  std::vector<const MEDFileIntFieldMultiTS *> fmtss2;
  for(std::vector<const MEDFileAnyTypeFieldMultiTS *>::const_iterator it=fmtss.begin();it!=fmtss.end();it++)
    {
      if(!(*it))
        throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::Aggregate : presence of null instance in input vector !");
      const MEDFileFieldMultiTS *elt1(dynamic_cast<const MEDFileFieldMultiTS *>(*it));
      if(elt1)
        {
          fmtss1.push_back(elt1);
          continue;
        }
      const MEDFileIntFieldMultiTS *elt2(dynamic_cast<const MEDFileIntFieldMultiTS *>(*it));
      if(elt2)
        {
          fmtss2.push_back(elt2);
          continue;
        }
      throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::Aggregate : not recognized type !");
    }
  if(fmtss1.size()!=sz && fmtss2.size()!=sz)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::Aggregate : type of data is not homogeneous !");
  if(fmtss1.size()==sz)
    return AggregateHelperFMTS<double>(fmtss1,dts);
  if(fmtss2.size()!=sz)
    return AggregateHelperFMTS<int>(fmtss2,dts);
  throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::Aggregate : not implemented yet !");
}

MEDFileAnyTypeFieldMultiTSIterator *MEDFileAnyTypeFieldMultiTS::iterator()
{
  return new MEDFileAnyTypeFieldMultiTSIterator(this);
}

//= MEDFileFieldMultiTS

MEDFileAnyTypeFieldMultiTS *MEDFileFieldMultiTS::shallowCpy() const
{
  return new MEDFileFieldMultiTS(*this);
}

/*!
 * This method performs a copy with datatype modification ( float64->int32 ) of \a this. The globals information are copied
 * following the given input policy.
 *
 * \param [in] isDeepCpyGlobs - a boolean that indicates the behaviour concerning globals (profiles and localizations)
 *                            By default (true) the globals are deeply copied.
 * \return MEDFileIntFieldMultiTS * - a new object that is the result of the conversion of \a this to int32 field.
 */
MEDFileIntFieldMultiTS *MEDFileFieldMultiTS::convertToInt(bool isDeepCpyGlobs) const
{
  MCAuto<MEDFileIntFieldMultiTS> ret;
  const MEDFileAnyTypeFieldMultiTSWithoutSDA *content(_content);
  if(content)
    {
      const MEDFileFieldMultiTSWithoutSDA *contc=dynamic_cast<const MEDFileFieldMultiTSWithoutSDA *>(content);
      if(!contc)
        throw INTERP_KERNEL::Exception("MEDFileFieldMultiTS::convertToInt : the content inside this is not FLOAT64 ! This is incoherent !");
      MCAuto<MEDFileIntFieldMultiTSWithoutSDA> newc(contc->convertToInt());
      ret=static_cast<MEDFileIntFieldMultiTS *>(MEDFileAnyTypeFieldMultiTS::BuildNewInstanceFromContent((MEDFileIntFieldMultiTSWithoutSDA *)newc));
    }
  else
    ret=MEDFileIntFieldMultiTS::New();
  if(isDeepCpyGlobs)
    ret->deepCpyGlobs(*this);
  else
    ret->shallowCpyGlobs(*this);
  return ret.retn();
}

MEDFileFieldMultiTS::MEDFileFieldMultiTS(med_idt fid, bool loadAll, const MEDFileMeshes *ms)
try:MEDFileTemplateFieldMultiTS<double>(fid,loadAll,ms)
{
}
catch(INTERP_KERNEL::Exception& e)
{ throw e; }

MEDFileFieldMultiTS::MEDFileFieldMultiTS(med_idt fid, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities)
try:MEDFileTemplateFieldMultiTS<double>(fid,fieldName,loadAll,ms,entities)
{
}
catch(INTERP_KERNEL::Exception& e)
{ throw e; }

MEDFileFieldMultiTS::MEDFileFieldMultiTS(const MEDFileFieldMultiTSWithoutSDA& other, bool shallowCopyOfContent):MEDFileTemplateFieldMultiTS<double>(other,shallowCopyOfContent)
{
}

std::vector< std::vector<DataArrayDouble *> > MEDFileFieldMultiTS::getFieldSplitedByType2(int iteration, int order, const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const
{
  return contentNotNull()->getFieldSplitedByType2(iteration,order,mname,types,typesF,pfls,locs);
}

//= MEDFileAnyTypeFieldMultiTSIterator

MEDFileAnyTypeFieldMultiTSIterator::MEDFileAnyTypeFieldMultiTSIterator(MEDFileAnyTypeFieldMultiTS *fmts):_fmts(fmts),_iter_id(0),_nb_iter(0)
{
  if(fmts)
    {
      fmts->incrRef();
      _nb_iter=fmts->getNumberOfTS();
    }
}

MEDFileAnyTypeFieldMultiTSIterator::~MEDFileAnyTypeFieldMultiTSIterator() 
{
}

MEDFileAnyTypeField1TS *MEDFileAnyTypeFieldMultiTSIterator::nextt()
{
  if(_iter_id<_nb_iter)
    {
      MEDFileAnyTypeFieldMultiTS *fmts(_fmts);
      if(fmts)
        return fmts->getTimeStepAtPos(_iter_id++);
      else
        return 0;
    }
  else
    return 0;
}

//= MEDFileIntFieldMultiTS
