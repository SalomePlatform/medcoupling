//  Copyright (C) 2007-2011  CEA/DEN, EDF R&D
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

#include "MEDFileField.hxx"
#include "MEDLoaderBase.hxx"
#include "MEDFileUtilities.hxx"

#include "InterpKernelAutoPtr.hxx"

#include <algorithm>

extern med_geometrie_element typmai[MED_NBR_GEOMETRIE_MAILLE+2];
extern INTERP_KERNEL::NormalizedCellType typmai2[MED_NBR_GEOMETRIE_MAILLE+2];
extern med_geometrie_element typmainoeud[1];
extern med_geometrie_element typmai3[32];

using namespace ParaMEDMEM;

MEDFileFieldPerMeshPerTypePerDisc *MEDFileFieldPerMeshPerTypePerDisc::New(MEDFileFieldPerMeshPerType *fath, med_idt fid) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldPerMeshPerTypePerDisc(fath,fid);
}

MEDFileFieldPerMeshPerTypePerDisc::MEDFileFieldPerMeshPerTypePerDisc(MEDFileFieldPerMeshPerType *fath, med_idt fid) throw(INTERP_KERNEL::Exception)
try:_father(fath)
{
  INTERP_KERNEL::AutoPtr<char> locname=MEDLoaderBase::buildEmptyString(MED_TAILLE_NOM);
  INTERP_KERNEL::AutoPtr<char> pflname=MEDLoaderBase::buildEmptyString(MED_TAILLE_NOM);
  std::string fieldName=getName();
  std::string meshName=getMeshName();
  int iteration=getIteration();
  int order=getOrder();
  const std::vector<std::string>& infos=getInfos();
  TypeOfField type=getType();
  INTERP_KERNEL::NormalizedCellType geoType=getGeoType();
  int nval=MEDnVal(fid,(char *)fieldName.c_str(),type==ON_CELLS?MED_MAILLE:MED_NOEUD,type==ON_CELLS?typmai3[(int)geoType]:MED_NONE,
                   iteration,order,(char *)meshName.c_str(),MED_COMPACT);
  _arr=DataArrayDouble::New();
  _arr->alloc(nval,infos.size());
  MEDchampLire(fid,(char *)getMeshName().c_str(),(char *)fieldName.c_str(),(unsigned char*)_arr->getPointer(),MED_FULL_INTERLACE,MED_ALL,locname,
               pflname,MED_COMPACT,type==ON_CELLS?MED_MAILLE:MED_NOEUD,
               type==ON_CELLS?typmai3[(int)geoType]:MED_NONE,
               iteration,order);
  _profile=pflname;
  _localization=locname;
}
catch(INTERP_KERNEL::Exception& e)
{
  throw e;
}

const MEDFileFieldPerMeshPerType *MEDFileFieldPerMeshPerTypePerDisc::getFather() const
{
  return _father;
}

int MEDFileFieldPerMeshPerTypePerDisc::getIteration() const
{
  return _father->getIteration();
}

int MEDFileFieldPerMeshPerTypePerDisc::getOrder() const
{
  return _father->getOrder();
}

std::string MEDFileFieldPerMeshPerTypePerDisc::getName() const
{
  return _father->getName();
}

TypeOfField MEDFileFieldPerMeshPerTypePerDisc::getType() const
{
  return _father->getType();
}

INTERP_KERNEL::NormalizedCellType MEDFileFieldPerMeshPerTypePerDisc::getGeoType() const
{
  return _father->getGeoType();
}

std::string MEDFileFieldPerMeshPerTypePerDisc::getMeshName() const
{
  return _father->getMeshName();
}

int MEDFileFieldPerMeshPerTypePerDisc::getNumberOfComponents() const
{
  return _father->getNumberOfComponents();
}

const std::vector<std::string>& MEDFileFieldPerMeshPerTypePerDisc::getInfos() const
{
  return _father->getInfos();
}

std::string MEDFileFieldPerMeshPerTypePerDisc::getProfile() const
{
  return _profile;
}

std::string MEDFileFieldPerMeshPerTypePerDisc::getLocalization() const
{
  return _localization;
}

MEDFileFieldPerMeshPerType *MEDFileFieldPerMeshPerType::New(MEDFileFieldPerMesh *fath, TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldPerMeshPerType(fath,type,geoType);
}

const MEDFileFieldPerMesh *MEDFileFieldPerMeshPerType::getFather() const
{
  return _father;
}

int MEDFileFieldPerMeshPerType::getIteration() const
{
  return _father->getIteration();
}

int MEDFileFieldPerMeshPerType::getOrder() const
{
  return _father->getOrder();
}

std::string MEDFileFieldPerMeshPerType::getName() const
{
  return _father->getName();
}

std::string MEDFileFieldPerMeshPerType::getMeshName() const
{
  return _father->getMeshName();
}

TypeOfField MEDFileFieldPerMeshPerType::getType() const
{
  return _type;
}

INTERP_KERNEL::NormalizedCellType MEDFileFieldPerMeshPerType::getGeoType() const
{
  return _geo_type;
}


int MEDFileFieldPerMeshPerType::getNumberOfComponents() const
{
  return _father->getNumberOfComponents();
}

const std::vector<std::string>& MEDFileFieldPerMeshPerType::getInfos() const
{
  return _father->getInfos();
}

std::vector<std::string> MEDFileFieldPerMeshPerType::getPflsReallyUsed() const
{
  std::vector<std::string> ret;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it1=_field_pm_pt_pd.begin();it1!=_field_pm_pt_pd.end();it1++)
    {
      std::string tmp=(*it1)->getProfile();
      if(!tmp.empty())
        ret.push_back(tmp);
    }
  return ret;
}

std::vector<std::string> MEDFileFieldPerMeshPerType::getLocsReallyUsed() const
{
  std::vector<std::string> ret;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it1=_field_pm_pt_pd.begin();it1!=_field_pm_pt_pd.end();it1++)
    {
      std::string tmp=(*it1)->getLocalization();
      if(!tmp.empty())
        ret.push_back(tmp);
    }
  return ret;
}

MEDFileFieldPerMeshPerType::MEDFileFieldPerMeshPerType(MEDFileFieldPerMesh *fath, TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception):_father(fath),_type(type),_geo_type(geoType)
{
}

void MEDFileFieldPerMeshPerType::finishLoading(med_idt fid) throw(INTERP_KERNEL::Exception)
{
  //for MED3 porting this limitation will be killed !
  _field_pm_pt_pd.resize(1);
  _field_pm_pt_pd[0]=MEDFileFieldPerMeshPerTypePerDisc::New(this,fid);
}

MEDFileFieldPerMesh *MEDFileFieldPerMesh::New(MEDFileField1TSWithoutDAS *fath, const char *meshName, double time)
{
  return new MEDFileFieldPerMesh(fath,meshName,time);
}

void MEDFileFieldPerMesh::pushBack(TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType)
{
  _types.push_back(type);
  _geo_types.push_back(geoType);
}

void MEDFileFieldPerMesh::finishLoading(med_idt fid) throw(INTERP_KERNEL::Exception)
{
  int sz=_types.size();
  _field_pm_pt.resize(sz);
  for(int i=0;i<sz;i++)
    {
      _field_pm_pt[i]=MEDFileFieldPerMeshPerType::New(this,_types[i],_geo_types[i]);
      _field_pm_pt[i]->finishLoading(fid);
    }
}

double MEDFileFieldPerMesh::getTime() const
{
  return _time;
}

int MEDFileFieldPerMesh::getIteration() const
{
  return _father->getIteration();
}

int MEDFileFieldPerMesh::getOrder() const
{
  return _father->getOrder();
}

std::string MEDFileFieldPerMesh::getName() const
{
  return _father->getName();
}

std::string MEDFileFieldPerMesh::getMeshName() const
{
  return _mesh_name;
}

int MEDFileFieldPerMesh::getNumberOfComponents() const
{
  return _father->getNumberOfComponents();
}

const std::vector<std::string>& MEDFileFieldPerMesh::getInfos() const
{
  return _father->getInfos();
}

std::vector<std::string> MEDFileFieldPerMesh::getPflsReallyUsed() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getPflsReallyUsed();
      for(std::vector<std::string>::const_iterator it2=tmp.begin();it2!=tmp.end();it2++)
        if(ret2.find(*it2)==ret2.end())
          {
            ret.push_back(*it2);
            ret2.insert(*it2);
          }
    }
  return ret;
}

std::vector<std::string> MEDFileFieldPerMesh::getLocsReallyUsed() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getLocsReallyUsed();
      for(std::vector<std::string>::const_iterator it2=tmp.begin();it2!=tmp.end();it2++)
        if(ret2.find(*it2)==ret2.end())
          {
            ret.push_back(*it2);
            ret2.insert(*it2);
          }
    }
  return ret;
}

MEDFileFieldPerMesh::MEDFileFieldPerMesh(MEDFileField1TSWithoutDAS *fath, const char *meshName, double time):_father(fath),_mesh_name(meshName),_time(time)
{
}

void MEDFieldFieldGlobs::loadProfileInFile(med_idt fid, int id, const char *pflName, int lgth) throw(INTERP_KERNEL::Exception)
{
  _pfls[id]=DataArrayInt::New();
  _pfls[id]->alloc(lgth,1);
  MEDprofilLire(fid,_pfls[id]->getPointer(),(char *)pflName);
}

void MEDFieldFieldGlobs::loadProfileInFile(med_idt fid, int i)
{
  INTERP_KERNEL::AutoPtr<char> pflName=MEDLoaderBase::buildEmptyString(MED_TAILLE_NOM);
  int sz;
  MEDprofilInfo(fid,i,pflName,&sz);
  _pfls[i]=DataArrayInt::New();
  _pfls[i]->alloc(sz,1);
  MEDprofilLire(fid,_pfls[i]->getPointer(),(char *)pflName);
}

std::vector<std::string> MEDFieldFieldGlobs::getPfls() const
{
  int sz=_pfls.size();
  std::vector<std::string> ret(sz);
  for(int i=0;i<sz;i++)
    ret[i]=_pfls[i]->getName();
  return ret;
}

std::vector<std::string> MEDFieldFieldGlobs::getLocs() const
{
  int sz=_locs.size();
  std::vector<std::string> ret(sz);
  for(int i=0;i<sz;i++)
    ret[i]=_locs[i]->getName();
  return ret;
}

MEDFileField1TSWithoutDAS *MEDFileField1TSWithoutDAS::New(const char *fieldName, int iteration, int order, const std::vector<std::string>& infos)
{
  return new MEDFileField1TSWithoutDAS(fieldName,iteration,order,infos);
}

void MEDFileField1TSWithoutDAS::pushBack(TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType, double time, const char *meshName) const
{
  _types.push_back(type);
  _geo_types.push_back(geoType);
  _times.push_back(time);
  _meshes.push_back(meshName);
}

void MEDFileField1TSWithoutDAS::finishLoading(med_idt fid) throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> meshesOnlyOnce;
  int nbOfTurn=_meshes.size();
  for(int i=0;i<nbOfTurn;i++)
    {
      std::vector<std::string>::iterator it=std::find(meshesOnlyOnce.begin(),meshesOnlyOnce.end(),_meshes[i]);
      if(it==meshesOnlyOnce.end())
        {
          _field_per_mesh.push_back(MEDFileFieldPerMesh::New(this,_meshes[i].c_str(),_times[i]));
          meshesOnlyOnce.push_back(_meshes[i]);
          _field_per_mesh.back()->pushBack(_types[i],_geo_types[i]);
        }
      else
        {
          int w=std::distance(meshesOnlyOnce.begin(),it);
          _field_per_mesh[w]->pushBack(_types[i],_geo_types[i]);
        }
    }
  int nbOfMeshes=_field_per_mesh.size();
  for(int i=0;i<nbOfMeshes;i++)
    _field_per_mesh[i]->finishLoading(fid);
}

std::vector<std::string> MEDFileField1TSWithoutDAS::getPflsReallyUsed() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getPflsReallyUsed();
      for(std::vector<std::string>::const_iterator it2=tmp.begin();it2!=tmp.end();it2++)
        if(ret2.find(*it2)==ret2.end())
          {
            ret.push_back(*it2);
            ret2.insert(*it2);
          }
    }
  return ret;
}

std::vector<std::string> MEDFileField1TSWithoutDAS::getLocsReallyUsed() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getLocsReallyUsed();
      for(std::vector<std::string>::const_iterator it2=tmp.begin();it2!=tmp.end();it2++)
        if(ret2.find(*it2)==ret2.end())
          {
            ret.push_back(*it2);
            ret2.insert(*it2);
          }
    }
  return ret;
}

MEDFileField1TSWithoutDAS::MEDFileField1TSWithoutDAS(const char *fieldName, int iteration, int order, const std::vector<std::string>& infos):_name(fieldName),_infos(infos),_iteration(iteration),_order(order)
{
}

MEDFileField1TS *MEDFileField1TS::New(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileField1TS(fileName,fieldName,iteration,order);
}

MEDFileField1TS::MEDFileField1TS(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
try:MEDFileField1TSWithoutDAS(fieldName,iteration,order,std::vector<std::string>())
{
  MEDFileUtilities::CheckFileForRead(fileName);
  med_idt fid=MEDouvrir((char *)fileName,MED_LECTURE);
  int nbFields=MEDnChamp(fid,0);
  med_type_champ typcha;
  bool found=false;
  std::vector<std::string> fns(nbFields);
  for(int i=0;i<nbFields && !found;i++)
    {
      int ncomp=MEDnChamp(fid,i+1);
      INTERP_KERNEL::AutoPtr<char> comp=new char[ncomp*MED_TAILLE_PNOM+1];
      INTERP_KERNEL::AutoPtr<char> unit=new char[ncomp*MED_TAILLE_PNOM+1];
      INTERP_KERNEL::AutoPtr<char> nomcha=MEDLoaderBase::buildEmptyString(MED_TAILLE_NOM);
      MEDchampInfo(fid,i+1,nomcha,&typcha,comp,unit,ncomp);
      std::string tmp(nomcha);
      fns[i]=tmp;
      found=(tmp==fieldName);
      if(found)
        {
          _infos.resize(ncomp);
          for(int j=0;j<ncomp;j++)
            _infos[j]=MEDLoaderBase::buildUnionUnit((char *)comp+i*MED_TAILLE_PNOM,MED_TAILLE_PNOM,(char *)unit+i*MED_TAILLE_PNOM,MED_TAILLE_PNOM);
        }
    }
  if(!found)
    {
      std::ostringstream oss; oss << "No such field '" << fieldName << "' in file '" << fileName << "' ! Available fields are : ";
      std::copy(fns.begin(),fns.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  finishLoading(fid);
  //
  std::vector<std::string> profiles=getPflsReallyUsed();
  int sz=profiles.size();
  _pfls.resize(sz);
  for(int i=0;i<sz;i++)
    loadProfileInFile(fid,i,profiles[i].c_str(),37);//tony
  MEDfermer(fid);
}
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

MEDFileFieldMultiTSWithoutDAS *MEDFileFieldMultiTSWithoutDAS::New(med_idt fid, const char *fieldName, int id, const std::vector<std::string>& infos) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldMultiTSWithoutDAS(fid,fieldName,id,infos);
}

MEDFileFieldMultiTSWithoutDAS::MEDFileFieldMultiTSWithoutDAS(const char *fieldName):_name(fieldName)
{
}

MEDFileFieldMultiTSWithoutDAS::MEDFileFieldMultiTSWithoutDAS(med_idt fid, const char *fieldName, int id,  const std::vector<std::string>& infos) throw(INTERP_KERNEL::Exception)
try:_name(fieldName),_infos(infos)
{
  finishLoading(fid);
}
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

void MEDFileFieldMultiTSWithoutDAS::finishLoading(med_idt fid) throw(INTERP_KERNEL::Exception)
{
  for(int i=0;i<MED_NBR_GEOMETRIE_MAILLE+2;i++)
    {
      int nbPdt=MEDnPasdetemps(fid,(char *)_name.c_str(),MED_MAILLE,typmai[i]);
      for(int j=0;j<nbPdt;j++)
        appendTimeStepEntry(fid,MED_MAILLE,i,j);
    }
  int nbPdt=MEDnPasdetemps(fid,(char *)_name.c_str(),MED_NOEUD,MED_NONE);
  for(int j=0;j<nbPdt;j++)
    appendTimeStepEntry(fid,MED_NOEUD,0,j);
  int nbOfTimeSteps=_time_steps.size();
  for(int i=0;i<nbOfTimeSteps;i++)
    _time_steps[i]->finishLoading(fid);
}

void MEDFileFieldMultiTSWithoutDAS::appendTimeStepEntry(med_idt fid, med_entite_maillage entity, int i, int j) throw(INTERP_KERNEL::Exception)
{
  std::vector< std::pair<int,int> > ts;
  med_int ngauss=0;
  med_int numdt=0,numo=0,nbrefmaa;
  med_float dt=0.0;
  med_booleen local;
  INTERP_KERNEL::AutoPtr<char> maa_ass=MEDLoaderBase::buildEmptyString(MED_TAILLE_NOM);
  INTERP_KERNEL::AutoPtr<char> dt_unit=MEDLoaderBase::buildEmptyString(MED_TAILLE_PNOM);
  MEDpasdetempsInfo(fid,(char *)_name.c_str(),MED_MAILLE,typmai[i],j+1,&ngauss,&numdt,&numo,dt_unit,&dt,maa_ass,&local,&nbrefmaa);
  std::pair<int,int> p(numdt,numo);
  std::vector< std::pair<int,int> >::iterator where=std::find(ts.begin(),ts.end(),p);
  if(where==ts.end())
    {
      ts.push_back(p);
      _time_steps.push_back(MEDFileField1TSWithoutDAS::New(_name.c_str(),numdt,numo,_infos));
      _time_steps.back()->pushBack(entity==MED_MAILLE?ON_CELLS:ON_NODES,entity==MED_MAILLE?typmai2[i]:INTERP_KERNEL::NORM_ERROR,dt,maa_ass);
    }
  else
    {
      int w=std::distance(ts.begin(),where);
      _time_steps[w]->pushBack(entity==MED_MAILLE?ON_CELLS:ON_NODES,entity==MED_MAILLE?typmai2[i]:INTERP_KERNEL::NORM_ERROR,dt,maa_ass);
    }
}

std::vector<std::string> MEDFileFieldMultiTSWithoutDAS::getPflsReallyUsed() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileField1TSWithoutDAS > >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getPflsReallyUsed();
      for(std::vector<std::string>::const_iterator it2=tmp.begin();it2!=tmp.end();it2++)
        if(ret2.find(*it2)==ret2.end())
          {
            ret.push_back(*it2);
            ret2.insert(*it2);
          }
    }
  return ret;
}

std::vector<std::string> MEDFileFieldMultiTSWithoutDAS::getLocsReallyUsed() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileField1TSWithoutDAS > >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getLocsReallyUsed();
      for(std::vector<std::string>::const_iterator it2=tmp.begin();it2!=tmp.end();it2++)
        if(ret2.find(*it2)==ret2.end())
          {
            ret.push_back(*it2);
            ret2.insert(*it2);
          }
    }
  return ret;
}

MEDFileFieldMultiTS *MEDFileFieldMultiTS::New(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldMultiTS(fileName,fieldName);
}

MEDFileFieldMultiTS::MEDFileFieldMultiTS(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
try:MEDFileFieldMultiTSWithoutDAS(fieldName)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  med_idt fid=MEDouvrir((char *)fileName,MED_LECTURE);
  int nbFields=MEDnChamp(fid,0);
  med_type_champ typcha;
  bool found=false;
  std::vector<std::string> fns(nbFields);
  for(int i=0;i<nbFields && !found;i++)
    {
      int ncomp=MEDnChamp(fid,i+1);
      INTERP_KERNEL::AutoPtr<char> comp=new char[ncomp*MED_TAILLE_PNOM+1];
      INTERP_KERNEL::AutoPtr<char> unit=new char[ncomp*MED_TAILLE_PNOM+1];
      INTERP_KERNEL::AutoPtr<char> nomcha=MEDLoaderBase::buildEmptyString(MED_TAILLE_NOM);
      MEDchampInfo(fid,i+1,nomcha,&typcha,comp,unit,ncomp);
      std::string tmp(nomcha);
      fns[i]=tmp;
      found=(tmp==fieldName);
      if(found)
        {
          _infos.resize(ncomp);
          for(int j=0;j<ncomp;j++)
            _infos[j]=MEDLoaderBase::buildUnionUnit((char *)comp+i*MED_TAILLE_PNOM,MED_TAILLE_PNOM,(char *)unit+i*MED_TAILLE_PNOM,MED_TAILLE_PNOM);
        }
    }
  if(!found)
    {
      std::ostringstream oss; oss << "No such field '" << fieldName << "' in file '" << fileName << "' ! Available fields are : ";
      std::copy(fns.begin(),fns.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  finishLoading(fid);
  MEDfermer(fid);
}
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

MEDFileFields *MEDFileFields::New(const char *fileName) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFields(fileName);
}

int MEDFileFields::getNumberOfFields() const
{
  return _fields.size();
}

MEDFileFields::MEDFileFields(const char *fileName) throw(INTERP_KERNEL::Exception)
try
  {
    MEDFileUtilities::CheckFileForRead(fileName);
    med_idt fid=MEDouvrir((char *)fileName,MED_LECTURE);
    int nbFields=MEDnChamp(fid,0);
    _fields.resize(nbFields);
    med_type_champ typcha;
    for(int i=0;i<nbFields;i++)
      {
        int ncomp=MEDnChamp(fid,i+1);
        INTERP_KERNEL::AutoPtr<char> comp=new char[ncomp*MED_TAILLE_PNOM+1];
        INTERP_KERNEL::AutoPtr<char> unit=new char[ncomp*MED_TAILLE_PNOM+1];
        INTERP_KERNEL::AutoPtr<char> nomcha=MEDLoaderBase::buildEmptyString(MED_TAILLE_NOM);
        MEDchampInfo(fid,i+1,nomcha,&typcha,comp,unit,ncomp);
        std::vector<std::string> infos(ncomp);
        for(int j=0;j<ncomp;j++)
          infos[j]=MEDLoaderBase::buildUnionUnit((char *)comp+i*MED_TAILLE_PNOM,MED_TAILLE_PNOM,(char *)unit+i*MED_TAILLE_PNOM,MED_TAILLE_PNOM);
        _fields[i]=MEDFileFieldMultiTSWithoutDAS::New(fid,nomcha,i,infos);
      }
    int nProfil=MEDnProfil(fid);
    _pfls.resize(nProfil);
    for(int i=0;i<nProfil;i++)
      loadProfileInFile(fid,i);
    MEDfermer(fid);
  }
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }
