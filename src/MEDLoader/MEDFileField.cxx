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

#include "MEDFileField.hxx"
#include "MEDFileMesh.hxx"
#include "MEDLoaderBase.hxx"
#include "MEDFileUtilities.hxx"

#include "InterpKernelAutoPtr.hxx"

#include <algorithm>

extern med_geometry_type typmai[MED_N_CELL_FIXED_GEO];
extern INTERP_KERNEL::NormalizedCellType typmai2[MED_N_CELL_FIXED_GEO];
extern med_geometry_type typmainoeud[1];
extern med_geometry_type typmai3[32];

using namespace ParaMEDMEM;

MEDFileFieldPerMeshPerTypePerDisc *MEDFileFieldPerMeshPerTypePerDisc::New(MEDFileFieldPerMeshPerType *fath, med_idt fid) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldPerMeshPerTypePerDisc(fath,fid);
}

MEDFileFieldPerMeshPerTypePerDisc::MEDFileFieldPerMeshPerTypePerDisc(MEDFileFieldPerMeshPerType *fath, med_idt fid) throw(INTERP_KERNEL::Exception)
try:_father(fath)
{
  INTERP_KERNEL::AutoPtr<char> locname=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> pflname=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  std::string fieldName=getName();
  std::string meshName=getMeshName();
  int iteration=getIteration();
  int order=getOrder();
  const std::vector<std::string>& infos=getInfos();
  TypeOfField type=getType();
  INTERP_KERNEL::NormalizedCellType geoType=getGeoType();
  int profilesize,nbi;
  int nval=MEDfieldnValueWithProfile(fid,fieldName.c_str(),iteration,order,type==ON_CELLS?MED_CELL:MED_NODE,type==ON_CELLS?typmai3[(int)geoType]:MED_NONE,1,MED_COMPACT_PFLMODE,
                                     pflname,&profilesize,locname,&nbi);//to generalize
  _arr=DataArrayDouble::New();
  _arr->alloc(nval,infos.size());
  MEDfieldValueWithProfileRd(fid,fieldName.c_str(),iteration,order,type==ON_CELLS?MED_CELL:MED_NODE,type==ON_CELLS?typmai3[(int)geoType]:MED_NONE,MED_COMPACT_PFLMODE,
                             pflname,MED_FULL_INTERLACE,MED_ALL_CONSTITUENT,reinterpret_cast<unsigned char*>(_arr->getPointer()));
  _profile=MEDLoaderBase::buildStringFromFortran(pflname,MED_NAME_SIZE);
  _localization=MEDLoaderBase::buildStringFromFortran(locname,MED_NAME_SIZE);
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

int MEDFileFieldPerMeshPerTypePerDisc::getNumberOfTuples() const
{
  return _arr->getNumberOfTuples();
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

bool MEDFileFieldPerMeshPerType::isOnNode(int& type, int& number, const DataArrayInt* &arrs) const throw(INTERP_KERNEL::Exception)
{
  if(_type!=ON_NODES)
    return false;
  type=INTERP_KERNEL::NORM_ERROR;
  number=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++)
    number+=(*it)->getNumberOfTuples();
  return true;
}

bool MEDFileFieldPerMeshPerType::isOnCell(int dimDimReq, int& type, int& number, const DataArrayInt* &arrs) const throw(INTERP_KERNEL::Exception)
{
  if(_type!=ON_CELLS)
    return false;
  type=(int)_geo_type;
  number=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++)
    number+=(*it)->getNumberOfTuples();
  return true;
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

MEDFileFieldPerMesh *MEDFileFieldPerMesh::New(MEDFileField1TSWithoutDAS *fath, double time)
{
  return new MEDFileFieldPerMesh(fath,time);
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
  return _father->getMeshName();
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

MEDCouplingFieldDouble *MEDFileFieldPerMesh::getFieldOnMeshAtLevel(int meshDimRelToMaxExt, const MEDFieldFieldGlobs *glob, const MEDCouplingMesh *mesh) const throw(INTERP_KERNEL::Exception)
{
  std::vector<const DataArrayInt *> arrs;
  std::vector<int> distrib=getDistributionOfTypes(meshDimRelToMaxExt,mesh->getMeshDimension(),arrs);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da=mesh->checkTypeConsistencyAndContig(distrib,arrs);
  return 0;
}

std::vector<int> MEDFileFieldPerMesh::getDistributionOfTypes(int meshDimRelToMaxExt, int mdim, std::vector<const DataArrayInt *>& arrs) const throw(INTERP_KERNEL::Exception)
{
  if(meshDimRelToMaxExt>1)
    throw INTERP_KERNEL::Exception("Invalid meshDimRelToMaxExt parameter passed ! must be 1 for node and 0,-1,-2,-3 for cells !");
  std::vector<int> ret;
  int tmp2,tmp3;
  arrs.clear();
  const DataArrayInt *arrTmp=0;
  if(meshDimRelToMaxExt==1)
    {//On node
      for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerType> >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
        if((*it)->isOnNode(tmp2,tmp3,arrTmp))
          {
            ret.push_back(-1);
            ret.push_back(tmp2);
            ret.push_back(tmp3);
            if(arrTmp)
              arrs.push_back(arrTmp);
          }
      if(ret.size()>2)
        throw INTERP_KERNEL::Exception("Detected on same time step two different discretization ... Should never happen !");
    }
  else
    {
      int dimTarget=mdim+meshDimRelToMaxExt;
      for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerType> >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
        if((*it)->isOnCell(dimTarget,tmp2,tmp3,arrTmp))
          {
            ret.push_back((int)(*it)->getGeoType());
            ret.push_back(tmp2);
            ret.push_back(tmp3);
            if(arrTmp)
              arrs.push_back(arrTmp);
          }
    }
  if(ret.empty())
    throw INTERP_KERNEL::Exception("No field part correspond to requested meshdimRel !");
  return ret;
}

MEDFileFieldPerMesh::MEDFileFieldPerMesh(MEDFileField1TSWithoutDAS *fath, double time):_father(fath),_time(time)
{
}

void MEDFieldFieldGlobs::loadProfileInFile(med_idt fid, int id, const char *pflName, int lgth) throw(INTERP_KERNEL::Exception)
{
  _pfls[id]=DataArrayInt::New();
  _pfls[id]->alloc(lgth,1);
  MEDprofileRd(fid,pflName,_pfls[id]->getPointer());
  _pfls[id]->applyLin(1,-1,0);//Converting into C format
}

void MEDFieldFieldGlobs::loadProfileInFile(med_idt fid, int i)
{
  INTERP_KERNEL::AutoPtr<char> pflName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  int sz;
  MEDprofileInfo(fid,i,pflName,&sz);
  _pfls[i]=DataArrayInt::New();
  _pfls[i]->alloc(sz,1);
  MEDprofileRd(fid,pflName,_pfls[i]->getPointer());
}

MEDFieldFieldGlobs::MEDFieldFieldGlobs(const char *fname):_file_name(fname)
{
}

void MEDFieldFieldGlobs::setFileName(const char *fileName)
{
  _file_name=fileName;
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

MEDFileField1TSWithoutDAS *MEDFileField1TSWithoutDAS::New(const char *fieldName, const char *meshName, int iteration, int order, int meshIt, int meshOrder, const std::vector<std::string>& infos)
{
  return new MEDFileField1TSWithoutDAS(fieldName,meshName,iteration,order,meshIt,meshOrder,infos);
}

void MEDFileField1TSWithoutDAS::pushBack(TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType, double time) const
{
  _types.push_back(type);
  _geo_types.push_back(geoType);
  _times.push_back(time);
}

void MEDFileField1TSWithoutDAS::finishLoading(med_idt fid) throw(INTERP_KERNEL::Exception)
{
  _field_per_mesh.push_back(MEDFileFieldPerMesh::New(this,_times[0]));
  _field_per_mesh[0]->finishLoading(fid);
}

MEDCouplingFieldDouble *MEDFileField1TSWithoutDAS::getFieldAtLevel(int meshDimRelToMaxExt, const MEDFieldFieldGlobs *glob) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> mm=MEDFileMesh::New(glob->getFileName(),_mesh_name.c_str(),_mesh_iteration,_mesh_order);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> m=mm->getGenMeshAtLevel(meshDimRelToMaxExt);
  return _field_per_mesh[0]->getFieldOnMeshAtLevel(meshDimRelToMaxExt,glob,m);
}

MEDCouplingFieldDouble *MEDFileField1TSWithoutDAS::getFieldOnMeshAtLevel(int meshDimRelToMaxExt, const MEDFieldFieldGlobs *glob, const MEDCouplingMesh *mesh) const throw(INTERP_KERNEL::Exception)
{
  return _field_per_mesh[0]->getFieldOnMeshAtLevel(meshDimRelToMaxExt,glob,mesh);
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

MEDFileField1TSWithoutDAS::MEDFileField1TSWithoutDAS(const char *fieldName, const char *meshName, int iteration, int order, int meshIt, int meshOrder,
                                                     const std::vector<std::string>& infos):_name(fieldName),_mesh_name(meshName),_infos(infos),_iteration(iteration),_order(order),
                                                                                            _mesh_iteration(meshIt),_mesh_order(meshOrder)
{
}

MEDFileField1TS *MEDFileField1TS::New(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileField1TS(fileName,fieldName,iteration,order);
}

MEDFileField1TS::MEDFileField1TS(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
try:MEDFileField1TSWithoutDAS(fieldName,"",iteration,order,-1,-1,std::vector<std::string>()),MEDFieldFieldGlobs(fileName)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  int nbFields=MEDnField(fid);
  med_field_type typcha;
  bool found=false;
  std::vector<std::string> fns(nbFields);
  for(int i=0;i<nbFields && !found;i++)
    {
      int ncomp=MEDfieldnComponent(fid,i+1);
      INTERP_KERNEL::AutoPtr<char> comp=MEDLoaderBase::buildEmptyString(ncomp*MED_SNAME_SIZE);
      INTERP_KERNEL::AutoPtr<char> unit=MEDLoaderBase::buildEmptyString(ncomp*MED_SNAME_SIZE);
      INTERP_KERNEL::AutoPtr<char> dtunit=MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE);
      INTERP_KERNEL::AutoPtr<char> nomcha=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
      INTERP_KERNEL::AutoPtr<char> nomMaa=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
      med_bool localMesh;
      int nbOfStep;
      MEDfieldInfo(fid,i+1,nomcha,nomMaa,&localMesh,&typcha,comp,unit,dtunit,&nbOfStep);
      std::string tmp(nomcha);
      fns[i]=tmp;
      found=(tmp==fieldName);
      if(found)
        {
          std::string mname=MEDLoaderBase::buildStringFromFortran(nomMaa,MED_NAME_SIZE);
          _mesh_name=mname;
          _infos.resize(ncomp);
          for(int j=0;j<ncomp;j++)
            _infos[j]=MEDLoaderBase::buildUnionUnit((char *)comp+i*MED_SNAME_SIZE,MED_SNAME_SIZE,(char *)unit+i*MED_SNAME_SIZE,MED_SNAME_SIZE);
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
}
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

MEDCouplingFieldDouble *MEDFileField1TS::getFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception)
{
  return MEDFileField1TSWithoutDAS::getFieldAtLevel(meshDimRelToMaxExt,this);
}

MEDCouplingFieldDouble *MEDFileField1TS::getFieldOnMeshAtLevel(int meshDimRelToMaxExt, const MEDCouplingMesh *mesh) const throw(INTERP_KERNEL::Exception)
{
  return MEDFileField1TSWithoutDAS::getFieldOnMeshAtLevel(meshDimRelToMaxExt,this,mesh);
}

MEDFileFieldMultiTSWithoutDAS *MEDFileFieldMultiTSWithoutDAS::New(med_idt fid, const char *fieldName, int id, const std::vector<std::string>& infos, int nbOfStep) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldMultiTSWithoutDAS(fid,fieldName,id,infos,nbOfStep);
}

MEDFileFieldMultiTSWithoutDAS::MEDFileFieldMultiTSWithoutDAS(const char *fieldName):_name(fieldName)
{
}

MEDFileFieldMultiTSWithoutDAS::MEDFileFieldMultiTSWithoutDAS(med_idt fid, const char *fieldName, int id, const std::vector<std::string>& infos, int nbOfStep) throw(INTERP_KERNEL::Exception)
try:_name(fieldName),_infos(infos)
{
  finishLoading(fid,nbOfStep);
}
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

void MEDFileFieldMultiTSWithoutDAS::finishLoading(med_idt fid, int nbPdt) throw(INTERP_KERNEL::Exception)
{
  for(int i=0;i<MED_N_CELL_FIXED_GEO;i++)
    for(int j=0;j<nbPdt;j++)
      appendTimeStepEntry(fid,MED_CELL,i,j);
  for(int j=0;j<nbPdt;j++)
    appendTimeStepEntry(fid,MED_NODE,0,j);
  int nbOfTimeSteps=_time_steps.size();
  for(int i=0;i<nbOfTimeSteps;i++)
    _time_steps[i]->finishLoading(fid);
}

void MEDFileFieldMultiTSWithoutDAS::appendTimeStepEntry(med_idt fid, med_entity_type entity, int i, int j) throw(INTERP_KERNEL::Exception)
{
  std::vector< std::pair<int,int> > ts;
  med_int numdt=0,numo=0;
  med_int meshIt=0,meshOrder=0;
  med_float dt=0.0;
  MEDfieldComputingStepMeshInfo(fid,_name.c_str(),j+1,&numdt,&numo,&dt,&meshIt,&meshOrder);
  std::pair<int,int> p(numdt,numo);
  std::vector< std::pair<int,int> >::iterator where=std::find(ts.begin(),ts.end(),p);
  if(where==ts.end())
    {
      ts.push_back(p);
      _time_steps.push_back(MEDFileField1TSWithoutDAS::New(_name.c_str(),_mesh_name.c_str(),numdt,numo,meshIt,meshOrder,_infos));
      _time_steps.back()->pushBack(entity==MED_CELL?ON_CELLS:ON_NODES,entity==MED_CELL?typmai2[i]:INTERP_KERNEL::NORM_ERROR,dt);
    }
  else
    {
      int w=std::distance(ts.begin(),where);
      _time_steps[w]->pushBack(entity==MED_CELL?ON_CELLS:ON_NODES,entity==MED_CELL?typmai2[i]:INTERP_KERNEL::NORM_ERROR,dt);
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

int MEDFileFieldMultiTSWithoutDAS::getNumberOfTS() const
{
  return _time_steps.size();
}

MEDFileFieldMultiTS *MEDFileFieldMultiTS::New(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldMultiTS(fileName,fieldName);
}

MEDFileFieldMultiTS::MEDFileFieldMultiTS(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
try:MEDFileFieldMultiTSWithoutDAS(fieldName),MEDFieldFieldGlobs(fileName)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  int nbFields=MEDnField(fid);
  med_field_type typcha;
  bool found=false;
  std::vector<std::string> fns(nbFields);
  int nbstep2=-1;
  for(int i=0;i<nbFields && !found;i++)
    {
      int ncomp=MEDfieldnComponent(fid,i+1);
      INTERP_KERNEL::AutoPtr<char> comp=MEDLoaderBase::buildEmptyString(ncomp*MED_SNAME_SIZE);
      INTERP_KERNEL::AutoPtr<char> unit=MEDLoaderBase::buildEmptyString(ncomp*MED_SNAME_SIZE);
      INTERP_KERNEL::AutoPtr<char> dtunit=MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE);
      INTERP_KERNEL::AutoPtr<char> nomcha=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
      INTERP_KERNEL::AutoPtr<char> nomMaa=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
      med_bool localMesh;
      int nbOfStep;
      MEDfieldInfo(fid,i+1,nomcha,nomMaa,&localMesh,&typcha,comp,unit,dtunit,&nbOfStep);
      std::string tmp(nomcha);
      fns[i]=tmp;
      found=(tmp==fieldName);
      if(found)
        {
          nbstep2=nbOfStep;
          _infos.resize(ncomp);
          for(int j=0;j<ncomp;j++)
            _infos[j]=MEDLoaderBase::buildUnionUnit((char *)comp+i*MED_SNAME_SIZE,MED_SNAME_SIZE,(char *)unit+i*MED_SNAME_SIZE,MED_SNAME_SIZE);
        }
    }
  if(!found)
    {
      std::ostringstream oss; oss << "No such field '" << fieldName << "' in file '" << fileName << "' ! Available fields are : ";
      std::copy(fns.begin(),fns.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  finishLoading(fid,nbstep2);
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
try:MEDFieldFieldGlobs(fileName)
  {
    MEDFileUtilities::CheckFileForRead(fileName);
    MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
    int nbFields=MEDnField(fid);
    _fields.resize(nbFields);
    med_field_type typcha;
    for(int i=0;i<nbFields;i++)
      {
        int ncomp=MEDfieldnComponent(fid,i+1);
        INTERP_KERNEL::AutoPtr<char> comp=MEDLoaderBase::buildEmptyString(ncomp*MED_SNAME_SIZE);
        INTERP_KERNEL::AutoPtr<char> unit=MEDLoaderBase::buildEmptyString(ncomp*MED_SNAME_SIZE);
        INTERP_KERNEL::AutoPtr<char> dtunit=MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE);
        INTERP_KERNEL::AutoPtr<char> nomcha=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
        INTERP_KERNEL::AutoPtr<char> nomMaa=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
        med_bool localMesh;
        int nbOfStep;
        MEDfieldInfo(fid,i+1,nomcha,nomMaa,&localMesh,&typcha,comp,unit,dtunit,&nbOfStep);
        std::vector<std::string> infos(ncomp);
        for(int j=0;j<ncomp;j++)
          infos[j]=MEDLoaderBase::buildUnionUnit((char *)comp+i*MED_SNAME_SIZE,MED_SNAME_SIZE,(char *)unit+i*MED_SNAME_SIZE,MED_SNAME_SIZE);
        _fields[i]=MEDFileFieldMultiTSWithoutDAS::New(fid,nomcha,i,infos,nbOfStep);
      }
    int nProfil=MEDnProfile(fid);
    _pfls.resize(nProfil);
    for(int i=0;i<nProfil;i++)
      loadProfileInFile(fid,i);
  }
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }
