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
#include "CellModel.hxx"

#include <algorithm>

extern med_geometry_type typmai[MED_N_CELL_FIXED_GEO];
extern INTERP_KERNEL::NormalizedCellType typmai2[MED_N_CELL_FIXED_GEO];
extern med_geometry_type typmainoeud[1];
extern med_geometry_type typmai3[32];

using namespace ParaMEDMEM;

MEDFileFieldLoc *MEDFileFieldLoc::New(med_idt fid, const char *locName)
{
  return new MEDFileFieldLoc(fid,locName);
}

MEDFileFieldLoc::MEDFileFieldLoc(med_idt fid, const char *locName):_name(locName)
{
  med_geometry_type geotype;
  med_geometry_type sectiongeotype;
  int nsectionmeshcell;
  INTERP_KERNEL::AutoPtr<char> geointerpname=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> sectionmeshname=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  MEDlocalizationInfoByName(fid,locName,&geotype,&_dim,&_nb_gauss_pt,geointerpname,sectionmeshname,&nsectionmeshcell,&sectiongeotype);
  _geo_type=(INTERP_KERNEL::NormalizedCellType)(std::distance(typmai3,std::find(typmai3,typmai3+32,geotype)));
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::getCellModel(_geo_type);
  _nb_node_per_cell=cm.getNumberOfNodes();
  _ref_coo.resize(_dim*_nb_node_per_cell);
  _gs_coo.resize(_dim*_nb_gauss_pt);
  _w.resize(_nb_gauss_pt);
  MEDlocalizationRd(fid,locName,MED_FULL_INTERLACE,&_ref_coo[0],&_gs_coo[0],&_w[0]);
}

void MEDFileFieldLoc::writeLL(med_idt fid) const
{
  MEDlocalizationWr(fid,_name.c_str(),typmai3[(int)_geo_type],_dim,&_ref_coo[0],MED_FULL_INTERLACE,_nb_gauss_pt,&_gs_coo[0],&_w[0],MED_NO_INTERPOLATION,MED_NO_MESH_SUPPORT);
}

MEDFileFieldPerMeshPerTypePerDisc *MEDFileFieldPerMeshPerTypePerDisc::New(MEDFileFieldPerMeshPerType *fath, med_idt fid, int profileIt) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldPerMeshPerTypePerDisc(fath,fid,profileIt);
}

MEDFileFieldPerMeshPerTypePerDisc::MEDFileFieldPerMeshPerTypePerDisc(MEDFileFieldPerMeshPerType *fath, med_idt fid, int profileIt) throw(INTERP_KERNEL::Exception)
try:_father(fath),_profile_it(profileIt)
{
  INTERP_KERNEL::AutoPtr<char> locname=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> pflname=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  std::string fieldName=getName();
  std::string meshName=getMeshName();
  int iteration=getIteration();
  int order=getOrder();
  const std::vector<std::string>& infos=getInfo();
  TypeOfField type=getType();
  INTERP_KERNEL::NormalizedCellType geoType=getGeoType();
  int profilesize,nbi;
  med_geometry_type mgeoti;
  med_entity_type menti=MEDFileFieldPerMeshPerType::ConvertIntoMEDFileType(type,geoType,mgeoti);
  _nval=MEDfieldnValueWithProfile(fid,fieldName.c_str(),iteration,order,menti,mgeoti,profileIt,MED_COMPACT_PFLMODE,
                                  pflname,&profilesize,locname,&nbi);
  _arr=DataArrayDouble::New();
  _arr->alloc(_nval*nbi,infos.size());
  MEDfieldValueWithProfileRd(fid,fieldName.c_str(),iteration,order,menti,mgeoti,MED_COMPACT_PFLMODE,
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

double MEDFileFieldPerMeshPerTypePerDisc::getTime() const
{
  return _father->getTime();
}

std::string MEDFileFieldPerMeshPerTypePerDisc::getName() const
{
  return _father->getName();
}

std::string MEDFileFieldPerMeshPerTypePerDisc::getMeshName() const
{
  return _father->getMeshName();
}

TypeOfField MEDFileFieldPerMeshPerTypePerDisc::getType() const
{
  return _father->getType();
}

INTERP_KERNEL::NormalizedCellType MEDFileFieldPerMeshPerTypePerDisc::getGeoType() const
{
  return _father->getGeoType();
}

int MEDFileFieldPerMeshPerTypePerDisc::getNumberOfComponents() const
{
  return _father->getNumberOfComponents();
}

int MEDFileFieldPerMeshPerTypePerDisc::getNumberOfTuples() const
{
  return _arr->getNumberOfTuples();
}

const std::vector<std::string>& MEDFileFieldPerMeshPerTypePerDisc::getInfo() const
{
  return _father->getInfo();
}

std::string MEDFileFieldPerMeshPerTypePerDisc::getProfile() const
{
  return _profile;
}

std::string MEDFileFieldPerMeshPerTypePerDisc::getLocalization() const
{
  return _localization;
}

void MEDFileFieldPerMeshPerTypePerDisc::writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception)
{
  TypeOfField type=getType();
  INTERP_KERNEL::NormalizedCellType geoType=getGeoType();
  med_geometry_type mgeoti;
  med_entity_type menti=MEDFileFieldPerMeshPerType::ConvertIntoMEDFileType(type,geoType,mgeoti);
  MEDfieldValueWithProfileWr(fid,getName().c_str(),getIteration(),getOrder(),getTime(),menti,mgeoti,
                             MED_COMPACT_PFLMODE,_profile.c_str(),_localization.c_str(),MED_FULL_INTERLACE,MED_ALL_CONSTITUENT,_nval,
                             reinterpret_cast<const unsigned char*>(_arr->getConstPointer()));
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

double MEDFileFieldPerMeshPerType::getTime() const
{
  return _father->getTime();
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

const std::vector<std::string>& MEDFileFieldPerMeshPerType::getInfo() const
{
  return _father->getInfo();
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
      if(!tmp.empty() && tmp!=MED_GAUSS_ELNO)
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
  INTERP_KERNEL::NormalizedCellType geoType=getGeoType();
  INTERP_KERNEL::AutoPtr<char> pflName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> locName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  med_geometry_type mgeoti;
  med_entity_type menti=ConvertIntoMEDFileType(getType(),geoType,mgeoti);
  int nbProfiles=MEDfieldnProfile(fid,getName().c_str(),getIteration(),getOrder(),menti,mgeoti,pflName,locName);
  _field_pm_pt_pd.resize(nbProfiles);
  for(int i=0;i<nbProfiles;i++)
    {
      _field_pm_pt_pd[i]=MEDFileFieldPerMeshPerTypePerDisc::New(this,fid,i+1);
    }
}

void MEDFileFieldPerMeshPerType::writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception)
{
  _field_pm_pt_pd[0]->writeLL(fid);
}

med_entity_type MEDFileFieldPerMeshPerType::ConvertIntoMEDFileType(TypeOfField ikType, INTERP_KERNEL::NormalizedCellType ikGeoType, med_geometry_type& medfGeoType)
{
  switch(ikType)
    {
    case ON_CELLS:
      medfGeoType=typmai3[(int)ikGeoType];
      return MED_CELL;
    case ON_NODES:
      medfGeoType=MED_NONE;
      return MED_NODE;
    case ON_GAUSS_NE:
      medfGeoType=typmai3[(int)ikGeoType];
      return MED_NODE_ELEMENT;
    default:
      throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerType::ConvertIntoMEDFileType : unexpected entity type ! internal error");
    }
  return MED_UNDEF_ENTITY_TYPE;
}

MEDFileFieldPerMesh *MEDFileFieldPerMesh::New(MEDFileField1TSWithoutDAS *fath, int meshCsit, int meshIteration, int meshOrder)
{
  return new MEDFileFieldPerMesh(fath,meshCsit,meshIteration,meshOrder);
}

void MEDFileFieldPerMesh::pushBack(TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType)
{
  _types.push_back(type);
  _geo_types.push_back(geoType);
}

void MEDFileFieldPerMesh::finishLoading(med_idt fid) throw(INTERP_KERNEL::Exception)
{
  INTERP_KERNEL::AutoPtr<char> meshName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> pflName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> locName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  _field_pm_pt.clear();
  for(int i=0;i<MED_N_CELL_FIXED_GEO;i++)
    {
      int nbProfile=MEDfield23nProfile(fid,getName().c_str(),getIteration(),getOrder(),MED_CELL,typmai[i],_mesh_csit,meshName,pflName,locName);
      if(nbProfile>0)
        {
          _field_pm_pt.resize(_field_pm_pt.size()+1);
          _field_pm_pt.back()=MEDFileFieldPerMeshPerType::New(this,ON_CELLS,typmai2[i]);
          _mesh_name=MEDLoaderBase::buildStringFromFortran(meshName,MED_NAME_SIZE+1);
          _field_pm_pt.back()->finishLoading(fid);
        }
      nbProfile=MEDfield23nProfile(fid,getName().c_str(),getIteration(),getOrder(),MED_NODE_ELEMENT,typmai[i],_mesh_csit,meshName,pflName,locName);
      if(nbProfile>0)
        {
          _field_pm_pt.resize(_field_pm_pt.size()+1);
          _field_pm_pt.back()=MEDFileFieldPerMeshPerType::New(this,ON_GAUSS_NE,typmai2[i]);
          _mesh_name=MEDLoaderBase::buildStringFromFortran(meshName,MED_NAME_SIZE+1);
          _field_pm_pt.back()->finishLoading(fid);
        }
    }
  int nbProfile=MEDfield23nProfile(fid,getName().c_str(),getIteration(),getOrder(),MED_NODE,MED_NONE,_mesh_csit,meshName,pflName,locName);
  if(nbProfile>0)
    {
      _field_pm_pt.resize(_field_pm_pt.size()+1);
      _field_pm_pt.back()=MEDFileFieldPerMeshPerType::New(this,ON_NODES,INTERP_KERNEL::NORM_ERROR);
      _mesh_name=MEDLoaderBase::buildStringFromFortran(meshName,MED_NAME_SIZE+1);
      _field_pm_pt.back()->finishLoading(fid);
    }
}

void MEDFileFieldPerMesh::writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception)
{
  int nbOfTypes=_field_pm_pt.size();
  for(int i=0;i<nbOfTypes;i++)
    _field_pm_pt[i]->writeLL(fid);
}

double MEDFileFieldPerMesh::getTime() const
{
  return _father->getTime();
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

int MEDFileFieldPerMesh::getNumberOfComponents() const
{
  return _father->getNumberOfComponents();
}

const std::vector<std::string>& MEDFileFieldPerMesh::getInfo() const
{
  return _father->getInfo();
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

MEDFileFieldPerMesh::MEDFileFieldPerMesh(MEDFileField1TSWithoutDAS *fath, int meshCsit, int meshIteration, int meshOrder):_mesh_iteration(meshIteration),_mesh_order(meshOrder),
                                                                                                                          _mesh_csit(meshCsit),_father(fath)
{
}

void MEDFieldFieldGlobs::loadProfileInFile(med_idt fid, int id, const char *pflName) throw(INTERP_KERNEL::Exception)
{
  _pfls[id]=DataArrayInt::New();
  int lgth=MEDprofileSizeByName(fid,pflName);
  _pfls[id]->setName(pflName);
  _pfls[id]->alloc(lgth,1);
  MEDprofileRd(fid,pflName,_pfls[id]->getPointer());
  _pfls[id]->applyLin(1,-1,0);//Converting into C format
}

void MEDFieldFieldGlobs::loadProfileInFile(med_idt fid, int i)
{
  INTERP_KERNEL::AutoPtr<char> pflName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  int sz;
  MEDprofileInfo(fid,i+1,pflName,&sz);
  std::string pflCpp=MEDLoaderBase::buildStringFromFortran(pflName,MED_NAME_SIZE);
  _pfls[i]=DataArrayInt::New();
  _pfls[i]->alloc(sz,1);
  _pfls[i]->setName(pflCpp.c_str());
  MEDprofileRd(fid,pflName,_pfls[i]->getPointer());
}

void MEDFieldFieldGlobs::writeGlobals(med_idt fid, const MEDFileWritable& opt) const throw(INTERP_KERNEL::Exception)
{
  int nbOfPfls=_pfls.size();
  for(int i=0;i<nbOfPfls;i++)
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> cpy=_pfls[i]->deepCpy();
      cpy->applyLin(1,1,0);
      INTERP_KERNEL::AutoPtr<char> pflName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
      MEDLoaderBase::safeStrCpy(_pfls[i]->getName().c_str(),MED_NAME_SIZE,pflName,opt.getTooLongStrPolicy());
      MEDprofileWr(fid,pflName,_pfls[i]->getNumberOfTuples(),cpy->getConstPointer());
    }
  //
  int nbOfLocs=_locs.size();
  for(int i=0;i<nbOfLocs;i++)
    _locs[i]->writeLL(fid);
}

void MEDFieldFieldGlobs::loadGlobals(med_idt fid) throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> profiles=getPflsReallyUsed();
  int sz=profiles.size();
  _pfls.resize(sz);
  for(int i=0;i<sz;i++)
    loadProfileInFile(fid,i,profiles[i].c_str());
  //
  std::vector<std::string> locs=getLocsReallyUsed();
  sz=locs.size();
  _locs.resize(sz);
  for(int i=0;i<sz;i++)
    _locs[i]=MEDFileFieldLoc::New(fid,locs[i].c_str());
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

MEDFileField1TSWithoutDAS *MEDFileField1TSWithoutDAS::New(const char *fieldName, int csit, int iteration, int order, const std::vector<std::string>& infos)
{
  return new MEDFileField1TSWithoutDAS(fieldName,csit,iteration,order,infos);
}

std::string MEDFileField1TSWithoutDAS::getMeshName() const throw(INTERP_KERNEL::Exception)
{
  if(_field_per_mesh.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::getMeshName : No field set !");
  return _field_per_mesh[0]->getMeshName();
}

void MEDFileField1TSWithoutDAS::pushBack(TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType, double time) const
{
  _types.push_back(type);
  _geo_types.push_back(geoType);
  _times.push_back(time);
}

void MEDFileField1TSWithoutDAS::finishLoading(med_idt fid) throw(INTERP_KERNEL::Exception)
{
  med_int numdt,numit;
  med_float dt;
  med_int nmesh;
  med_bool localMesh;
  med_int meshnumdt,meshnumit;
  INTERP_KERNEL::AutoPtr<char> meshName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  MEDfieldComputingStepInfo(fid,_name.c_str(),_csit,&numdt,&numit,&_dt);
  MEDfield23ComputingStepMeshInfo(fid,_name.c_str(),_csit,&numdt,&numit,&dt,&nmesh,meshName,&localMesh,&meshnumdt,&meshnumit);
  if(_iteration!=numdt || _order!=numit)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutDAS::finishLoading : unexpected exception internal error !");
  _field_per_mesh.resize(nmesh);
  for(int i=0;i<nmesh;i++)
    _field_per_mesh[i]=MEDFileFieldPerMesh::New(this,i+1,meshnumdt,meshnumit);
  for(int i=0;i<nmesh;i++)
    _field_per_mesh[i]->finishLoading(fid);
}

std::vector<std::string> MEDFileField1TSWithoutDAS::getPflsReallyUsed2() const
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

std::vector<std::string> MEDFileField1TSWithoutDAS::getLocsReallyUsed2() const
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

void MEDFileField1TSWithoutDAS::writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception)
{
  if(_field_per_mesh.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutDAS::writeLL : empty field !");
  if(_field_per_mesh.size()>1)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutDAS::writeLL : In MED3.0 mode in writting mode only ONE underlying mesh supported !");
  _field_per_mesh[0]->writeLL(fid);
}

MEDCouplingFieldDouble *MEDFileField1TSWithoutDAS::getFieldAtLevel(int meshDimRelToMaxExt, const MEDFieldFieldGlobs *glob) const throw(INTERP_KERNEL::Exception)
{
  //MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> mm=MEDFileMesh::New(glob->getFileName(),getMeshName().c_str(),_mesh_iteration,_mesh_order);
  //MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> m=mm->getGenMeshAtLevel(meshDimRelToMaxExt);
  return 0;//_field_per_mesh[0]->getFieldOnMeshAtLevel(meshDimRelToMaxExt,glob,m);
}

MEDCouplingFieldDouble *MEDFileField1TSWithoutDAS::getFieldOnMeshAtLevel(int meshDimRelToMaxExt, const MEDFieldFieldGlobs *glob, const MEDCouplingMesh *mesh) const throw(INTERP_KERNEL::Exception)
{
  return _field_per_mesh[0]->getFieldOnMeshAtLevel(meshDimRelToMaxExt,glob,mesh);
}

MEDFileField1TSWithoutDAS::MEDFileField1TSWithoutDAS(const char *fieldName, int csit, int iteration, int order,
                                                     const std::vector<std::string>& infos):_name(fieldName),_infos(infos),_csit(csit),_iteration(iteration),_order(order)
{
}

MEDFileField1TS *MEDFileField1TS::New(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileField1TS(fileName,fieldName,iteration,order);
}

void MEDFileField1TS::write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception)
{
  med_access_mode medmod=MEDFileUtilities::TraduceWriteMode(mode);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,medmod);
  int nbComp=_infos.size();
  INTERP_KERNEL::AutoPtr<char> comp=MEDLoaderBase::buildEmptyString(nbComp*MED_SNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> unit=MEDLoaderBase::buildEmptyString(nbComp*MED_SNAME_SIZE);
  for(int i=0;i<nbComp;i++)
    {
      std::string info=_infos[i];
      std::string c,u;
      MEDLoaderBase::splitIntoNameAndUnit(info,c,u);
      MEDLoaderBase::safeStrCpy2(c.c_str(),MED_SNAME_SIZE-1,comp+i*MED_SNAME_SIZE,_too_long_str);
      MEDLoaderBase::safeStrCpy2(u.c_str(),MED_SNAME_SIZE-1,unit+i*MED_SNAME_SIZE,_too_long_str);
    }
  MEDfieldCr(fid,_name.c_str(),MED_FLOAT64,nbComp,comp,unit,getDtUnit().c_str(),getMeshName().c_str());
  writeGlobals(fid,*this);
  writeLL(fid);
}

MEDFileField1TS::MEDFileField1TS(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
try:MEDFileField1TSWithoutDAS(fieldName,-1,iteration,order,std::vector<std::string>()),MEDFieldFieldGlobs(fileName)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  int nbFields=MEDnField(fid);
  med_field_type typcha;
  bool found=false;
  std::vector<std::string> fns(nbFields);
  int nbOfStep2=-1;
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
          nbOfStep2=nbOfStep;
          std::string mname=MEDLoaderBase::buildStringFromFortran(nomMaa,MED_NAME_SIZE);
          _infos.resize(ncomp);
          for(int j=0;j<ncomp;j++)
            _infos[j]=MEDLoaderBase::buildUnionUnit((char *)comp+j*MED_SNAME_SIZE,MED_SNAME_SIZE,(char *)unit+j*MED_SNAME_SIZE,MED_SNAME_SIZE);
        }
    }
  if(!found)
    {
      std::ostringstream oss; oss << "No such field '" << fieldName << "' in file '" << fileName << "' ! Available fields are : ";
      std::copy(fns.begin(),fns.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  found=false;
  std::vector< std::pair<int,int> > dtits(nbOfStep2);
  for(int i=0;i<nbOfStep2 && !found;i++)
    {
      med_int numdt,numit;
      med_float dt;
      MEDfieldComputingStepInfo(fid,fieldName,i+1,&numdt,&numit,&dt);
      if(numdt==iteration && numit==order)
        {
          found=true;
          _csit=i+1;
        }
      else
        dtits[i]=std::pair<int,int>(numdt,numit);
    }
  if(!found)
    {
      std::ostringstream oss; oss << "No such iteration (" << iteration << "," << order << ") in existing field '" << fieldName << "' in file '" << fileName << "' ! Available iterations are : ";
      for(std::vector< std::pair<int,int> >::const_iterator iter=dtits.begin();iter!=dtits.end();iter++)
        oss << "(" << (*iter).first << "," << (*iter).second << "), ";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  finishLoading(fid);
  //
  loadGlobals(fid);
}
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

std::vector<std::string> MEDFileField1TS::getPflsReallyUsed() const
{
  return getPflsReallyUsed2();
}

std::vector<std::string> MEDFileField1TS::getLocsReallyUsed() const
{
  return getLocsReallyUsed2();
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
try:_name(fieldName)
{
  finishLoading(fid,nbOfStep);
}
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

const std::vector<std::string>& MEDFileFieldMultiTSWithoutDAS::getInfo() const throw(INTERP_KERNEL::Exception)
{
  if(_time_steps.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutDAS::getInfos : not time steps !");
  return _time_steps[0]->getInfo();
}

std::string MEDFileFieldMultiTSWithoutDAS::getMeshName() const throw(INTERP_KERNEL::Exception)
{
  if(_time_steps.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutDAS::getMeshName : not time steps !");
  return _time_steps[0]->getMeshName();
}

std::string MEDFileFieldMultiTSWithoutDAS::getDtUnit() const throw(INTERP_KERNEL::Exception)
{
  if(_time_steps.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutDAS::getMeshName : not time steps !");
  return _time_steps[0]->getDtUnit();
}

void MEDFileFieldMultiTSWithoutDAS::finishLoading(med_idt fid, int nbPdt) throw(INTERP_KERNEL::Exception)
{
  _time_steps.resize(nbPdt);
  for(int i=0;i<nbPdt;i++)
    {
      std::vector< std::pair<int,int> > ts;
      med_int numdt=0,numo=0;
      med_int meshIt=0,meshOrder=0;
      med_float dt=0.0;
      MEDfieldComputingStepMeshInfo(fid,_name.c_str(),i+1,&numdt,&numo,&dt,&meshIt,&meshOrder);
      _time_steps[i]=MEDFileField1TSWithoutDAS::New(_name.c_str(),i+1,numdt,numo,_infos);
      _time_steps[i]->finishLoading(fid);
    }
}

void MEDFileFieldMultiTSWithoutDAS::writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception)
{
  if(_time_steps.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutDAS::writeLL : no time steps set !");
  std::vector<std::string> infos(getInfo());
  int nbComp=infos.size();
  INTERP_KERNEL::AutoPtr<char> comp=MEDLoaderBase::buildEmptyString(nbComp*MED_SNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> unit=MEDLoaderBase::buildEmptyString(nbComp*MED_SNAME_SIZE);
  for(int i=0;i<nbComp;i++)
    {
      std::string info=infos[i];
      std::string c,u;
      MEDLoaderBase::splitIntoNameAndUnit(info,c,u);
      MEDLoaderBase::safeStrCpy2(c.c_str(),MED_SNAME_SIZE-1,comp+i*MED_SNAME_SIZE,_too_long_str);
      MEDLoaderBase::safeStrCpy2(u.c_str(),MED_SNAME_SIZE-1,unit+i*MED_SNAME_SIZE,_too_long_str);
    }
  MEDfieldCr(fid,_name.c_str(),MED_FLOAT64,nbComp,comp,unit,getDtUnit().c_str(),getMeshName().c_str());
  int nbOfTS=_time_steps.size();
  for(int i=0;i<nbOfTS;i++)
    {
      _time_steps[i]->copyOptionsFrom(*this);
      _time_steps[i]->writeLL(fid);
    }
}

int MEDFileFieldMultiTSWithoutDAS::getNumberOfTS() const
{
  return _time_steps.size();
}

MEDFileFieldMultiTS *MEDFileFieldMultiTS::New(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldMultiTS(fileName,fieldName);
}

void MEDFileFieldMultiTS::write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception)
{
  med_access_mode medmod=MEDFileUtilities::TraduceWriteMode(mode);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,medmod);
  writeGlobals(fid,*this);
  writeLL(fid);
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
            _infos[j]=MEDLoaderBase::buildUnionUnit((char *)comp+j*MED_SNAME_SIZE,MED_SNAME_SIZE,(char *)unit+j*MED_SNAME_SIZE,MED_SNAME_SIZE);
        }
    }
  if(!found)
    {
      std::ostringstream oss; oss << "No such field '" << fieldName << "' in file '" << fileName << "' ! Available fields are : ";
      std::copy(fns.begin(),fns.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  finishLoading(fid,nbstep2);
  loadGlobals(fid);
}
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

std::vector<std::string> MEDFileFieldMultiTS::getPflsReallyUsed() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileField1TSWithoutDAS > >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
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

std::vector<std::string> MEDFileFieldMultiTS::getLocsReallyUsed() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileField1TSWithoutDAS > >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
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

std::vector<std::string> MEDFileFields::getPflsReallyUsed() const
{
  return std::vector<std::string>();//tony
}

std::vector<std::string> MEDFileFields::getLocsReallyUsed() const
{
  return std::vector<std::string>();//tony
}

