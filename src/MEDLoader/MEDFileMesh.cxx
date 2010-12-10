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

#include "MEDFileMesh.hxx"
#include "MEDFileUtilities.hxx"
#include "MEDLoader.hxx"
#include "MEDLoaderBase.hxx"

#include "MEDCouplingUMesh.hxx"

#include "InterpKernelAutoPtr.hxx"

using namespace ParaMEDMEM;

MEDFileMesh *MEDFileMesh::New(const char *fileName, const char *mName) throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("Not implemented yet !");
}

MEDFileUMesh *MEDFileUMesh::New(const char *fileName, const char *mName) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileUMesh(fileName,mName);
}

MEDFileUMesh *MEDFileUMesh::New(const char *fileName) throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> ms=MEDLoader::GetMeshNames(fileName);
  if(ms.empty())
    {
      std::ostringstream oss; oss << "MEDFileUMesh::New : no meshes in file \"" << fileName << "\" !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return new MEDFileUMesh(fileName,ms.front().c_str());
}

MEDFileUMesh *MEDFileUMesh::New()
{
  return new MEDFileUMesh;
}

MEDFileUMesh::MEDFileUMesh():_too_long_str(0)
{
}

MEDFileUMesh::MEDFileUMesh(const char *fileName, const char *mName) throw(INTERP_KERNEL::Exception)
  try:_too_long_str(0)
  {
    MEDFileUtilities::CheckFileForRead(fileName);
    med_idt fid=MEDouvrir((char *)fileName,MED_LECTURE);
    MEDFileUMeshL2 loaderl2;
    loaderl2.loadAll(fid,MEDFileUMeshL2::getMeshIdFromName(fid,mName),mName);
    int lev=loaderl2.getNumberOfLevels();
     _ms.resize(lev);
     for(int i=0;i<lev;i++)
       {
         if(!loaderl2.emptyLev(i))
           _ms[i]=new MEDFileUMeshSplitL1(loaderl2,mName,i);
         else
           _ms[i]=0;
       }
    MEDFileUMeshL2::readFamiliesAndGrps(fid,mName,_families,_groups);
    MEDfermer(fid);
    //
    setName(loaderl2.getName());
    setDescription(loaderl2.getDescription());
    _coords=loaderl2.getCoords();
    _fam_coords=loaderl2.getCoordsFamily();
    _num_coords=loaderl2.getCoordsNum();
  }
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

MEDFileUMesh::~MEDFileUMesh()
{
}

void MEDFileUMesh::write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception)
{
  if(_name.empty())
    throw INTERP_KERNEL::Exception("MEDFileUMesh : name is empty. MED file ask for a NON EMPTY name !");
  med_mode_acces medmod=MEDFileUtilities::TraduceWriteMode(mode);
  med_idt fid=MEDouvrir((char *)fileName,medmod);
  std::ostringstream oss; oss << "MEDFileUMesh : error on attempt to write in file : \"" << fileName << "\""; 
  MEDFileUtilities::CheckMEDCode(fid,fid,oss.str().c_str());
  const DataArrayDouble *coo=_coords;
  INTERP_KERNEL::AutoPtr<char> maa=MEDLoaderBase::buildEmptyString(MED_TAILLE_NOM);
  INTERP_KERNEL::AutoPtr<char> desc=MEDLoaderBase::buildEmptyString(MED_TAILLE_DESC);
  MEDLoaderBase::safeStrCpy(_name.c_str(),MED_TAILLE_NOM,maa,_too_long_str);
  MEDLoaderBase::safeStrCpy(_desc_name.c_str(),MED_TAILLE_DESC,desc,_too_long_str);
  int spaceDim=coo?coo->getNumberOfComponents():0;
  MEDmaaCr(fid,maa,spaceDim,MED_NON_STRUCTURE,desc);
  MEDdimEspaceCr(fid,maa,spaceDim);
  MEDFileUMeshL2::writeCoords(fid,maa,_coords,_fam_coords,_num_coords);
  int mdim=getMeshDimension();
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::const_iterator it=_ms.begin();it!=_ms.end();it++)
    if((const MEDFileUMeshSplitL1 *)(*it)!=0)
      (*it)->write(fid,maa,mdim);
  MEDFileUMeshL2::writeFamiliesAndGrps(fid,maa,_families,_groups);
  MEDfermer(fid);
}

int MEDFileUMesh::getNumberOfNonEmptyLevels() const
{
  int ret=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::const_iterator it=_ms.begin();it!=_ms.end();it++)
    if((const MEDFileUMeshSplitL1 *)(*it)!=0)
      if(!(*it)->empty())
        ret++;
  return ret;
}

std::vector<int> MEDFileUMesh::getNonEmptyLevels() const
{
  std::vector<int> ret;
  int lev=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::const_iterator it=_ms.begin();it!=_ms.end();it++,lev--)
    if((const MEDFileUMeshSplitL1 *)(*it)!=0)
      if(!(*it)->empty())
        ret.push_back(lev);
  return ret;
}

int MEDFileUMesh::getFamilyId(const char *name) const throw(INTERP_KERNEL::Exception)
{
  std::string oname(name);
  std::map<std::string, int>::const_iterator it=_families.find(oname);
  std::vector<std::string> fams=getFamiliesNames();
  if(it==_families.end())
    {
      std::ostringstream oss; oss << "No such familyname \"" << name << "\" !\nAvailable families are :";
      std::copy(fams.begin(),fams.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return (*it).second;
}

std::vector<int> MEDFileUMesh::getFamiliesIds(const std::vector<std::string>& famNames) const throw(INTERP_KERNEL::Exception)
{
  std::vector<int> famIds;
  for(std::vector<std::string>::const_iterator it=famNames.begin();it!=famNames.end();it++)
    {
      std::map<std::string,int>::const_iterator it2=_families.find(*it);
      if(it2==_families.end())
        {
          std::ostringstream oss; oss << "No such family in mesh \"" << _name << "\" : " << *it; 
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      famIds.push_back((*it2).second);
    }
  return famIds;
}

std::string MEDFileUMesh::getFamilyNameGivenId(int id) const throw(INTERP_KERNEL::Exception)
{
  for(std::map<std::string,int>::const_iterator it=_families.begin();it!=_families.end();it++)
    if((*it).second==id)
      return (*it).first;
  std::ostringstream oss; oss << "MEDFileUMesh::getFamilyNameGivenId : no such family id : " << id;
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

int MEDFileUMesh::getMeshDimension() const
{
  int lev=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::const_iterator it=_ms.begin();it!=_ms.end();it++,lev++)
    if((const MEDFileUMeshSplitL1 *)(*it)!=0)
      return (*it)->getMeshDimension()+lev;
  throw INTERP_KERNEL::Exception("MEDFileUMesh::getMeshDimension : impossible to find a mesh dimension !");
}

std::vector<std::string> MEDFileUMesh::getFamiliesOnGroup(const char *name) const throw(INTERP_KERNEL::Exception)
{
  std::string oname(name);
  std::map<std::string, std::vector<std::string> >::const_iterator it=_groups.find(oname);
  std::vector<std::string> grps=getGroupsNames();
  if(it==_groups.end())
    {
      std::ostringstream oss; oss << "No such groupname \"" << name << "\" !\nAvailable groups are :";
      std::copy(grps.begin(),grps.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return (*it).second;
}

std::vector<std::string> MEDFileUMesh::getGroupsNames() const
{
  std::vector<std::string> ret(_groups.size());
  int i=0;
  for(std::map<std::string, std::vector<std::string> >::const_iterator it=_groups.begin();it!=_groups.end();it++,i++)
    ret[i]=(*it).first;
  return ret;
}

std::vector<std::string> MEDFileUMesh::getFamiliesNames() const
{
  std::vector<std::string> ret(_families.size());
  int i=0;
  for(std::map<std::string, int >::const_iterator it=_families.begin();it!=_families.end();it++,i++)
    ret[i]=(*it).first;
  return ret;
}

void MEDFileUMesh::removeGroup(const char *name) throw(INTERP_KERNEL::Exception)
{
  std::string oname(name);
  std::map<std::string, std::vector<std::string> >::iterator it=_groups.find(oname);
  std::vector<std::string> grps=getGroupsNames();
  if(it==_groups.end())
    {
      std::ostringstream oss; oss << "No such groupname \"" << name << "\" !\nAvailable groups are :";
      std::copy(grps.begin(),grps.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  _groups.erase(it);
}

void MEDFileUMesh::removeFamily(const char *name) throw(INTERP_KERNEL::Exception)
{
  std::string oname(name);
  std::map<std::string, int >::iterator it=_families.find(oname);
  std::vector<std::string> fams=getFamiliesNames();
  if(it==_families.end())
    {
      std::ostringstream oss; oss << "No such familyname \"" << name << "\" !\nAvailable families are :";
      std::copy(fams.begin(),fams.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  _families.erase(it);
}

void MEDFileUMesh::changeGroupName(const char *oldName, const char *newName) throw(INTERP_KERNEL::Exception)
{
  std::string oname(oldName);
  std::map<std::string, std::vector<std::string> >::iterator it=_groups.find(oname);
  std::vector<std::string> grps=getGroupsNames();
  if(it==_groups.end())
    {
      std::ostringstream oss; oss << "No such groupname \"" << oldName << "\" !\nAvailable groups are :";
      std::copy(grps.begin(),grps.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  std::string nname(newName);
  it=_groups.find(nname);
  if(it!=_groups.end())
    {
      std::ostringstream oss; oss << "Such groupname \"" << newName << " already exists ! Kill it before !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  std::vector<std::string> cpy=(*it).second;
  _groups.erase(it);
  _groups[newName]=cpy;
}

void MEDFileUMesh::changeFamilyName(const char *oldName, const char *newName) throw(INTERP_KERNEL::Exception)
{
  std::string oname(oldName);
  std::map<std::string, int >::iterator it=_families.find(oname);
  std::vector<std::string> fams=getFamiliesNames();
  if(it==_families.end())
    {
      std::ostringstream oss; oss << "No such familyname \"" << oldName << "\" !\nAvailable families are :";
      std::copy(fams.begin(),fams.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  std::string nname(newName);
  it=_families.find(nname);
  if(it!=_families.end())
    {
      std::ostringstream oss; oss << "Such familyname \"" << newName << " already exists ! Kill it before !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  int cpy=(*it).second;
  _families.erase(it);
  _families[newName]=cpy;
}

DataArrayDouble *MEDFileUMesh::getCoords() const
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> tmp(_coords);
  if((DataArrayDouble *)tmp)
    {
      tmp->incrRef();
      return tmp;
    }
  return 0;
}

MEDCouplingUMesh *MEDFileUMesh::getGroup(int meshDimRelToMax, const char *grp, bool renum) const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> tmp(1);
  tmp[0]=grp;
  MEDCouplingUMesh *ret=getGroups(meshDimRelToMax,tmp,renum);
  ret->setName(grp);
  return ret;
}

DataArrayInt *MEDFileUMesh::getGroupArr(int meshDimRelToMax, const char *grp, bool renum) const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> tmp(1);
  tmp[0]=grp;
  DataArrayInt *ret=getGroupsArr(meshDimRelToMax,tmp,renum);
  ret->setName(grp);
  return ret;
}

MEDCouplingUMesh *MEDFileUMesh::getGroups(int meshDimRelToMax, const std::vector<std::string>& grps, bool renum) const throw(INTERP_KERNEL::Exception)
{
  std::set<std::string> fams;
  for(std::vector<std::string>::const_iterator it=grps.begin();it!=grps.end();it++)
    {
      std::map<std::string, std::vector<std::string> >::const_iterator it2=_groups.find(*it);
      if(it2==_groups.end())
        {
          std::ostringstream oss; oss << "No such group in mesh \"" << _name << "\" : " << *it; 
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      fams.insert((*it2).second.begin(),(*it2).second.end());
    }
  std::vector<std::string> fams2(fams.begin(),fams.end());
  return getFamilies(meshDimRelToMax,fams2,renum);
}

DataArrayInt *MEDFileUMesh::getGroupsArr(int meshDimRelToMax, const std::vector<std::string>& grps, bool renum) const throw(INTERP_KERNEL::Exception)
{
  std::set<std::string> fams;
  for(std::vector<std::string>::const_iterator it=grps.begin();it!=grps.end();it++)
    {
      std::map<std::string, std::vector<std::string> >::const_iterator it2=_groups.find(*it);
      if(it2==_groups.end())
        {
          std::ostringstream oss; oss << "No such group in mesh \"" << _name << "\" : " << *it; 
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      fams.insert((*it2).second.begin(),(*it2).second.end());
    }
  std::vector<std::string> fams2(fams.begin(),fams.end());
  return getFamiliesArr(meshDimRelToMax,fams2,renum);
}

MEDCouplingUMesh *MEDFileUMesh::getFamily(int meshDimRelToMax, const char *fam, bool renum) const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> tmp(1);
  tmp[0]=fam;
  MEDCouplingUMesh *ret=getFamilies(meshDimRelToMax,tmp,renum);
  ret->setName(fam);
  return ret;
}

DataArrayInt *MEDFileUMesh::getFamilyArr(int meshDimRelToMax, const char *fam, bool renum) const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> tmp(1);
  tmp[0]=fam;
  DataArrayInt *ret=getFamiliesArr(meshDimRelToMax,tmp,renum);
  ret->setName(fam);
  return ret;
}

MEDCouplingUMesh *MEDFileUMesh::getFamilies(int meshDimRelToMax, const std::vector<std::string>& fams, bool renum) const throw(INTERP_KERNEL::Exception)
{
  std::vector<int> famIds=getFamiliesIds(fams);
  const MEDFileUMeshSplitL1 *l1=getMeshAtLevSafe(meshDimRelToMax);
  return l1->getFamilyPart(famIds,renum);
}

DataArrayInt *MEDFileUMesh::getFamiliesArr(int meshDimRelToMax, const std::vector<std::string>& fams, bool renum) const throw(INTERP_KERNEL::Exception)
{
  std::vector<int> famIds=getFamiliesIds(fams);
  const MEDFileUMeshSplitL1 *l1=getMeshAtLevSafe(meshDimRelToMax);
  return l1->getFamilyPartArr(famIds,renum);
}

DataArrayInt *MEDFileUMesh::getNodeGroupArr(const char *grp, bool renum) const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> tmp(1);
  tmp[0]=grp;
  DataArrayInt *ret=getNodeGroupsArr(tmp,renum);
  ret->setName(grp);
  return ret;
}

DataArrayInt *MEDFileUMesh::getNodeGroupsArr(const std::vector<std::string>& grps, bool renum) const throw(INTERP_KERNEL::Exception)
{
  std::set<std::string> fams;
  for(std::vector<std::string>::const_iterator it=grps.begin();it!=grps.end();it++)
    {
      std::map<std::string, std::vector<std::string> >::const_iterator it2=_groups.find(*it);
      if(it2==_groups.end())
        {
          std::ostringstream oss; oss << "No such group in mesh \"" << _name << "\" : " << *it; 
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      fams.insert((*it2).second.begin(),(*it2).second.end());
    }
  std::vector<std::string> fams2(fams.begin(),fams.end());
  return getNodeFamiliesArr(fams2,renum);
}

DataArrayInt *MEDFileUMesh::getNodeFamilyArr(const char *fam, bool renum) const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> tmp(1);
  tmp[0]=fam;
  DataArrayInt *ret=getNodeFamiliesArr(tmp,renum);
  ret->setName(fam);
  return ret;
}

DataArrayInt *MEDFileUMesh::getNodeFamiliesArr(const std::vector<std::string>& fams, bool renum) const throw(INTERP_KERNEL::Exception)
{
  std::vector<int> famIds=getFamiliesIds(fams);
  DataArrayInt *da=_fam_coords->getIdsEqualList(famIds);
  if(renum)
    return MEDFileUMeshSplitL1::renumber(_num_coords,da);
  return da;
}

MEDCouplingUMesh *MEDFileUMesh::getMeshAtRank(int meshDimRelToMax, bool renum) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileUMeshSplitL1 *l1=getMeshAtLevSafe(meshDimRelToMax);
  return l1->getWholeMesh(renum);
}

MEDCouplingUMesh *MEDFileUMesh::getRank0Mesh(bool renum) const throw(INTERP_KERNEL::Exception)
{
  return getMeshAtRank(0,renum);
}

MEDCouplingUMesh *MEDFileUMesh::getRankM1Mesh(bool renum) const throw(INTERP_KERNEL::Exception)
{
  return getMeshAtRank(-1,renum);
}

MEDCouplingUMesh *MEDFileUMesh::getRankM2Mesh(bool renum) const throw(INTERP_KERNEL::Exception)
{
  return getMeshAtRank(-2,renum);
}

MEDCouplingUMesh *MEDFileUMesh::getRankM3Mesh(bool renum) const throw(INTERP_KERNEL::Exception)
{
  return getMeshAtRank(-3,renum);
}

const MEDFileUMeshSplitL1 *MEDFileUMesh::getMeshAtLevSafe(int meshDimRelToMax) const throw(INTERP_KERNEL::Exception)
{
  int tracucedRk=-meshDimRelToMax;
  if(tracucedRk>=(int)_ms.size())
    throw INTERP_KERNEL::Exception("Invalid mesh dim relative to max given ! To low !");
  if((const MEDFileUMeshSplitL1 *)_ms[tracucedRk]==0)
    throw INTERP_KERNEL::Exception("On specified lev (or entity) no cells exists !");
  return _ms[tracucedRk];
}

void MEDFileUMesh::checkMeshDimCoherency(int meshDim, int meshDimRelToMax) const throw(INTERP_KERNEL::Exception)
{
  if(-meshDimRelToMax>=(int)_ms.size())
    throw INTERP_KERNEL::Exception("MEDFileUMesh::checkMeshDimCoherency : The meshdim of mesh is not managed by 'this' !");
  int i=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::const_iterator it=_ms.begin();it!=_ms.end();it++,i++)
    {
      if(((const MEDFileUMeshSplitL1*) (*it))!=0)
        {
          int ref=(*it)->getMeshDimension();
          if(ref+i!=meshDim+meshDimRelToMax)
            throw INTERP_KERNEL::Exception("MEDFileUMesh::checkMeshDimCoherency : no coherency between levels !");
        }
    }
}

void MEDFileUMesh::setCoords(DataArrayDouble *coords) throw(INTERP_KERNEL::Exception)
{
  coords->checkAllocated();
  int nbOfTuples=coords->getNumberOfTuples();
  _coords=coords;
  coords->incrRef();
  _fam_coords=DataArrayInt::New();
  _fam_coords->alloc(nbOfTuples,1);
  _fam_coords->fillWithZero();
}

void MEDFileUMesh::addNodeGroup(const std::vector<int>& ids) throw(INTERP_KERNEL::Exception)
{
  const DataArrayDouble *coords=_coords;
  if(!coords)
    throw INTERP_KERNEL::Exception("addNodeGroup : no coords set !");
  DataArrayInt *sub=_fam_coords->selectByTupleIdSafe(&ids[0],&ids[0]+ids.size());
  std::set<int> ssub(sub->getConstPointer(),sub->getConstPointer()+sub->getNumberOfTuples());
  
}

void MEDFileUMesh::setFamilyNameAttachedOnId(int id, const std::string& newFamName) throw(INTERP_KERNEL::Exception)
{
  std::string oldName=getFamilyNameGivenId(id);
  _families.erase(oldName);
  _families[newFamName]=id;
}

void MEDFileUMesh::setMeshAtRank(int meshDimRelToMax, MEDCouplingUMesh *m) throw(INTERP_KERNEL::Exception)
{
  std::vector<int> levSet=getNonEmptyLevels();
  if(std::find(levSet.begin(),levSet.end(),meshDimRelToMax)==levSet.end())
    {
      if((DataArrayDouble *)_coords==0)
        {
          DataArrayDouble *c=m->getCoords();
          if(c)
            c->incrRef();
          _coords=c;
        }
      else
        {
          if(m->getCoords()!=_coords)
            throw INTERP_KERNEL::Exception("MEDFileUMesh::setMeshAtRank : Invalid Given Mesh ! The coordinates are not the same ! try to use tryToShareSameCoords !");
          int sz=(-meshDimRelToMax)+1;
          if(sz>=(int)_ms.size())
            _ms.resize(sz);
          checkMeshDimCoherency(m->getMeshDimension(),sz);
          _ms[sz]=new MEDFileUMeshSplitL1(m);
        }
    }
}

void MEDFileUMesh::setGroupsFromScratch(int meshDimRelToMax, const std::vector<MEDCouplingUMesh *>& ms) throw(INTERP_KERNEL::Exception)
{
  if(ms.empty())
    throw INTERP_KERNEL::Exception("MEDFileUMesh::setGroupsFromScratch : expecting a non empty vector !");
  int sz=(-meshDimRelToMax)+1;
  if(sz>=(int)_ms.size())
    _ms.resize(sz);
  checkMeshDimCoherency(ms[0]->getMeshDimension(),meshDimRelToMax);
  DataArrayDouble *coo=checkMultiMesh(ms);
  if((DataArrayDouble *)_coords==0)
    {
      coo->incrRef();
      _coords=coo;
    }
  else
    if((DataArrayDouble *)_coords!=coo)
      throw INTERP_KERNEL::Exception("MEDFileUMesh::setGroupsFromScratch : coordinates mismatches !");
  _ms[-meshDimRelToMax]->setGroupsFromScratch(ms,_families,_groups);
}

void MEDFileUMesh::setGroupsOnSetMesh(int meshDimRelToMax, const std::vector<MEDCouplingUMesh *>& ms) throw(INTERP_KERNEL::Exception)
{
  if(ms.empty())
    throw INTERP_KERNEL::Exception("MEDFileUMesh::setGroupsOnSetMesh : expecting a non empty vector !");
  
}

DataArrayDouble *MEDFileUMesh::checkMultiMesh(const std::vector<MEDCouplingUMesh *>& ms) const throw(INTERP_KERNEL::Exception)
{
  return 0;
}
