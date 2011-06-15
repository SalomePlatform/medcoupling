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

#include <limits>

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

MEDFileUMesh::MEDFileUMesh():_too_long_str(0),_zipconn_pol(2)
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
    computeRevNum();
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
  if(!existsFamily(0))
    (const_cast<MEDFileUMesh *>(this))->addFamily(0,"FAMILLE_ZERO");
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
  int mdim=getMeshDimension();
  MEDmaaCr(fid,maa,spaceDim,MED_NON_STRUCTURE,desc);//spaceDim is false but to make reader happy for 3DSurf and 2DCurve meshes !
  MEDdimEspaceCr(fid,maa,spaceDim);
  MEDFileUMeshL2::writeCoords(fid,maa,_coords,_fam_coords,_num_coords);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::const_iterator it=_ms.begin();it!=_ms.end();it++)
    if((const MEDFileUMeshSplitL1 *)(*it)!=0)
      (*it)->write(fid,maa,mdim);
  MEDFileUMeshL2::writeFamiliesAndGrps(fid,maa,_families,_groups,_too_long_str);
  MEDfermer(fid);
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

std::vector<int> MEDFileUMesh::getNonEmptyLevelsExt() const
{
  std::vector<int> ret0=getNonEmptyLevels();
  if((const DataArrayDouble *) _coords)
    {
      std::vector<int> ret(ret0.size()+1);
      ret[0]=1;
      std::copy(ret0.begin(),ret0.end(),ret.begin()+1);
      return ret;
    }
  return ret0;
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

int MEDFileUMesh::getMaxFamilyId() const throw(INTERP_KERNEL::Exception)
{
  if(_families.empty())
    throw INTERP_KERNEL::Exception("MEDFileUMesh::getMaxFamilyId : no families set !");
  int ret=-std::numeric_limits<int>::max();
  for(std::map<std::string,int>::const_iterator it=_families.begin();it!=_families.end();it++)
    {
      ret=std::max((*it).second,ret);
    }
  return ret;
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

int MEDFileUMesh::getSizeAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception)
{
  if(meshDimRelToMaxExt==1)
    {
      if(!((const DataArrayDouble *)_coords))
        throw INTERP_KERNEL::Exception("MEDFileUMesh::getSizeAtLevel : no coordinates specified !");
      return _coords->getNumberOfTuples();
    }
  return getMeshAtLevSafe(meshDimRelToMaxExt)->getSize();
}

const DataArrayInt *MEDFileUMesh::getFamilyFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception)
{
  if(meshDimRelToMaxExt==1)
    {
      if(!((const DataArrayInt *)_fam_coords))
        throw INTERP_KERNEL::Exception("MEDFileUMesh::getFamilyFieldAtLevel : no coordinates specified !");
      return _fam_coords;
    }
  const MEDFileUMeshSplitL1 *l1=getMeshAtLevSafe(meshDimRelToMaxExt);
  return l1->getFamilyField();
}

const DataArrayInt *MEDFileUMesh::getNumberFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception)
{
  if(meshDimRelToMaxExt==1)
    {
      if(!((const DataArrayInt *)_num_coords))
        throw INTERP_KERNEL::Exception("MEDFileUMesh::getNumberFieldAtLevel : no coordinates renum specified !");
      return _num_coords;
    }
  const MEDFileUMeshSplitL1 *l1=getMeshAtLevSafe(meshDimRelToMaxExt);
  return l1->getNumberField();
}

const DataArrayInt *MEDFileUMesh::getRevNumberFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception)
{
  if(meshDimRelToMaxExt==1)
    {
      if(!((const DataArrayInt *)_num_coords))
        throw INTERP_KERNEL::Exception("MEDFileUMesh::getRevNumberFieldAtLevel : no coordinates renum specified !");
      return _rev_num_coords;
    }
  const MEDFileUMeshSplitL1 *l1=getMeshAtLevSafe(meshDimRelToMaxExt);
  return l1->getRevNumberField();
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

std::vector<std::string> MEDFileUMesh::getGroupsOnFamily(const char *name) const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> ret;
  for(std::map<std::string, std::vector<std::string> >::const_iterator it1=_groups.begin();it1!=_groups.end();it1++)
    {
      for(std::vector<std::string>::const_iterator it2=(*it1).second.begin();it2!=(*it1).second.end();it2++)
        if((*it2)==name)
          {
            ret.push_back((*it1).first);
            break;
          }
    }
  return ret;
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

MEDCouplingUMesh *MEDFileUMesh::getGroup(int meshDimRelToMaxExt, const char *grp, bool renum) const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> tmp(1);
  tmp[0]=grp;
  MEDCouplingUMesh *ret=getGroups(meshDimRelToMaxExt,tmp,renum);
  ret->setName(grp);
  return ret;
}

DataArrayInt *MEDFileUMesh::getGroupArr(int meshDimRelToMaxExt, const char *grp, bool renum) const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> tmp(1);
  tmp[0]=grp;
  DataArrayInt *ret=getGroupsArr(meshDimRelToMaxExt,tmp,renum);
  ret->setName(grp);
  return ret;
}

MEDCouplingUMesh *MEDFileUMesh::getGroups(int meshDimRelToMaxExt, const std::vector<std::string>& grps, bool renum) const throw(INTERP_KERNEL::Exception)
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
  return getFamilies(meshDimRelToMaxExt,fams2,renum);
}

DataArrayInt *MEDFileUMesh::getGroupsArr(int meshDimRelToMaxExt, const std::vector<std::string>& grps, bool renum) const throw(INTERP_KERNEL::Exception)
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
  return getFamiliesArr(meshDimRelToMaxExt,fams2,renum);
}

MEDCouplingUMesh *MEDFileUMesh::getFamily(int meshDimRelToMaxExt, const char *fam, bool renum) const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> tmp(1);
  tmp[0]=fam;
  MEDCouplingUMesh *ret=getFamilies(meshDimRelToMaxExt,tmp,renum);
  ret->setName(fam);
  return ret;
}

DataArrayInt *MEDFileUMesh::getFamilyArr(int meshDimRelToMaxExt, const char *fam, bool renum) const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> tmp(1);
  tmp[0]=fam;
  DataArrayInt *ret=getFamiliesArr(meshDimRelToMaxExt,tmp,renum);
  ret->setName(fam);
  return ret;
}

MEDCouplingUMesh *MEDFileUMesh::getFamilies(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum) const throw(INTERP_KERNEL::Exception)
{
  if(meshDimRelToMaxExt==1)
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr=getFamiliesArr(1,fams,renum);
      MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> ret=MEDCouplingUMesh::New();
      MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> c=_coords->selectByTupleId(arr->getConstPointer(),arr->getConstPointer()+arr->getNbOfElems());
      ret->setCoords(c);
      ret->incrRef();
      return ret;
    }
  std::vector<int> famIds=getFamiliesIds(fams);
  const MEDFileUMeshSplitL1 *l1=getMeshAtLevSafe(meshDimRelToMaxExt);
  return l1->getFamilyPart(famIds,renum);
}

DataArrayInt *MEDFileUMesh::getFamiliesArr(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum) const throw(INTERP_KERNEL::Exception)
{
  std::vector<int> famIds=getFamiliesIds(fams);
  if(meshDimRelToMaxExt==1)
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da=_fam_coords->getIdsEqualList(famIds);
      return MEDFileUMeshSplitL1::Renumber(_num_coords,da);
    }
  const MEDFileUMeshSplitL1 *l1=getMeshAtLevSafe(meshDimRelToMaxExt);
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
  return getGroupsArr(1,grps,renum);
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
  return getFamiliesArr(1,fams,renum);
}

MEDCouplingUMesh *MEDFileUMesh::getMeshAtLevel(int meshDimRelToMaxExt, bool renum) const throw(INTERP_KERNEL::Exception)
{
  if(meshDimRelToMaxExt==1)
    {
      if(!renum)
        {
          MEDCouplingUMesh *umesh=MEDCouplingUMesh::New();
          MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> cc=_coords->deepCpy();
          umesh->setCoords(cc);
          return umesh;
        }
    }
  const MEDFileUMeshSplitL1 *l1=getMeshAtLevSafe(meshDimRelToMaxExt);
  return l1->getWholeMesh(renum);
}

MEDCouplingUMesh *MEDFileUMesh::getLevel0Mesh(bool renum) const throw(INTERP_KERNEL::Exception)
{
  return getMeshAtLevel(0,renum);
}

MEDCouplingUMesh *MEDFileUMesh::getLevelM1Mesh(bool renum) const throw(INTERP_KERNEL::Exception)
{
  return getMeshAtLevel(-1,renum);
}

MEDCouplingUMesh *MEDFileUMesh::getLevelM2Mesh(bool renum) const throw(INTERP_KERNEL::Exception)
{
  return getMeshAtLevel(-2,renum);
}

MEDCouplingUMesh *MEDFileUMesh::getLevelM3Mesh(bool renum) const throw(INTERP_KERNEL::Exception)
{
  return getMeshAtLevel(-3,renum);
}

bool MEDFileUMesh::existsFamily(int famId) const
{
  for(std::map<std::string,int>::const_iterator it2=_families.begin();it2!=_families.end();it2++)
    if((*it2).second==famId)
      return true;
  return false;
}

bool MEDFileUMesh::existsFamily(const char *familyName) const
{
  std::string fname(familyName);
  return _families.find(fname)!=_families.end();
}

const MEDFileUMeshSplitL1 *MEDFileUMesh::getMeshAtLevSafe(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception)
{
  if(meshDimRelToMaxExt==1)
    throw INTERP_KERNEL::Exception("Dimension request is invalid : asking for node level (1) !");
  if(meshDimRelToMaxExt>1)
    throw INTERP_KERNEL::Exception("Dimension request is invalid (>1) !");
  int tracucedRk=-meshDimRelToMaxExt;
  if(tracucedRk>=(int)_ms.size())
    throw INTERP_KERNEL::Exception("Invalid mesh dim relative to max given ! To low !");
  if((const MEDFileUMeshSplitL1 *)_ms[tracucedRk]==0)
    throw INTERP_KERNEL::Exception("On specified lev (or entity) no cells exists !");
  return _ms[tracucedRk];
}

MEDFileUMeshSplitL1 *MEDFileUMesh::getMeshAtLevSafe(int meshDimRelToMaxExt) throw(INTERP_KERNEL::Exception)
{
   if(meshDimRelToMaxExt==1)
    throw INTERP_KERNEL::Exception("Dimension request is invalid : asking for node level (1) !");
  if(meshDimRelToMaxExt>1)
    throw INTERP_KERNEL::Exception("Dimension request is invalid (>1) !");
  int tracucedRk=-meshDimRelToMaxExt;
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
          if(ref+i!=meshDim-meshDimRelToMax)
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

void MEDFileUMesh::setGroupsAtLevel(int meshDimRelToMaxExt, const std::vector<const DataArrayInt *>& grps, bool renum) throw(INTERP_KERNEL::Exception)
{
  if(grps.empty())
    return ;
  std::set<std::string> grpsName;
  std::vector<std::string> grpsName2(grps.size());
  int i=0;

  for(std::vector<const DataArrayInt *>::const_iterator it=grps.begin();it!=grps.end();it++,i++)
    {
      grpsName.insert((*it)->getName());
      grpsName2[i]=(*it)->getName();
    }
  if(grpsName.size()!=grps.size())
    throw INTERP_KERNEL::Exception("MEDFileUMesh::setGroupsAtLevel : groups name must be different each other !");
  if(grpsName.find(std::string(""))!=grpsName.end())
    throw INTERP_KERNEL::Exception("MEDFileUMesh::setGroupsAtLevel : groups name must be different empty string !");
  int sz=getSizeAtLevel(meshDimRelToMaxExt);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> fam;
  std::vector< std::vector<int> > fidsOfGroups;
  if(!renum)
    {
      fam=DataArrayInt::MakePartition(grps,sz,fidsOfGroups);
    }
  else
    {
      std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > grps2(grps.size());
      for(unsigned int i=0;i<grps.size();i++)
        {
          grps2[i]=MEDFileUMeshSplitL1::Renumber(getRevNumberFieldAtLevel(meshDimRelToMaxExt),grps[i]);
          grps2[i]->setName(grps[i]->getName().c_str());
        }
      std::vector<const DataArrayInt *> grps3(grps2.begin(),grps2.end());
      fam=DataArrayInt::MakePartition(grps3,sz,fidsOfGroups);
    }
  int offset=1;
  if(!_families.empty())
    offset=getMaxFamilyId()+1;
  TranslateFamilyIds(offset,fam,fidsOfGroups);
  std::set<int> ids=fam->getDifferentValues();
  appendFamilyEntries(ids,fidsOfGroups,grpsName2);
  setFamilyArr(meshDimRelToMaxExt,fam);
}

void MEDFileUMesh::eraseGroupsAtLevel(int meshDimRelToMaxExt) throw(INTERP_KERNEL::Exception)
{
  if(meshDimRelToMaxExt==1)
    {
      if((DataArrayInt *)_fam_coords)
        _fam_coords->fillWithZero();
      return ;
    }
  MEDFileUMeshSplitL1 *l1=getMeshAtLevSafe(meshDimRelToMaxExt);
  l1->eraseFamilyField();
  optimizeFamilies();
}

void MEDFileUMesh::optimizeFamilies() throw(INTERP_KERNEL::Exception)
{
  std::vector<int> levs=getNonEmptyLevelsExt();
  std::set<int> allFamsIds;
  for(std::vector<int>::const_iterator it=levs.begin();it!=levs.end();it++)
    {
      const DataArrayInt *ffield=getFamilyFieldAtLevel(*it);
      std::set<int> ids=ffield->getDifferentValues();
      std::set<int> res;
      std::set_union(ids.begin(),ids.end(),allFamsIds.begin(),allFamsIds.end(),std::inserter(res,res.begin()));
      allFamsIds=res;
    }
  std::set<std::string> famNamesToKill;
  for(std::map<std::string,int>::const_iterator it=_families.begin();it!=_families.end();it++)
    {
      if(allFamsIds.find((*it).second)!=allFamsIds.end())
        famNamesToKill.insert((*it).first);
    }
  for(std::set<std::string>::const_iterator it=famNamesToKill.begin();it!=famNamesToKill.end();it++)
    _families.erase(*it);
  std::vector<std::string> grpNamesToKill;
  for(std::map<std::string, std::vector<std::string> >::iterator it=_groups.begin();it!=_groups.end();it++)
    {
      std::vector<std::string> tmp;
      for(std::vector<std::string>::const_iterator it2=(*it).second.begin();it2!=(*it).second.end();it2++)
        {
          if(famNamesToKill.find(*it2)==famNamesToKill.end())
            tmp.push_back(*it2);
        }
      if(!tmp.empty())
        (*it).second=tmp;
      else
        tmp.push_back((*it).first);
    }
  for(std::vector<std::string>::const_iterator it=grpNamesToKill.begin();it!=grpNamesToKill.end();it++)
    _groups.erase(*it);
}

void MEDFileUMesh::setFamilyField(DataArrayInt *arr, const std::vector< std::vector< int > > &userfids, const std::vector<std::string>& grpNames, bool renum) throw(INTERP_KERNEL::Exception)
{
  
}

void MEDFileUMesh::addNodeGroup(const std::string& name, const std::vector<int>& ids) throw(INTERP_KERNEL::Exception)
{
  const DataArrayDouble *coords=_coords;
  if(!coords)
    throw INTERP_KERNEL::Exception("addNodeGroup : no coords set !");
  DataArrayInt *sub=_fam_coords->selectByTupleIdSafe(&ids[0],&ids[0]+ids.size());
  std::set<int> ssub(sub->getConstPointer(),sub->getConstPointer()+sub->getNumberOfTuples());
  
}

void MEDFileUMesh::copyFamGrpMapsFrom(const MEDFileUMesh& other)
{
  _groups=other._groups;
  _families=other._families;
}

void MEDFileUMesh::setFamilyInfo(const std::map<std::string,int>& info)
{
  _families=info;
}

void MEDFileUMesh::setGroupInfo(const std::map<std::string, std::vector<std::string> >&info)
{
  _groups=info;
}

void MEDFileUMesh::setFamilyNameAttachedOnId(int id, const std::string& newFamName) throw(INTERP_KERNEL::Exception)
{
  std::string oldName=getFamilyNameGivenId(id);
  _families.erase(oldName);
  _families[newFamName]=id;
}

void MEDFileUMesh::setMeshAtLevel(int meshDimRelToMax, MEDCouplingUMesh *m) throw(INTERP_KERNEL::Exception)
{
  setMeshAtLevelGen(meshDimRelToMax,m,true);
}

void MEDFileUMesh::setMeshAtLevelOld(int meshDimRelToMax, MEDCouplingUMesh *m) throw(INTERP_KERNEL::Exception)
{
  setMeshAtLevelGen(meshDimRelToMax,m,false);
}

void MEDFileUMesh::setMeshAtLevelGen(int meshDimRelToMax, MEDCouplingUMesh *m, bool newOrOld) throw(INTERP_KERNEL::Exception)
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
            throw INTERP_KERNEL::Exception("MEDFileUMesh::setMeshAtLevel : Invalid Given Mesh ! The coordinates are not the same ! try to use tryToShareSameCoords !");
          int sz=(-meshDimRelToMax)+1;
          if(sz>=(int)_ms.size())
            _ms.resize(sz);
          checkMeshDimCoherency(m->getMeshDimension(),meshDimRelToMax);
          _ms[sz-1]=new MEDFileUMeshSplitL1(m,newOrOld);
        }
    }
}

void MEDFileUMesh::setGroupsFromScratch(int meshDimRelToMax, const std::vector<const MEDCouplingUMesh *>& ms) throw(INTERP_KERNEL::Exception)
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
  std::vector<DataArrayInt *> corr;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m=MEDCouplingUMesh::FuseUMeshesOnSameCoords(ms,_zipconn_pol,corr);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > corr3(corr.begin(),corr.end());
  setMeshAtLevel(meshDimRelToMax,m);
  std::vector<const DataArrayInt *> corr2(corr.begin(),corr.end());
  setGroupsAtLevel(meshDimRelToMax,corr2,true);
}

void MEDFileUMesh::setGroupsOnSetMesh(int meshDimRelToMax, const std::vector<const MEDCouplingUMesh *>& ms, bool renum) throw(INTERP_KERNEL::Exception)
{
  if(ms.empty())
    throw INTERP_KERNEL::Exception("MEDFileUMesh::setGroupsOnSetMesh : expecting a non empty vector !");
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
      throw INTERP_KERNEL::Exception("MEDFileUMesh::setGroupsOnSetMesh : coordinates mismatches !");
  MEDCouplingUMesh *m=getMeshAtLevel(meshDimRelToMax,renum);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > corr(ms.size());
  int i=0;
  for(std::vector<const MEDCouplingUMesh *>::const_iterator it=ms.begin();it!=ms.end();it++,i++)
    {
      DataArrayInt *arr=0;
      bool test=m->areCellsIncludedIn(*it,_zipconn_pol,arr);
      corr[i]=arr;
      if(!test)
        {
          std::ostringstream oss; oss << "MEDFileUMesh::setGroupsOnSetMesh : mesh #" << i << " is not part of whole mesh !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  std::vector<const DataArrayInt *> corr2(corr.begin(),corr.end());
  setGroupsAtLevel(meshDimRelToMax,corr2,renum);
}

DataArrayDouble *MEDFileUMesh::checkMultiMesh(const std::vector<const MEDCouplingUMesh *>& ms) const throw(INTERP_KERNEL::Exception)
{
  const DataArrayDouble *ret=ms[0]->getCoords();
  int mdim=ms[0]->getMeshDimension();
  for(unsigned int i=1;i<ms.size();i++)
    {
      ms[i]->checkCoherency();
      if(ms[i]->getCoords()!=ret)
        throw INTERP_KERNEL::Exception("MEDFileUMesh::checkMultiMesh : meshes must share the same coords !");
      if(ms[i]->getMeshDimension()!=mdim)
        throw INTERP_KERNEL::Exception("MEDFileUMesh::checkMultiMesh : meshes have not same mesh dimension !");
    }
  return const_cast<DataArrayDouble *>(ret);
}

void MEDFileUMesh::setFamilyArr(int meshDimRelToMaxExt, DataArrayInt *famArr)
{
  if(meshDimRelToMaxExt==1)
    {
      famArr->incrRef();
      _fam_coords=famArr;
      return ;
    }
  if(meshDimRelToMaxExt>1)
    throw INTERP_KERNEL::Exception("MEDFileUMesh::setFamilyArr : Dimension request is invalid (>1) !");
  int traducedRk=-meshDimRelToMaxExt;
  if(traducedRk>=(int)_ms.size())
    throw INTERP_KERNEL::Exception("Invalid mesh dim relative to max given ! To low !");
  if((MEDFileUMeshSplitL1 *)_ms[traducedRk]==0)
    throw INTERP_KERNEL::Exception("On specified lev (or entity) no cells exists !");
  return _ms[traducedRk]->setFamilyArr(famArr);
}

void MEDFileUMesh::setRenumArr(int meshDimRelToMaxExt, DataArrayInt *renumArr)
{
  if(meshDimRelToMaxExt==1)
    {
      if(renumArr)
        renumArr->incrRef();
      _num_coords=renumArr;
      computeRevNum();
      return ;
    }
  if(meshDimRelToMaxExt>1)
    throw INTERP_KERNEL::Exception("MEDFileUMesh::setRenumArr : Dimension request is invalid (>1) !");
  int traducedRk=-meshDimRelToMaxExt;
  if(traducedRk>=(int)_ms.size())
    throw INTERP_KERNEL::Exception("Invalid mesh dim relative to max given ! To low !");
  if((MEDFileUMeshSplitL1 *)_ms[traducedRk]==0)
    throw INTERP_KERNEL::Exception("On specified lev (or entity) no cells exists !");
  return _ms[traducedRk]->setRenumArr(renumArr);
}

void MEDFileUMesh::TranslateFamilyIds(int offset, DataArrayInt *famArr, std::vector< std::vector<int> >& famIdsPerGrp)
{
  famArr->applyLin(1,offset,0);
  for(std::vector< std::vector<int> >::iterator it1=famIdsPerGrp.begin();it1!=famIdsPerGrp.end();it1++)
    std::transform((*it1).begin(),(*it1).end(),(*it1).begin(),std::bind2nd(std::plus<int>(),offset));
}

/*!
 * This method append into '_families' attribute the families whose ids are in 'famIds'. Warning 'famIds' are expected to be ids
 * not in '_families'. Groups information are given in parameters in order to give to families representative names.
 * For the moment, the two last input parameters are not taken into account.
 */
void MEDFileUMesh::appendFamilyEntries(const std::set<int>& famIds, const std::vector< std::vector<int> >& fidsOfGrps, const std::vector<std::string>& grpNames)
{
  std::map<int,std::string> famInv;
  for(std::set<int>::const_iterator it=famIds.begin();it!=famIds.end();it++)
    {
      std::ostringstream oss;
      oss << "Family_" << (*it);
      _families[oss.str()]=(*it);
      famInv[*it]=oss.str();
    }
  int i=0;
  for(std::vector< std::vector<int> >::const_iterator it1=fidsOfGrps.begin();it1!=fidsOfGrps.end();it1++,i++)
    {
      for(std::vector<int>::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
        {
          _groups[grpNames[i]].push_back(famInv[*it2]);
        }
    }
}

/*!
 * This method appends a new entry in _families attribute. An exception is thrown if either the famId is already
 * kept by an another familyName. An exception is thrown if name 'familyName' is alreadyset with a different 'famId'.
 */
void MEDFileUMesh::addFamily(int famId, const char *familyName) throw(INTERP_KERNEL::Exception)
{
  std::string fname(familyName);
  std::map<std::string,int>::const_iterator it=_families.find(fname);
  if(it==_families.end())
    {
       for(std::map<std::string,int>::const_iterator it2=_families.begin();it2!=_families.end();it2++)
         if((*it2).second==famId)
           {
             std::ostringstream oss;
             oss << "MEDFileUMesh::addFamily : Family \"" << (*it2).first << "\" already exists with specified id : " << famId << " !";
             throw INTERP_KERNEL::Exception(oss.str().c_str());
           }
       _families[fname]=famId;
    }
  else
    {
      if((*it).second!=famId)
        {
          std::ostringstream oss;
          oss << "MEDFileUMesh::addFamily : Family \"" << fname << "\" already exists but has id set to " << (*it).second << " different from asked famId " << famId << " !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
}

void MEDFileUMesh::computeRevNum() const
{
  if((const DataArrayInt *)_num_coords)
    {
      int pos;
      int maxValue=_num_coords->getMaxValue(pos);
      _rev_num_coords=_num_coords->invertArrayN2O2O2N(maxValue+1);
    }
}
