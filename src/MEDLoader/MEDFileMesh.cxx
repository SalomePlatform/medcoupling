// Copyright (C) 2007-2011  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

#include "MEDFileMesh.hxx"
#include "MEDFileUtilities.hxx"
#include "MEDLoader.hxx"
#include "MEDLoaderBase.hxx"

#include "MEDCouplingUMesh.hxx"

#include "InterpKernelAutoPtr.hxx"

#include <limits>
#include <cmath>

using namespace ParaMEDMEM;

const char MEDFileMesh::DFT_FAM_NAME[]="FAMILLE_ZERO";

MEDFileMesh::MEDFileMesh():_order(-1),_iteration(-1),_time(0.)
{
}

MEDFileMesh *MEDFileMesh::New(const char *fileName) throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> ms=MEDLoader::GetMeshNames(fileName);
  if(ms.empty())
    {
      std::ostringstream oss; oss << "MEDFileMesh::New : no meshes in file \"" << fileName << "\" !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  MEDFileUtilities::CheckFileForRead(fileName);
  ParaMEDMEM::MEDCouplingMeshType meshType;
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  int dt,it;
  std::string dummy2;
  MEDFileMeshL2::GetMeshIdFromName(fid,ms.front().c_str(),meshType,dt,it,dummy2);
  switch(meshType)
    {
    case UNSTRUCTURED:
      {
        MEDCouplingAutoRefCountObjectPtr<MEDFileUMesh> ret=MEDFileUMesh::New();
        ret->loadUMeshFromFile(fid,ms.front().c_str(),dt,it);
        ret->incrRef();
        return (MEDFileUMesh *)ret;
      }
    case CARTESIAN:
      {
        MEDCouplingAutoRefCountObjectPtr<MEDFileCMesh> ret=MEDFileCMesh::New();
        ret->loadCMeshFromFile(fid,ms.front().c_str(),dt,it);
        ret->incrRef();
        return (MEDFileCMesh *)ret;
      }
    default:
      {
        std::ostringstream oss; oss << "MEDFileMesh::New : MED file exists and has mesh '" << ms.front() << "' exists but unsupported type yet !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    }
}

MEDFileMesh *MEDFileMesh::New(const char *fileName, const char *mName, int dt, int it) throw(INTERP_KERNEL::Exception)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  ParaMEDMEM::MEDCouplingMeshType meshType;
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  int dummy0,dummy1;
  std::string dummy2;
  MEDFileMeshL2::GetMeshIdFromName(fid,mName,meshType,dummy0,dummy1,dummy2);
  switch(meshType)
    {
    case UNSTRUCTURED:
      {
        MEDCouplingAutoRefCountObjectPtr<MEDFileUMesh> ret=MEDFileUMesh::New();
        ret->loadUMeshFromFile(fid,mName,dt,it);
        ret->incrRef();
        return (MEDFileUMesh *)ret;
      }
    case CARTESIAN:
      {
        MEDCouplingAutoRefCountObjectPtr<MEDFileCMesh> ret=MEDFileCMesh::New();
        ret->loadCMeshFromFile(fid,mName,dt,it);
        ret->incrRef();
        return (MEDFileCMesh *)ret;
      }
    default:
      {
        std::ostringstream oss; oss << "MEDFileMesh::New : MED file exists and has mesh '" << mName << "' exists but unsupported type yet !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    }
}

bool MEDFileMesh::isEqual(const MEDFileMesh *other, double eps, std::string& what) const
{
  if(_order!=other->_order)
    {
      what="Orders differ !";
      return false;
    }
  if(_iteration!=other->_iteration)
    {
      what="Iterations differ !";
      return false;
    }
  if(fabs(_time-other->_time)>eps)
    {
      what="Time values differ !";
      return false;
    }
  if(_dt_unit!=other->_dt_unit)
    {
      what="Time units differ !";
      return false;
    }
  if(_name!=other->_name)
    {
      what="Names differ !";
      return false;
    }
  if(_univ_name!=other->_univ_name)
    {
      what="Univ names differ !";
      return false;
    }
  if(_desc_name!=other->_desc_name)
    {
      what="Description names differ !";
      return false;
    }
  if(!areGrpsEqual(other,what))
    return false;
  if(!areFamsEqual(other,what))
    return false;
  return true;
}

void MEDFileMesh::clearNonDiscrAttributes() const
{
  
}

void MEDFileMesh::copyFamGrpMapsFrom(const MEDFileMesh& other)
{
  _groups=other._groups;
  _families=other._families;
}

std::vector<std::string> MEDFileMesh::getFamiliesOnGroup(const char *name) const throw(INTERP_KERNEL::Exception)
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

std::vector<std::string> MEDFileMesh::getFamiliesOnGroups(const std::vector<std::string>& grps) const throw(INTERP_KERNEL::Exception)
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
  return fams2;
}

std::vector<int> MEDFileMesh::getFamiliesIdsOnGroup(const char *name) const throw(INTERP_KERNEL::Exception)
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
  return getFamiliesIds((*it).second);
}

/*!
 * This method sets families at a corresponding groups existing or not. If it existed, it is replaced by new 'fams'.
 * Each entry in 'fams' is checked if it is not still existing default id 0 is set.
 */
void MEDFileMesh::setFamiliesOnGroup(const char *name, const std::vector<std::string>& fams) throw(INTERP_KERNEL::Exception)
{
  std::string oname(name);
  _groups[oname]=fams;
  for(std::vector<std::string>::const_iterator it1=fams.begin();it1!=fams.end();it1++)
    {
      std::map<std::string,int>::iterator it2=_families.find(*it1);
      if(it2==_families.end())
        _families[*it1]=0;
    }
}

/*!
 * Behaves as MEDFileMesh::setFamiliesOnGroup, except that if there is presence of a family id in 'famIds' not existing an exception is thrown.
 * If several families have same id the first one in lexical order is taken into account.
 */
void MEDFileMesh::setFamiliesIdsOnGroup(const char *name, const std::vector<int>& famIds) throw(INTERP_KERNEL::Exception)
{
  std::string oname(name);
  std::vector<std::string> fams(famIds.size());
  int i=0;
  for(std::vector<int>::const_iterator it1=famIds.begin();it1!=famIds.end();it1++,i++)
    {
      std::string name=getFamilyNameGivenId(*it1);
      fams[i]=name;
    }
  _groups[oname]=fams;
}

std::vector<std::string> MEDFileMesh::getGroupsOnFamily(const char *name) const throw(INTERP_KERNEL::Exception)
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

/*!
 * This method expects that family 'famName' is already existing. If not an exception will be thrown.
 */
void MEDFileMesh::setGroupsOnFamily(const char *famName, const std::vector<std::string>& grps) throw(INTERP_KERNEL::Exception)
{
  std::string fName(famName);
  const std::map<std::string,int>::const_iterator it=_families.find(fName);
  if(it==_families.end())
    {
      std::vector<std::string> fams=getFamiliesNames();
      std::ostringstream oss; oss << "No such familyname \"" << fName << "\" !\nAvailable families are :";
      std::copy(fams.begin(),fams.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  for(std::vector<std::string>::const_iterator it=grps.begin();it!=grps.end();it++)
    {
      std::map< std::string, std::vector<std::string> >::iterator it2=_groups.find(*it);
      if(it2!=_groups.end())
        (*it2).second.push_back(fName);
      else
        {
          std::vector<std::string> grps(1,fName);
          _groups[*it]=grps;
        }
    }
}

std::vector<std::string> MEDFileMesh::getGroupsNames() const
{
  std::vector<std::string> ret(_groups.size());
  int i=0;
  for(std::map<std::string, std::vector<std::string> >::const_iterator it=_groups.begin();it!=_groups.end();it++,i++)
    ret[i]=(*it).first;
  return ret;
}

std::vector<std::string> MEDFileMesh::getFamiliesNames() const
{
  std::vector<std::string> ret(_families.size());
  int i=0;
  for(std::map<std::string, int >::const_iterator it=_families.begin();it!=_families.end();it++,i++)
    ret[i]=(*it).first;
  return ret;
}

/*!
 * This method scans every families and for each families shared by only one group, the corresponding family takes the same name than the group.
 */
void MEDFileMesh::assignFamilyNameWithGroupName() throw(INTERP_KERNEL::Exception)
{
  std::map<std::string, std::vector<std::string> > groups(_groups);
  std::map<std::string,int> newFams;
  for(std::map<std::string,int>::const_iterator it=_families.begin();it!=_families.end();it++)
    {
      std::vector<std::string> grps=getGroupsOnFamily((*it).first.c_str());
      if(grps.size()==1 && groups[grps[0]].size()==1)
        {
          if(newFams.find(grps[0])!=newFams.end())
            {
              std::ostringstream oss; oss << "MEDFileMesh::assignFamilyNameWithGroupName : Family \"" << grps[0] << "\" already exists !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
          newFams[grps[0]]=(*it).second;
          std::vector<std::string>& grps2=groups[grps[0]];
          std::size_t pos=std::distance(grps2.begin(),std::find(grps2.begin(),grps2.end(),(*it).first));
          grps2[pos]=grps[0];
        }
      else
        {
          if(newFams.find((*it).first)!=newFams.end())
            {
              std::ostringstream oss; oss << "MEDFileMesh::assignFamilyNameWithGroupName : Family \"" << (*it).first << "\" already exists !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
          newFams[(*it).first]=(*it).second;
        }
    }
  _families=newFams;
  _groups=groups;
}

void MEDFileMesh::removeGroup(const char *name) throw(INTERP_KERNEL::Exception)
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

void MEDFileMesh::removeFamily(const char *name) throw(INTERP_KERNEL::Exception)
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

void MEDFileMesh::changeGroupName(const char *oldName, const char *newName) throw(INTERP_KERNEL::Exception)
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
      std::ostringstream oss; oss << "Such groupname \"" << newName << "\" already exists ! Kill it before !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  std::vector<std::string> cpy=(*it).second;
  _groups.erase(it);
  _groups[newName]=cpy;
}

void MEDFileMesh::changeFamilyId(int oldId, int newId) throw(INTERP_KERNEL::Exception)
{
  changeFamilyIdArr(oldId,newId);
  std::map<std::string,int> fam2;
  for(std::map<std::string,int>::const_iterator it=_families.begin();it!=_families.end();it++)
    {
      if((*it).second==oldId)
        fam2[(*it).first]=newId;
      else
        fam2[(*it).first]=(*it).second;
    }
  _families=fam2;
}

void MEDFileMesh::changeFamilyName(const char *oldName, const char *newName) throw(INTERP_KERNEL::Exception)
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
  std::map<std::string, int >::iterator it2=_families.find(nname);
  if(it2!=_families.end())
    {
      std::ostringstream oss; oss << "Such familyname \"" << newName << " already exists ! Kill it before !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  int cpy=(*it).second;
  _families.erase(it);
  _families[newName]=cpy;
}

bool MEDFileMesh::areFamsEqual(const MEDFileMesh *other, std::string& what) const
{
  if(_families==other->_families)
    return true;
  std::map<std::string,int> fam0;
  std::map<std::string,int> fam1;
  for(std::map<std::string,int>::const_iterator it=_families.begin();it!=_families.end();it++)
    if((*it).second!=0)
      fam0[(*it).first]=(*it).second;
  for(std::map<std::string,int>::const_iterator it=other->_families.begin();it!=other->_families.end();it++)
    if((*it).second!=0)
      fam1[(*it).first]=(*it).second;
  return fam0==fam1;
}

bool MEDFileMesh::areGrpsEqual(const MEDFileMesh *other, std::string& what) const
{
  if(_groups==other->_groups)
    return true;
  bool ret=true;
  std::size_t sz=_groups.size();
  if(sz!=other->_groups.size())
    {
      what="Groups differ because not same number !\n";
      ret=false;
    }
  if(ret)
    {
      std::map<std::string, std::vector<std::string> >::const_iterator it1=_groups.begin();
      for(std::size_t i=0;i<sz && ret;i++,it1++)
        {
          std::map<std::string, std::vector<std::string> >::const_iterator it2=other->_groups.find((*it1).first);
          if(it2!=other->_groups.end())
            {
              std::set<std::string> s1((*it1).second.begin(),(*it1).second.end());
              std::set<std::string> s2((*it2).second.begin(),(*it2).second.end());
              ret=(s1==s2);
            }
          else
            {
              ret=false;
              what="A group in first mesh exists not in other !\n";
            }
        }
    }
  if(!ret)
    {
      std::ostringstream oss; oss << "Groups description differs :\n";
      oss << "First group description :\n";
      for(std::map<std::string, std::vector<std::string> >::const_iterator it=_groups.begin();it!=_groups.end();it++)
        {
          oss << " Group \"" << (*it).first << "\" on following families :\n";
          for(std::vector<std::string>::const_iterator it2=(*it).second.begin();it2!=(*it).second.end();it2++)
        oss << "    \"" << *it2 << "\n";
        }
      oss << "Second group description :\n";
      for(std::map<std::string, std::vector<std::string> >::const_iterator it=other->_groups.begin();it!=other->_groups.end();it++)
        {
          oss << " Group \"" << (*it).first << "\" on following families :\n";
          for(std::vector<std::string>::const_iterator it2=(*it).second.begin();it2!=(*it).second.end();it2++)
            oss << "    \"" << *it2 << "\n";
        }
      what+=oss.str();
    }
  return ret;
}

bool MEDFileMesh::existsFamily(int famId) const
{
  for(std::map<std::string,int>::const_iterator it2=_families.begin();it2!=_families.end();it2++)
    if((*it2).second==famId)
      return true;
  return false;
}

bool MEDFileMesh::existsFamily(const char *familyName) const
{
  std::string fname(familyName);
  return _families.find(fname)!=_families.end();
}

void MEDFileMesh::setFamilyId(const char *familyName, int id)
{
  std::string fname(familyName);
  _families[fname]=id;
}

/*!
 * This method appends a new entry in _families attribute. An exception is thrown if either the famId is already
 * kept by an another familyName. An exception is thrown if name 'familyName' is alreadyset with a different 'famId'.
 */
void MEDFileMesh::addFamily(const char *familyName, int famId) throw(INTERP_KERNEL::Exception)
{
  std::string fname(familyName);
  std::map<std::string,int>::const_iterator it=_families.find(fname);
  if(it==_families.end())
    {
       for(std::map<std::string,int>::const_iterator it2=_families.begin();it2!=_families.end();it2++)
         if((*it2).second==famId)
           {
             std::ostringstream oss;
             oss << "MEDFileMesh::addFamily : Family \"" << (*it2).first << "\" already exists with specified id : " << famId << " !";
             throw INTERP_KERNEL::Exception(oss.str().c_str());
           }
       _families[fname]=famId;
    }
  else
    {
      if((*it).second!=famId)
        {
          std::ostringstream oss;
          oss << "MEDFileMesh::addFamily : Family \"" << fname << "\" already exists but has id set to " << (*it).second << " different from asked famId " << famId << " !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
}

void MEDFileMesh::addGrpOnFamily(const char *grpName, const char *famName) throw(INTERP_KERNEL::Exception)
{
  std::string grpn(grpName);
  std::string famn(famName);
  if(grpn.empty() || famn.empty())
    throw INTERP_KERNEL::Exception("MEDFileMesh::addGrpOnFamily : input strings must be non null !");
  std::map<std::string, std::vector<std::string> >::iterator it=_groups.find(grpn);
  if(it==_groups.end())
    {
      _groups[grpn].push_back(famn);
    }
  else
    {
      std::vector<std::string>::iterator it2=std::find((*it).second.begin(),(*it).second.end(),famn);
      if(it2==(*it).second.end())
        (*it).second.push_back(famn);
    }
}

void MEDFileMesh::setFamilyInfo(const std::map<std::string,int>& info)
{
  _families=info;
}

void MEDFileMesh::setGroupInfo(const std::map<std::string, std::vector<std::string> >&info)
{
  _groups=info;
}

int MEDFileMesh::getFamilyId(const char *name) const throw(INTERP_KERNEL::Exception)
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

std::vector<int> MEDFileMesh::getFamiliesIds(const std::vector<std::string>& fams) const throw(INTERP_KERNEL::Exception)
{
  std::vector<int> ret(fams.size());
  int i=0;
  for(std::vector<std::string>::const_iterator it=fams.begin();it!=fams.end();it++,i++)
    {
      std::map<std::string, int>::const_iterator it2=_families.find(*it);
      if(it2==_families.end())
        {
          std::vector<std::string> fams2=getFamiliesNames();
          std::ostringstream oss; oss << "No such familyname \"" << *it << "\" in input list !\nAvailable families are :";
          std::copy(fams2.begin(),fams2.end(),std::ostream_iterator<std::string>(oss," "));
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      ret[i]=(*it2).second;
    }
  return ret;
}

int MEDFileMesh::getMaxFamilyId() const throw(INTERP_KERNEL::Exception)
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

/*!
 * Returns the first (in lexical order) family name having family id equal to 'id'.
 */
std::string MEDFileMesh::getFamilyNameGivenId(int id) const throw(INTERP_KERNEL::Exception)
{
  for(std::map<std::string,int>::const_iterator it=_families.begin();it!=_families.end();it++)
    if((*it).second==id)
      return (*it).first;
  std::ostringstream oss; oss << "MEDFileUMesh::getFamilyNameGivenId : no such family id : " << id;
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

std::string MEDFileMesh::simpleRepr() const
{
  std::ostringstream oss;
  oss << "(*************************************)\n(* GENERAL INFORMATION ON THE MESH : *)\n(*************************************)\n";
  oss << "- Name of the mesh : <<" << getName() << ">>\n";
  oss << "- Description associated to the mesh : " << getDescription() << std::endl;
  return oss.str();
}

DataArrayInt *MEDFileMesh::getGroupArr(int meshDimRelToMaxExt, const char *grp, bool renum) const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> tmp(1);
  tmp[0]=grp;
  DataArrayInt *ret=getGroupsArr(meshDimRelToMaxExt,tmp,renum);
  ret->setName(grp);
  return ret;
}

DataArrayInt *MEDFileMesh::getGroupsArr(int meshDimRelToMaxExt, const std::vector<std::string>& grps, bool renum) const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> fams2=getFamiliesOnGroups(grps);
  return getFamiliesArr(meshDimRelToMaxExt,fams2,renum);
}

DataArrayInt *MEDFileMesh::getFamilyArr(int meshDimRelToMaxExt, const char *fam, bool renum) const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> tmp(1);
  tmp[0]=fam;
  DataArrayInt *ret=getFamiliesArr(meshDimRelToMaxExt,tmp,renum);
  ret->setName(fam);
  return ret;
}

DataArrayInt *MEDFileMesh::getNodeGroupArr(const char *grp, bool renum) const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> tmp(1);
  tmp[0]=grp;
  DataArrayInt *ret=getNodeGroupsArr(tmp,renum);
  ret->setName(grp);
  return ret;
}

DataArrayInt *MEDFileMesh::getNodeGroupsArr(const std::vector<std::string>& grps, bool renum) const throw(INTERP_KERNEL::Exception)
{
  return getGroupsArr(1,grps,renum);
}

DataArrayInt *MEDFileMesh::getNodeFamilyArr(const char *fam, bool renum) const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> tmp(1);
  tmp[0]=fam;
  DataArrayInt *ret=getNodeFamiliesArr(tmp,renum);
  ret->setName(fam);
  return ret;
}

DataArrayInt *MEDFileMesh::getNodeFamiliesArr(const std::vector<std::string>& fams, bool renum) const throw(INTERP_KERNEL::Exception)
{
  return getFamiliesArr(1,fams,renum);
}

void MEDFileMesh::setGroupsAtLevel(int meshDimRelToMaxExt, const std::vector<const DataArrayInt *>& grps, bool renum) throw(INTERP_KERNEL::Exception)
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
  setFamilyFieldArr(meshDimRelToMaxExt,fam);
}

/*!
 * This method append into '_families' attribute the families whose ids are in 'famIds'. Warning 'famIds' are expected to be ids
 * not in '_families'. Groups information are given in parameters in order to give to families representative names.
 * For the moment, the two last input parameters are not taken into account.
 */
void MEDFileMesh::appendFamilyEntries(const std::set<int>& famIds, const std::vector< std::vector<int> >& fidsOfGrps, const std::vector<std::string>& grpNames)
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

void MEDFileMesh::TranslateFamilyIds(int offset, DataArrayInt *famArr, std::vector< std::vector<int> >& famIdsPerGrp)
{
  famArr->applyLin(1,offset,0);
  for(std::vector< std::vector<int> >::iterator it1=famIdsPerGrp.begin();it1!=famIdsPerGrp.end();it1++)
    std::transform((*it1).begin(),(*it1).end(),(*it1).begin(),std::bind2nd(std::plus<int>(),offset));
}

/*!
 * This method should be called by any set* method of subclasses to deal automatically with _name attribute.
 * If _name attribute is empty the name of 'm' if taken as _name attribute.
 * If _name is not empty and that 'm' has the same name nothing is done.
 * If _name is not emplt and that 'm' has \b NOT the same name an exception is thrown.
 */
void MEDFileMesh::dealWithTinyInfo(const MEDCouplingMesh *m) throw(INTERP_KERNEL::Exception)
{
  if(_name.empty())
    _name=m->getName();
  else
    {
      std::string name(m->getName());
      if(!name.empty())
        {
          if(_name!=name)
            {
              std::ostringstream oss; oss << "MEDFileMesh::dealWithTinyInfo : name of current MEDfile mesh is '" << _name << "' whereas name of input mesh is : '";
              oss << name << "' ! Names must match !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
    }
}

void MEDFileMesh::getFamilyRepr(std::ostream& oss) const
{
  oss << "(**************************)\n(* FAMILIES OF THE MESH : *)\n(**************************)\n";
  for(std::map<std::string,int>::const_iterator it=_families.begin();it!=_families.end();it++)
    {
      oss << "- Family with name \"" << (*it).first << "\" with number " << (*it).second << std::endl;
      oss << "  - Groups lying on this family : ";
      std::vector<std::string> grps=getGroupsOnFamily((*it).first.c_str());
      std::copy(grps.begin(),grps.end(),std::ostream_iterator<std::string>(oss," "));
      oss << std::endl << std::endl;
    }
}

MEDFileUMesh *MEDFileUMesh::New(const char *fileName, const char *mName, int dt, int it) throw(INTERP_KERNEL::Exception)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  return new MEDFileUMesh(fid,mName,dt,it);
}

MEDFileUMesh *MEDFileUMesh::New(const char *fileName) throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> ms=MEDLoader::GetMeshNames(fileName);
  if(ms.empty())
    {
      std::ostringstream oss; oss << "MEDFileUMesh::New : no meshes in file \"" << fileName << "\" !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  int dt,it;
  ParaMEDMEM::MEDCouplingMeshType meshType;
  std::string dummy2;
  MEDFileMeshL2::GetMeshIdFromName(fid,ms.front().c_str(),meshType,dt,it,dummy2);
  return new MEDFileUMesh(fid,ms.front().c_str(),dt,it);
}

MEDFileUMesh *MEDFileUMesh::New()
{
  return new MEDFileUMesh;
}

bool MEDFileUMesh::isEqual(const MEDFileMesh *other, double eps, std::string& what) const
{
  if(!MEDFileMesh::isEqual(other,eps,what))
    return false;
  const MEDFileUMesh *otherC=dynamic_cast<const MEDFileUMesh *>(other);
  if(!otherC)
    {
      what="Mesh types differ ! This is unstructured and other is NOT !";
      return false;
    }
  clearNonDiscrAttributes();
  otherC->clearNonDiscrAttributes();
  const DataArrayDouble *coo1=_coords;
  const DataArrayDouble *coo2=otherC->_coords;
  if((coo1==0 && coo2!=0) || (coo1!=0 && coo2==0))
    {
      what="Mismatch of coordinates ! One is defined and not other !";
      return false;
    }
  if(coo1)
    {
      bool ret=coo1->isEqual(*coo2,eps);
      if(!ret)
        {
          what="Coords differ !";
          return false;
        }
    }
  const DataArrayInt *famc1=_fam_coords;
  const DataArrayInt *famc2=otherC->_fam_coords;
  if((famc1==0 && famc2!=0) || (famc1!=0 && famc2==0))
    {
      what="Mismatch of families arr on nodes ! One is defined and not other !";
      return false;
    }
  if(famc1)
    {
      bool ret=famc1->isEqual(*famc2);
      if(!ret)
        {
          what="Families arr on node differ !";
          return false;
        }
    }
  const DataArrayInt *numc1=_num_coords;
  const DataArrayInt *numc2=otherC->_num_coords;
  if((numc1==0 && numc2!=0) || (numc1!=0 && numc2==0))
    {
      what="Mismatch of numbering arr on nodes ! One is defined and not other !";
      return false;
    }
  if(numc1)
    {
      bool ret=numc1->isEqual(*numc2);
      if(!ret)
        {
          what="Numbering arr on node differ !";
          return false;
        }
    }
  if(_ms.size()!=otherC->_ms.size())
    {
      what="Number of levels differs !";
      return false;
    }
  std::size_t sz=_ms.size();
  for(std::size_t i=0;i<sz;i++)
    {
      const MEDFileUMeshSplitL1 *s1=_ms[i];
      const MEDFileUMeshSplitL1 *s2=otherC->_ms[i];
      if((s1==0 && s2!=0) || (s1!=0 && s2==0))
        {
          what="Mismatch of presence of sub levels !";
          return false;
        }
      if(s1)
        {
          bool ret=s1->isEqual(s2,eps,what);
          if(!ret)
            return false;
        }
    }
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> > _ms;
  return true;
}

void MEDFileUMesh::clearNonDiscrAttributes() const
{
  MEDFileMesh::clearNonDiscrAttributes();
  const DataArrayDouble *coo1=_coords;
  if(coo1)
    ((DataArrayDouble *)coo1)->setName("");//This parameter is not discriminant for comparison
  const DataArrayInt *famc1=_fam_coords;
  if(famc1)
    ((DataArrayInt *)famc1)->setName("");//This parameter is not discriminant for comparison
  const DataArrayInt *numc1=_num_coords;
  if(numc1)
    ((DataArrayInt *)numc1)->setName("");//This parameter is not discriminant for comparison
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::const_iterator it=_ms.begin();it!=_ms.end();it++)
    {
      const MEDFileUMeshSplitL1 *tmp=(*it);
      if(tmp)
        tmp->clearNonDiscrAttributes();
    }
}

MEDFileUMesh::MEDFileUMesh()
{
}

MEDFileUMesh::MEDFileUMesh(med_idt fid, const char *mName, int dt, int it) throw(INTERP_KERNEL::Exception)
try
  {
    loadUMeshFromFile(fid,mName,dt,it);
  }
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

void MEDFileUMesh::loadUMeshFromFile(med_idt fid, const char *mName, int dt, int it) throw(INTERP_KERNEL::Exception)
{
  MEDFileUMeshL2 loaderl2;
  ParaMEDMEM::MEDCouplingMeshType meshType;
  int dummy0,dummy1;
  std::string dummy2;
  int mid=MEDFileUMeshL2::GetMeshIdFromName(fid,mName,meshType,dummy0,dummy1,dummy2);
  if(meshType!=UNSTRUCTURED)
    {
         std::ostringstream oss; oss << "Trying to load as unstructured an existing mesh with name '" << mName << "' !";
         throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  loaderl2.loadAll(fid,mid,mName,dt,it);
  int lev=loaderl2.getNumberOfLevels();
  _ms.resize(lev);
  for(int i=0;i<lev;i++)
    {
      if(!loaderl2.emptyLev(i))
        _ms[i]=new MEDFileUMeshSplitL1(loaderl2,mName,i);
      else
        _ms[i]=0;
    }
  MEDFileMeshL2::ReadFamiliesAndGrps(fid,mName,_families,_groups);
  //
  setName(loaderl2.getName());
  setDescription(loaderl2.getDescription());
  setIteration(loaderl2.getIteration());
  setOrder(loaderl2.getOrder());
  setTimeValue(loaderl2.getTime());
  setTimeUnit(loaderl2.getTimeUnit());
  _coords=loaderl2.getCoords();
  _fam_coords=loaderl2.getCoordsFamily();
  _num_coords=loaderl2.getCoordsNum();
  computeRevNum();
}

MEDFileUMesh::~MEDFileUMesh()
{
}

void MEDFileUMesh::write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception)
{
  if(_name.empty())
    throw INTERP_KERNEL::Exception("MEDFileUMesh : name is empty. MED file ask for a NON EMPTY name !");
  if(!existsFamily(0))
    (const_cast<MEDFileUMesh *>(this))->addFamily(DFT_FAM_NAME,0);
  med_access_mode medmod=MEDFileUtilities::TraduceWriteMode(mode);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,medmod);
  std::ostringstream oss; oss << "MEDFileUMesh : error on attempt to write in file : \"" << fileName << "\""; 
  MEDFileUtilities::CheckMEDCode(fid,fid,oss.str().c_str());
  const DataArrayDouble *coo=_coords;
  INTERP_KERNEL::AutoPtr<char> maa=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> desc=MEDLoaderBase::buildEmptyString(MED_COMMENT_SIZE);
  MEDLoaderBase::safeStrCpy(_name.c_str(),MED_NAME_SIZE,maa,_too_long_str);
  MEDLoaderBase::safeStrCpy(_desc_name.c_str(),MED_COMMENT_SIZE,desc,_too_long_str);
  int spaceDim=coo?coo->getNumberOfComponents():0;
  int mdim=getMeshDimension();
  INTERP_KERNEL::AutoPtr<char> comp=MEDLoaderBase::buildEmptyString(spaceDim*MED_SNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> unit=MEDLoaderBase::buildEmptyString(spaceDim*MED_SNAME_SIZE);
  for(int i=0;i<spaceDim;i++)
    {
      std::string info=coo->getInfoOnComponent(i);
      std::string c,u;
      MEDLoaderBase::splitIntoNameAndUnit(info,c,u);
      MEDLoaderBase::safeStrCpy2(c.c_str(),MED_SNAME_SIZE-1,comp+i*MED_SNAME_SIZE,_too_long_str);//MED_TAILLE_PNOM-1 to avoid to write '\0' on next compo
      MEDLoaderBase::safeStrCpy2(u.c_str(),MED_SNAME_SIZE-1,unit+i*MED_SNAME_SIZE,_too_long_str);//MED_TAILLE_PNOM-1 to avoid to write '\0' on next compo
    }
  MEDmeshCr(fid,maa,spaceDim,mdim,MED_UNSTRUCTURED_MESH,desc,"",MED_SORT_DTIT,MED_CARTESIAN,comp,unit);
  MEDFileUMeshL2::WriteCoords(fid,maa,_iteration,_order,_time,_coords,_fam_coords,_num_coords);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::const_iterator it=_ms.begin();it!=_ms.end();it++)
    if((const MEDFileUMeshSplitL1 *)(*it)!=0)
      (*it)->write(fid,maa,mdim);
  MEDFileUMeshL2::WriteFamiliesAndGrps(fid,maa,_families,_groups,_too_long_str);
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

/*!
 * This methods returns all relative mesh levels where group 'grp' is defined \b excluded \b nodes.
 * To include nodes call MEDFileUMesh::getGrpNonEmptyLevelsExt method.
 */
std::vector<int> MEDFileUMesh::getGrpNonEmptyLevels(const char *grp) const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> fams=getFamiliesOnGroup(grp);
  return getFamsNonEmptyLevels(fams);
}

/*!
 * This method is a generalization of MEDFileUMesh::getGrpNonEmptyLevelsExt. It looks at the node level to state if the group 'grp' has a part lying on node.
 */
std::vector<int> MEDFileUMesh::getGrpNonEmptyLevelsExt(const char *grp) const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> fams=getFamiliesOnGroup(grp);
  return getFamsNonEmptyLevelsExt(fams);
}

/*!
 * This methods returns all relative mesh levels where family 'fam' is defined \b excluded \b nodes.
 * To include nodes call MEDFileUMesh::getFamNonEmptyLevelsExt method.
 */
std::vector<int> MEDFileUMesh::getFamNonEmptyLevels(const char *fam) const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> fams(1,std::string(fam));
  return getFamsNonEmptyLevels(fams);
}

/*!
 * This method is a generalization of MEDFileUMesh::getFamNonEmptyLevels. It looks at the node level to state if the family 'fam' has a part lying on node.
 */
std::vector<int> MEDFileUMesh::getFamNonEmptyLevelsExt(const char *fam) const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> fams(1,std::string(fam));
  return getFamsNonEmptyLevelsExt(fams);
}

/*!
 * This methods returns all relative mesh levels where groups 'grps' are defined \b excluded \b nodes.
 * To include nodes call MEDFileUMesh::getGrpsNonEmptyLevelsExt method.
 */
std::vector<int> MEDFileUMesh::getGrpsNonEmptyLevels(const std::vector<std::string>& grps) const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> fams=getFamiliesOnGroups(grps);
  return getFamsNonEmptyLevels(fams);
}

/*!
 * This method is a generalization of MEDFileUMesh::getGrpsNonEmptyLevels. It looks at the node level to state if the families 'fams' has a part lying on node.
 */
std::vector<int> MEDFileUMesh::getGrpsNonEmptyLevelsExt(const std::vector<std::string>& grps) const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> fams=getFamiliesOnGroups(grps);
  return getFamsNonEmptyLevelsExt(fams);
}

/*!
 * This methods returns all relative mesh levels where families 'fams' are defined \b excluded \b nodes.
 * To include nodes call MEDFileUMesh::getFamsNonEmptyLevelsExt method.
 */
std::vector<int> MEDFileUMesh::getFamsNonEmptyLevels(const std::vector<std::string>& fams) const throw(INTERP_KERNEL::Exception)
{
  std::vector<int> ret;
  std::vector<int> levs=getNonEmptyLevels();
  std::vector<int> famIds=getFamiliesIds(fams);
  for(std::vector<int>::const_iterator it=levs.begin();it!=levs.end();it++)
    if(_ms[-(*it)]->presenceOfOneFams(famIds))
      ret.push_back(*it);
  return ret;
}

/*!
 * This method is a generalization of MEDFileUMesh::getFamsNonEmptyLevels. It looks at the node level to state if the families 'fams' has a part lying on node.
 */
std::vector<int> MEDFileUMesh::getFamsNonEmptyLevelsExt(const std::vector<std::string>& fams) const throw(INTERP_KERNEL::Exception)
{
  std::vector<int> ret0=getFamsNonEmptyLevels(fams);
  const DataArrayInt *famCoords=_fam_coords;
  if(!famCoords)
    return ret0;
  std::vector<int> famIds=getFamiliesIds(fams);
  if(famCoords->presenceOfValue(famIds))
    {
      std::vector<int> ret(ret0.size()+1);
      ret[0]=1;
      std::copy(ret0.begin(),ret0.end(),ret.begin()+1);
      return ret;
    }
  else
    return ret0;
}

int MEDFileUMesh::getMeshDimension() const throw(INTERP_KERNEL::Exception)
{
  int lev=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::const_iterator it=_ms.begin();it!=_ms.end();it++,lev++)
    if((const MEDFileUMeshSplitL1 *)(*it)!=0)
      return (*it)->getMeshDimension()+lev;
  throw INTERP_KERNEL::Exception("MEDFileUMesh::getMeshDimension : impossible to find a mesh dimension !");
}

std::string MEDFileUMesh::simpleRepr() const
{
  std::ostringstream oss;
  oss << MEDFileMesh::simpleRepr();
  const DataArrayDouble *coo=_coords;
  oss << "- The dimension of the space is ";
  static const char MSG1[]= "*** NO COORDS SET ***";
  static const char MSG2[]= "*** NO CONNECTIVITY SET FOR THIS LEVEL***";
  if(coo)
    oss << _coords->getNumberOfComponents() << std::endl;
  else
    oss << MSG1 << std::endl;
  oss << "- Type of the mesh : UNSTRUCTURED\n";
  oss << "- Number of nodes : ";
  if(coo)
    oss << _coords->getNumberOfTuples() << std::endl;
  else
    oss << MSG1 << std::endl;
  std::size_t nbOfLev=_ms.size();
  oss << "- Number of levels allocated : " << nbOfLev << std::endl;
  for(std::size_t i=0;i<nbOfLev;i++)
    {
      const MEDFileUMeshSplitL1 *lev=_ms[i];
      oss << "  - Level #" << -((int) i) << " has dimension : ";
      if(lev)
        {
          oss << lev->getMeshDimension() << std::endl;
          lev->simpleRepr(oss);
        }
      else
        oss << MSG2 << std::endl;
    }
  oss << "- Number of families : " << _families.size() << std::endl << std::endl;
  if(coo)
    {
      oss << "(***********************)\n(* NODES OF THE MESH : *)\n(***********************)\n";
      oss << "- Names of coordinates :" << std::endl;
      std::vector<std::string> vars=coo->getVarsOnComponent();
      std::copy(vars.begin(),vars.end(),std::ostream_iterator<std::string>(oss," "));
      oss << std::endl << "- Units of coordinates : " << std::endl;
      std::vector<std::string> units=coo->getUnitsOnComponent();
      std::copy(units.begin(),units.end(),std::ostream_iterator<std::string>(oss," "));
    }
  oss << std::endl << std::endl;
  getFamilyRepr(oss);
  return oss.str();
}

std::string MEDFileUMesh::advancedRepr() const
{
  return simpleRepr();
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
    return _num_coords;
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
  synchronizeTinyInfoOnLeaves();
  std::vector<std::string> tmp(1);
  tmp[0]=grp;
  MEDCouplingUMesh *ret=getGroups(meshDimRelToMaxExt,tmp,renum);
  ret->setName(grp);
  return ret;
}

MEDCouplingUMesh *MEDFileUMesh::getGroups(int meshDimRelToMaxExt, const std::vector<std::string>& grps, bool renum) const throw(INTERP_KERNEL::Exception)
{
  synchronizeTinyInfoOnLeaves();
  std::vector<std::string> fams2=getFamiliesOnGroups(grps);
  return getFamilies(meshDimRelToMaxExt,fams2,renum);
}

MEDCouplingUMesh *MEDFileUMesh::getFamily(int meshDimRelToMaxExt, const char *fam, bool renum) const throw(INTERP_KERNEL::Exception)
{
  synchronizeTinyInfoOnLeaves();
  std::vector<std::string> tmp(1);
  tmp[0]=fam;
  MEDCouplingUMesh *ret=getFamilies(meshDimRelToMaxExt,tmp,renum);
  ret->setName(fam);
  return ret;
}

MEDCouplingUMesh *MEDFileUMesh::getFamilies(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum) const throw(INTERP_KERNEL::Exception)
{
  synchronizeTinyInfoOnLeaves();
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
      if((const DataArrayInt *)_fam_coords)
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da=_fam_coords->getIdsEqualList(famIds);
          if(renum)
            return MEDFileUMeshSplitL1::Renumber(_num_coords,da);
          else
            {
              da->incrRef();
              return da;
            }
        }
      else
        throw INTERP_KERNEL::Exception("MEDFileUMesh::getFamiliesArr : no family array specified on nodes !");
    }
  const MEDFileUMeshSplitL1 *l1=getMeshAtLevSafe(meshDimRelToMaxExt);
  return l1->getFamilyPartArr(famIds,renum);
}

MEDCouplingUMesh *MEDFileUMesh::getMeshAtLevel(int meshDimRelToMaxExt, bool renum) const throw(INTERP_KERNEL::Exception)
{
  synchronizeTinyInfoOnLeaves();
  if(meshDimRelToMaxExt==1)
    {
      if(!renum)
        {
          MEDCouplingUMesh *umesh=MEDCouplingUMesh::New();
          MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> cc=_coords->deepCpy();
          umesh->setCoords(cc);
          MEDFileUMeshSplitL1::ClearNonDiscrAttributes(umesh);
          return umesh;
        }
    }
  const MEDFileUMeshSplitL1 *l1=getMeshAtLevSafe(meshDimRelToMaxExt);
  return l1->getWholeMesh(renum);
}

MEDCouplingMesh *MEDFileUMesh::getGenMeshAtLevel(int meshDimRelToMax, bool renum) const throw(INTERP_KERNEL::Exception)
{
  return getMeshAtLevel(meshDimRelToMax,renum);
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

void MEDFileUMesh::setFamilyNameAttachedOnId(int id, const std::string& newFamName) throw(INTERP_KERNEL::Exception)
{
  std::string oldName=getFamilyNameGivenId(id);
  _families.erase(oldName);
  _families[newFamName]=id;
}

void MEDFileUMesh::removeMeshAtLevel(int meshDimRelToMax) throw(INTERP_KERNEL::Exception)
{
  std::vector<int> levSet=getNonEmptyLevels();
  std::vector<int>::const_iterator it=std::find(levSet.begin(),levSet.end(),meshDimRelToMax);
  if(it==levSet.end())
    throw INTERP_KERNEL::Exception("MEDFileUMesh::removeMeshAtLevel : the requested level is not existing !");
  int pos=(-meshDimRelToMax);
  _ms[pos]=0;
}

void MEDFileUMesh::setMeshAtLevel(int meshDimRelToMax, MEDCouplingUMesh *m, bool newOrOld) throw(INTERP_KERNEL::Exception)
{
  setMeshAtLevelGen(meshDimRelToMax,m,newOrOld);
}

void MEDFileUMesh::setMeshAtLevelGen(int meshDimRelToMax, MEDCouplingUMesh *m, bool newOrOld) throw(INTERP_KERNEL::Exception)
{
  dealWithTinyInfo(m);
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
      if(m->getCoords()!=_coords)
        throw INTERP_KERNEL::Exception("MEDFileUMesh::setMeshAtLevel : Invalid Given Mesh ! The coordinates are not the same ! try to use tryToShareSameCoords !");
      int sz=(-meshDimRelToMax)+1;
      if(sz>=(int)_ms.size())
        _ms.resize(sz);
      checkMeshDimCoherency(m->getMeshDimension(),meshDimRelToMax);
      _ms[sz-1]=new MEDFileUMeshSplitL1(m,newOrOld);
    }
  else
    _ms[-meshDimRelToMax]=new MEDFileUMeshSplitL1(m,newOrOld);
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

void MEDFileUMesh::setFamilyFieldArr(int meshDimRelToMaxExt, DataArrayInt *famArr) throw(INTERP_KERNEL::Exception)
{
  if(meshDimRelToMaxExt==1)
    {
      famArr->incrRef();
      _fam_coords=famArr;
      return ;
    }
  if(meshDimRelToMaxExt>1)
    throw INTERP_KERNEL::Exception("MEDFileUMesh::setFamilyFieldArr : Dimension request is invalid (>1) !");
  int traducedRk=-meshDimRelToMaxExt;
  if(traducedRk>=(int)_ms.size())
    throw INTERP_KERNEL::Exception("Invalid mesh dim relative to max given ! To low !");
  if((MEDFileUMeshSplitL1 *)_ms[traducedRk]==0)
    throw INTERP_KERNEL::Exception("On specified lev (or entity) no cells exists !");
  return _ms[traducedRk]->setFamilyArr(famArr);
}

void MEDFileUMesh::setRenumFieldArr(int meshDimRelToMaxExt, DataArrayInt *renumArr) throw(INTERP_KERNEL::Exception)
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

void MEDFileUMesh::synchronizeTinyInfoOnLeaves() const
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::const_iterator it=_ms.begin();it!=_ms.end();it++)
    if((const MEDFileUMeshSplitL1 *)(*it))
      (*it)->synchronizeTinyInfo(*this);
}

void MEDFileUMesh::changeFamilyIdArr(int oldId, int newId) throw(INTERP_KERNEL::Exception)
{
  DataArrayInt *arr=_fam_coords;
  if(arr)
    arr->changeValue(oldId,newId);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::iterator it=_ms.begin();it!=_ms.end();it++)
    {
      MEDFileUMeshSplitL1 *sp=(*it);
      if(sp)
        {
          sp->changeFamilyIdArr(oldId,newId);
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

MEDFileCMesh *MEDFileCMesh::New()
{
  return new MEDFileCMesh;
}

MEDFileCMesh *MEDFileCMesh::New(const char *fileName) throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> ms=MEDLoader::GetMeshNames(fileName);
  if(ms.empty())
    {
      std::ostringstream oss; oss << "MEDFileUMesh::New : no meshes in file \"" << fileName << "\" !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  int dt,it;
  ParaMEDMEM::MEDCouplingMeshType meshType;
  std::string dummy2;
  MEDFileMeshL2::GetMeshIdFromName(fid,ms.front().c_str(),meshType,dt,it,dummy2);
  return new MEDFileCMesh(fid,ms.front().c_str(),dt,it);
}

MEDFileCMesh *MEDFileCMesh::New(const char *fileName, const char *mName, int dt, int it) throw(INTERP_KERNEL::Exception)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  return new MEDFileCMesh(fid,mName,dt,it);
}

int MEDFileCMesh::getMeshDimension() const throw(INTERP_KERNEL::Exception)
{
  if(!((const MEDCouplingCMesh*)_cmesh))
    throw INTERP_KERNEL::Exception("MEDFileCMesh::getMeshDimension : unable to get meshdimension because no mesh set !");
  return _cmesh->getMeshDimension();
}

std::string MEDFileCMesh::simpleRepr() const
{
  return MEDFileMesh::simpleRepr();
}

std::string MEDFileCMesh::advancedRepr() const
{
  return simpleRepr();
}

bool MEDFileCMesh::isEqual(const MEDFileMesh *other, double eps, std::string& what) const
{
  if(!MEDFileMesh::isEqual(other,eps,what))
    return false;
  const MEDFileCMesh *otherC=dynamic_cast<const MEDFileCMesh *>(other);
  if(!otherC)
    {
      what="Mesh types differ ! This is cartesian and other is NOT !";
      return false;
    }
  clearNonDiscrAttributes();
  otherC->clearNonDiscrAttributes();
  const MEDCouplingCMesh *coo1=_cmesh;
  const MEDCouplingCMesh *coo2=otherC->_cmesh;
  if((coo1==0 && coo2!=0) || (coo1!=0 && coo2==0))
    {
      what="Mismatch of cartesian meshes ! One is defined and not other !";
      return false;
    }
  if(coo1)
    {
      bool ret=coo1->isEqual(coo2,eps);
      if(!ret)
        {
          what="cartesian meshes differ !";
          return false;
        }
    }
  const DataArrayInt *famc1=_fam_nodes;
  const DataArrayInt *famc2=otherC->_fam_nodes;
  if((famc1==0 && famc2!=0) || (famc1!=0 && famc2==0))
    {
      what="Mismatch of families arr on nodes ! One is defined and not other !";
      return false;
    }
  if(famc1)
    {
      bool ret=famc1->isEqual(*famc2);
      if(!ret)
        {
          what="Families arr on nodes differ !";
          return false;
        }
    }
  famc1=_fam_cells;
  famc2=otherC->_fam_cells;
  if((famc1==0 && famc2!=0) || (famc1!=0 && famc2==0))
    {
      what="Mismatch of families arr on cells ! One is defined and not other !";
      return false;
    }
  if(famc1)
    {
      bool ret=famc1->isEqual(*famc2);
      if(!ret)
        {
          what="Families arr on cells differ !";
          return false;
        }
    }
  famc1=_num_nodes;
  famc2=otherC->_num_nodes;
  if((famc1==0 && famc2!=0) || (famc1!=0 && famc2==0))
    {
      what="Mismatch of numbering arr on nodes ! One is defined and not other !";
      return false;
    }
  if(famc1)
    {
      bool ret=famc1->isEqual(*famc2);
      if(!ret)
        {
          what="Numbering arr on nodes differ !";
          return false;
        }
    }
  famc1=_num_cells;
  famc2=otherC->_num_cells;
  if((famc1==0 && famc2!=0) || (famc1!=0 && famc2==0))
    {
      what="Mismatch of numbering arr on cells ! One is defined and not other !";
      return false;
    }
  if(famc1)
    {
      bool ret=famc1->isEqual(*famc2);
      if(!ret)
        {
          what="Numbering arr on cells differ !";
          return false;
        }
    }
  return true;
}

void MEDFileCMesh::clearNonDiscrAttributes() const
{
  MEDFileMesh::clearNonDiscrAttributes();
  MEDFileUMeshSplitL1::ClearNonDiscrAttributes(_cmesh);
  const DataArrayInt *tmp=_fam_nodes;
  if(tmp)
    (const_cast<DataArrayInt *>(tmp))->setName("");
  tmp=_num_nodes;
  if(tmp)
    (const_cast<DataArrayInt *>(tmp))->setName("");
  tmp=_fam_cells;
  if(tmp)
    (const_cast<DataArrayInt *>(tmp))->setName("");
  tmp=_num_cells;
  if(tmp)
    (const_cast<DataArrayInt *>(tmp))->setName("");
}

MEDFileCMesh::MEDFileCMesh()
{
}

MEDFileCMesh::MEDFileCMesh(med_idt fid, const char *mName, int dt, int it) throw(INTERP_KERNEL::Exception)
try
  {
    loadCMeshFromFile(fid,mName,dt,it);
  }
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

void MEDFileCMesh::loadCMeshFromFile(med_idt fid, const char *mName, int dt, int it) throw(INTERP_KERNEL::Exception)
{
  MEDFileCMeshL2 loaderl2;
  ParaMEDMEM::MEDCouplingMeshType meshType;
  int dummy0,dummy1;
  std::string dtunit;
  int mid=MEDFileMeshL2::GetMeshIdFromName(fid,mName,meshType,dummy0,dummy1,dtunit);
  if(meshType!=CARTESIAN)
    {
      std::ostringstream oss; oss << "Trying to load as cartesian an existing mesh with name '" << mName << "' that is NOT cartesian !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  loaderl2.loadAll(fid,mid,mName,dt,it);
  MEDCouplingCMesh *mesh=loaderl2.getMesh();
  mesh->incrRef();
  _cmesh=mesh;
  setName(loaderl2.getName());
  setDescription(loaderl2.getDescription());
  setIteration(loaderl2.getIteration());
  setOrder(loaderl2.getOrder());
  setTimeValue(loaderl2.getTime());
  setTimeUnit(loaderl2.getTimeUnit());
  MEDFileMeshL2::ReadFamiliesAndGrps(fid,mName,_families,_groups);
  med_bool chgt=MED_FALSE,trsf=MED_FALSE;
  int nbOfElt=MEDmeshnEntity(fid,mName,dt,it,MED_NODE,MED_NONE,MED_FAMILY_NUMBER,MED_NODAL,&chgt,&trsf);
  if(nbOfElt>0)
    {
      _fam_nodes=DataArrayInt::New();
      _fam_nodes->alloc(nbOfElt,1);
      MEDmeshEntityFamilyNumberRd(fid,mName,dt,it,MED_NODE,MED_NONE,_fam_nodes->getPointer());
    }
  nbOfElt=MEDmeshnEntity(fid,mName,dt,it,MED_NODE,MED_NONE,MED_NUMBER,MED_NODAL,&chgt,&trsf);
  if(nbOfElt>0)
    {
      _num_nodes=DataArrayInt::New();
      _num_nodes->alloc(nbOfElt,1);
      MEDmeshEntityNumberRd(fid,mName,dt,it,MED_NODE,MED_NONE,_num_nodes->getPointer());
    }
  int spaceDim=mesh->getSpaceDimension();
  med_geometry_type geoTypeReq=MED_NONE;
  switch(spaceDim)
    {
    case 3:
      geoTypeReq=MED_HEXA8;
      break;
    case 2:
      geoTypeReq=MED_QUAD4;
      break;
    case 1:
      geoTypeReq=MED_SEG2;
      break;
    case 0:
      geoTypeReq=MED_POINT1;
      break;
    default:
      throw INTERP_KERNEL::Exception("Invalid spacedim detected for cartesian mesh ! Must be in (1,2,3) !");
    }
  nbOfElt=MEDmeshnEntity(fid,mName,dt,it,MED_CELL,geoTypeReq,MED_FAMILY_NUMBER,MED_NODAL,&chgt,&trsf);
  if(nbOfElt>0)
    {
      _fam_cells=DataArrayInt::New();
      _fam_cells->alloc(nbOfElt,1);
      MEDmeshEntityFamilyNumberRd(fid,mName,dt,it,MED_CELL,geoTypeReq,_fam_cells->getPointer());
    }
  nbOfElt=MEDmeshnEntity(fid,mName,dt,it,MED_CELL,geoTypeReq,MED_NUMBER,MED_NODAL,&chgt,&trsf);
  if(nbOfElt>0)
    {
      _num_cells=DataArrayInt::New();
      _num_cells->alloc(nbOfElt,1);
      MEDmeshEntityNumberRd(fid,mName,dt,it,MED_CELL,geoTypeReq,_num_cells->getPointer());
    }
}

void MEDFileCMesh::changeFamilyIdArr(int oldId, int newId) throw(INTERP_KERNEL::Exception)
{
  DataArrayInt *arr=_fam_nodes;
  if(arr)
    arr->changeValue(oldId,newId);
  arr=_fam_cells;
  if(arr)
    arr->changeValue(oldId,newId);
}

const MEDCouplingCMesh *MEDFileCMesh::getMesh() const
{
  synchronizeTinyInfoOnLeaves();
  return _cmesh;
}

MEDCouplingMesh *MEDFileCMesh::getGenMeshAtLevel(int meshDimRelToMax, bool renum) const throw(INTERP_KERNEL::Exception)
{
  if(renum)
    throw INTERP_KERNEL::Exception("MEDFileCMesh does not support renumbering ! To do it perform request of renum array directly !");
  if(meshDimRelToMax!=0)
    throw INTERP_KERNEL::Exception("MEDFileCMesh does not support multi level for mesh 0 expected as input !");
  const MEDCouplingCMesh *m=getMesh();
  if(m)
    m->incrRef();
  return const_cast<MEDCouplingCMesh *>(m);
}

void MEDFileCMesh::setMesh(MEDCouplingCMesh *m) throw(INTERP_KERNEL::Exception)
{
  dealWithTinyInfo(m);
  if(m)
    m->incrRef();
  _cmesh=m;
}

void MEDFileCMesh::write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception)
{
  if(_name.empty())
    throw INTERP_KERNEL::Exception("MEDFileCMesh : name is empty. MED file ask for a NON EMPTY name !");
  if(!existsFamily(0))
    (const_cast<MEDFileCMesh *>(this))->addFamily(DFT_FAM_NAME,0);
  med_access_mode medmod=MEDFileUtilities::TraduceWriteMode(mode);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,medmod);
  std::ostringstream oss; oss << "MEDFileCMesh : error on attempt to write in file : \"" << fileName << "\""; 
  MEDFileUtilities::CheckMEDCode(fid,fid,oss.str().c_str());
  INTERP_KERNEL::AutoPtr<char> maa=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> desc=MEDLoaderBase::buildEmptyString(MED_COMMENT_SIZE);
  INTERP_KERNEL::AutoPtr<char> dtunit=MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE);
  MEDLoaderBase::safeStrCpy(_name.c_str(),MED_NAME_SIZE,maa,_too_long_str);
  MEDLoaderBase::safeStrCpy(_desc_name.c_str(),MED_COMMENT_SIZE,desc,_too_long_str);
  MEDLoaderBase::safeStrCpy(_dt_unit.c_str(),MED_LNAME_SIZE,dtunit,_too_long_str);
  int spaceDim=_cmesh->getSpaceDimension();
  INTERP_KERNEL::AutoPtr<char> comp=MEDLoaderBase::buildEmptyString(spaceDim*MED_SNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> unit=MEDLoaderBase::buildEmptyString(spaceDim*MED_SNAME_SIZE);
  for(int i=0;i<spaceDim;i++)
    {
      std::string info(_cmesh->getCoordsAt(i)->getInfoOnComponent(0));
      std::string c,u;
      MEDLoaderBase::splitIntoNameAndUnit(info,c,u);
      MEDLoaderBase::safeStrCpy2(c.c_str(),MED_SNAME_SIZE-1,comp+i*MED_SNAME_SIZE,_too_long_str);//MED_TAILLE_PNOM-1 to avoid to write '\0' on next compo
      MEDLoaderBase::safeStrCpy2(u.c_str(),MED_SNAME_SIZE-1,unit+i*MED_SNAME_SIZE,_too_long_str);//MED_TAILLE_PNOM-1 to avoid to write '\0' on next compo
    }
  MEDmeshCr(fid,maa,spaceDim,spaceDim,MED_STRUCTURED_MESH,desc,dtunit,MED_SORT_DTIT,MED_CARTESIAN,comp,unit);
  MEDmeshGridTypeWr(fid,maa,MED_CARTESIAN_GRID);
  for(int i=0;i<spaceDim;i++)
    {
      const DataArrayDouble *da=_cmesh->getCoordsAt(i);
      MEDmeshGridIndexCoordinateWr(fid,maa,_iteration,_order,_time,i+1,da->getNumberOfTuples(),da->getConstPointer());
    }
  //
  med_geometry_type geoTypeReq=MED_NONE;
  switch(spaceDim)
    {
    case 3:
      geoTypeReq=MED_HEXA8;
      break;
    case 2:
      geoTypeReq=MED_QUAD4;
      break;
    case 1:
      geoTypeReq=MED_SEG2;
      break;
    case 0:
      geoTypeReq=MED_POINT1;
      break;
    default:
      throw INTERP_KERNEL::Exception("Invalid spacedim detected for cartesian mesh ! Must be in (1,2,3) !");
    }
  //
  if((const DataArrayInt *)_fam_cells)
    MEDmeshEntityFamilyNumberWr(fid,maa,_iteration,_order,MED_CELL,geoTypeReq,_fam_cells->getNumberOfTuples(),_fam_cells->getConstPointer());
  if((const DataArrayInt *)_fam_nodes)
    MEDmeshEntityFamilyNumberWr(fid,maa,_iteration,_order,MED_NODE,MED_NONE,_fam_nodes->getNumberOfTuples(),_fam_nodes->getConstPointer());
  //
  MEDFileUMeshL2::WriteFamiliesAndGrps(fid,maa,_families,_groups,_too_long_str);
}

int MEDFileCMesh::getSizeAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception)
{
  if(meshDimRelToMaxExt!=0 && meshDimRelToMaxExt!=1)
    throw INTERP_KERNEL::Exception("MEDFileCMesh::getSizeAtLevel : Only available for levels 0 or 1 !");
  if(!((const MEDCouplingCMesh *)_cmesh))
    throw INTERP_KERNEL::Exception("MEDFileCMesh::getSizeAtLevel : No cartesian mesh set !");
  if(meshDimRelToMaxExt==0)
    return _cmesh->getNumberOfCells();
  else
    return _cmesh->getNumberOfNodes();
}

void MEDFileCMesh::synchronizeTinyInfoOnLeaves() const
{
  const MEDCouplingCMesh *cmesh=_cmesh;
  (const_cast<MEDCouplingCMesh *>(cmesh))->setName(_name.c_str());
  (const_cast<MEDCouplingCMesh *>(cmesh))->setDescription(_desc_name.c_str());
  (const_cast<MEDCouplingCMesh *>(cmesh))->setTime(_time,_iteration,_order);
  (const_cast<MEDCouplingCMesh *>(cmesh))->setTimeUnit(_dt_unit.c_str());
}

std::vector<int> MEDFileCMesh::getNonEmptyLevels() const
{
  std::vector<int> ret(1);
  return ret;
}

std::vector<int> MEDFileCMesh::getNonEmptyLevelsExt() const
{
  std::vector<int> ret(2);
  ret[0]=1;
  return ret;
}

DataArrayInt *MEDFileCMesh::getFamiliesArr(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum) const throw(INTERP_KERNEL::Exception)
{
  if(meshDimRelToMaxExt!=0 && meshDimRelToMaxExt!=1)
    throw INTERP_KERNEL::Exception("MEDFileCMesh::getFamiliesArr : Only available for levels 0 or 1 !");
  std::vector<int> famIds=getFamiliesIds(fams);
  if(meshDimRelToMaxExt==1)
    {
      if((const DataArrayInt *)_fam_nodes)
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da=_fam_nodes->getIdsEqualList(famIds);
          if(renum)
            return MEDFileUMeshSplitL1::Renumber(_num_nodes,da);
          else
            {
              da->incrRef();
              return da;
            }
        }
      else
        throw INTERP_KERNEL::Exception("MEDFileCMesh::getFamiliesArr : no family array specified on nodes !");
    }
  else
    {
      if((const DataArrayInt *)_fam_cells)
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da=_fam_cells->getIdsEqualList(famIds);
          if(renum)
            return MEDFileUMeshSplitL1::Renumber(_num_cells,da);
          else
            {
              da->incrRef();
              return da;
            }
        }
      else
        throw INTERP_KERNEL::Exception("MEDFileCMesh::getFamiliesArr : no family array specified on cells !");
    }
}

void MEDFileCMesh::setFamilyFieldArr(int meshDimRelToMaxExt, DataArrayInt *famArr) throw(INTERP_KERNEL::Exception)
{
  if(meshDimRelToMaxExt!=0 && meshDimRelToMaxExt!=1)
    throw INTERP_KERNEL::Exception("MEDFileCMesh::setRenumFieldArr : Only available for levels 0 or 1 !");
  if(famArr)
    famArr->incrRef();
  if(meshDimRelToMaxExt==0)
    _fam_cells=famArr;
  else
    _fam_nodes=famArr;
}

void MEDFileCMesh::setRenumFieldArr(int meshDimRelToMaxExt, DataArrayInt *renumArr) throw(INTERP_KERNEL::Exception)
{
  if(meshDimRelToMaxExt!=0 && meshDimRelToMaxExt!=1)
    throw INTERP_KERNEL::Exception("MEDFileCMesh::setRenumFieldArr : Only available for levels 0 or 1 !");
  if(renumArr)
    renumArr->incrRef();
  if(meshDimRelToMaxExt==0)
    _num_cells=renumArr;
  else
    _num_nodes=renumArr;
}

const DataArrayInt *MEDFileCMesh::getFamilyFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception)
{
  if(meshDimRelToMaxExt!=0 && meshDimRelToMaxExt!=1)
    throw INTERP_KERNEL::Exception("MEDFileCMesh::getFamilyFieldAtLevel : Only available for levels 0 or 1 !");
  if(meshDimRelToMaxExt==0)
    return _fam_cells;
  else
    return _fam_nodes;
}

const DataArrayInt *MEDFileCMesh::getNumberFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception)
{
  if(meshDimRelToMaxExt!=0 && meshDimRelToMaxExt!=1)
    throw INTERP_KERNEL::Exception("MEDFileCMesh::getNumberFieldAtLevel : Only available for levels 0 or 1 !");
  if(meshDimRelToMaxExt==0)
    return _num_cells;
  else
    return _num_nodes;
}

const DataArrayInt *MEDFileCMesh::getRevNumberFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception)
{
  if(meshDimRelToMaxExt!=0 && meshDimRelToMaxExt!=1)
    throw INTERP_KERNEL::Exception("MEDFileCMesh::getRevNumberFieldAtLevel : Only available for levels 0 or 1 !");
  if(meshDimRelToMaxExt==0)
    {
      if((const DataArrayInt *)_num_cells)
        {
          int pos;
          int maxValue=_num_cells->getMaxValue(pos);
          _rev_num_cells=_num_cells->invertArrayN2O2O2N(maxValue+1);
          return _rev_num_cells;
        }
      else
        throw INTERP_KERNEL::Exception("MEDFileCMesh::getRevNumberFieldAtLevel : no cell renumbering for a request on reverse numbering !");
    }
  else
    {
      if((const DataArrayInt *)_num_nodes)
        {
          int pos;
          int maxValue=_num_nodes->getMaxValue(pos);
          _rev_num_nodes=_num_nodes->invertArrayN2O2O2N(maxValue+1);
          return _rev_num_nodes;
        }
      else
        throw INTERP_KERNEL::Exception("MEDFileCMesh::getRevNumberFieldAtLevel : no node renumbering for a request on reverse numbering !");
    }
}


MEDFileMeshMultiTS *MEDFileMeshMultiTS::New()
{
  return new MEDFileMeshMultiTS;
}

MEDFileMeshMultiTS *MEDFileMeshMultiTS::New(const char *fileName) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileMeshMultiTS(fileName);
}

MEDFileMeshMultiTS *MEDFileMeshMultiTS::New(const char *fileName, const char *mName) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileMeshMultiTS(fileName,mName);
}

const char *MEDFileMeshMultiTS::getName() const throw(INTERP_KERNEL::Exception)
{
  if(_mesh_one_ts.empty())
    throw INTERP_KERNEL::Exception("MEDFileMeshMultiTS::getName : no time steps set !");
  return _mesh_one_ts[0]->getName();
}

MEDFileMesh *MEDFileMeshMultiTS::getOneTimeStep() const throw(INTERP_KERNEL::Exception)
{
  if(_mesh_one_ts.empty())
    throw INTERP_KERNEL::Exception("MEDFileMeshMultiTS::getOneTimeStep : empty time step set !");
  return const_cast<MEDFileMesh *>(static_cast<const MEDFileMesh *>(_mesh_one_ts[0]));
}

void MEDFileMeshMultiTS::setOneTimeStep(MEDFileMesh *mesh1TimeStep) throw(INTERP_KERNEL::Exception)
{
  if(!mesh1TimeStep)
    throw INTERP_KERNEL::Exception("MEDFileMeshMultiTS::setOneTimeStep : input pointer should be different from 0 !");
  _mesh_one_ts.resize(1);
  mesh1TimeStep->incrRef();
  //MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> toto=mesh1TimeStep;
  _mesh_one_ts[0]=mesh1TimeStep;
}

void MEDFileMeshMultiTS::write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> >::const_iterator it=_mesh_one_ts.begin();it!=_mesh_one_ts.end();it++)
    {
      (*it)->copyOptionsFrom(*this);
      (*it)->write(fileName,mode);
    }
}

void MEDFileMeshMultiTS::loadFromFile(const char *fileName, const char *mName) throw(INTERP_KERNEL::Exception)
{//for the moment to be improved
  _mesh_one_ts.resize(1);
  _mesh_one_ts[0]=MEDFileMesh::New(fileName,mName,-1,-1);
}

MEDFileMeshMultiTS::MEDFileMeshMultiTS()
{
}

MEDFileMeshMultiTS::MEDFileMeshMultiTS(const char *fileName) throw(INTERP_KERNEL::Exception)
try
  {
    std::vector<std::string> ms=MEDLoader::GetMeshNames(fileName);
    if(ms.empty())
    {
      std::ostringstream oss; oss << "MEDFileUMesh::New : no meshes in file \"" << fileName << "\" !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
    MEDFileUtilities::CheckFileForRead(fileName);
    MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
    int dt,it;
    ParaMEDMEM::MEDCouplingMeshType meshType;
    std::string dummy2;
    MEDFileMeshL2::GetMeshIdFromName(fid,ms.front().c_str(),meshType,dt,it,dummy2);
    loadFromFile(fileName,ms.front().c_str());
  }
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

MEDFileMeshMultiTS::MEDFileMeshMultiTS(const char *fileName, const char *mName) throw(INTERP_KERNEL::Exception)
try
  {
    loadFromFile(fileName,mName);
  }
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

MEDFileMeshes *MEDFileMeshes::New()
{
  return new MEDFileMeshes;
}

MEDFileMeshes *MEDFileMeshes::New(const char *fileName) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileMeshes(fileName);
}

void MEDFileMeshes::write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception)
{
  checkCoherency();
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileMeshMultiTS> >::const_iterator it=_meshes.begin();it!=_meshes.end();it++)
    {
      (*it)->copyOptionsFrom(*this);
      (*it)->write(fileName,mode);
    }
}

int MEDFileMeshes::getNumberOfMeshes() const throw(INTERP_KERNEL::Exception)
{
  return _meshes.size();
}

MEDFileMesh *MEDFileMeshes::getMeshAtPos(int i) const throw(INTERP_KERNEL::Exception)
{
  if(i<0 || i>=(int)_meshes.size())
    {
      std::ostringstream oss; oss << "MEDFileMeshes::getMeshAtPos : invalid mesh id given in parameter ! Should be in [0;" << _meshes.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return _meshes[i]->getOneTimeStep();
}

MEDFileMesh *MEDFileMeshes::getMeshWithName(const char *mname) const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> ms=getMeshesNames();
  std::vector<std::string>::iterator it=std::find(ms.begin(),ms.end(),mname);
  if(it==ms.end())
    {
      std::ostringstream oss; oss << "MEDFileMeshes::getMeshWithName : Mesh  \"" << mname << "\" does not exist in this ! Existing are : ";
      std::copy(ms.begin(),ms.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return getMeshAtPos((int)std::distance(ms.begin(),it));
}

std::vector<std::string> MEDFileMeshes::getMeshesNames() const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> ret(_meshes.size());
  int i=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileMeshMultiTS> >::const_iterator it=_meshes.begin();it!=_meshes.end();it++,i++)
    {
      const MEDFileMeshMultiTS *f=(*it);
      if(f)
        {
          ret[i]=f->getName();
        }
      else
        {
          std::ostringstream oss; oss << "MEDFileMeshes::getMeshesNames : At rank #" << i << " mesh is not defined !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  return ret;
}

void MEDFileMeshes::resize(int newSize) throw(INTERP_KERNEL::Exception)
{
  _meshes.resize(newSize);
}

void MEDFileMeshes::pushMesh(MEDFileMesh *mesh) throw(INTERP_KERNEL::Exception)
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDFileMeshes::pushMesh : invalid input pointer ! should be different from 0 !");
  MEDFileMeshMultiTS *elt=MEDFileMeshMultiTS::New();
  elt->setOneTimeStep(mesh);
  _meshes.push_back(elt);
}

void MEDFileMeshes::setMeshAtPos(int i, MEDFileMesh *mesh) throw(INTERP_KERNEL::Exception)
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDFileMeshes::setMeshAtPos : invalid input pointer ! should be different from 0 !");
  if(i>=(int)_meshes.size())
    _meshes.resize(i+1);
  MEDFileMeshMultiTS *elt=MEDFileMeshMultiTS::New();
  elt->setOneTimeStep(mesh);
  _meshes[i]=elt;
}

void MEDFileMeshes::destroyMeshAtPos(int i) throw(INTERP_KERNEL::Exception)
{
  if(i<0 || i>=(int)_meshes.size())
    {
      std::ostringstream oss; oss << "MEDFileMeshes::destroyMeshAtPos : Invalid given id in input (" << i << ") should be in [0," << _meshes.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  _meshes.erase(_meshes.begin()+i);
}

void MEDFileMeshes::loadFromFile(const char *fileName) throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> ms=MEDLoader::GetMeshNames(fileName);
  int i=0;
  _meshes.resize(ms.size());
  for(std::vector<std::string>::const_iterator it=ms.begin();it!=ms.end();it++,i++)
    _meshes[i]=MEDFileMeshMultiTS::New(fileName,(*it).c_str());
}

MEDFileMeshes::MEDFileMeshes()
{
}

MEDFileMeshes::MEDFileMeshes(const char *fileName) throw(INTERP_KERNEL::Exception)
try
  {
    loadFromFile(fileName);
  }
catch(INTERP_KERNEL::Exception& e)
  {
  }

std::string MEDFileMeshes::simpleRepr() const
{
  std::ostringstream oss;
  oss << "(*****************)\n(* MEDFileMeshes *)\n(*****************)\n\n";
  simpleReprWithoutHeader(oss);
  return oss.str();
}

void MEDFileMeshes::simpleReprWithoutHeader(std::ostream& oss) const
{
  int nbOfMeshes=getNumberOfMeshes();
  oss << "There are " << nbOfMeshes << " meshes with the following names : \n";
  std::vector<std::string> mns=getMeshesNames();
  for(int i=0;i<nbOfMeshes;i++)
    oss << "  - #" << i << " \"" << mns[i] << "\"\n";
}

void MEDFileMeshes::checkCoherency() const throw(INTERP_KERNEL::Exception)
{
  static const char MSG[]="MEDFileMeshes::checkCoherency : mesh at rank ";
  int i=0;
  std::set<std::string> s;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileMeshMultiTS> >::const_iterator it=_meshes.begin();it!=_meshes.end();it++,i++)
    {
      const MEDFileMeshMultiTS *elt=(*it);
      if(!elt)
        {
          std::ostringstream oss; oss << MSG << i << "/" << _meshes.size() << " is empty !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      std::size_t sz=s.size();
      s.insert(std::string((*it)->getName()));
      if(s.size()==sz)
        {
          std::ostringstream oss; oss << MSG << i << " has a name (\"" << (*it)->getName() << "\") already used by an another mesh in list !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
}
