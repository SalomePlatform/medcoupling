// Copyright (C) 2007-2015  CEA/DEN, EDF R&D
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

#include "MEDFileMesh.hxx"
#include "MEDFileUtilities.hxx"
#include "MEDFileFieldOverView.hxx"
#include "MEDFileField.hxx"
#include "MEDLoader.hxx"
#include "MEDFileSafeCaller.txx"
#include "MEDLoaderBase.hxx"

#include "MEDCouplingUMesh.hxx"

#include "InterpKernelAutoPtr.hxx"

#include <limits>
#include <cmath>

extern med_geometry_type typmai3[34];

using namespace ParaMEDMEM;

const char MEDFileMesh::DFT_FAM_NAME[]="FAMILLE_ZERO";

MEDFileMesh::MEDFileMesh():_order(-1),_iteration(-1),_time(0.),_univ_wr_status(true)
{
}

std::size_t MEDFileMesh::getHeapMemorySizeWithoutChildren() const
{
  std::size_t ret(_dt_unit.capacity()+_name.capacity()+_univ_name.capacity()+_desc_name.capacity());
  for(std::map<std::string, std::vector<std::string> >::const_iterator it=_groups.begin();it!=_groups.end();it++)
    {
      ret+=(*it).first.capacity()+(*it).second.capacity()*sizeof(std::string);
      for(std::vector<std::string>::const_iterator it2=(*it).second.begin();it2!=(*it).second.end();it2++)
        ret+=(*it2).capacity();
    }
  for(std::map<std::string,int>::const_iterator it=_families.begin();it!=_families.end();it++)
    ret+=(*it).first.capacity()+sizeof(int);
  return ret;
}

std::vector<const BigMemoryObject *> MEDFileMesh::getDirectChildrenWithNull() const
{
  return std::vector<const BigMemoryObject *>();
}

/*!
 * Returns a new MEDFileMesh holding the mesh data that has been read from a given MED
 * file. The first mesh in the file is loaded.
 *  \param [in] fileName - the name of MED file to read.
 *  \return MEDFileMesh * - a new instance of MEDFileMesh. The caller is to delete this
 *          mesh using decrRef() as it is no more needed. 
 *  \throw If the file is not readable.
 *  \throw If there is no meshes in the file.
 *  \throw If the mesh in the file is of a not supported type.
 */
MEDFileMesh *MEDFileMesh::New(const std::string& fileName, MEDFileMeshReadSelector *mrs)
{
  std::vector<std::string> ms=MEDLoader::GetMeshNames(fileName);
  if(ms.empty())
    {
      std::ostringstream oss; oss << "MEDFileMesh::New : no meshes in file \"" << fileName << "\" !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  MEDFileUtilities::CheckFileForRead(fileName);
  ParaMEDMEM::MEDCouplingMeshType meshType;
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName.c_str(),MED_ACC_RDONLY);
  int dt,it;
  std::string dummy2;
  MEDFileMeshL2::GetMeshIdFromName(fid,ms.front(),meshType,dt,it,dummy2);
  switch(meshType)
  {
    case UNSTRUCTURED:
      {
        MEDCouplingAutoRefCountObjectPtr<MEDFileUMesh> ret=MEDFileUMesh::New();
        ret->loadUMeshFromFile(fid,ms.front(),dt,it,mrs);
        ret->loadJointsFromFile(fid);
        return (MEDFileUMesh *)ret.retn();
      }
    case CARTESIAN:
      {
        MEDCouplingAutoRefCountObjectPtr<MEDFileCMesh> ret=MEDFileCMesh::New();
        ret->loadCMeshFromFile(fid,ms.front(),dt,it,mrs);
        ret->loadJointsFromFile(fid);
        return (MEDFileCMesh *)ret.retn();
      }
    case CURVE_LINEAR:
      {
        MEDCouplingAutoRefCountObjectPtr<MEDFileCurveLinearMesh> ret=MEDFileCurveLinearMesh::New();
        ret->loadCLMeshFromFile(fid,ms.front(),dt,it,mrs);
        ret->loadJointsFromFile(fid);
        return (MEDFileCurveLinearMesh *)ret.retn();
      }
    default:
      {
        std::ostringstream oss; oss << "MEDFileMesh::New : MED file exists and has mesh '" << ms.front() << "' exists but unsupported type yet !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
  }
}

/*!
 * Returns a new MEDFileMesh holding the mesh data that has been read from a given MED
 * file. The mesh to load is specified by its name and numbers of a time step and an
 * iteration.
 *  \param [in] fileName - the name of MED file to read.
 *  \param [in] mName - the name of the mesh to read.
 *  \param [in] dt - the number of a time step.
 *  \param [in] it - the number of an iteration.
 *  \param [in] joints - the sub-domain joints to use instead of those that can be read
 *          from the MED file. Usually this joints are those just read by another iteration
 *          of mName mesh, when this method is called by MEDFileMeshMultiTS::New()
 *  \return MEDFileMesh * - a new instance of MEDFileMesh. The caller is to delete this
 *          mesh using decrRef() as it is no more needed. 
 *  \throw If the file is not readable.
 *  \throw If there is no mesh with given attributes in the file.
 *  \throw If the mesh in the file is of a not supported type.
 */
MEDFileMesh *MEDFileMesh::New(const std::string& fileName, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs, MEDFileJoints* joints)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  ParaMEDMEM::MEDCouplingMeshType meshType;
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName.c_str(),MED_ACC_RDONLY);
  int dummy0,dummy1;
  std::string dummy2;
  MEDFileMeshL2::GetMeshIdFromName(fid,mName,meshType,dummy0,dummy1,dummy2);
  switch(meshType)
  {
    case UNSTRUCTURED:
      {
        MEDCouplingAutoRefCountObjectPtr<MEDFileUMesh> ret=MEDFileUMesh::New();
        ret->loadUMeshFromFile(fid,mName,dt,it,mrs);
        ret->loadJointsFromFile(fid,joints);
        return (MEDFileUMesh *)ret.retn();
      }
    case CARTESIAN:
      {
        MEDCouplingAutoRefCountObjectPtr<MEDFileCMesh> ret=MEDFileCMesh::New();
        ret->loadCMeshFromFile(fid,mName,dt,it,mrs);
        ret->loadJointsFromFile(fid,joints);
        return (MEDFileCMesh *)ret.retn();
      }
    case CURVE_LINEAR:
      {
        MEDCouplingAutoRefCountObjectPtr<MEDFileCurveLinearMesh> ret=MEDFileCurveLinearMesh::New();
        ret->loadCLMeshFromFile(fid,mName,dt,it,mrs);
        ret->loadJointsFromFile(fid,joints);
        return (MEDFileCurveLinearMesh *)ret.retn();
      }
    default:
      {
        std::ostringstream oss; oss << "MEDFileMesh::New : MED file exists and has mesh '" << mName << "' exists but unsupported type yet !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
  }
}

/*!
 * Writes \a this mesh into an open MED file specified by its descriptor.
 *  \param [in] fid - the MED file descriptor.
 *  \throw If the mesh name is not set.
 *  \throw If the file is open for reading only.
 *  \throw If the writing mode == 1 and the same data is present in an existing file.
 */
void MEDFileMesh::write(med_idt fid) const
{
  if(!existsFamily(0))
    const_cast<MEDFileMesh *>(this)->addFamily(DFT_FAM_NAME,0);
  if(_name.empty())
    throw INTERP_KERNEL::Exception("MEDFileMesh : name is empty. MED file ask for a NON EMPTY name !");
  writeLL(fid);
}

/*!
 * Writes \a this mesh into a MED file specified by its name.
 *  \param [in] fileName - the MED file name.
 *  \param [in] mode - the writing mode. For more on \a mode, see \ref AdvMEDLoaderBasics.
 * - 2 - erase; an existing file is removed.
 * - 1 - append; same data should not be present in an existing file.
 * - 0 - overwrite; same data present in an existing file is overwritten.
 *  \throw If the mesh name is not set.
 *  \throw If \a mode == 1 and the same data is present in an existing file.
 */
void MEDFileMesh::write(const std::string& fileName, int mode) const
{
  med_access_mode medmod=MEDFileUtilities::TraduceWriteMode(mode);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName.c_str(),medmod);
  std::ostringstream oss; oss << "MEDFileMesh : error on attempt to write in file : \"" << fileName << "\""; 
  MEDFileUtilities::CheckMEDCode(fid,fid,oss.str());
  write(fid);
}

/*!
 * Checks if \a this and another mesh are equal.
 *  \param [in] other - the mesh to compare with.
 *  \param [in] eps - a precision used to compare real values.
 *  \param [in,out] what - the string returning description of unequal data.
 *  \return bool - \c true if the meshes are equal, \c false, else.
 */
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
  //univ_name has been ignored -> not a bug because it is a mutable attribute
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

void MEDFileMesh::setName(const std::string& name)
{
  _name=name;
}

/*!
 * Clears redundant attributes of incorporated data arrays.
 */
void MEDFileMesh::clearNonDiscrAttributes() const
{

}

bool MEDFileMesh::changeNames(const std::vector< std::pair<std::string,std::string> >& modifTab)
{
  for(std::vector< std::pair<std::string,std::string> >::const_iterator it=modifTab.begin();it!=modifTab.end();it++)
    {
      if((*it).first==_name)
        {
          _name=(*it).second;
          return true;
        }
    }
  return false;
}

/*!
 * Copies data on groups and families from another mesh.
 *  \param [in] other - the mesh to copy the data from.
 */
void MEDFileMesh::copyFamGrpMapsFrom(const MEDFileMesh& other)
{
  _groups=other._groups;
  _families=other._families;
}


/*!
 * This method clear all the groups in the map.
 * So this method does not operate at all on arrays.
 * So this method can lead to orphan families.
 * 
 * \sa MEDFileMesh::clearFamMap, MEDFileMesh::clearFamGrpMaps
 */
void MEDFileMesh::clearGrpMap()
{
  _groups.clear();
}

/*!
 * This method clear all the families in the map.
 * So this method does not operate at all on arrays.
 * WARNING ! if there are some groups lying on cleared families, those groups will be impacted !
 *
 * \sa MEDFileMesh::clearFamMap, MEDFileMesh::clearFamGrpMaps
 */
void MEDFileMesh::clearFamMap()
{
  _families.clear();
}

/*!
 * This method clear all the families and groups in the map.
 * So this method does not operate at all on arrays.
 * As all groups and families entry will be removed after 
 * the call of MEDFileMesh::setFamilyFieldArr method with 0 or None (python) in the 2nd parameter can be useful to reduce the size of the object.
 *
 * \sa MEDFileMesh::clearFamMap, MEDFileMesh::clearFamMap
 */
void MEDFileMesh::clearFamGrpMaps()
{
  clearGrpMap();
  clearFamMap();
}

/*!
 * Returns names of families constituting a group.
 *  \param [in] name - the name of the group of interest.
 *  \return std::vector<std::string> - a sequence of names of the families.
 *  \throw If the name of a nonexistent group is specified.
 */
std::vector<std::string> MEDFileMesh::getFamiliesOnGroup(const std::string& name) const
{
  std::string oname(name);
  std::map<std::string, std::vector<std::string> >::const_iterator it=_groups.find(oname);
  if(it==_groups.end())
    {
      std::vector<std::string> grps=getGroupsNames();
      std::ostringstream oss; oss << "No such groupname \"" << name << "\" !\nAvailable groups are :";
      std::copy(grps.begin(),grps.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return (*it).second;
}

/*!
 * Returns names of families constituting some groups.
 *  \param [in] grps - a sequence of names of groups of interest.
 *  \return std::vector<std::string> - a sequence of names of the families.
 *  \throw If a name of a nonexistent group is present in \a grps.
 */
std::vector<std::string> MEDFileMesh::getFamiliesOnGroups(const std::vector<std::string>& grps) const
{
  std::set<std::string> fams;
  for(std::vector<std::string>::const_iterator it=grps.begin();it!=grps.end();it++)
    {
      std::map<std::string, std::vector<std::string> >::const_iterator it2=_groups.find(*it);
      if(it2==_groups.end())
        {
          std::ostringstream oss; oss << "No such group in mesh \"" << _name << "\" : " << *it; 
          std::vector<std::string> grps2=getGroupsNames(); oss << "\" !\nAvailable groups are :";
          std::copy(grps2.begin(),grps2.end(),std::ostream_iterator<std::string>(oss," "));
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      fams.insert((*it2).second.begin(),(*it2).second.end());
    }
  std::vector<std::string> fams2(fams.begin(),fams.end());
  return fams2;
}

/*!
 * Returns ids of families constituting a group.
 *  \param [in] name - the name of the group of interest.
 *  \return std::vector<int> - sequence of ids of the families.
 *  \throw If the name of a nonexistent group is specified.
 */
std::vector<int> MEDFileMesh::getFamiliesIdsOnGroup(const std::string& name) const
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
 * Sets names of families constituting a group. If data on families of this group is
 * already present, it is overwritten. Every family in \a fams is checked, and if a
 family is not yet in \a this mesh, the default group id \c 0 is assigned to it.
 *  \param [in] name - the name of the group of interest.
 *  \param [in] fams - a sequence of names of families constituting the group.
 */
void MEDFileMesh::setFamiliesOnGroup(const std::string& name, const std::vector<std::string>& fams)
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
 * Sets families constituting a group. The families are specified by their ids.
 * If a family name is not found by its id, an exception is thrown.
 * If several families have same id, the first one in lexical order is taken.
 *  \param [in] name - the name of the group of interest.
 *  \param [in] famIds - a sequence of ids of families constituting the group.
 *  \throw If a family name is not found by its id.
 */
void MEDFileMesh::setFamiliesIdsOnGroup(const std::string& name, const std::vector<int>& famIds)
{
  std::string oname(name);
  std::vector<std::string> fams(famIds.size());
  int i=0;
  for(std::vector<int>::const_iterator it1=famIds.begin();it1!=famIds.end();it1++,i++)
    {
      std::string name2=getFamilyNameGivenId(*it1);
      fams[i]=name2;
    }
  _groups[oname]=fams;
}

/*!
 * Returns names of groups including a given family.
 *  \param [in] name - the name of the family of interest.
 *  \return std::vector<std::string> - a sequence of names of groups including the family.
 */
std::vector<std::string> MEDFileMesh::getGroupsOnFamily(const std::string& name) const
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
 * Adds an existing family to groups.
 *  \param [in] famName - a name of family to add to \a grps.
 *  \param [in] grps - a sequence of group names to add the family in.
 *  \throw If a family named \a famName not yet exists.
 */
void MEDFileMesh::setGroupsOnFamily(const std::string& famName, const std::vector<std::string>& grps)
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
  for(std::vector<std::string>::const_iterator it3=grps.begin();it3!=grps.end();it3++)
    {
      std::map< std::string, std::vector<std::string> >::iterator it2=_groups.find(*it3);
      if(it2!=_groups.end())
        (*it2).second.push_back(fName);
      else
        {
          std::vector<std::string> grps2(1,fName);
          _groups[*it3]=grps2;
        }
    }
}

/*!
 * Returns names of all groups of \a this mesh.
 *  \return std::vector<std::string> - a sequence of group names.
 */
std::vector<std::string> MEDFileMesh::getGroupsNames() const
{
  std::vector<std::string> ret(_groups.size());
  int i=0;
  for(std::map<std::string, std::vector<std::string> >::const_iterator it=_groups.begin();it!=_groups.end();it++,i++)
    ret[i]=(*it).first;
  return ret;
}

/*!
 * Returns names of all families of \a this mesh.
 *  \return std::vector<std::string> - a sequence of family names.
 */
std::vector<std::string> MEDFileMesh::getFamiliesNames() const
{
  std::vector<std::string> ret(_families.size());
  int i=0;
  for(std::map<std::string, int >::const_iterator it=_families.begin();it!=_families.end();it++,i++)
    ret[i]=(*it).first;
  return ret;
}

/*!
 * Changes a name of every family, included in one group only, to be same as the group name.
 *  \throw If there are families with equal names in \a this mesh.
 */
void MEDFileMesh::assignFamilyNameWithGroupName()
{
  std::map<std::string, std::vector<std::string> > groups(_groups);
  std::map<std::string,int> newFams;
  for(std::map<std::string,int>::const_iterator it=_families.begin();it!=_families.end();it++)
    {
      std::vector<std::string> grps=getGroupsOnFamily((*it).first);
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

/*!
 * Removes all groups lying on no family. If there is no empty groups, \a this is let untouched.
 * 
 * \return the removed groups.
 */
std::vector<std::string> MEDFileMesh::removeEmptyGroups()
{
  std::vector<std::string> ret;
  std::map<std::string, std::vector<std::string> > newGrps;
  for(std::map<std::string, std::vector<std::string> >::const_iterator it=_groups.begin();it!=_groups.end();it++)
    {
      if((*it).second.empty())
        ret.push_back((*it).first);
      else
        newGrps[(*it).first]=(*it).second;
    }
  if(!ret.empty())
    _groups=newGrps;
  return ret;
}

/*!
 * Removes a group from \a this mesh.
 *  \param [in] name - the name of the group to remove.
 *  \throw If no group with such a \a name exists.
 */
void MEDFileMesh::removeGroup(const std::string& name)
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

/*!
 * Removes a family from \a this mesh.
 *  \param [in] name - the name of the family to remove.
 *  \throw If no family with such a \a name exists.
 */
void MEDFileMesh::removeFamily(const std::string& name)
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
  for(std::map<std::string, std::vector<std::string> >::iterator it3=_groups.begin();it3!=_groups.end();it3++)
    {
      std::vector<std::string>& v=(*it3).second;
      std::vector<std::string>::iterator it4=std::find(v.begin(),v.end(),oname);
      if(it4!=v.end())
        v.erase(it4);
    }
}

/*!
 * Removes all groups in \a this that are orphan. A group is orphan if this group lies on
 * a set of families, themselves orphan. A family is said orphan if its id appears nowhere in
 * family field whatever its level. This method also suppresses the orphan families.
 * 
 * \return - The list of removed groups names. 
 *
 * \sa MEDFileMesh::removeOrphanFamilies.
 */
std::vector<std::string> MEDFileMesh::removeOrphanGroups()
{
  removeOrphanFamilies();
  return removeEmptyGroups();
}

/*!
 * Removes all families in \a this that are orphan. A family is said orphan if its id appears nowhere in
 * family field whatever its level. Groups are updated in consequence, that is to say all groups lying on orphan family, will see their families list modified.
 * 
 * \return - The list of removed families names.
 * \sa MEDFileMesh::removeOrphanGroups.
 */
std::vector<std::string> MEDFileMesh::removeOrphanFamilies()
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> allFamIdsInUse=computeAllFamilyIdsInUse();
  std::vector<std::string> ret;
  if(!((DataArrayInt*)allFamIdsInUse))
    {
      ret=getFamiliesNames();
      _families.clear(); _groups.clear();
      return ret;
    }
  std::map<std::string,int> famMap;
  std::map<std::string, std::vector<std::string> > grps(_groups);
  for(std::map<std::string,int>::const_iterator it=_families.begin();it!=_families.end();it++)
    {
      if(allFamIdsInUse->presenceOfValue((*it).second))
        famMap[(*it).first]=(*it).second;
      else
        {
          ret.push_back((*it).first);
          std::vector<std::string> grpsOnEraseFam=getGroupsOnFamily((*it).first);
          for(std::vector<std::string>::const_iterator it2=grpsOnEraseFam.begin();it2!=grpsOnEraseFam.end();it2++)
            {
              std::map<std::string, std::vector<std::string> >::iterator it3=grps.find(*it2);//it3!=grps.empty() thanks to copy
              std::vector<std::string>& famv=(*it3).second;
              std::vector<std::string>::iterator it4=std::find(famv.begin(),famv.end(),(*it).first);//it4!=famv.end() thanks to copy
              famv.erase(it4);
            }
        }
    }
  if(!ret.empty())
    { _families=famMap; _groups=grps; }
  return ret;
}

/*!
 * This method operates only on maps in \a this. The arrays are not considered here. So this method will remove a family (except "FAMILLE_ZERO" family) if no group lies on it whatever
 * this family is orphan or not.
 *
 * \warning this method is different from removeOrphanFamilies that scans family field array to find orphan families.
 */
void MEDFileMesh::removeFamiliesReferedByNoGroups()
{
  std::map<std::string,int> fams;
  std::set<std::string> sfams;
  for(std::map<std::string,int>::const_iterator it=_families.begin();it!=_families.end();it++)
    sfams.insert((*it).first);
  for(std::map<std::string, std::vector<std::string> >::const_iterator it0=_groups.begin();it0!=_groups.end();it0++)
    for(std::vector<std::string>::const_iterator it1=(*it0).second.begin();it1!=(*it0).second.end();it1++)
      sfams.erase(*it1);
  for(std::set<std::string>::const_iterator it=sfams.begin();it!=sfams.end();it++)
    if(*it!=DFT_FAM_NAME)
      _families.erase(*it);
}

/*!
 * This method has no impact on groups. This method only works on families. This method firstly removes families not refered by any groups in \a this, then all unused entities
 * are put as belonging to family 0 ("FAMILLE_ZERO"). Finally, all orphanFamilies are killed.
 * This method raises an exception if "FAMILLE_ZERO" is already belonging to a group.
 *
 * \sa MEDFileMesh::removeOrphanFamilies
 */
void MEDFileMesh::rearrangeFamilies()
{
  checkOrphanFamilyZero();
  removeFamiliesReferedByNoGroups();
  //
  std::vector<int> levels(getNonEmptyLevelsExt());
  std::set<int> idsRefed;
  for(std::map<std::string,int>::const_iterator it=_families.begin();it!=_families.end();it++)
    idsRefed.insert((*it).second);
  for(std::vector<int>::const_iterator it=levels.begin();it!=levels.end();it++)
    {
      const DataArrayInt *fams(0);
      try
      {
          fams=getFamilyFieldAtLevel(*it);
      }
      catch(INTERP_KERNEL::Exception& e) { }
      if(!fams)
        continue;
      std::vector<bool> v(fams->getNumberOfTuples(),false);
      for(std::set<int>::const_iterator pt=idsRefed.begin();pt!=idsRefed.end();pt++)
        fams->switchOnTupleEqualTo(*pt,v);
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> unfetchedIds(DataArrayInt::BuildListOfSwitchedOff(v));
      if(!unfetchedIds->empty())
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> newFams(fams->deepCpy());
          newFams->setPartOfValuesSimple3(0,unfetchedIds->begin(),unfetchedIds->end(),0,1,1);
          setFamilyFieldArr(*it,newFams);
        }
    }
  removeOrphanFamilies();
}

/*!
 * This method only checks that "FAMILLE_ZERO" is orphan (not belonging to a group).
 */
void MEDFileMesh::checkOrphanFamilyZero() const
{
  for(std::map<std::string, std::vector<std::string> >::const_iterator it=_groups.begin();it!=_groups.end();it++)
    {
      if(std::find((*it).second.begin(),(*it).second.end(),DFT_FAM_NAME)!=(*it).second.end())
        {
          std::ostringstream oss; oss << "MEDFileMesh::rearrangeFamilies : Groups \"" << (*it).first << "\" is lying on family \"" << DFT_FAM_NAME << "\" !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
}

/*!
 * Renames a group in \a this mesh.
 *  \param [in] oldName - a current name of the group to rename.
 *  \param [in] newName - a new group name.
 *  \throw If no group named \a oldName exists in \a this mesh.
 *  \throw If a group named \a newName already exists.
 */
void MEDFileMesh::changeGroupName(const std::string& oldName, const std::string& newName)
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
  std::map<std::string, std::vector<std::string> >::iterator it2=_groups.find(nname);
  if(it2!=_groups.end())
    {
      std::ostringstream oss; oss << "Such groupname \"" << newName << "\" already exists ! Kill it before !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  std::vector<std::string> cpy=(*it).second;
  _groups.erase(it);
  _groups[newName]=cpy;
}

/*!
 * Changes an id of a family in \a this mesh. 
 * This method calls changeFamilyIdArr().
 *  \param [in] oldId - a current id of the family.
 *  \param [in] newId - a new family id.
 */
void MEDFileMesh::changeFamilyId(int oldId, int newId)
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

/*!
 * Renames a family in \a this mesh.
 *  \param [in] oldName - a current name of the family to rename.
 *  \param [in] newName - a new family name.
 *  \throw If no family named \a oldName exists in \a this mesh.
 *  \throw If a family named \a newName already exists.
 */
void MEDFileMesh::changeFamilyName(const std::string& oldName, const std::string& newName)
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
  for(std::map<std::string, std::vector<std::string> >::iterator it3=_groups.begin();it3!=_groups.end();it3++)
    {
      std::vector<std::string>& v=(*it3).second;
      std::vector<std::string>::iterator it4=std::find(v.begin(),v.end(),oname);
      if(it4!=v.end())
        (*it4)=nname;
    }
}

/*!
 * Checks if \a this and another mesh contains the same families.
 *  \param [in] other - the mesh to compare with \a this one.
 *  \param [in,out] what - an unused parameter.
 *  \return bool - \c true if number of families and their ids are the same in the two
 *          meshes. Families with the id == \c 0 are not considered.
 */
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

/*!
 * Checks if \a this and another mesh contains the same groups.
 *  \param [in] other - the mesh to compare with \a this one.
 *  \param [in,out] what - a string describing a difference of groups of the two meshes
 *          in case if this method returns \c false.
 *  \return bool - \c true if number of groups and families constituting them are the
 *          same in the two meshes.
 */
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

/*!
 * Checks if a group with a given name exists in \a this mesh.
 *  \param [in] groupName - the group name.
 *  \return bool - \c true the group \a groupName exists in \a this mesh.
 */
bool MEDFileMesh::existsGroup(const std::string& groupName) const
{
  std::string grpName(groupName);
  return _groups.find(grpName)!=_groups.end();
}

/*!
 * Checks if a family with a given id exists in \a this mesh.
 *  \param [in] famId - the family id.
 *  \return bool - \c true the family with the id \a famId exists in \a this mesh.
 */
bool MEDFileMesh::existsFamily(int famId) const
{
  for(std::map<std::string,int>::const_iterator it2=_families.begin();it2!=_families.end();it2++)
    if((*it2).second==famId)
      return true;
  return false;
}

/*!
 * Checks if a family with a given name exists in \a this mesh.
 *  \param [in] familyName - the family name.
 *  \return bool - \c true the family \a familyName exists in \a this mesh.
 */
bool MEDFileMesh::existsFamily(const std::string& familyName) const
{
  std::string fname(familyName);
  return _families.find(fname)!=_families.end();
}

/*!
 * Sets an id of a family.
 *  \param [in] familyName - the family name.
 *  \param [in] id - a new id of the family.
 */
void MEDFileMesh::setFamilyId(const std::string& familyName, int id)
{
  std::string fname(familyName);
  _families[fname]=id;
}

void MEDFileMesh::setFamilyIdUnique(const std::string& familyName, int id)
{
  std::string fname(familyName);
  for(std::map<std::string,int>::const_iterator it=_families.begin();it!=_families.end();it++)
    if((*it).second==id)
      {
        if((*it).first!=familyName)
          {
            std::ostringstream oss; oss << "MEDFileMesh::setFamilyIdUnique : Family id #" << id << " is already belonging to family with name \"" << (*it).first << "\" !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
  _families[fname]=id;
}

/*!
 * Adds a family to \a this mesh.
 *  \param [in] familyName - a name of the family.
 *  \param [in] famId - an id of the family.
 *  \throw If a family with the same name or id already exists in \a this mesh.
 */
void MEDFileMesh::addFamily(const std::string& familyName, int famId)
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

/*!
 * Creates a group including all mesh entities of given dimension.
 * \warning This method does \b not guarantee that the created group includes mesh
 * entities of only \a meshDimRelToMaxExt dimension in the case if some family id is
 * present in family fields of different dimensions. To assure this, call
 * ensureDifferentFamIdsPerLevel() \b before calling this method.
 *  \param [in] meshDimRelToMaxExt - a relative dimension of mesh entities to include to
 *          the group.
 *  \param [in] groupName - a name of the new group.
 *  \throw If a group named \a groupName already exists.
 *  \throw If no mesh entities of dimension \a meshDimRelToMaxExt exist in \a this mesh.
 *  \throw If no family field of dimension \a meshDimRelToMaxExt is present in \a this mesh.
 */
void MEDFileMesh::createGroupOnAll(int meshDimRelToMaxExt, const std::string& groupName)
{
  std::string grpName(groupName);
  std::vector<int> levs=getNonEmptyLevelsExt();
  if(std::find(levs.begin(),levs.end(),meshDimRelToMaxExt)==levs.end())
    {
      std::ostringstream oss; oss << "MEDFileMesh::createGroupOnAll : The relative ext dimension " << meshDimRelToMaxExt << " is not available !" << std::endl;
      oss << "Available relative ext levels are : ";
      std::copy(levs.begin(),levs.end(),std::ostream_iterator<int>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  if(existsGroup(groupName))
    {
      std::ostringstream oss; oss << "MEDFileMesh::createGroupOnAll : The groups \"" << grpName << "\" already exists in this !" << std::endl;
      oss << "Already existing groups are : ";
      std::copy(levs.begin(),levs.end(),std::ostream_iterator<int>(oss," "));
      oss << std::endl << "Please choose an another group name or call removeGroup(\"" << grpName << "\") method !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  const DataArrayInt *fieldFamIds=getFamilyFieldAtLevel(meshDimRelToMaxExt);
  if(fieldFamIds==0)
    throw INTERP_KERNEL::Exception("MEDFileMesh::createGroupOnAll : Family field arr ids is not defined for this level !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> famIds=fieldFamIds->getDifferentValues();
  std::vector<std::string> familiesOnWholeGroup;
  for(const int *it=famIds->begin();it!=famIds->end();it++)
    {
      bool tmp;
      familiesOnWholeGroup.push_back(findOrCreateAndGiveFamilyWithId(*it,tmp));
    }
  _groups[grpName]=familiesOnWholeGroup;
}

/*!
 * Ensures that given family ids do not present in family fields of dimensions different
 * than given ones. If a family id is present in the family fields of dimensions different
 * than the given ones, a new family is created and the whole data is updated accordingly.
 *  \param [in] famIds - a sequence of family ids to check.
 *  \param [in] vMeshDimRelToMaxExt - a sequence of relative dimensions to which the \a
 *          famIds should exclusively belong.
 *  \return bool - \c true if no modification is done in \a this mesh by this method.
 */
bool MEDFileMesh::keepFamIdsOnlyOnLevs(const std::vector<int>& famIds, const std::vector<int>& vMeshDimRelToMaxExt)
{
  std::set<int> levsInput(vMeshDimRelToMaxExt.begin(),vMeshDimRelToMaxExt.end());
  std::vector<int> levs=getNonEmptyLevelsExt();
  std::set<int> levs2(levs.begin(),levs.end());
  std::vector<int> levsToTest;
  std::set_difference(levs2.begin(),levs2.end(),levsInput.begin(),levsInput.end(),std::back_insert_iterator< std::vector<int> >(levsToTest));
  std::set<int> famIds2(famIds.begin(),famIds.end());
  bool ret=true;
  int maxFamId=1;
  if(!_families.empty())
    maxFamId=getMaxFamilyId()+1;
  std::vector<std::string> allFams=getFamiliesNames();
  for(std::vector<int>::const_iterator it=levsToTest.begin();it!=levsToTest.end();it++)
    {
      const DataArrayInt *fieldFamIds=getFamilyFieldAtLevel(*it);
      if(fieldFamIds)
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> famIds3=fieldFamIds->getDifferentValues();
          std::vector<int> tmp;
          std::set_intersection(famIds3->begin(),famIds3->end(),famIds2.begin(),famIds2.end(),std::back_insert_iterator< std::vector<int> >(tmp));
          for(std::vector<int>::const_iterator it2=tmp.begin();it2!=tmp.end();it2++)
            {
              ret=false;
              std::string famName=getFamilyNameGivenId(*it2);
              std::ostringstream oss; oss << "Family_" << maxFamId;
              std::string zeName=CreateNameNotIn(oss.str(),allFams);
              addFamilyOnAllGroupsHaving(famName,zeName);
              _families[zeName]=maxFamId;
              (const_cast<DataArrayInt *>(fieldFamIds))->changeValue(*it2,maxFamId);
              maxFamId++;
            }
        }
    }
  return ret;
}

/*!
 * Adds a family to a given group in \a this mesh. If the group with a given name does
 * not exist, it is created.
 *  \param [in] grpName - the name of the group to add the family in.
 *  \param [in] famName - the name of the family to add to the group named \a grpName.
 *  \throw If \a grpName or \a famName is an empty string.
 *  \throw If no family named \a famName is present in \a this mesh.
 */
void MEDFileMesh::addFamilyOnGrp(const std::string& grpName, const std::string& famName)
{
  std::string grpn(grpName);
  std::string famn(famName);
  if(grpn.empty() || famn.empty())
    throw INTERP_KERNEL::Exception("MEDFileMesh::addFamilyOnGrp : input strings must be non null !");
  std::vector<std::string> fams=getFamiliesNames();
  if(std::find(fams.begin(),fams.end(),famn)==fams.end())
    {
      std::ostringstream oss; oss << "MEDFileMesh::addFamilyOnGrp : Family \"" << famn << "\" does not exist !" << std::endl;
      oss << "Create this family or choose an existing one ! Existing fams are : ";
      std::copy(fams.begin(),fams.end(),std::ostream_iterator<std::string>(oss," ")); oss << ".";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
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

/*!
 * This method adds to all groups lying on family with name 'famName' the other family name 'otherFamName'.
 * This method is quite underground because it can lead to unconsistency because family 'otherFamName' is \b not added into _families.
 * This method is used by MEDFileMesh::keepFamIdsOnlyOnLevs method.
 */
void MEDFileMesh::addFamilyOnAllGroupsHaving(const std::string& famName, const std::string& otherFamName)
{
  std::string famNameCpp(famName);
  std::string otherCpp(otherFamName);
  for(std::map<std::string, std::vector<std::string> >::iterator it=_groups.begin();it!=_groups.end();it++)
    {
      std::vector<std::string>& v=(*it).second;
      if(std::find(v.begin(),v.end(),famNameCpp)!=v.end())
        {
          v.push_back(otherCpp);
        }
    }
}

/*!
 * \param [in] ids ids and group name of the new group to add. The ids should be sorted and different each other (MED file norm).
 * \parma [in,out] famArr family array on level of interest to be renumbered. The input pointer should be not \c NULL (no check of that will be performed)
 */
void MEDFileMesh::addGroupUnderground(bool isNodeGroup, const DataArrayInt *ids, DataArrayInt *famArr)
{
  if(!ids)
    throw INTERP_KERNEL::Exception("MEDFileUMesh::addGroup : NULL pointer in input !");
  std::string grpName(ids->getName());
  if(grpName.empty())
    throw INTERP_KERNEL::Exception("MEDFileUMesh::addGroup : empty group name ! MED file format do not accept empty group name !");
  ids->checkStrictlyMonotonic(true);
  famArr->incrRef(); MEDCouplingAutoRefCountObjectPtr<DataArrayInt> famArrTmp(famArr);
  std::vector<std::string> grpsNames=getGroupsNames();
  if(std::find(grpsNames.begin(),grpsNames.end(),grpName)!=grpsNames.end())
    {
      std::ostringstream oss; oss << "MEDFileUMesh::addGroup : Group with name \"" << grpName << "\" already exists ! Destroy it before calling this method !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  std::list< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > allFamIds(getAllNonNullFamilyIds());
  allFamIds.erase(std::find(allFamIds.begin(),allFamIds.end(),famArrTmp));
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> famIds=famArr->selectByTupleIdSafe(ids->begin(),ids->end());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> diffFamIds=famIds->getDifferentValues();
  std::vector<int> familyIds;
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > idsPerfamiliyIds;
  int maxVal=getTheMaxAbsFamilyId()+1;
  std::map<std::string,int> families(_families);
  std::map<std::string, std::vector<std::string> > groups(_groups);
  std::vector<std::string> fams;
  bool created(false);
  for(const int *famId=diffFamIds->begin();famId!=diffFamIds->end();famId++)
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ids2Tmp=famIds->getIdsEqual(*famId);
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ids2=ids->selectByTupleId(ids2Tmp->begin(),ids2Tmp->end());
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ids1=famArr->getIdsEqual(*famId);
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret0(ids1->buildSubstractionOptimized(ids2));
      if(ret0->empty())
        {
          bool isFamPresent=false;
          for(std::list< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> >::const_iterator itl=allFamIds.begin();itl!=allFamIds.end() && !isFamPresent;itl++)
            isFamPresent=(*itl)->presenceOfValue(*famId);
          if(!isFamPresent)
            { familyIds.push_back(*famId); idsPerfamiliyIds.push_back(ret0); fams.push_back(FindOrCreateAndGiveFamilyWithId(families,*famId,created)); } // adding *famId in grp
          else
            {
              familyIds.push_back(isNodeGroup?maxVal:-maxVal); idsPerfamiliyIds.push_back(ids2);
              std::string locFamName=FindOrCreateAndGiveFamilyWithId(families,isNodeGroup?maxVal:-maxVal,created);
              fams.push_back(locFamName);
              if(existsFamily(*famId))
                {
                  std::string locFamName2=getFamilyNameGivenId(*famId); std::vector<std::string> v(2); v[0]=locFamName2; v[1]=locFamName;
                  ChangeAllGroupsContainingFamily(groups,getFamilyNameGivenId(*famId),v);
                }
              maxVal++;
            } // modifying all other groups on *famId to lie on maxVal and lie the grp on maxVal
        }
      else
        {
          familyIds.push_back(isNodeGroup?maxVal:-maxVal); idsPerfamiliyIds.push_back(ret0); // modifying all other groups on *famId to lie on maxVal and on maxVal+1
          familyIds.push_back(isNodeGroup?maxVal+1:-maxVal-1); idsPerfamiliyIds.push_back(ids2);//grp lie only on maxVal+1
          std::string n2(FindOrCreateAndGiveFamilyWithId(families,isNodeGroup?maxVal+1:-maxVal-1,created)); fams.push_back(n2);
          if(existsFamily(*famId))
            {
              std::string n1(FindOrCreateAndGiveFamilyWithId(families,isNodeGroup?maxVal:-maxVal,created)); std::vector<std::string> v(2); v[0]=n1; v[1]=n2;
              ChangeAllGroupsContainingFamily(groups,getFamilyNameGivenId(*famId),v);
            }
          maxVal+=2;
        }
    }
  for(std::size_t i=0;i<familyIds.size();i++)
    {
      DataArrayInt *da=idsPerfamiliyIds[i];
      famArr->setPartOfValuesSimple3(familyIds[i],da->begin(),da->end(),0,1,1);
    }
  _families=families;
  _groups=groups;
  _groups[grpName]=fams;
}

void MEDFileMesh::changeAllGroupsContainingFamily(const std::string& familyNameToChange, const std::vector<std::string>& newFamiliesNames)
{
  ChangeAllGroupsContainingFamily(_groups,familyNameToChange,newFamiliesNames);
}

void MEDFileMesh::ChangeAllGroupsContainingFamily(std::map<std::string, std::vector<std::string> >& groups, const std::string& familyNameToChange, const std::vector<std::string>& newFamiliesNames)
{
  std::string fam(familyNameToChange);
  for(std::map<std::string, std::vector<std::string> >::iterator it=groups.begin();it!=groups.end();it++)
    {
      std::vector<std::string>& fams((*it).second);
      std::vector<std::string>::iterator it2=std::find(fams.begin(),fams.end(),fam);
      if(it2!=fams.end())
        {
          fams.erase(it2);
          fams.insert(fams.end(),newFamiliesNames.begin(),newFamiliesNames.end());
        }
    }
}

/*!
 * Returns a name of the family having a given id or, if no such a family exists, creates
 * a new uniquely named family and returns its name.
 *  \param [in] id - the id of the family whose name is required.
 *  \param [out] created - returns \c true if the new family has been created, \c false, else.
 *  \return std::string - the name of the existing or the created family.
 *  \throw If it is not possible to create a unique family name.
 */
std::string MEDFileMesh::findOrCreateAndGiveFamilyWithId(int id, bool& created)
{
  return FindOrCreateAndGiveFamilyWithId(_families,id,created);
}

/*!
 * If it exists a family whose family id is equal to 'id' this method behaves as MEDFileMesh::getFamilyNameGivenId.
 * In this case, 'this' internal states remains unchanged and 'created' out parameter will be set to false.
 * If there is no family whose family id is equal to 'id' a family is created with a name different from those
 * already existing. In this case 'created' will be returned with a value set to true, and internal state
 * will be modified.
 * This method will throws an exception if it is not possible to create a unique family name.
 */
std::string MEDFileMesh::FindOrCreateAndGiveFamilyWithId(std::map<std::string,int>& families, int id, bool& created)
{
  std::vector<std::string> famAlreadyExisting(families.size());
  int ii=0;
  for(std::map<std::string,int>::const_iterator it=families.begin();it!=families.end();it++,ii++)
    {
      if((*it).second!=id)
        {
          famAlreadyExisting[ii]=(*it).first;
        }
      else
        {
          created=false;
          return (*it).first;
        }
    }
  created=true;
  std::ostringstream oss; oss << "Family_" << id;
  std::string ret=CreateNameNotIn(oss.str(),famAlreadyExisting);
  families[ret]=id;
  return ret;
}

/*!
 * Sets names and ids of all families in \a this mesh.
 *  \param [in] info - a map of a family name to a family id.
 */
void MEDFileMesh::setFamilyInfo(const std::map<std::string,int>& info)
{
  _families=info;
}

/*!
 * Sets names of all groups and families constituting them in \a this mesh.
 *  \param [in] info - a map of a group name to a vector of names of families
 *          constituting the group.
 */
void MEDFileMesh::setGroupInfo(const std::map<std::string, std::vector<std::string> >&info)
{
  _groups=info;
}

/*!
 * Returns an id of the family having a given name.
 *  \param [in] name - the name of the family of interest.
 *  \return int - the id of the family of interest.
 *  \throw If no family with such a \a name exists.
 */
int MEDFileMesh::getFamilyId(const std::string& name) const
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

/*!
 * Returns ids of the families having given names.
 *  \param [in] fams - a sequence of the names of families of interest.
 *  \return std::vector<int> - a sequence of the ids of families of interest.
 *  \throw If \a fams contains a name of an inexistent family.
 */
std::vector<int> MEDFileMesh::getFamiliesIds(const std::vector<std::string>& fams) const
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

/*!
 * Returns a maximal abs(id) of families in \a this mesh.
 *  \return int - the maximal norm of family id.
 *  \throw If there are no families in \a this mesh.
 */
int MEDFileMesh::getMaxAbsFamilyId() const
{
  if(_families.empty())
    throw INTERP_KERNEL::Exception("MEDFileMesh::getMaxFamilyId : no families set !");
  int ret=-std::numeric_limits<int>::max();
  for(std::map<std::string,int>::const_iterator it=_families.begin();it!=_families.end();it++)
    {
      ret=std::max(std::abs((*it).second),ret);
    }
  return ret;
}

/*!
 * Returns a maximal id of families in \a this mesh.
 *  \return int - the maximal family id.
 *  \throw If there are no families in \a this mesh.
 */
int MEDFileMesh::getMaxFamilyId() const
{
  if(_families.empty())
    throw INTERP_KERNEL::Exception("MEDFileMesh::getMaxFamilyId : no families set !");
  int ret=-std::numeric_limits<int>::max();
  for(std::map<std::string,int>::const_iterator it=_families.begin();it!=_families.end();it++)
    {
      ret=std::max((*it).second,ret);
    }
  return ret;
}

/*!
 * Returns a minimal id of families in \a this mesh.
 *  \return int - the minimal family id.
 *  \throw If there are no families in \a this mesh.
 */
int MEDFileMesh::getMinFamilyId() const
{
  if(_families.empty())
    throw INTERP_KERNEL::Exception("MEDFileMesh::getMinFamilyId : no families set !");
  int ret=std::numeric_limits<int>::max();
  for(std::map<std::string,int>::const_iterator it=_families.begin();it!=_families.end();it++)
    {
      ret=std::min((*it).second,ret);
    }
  return ret;
}

/*!
 * Returns a maximal id of families in \a this mesh. Not only named families are
 * considered but all family fields as well.
 *  \return int - the maximal family id.
 */
int MEDFileMesh::getTheMaxAbsFamilyId() const
{
  int m1=-std::numeric_limits<int>::max();
  for(std::map<std::string,int>::const_iterator it=_families.begin();it!=_families.end();it++)
    m1=std::max(std::abs((*it).second),m1);
  int m2=getMaxAbsFamilyIdInArrays();
  return std::max(m1,m2);
}

/*!
 * Returns a maximal id of families in \a this mesh. Not only named families are
 * considered but all family fields as well.
 *  \return int - the maximal family id.
 */
int MEDFileMesh::getTheMaxFamilyId() const
{
  int m1=-std::numeric_limits<int>::max();
  for(std::map<std::string,int>::const_iterator it=_families.begin();it!=_families.end();it++)
    m1=std::max((*it).second,m1);
  int m2=getMaxFamilyIdInArrays();
  return std::max(m1,m2);
}

/*!
 * Returns a minimal id of families in \a this mesh. Not only named families are
 * considered but all family fields as well.
 *  \return int - the minimal family id.
 */
int MEDFileMesh::getTheMinFamilyId() const
{
  int m1=std::numeric_limits<int>::max();
  for(std::map<std::string,int>::const_iterator it=_families.begin();it!=_families.end();it++)
    m1=std::min((*it).second,m1);
  int m2=getMinFamilyIdInArrays();
  return std::min(m1,m2);
}

/*!
 * This method only considers the maps. The contain of family array is ignored here.
 * 
 * \sa MEDFileMesh::computeAllFamilyIdsInUse
 */
DataArrayInt *MEDFileMesh::getAllFamiliesIdsReferenced() const
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  std::set<int> v;
  for(std::map<std::string,int>::const_iterator it=_families.begin();it!=_families.end();it++)
    v.insert((*it).second);
  ret->alloc((int)v.size(),1);
  std::copy(v.begin(),v.end(),ret->getPointer());
  return ret.retn();
}

/*!
 * This method does not consider map of family name, family id. Only family field array on different levels is considered.
 * 
 * \sa MEDFileMesh::getAllFamiliesIdsReferenced
 */
DataArrayInt *MEDFileMesh::computeAllFamilyIdsInUse() const
{
  std::vector<int> famLevs=getFamArrNonEmptyLevelsExt();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret;
  for(std::vector<int>::const_iterator it=famLevs.begin();it!=famLevs.end();it++)
    {
      const DataArrayInt *arr=getFamilyFieldAtLevel(*it);//arr not null due to spec of getFamArrNonEmptyLevelsExt
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> dv=arr->getDifferentValues();
      if((DataArrayInt *) ret)
        ret=dv->buildUnion(ret);
      else
        ret=dv;
    }
  return ret.retn();
}

/*!
 * true is returned if no modification has been needed. false if family
 * renumbering has been needed.       
 */
bool MEDFileMesh::ensureDifferentFamIdsPerLevel()
{
  std::vector<int> levs=getNonEmptyLevelsExt();
  std::set<int> allFamIds;
  int maxId=getMaxFamilyId()+1;
  std::map<int,std::vector<int> > famIdsToRenum;
  for(std::vector<int>::const_iterator it=levs.begin();it!=levs.end();it++)
    {
      const DataArrayInt *fam=getFamilyFieldAtLevel(*it);
      if(fam)
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp=fam->getDifferentValues();
          std::set<int> r2;
          std::set_intersection(tmp->begin(),tmp->end(),allFamIds.begin(),allFamIds.end(),std::inserter(r2,r2.end()));
          if(!r2.empty())
            famIdsToRenum[*it].insert(famIdsToRenum[*it].end(),r2.begin(),r2.end());
          std::set<int> r3;
          std::set_union(tmp->begin(),tmp->end(),allFamIds.begin(),allFamIds.end(),std::inserter(r3,r3.end()));
        }
    }
  if(famIdsToRenum.empty())
    return true;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> allIds=getAllFamiliesIdsReferenced();
  for(std::map<int,std::vector<int> >::const_iterator it2=famIdsToRenum.begin();it2!=famIdsToRenum.end();it2++)
    {
      DataArrayInt *fam=const_cast<DataArrayInt *>(getFamilyFieldAtLevel((*it2).first));
      int *famIdsToChange=fam->getPointer();
      std::map<int,int> ren;
      for(std::vector<int>::const_iterator it3=(*it2).second.begin();it3!=(*it2).second.end();it3++,maxId++)
        {
          if(allIds->presenceOfValue(*it3))
            {
              std::string famName=getFamilyNameGivenId(*it3);
              std::vector<std::string> grps=getGroupsOnFamily(famName);
              ren[*it3]=maxId;
              bool dummy;
              std::string newFam=findOrCreateAndGiveFamilyWithId(maxId,dummy);
              for(std::vector<std::string>::const_iterator it4=grps.begin();it4!=grps.end();it4++)
                addFamilyOnGrp((*it4),newFam);
            }
        }
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ids=fam->getIdsEqualList(&(*it2).second[0],&(*it2).second[0]+(*it2).second.size());
      for(const int *id=ids->begin();id!=ids->end();id++)
        famIdsToChange[*id]=ren[famIdsToChange[*id]];
    }
  return false;
}

/*!
 * This method normalizes fam id with the policy explained underneath. This policy is close to those implemented in SMESH.
 * Level #0 famids > 0, Level #-1 famids < 0, Level #-2 famids=0, Level #1 famids=0
 * This policy is those used by SMESH and Trio and that is the opposite of those in MED file.
 * This method will throw an exception if a same family id is detected in different level.
 * \warning This policy is the opposite of those in MED file documentation ...
 */
void MEDFileMesh::normalizeFamIdsTrio()
{
  ensureDifferentFamIdsPerLevel();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> allIds=getAllFamiliesIdsReferenced();
  std::vector<int> levs=getNonEmptyLevelsExt();
  std::set<int> levsS(levs.begin(),levs.end());
  std::set<std::string> famsFetched;
  std::map<std::string,int> families;
  if(std::find(levs.begin(),levs.end(),0)!=levs.end())
    {
      levsS.erase(0);
      const DataArrayInt *fam=getFamilyFieldAtLevel(0);
      if(fam)
        {
          int refId=1;
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp=fam->getDifferentValues();
          std::map<int,int> ren;
          for(const int *it=tmp->begin();it!=tmp->end();it++,refId++)
            ren[*it]=refId;
          int nbOfTuples=fam->getNumberOfTuples();
          int *start=const_cast<DataArrayInt *>(fam)->getPointer();
          for(int *w=start;w!=start+nbOfTuples;w++)
            *w=ren[*w];
          for(const int *it=tmp->begin();it!=tmp->end();it++)
            {
              if(allIds->presenceOfValue(*it))
                {
                  std::string famName=getFamilyNameGivenId(*it);
                  families[famName]=ren[*it];
                  famsFetched.insert(famName);
                }
            }
        }
    }
  if(std::find(levs.begin(),levs.end(),-1)!=levs.end())
    {
      levsS.erase(-1);
      const DataArrayInt *fam=getFamilyFieldAtLevel(-1);
      if(fam)
        {
          int refId=-1;
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp=fam->getDifferentValues();
          std::map<int,int> ren;
          for(const int *it=tmp->begin();it!=tmp->end();it++,refId--)
            ren[*it]=refId;
          int nbOfTuples=fam->getNumberOfTuples();
          int *start=const_cast<DataArrayInt *>(fam)->getPointer();
          for(int *w=start;w!=start+nbOfTuples;w++)
            *w=ren[*w];
          for(const int *it=tmp->begin();it!=tmp->end();it++)
            {
              if(allIds->presenceOfValue(*it))
                {
                  std::string famName=getFamilyNameGivenId(*it);
                  families[famName]=ren[*it];
                  famsFetched.insert(famName);
                }
            }
        }
    }
  for(std::set<int>::const_iterator it2=levsS.begin();it2!=levsS.end();it2++)
    {
      DataArrayInt *fam=const_cast<DataArrayInt*>(getFamilyFieldAtLevel(*it2));
      if(fam)
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp=fam->getDifferentValues();
          fam->fillWithZero();
          for(const int *it3=tmp->begin();it3!=tmp->end();it3++)
            if(allIds->presenceOfValue(*it3))
              {
                std::string famName=getFamilyNameGivenId(*it3);
                families[famName]=0;
                famsFetched.insert(famName);
              }
        }
    }
  //
  std::vector<std::string> allFams=getFamiliesNames();
  std::set<std::string> allFamsS(allFams.begin(),allFams.end());
  std::set<std::string> unFetchedIds;
  std::set_difference(allFamsS.begin(),allFamsS.end(),famsFetched.begin(),famsFetched.end(),std::inserter(unFetchedIds,unFetchedIds.end()));
  for(std::set<std::string>::const_iterator it4=unFetchedIds.begin();it4!=unFetchedIds.end();it4++)
    families[*it4]=_families[*it4];
  _families=families;
}

/*!
 * This method normalizes fam id with the following policy.
 * Level #0 famids < 0, Level #-1 famids < 0 and for Level #1 famids >= 0
 * This policy is those defined in the MED file format but is the opposite of those implemented in SMESH and Trio.
 * This method will throw an exception if a same family id is detected in different level.
 */
void MEDFileMesh::normalizeFamIdsMEDFile()
{
  ensureDifferentFamIdsPerLevel();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> allIds=getAllFamiliesIdsReferenced();
  std::vector<int> levs=getNonEmptyLevelsExt();
  std::set<int> levsS(levs.begin(),levs.end());
  std::set<std::string> famsFetched;
  std::map<std::string,int> families;
  int refId=1;
  if(std::find(levs.begin(),levs.end(),1)!=levs.end())
    {
      levsS.erase(1);
      const DataArrayInt *fam=getFamilyFieldAtLevel(1);
      if(fam)
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp=fam->getDifferentValues();
          std::map<int,int> ren;
          for(const int *it=tmp->begin();it!=tmp->end();it++,refId++)
            ren[*it]=refId;
          int nbOfTuples=fam->getNumberOfTuples();
          int *start=const_cast<DataArrayInt *>(fam)->getPointer();
          for(int *w=start;w!=start+nbOfTuples;w++)
            *w=ren[*w];
          for(const int *it=tmp->begin();it!=tmp->end();it++)
            {
              if(allIds->presenceOfValue(*it))
                {
                  std::string famName=getFamilyNameGivenId(*it);
                  families[famName]=ren[*it];
                  famsFetched.insert(famName);
                }
            }
        }
    }
  refId=-1;
  for(std::set<int>::const_reverse_iterator it2=levsS.rbegin();it2!=levsS.rend();it2++)
    {
      const DataArrayInt *fam=getFamilyFieldAtLevel(*it2);
      if(fam)
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp=fam->getDifferentValues();
          std::map<int,int> ren;
          for(const int *it=tmp->begin();it!=tmp->end();it++,refId--)
            ren[*it]=refId;
          int nbOfTuples=fam->getNumberOfTuples();
          int *start=const_cast<DataArrayInt *>(fam)->getPointer();
          for(int *w=start;w!=start+nbOfTuples;w++)
            *w=ren[*w];
          for(const int *it=tmp->begin();it!=tmp->end();it++)
            {
              if(allIds->presenceOfValue(*it))
                {
                  std::string famName=getFamilyNameGivenId(*it);
                  families[famName]=ren[*it];
                  famsFetched.insert(famName);
                }
            }
        }
    }
  //
  std::vector<std::string> allFams=getFamiliesNames();
  std::set<std::string> allFamsS(allFams.begin(),allFams.end());
  std::set<std::string> unFetchedIds;
  std::set_difference(allFamsS.begin(),allFamsS.end(),famsFetched.begin(),famsFetched.end(),std::inserter(unFetchedIds,unFetchedIds.end()));
  for(std::set<std::string>::const_iterator it4=unFetchedIds.begin();it4!=unFetchedIds.end();it4++)
    families[*it4]=_families[*it4];
  _families=families;
}

/*!
 * Returns a name of the family by its id. If there are several families having the given
 * id, the name first in lexical order is returned.
 *  \param [in] id - the id of the family whose name is required.
 *  \return std::string - the name of the found family.
 *  \throw If no family with the given \a id exists.
 */
std::string MEDFileMesh::getFamilyNameGivenId(int id) const
{
  for(std::map<std::string,int>::const_iterator it=_families.begin();it!=_families.end();it++)
    if((*it).second==id)
      return (*it).first;
  std::ostringstream oss; oss << "MEDFileUMesh::getFamilyNameGivenId : no such family id : " << id;
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

/*!
 * Returns a string describing \a this mesh. This description includes the mesh name and
 * the mesh description string.
 *  \return std::string - the mesh information string.
 */
std::string MEDFileMesh::simpleRepr() const
{
  std::ostringstream oss;
  oss << "(*************************************)\n(* GENERAL INFORMATION ON THE MESH : *)\n(*************************************)\n";
  oss << "- Name of the mesh : <<" << getName() << ">>\n";
  oss << "- Description associated to the mesh : " << getDescription() << std::endl;
  return oss.str();
}

/*!
 * This method is nearly like getFamilyFieldAtLevel method. Except that if the array does not exist at the specified level \a meshDimRelToMaxExt
 * an empty one is created.
 */
DataArrayInt *MEDFileMesh::getOrCreateAndGetFamilyFieldAtLevel(int meshDimRelToMaxExt)
{
  DataArrayInt *ret(getFamilyFieldAtLevel(meshDimRelToMaxExt));
  if(ret)
    return ret;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr(DataArrayInt::New());
  arr->alloc(getSizeAtLevel(meshDimRelToMaxExt),1);
  arr->fillWithZero();
  setFamilyFieldArr(meshDimRelToMaxExt,arr);
  return getFamilyFieldAtLevel(meshDimRelToMaxExt);
}

/*!
 * Returns ids of mesh entities contained in a given group of a given dimension.
 *  \param [in] meshDimRelToMaxExt - a relative dimension of the mesh entities whose ids
 *          are required.
 *  \param [in] grp - the name of the group of interest.
 *  \param [in] renum - if \c true, the optional numbers of entities, if available, are
 *          returned instead of ids. 
 *  \return DataArrayInt * - a new instance of DataArrayInt holding either ids or
 *          numbers, if available and required, of mesh entities of the group. The caller
 *          is to delete this array using decrRef() as it is no more needed. 
 *  \throw If the name of a nonexistent group is specified.
 *  \throw If the family field is missing for \a meshDimRelToMaxExt.
 */
DataArrayInt *MEDFileMesh::getGroupArr(int meshDimRelToMaxExt, const std::string& grp, bool renum) const
{
  std::vector<std::string> tmp(1);
  tmp[0]=grp;
  DataArrayInt *ret=getGroupsArr(meshDimRelToMaxExt,tmp,renum);
  ret->setName(grp);
  return ret;
}

/*!
 * Returns ids of mesh entities contained in given groups of a given dimension.
 *  \param [in] meshDimRelToMaxExt - a relative dimension of the mesh entities whose ids
 *          are required.
 *  \param [in] grps - the names of the groups of interest.
 *  \param [in] renum - if \c true, the optional numbers of entities, if available, are
 *          returned instead of ids.
 *  \return DataArrayInt * - a new instance of DataArrayInt holding either ids or
 *          numbers, if available and required, of mesh entities of the groups. The caller
 *          is to delete this array using decrRef() as it is no more needed. 
 *  \throw If the name of a nonexistent group is present in \a grps.
 *  \throw If the family field is missing for \a meshDimRelToMaxExt.
 */
DataArrayInt *MEDFileMesh::getGroupsArr(int meshDimRelToMaxExt, const std::vector<std::string>& grps, bool renum) const
{
  std::vector<std::string> fams2=getFamiliesOnGroups(grps);
  return getFamiliesArr(meshDimRelToMaxExt,fams2,renum);
}

/*!
 * Returns ids of mesh entities contained in a given family of a given dimension.
 *  \param [in] meshDimRelToMaxExt - a relative dimension of the mesh entities whose ids
 *          are required.
 *  \param [in] fam - the name of the family of interest.
 *  \param [in] renum - if \c true, the optional numbers of entities, if available, are
 *          returned instead of ids. 
 *  \return DataArrayInt * - a new instance of DataArrayInt holding either ids or
 *          numbers, if available and required, of mesh entities of the family. The caller
 *          is to delete this array using decrRef() as it is no more needed. 
 *  \throw If the family field is missing for \a meshDimRelToMaxExt.
 */
DataArrayInt *MEDFileMesh::getFamilyArr(int meshDimRelToMaxExt, const std::string& fam, bool renum) const
{
  std::vector<std::string> tmp(1);
  tmp[0]=fam;
  DataArrayInt *ret=getFamiliesArr(meshDimRelToMaxExt,tmp,renum);
  ret->setName(fam);
  return ret;
}

/*!
 * Returns ids of nodes contained in a given group.
 *  \param [in] grp - the name of the group of interest.
 *  \param [in] renum - if \c true, the optional numbers of nodes, if available, are
 *          returned instead of ids. 
 *  \return DataArrayInt * - a new instance of DataArrayInt holding either ids or
 *          numbers, if available and required, of nodes of the group. The caller
 *          is to delete this array using decrRef() as it is no more needed. 
 *  \throw If the name of a nonexistent group is specified.
 *  \throw If the family field is missing for nodes.
 */
DataArrayInt *MEDFileMesh::getNodeGroupArr(const std::string& grp, bool renum) const
{
  std::vector<std::string> tmp(1);
  tmp[0]=grp;
  DataArrayInt *ret=getNodeGroupsArr(tmp,renum);
  ret->setName(grp);
  return ret;
}

/*!
 * Returns ids of nodes contained in given groups.
 *  \param [in] grps - the names of the groups of interest.
 *  \param [in] renum - if \c true, the optional numbers of nodes, if available, are
 *          returned instead of ids. 
 *  \return DataArrayInt * - a new instance of DataArrayInt holding either ids or
 *          numbers, if available and required, of nodes of the groups. The caller
 *          is to delete this array using decrRef() as it is no more needed. 
 *  \throw If the name of a nonexistent group is present in \a grps.
 *  \throw If the family field is missing for nodes.
 */
DataArrayInt *MEDFileMesh::getNodeGroupsArr(const std::vector<std::string>& grps, bool renum) const
{
  return getGroupsArr(1,grps,renum);
}

/*!
 * Returns ids of nodes contained in a given group.
 *  \param [in] grp - the name of the group of interest.
 *  \param [in] renum - if \c true, the optional numbers of nodes, if available, are
 *          returned instead of ids. 
 *  \return DataArrayInt * - a new instance of DataArrayInt holding either ids or
 *          numbers, if available and required, of nodes of the group. The caller
 *          is to delete this array using decrRef() as it is no more needed. 
 *  \throw If the name of a nonexistent group is specified.
 *  \throw If the family field is missing for nodes.
 */
DataArrayInt *MEDFileMesh::getNodeFamilyArr(const std::string& fam, bool renum) const
{
  std::vector<std::string> tmp(1);
  tmp[0]=fam;
  DataArrayInt *ret=getNodeFamiliesArr(tmp,renum);
  ret->setName(fam);
  return ret;
}

/*!
 * Returns ids of nodes contained in given families.
 *  \param [in] fams - the names of the families of interest.
 *  \param [in] renum - if \c true, the optional numbers of nodes, if available, are
 *          returned instead of ids. 
 *  \return DataArrayInt * - a new instance of DataArrayInt holding either ids or
 *          numbers, if available and required, of nodes of the families. The caller
 *          is to delete this array using decrRef() as it is no more needed. 
 *  \throw If the family field is missing for nodes.
 */
DataArrayInt *MEDFileMesh::getNodeFamiliesArr(const std::vector<std::string>& fams, bool renum) const
{
  return getFamiliesArr(1,fams,renum);
}

/*!
 * Adds groups of given dimension and creates corresponding families and family fields
 * given ids of mesh entities of each group.
 *  \param [in] meshDimRelToMaxExt - the relative mesh dimension of given mesh entities.
 *  \param [in] grps - a sequence of arrays of ids each describing a group.
 *  \param [in] renum - \c true means that \a grps contains not ids but optional numbers
 *          of mesh entities.
 *  \throw If names of some groups in \a grps are equal.
 *  \throw If \a grps includes a group with an empty name.
 *  \throw If \a grps includes invalid ids (or numbers if \a renum == \c true ).
 *  \throw If there are no mesh entities of \a meshDimRelToMaxExt dimension in \a this mesh.
 */
void MEDFileMesh::setGroupsAtLevel(int meshDimRelToMaxExt, const std::vector<const DataArrayInt *>& grps, bool renum)
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
      for(unsigned int ii=0;ii<grps.size();ii++)
        {
          grps2[ii]=MEDFileUMeshSplitL1::Renumber(getRevNumberFieldAtLevel(meshDimRelToMaxExt),grps[ii]);
          grps2[ii]->setName(grps[ii]->getName());
        }
      std::vector<const DataArrayInt *> grps3(grps2.begin(),grps2.end());
      fam=DataArrayInt::MakePartition(grps3,sz,fidsOfGroups);
    }
  int offset=1;
  if(!_families.empty())
    offset=getMaxAbsFamilyId()+1;
  TranslateFamilyIds(meshDimRelToMaxExt==1?offset:-offset,fam,fidsOfGroups);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ids=fam->getDifferentValues();
  appendFamilyEntries(ids,fidsOfGroups,grpsName2);
  setFamilyFieldArr(meshDimRelToMaxExt,fam);
}

/*!
 * This method append into '_families' attribute the families whose ids are in 'famIds'. Warning 'famIds' are expected to be ids
 * not in '_families'. Groups information are given in parameters in order to give to families representative names.
 * For the moment, the two last input parameters are not taken into account.
 */
void MEDFileMesh::appendFamilyEntries(const DataArrayInt *famIds, const std::vector< std::vector<int> >& fidsOfGrps, const std::vector<std::string>& grpNames)
{
  std::map<int,std::string> famInv;
  for(const int *it=famIds->begin();it!=famIds->end();it++)
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

std::vector<INTERP_KERNEL::NormalizedCellType> MEDFileMesh::getAllGeoTypes() const
{
  std::vector<int> levs(getNonEmptyLevels());
  std::vector<INTERP_KERNEL::NormalizedCellType> ret;
  for(std::vector<int>::const_iterator it=levs.begin();it!=levs.end();it++)
    {
      std::vector<INTERP_KERNEL::NormalizedCellType> elts(getGeoTypesAtLevel(*it));
      ret.insert(ret.end(),elts.begin(),elts.end());
    }
  return ret;
}

std::vector<int> MEDFileMesh::getDistributionOfTypes(int meshDimRelToMax) const
{
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> mLev(getGenMeshAtLevel(meshDimRelToMax));
  return mLev->getDistributionOfTypes();
}

void MEDFileMesh::TranslateFamilyIds(int offset, DataArrayInt *famArr, std::vector< std::vector<int> >& famIdsPerGrp)
{
  famArr->applyLin(offset>0?1:-1,offset,0);
  for(std::vector< std::vector<int> >::iterator it1=famIdsPerGrp.begin();it1!=famIdsPerGrp.end();it1++)
    {
      if(offset<0)
        std::transform((*it1).begin(),(*it1).end(),(*it1).begin(),std::negate<int>());
      std::transform((*it1).begin(),(*it1).end(),(*it1).begin(),std::bind2nd(std::plus<int>(),offset));
    }
}

/*!
 * Warning no check is done on 'nameTry' in parameter. It should be non empty.
 * This method returns a name close to 'nameTry' so that it is not already into 'namesToAvoid'.
 * If this method fails to find such a name it will throw an exception.
 */
std::string MEDFileMesh::CreateNameNotIn(const std::string& nameTry, const std::vector<std::string>& namesToAvoid)
{
  //attempt #0
  if(std::find(namesToAvoid.begin(),namesToAvoid.end(),nameTry)==namesToAvoid.end())
    return nameTry;
  //attempt #1
  std::size_t len=nameTry.length();
  for(std::size_t ii=1;ii<len;ii++)
    {
      std::string tmp=nameTry.substr(ii,len-ii);
      if(std::find(namesToAvoid.begin(),namesToAvoid.end(),tmp)==namesToAvoid.end())
        return tmp;
    }
  //attempt #2
  if(len>=1)
    {
      for(std::size_t i=1;i<30;i++)
        {
          std::string tmp1(nameTry.at(0),i);
          tmp1+=nameTry;
          if(std::find(namesToAvoid.begin(),namesToAvoid.end(),tmp1)==namesToAvoid.end())
            return tmp1;
        }
    }
  //attempt #3
  std::string tmp2;
  for(std::vector<std::string>::const_iterator it2=namesToAvoid.begin();it2!=namesToAvoid.end();it2++)
    tmp2+=(*it2);
  if(std::find(namesToAvoid.begin(),namesToAvoid.end(),tmp2)==namesToAvoid.end())
    return tmp2;
  throw INTERP_KERNEL::Exception("MEDFileMesh::CreateNameNotIn : impossible to find a not already used name !");
}

int MEDFileMesh::PutInThirdComponentOfCodeOffset(std::vector<int>& code, int strt)
{
  std::size_t nbOfChunks=code.size()/3;
  if(code.size()%3!=0)
    throw INTERP_KERNEL::Exception("MEDFileMesh::PutInThirdComponentOfCodeOffset : code has invalid size : should be of size 3*x !");
  int ret=strt;
  for(std::size_t i=0;i<nbOfChunks;i++)
    {
      code[3*i+2]=ret;
      ret+=code[3*i+1];
    }
  return ret;
}

/*!
 * This method should be called by any set* method of subclasses to deal automatically with _name attribute.
 * If _name attribute is empty the name of 'm' if taken as _name attribute.
 * If _name is not empty and that 'm' has the same name nothing is done.
 * If _name is not emplt and that 'm' has \b NOT the same name an exception is thrown.
 */
void MEDFileMesh::dealWithTinyInfo(const MEDCouplingMesh *m)
{
  if(!m)
    throw INTERP_KERNEL::Exception("MEDFileMesh::dealWithTinyInfo : input mesh in NULL !");
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
  if(_desc_name.empty())
    _desc_name=m->getDescription();
  else
    {
      std::string name(m->getDescription());
      if(!name.empty())
        {
          if(_desc_name!=name)
            {
              std::ostringstream oss; oss << "MEDFileMesh::dealWithTinyInfo : description of current MEDfile mesh is '" << _desc_name << "' whereas name of input mesh is : '";
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
      std::vector<std::string> grps=getGroupsOnFamily((*it).first);
      std::copy(grps.begin(),grps.end(),std::ostream_iterator<std::string>(oss," "));
      oss << std::endl << std::endl;
    }
}

/*!
 * Returns a new MEDFileUMesh holding the mesh data that has been read from a given MED
 * file. The mesh to load is specified by its name and numbers of a time step and an
 * iteration.
 *  \param [in] fileName - the name of MED file to read.
 *  \param [in] mName - the name of the mesh to read.
 *  \param [in] dt - the number of a time step.
 *  \param [in] it - the number of an iteration.
 *  \return MEDFileUMesh * - a new instance of MEDFileUMesh. The caller is to delete this
 *          mesh using decrRef() as it is no more needed. 
 *  \throw If the file is not readable.
 *  \throw If there is no mesh with given attributes in the file.
 *  \throw If the mesh in the file is not an unstructured one.
 */
MEDFileUMesh *MEDFileUMesh::New(const std::string& fileName, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName.c_str(),MED_ACC_RDONLY);
  return new MEDFileUMesh(fid,mName,dt,it,mrs);
}

/*!
 * Returns a new MEDFileUMesh holding the mesh data that has been read from a given MED
 * file. The first mesh in the file is loaded.
 *  \param [in] fileName - the name of MED file to read.
 *  \return MEDFileUMesh * - a new instance of MEDFileUMesh. The caller is to delete this
 *          mesh using decrRef() as it is no more needed. 
 *  \throw If the file is not readable.
 *  \throw If there is no meshes in the file.
 *  \throw If the mesh in the file is not an unstructured one.
 */
MEDFileUMesh *MEDFileUMesh::New(const std::string& fileName, MEDFileMeshReadSelector *mrs)
{
  std::vector<std::string> ms=MEDLoader::GetMeshNames(fileName);
  if(ms.empty())
    {
      std::ostringstream oss; oss << "MEDFileUMesh::New : no meshes in file \"" << fileName << "\" !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName.c_str(),MED_ACC_RDONLY);
  int dt,it;
  ParaMEDMEM::MEDCouplingMeshType meshType;
  std::string dummy2;
  MEDFileMeshL2::GetMeshIdFromName(fid,ms.front(),meshType,dt,it,dummy2);
  return new MEDFileUMesh(fid,ms.front(),dt,it,mrs);
}

/*!
 * Returns an empty instance of MEDFileUMesh.
 *  \return MEDFileUMesh * - a new instance of MEDFileUMesh. The caller is to delete this
 *          mesh using decrRef() as it is no more needed. 
 */
MEDFileUMesh *MEDFileUMesh::New()
{
  return new MEDFileUMesh;
}

/*!
 * This method loads from file with name \a fileName the mesh called \a mName as New does. The difference is that
 * here only a part of cells contained in the file will be loaded. The selection of cell is specified using the two consecutive parameters
 * \a types and \a slicPerTyp. This method allows to load from a mesh (typically huge) in a MED file a part of cells of that mesh.
 * The part of cells is specified using triplet (start,stop,step) for each geometric type. Only nodes lying on selected cells will be loaded to reduce
 * at most the memory consumtion.
 *
 * \param [in] fileName - the name of the file.
 * \param [in] mName - the name of the mesh to be read.
 * \param [in] types - the list of the geo types of which some part will be taken. A geometric type in \a types must appear only once at most.
 * \param [in] slicPerType - an array of size 3 times larger than \a types that specifies for each type in \a types (in the same order) resp the start, the stop and the step.
 * \param [in] dt - the iteration, that is to say the first element of the pair that locates the asked time step.
 * \param [in] it - the order, that is to say the second element of the pair that locates the asked time step.
 * \param [in] mrs - the request for what to be loaded.
 * \return MEDFileUMesh * - a new instance of MEDFileUMesh. The caller is to delete this mesh using decrRef() as it is no more needed.
 */
MEDFileUMesh *MEDFileUMesh::LoadPartOf(const std::string& fileName, const std::string& mName, const std::vector<INTERP_KERNEL::NormalizedCellType>& types, const std::vector<int>& slicPerTyp, int dt, int it, MEDFileMeshReadSelector *mrs)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid(MEDfileOpen(fileName.c_str(),MED_ACC_RDONLY));
  return MEDFileUMesh::LoadPartOf(fid,mName,types,slicPerTyp,dt,it,mrs);
}

/*!
 * Please refer to the other MEDFileUMesh::LoadPartOf method that has the same semantic and the same parameter (excepted the first).
 * This method is \b NOT wrapped into python.
 */
MEDFileUMesh *MEDFileUMesh::LoadPartOf(med_idt fid, const std::string& mName, const std::vector<INTERP_KERNEL::NormalizedCellType>& types, const std::vector<int>& slicPerTyp, int dt, int it, MEDFileMeshReadSelector *mrs)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileUMesh> ret(MEDFileUMesh::New());
  ret->loadPartUMeshFromFile(fid,mName,types,slicPerTyp,dt,it,mrs);
  return ret.retn();
}

std::size_t MEDFileUMesh::getHeapMemorySizeWithoutChildren() const
{
  std::size_t ret(MEDFileMesh::getHeapMemorySizeWithoutChildren());
  ret+=_ms.capacity()*(sizeof(MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1>));
  return ret;
}

std::vector<const BigMemoryObject *> MEDFileUMesh::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret(MEDFileMesh::getDirectChildrenWithNull());
  ret.push_back((const DataArrayDouble*)_coords);
  ret.push_back((const DataArrayInt *)_fam_coords);
  ret.push_back((const DataArrayInt *)_num_coords);
  ret.push_back((const DataArrayInt *)_rev_num_coords);
  ret.push_back((const DataArrayAsciiChar *)_name_coords);
  ret.push_back((const PartDefinition *)_part_coords);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::const_iterator it=_ms.begin();it!=_ms.end();it++)
    ret.push_back((const MEDFileUMeshSplitL1*) *it);
  return ret;
}

MEDFileMesh *MEDFileUMesh::shallowCpy() const
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileUMesh> ret=new MEDFileUMesh(*this);
  return ret.retn();
}

MEDFileMesh *MEDFileUMesh::createNewEmpty() const
{
  return new MEDFileUMesh;
}

MEDFileMesh *MEDFileUMesh::deepCpy() const
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileUMesh> ret=new MEDFileUMesh(*this);
  if((const DataArrayDouble*)_coords)
    ret->_coords=_coords->deepCpy();
  if((const DataArrayInt*)_fam_coords)
    ret->_fam_coords=_fam_coords->deepCpy();
  if((const DataArrayInt*)_num_coords)
    ret->_num_coords=_num_coords->deepCpy();
  if((const DataArrayInt*)_rev_num_coords)
    ret->_rev_num_coords=_rev_num_coords->deepCpy();
  if((const DataArrayAsciiChar*)_name_coords)
    ret->_name_coords=_name_coords->deepCpy();
  std::size_t i=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::const_iterator it=_ms.begin();it!=_ms.end();it++,i++)
    {
      if((const MEDFileUMeshSplitL1 *)(*it))
        ret->_ms[i]=(*it)->deepCpy(ret->_coords);
    }
  if((const PartDefinition*)_part_coords)
    ret->_part_coords=_part_coords->deepCpy();
  return ret.retn();
}

/*!
 * Checks if \a this and another mesh are equal.
 *  \param [in] other - the mesh to compare with.
 *  \param [in] eps - a precision used to compare real values.
 *  \param [in,out] what - the string returning description of unequal data.
 *  \return bool - \c true if the meshes are equal, \c false, else.
 */
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
  const DataArrayAsciiChar *namec1=_name_coords;
  const DataArrayAsciiChar *namec2=otherC->_name_coords;
  if((namec1==0 && namec2!=0) || (namec1!=0 && namec2==0))
    {
      what="Mismatch of naming arr on nodes ! One is defined and not other !";
      return false;
    }
  if(namec1)
    {
      bool ret=namec1->isEqual(*namec2);
      if(!ret)
        {
          what="Names arr on node differ !";
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
  const PartDefinition *pd0(_part_coords),*pd1(otherC->_part_coords);
  if(!pd0 && !pd1)
    return true;
  if((!pd0 && pd1) || (pd0 && !pd1))
    {
      what=std::string("node part def is defined only for one among this or other !");
      return false;
    }
  return pd0->isEqual(pd1,what);
}

/*!
 * Clears redundant attributes of incorporated data arrays.
 */
void MEDFileUMesh::clearNonDiscrAttributes() const
{
  MEDFileMesh::clearNonDiscrAttributes();
  const DataArrayDouble *coo1=_coords;
  if(coo1)
    (const_cast<DataArrayDouble *>(coo1))->setName("");//This parameter is not discriminant for comparison
  const DataArrayInt *famc1=_fam_coords;
  if(famc1)
    (const_cast<DataArrayInt *>(famc1))->setName("");//This parameter is not discriminant for comparison
  const DataArrayInt *numc1=_num_coords;
  if(numc1)
    (const_cast<DataArrayInt *>(numc1))->setName("");//This parameter is not discriminant for comparison
  const DataArrayAsciiChar *namc1=_name_coords;
  if(namc1)
    (const_cast<DataArrayAsciiChar *>(namc1))->setName("");//This parameter is not discriminant for comparison
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::const_iterator it=_ms.begin();it!=_ms.end();it++)
    {
      const MEDFileUMeshSplitL1 *tmp=(*it);
      if(tmp)
        tmp->clearNonDiscrAttributes();
    }
}

void MEDFileUMesh::setName(const std::string& name)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::iterator it=_ms.begin();it!=_ms.end();it++)
    if((MEDFileUMeshSplitL1 *)(*it)!=0)
      (*it)->setName(name);
  MEDFileMesh::setName(name);
}

MEDFileUMesh::MEDFileUMesh()
{
}

MEDFileUMesh::MEDFileUMesh(med_idt fid, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs)
try
{
    loadUMeshFromFile(fid,mName,dt,it,mrs);
    loadJointsFromFile(fid);
  }
catch(INTERP_KERNEL::Exception& e)
{
    throw e;
}

/*!
 * This method loads only a part of specified cells (given by range of cell ID per geometric type)
 * See MEDFileUMesh::LoadPartOf for detailed description.
 *
 * \sa loadUMeshFromFile
 */
void MEDFileUMesh::loadPartUMeshFromFile(med_idt fid, const std::string& mName, const std::vector<INTERP_KERNEL::NormalizedCellType>& types, const std::vector<int>& slicPerTyp, int dt, int it, MEDFileMeshReadSelector *mrs)
{
  MEDFileUMeshL2 loaderl2;
  ParaMEDMEM::MEDCouplingMeshType meshType;
  int dummy0,dummy1;
  std::string dummy2;
  int mid(MEDFileUMeshL2::GetMeshIdFromName(fid,mName,meshType,dummy0,dummy1,dummy2));
  if(meshType!=UNSTRUCTURED)
    {
      std::ostringstream oss; oss << "loadPartUMeshFromFile : Trying to load as unstructured an existing mesh with name '" << mName << "' !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  loaderl2.loadPart(fid,mid,mName,types,slicPerTyp,dt,it,mrs);
  dispatchLoadedPart(fid,loaderl2,mName,mrs);
}

/*!
 * \brief Write joints in a file
 */
void MEDFileMesh::writeJoints(med_idt fid) const
{
  if ( (const MEDFileJoints*) _joints )
    _joints->write(fid);
}

/*!
 * \brief Load joints in a file or use provided ones
 */
//================================================================================
/*!
 * \brief Load joints in a file or use provided ones
 *  \param [in] fid - MED file descriptor
 *  \param [in] toUseInstedOfReading - optional joints to use instead of reading,
 *          Usually this joints are those just read by another iteration
 *          of namesake mesh, when this method is called by MEDFileMeshMultiTS::New()
 */
//================================================================================

void MEDFileMesh::loadJointsFromFile(med_idt fid, MEDFileJoints* toUseInstedOfReading)
{
  if ( toUseInstedOfReading )
    setJoints( toUseInstedOfReading );
  else
    _joints = MEDFileJoints::New( fid, _name );
}

/*!
 * \brief Return number of joints, which is equal to number of adjacent mesh domains
 */
int MEDFileMesh::getNumberOfJoints()
{
  return ( (MEDFileJoints*) _joints ) ? _joints->getNumberOfJoints() : 0;
}

/*!
 * \brief Return joints with all adjacent mesh domains
 */
MEDFileJoints * MEDFileMesh::getJoints() const
{
  return const_cast<MEDFileJoints*>(& (*_joints));
}

void MEDFileMesh::setJoints( MEDFileJoints* joints )
{
  if ( joints != _joints )
    {
      _joints = joints;
      if ( joints )
        joints->incrRef();
    }
}

/*!
 * This method loads \b all \b the \b mesh \a mName in the file with \a fid descriptor.
 *
 * \sa loadPartUMeshFromFile
 */
void MEDFileUMesh::loadUMeshFromFile(med_idt fid, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs)
{
  MEDFileUMeshL2 loaderl2;
  ParaMEDMEM::MEDCouplingMeshType meshType;
  int dummy0,dummy1;
  std::string dummy2;
  int mid(MEDFileUMeshL2::GetMeshIdFromName(fid,mName,meshType,dummy0,dummy1,dummy2));
  if(meshType!=UNSTRUCTURED)
    {
      std::ostringstream oss; oss << "Trying to load as unstructured an existing mesh with name '" << mName << "' !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  loaderl2.loadAll(fid,mid,mName,dt,it,mrs);
  dispatchLoadedPart(fid,loaderl2,mName,mrs);
}

void MEDFileUMesh::dispatchLoadedPart(med_idt fid, const MEDFileUMeshL2& loaderl2, const std::string& mName, MEDFileMeshReadSelector *mrs)
{
  int lev=loaderl2.getNumberOfLevels();
  _ms.resize(lev);
  for(int i=0;i<lev;i++)
    {
      if(!loaderl2.emptyLev(i))
        _ms[i]=new MEDFileUMeshSplitL1(loaderl2,mName,i);
      else
        _ms[i]=0;
    }
  MEDFileMeshL2::ReadFamiliesAndGrps(fid,mName,_families,_groups,mrs);
  //
  setName(loaderl2.getName());
  setDescription(loaderl2.getDescription());
  setUnivName(loaderl2.getUnivName());
  setIteration(loaderl2.getIteration());
  setOrder(loaderl2.getOrder());
  setTimeValue(loaderl2.getTime());
  setTimeUnit(loaderl2.getTimeUnit());
  _coords=loaderl2.getCoords();
  if(!mrs || mrs->isNodeFamilyFieldReading())
    _fam_coords=loaderl2.getCoordsFamily();
  if(!mrs || mrs->isNodeNumFieldReading())
    _num_coords=loaderl2.getCoordsNum();
  if(!mrs || mrs->isNodeNameFieldReading())
    _name_coords=loaderl2.getCoordsName();
  _part_coords=loaderl2.getPartDefOfCoo();
  computeRevNum();
}

MEDFileUMesh::~MEDFileUMesh()
{
}

void MEDFileUMesh::writeLL(med_idt fid) const
{
  const DataArrayDouble *coo=_coords;
  INTERP_KERNEL::AutoPtr<char> maa=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> desc=MEDLoaderBase::buildEmptyString(MED_COMMENT_SIZE);
  MEDLoaderBase::safeStrCpy(_name.c_str(),MED_NAME_SIZE,maa,_too_long_str);
  MEDLoaderBase::safeStrCpy(_desc_name.c_str(),MED_COMMENT_SIZE,desc,_too_long_str);
  int spaceDim=coo?coo->getNumberOfComponents():0;
  int mdim(0);
  if(!_ms.empty())
    mdim=getMeshDimension();
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
  MEDFILESAFECALLERWR0(MEDmeshCr,(fid,maa,spaceDim,mdim,MED_UNSTRUCTURED_MESH,desc,"",MED_SORT_DTIT,MED_CARTESIAN,comp,unit));
  if(_univ_wr_status)
    MEDFILESAFECALLERWR0(MEDmeshUniversalNameWr,(fid,maa));
  std::string meshName(MEDLoaderBase::buildStringFromFortran(maa,MED_NAME_SIZE));
  MEDFileUMeshL2::WriteCoords(fid,meshName,_iteration,_order,_time,_coords,_fam_coords,_num_coords,_name_coords);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::const_iterator it=_ms.begin();it!=_ms.end();it++)
    if((const MEDFileUMeshSplitL1 *)(*it)!=0)
      (*it)->write(fid,meshName,mdim);
  MEDFileUMeshL2::WriteFamiliesAndGrps(fid,meshName,_families,_groups,_too_long_str);

  writeJoints(fid);
}

/*!
 * Returns relative dimensions of mesh entities (excluding nodes) present in \a this mesh.
 *  \return std::vector<int> - a sequence of the relative dimensions.
 */
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

/*!
 * Returns relative dimensions of mesh entities (including nodes) present in \a this mesh.
 *  \return std::vector<int> - a sequence of the relative dimensions.
 */
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

std::vector<int> MEDFileUMesh::getFamArrNonEmptyLevelsExt() const
{
  std::vector<int> ret;
  const DataArrayInt *famCoo(_fam_coords);
  if(famCoo)
    ret.push_back(1);
  int lev=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::const_iterator it=_ms.begin();it!=_ms.end();it++,lev--)
    {
      const MEDFileUMeshSplitL1 *cur(*it);
      if(cur)
        if(cur->getFamilyField())
          ret.push_back(lev);
    }
  return ret;
}

std::vector<int> MEDFileUMesh::getNumArrNonEmptyLevelsExt() const
{
  std::vector<int> ret;
  const DataArrayInt *numCoo(_num_coords);
  if(numCoo)
    ret.push_back(1);
  int lev=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::const_iterator it=_ms.begin();it!=_ms.end();it++,lev--)
    {
      const MEDFileUMeshSplitL1 *cur(*it);
      if(cur)
        if(cur->getNumberField())
          ret.push_back(lev);
    }
  return ret;
}

std::vector<int> MEDFileUMesh::getNameArrNonEmptyLevelsExt() const
{
  std::vector<int> ret;
  const DataArrayAsciiChar *nameCoo(_name_coords);
  if(nameCoo)
    ret.push_back(1);
  int lev=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::const_iterator it=_ms.begin();it!=_ms.end();it++,lev--)
    {
      const MEDFileUMeshSplitL1 *cur(*it);
      if(cur)
        if(cur->getNameField())
          ret.push_back(lev);
    }
  return ret;
}

/*!
 * Returns all relative mesh levels (**excluding nodes**) where a given group is defined.
 * To include nodes, call getGrpNonEmptyLevelsExt() method.
 *  \param [in] grp - the name of the group of interest.
 *  \return std::vector<int> - a sequence of the relative dimensions.
 */
std::vector<int> MEDFileUMesh::getGrpNonEmptyLevels(const std::string& grp) const
{
  std::vector<std::string> fams=getFamiliesOnGroup(grp);
  return getFamsNonEmptyLevels(fams);
}

/*!
 * Returns all relative mesh levels (including nodes) where a given group is defined.
 *  \param [in] grp - the name of the group of interest.
 *  \return std::vector<int> - a sequence of the relative dimensions.
 */
std::vector<int> MEDFileUMesh::getGrpNonEmptyLevelsExt(const std::string& grp) const
{
  std::vector<std::string> fams=getFamiliesOnGroup(grp);
  return getFamsNonEmptyLevelsExt(fams);
}

/*!
 * Returns all relative mesh levels (**excluding nodes**) where a given family is defined.
 * To include nodes, call getFamNonEmptyLevelsExt() method.
 *  \param [in] fam - the name of the family of interest.
 *  \return std::vector<int> - a sequence of the relative dimensions.
 */
std::vector<int> MEDFileUMesh::getFamNonEmptyLevels(const std::string& fam) const
{
  std::vector<std::string> fams(1,std::string(fam));
  return getFamsNonEmptyLevels(fams);
}

/*!
 * Returns all relative mesh levels (including nodes) where a given family is defined.
 *  \param [in] fam - the name of the family of interest.
 *  \return std::vector<int> - a sequence of the relative dimensions.
 */
std::vector<int> MEDFileUMesh::getFamNonEmptyLevelsExt(const std::string& fam) const
{
  std::vector<std::string> fams(1,std::string(fam));
  return getFamsNonEmptyLevelsExt(fams);
}

/*!
 * Returns all relative mesh levels (**excluding nodes**) where given groups are defined.
 * To include nodes, call getGrpsNonEmptyLevelsExt() method.
 *  \param [in] grps - a sequence of names of the groups of interest.
 *  \return std::vector<int> - a sequence of the relative dimensions.
 */
std::vector<int> MEDFileUMesh::getGrpsNonEmptyLevels(const std::vector<std::string>& grps) const
{
  std::vector<std::string> fams=getFamiliesOnGroups(grps);
  return getFamsNonEmptyLevels(fams);
}

/*!
 * Returns all relative mesh levels (including nodes) where given groups are defined.
 *  \param [in] grps - a sequence of names of the groups of interest.
 *  \return std::vector<int> - a sequence of the relative dimensions.
 */
std::vector<int> MEDFileUMesh::getGrpsNonEmptyLevelsExt(const std::vector<std::string>& grps) const
{
  std::vector<std::string> fams=getFamiliesOnGroups(grps);
  return getFamsNonEmptyLevelsExt(fams);
}

/*!
 * Returns all relative mesh levels (**excluding nodes**) where given families are defined.
 * To include nodes, call getFamsNonEmptyLevelsExt() method.
 *  \param [in] fams - the name of the family of interest.
 *  \return std::vector<int> - a sequence of the relative dimensions.
 */
std::vector<int> MEDFileUMesh::getFamsNonEmptyLevels(const std::vector<std::string>& fams) const
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
 * Returns all relative mesh levels (including nodes) where given families are defined.
 *  \param [in] fams - the names of the families of interest.
 *  \return std::vector<int> - a sequence of the relative dimensions.
 */
std::vector<int> MEDFileUMesh::getFamsNonEmptyLevelsExt(const std::vector<std::string>& fams) const
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

/*!
 * Returns names of groups that partly or fully appear on the level \a meshDimRelToMaxExt.
 *  \param [in] meshDimRelToMaxExt - a relative dimension of interest.
 *  \return std::vector<std::string> - a sequence of group names at \a meshDimRelToMaxExt
 *          level. 
 */
std::vector<std::string> MEDFileUMesh::getGroupsOnSpecifiedLev(int meshDimRelToMaxExt) const
{
  std::vector<std::string> ret;
  std::vector<std::string> allGrps=getGroupsNames();
  for(std::vector<std::string>::const_iterator it=allGrps.begin();it!=allGrps.end();it++)
    {
      std::vector<int> levs=getGrpNonEmptyLevelsExt((*it));
      if(std::find(levs.begin(),levs.end(),meshDimRelToMaxExt)!=levs.end())
        ret.push_back(*it);
    }
  return ret;
}

int MEDFileUMesh::getMaxAbsFamilyIdInArrays() const
{
  int ret=-std::numeric_limits<int>::max(),tmp=-1;
  if((const DataArrayInt *)_fam_coords)
    {
      int val=_fam_coords->getMaxValue(tmp);
      ret=std::max(ret,std::abs(val));
    }
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::const_iterator it=_ms.begin();it!=_ms.end();it++)
    {
      if((const MEDFileUMeshSplitL1 *)(*it))
        {
          const DataArrayInt *da=(*it)->getFamilyField();
          if(da)
            {
              int val=da->getMaxValue(tmp);
              ret=std::max(ret,std::abs(val));
            }
        }
    }
  return ret;
}

int MEDFileUMesh::getMaxFamilyIdInArrays() const
{
  int ret=-std::numeric_limits<int>::max(),tmp=-1;
  if((const DataArrayInt *)_fam_coords)
    {
      int val=_fam_coords->getMaxValue(tmp);
      ret=std::max(ret,val);
    }
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::const_iterator it=_ms.begin();it!=_ms.end();it++)
    {
      if((const MEDFileUMeshSplitL1 *)(*it))
        {
          const DataArrayInt *da=(*it)->getFamilyField();
          if(da)
            {
              int val=da->getMaxValue(tmp);
              ret=std::max(ret,val);
            }
        }
    }
  return ret;
}

int MEDFileUMesh::getMinFamilyIdInArrays() const
{
  int ret=std::numeric_limits<int>::max(),tmp=-1;
  if((const DataArrayInt *)_fam_coords)
    {
      int val=_fam_coords->getMinValue(tmp);
      ret=std::min(ret,val);
    }
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::const_iterator it=_ms.begin();it!=_ms.end();it++)
    {
      if((const MEDFileUMeshSplitL1 *)(*it))
        {
          const DataArrayInt *da=(*it)->getFamilyField();
          if(da)
            {
              int val=da->getMinValue(tmp);
              ret=std::min(ret,val);
            }
        }
    }
  return ret;
}

/*!
 * Returns the dimension on cells in \a this mesh.
 *  \return int - the mesh dimension.
 *  \throw If there are no cells in this mesh.
 */
int MEDFileUMesh::getMeshDimension() const
{
  int lev=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::const_iterator it=_ms.begin();it!=_ms.end();it++,lev++)
    if((const MEDFileUMeshSplitL1 *)(*it)!=0)
      return (*it)->getMeshDimension()+lev;
  throw INTERP_KERNEL::Exception("MEDFileUMesh::getMeshDimension : impossible to find a mesh dimension !");
}

/*!
 * Returns the space dimension of \a this mesh that is equal to number of components in
 * the node coordinates array.
 *  \return int - the space dimension of \a this mesh.
 *  \throw If the node coordinates array is not available.
 */
int MEDFileUMesh::getSpaceDimension() const
{
  const DataArrayDouble *coo=_coords;
  if(!coo)
    throw INTERP_KERNEL::Exception(" MEDFileUMesh::getSpaceDimension : no coords set !");
  return coo->getNumberOfComponents();
}

/*!
 * Returns a string describing \a this mesh.
 *  \return std::string - the mesh information string.
 */
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

/*!
 * Returns a full textual description of \a this mesh.
 *  \return std::string - the string holding the mesh description.
 */
std::string MEDFileUMesh::advancedRepr() const
{
  return simpleRepr();
}

/*!
 * Returns number of mesh entities of a given relative dimension in \a this mesh.
 *  \param [in] meshDimRelToMaxExt - the relative dimension of interest.
 *  \return int - the number of entities.
 *  \throw If no mesh entities of dimension \a meshDimRelToMaxExt are available in \a this mesh.
 */
int MEDFileUMesh::getSizeAtLevel(int meshDimRelToMaxExt) const
{
  if(meshDimRelToMaxExt==1)
    {
      if(!((const DataArrayDouble *)_coords))
        throw INTERP_KERNEL::Exception("MEDFileUMesh::getSizeAtLevel : no coordinates specified !");
      return _coords->getNumberOfTuples();
    }
  return getMeshAtLevSafe(meshDimRelToMaxExt)->getSize();
}

/*!
 * Returns the family field for mesh entities of a given dimension.
 *  \param [in] meshDimRelToMaxExt - the relative dimension of mesh entities.
 *  \return const DataArrayInt * - the family field. It is an array of ids of families
 *          each mesh entity belongs to. It can be \c NULL.
 */
const DataArrayInt *MEDFileUMesh::getFamilyFieldAtLevel(int meshDimRelToMaxExt) const
{
  if(meshDimRelToMaxExt==1)
    return _fam_coords;
  const MEDFileUMeshSplitL1 *l1(getMeshAtLevSafe(meshDimRelToMaxExt));
  return l1->getFamilyField();
}

DataArrayInt *MEDFileUMesh::getFamilyFieldAtLevel(int meshDimRelToMaxExt)
{
  if(meshDimRelToMaxExt==1)
    return _fam_coords;
  MEDFileUMeshSplitL1 *l1(getMeshAtLevSafe(meshDimRelToMaxExt));
  return l1->getFamilyField();
}

/*!
 * Returns the optional numbers of mesh entities of a given dimension.
 *  \param [in] meshDimRelToMaxExt - the relative dimension of mesh entities.
 *  \return const DataArrayInt * - the array of the entity numbers.
 *  \throw If there are no mesh entities of \a meshDimRelToMaxExt dimension in \a this mesh.
 */
const DataArrayInt *MEDFileUMesh::getNumberFieldAtLevel(int meshDimRelToMaxExt) const
{
  if(meshDimRelToMaxExt==1)
    return _num_coords;
  const MEDFileUMeshSplitL1 *l1=getMeshAtLevSafe(meshDimRelToMaxExt);
  return l1->getNumberField();
}

const DataArrayAsciiChar *MEDFileUMesh::getNameFieldAtLevel(int meshDimRelToMaxExt) const
{
  if(meshDimRelToMaxExt==1)
    return _name_coords;
  const MEDFileUMeshSplitL1 *l1=getMeshAtLevSafe(meshDimRelToMaxExt);
  return l1->getNameField();
}

/*!
 * This method returns for a specified relative level \a meshDimRelToMaxExt the part effectively read (if the instance is the result of the read of a file).
 *
 * \param [in] meshDimRelToMaxExt - the extended relative level for which the part definition is requested.
 * \param [in] gt - The input geometric type for which the part definition is requested.
 * \return the part definition owned by \a this. So no need to deallocate the returned instance.
 */
const PartDefinition *MEDFileUMesh::getPartDefAtLevel(int meshDimRelToMaxExt, INTERP_KERNEL::NormalizedCellType gt) const
{
  if(meshDimRelToMaxExt==1)
    return _part_coords;
  const MEDFileUMeshSplitL1 *l1(getMeshAtLevSafe(meshDimRelToMaxExt));
  return l1->getPartDef(gt);
}

int MEDFileUMesh::getNumberOfNodes() const
{
  const DataArrayDouble *coo(_coords);
  if(!coo)
    throw INTERP_KERNEL::Exception(" MEDFileUMesh::getNumberOfNodes : no coords set !");
  return coo->getNumberOfTuples();
}

int MEDFileUMesh::getNumberOfCellsAtLevel(int meshDimRelToMaxExt) const
{
  const MEDFileUMeshSplitL1 *l1(getMeshAtLevSafe(meshDimRelToMaxExt));
  return l1->getNumberOfCells();
}

bool MEDFileUMesh::hasImplicitPart() const
{
  return false;
}

int MEDFileUMesh::buildImplicitPartIfAny(INTERP_KERNEL::NormalizedCellType gt) const
{
  throw INTERP_KERNEL::Exception("MEDFileUMesh::buildImplicitPartIfAny : unstructured meshes do not have implicit part !");
}

void MEDFileUMesh::releaseImplicitPartIfAny() const
{
}

void MEDFileUMesh::whichAreNodesFetched(const MEDFileField1TSStructItem& st, const MEDFileFieldGlobsReal *globs, std::vector<bool>& nodesFetched) const
{
  std::size_t sz(st.getNumberOfItems());
  for(std::size_t i=0;i<sz;i++)
    {
      INTERP_KERNEL::NormalizedCellType curGt(st[i].getGeo());
      const MEDCoupling1GTUMesh *m(getDirectUndergroundSingleGeoTypeMesh(curGt));
      if(st[i].getPflName().empty())
        m->computeNodeIdsAlg(nodesFetched);
      else
        {
          const DataArrayInt *arr(globs->getProfile(st[i].getPflName()));
          MEDCouplingAutoRefCountObjectPtr<MEDCoupling1GTUMesh> m2(dynamic_cast<MEDCoupling1GTUMesh *>(m->buildPartOfMySelf(arr->begin(),arr->end(),true)));
          m2->computeNodeIdsAlg(nodesFetched);
        }
    }
}

/*!
 * Returns the optional numbers of mesh entities of a given dimension transformed using
 * DataArrayInt::invertArrayN2O2O2N().
 *  \param [in] meshDimRelToMaxExt - the relative dimension of mesh entities.
 *  \return const DataArrayInt * - the array of the entity numbers transformed using
 *          DataArrayInt::invertArrayN2O2O2N().
 *  \throw If there are no mesh entities of \a meshDimRelToMaxExt dimension in \a this mesh.
 */
const DataArrayInt *MEDFileUMesh::getRevNumberFieldAtLevel(int meshDimRelToMaxExt) const
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

/*!
 * Returns a pointer to the node coordinates array of \a this mesh \b without
 * incrementing its reference counter, thus there is no need to decrRef() it by the caller.
 */
DataArrayDouble *MEDFileUMesh::getCoords() const
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> tmp(_coords);
  if((DataArrayDouble *)tmp)
    {
      return tmp;
    }
  return 0;
}

/*!
 * Returns a new MEDCouplingUMesh corresponding to mesh entities included in a given
 * group of \a this mesh. Only mesh entities of a given dimension are included in the
 * new mesh.
 *  \param [in] meshDimRelToMaxExt - the relative dimension of mesh entities of interest.
 *  \param [in] grp - the name of the group whose mesh entities are included in the
 *          new mesh.
 *  \param [in] renum - if \c true, cells and nodes of the result mesh are permuted
 *          according to the optional numbers of entities, if available.
 *  \return MEDCouplingUMesh * - a new instance of MEDCouplingUMesh. The caller is to
 *          delete this mesh using decrRef() as it is no more needed.
 *  \throw If the name of a nonexistent group is specified.
 *  \throw If there are no mesh entities of \a meshDimRelToMaxExt dimension in \a this mesh.
 */
MEDCouplingUMesh *MEDFileUMesh::getGroup(int meshDimRelToMaxExt, const std::string& grp, bool renum) const
{
  synchronizeTinyInfoOnLeaves();
  std::vector<std::string> tmp(1);
  tmp[0]=grp;
  return getGroups(meshDimRelToMaxExt,tmp,renum);
}

/*!
 * Returns a new MEDCouplingUMesh corresponding to mesh entities included in given
 * groups of \a this mesh. Only mesh entities of a given dimension are included in the
 * new mesh.
 *  \param [in] meshDimRelToMaxExt - the relative dimension of mesh entities of interest.
 *  \param [in] grps - a sequence of group names whose mesh entities are included in the
 *          new mesh.
 *  \param [in] renum - if \c true, cells and nodes of the result mesh are permuted
 *          according to the optional numbers of entities, if available.
 *  \return MEDCouplingUMesh * - a new instance of MEDCouplingUMesh. The caller is to
 *          delete this mesh using decrRef() as it is no more needed.
 *  \throw If a name of a nonexistent group is present in \a grps.
 *  \throw If there are no mesh entities of \a meshDimRelToMaxExt dimension in \a this mesh.
 */
MEDCouplingUMesh *MEDFileUMesh::getGroups(int meshDimRelToMaxExt, const std::vector<std::string>& grps, bool renum) const
{
  synchronizeTinyInfoOnLeaves();
  std::vector<std::string> fams2=getFamiliesOnGroups(grps);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> zeRet=getFamilies(meshDimRelToMaxExt,fams2,renum);
  if(grps.size()==1 && ((MEDCouplingUMesh *)zeRet))
    zeRet->setName(grps[0]);
  return zeRet.retn();
}

/*!
 * Returns a new MEDCouplingUMesh corresponding to mesh entities included in a given
 * family of \a this mesh. Only mesh entities of a given dimension are included in the
 * new mesh.
 *  \param [in] meshDimRelToMaxExt - the relative dimension of mesh entities of interest.
 *  \param [in] fam - the name of the family whose mesh entities are included in the
 *          new mesh.
 *  \param [in] renum - if \c true, cells and nodes of the result mesh are permuted
 *          according to the optional numbers of entities, if available.
 *  \return MEDCouplingUMesh * - a new instance of MEDCouplingUMesh. The caller is to
 *          delete this mesh using decrRef() as it is no more needed.
 *  \throw If a name of a nonexistent family is present in \a grps.
 *  \throw If there are no mesh entities of \a meshDimRelToMaxExt dimension in \a this mesh.
 */
MEDCouplingUMesh *MEDFileUMesh::getFamily(int meshDimRelToMaxExt, const std::string& fam, bool renum) const
{
  synchronizeTinyInfoOnLeaves();
  std::vector<std::string> tmp(1);
  tmp[0]=fam;
  return getFamilies(meshDimRelToMaxExt,tmp,renum);
}

/*!
 * Returns a new MEDCouplingUMesh corresponding to mesh entities included in given
 * families of \a this mesh. Only mesh entities of a given dimension are included in the
 * new mesh.
 *  \param [in] meshDimRelToMaxExt - the relative dimension of mesh entities of interest.
 *  \param [in] fams - a sequence of family names whose mesh entities are included in the
 *          new mesh.
 *  \param [in] renum - if \c true, cells and nodes of the result mesh are permuted
 *          according to the optional numbers of entities, if available.
 *  \return MEDCouplingUMesh * - a new instance of MEDCouplingUMesh. The caller is to
 *          delete this mesh using decrRef() as it is no more needed.
 *  \throw If a name of a nonexistent family is present in \a fams.
 *  \throw If there are no mesh entities of \a meshDimRelToMaxExt dimension in \a this mesh.
 */
MEDCouplingUMesh *MEDFileUMesh::getFamilies(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum) const
{
  synchronizeTinyInfoOnLeaves();
  if(meshDimRelToMaxExt==1)
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr=getFamiliesArr(1,fams,renum);
      MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> ret=MEDCouplingUMesh::New();
      MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> c=_coords->selectByTupleId(arr->getConstPointer(),arr->getConstPointer()+arr->getNbOfElems());
      ret->setCoords(c);
      return ret.retn();
    }
  std::vector<int> famIds=getFamiliesIds(fams);
  const MEDFileUMeshSplitL1 *l1=getMeshAtLevSafe(meshDimRelToMaxExt);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> zeRet;
  if(!famIds.empty())
    zeRet=l1->getFamilyPart(&famIds[0],&famIds[0]+famIds.size(),renum);
  else
    zeRet=l1->getFamilyPart(0,0,renum);
  if(fams.size()==1 && ((MEDCouplingUMesh *)zeRet))
    zeRet->setName(fams[0]);
  return zeRet.retn();
}

/*!
 * Returns ids of mesh entities contained in given families of a given dimension.
 *  \param [in] meshDimRelToMaxExt - a relative dimension of the mesh entities whose ids
 *          are required.
 *  \param [in] fams - the names of the families of interest.
 *  \param [in] renum - if \c true, the optional numbers of entities, if available, are
 *          returned instead of ids.
 *  \return DataArrayInt * - a new instance of DataArrayInt holding either ids or
 *          numbers, if available and required, of mesh entities of the families. The caller
 *          is to delete this array using decrRef() as it is no more needed. 
 *  \throw If the family field is missing for \a meshDimRelToMaxExt.
 */
DataArrayInt *MEDFileUMesh::getFamiliesArr(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum) const
{
  std::vector<int> famIds=getFamiliesIds(fams);
  if(meshDimRelToMaxExt==1)
    {
      if((const DataArrayInt *)_fam_coords)
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da;
          if(!famIds.empty())
            da=_fam_coords->getIdsEqualList(&famIds[0],&famIds[0]+famIds.size());
          else
            da=_fam_coords->getIdsEqualList(0,0);
          if(renum)
            return MEDFileUMeshSplitL1::Renumber(_num_coords,da);
          else
            return da.retn();
        }
      else
        throw INTERP_KERNEL::Exception("MEDFileUMesh::getFamiliesArr : no family array specified on nodes !");
    }
  const MEDFileUMeshSplitL1 *l1=getMeshAtLevSafe(meshDimRelToMaxExt);
  if(!famIds.empty())
    return l1->getFamilyPartArr(&famIds[0],&famIds[0]+famIds.size(),renum);
  else
    return l1->getFamilyPartArr(0,0,renum);
}

/*!
 * Returns a MEDCouplingUMesh of a given relative dimension.
 * \warning If \a meshDimRelToMaxExt == 1 (which means nodes), the returned mesh **is not
 * valid**. This is a feature, because MEDLoader does not create cells that do not exist! 
 * To build a valid MEDCouplingUMesh from the returned one in this case,
 * call MEDCouplingUMesh::Build0DMeshFromCoords().
 *  \param [in] meshDimRelToMax - the relative dimension of interest.
 *  \param [in] renum - if \c true, the returned mesh is permuted according to the
 *          optional numbers of mesh entities.
 *  \return MEDCouplingUMesh * - a pointer to MEDCouplingUMesh that the caller is to
 *          delete using decrRef() as it is no more needed. 
 *  \throw If there are no mesh entities of \a meshDimRelToMaxExt dimension in \a this mesh.
 * \sa getGenMeshAtLevel()
 */
MEDCouplingUMesh *MEDFileUMesh::getMeshAtLevel(int meshDimRelToMaxExt, bool renum) const
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
          umesh->setName(getName());
          return umesh;
        }
    }
  const MEDFileUMeshSplitL1 *l1=getMeshAtLevSafe(meshDimRelToMaxExt);
  return l1->getWholeMesh(renum);
}

/*!
 * Returns a MEDCouplingUMesh of a given relative dimension.
 * \warning If \a meshDimRelToMaxExt == 1 (which means nodes), the returned mesh **is not
 * valid**. This is a feature, because MEDLoader does not create cells that do not exist! 
 * To build a valid MEDCouplingUMesh from the returned one in this case,
 * call MEDCouplingUMesh::Build0DMeshFromCoords().
 *  \param [in] meshDimRelToMax - the relative dimension of interest.
 *  \param [in] renum - if \c true, the returned mesh is permuted according to the
 *          optional numbers of mesh entities.
 *  \return MEDCouplingMesh * - a pointer to MEDCouplingUMesh that the caller is to
 *          delete using decrRef() as it is no more needed. 
 *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in \a this mesh.
 * \sa getMeshAtLevel()
 */
MEDCouplingMesh *MEDFileUMesh::getGenMeshAtLevel(int meshDimRelToMax, bool renum) const
{
  return getMeshAtLevel(meshDimRelToMax,renum);
}

std::vector<int> MEDFileUMesh::getDistributionOfTypes(int meshDimRelToMax) const
{
  const MEDFileUMeshSplitL1 *l1(getMeshAtLevSafe(meshDimRelToMax));
  return l1->getDistributionOfTypes();
}

/*!
 * Returns a MEDCouplingUMesh of a relative dimension == 0.
 *  \param [in] renum - if \c true, the returned mesh is permuted according to the
 *          optional numbers of mesh entities.
 *  \return MEDCouplingUMesh * - a pointer to MEDCouplingUMesh that the caller is to
 *          delete using decrRef() as it is no more needed. 
 *  \throw If there are no mesh entities of the relative dimension == 0 in \a this mesh.
 */
MEDCouplingUMesh *MEDFileUMesh::getLevel0Mesh(bool renum) const
{
  return getMeshAtLevel(0,renum);
}

/*!
 * Returns a MEDCouplingUMesh of a relative dimension == -1.
 *  \param [in] renum - if \c true, the returned mesh is permuted according to the
 *          optional numbers of mesh entities.
 *  \return MEDCouplingUMesh * - a pointer to MEDCouplingUMesh that the caller is to
 *          delete using decrRef() as it is no more needed. 
 *  \throw If there are no mesh entities of the relative dimension == -1 in \a this mesh.
 */
MEDCouplingUMesh *MEDFileUMesh::getLevelM1Mesh(bool renum) const
{
  return getMeshAtLevel(-1,renum);
}

/*!
 * Returns a MEDCouplingUMesh of a relative dimension == -2.
 *  \param [in] renum - if \c true, the returned mesh is permuted according to the
 *          optional numbers of mesh entities.
 *  \return MEDCouplingUMesh * - a pointer to MEDCouplingUMesh that the caller is to
 *          delete using decrRef() as it is no more needed. 
 *  \throw If there are no mesh entities of the relative dimension == -2 in \a this mesh.
 */
MEDCouplingUMesh *MEDFileUMesh::getLevelM2Mesh(bool renum) const
{
  return getMeshAtLevel(-2,renum);
}

/*!
 * Returns a MEDCouplingUMesh of a relative dimension == -3.
 *  \param [in] renum - if \c true, the returned mesh is permuted according to the
 *          optional numbers of mesh entities.
 *  \return MEDCouplingUMesh * - a pointer to MEDCouplingUMesh that the caller is to
 *          delete using decrRef() as it is no more needed. 
 *  \throw If there are no mesh entities of the relative dimension == -3 in \a this mesh.
 */
MEDCouplingUMesh *MEDFileUMesh::getLevelM3Mesh(bool renum) const
{
  return getMeshAtLevel(-3,renum);
}

/*!
 * This method is for advanced users. There is two storing strategy of mesh in \a this.
 * Either MEDCouplingUMesh, or vector of MEDCoupling1GTUMesh instances.
 * When assignement is done the first one is done, which is not optimal in write mode for MED file.
 * This method allows to switch from MEDCouplingUMesh mode to MEDCoupling1GTUMesh mode.
 */
void MEDFileUMesh::forceComputationOfParts() const
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::const_iterator it=_ms.begin();it!=_ms.end();it++)
    {
      const MEDFileUMeshSplitL1 *elt(*it);
      if(elt)
        elt->forceComputationOfParts();
    }
}

/*!
 * This method returns a vector of mesh parts containing each exactly one geometric type.
 * This method will never launch an automatic computation of split by type (an INTERP_KERNEL::Exception will be then thrown).
 * This method is only for memory aware users.
 * The returned pointers are **NOT** new object pointer. No need to mange them.
 */
std::vector<MEDCoupling1GTUMesh *> MEDFileUMesh::getDirectUndergroundSingleGeoTypeMeshes(int meshDimRelToMax) const
{
  const MEDFileUMeshSplitL1 *sp(getMeshAtLevSafe(meshDimRelToMax));
  return sp->getDirectUndergroundSingleGeoTypeMeshes();
}

/*!
 * This method returns the part of \a this having the geometric type \a gt.
 * If such part is not existing an exception will be thrown.
 * The returned pointer is **NOT** new object pointer. No need to mange it.
 */
MEDCoupling1GTUMesh *MEDFileUMesh::getDirectUndergroundSingleGeoTypeMesh(INTERP_KERNEL::NormalizedCellType gt) const
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(gt);
  int lev=(int)cm.getDimension()-getMeshDimension();
  const MEDFileUMeshSplitL1 *sp(getMeshAtLevSafe(lev));
  return sp->getDirectUndergroundSingleGeoTypeMesh(gt);
}

/*!
 * Given a relative level \a meshDimRelToMax it returns the sorted vector of geometric types present in \a this.
 * \throw if the reqsuested \a meshDimRelToMax does not exist.
 */
std::vector<INTERP_KERNEL::NormalizedCellType> MEDFileUMesh::getGeoTypesAtLevel(int meshDimRelToMax) const
{
  const MEDFileUMeshSplitL1 *sp(getMeshAtLevSafe(meshDimRelToMax));
  return sp->getGeoTypes();
}

/*!
 * This method extracts from whole family field ids the part relative to the input parameter \a gt.
 * \param [in] gt - the geometric type for which the family field is asked.
 * \return DataArrayInt * - a pointer to DataArrayInt that the caller is to
 *          delete using decrRef() as it is no more needed.
 * \sa MEDFileUMesh::extractNumberFieldOnGeoType
 */
DataArrayInt *MEDFileUMesh::extractFamilyFieldOnGeoType(INTERP_KERNEL::NormalizedCellType gt) const
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(gt);
  int lev=(int)cm.getDimension()-getMeshDimension();
  const MEDFileUMeshSplitL1 *sp(getMeshAtLevSafe(lev));
  return sp->extractFamilyFieldOnGeoType(gt);
}

/*!
 * This method extracts from whole number field ids the part relative to the input parameter \a gt.
 * \param [in] gt - the geometric type for which the number field is asked.
 * \return DataArrayInt * - a pointer to DataArrayInt that the caller is to
 *          delete using decrRef() as it is no more needed.
 * \sa MEDFileUMesh::extractFamilyFieldOnGeoType
 */
DataArrayInt *MEDFileUMesh::extractNumberFieldOnGeoType(INTERP_KERNEL::NormalizedCellType gt) const
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(gt);
  int lev=(int)cm.getDimension()-getMeshDimension();
  const MEDFileUMeshSplitL1 *sp(getMeshAtLevSafe(lev));
  return sp->extractNumberFieldOnGeoType(gt);
}

/*!
 * This method returns for specified geometric type \a gt the relative level to \a this.
 * If the relative level is empty an exception will be thrown.
 */
int MEDFileUMesh::getRelativeLevOnGeoType(INTERP_KERNEL::NormalizedCellType gt) const
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(gt);
  int ret((int)cm.getDimension()-getMeshDimension());
  getMeshAtLevSafe(ret);//To test that returned value corresponds to a valid level.
  return ret;
}

const MEDFileUMeshSplitL1 *MEDFileUMesh::getMeshAtLevSafe(int meshDimRelToMaxExt) const
{
  if(meshDimRelToMaxExt==1)
    throw INTERP_KERNEL::Exception("Dimension request is invalid : asking for node level (1) !");
  if(meshDimRelToMaxExt>1)
    throw INTERP_KERNEL::Exception("Dimension request is invalid (>1) !");
  int tracucedRk=-meshDimRelToMaxExt;
  if(tracucedRk>=(int)_ms.size())
    throw INTERP_KERNEL::Exception("Invalid mesh dim relative to max given ! Too low !");
  if((const MEDFileUMeshSplitL1 *)_ms[tracucedRk]==0)
    throw INTERP_KERNEL::Exception("On specified lev (or entity) no cells exists !");
  return _ms[tracucedRk];
}

MEDFileUMeshSplitL1 *MEDFileUMesh::getMeshAtLevSafe(int meshDimRelToMaxExt)
{
  if(meshDimRelToMaxExt==1)
    throw INTERP_KERNEL::Exception("Dimension request is invalid : asking for node level (1) !");
  if(meshDimRelToMaxExt>1)
    throw INTERP_KERNEL::Exception("Dimension request is invalid (>1) !");
  int tracucedRk=-meshDimRelToMaxExt;
  if(tracucedRk>=(int)_ms.size())
    throw INTERP_KERNEL::Exception("Invalid mesh dim relative to max given ! Too low !");
  if((const MEDFileUMeshSplitL1 *)_ms[tracucedRk]==0)
    throw INTERP_KERNEL::Exception("On specified lev (or entity) no cells exists !");
  return _ms[tracucedRk];
}

void MEDFileUMesh::checkMeshDimCoherency(int meshDim, int meshDimRelToMax) const
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

/*!
 * Sets the node coordinates array of \a this mesh.
 *  \param [in] coords - the new node coordinates array.
 *  \throw If \a coords == \c NULL.
 */
void MEDFileUMesh::setCoords(DataArrayDouble *coords)
{
  if(!coords)
    throw INTERP_KERNEL::Exception("MEDFileUMesh::setCoords : null pointer in input !");
  coords->checkAllocated();
  int nbOfTuples=coords->getNumberOfTuples();
  _coords=coords;
  coords->incrRef();
  _fam_coords=DataArrayInt::New();
  _fam_coords->alloc(nbOfTuples,1);
  _fam_coords->fillWithZero();
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::iterator it=_ms.begin();it!=_ms.end();it++)
    if((MEDFileUMeshSplitL1 *)(*it))
      (*it)->setCoords(coords);
}

/*!
 * Removes all groups of a given dimension in \a this mesh.
 *  \param [in] meshDimRelToMaxExt - the relative dimension of interest.
 *  \throw If there are no mesh entities of \a meshDimRelToMaxExt dimension in \a this mesh.
 */
void MEDFileUMesh::eraseGroupsAtLevel(int meshDimRelToMaxExt)
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

/*!
 * Removes all families with ids not present in the family fields of \a this mesh.
 */
void MEDFileUMesh::optimizeFamilies()
{
  std::vector<int> levs=getNonEmptyLevelsExt();
  std::set<int> allFamsIds;
  for(std::vector<int>::const_iterator it=levs.begin();it!=levs.end();it++)
    {
      const DataArrayInt *ffield=getFamilyFieldAtLevel(*it);
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ids=ffield->getDifferentValues();
      std::set<int> res;
      std::set_union(ids->begin(),ids->end(),allFamsIds.begin(),allFamsIds.end(),std::inserter(res,res.begin()));
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

/**
 * \b this must be filled at level 0 and -1, typically the -1 level being (part of) the descending connectivity
 * of the top level. This method build a "crack", or an inner boundary, in \b this along the group of level -1 named grpNameM1.
 * The boundary is built according to the following method:
 *  - all nodes along the boundary which are not lying on an internal extremity of the (-1)-level group are duplicated (so the
 * coordinates array is extended).
 *  - new (-1)-level cells are built lying on those new nodes. So the edges/faces along the group are duplicated.
 *  After this operation a top-level cell bordering the group will loose some neighbors (typically the cell which is  on the
 *  other side of the group is no more a neighbor)
 *   - finally, the connectivity of (part of) the top level-cells bordering the group is also modified so that some cells
 *  bordering the newly created boundary use the newly computed nodes.
 *
 *  \param[in] grpNameM1 name of the (-1)-level group defining the boundary
 *  \param[out] nodesDuplicated ids of the initial nodes which have been duplicated (and whose copy is put at the end of
 *  the coord array)
 *  \param[out] cellsModified ids of the cells whose connectivity has been modified (to use the newly created nodes)
 *  \param[out] cellsNotModified ids of the rest of cells bordering the new boundary whose connectivity remains unchanged.
 */
void MEDFileUMesh::buildInnerBoundaryAlongM1Group(const std::string& grpNameM1, DataArrayInt *&nodesDuplicated,
                                           DataArrayInt *&cellsModified, DataArrayInt *&cellsNotModified)
{
  std::vector<int> levs=getNonEmptyLevels();
  if(std::find(levs.begin(),levs.end(),0)==levs.end() || std::find(levs.begin(),levs.end(),-1)==levs.end())
    throw INTERP_KERNEL::Exception("MEDFileUMesh::buildInnerBoundaryAlongM1Group : This method works only for mesh definied on level 0 and -1 !");
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m0=getMeshAtLevel(0);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m1=getMeshAtLevel(-1);
  int nbNodes=m0->getNumberOfNodes();
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m11=getGroup(-1,grpNameM1);
  DataArrayInt *tmp00=0,*tmp11=0,*tmp22=0;
  m0->findNodesToDuplicate(*m11,tmp00,tmp11,tmp22);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> nodeIdsToDuplicate(tmp00);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> cellsToModifyConn0(tmp11);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> cellsToModifyConn1(tmp22);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> tmp0=static_cast<MEDCouplingUMesh *>(m0->buildPartOfMySelf(cellsToModifyConn0->begin(),cellsToModifyConn0->end(),true));
  // node renumbering of cells in m1 impacted by duplication of node but not in group 'grpNameM1' on level -1
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> descTmp0=DataArrayInt::New(),descITmp0=DataArrayInt::New(),revDescTmp0=DataArrayInt::New(),revDescITmp0=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> tmp0Desc=tmp0->buildDescendingConnectivity(descTmp0,descITmp0,revDescTmp0,revDescITmp0);
  descTmp0=0; descITmp0=0; revDescTmp0=0; revDescITmp0=0;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> cellsInM1ToRenumW2=tmp0Desc->getCellIdsLyingOnNodes(nodeIdsToDuplicate->begin(),nodeIdsToDuplicate->end(),false);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> cellsInM1ToRenumW3=static_cast<MEDCouplingUMesh *>(tmp0Desc->buildPartOfMySelf(cellsInM1ToRenumW2->begin(),cellsInM1ToRenumW2->end(),true));
  DataArrayInt *cellsInM1ToRenumW4Tmp=0;
  m1->areCellsIncludedIn(cellsInM1ToRenumW3,2,cellsInM1ToRenumW4Tmp);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> cellsInM1ToRenumW4(cellsInM1ToRenumW4Tmp);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> cellsInM1ToRenumW5=cellsInM1ToRenumW4->getIdsInRange(0,m1->getNumberOfCells());
  cellsInM1ToRenumW5->transformWithIndArr(cellsInM1ToRenumW4->begin(),cellsInM1ToRenumW4->end());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> grpIds=getGroupArr(-1,grpNameM1);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> cellsInM1ToRenum=cellsInM1ToRenumW5->buildSubstraction(grpIds);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m1Part=static_cast<MEDCouplingUMesh *>(m1->buildPartOfMySelf(cellsInM1ToRenum->begin(),cellsInM1ToRenum->end(),true));
  m1Part->duplicateNodesInConn(nodeIdsToDuplicate->begin(),nodeIdsToDuplicate->end(),nbNodes);
  m1->setPartOfMySelf(cellsInM1ToRenum->begin(),cellsInM1ToRenum->end(),*m1Part);
  // end of node renumbering of cells in m1 impacted by duplication of node but not in group of level -1 'grpNameM1'
  tmp0->duplicateNodes(nodeIdsToDuplicate->begin(),nodeIdsToDuplicate->end());
  m0->setCoords(tmp0->getCoords());
  m0->setPartOfMySelf(cellsToModifyConn0->begin(),cellsToModifyConn0->end(),*tmp0);
  m1->setCoords(m0->getCoords());
  _coords=m0->getCoords(); _coords->incrRef();
  // duplication of cells in group 'grpNameM1' on level -1
  m11->duplicateNodesInConn(nodeIdsToDuplicate->begin(),nodeIdsToDuplicate->end(),nbNodes); m11->setCoords(m0->getCoords());
  std::vector<const MEDCouplingUMesh *> v(2); v[0]=m1; v[1]=m11;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> newm1=MEDCouplingUMesh::AggregateSortedByTypeMeshesOnSameCoords(v,tmp00,tmp11);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> szOfCellGrpOfSameType(tmp00);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> idInMsOfCellGrpOfSameType(tmp11);
  //
  newm1->setName(getName());
  const DataArrayInt *fam=getFamilyFieldAtLevel(-1);
  if(!fam)
    throw INTERP_KERNEL::Exception("MEDFileUMesh::buildInnerBoundaryAlongM1Group : internal problem !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> newFam=DataArrayInt::New();
  newFam->alloc(newm1->getNumberOfCells(),1);
  // Get a new family ID: care must be taken if we need a positive ID or a negative one:
  // Positive ID for family of nodes, negative for all the rest.
  int idd;
  if (m1->getMeshDimension() == 0)
    idd=getMaxFamilyId()+1;
  else
    idd=getMinFamilyId()-1;
  int globStart=0,start=0,end,globEnd;
  int nbOfChunks=szOfCellGrpOfSameType->getNumberOfTuples();
  for(int i=0;i<nbOfChunks;i++)
    {
      globEnd=globStart+szOfCellGrpOfSameType->getIJ(i,0);
      if(idInMsOfCellGrpOfSameType->getIJ(i,0)==0)
        {
          end=start+szOfCellGrpOfSameType->getIJ(i,0);
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> part=fam->selectByTupleId2(start,end,1);
          newFam->setPartOfValues1(part,globStart,globEnd,1,0,1,1,true);
          start=end;
        }
      else
        {
          newFam->setPartOfValuesSimple1(idd,globStart,globEnd,1,0,1,1);
        }
      globStart=globEnd;
    }
  newm1->setCoords(getCoords());
  setMeshAtLevel(-1,newm1);
  setFamilyFieldArr(-1,newFam);
  std::string grpName2(grpNameM1); grpName2+="_dup";
  addFamily(grpName2,idd);
  addFamilyOnGrp(grpName2,grpName2);
  //
  fam=_fam_coords;
  if(fam)
    {
      int newNbOfNodes=getCoords()->getNumberOfTuples();
      newFam=DataArrayInt::New(); newFam->alloc(newNbOfNodes,1);
      newFam->setPartOfValues1(fam,0,nbNodes,1,0,1,1,true);
      newFam->setPartOfValuesSimple1(0,nbNodes,newNbOfNodes,1,0,1,1);
      _fam_coords=newFam;
    }
  nodesDuplicated=nodeIdsToDuplicate.retn();
  cellsModified=cellsToModifyConn0.retn();
  cellsNotModified=cellsToModifyConn1.retn();
}

/*!
 * \param [out] oldCode retrieves the distribution of types before the call if true is returned
 * \param [out] newCode etrieves the distribution of types after the call if true is returned
 * \param [out] o2nRenumCell tells for **all levels** the old 2 new renumbering of cells.
 * 
 * \return false if no modification has been performed linked to the unpolyzation. Neither cell type, not cell numbers. When false is returned no need of field on cells or on gauss renumbering.
 * Inversely, if true is returned, it means that distribution of cell by geometric type has changed and field on cell and field on gauss point must be renumbered.
 */
bool MEDFileUMesh::unPolyze(std::vector<int>& oldCode, std::vector<int>& newCode, DataArrayInt *& o2nRenumCell)
{
  o2nRenumCell=0; oldCode.clear(); newCode.clear();
  std::vector<int> levs=getNonEmptyLevels();
  bool ret=false;
  std::vector< const DataArrayInt* > renumCellsSplited;//same than memorySaverIfThrow
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > memorySaverIfThrow;//same than renumCellsSplited only in case of throw
  int start=0;
  int end=0;
  for(std::vector<int>::reverse_iterator it=levs.rbegin();it!=levs.rend();it++)
    {
      MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m=getMeshAtLevel(*it);
      std::vector<int> code1=m->getDistributionOfTypes();
      end=PutInThirdComponentOfCodeOffset(code1,start);
      oldCode.insert(oldCode.end(),code1.begin(),code1.end());
      bool hasChanged=m->unPolyze();
      DataArrayInt *fake=0;
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> o2nCellsPart=m->getLevArrPerCellTypes(MEDCouplingUMesh::MEDMEM_ORDER,
          MEDCouplingUMesh::MEDMEM_ORDER+MEDCouplingUMesh::N_MEDMEM_ORDER,fake);
      fake->decrRef();
      renumCellsSplited.push_back(o2nCellsPart); memorySaverIfThrow.push_back(o2nCellsPart);
      if(hasChanged)
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> o2nCellsPart2=o2nCellsPart->buildPermArrPerLevel();
          m->renumberCells(o2nCellsPart2->getConstPointer(),false);
          ret=true;
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> famField2,numField2;
          const DataArrayInt *famField=getFamilyFieldAtLevel(*it); if(famField) { famField->incrRef(); famField2=const_cast<DataArrayInt *>(famField); }
          const DataArrayInt *numField=getNumberFieldAtLevel(*it); if(numField) { numField->incrRef(); numField2=const_cast<DataArrayInt *>(numField); }
          setMeshAtLevel(*it,m);
          std::vector<int> code2=m->getDistributionOfTypes();
          end=PutInThirdComponentOfCodeOffset(code2,start);
          newCode.insert(newCode.end(),code2.begin(),code2.end());
          //
          if(o2nCellsPart2->isIdentity())
            continue;
          if(famField)
            {
              MEDCouplingAutoRefCountObjectPtr<DataArrayInt> newFamField=famField->renumber(o2nCellsPart2->getConstPointer());
              setFamilyFieldArr(*it,newFamField);
            }
          if(numField)
            {
              MEDCouplingAutoRefCountObjectPtr<DataArrayInt> newNumField=numField->renumber(o2nCellsPart2->getConstPointer());
              setRenumFieldArr(*it,newNumField);
            }
        }
      else
        {
          newCode.insert(newCode.end(),code1.begin(),code1.end());
        }
      start=end;
    }
  if(ret)
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> renumCells=DataArrayInt::Aggregate(renumCellsSplited);
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> o2nRenumCellRet=renumCells->buildPermArrPerLevel();
      o2nRenumCell=o2nRenumCellRet.retn();
    }
  return ret;
}

/*! \cond HIDDEN_ITEMS */
struct MEDLoaderAccVisit1
{
  MEDLoaderAccVisit1():_new_nb_of_nodes(0) { }
  int operator()(bool val) { return val?_new_nb_of_nodes++:-1; }
  int _new_nb_of_nodes;
};
/*! \endcond */

/*!
 * Array returned is the correspondance in \b old \b to \b new format. The returned array is newly created and should be dealt by the caller.
 * The maximum value stored in returned array is the number of nodes of \a this minus 1 after call of this method.
 * The size of returned array is the number of nodes of the old (previous to the call of this method) number of nodes.
 * -1 values in returned array means that the corresponding old node is no more used.
 *
 * \return newly allocated array containing correspondance in \b old \b to \b new format. If all nodes in \a this are fetched \c NULL pointer is returned and nothing
 *         is modified in \a this.
 * \throw If no coordinates are set in \a this or if there is in any available mesh in \a this a cell having a nodal connectivity containing a node id not in the range of
 *  set coordinates.
 */
DataArrayInt *MEDFileUMesh::zipCoords()
{
  const DataArrayDouble *coo(getCoords());
  if(!coo)
    throw INTERP_KERNEL::Exception("MEDFileUMesh::zipCoords : no coordinates set in this !");
  int nbOfNodes(coo->getNumberOfTuples());
  std::vector<bool> nodeIdsInUse(nbOfNodes,false);
  std::vector<int> neLevs(getNonEmptyLevels());
  for(std::vector<int>::const_iterator lev=neLevs.begin();lev!=neLevs.end();lev++)
    {
      const MEDFileUMeshSplitL1 *zeLev(getMeshAtLevSafe(*lev));
      if(zeLev->isMeshStoredSplitByType())
        {
          std::vector<MEDCoupling1GTUMesh *> ms(zeLev->getDirectUndergroundSingleGeoTypeMeshes());
          for(std::vector<MEDCoupling1GTUMesh *>::const_iterator it=ms.begin();it!=ms.end();it++)
            if(*it)
              (*it)->computeNodeIdsAlg(nodeIdsInUse);
        }
      else
        {
          MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> mesh(zeLev->getWholeMesh(false));
          mesh->computeNodeIdsAlg(nodeIdsInUse);
        }
    }
  int nbrOfNodesInUse((int)std::count(nodeIdsInUse.begin(),nodeIdsInUse.end(),true));
  if(nbrOfNodesInUse==nbOfNodes)
    return 0;//no need to update _part_coords
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret(DataArrayInt::New()); ret->alloc(nbOfNodes,1);
  std::transform(nodeIdsInUse.begin(),nodeIdsInUse.end(),ret->getPointer(),MEDLoaderAccVisit1());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret2(ret->invertArrayO2N2N2OBis(nbrOfNodesInUse));
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> newCoords(coo->selectByTupleIdSafe(ret2->begin(),ret2->end()));
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> newFamCoords;
  MEDCouplingAutoRefCountObjectPtr<DataArrayAsciiChar> newNameCoords;
  if((const DataArrayInt *)_fam_coords)
    newFamCoords=_fam_coords->selectByTupleIdSafe(ret2->begin(),ret2->end());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> newNumCoords;
  if((const DataArrayInt *)_num_coords)
    newNumCoords=_num_coords->selectByTupleIdSafe(ret2->begin(),ret2->end());
  if((const DataArrayAsciiChar *)_name_coords)
    newNameCoords=static_cast<DataArrayAsciiChar *>(_name_coords->selectByTupleIdSafe(ret2->begin(),ret2->end()));
  _coords=newCoords; _fam_coords=newFamCoords; _num_coords=newNumCoords; _name_coords=newNameCoords; _rev_num_coords=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::iterator it=_ms.begin();it!=_ms.end();it++)
    {
      if((MEDFileUMeshSplitL1*)*it)
        {
          (*it)->renumberNodesInConn(ret->begin());
          (*it)->setCoords(_coords);
        }
    }
  // updates _part_coords
  const PartDefinition *pc(_part_coords);
  if(pc)
    {
      MEDCouplingAutoRefCountObjectPtr<PartDefinition> tmpPD(DataArrayPartDefinition::New(ret2));
      _part_coords=tmpPD->composeWith(pc);
    }
  return ret.retn();
}

/*!
 * This method performs an extrusion along a path defined by \a m1D.
 * \a this is expected to be a mesh with max mesh dimension equal to 2.
 * \a m1D is expected to be a mesh with space dimesion equal to 3 and mesh dimension equal to 1.
 * Mesh dimensions of returned mesh is incremented by one compared to thoose in \a this.
 * This method scans all levels in \a this
 * and put them in the returned mesh. All groups in \a this are also put in the returned mesh.
 *
 * \param [in] m1D - the mesh defining the extrusion path.
 * \param [in] policy - defines the policy of extrusion (see MEDCouplingUMesh::buildExtrudedMesh for more details)
 * \return - a new reference on mesh (you have to deal with using decrRef). The returned mesh will have the same name than \a this.
 *
 * \sa MEDCouplingUMesh::buildExtrudedMesh
 */
MEDFileUMesh *MEDFileUMesh::buildExtrudedMesh(const MEDCouplingUMesh *m1D, int policy) const
{
  if(getMeshDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDFileUMesh::buildExtrudedMesh : this is expected to be with mesh dimension equal to 2 !");
  MEDCouplingAutoRefCountObjectPtr<MEDFileUMesh> ret(MEDFileUMesh::New());
  m1D->checkCoherency();
  if(m1D->getMeshDimension()!=1)
    throw INTERP_KERNEL::Exception("MEDFileUMesh::buildExtrudedMesh : input mesh must have a mesh dimension equal to one !");
  int nbRep(m1D->getNumberOfCells());
  std::vector<int> levs(getNonEmptyLevels());
  std::vector<std::string> grps(getGroupsNames());
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> > zeList;
  DataArrayDouble *coords(0);
  std::size_t nbOfLevsOut(levs.size()+1);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > o2ns(nbOfLevsOut);
  for(std::vector<int>::const_iterator lev=levs.begin();lev!=levs.end();lev++)
    {
      MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> item(getMeshAtLevel(*lev));
      item=item->clone(false);
      item->changeSpaceDimension(3+(*lev),0.);//no problem non const but change DataArrayDouble for coordinates do not alter data
      MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> tmp(static_cast<MEDCouplingUMesh *>(m1D->deepCpy()));
      tmp->changeSpaceDimension(3+(*lev),0.);
      MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> elt(item->buildExtrudedMesh(tmp,policy));
      zeList.push_back(elt);
      if(*lev==0)
        coords=elt->getCoords();
    }
  if(!coords)
    throw INTERP_KERNEL::Exception("MEDFileUMesh::buildExtrudedMesh : internal error !");
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> >::iterator it=zeList.begin();it!=zeList.end();it++)
    {
      (*it)->setName(getName());
      (*it)->setCoords(coords);
    }
  for(std::size_t ii=0;ii!=zeList.size();ii++)
    {
      int lev(levs[ii]);
      MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> elt(zeList[ii]);
      if(lev<=-1)
        {
          MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> elt1(getMeshAtLevel(lev+1));
          MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> elt2(elt1->clone(false));
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp(elt2->getNodalConnectivity()->deepCpy());
          elt2->setConnectivity(tmp,elt2->getNodalConnectivityIndex());
          elt2->shiftNodeNumbersInConn(nbRep*elt1->getNumberOfNodes());
          elt1->setCoords(elt->getCoords()); elt2->setCoords(elt->getCoords());
          std::vector<const MEDCouplingUMesh *> elts(3);
          elts[0]=elt; elts[1]=elt1; elts[2]=elt2;
          elt=MEDCouplingUMesh::MergeUMeshesOnSameCoords(elts);
          elt->setName(getName());
        }
      //
      o2ns[ii]=elt->sortCellsInMEDFileFrmt();
      ret->setMeshAtLevel(lev,elt);
    }
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> endLev(getMeshAtLevel(levs.back())),endLev2;
  endLev=endLev->clone(false); endLev->setCoords(coords);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp(endLev->getNodalConnectivity()->deepCpy());
  endLev2=endLev->clone(false); endLev2->setConnectivity(tmp,endLev->getNodalConnectivityIndex());
  endLev2->shiftNodeNumbersInConn(nbRep*getNumberOfNodes());
  endLev=MEDCouplingUMesh::MergeUMeshesOnSameCoords(endLev,endLev2);
  o2ns[levs.size()]=endLev->sortCellsInMEDFileFrmt();
  endLev->setName(getName());
  ret->setMeshAtLevel(levs.back()-1,endLev);
  //
  for(std::size_t ii=0;ii!=zeList.size();ii++)
    {
      int lev(levs[ii]);
      std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > outGrps;
      std::vector< const DataArrayInt * > outGrps2;
      if(lev<=-1)
        {
          for(std::vector<std::string>::const_iterator grp=grps.begin();grp!=grps.end();grp++)
            {
              MEDCouplingAutoRefCountObjectPtr<DataArrayInt> grpArr(getGroupArr(lev+1,*grp));
              if(!grpArr->empty())
                {
                  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> grpArr1(grpArr->deepCpy()),grpArr2(grpArr->deepCpy());
                  int offset0(zeList[ii]->getNumberOfCells());
                  int offset1(offset0+getNumberOfCellsAtLevel(lev+1));
                  grpArr1->applyLin(1,offset0); grpArr2->applyLin(1,offset1);
                  std::ostringstream oss; oss << grpArr2->getName() << "_top";
                  grpArr2->setName(oss.str());
                  grpArr1->transformWithIndArr(o2ns[ii]->begin(),o2ns[ii]->end());
                  grpArr2->transformWithIndArr(o2ns[ii]->begin(),o2ns[ii]->end());
                  outGrps.push_back(grpArr1); outGrps.push_back(grpArr2);
                  outGrps2.push_back(grpArr1); outGrps2.push_back(grpArr2);
                }
            }
        }
      //
      for(std::vector<std::string>::const_iterator grp=grps.begin();grp!=grps.end();grp++)
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> grpArr(getGroupArr(lev,*grp));
          if(!grpArr->empty())
            {
              int nbCellsB4Extrusion(getNumberOfCellsAtLevel(lev));
              std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > grpArrs(nbRep);
              std::vector< const DataArrayInt *> grpArrs2(nbRep);
              for(int iii=0;iii<nbRep;iii++)
                {
                  grpArrs[iii]=grpArr->deepCpy(); grpArrs[iii]->applyLin(1,iii*nbCellsB4Extrusion);
                  grpArrs2[iii]=grpArrs[iii];
                }
              MEDCouplingAutoRefCountObjectPtr<DataArrayInt> grpArrExt(DataArrayInt::Aggregate(grpArrs2));
              grpArrExt->transformWithIndArr(o2ns[ii]->begin(),o2ns[ii]->end());
              std::ostringstream grpName; grpName << *grp << "_extruded";
              grpArrExt->setName(grpName.str());
              outGrps.push_back(grpArrExt);
              outGrps2.push_back(grpArrExt);
            }
        }
      ret->setGroupsAtLevel(lev,outGrps2);
    }
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > outGrps;
  std::vector< const DataArrayInt * > outGrps2;
  for(std::vector<std::string>::const_iterator grp=grps.begin();grp!=grps.end();grp++)
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> grpArr1(getGroupArr(levs.back(),*grp));
      if(grpArr1->empty())
        continue;
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> grpArr2(grpArr1->deepCpy());
      std::ostringstream grpName; grpName << *grp << "_top";
      grpArr2->setName(grpName.str());
      grpArr2->applyLin(1,getNumberOfCellsAtLevel(levs.back()));
      outGrps.push_back(grpArr1); outGrps.push_back(grpArr2);
      outGrps2.push_back(grpArr1); outGrps2.push_back(grpArr2);
    }
  ret->setGroupsAtLevel(levs.back()-1,outGrps2);
  return ret.retn();
}

/*!
 * This method converts all linear cells in \a this into quadratic cells (following the \a conversionType policy).
 * All the cells converted are put in the returned instance. This method applies all the groups and families in \a this to returned instance.
 * Groups on nodes and families on nodes are copied directly to the returned instance without transformation.
 *
 * \param [in] conversionType - conversionType specifies the type of conversion expected. Only 0 (default) and 1 are supported presently. 0 those that creates the 'most' simple
 *             corresponding quadratic cells. 1 is those creating the 'most' complex.
 * \param [in] eps - detection threshold for coordinates.
 * \return A new instance that is the result of the conversion. The caller has the ownership of this returned instance.
 *
 * \sa MEDCouplingUMesh::convertLinearCellsToQuadratic , quadraticToLinear
 */
MEDFileUMesh *MEDFileUMesh::linearToQuadratic(int conversionType, double eps) const
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileUMesh> ret(MEDFileUMesh::New());
  int initialNbNodes(getNumberOfNodes());
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m0Tmp(getMeshAtLevel(0));
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m0(dynamic_cast<MEDCouplingUMesh *>(m0Tmp->deepCpy()));
  {
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> notUsed(m0->convertLinearCellsToQuadratic(conversionType));
  }
  DataArrayDouble *zeCoords(m0->getCoords());
  ret->setMeshAtLevel(0,m0);
  std::vector<int> levs(getNonEmptyLevels());
  const DataArrayInt *famField(getFamilyFieldAtLevel(0));
  if(famField)
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> famFieldCpy(famField->deepCpy());
      ret->setFamilyFieldArr(0,famFieldCpy);
    }
  famField=getFamilyFieldAtLevel(1);
  if(famField)
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> fam(DataArrayInt::New()); fam->alloc(zeCoords->getNumberOfTuples(),1);
      fam->fillWithZero();
      fam->setPartOfValues1(famField,0,initialNbNodes,1,0,1,1);
      ret->setFamilyFieldArr(1,fam);
    }
  ret->copyFamGrpMapsFrom(*this);
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> partZeCoords(zeCoords->selectByTupleId2(initialNbNodes,zeCoords->getNumberOfTuples(),1));
  for(std::vector<int>::const_iterator lev=levs.begin();lev!=levs.end();lev++)
    {
      if(*lev==0)
        continue;
      MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m1Tmp(getMeshAtLevel(*lev));
      MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m1(dynamic_cast<MEDCouplingUMesh *>(m1Tmp->deepCpy()));
      if(m1->getMeshDimension()!=0)
        {
          {
            MEDCouplingAutoRefCountObjectPtr<DataArrayInt> notUsed(m1->convertLinearCellsToQuadratic(conversionType));
          }//kill unused notUsed var
          MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> m1Coords(m1->getCoords()->selectByTupleId2(initialNbNodes,m1->getNumberOfNodes(),1));
          DataArrayInt *b(0);
          bool a(partZeCoords->areIncludedInMe(m1Coords,eps,b));
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> bSafe(b);
          if(!a)
            {
              std::ostringstream oss; oss << "MEDFileUMesh::linearCellsToQuadratic : for level " << *lev << " problem to identify nodes generated !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
          b->applyLin(1,initialNbNodes);
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> l0(DataArrayInt::New()); l0->alloc(initialNbNodes,1); l0->iota();
          std::vector<const DataArrayInt *> v(2); v[0]=l0; v[1]=b;
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> renum(DataArrayInt::Aggregate(v));
          m1->renumberNodesInConn(renum->begin());
        }
      m1->setCoords(zeCoords);
      ret->setMeshAtLevel(*lev,m1);
      famField=getFamilyFieldAtLevel(*lev);
      if(famField)
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> famFieldCpy(famField->deepCpy());
          ret->setFamilyFieldArr(*lev,famFieldCpy);
        }
    }
  return ret.retn();
}

/*!
 * This method converts all quadratic cells in \a this into linear cells.
 * All the cells converted are put in the returned instance. This method applies all the groups and families in \a this to returned instance.
 * Groups on nodes and families on nodes are copied directly to the returned instance without transformation.
 *
 * \param [in] eps - detection threshold for coordinates.
 * \return A new instance that is the result of the conversion. The caller has the ownership of this returned instance.
 *
 * \sa MEDCouplingUMesh::convertLinearCellsToQuadratic , linearToQuadratic
 */
MEDFileUMesh *MEDFileUMesh::quadraticToLinear(double eps) const
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileUMesh> ret(MEDFileUMesh::New());
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m0Tmp(getMeshAtLevel(0));
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m0(dynamic_cast<MEDCouplingUMesh *>(m0Tmp->deepCpy()));
  m0->convertQuadraticCellsToLinear();
  m0->zipCoords();
  DataArrayDouble *zeCoords(m0->getCoords());
  ret->setMeshAtLevel(0,m0);
  std::vector<int> levs(getNonEmptyLevels());
  const DataArrayInt *famField(getFamilyFieldAtLevel(0));
  if(famField)
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> famFieldCpy(famField->deepCpy());
      ret->setFamilyFieldArr(0,famFieldCpy);
    }
  famField=getFamilyFieldAtLevel(1);
  if(famField)
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> fam(famField->selectByTupleId2(0,zeCoords->getNumberOfTuples(),1));
      ret->setFamilyFieldArr(1,fam);
    }
  ret->copyFamGrpMapsFrom(*this);
  for(std::vector<int>::const_iterator lev=levs.begin();lev!=levs.end();lev++)
    {
      if(*lev==0)
        continue;
      MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m1Tmp(getMeshAtLevel(*lev));
      MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m1(dynamic_cast<MEDCouplingUMesh *>(m1Tmp->deepCpy()));
      m1->convertQuadraticCellsToLinear();
      m1->zipCoords();
      DataArrayInt *b(0);
      bool a(zeCoords->areIncludedInMe(m1->getCoords(),eps,b));
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> bSafe(b);
      if(!a)
        {
          std::ostringstream oss; oss << "MEDFileUMesh::quadraticToLinear : for level " << *lev << " problem to identify nodes generated !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      m1->renumberNodesInConn(b->begin());
      m1->setCoords(zeCoords);
      ret->setMeshAtLevel(*lev,m1);
      famField=getFamilyFieldAtLevel(*lev);
      if(famField)
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> famFieldCpy(famField->deepCpy());
          ret->setFamilyFieldArr(*lev,famFieldCpy);
        }
    }
  return ret.retn();
}

void MEDFileUMesh::serialize(std::vector<double>& tinyDouble, std::vector<int>& tinyInt, std::vector<std::string>& tinyStr, std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> >& bigArraysI, MEDCouplingAutoRefCountObjectPtr<DataArrayDouble>& bigArrayD)
{
  clearNonDiscrAttributes();
  forceComputationOfParts();
  tinyDouble.clear(); tinyInt.clear(); tinyStr.clear(); bigArraysI.clear(); bigArrayD=0;
  std::vector<int> layer0;
  layer0.push_back(_order); //0 i
  layer0.push_back(_iteration);//1 i
  layer0.push_back(getSpaceDimension());//2 i
  tinyDouble.push_back(_time);//0 d
  tinyStr.push_back(_name);//0 s
  tinyStr.push_back(_desc_name);//1 s
  for(int i=0;i<getSpaceDimension();i++)
    tinyStr.push_back(_coords->getInfoOnComponent(i));
  layer0.push_back((int)_families.size());//3 i <- key info aa layer#0
  for(std::map<std::string,int>::const_iterator it=_families.begin();it!=_families.end();it++)
    {
      tinyStr.push_back((*it).first);
      layer0.push_back((*it).second);
    }
  layer0.push_back((int)_groups.size());//3+aa i <- key info bb layer#0
  for(std::map<std::string, std::vector<std::string> >::const_iterator it0=_groups.begin();it0!=_groups.end();it0++)
    {
      layer0.push_back((int)(*it0).second.size());
      tinyStr.push_back((*it0).first);
      for(std::vector<std::string>::const_iterator it1=((*it0).second).begin();it1!=((*it0).second).end();it1++)
        tinyStr.push_back(*it1);
    }
  // sizeof(layer0)==3+aa+1+bb layer#0
  bigArrayD=_coords;// 0 bd
  bigArraysI.push_back(_fam_coords);// 0 bi
  bigArraysI.push_back(_num_coords);// 1 bi
  const PartDefinition *pd(_part_coords);
  if(!pd)
    layer0.push_back(-1);
  else
    {
      std::vector<int> tmp0;
      pd->serialize(tmp0,bigArraysI);
      tinyInt.push_back(tmp0.size());
      tinyInt.insert(tinyInt.end(),tmp0.begin(),tmp0.end());
    }
  //
  std::vector<int> layer1;
  std::vector<int> levs(getNonEmptyLevels());
  layer1.push_back((int)levs.size());// 0 i <- key
  layer1.insert(layer1.end(),levs.begin(),levs.end());
  for(std::vector<int>::const_iterator it=levs.begin();it!=levs.end();it++)
    {
      const MEDFileUMeshSplitL1 *lev(getMeshAtLevSafe(*it));
      lev->serialize(layer1,bigArraysI);
    }
  // put layers all together.
  tinyInt.push_back(layer0.size());
  tinyInt.insert(tinyInt.end(),layer0.begin(),layer0.end());
  tinyInt.push_back(layer1.size());
  tinyInt.insert(tinyInt.end(),layer1.begin(),layer1.end());
}

void MEDFileUMesh::unserialize(std::vector<double>& tinyDouble, std::vector<int>& tinyInt, std::vector<std::string>& tinyStr,
                               std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> >& bigArraysI, MEDCouplingAutoRefCountObjectPtr<DataArrayDouble>& bigArrayD)
{
  int sz0(tinyInt[0]);
  std::vector<int> layer0(tinyInt.begin()+1,tinyInt.begin()+1+sz0);
  int sz1(tinyInt[sz0+1]);
  std::vector<int> layer1(tinyInt.begin()+2+sz0,tinyInt.begin()+2+sz0+sz1);
  //
  std::reverse(layer0.begin(),layer0.end());
  std::reverse(layer1.begin(),layer1.end());
  std::reverse(tinyDouble.begin(),tinyDouble.end());
  std::reverse(tinyStr.begin(),tinyStr.end());
  std::reverse(bigArraysI.begin(),bigArraysI.end());
  //
  _order=layer0.back(); layer0.pop_back();
  _iteration=layer0.back(); layer0.pop_back();
  int spaceDim(layer0.back()); layer0.pop_back();
  _time=tinyDouble.back(); tinyDouble.pop_back();
  _name=tinyStr.back(); tinyStr.pop_back();
  _desc_name=tinyStr.back(); tinyStr.pop_back();
  _coords=bigArrayD; _coords->rearrange(spaceDim);
  for(int i=0;i<spaceDim;i++)
    {
      _coords->setInfoOnComponent(i,tinyStr.back());
      tinyStr.pop_back();
    }
  int nbOfFams(layer0.back()); layer0.pop_back();
  _families.clear();
  for(int i=0;i<nbOfFams;i++)
    {
      _families[tinyStr.back()]=layer0.back();
      tinyStr.pop_back(); layer0.pop_back();
    }
  int nbGroups(layer0.back()); layer0.pop_back();
  _groups.clear();
  for(int i=0;i<nbGroups;i++)
    {
      std::string grpName(tinyStr.back()); tinyStr.pop_back();
      int nbOfFamsOnGrp(layer0.back()); layer0.pop_back();
      std::vector<std::string> fams(nbOfFamsOnGrp);
      for(int j=0;j<nbOfFamsOnGrp;j++)
        {
          fams[j]=tinyStr.back(); tinyStr.pop_back();
        }
      _groups[grpName]=fams;
    }
  _fam_coords=bigArraysI.back(); bigArraysI.pop_back();
  _num_coords=bigArraysI.back(); bigArraysI.pop_back();
  _part_coords=0;
  int isPd(layer0.back()); layer0.pop_back();
  if(isPd!=-1)
    {
      std::vector<int> tmp0(layer0.begin(),layer0.begin()+isPd);
      layer0.erase(layer0.begin(),layer0.begin()+isPd);
      _part_coords=PartDefinition::Unserialize(tmp0,bigArraysI);
    }
  if(!layer0.empty())
    throw INTERP_KERNEL::Exception("MEDFileUMesh::unserialize : something wrong during unserialization #1 !");
  //
  int nbLevs(layer1.back()); layer1.pop_back();
  std::vector<int> levs(layer1.rbegin(),layer1.rbegin()+nbLevs); layer1.erase(layer1.end()-nbLevs,layer1.end());
  _ms.clear();
  int maxLev(-(*std::min_element(levs.begin(),levs.end())));
  _ms.resize(maxLev+1);
  for(int i=0;i<nbLevs;i++)
    {
      int lev(levs[i]);
      int pos(-lev);
      _ms[pos]=MEDFileUMeshSplitL1::Unserialize(_name,_coords,layer1,bigArraysI);
    }
}

/*!
 * Adds a group of nodes to \a this mesh.
 *  \param [in] ids - a DataArrayInt providing ids and a name of the group to add.
 *          The ids should be sorted and different each other (MED file norm).
 *
 *  \warning this method can alter default "FAMILLE_ZERO" family.
 *  For users sensitive to this a call to MEDFileMesh::rearrangeFamilies will be necessary after addGroup session.
 *
 *  \throw If the node coordinates array is not set.
 *  \throw If \a ids == \c NULL.
 *  \throw If \a ids->getName() == "".
 *  \throw If \a ids does not respect the MED file norm.
 *  \throw If a group with name \a ids->getName() already exists.
 */
void MEDFileUMesh::addNodeGroup(const DataArrayInt *ids)
{
  const DataArrayDouble *coords(_coords);
  if(!coords)
    throw INTERP_KERNEL::Exception("MEDFileUMesh::addNodeGroup : no coords set !");
  int nbOfNodes(coords->getNumberOfTuples());
  if(!((DataArrayInt *)_fam_coords))
    { _fam_coords=DataArrayInt::New(); _fam_coords->alloc(nbOfNodes,1); _fam_coords->fillWithZero(); }
  //
  addGroupUnderground(true,ids,_fam_coords);
}

/*!
 * Adds a group of nodes/cells/faces/edges to \a this mesh.
 *
 *  \param [in] ids - a DataArrayInt providing ids and a name of the group to add.
 *          The ids should be sorted and different each other (MED file norm).
 *
 * \warning this method can alter default "FAMILLE_ZERO" family.
 * For users sensitive to this a call to MEDFileMesh::rearrangeFamilies will be necessary after addGroup session.
 *
 *  \throw If the node coordinates array is not set.
 *  \throw If \a ids == \c NULL.
 *  \throw If \a ids->getName() == "".
 *  \throw If \a ids does not respect the MED file norm.
 *  \throw If a group with name \a ids->getName() already exists.
 */
void MEDFileUMesh::addGroup(int meshDimRelToMaxExt, const DataArrayInt *ids)
{
  std::vector<int> levs(getNonEmptyLevelsExt());
  if(std::find(levs.begin(),levs.end(),meshDimRelToMaxExt)==levs.end())
    { 
      std::ostringstream oss; oss << "MEDFileUMesh::addGroup : level " << meshDimRelToMaxExt << " not available ! Should be in ";
      std::copy(levs.begin(),levs.end(),std::ostream_iterator<int>(oss," ")); oss << " !"; throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  if(meshDimRelToMaxExt==1)
    { addNodeGroup(ids); return ; }
  MEDFileUMeshSplitL1 *lev(getMeshAtLevSafe(meshDimRelToMaxExt));
  DataArrayInt *fam(lev->getOrCreateAndGetFamilyField());
  addGroupUnderground(false,ids,fam);
}

/*!
 * Changes a name of a family specified by its id.
 *  \param [in] id - the id of the family of interest.
 *  \param [in] newFamName - the new family name.
 *  \throw If no family with the given \a id exists.
 */
void MEDFileUMesh::setFamilyNameAttachedOnId(int id, const std::string& newFamName)
{
  std::string oldName=getFamilyNameGivenId(id);
  _families.erase(oldName);
  _families[newFamName]=id;
}

/*!
 * Removes a mesh of a given dimension.
 *  \param [in] meshDimRelToMax - the relative dimension of interest.
 *  \throw If there is no mesh at level \a meshDimRelToMax in \a this mesh.
 */
void MEDFileUMesh::removeMeshAtLevel(int meshDimRelToMax)
{
  std::vector<int> levSet=getNonEmptyLevels();
  std::vector<int>::const_iterator it=std::find(levSet.begin(),levSet.end(),meshDimRelToMax);
  if(it==levSet.end())
    throw INTERP_KERNEL::Exception("MEDFileUMesh::removeMeshAtLevel : the requested level is not existing !");
  int pos=(-meshDimRelToMax);
  _ms[pos]=0;
}

/*!
 * Sets a new MEDCoupling1GTUMesh at a given level in \a this mesh.
 *  \param [in] meshDimRelToMax - a relative level to set the mesh at.
 *  \param [in] m - the new mesh to set.
 *  \throw If the name or the description of \a this mesh and \a m are not empty and are
 *         different. 
 *  \throw If the node coordinates array is set \a this in mesh and \a m refers to
 *         another node coordinates array.
 *  \throw If the mesh dimension of \a m does not correspond to \a meshDimRelToMax or
 *         to the existing meshes of other levels of \a this mesh.
 */
void MEDFileUMesh::setMeshAtLevel(int meshDimRelToMax, MEDCoupling1GTUMesh *m)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> elt(new MEDFileUMeshSplitL1(m));
  checkAndGiveEntryInSplitL1(meshDimRelToMax,m)=elt;
}

/*!
 * Sets a new MEDCouplingUMesh at a given level in \a this mesh.
 *  \param [in] meshDimRelToMax - a relative level to set the mesh at.
 *  \param [in] m - the new mesh to set.
 *  \param [in] newOrOld - if \c true, cells in \a m are sorted by type to be ready for 
 *         writing \a this mesh in a MED file.
 *  \throw If the name or the description of \a this mesh and \a m are not empty and are
 *         different. 
 *  \throw If the node coordinates array is set \a this in mesh and \a m refers to
 *         another node coordinates array.
 *  \throw If the mesh dimension of \a m does not correspond to \a meshDimRelToMax or
 *         to the existing meshes of other levels of \a this mesh.
 */
void MEDFileUMesh::setMeshAtLevel(int meshDimRelToMax, MEDCouplingUMesh *m, bool newOrOld)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> elt(new MEDFileUMeshSplitL1(m,newOrOld));
  checkAndGiveEntryInSplitL1(meshDimRelToMax,m)=elt;
}

MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1>& MEDFileUMesh::checkAndGiveEntryInSplitL1(int meshDimRelToMax, MEDCouplingPointSet *m)
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
      return _ms[sz-1];
    }
  else
    return _ms[-meshDimRelToMax];
}

/*!
 * This method allows to set at once the content of different levels in \a this.
 * This method is equivalent to a series of call to MEDFileUMesh::setMeshAtLevel.
 *
 * \param [in] ms - List of unstructured meshes lying on the same coordinates and having different mesh dimesnion.
 * \param [in] renum - the parameter (set to false by default) that tells the beheviour if there is a mesh on \a ms that is not geo type sorted.
 *                     If false, an exception ois thrown. If true the mesh is reordered automatically. It is highly recommanded to let this parameter to false.
 *
 * \throw If \a there is a null pointer in \a ms.
 * \sa MEDFileUMesh::setMeshAtLevel
 */
void MEDFileUMesh::setMeshes(const std::vector<const MEDCouplingUMesh *>& ms, bool renum)
{
  if(ms.empty())
    return ;
  const MEDCouplingUMesh *mRef=ms[0];
  if(!mRef)
    throw INTERP_KERNEL::Exception("MEDFileUMesh::setMeshes : null instance in the first element of input meshes !");
  std::string name(mRef->getName());
  const DataArrayDouble *coo(mRef->getCoords());
  std::set<int> s;
  int zeDim=-1;
  for(std::vector<const MEDCouplingUMesh *>::const_iterator it=ms.begin();it!=ms.end();it++)
    {
      const MEDCouplingUMesh *cur(*it);
      if(!cur)
        throw INTERP_KERNEL::Exception("MEDFileUMesh::setMeshes : null instance in input vector of meshes !");
      if(coo!=cur->getCoords())
        throw INTERP_KERNEL::Exception("MEDFileUMesh::setMeshes : The input meshes do not share the same coordinates !");
      int mdim=cur->getMeshDimension();
      zeDim=std::max(zeDim,mdim);
      if(s.find(mdim)!=s.end())
        throw INTERP_KERNEL::Exception("MEDFileUMesh::setMeshes : The input meshes must share the same coordinates pointer, and should have different mesh dimension each other !");
    }
  for(std::vector<const MEDCouplingUMesh *>::const_iterator it=ms.begin();it!=ms.end();it++)
    {
      int mdim=(*it)->getMeshDimension();
      setName((*it)->getName());
      setMeshAtLevel(mdim-zeDim,const_cast<MEDCouplingUMesh *>(*it),renum);
    }
  setName(name);
}

/*!
 * Creates one MEDCouplingUMesh at a given level in \a this mesh from a sequence of
 * meshes each representing a group, and creates corresponding groups in \a this mesh.
 * The given meshes must share the same node coordinates array.
 *  \param [in] meshDimRelToMax - the relative dimension to create the mesh and groups at.
 *  \param [in] ms - the sequence of meshes. Each mesh in \a ms represents a group to
 *          create in \a this mesh.
 *  \throw If \a ms is empty.
 *  \throw If dimension of meshes in \a ms does not correspond to \a meshDimRelToMax or
 *         to the existing meshes of other levels of \a this mesh.
 *  \throw If the meshes in \a ms do not share the same node coordinates array.
 *  \throw If the node coordinates array of \a this mesh (if any) is not the same as that
 *         of the given meshes.
 *  \throw If \a ms[ i ] is not well defined (MEDCouplingUMesh::checkCoherency()).
 *  \throw If names of some meshes in \a ms are equal.
 *  \throw If \a ms includes a mesh with an empty name.
 */
void MEDFileUMesh::setGroupsFromScratch(int meshDimRelToMax, const std::vector<const MEDCouplingUMesh *>& ms, bool renum)
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
  setMeshAtLevel(meshDimRelToMax,m,renum);
  std::vector<const DataArrayInt *> corr2(corr.begin(),corr.end());
  setGroupsAtLevel(meshDimRelToMax,corr2,true);
}

/*!
 * Creates groups at a given level in \a this mesh from a sequence of
 * meshes each representing a group.
 * The given meshes must share the same node coordinates array.
 *  \param [in] meshDimRelToMax - the relative dimension to create the groups at.
 *  \param [in] ms - the sequence of meshes. Each mesh in \a ms represents a group to
 *         create in \a this mesh.
 *  \param [in] renum - if \c true, then the optional numbers of entities are taken into
 *         account. 
 *  \throw If \a ms is empty.
 *  \throw If dimension of meshes in \a ms does not correspond to \a meshDimRelToMax or
 *         to the existing meshes of other levels of \a this mesh.
 *  \throw If the meshes in \a ms do not share the same node coordinates array.
 *  \throw If the node coordinates array of \a this mesh (if any) is not the same as that
 *         of the given meshes.
 *  \throw If \a ms[ i ] is not well defined (MEDCouplingUMesh::checkCoherency()).
 *  \throw If names of some meshes in \a ms are equal.
 *  \throw If \a ms includes a mesh with an empty name.
 */
void MEDFileUMesh::setGroupsOnSetMesh(int meshDimRelToMax, const std::vector<const MEDCouplingUMesh *>& ms, bool renum)
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

DataArrayDouble *MEDFileUMesh::checkMultiMesh(const std::vector<const MEDCouplingUMesh *>& ms) const
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

/*!
 * Sets the family field of a given relative dimension.
 *  \param [in] meshDimRelToMaxExt - the relative dimension of entities for which
 *          the family field is set.
 *  \param [in] famArr - the array of the family field.
 *  \throw If there are no mesh entities of \a meshDimRelToMaxExt dimension in \a this mesh.
 *  \throw If \a famArr has an invalid size.
 */
void MEDFileUMesh::setFamilyFieldArr(int meshDimRelToMaxExt, DataArrayInt *famArr)
{
  if(meshDimRelToMaxExt==1)
    {
      if(!famArr)
        {
          _fam_coords=0;
          return ;
        }
      DataArrayDouble *coo(_coords);
      if(!coo)
        throw INTERP_KERNEL::Exception("MEDFileUMesh::setFamilyFieldArr : the coordinates have not been set !");
      famArr->checkNbOfTuplesAndComp(coo->getNumberOfTuples(),1,"MEDFileUMesh::setFamilyFieldArr : Problem in size of node family arr ! ");
      famArr->incrRef();
      _fam_coords=famArr;
      return ;
    }
  if(meshDimRelToMaxExt>1)
    throw INTERP_KERNEL::Exception("MEDFileUMesh::setFamilyFieldArr : Dimension request is invalid (>1) !");
  int traducedRk=-meshDimRelToMaxExt;
  if(traducedRk>=(int)_ms.size())
    throw INTERP_KERNEL::Exception("Invalid mesh dim relative to max given ! Too low !");
  if((MEDFileUMeshSplitL1 *)_ms[traducedRk]==0)
    throw INTERP_KERNEL::Exception("On specified lev (or entity) no cells exists !");
  return _ms[traducedRk]->setFamilyArr(famArr);
}

/*!
 * Sets the optional numbers of mesh entities of a given dimension.
 *  \param [in] meshDimRelToMaxExt - the relative dimension of mesh entities.
 *  \param [in] renumArr - the array of the numbers.
 *  \throw If there are no mesh entities of \a meshDimRelToMaxExt dimension in \a this mesh.
 *  \throw If \a renumArr has an invalid size.
 */
void MEDFileUMesh::setRenumFieldArr(int meshDimRelToMaxExt, DataArrayInt *renumArr)
{
  if(meshDimRelToMaxExt==1)
    {
      if(!renumArr)
        {
          _num_coords=0;
          _rev_num_coords=0;
          return ;
        }
      DataArrayDouble *coo(_coords);
      if(!coo)
        throw INTERP_KERNEL::Exception("MEDFileUMesh::setRenumFieldArr : the coordinates have not been set !");
      renumArr->checkNbOfTuplesAndComp(coo->getNumberOfTuples(),1,"MEDFileUMesh::setRenumArr : Problem in size of node numbering arr ! ");
      renumArr->incrRef();
      _num_coords=renumArr;
      computeRevNum();
      return ;
    }
  if(meshDimRelToMaxExt>1)
    throw INTERP_KERNEL::Exception("MEDFileUMesh::setRenumArr : Dimension request is invalid (>1) !");
  int traducedRk=-meshDimRelToMaxExt;
  if(traducedRk>=(int)_ms.size())
    throw INTERP_KERNEL::Exception("Invalid mesh dim relative to max given ! Too low !");
  if((MEDFileUMeshSplitL1 *)_ms[traducedRk]==0)
    throw INTERP_KERNEL::Exception("On specified lev (or entity) no cells exists !");
  return _ms[traducedRk]->setRenumArr(renumArr);
}

/*!
 * Sets the optional names of mesh entities of a given dimension.
 *  \param [in] meshDimRelToMaxExt - the relative dimension of mesh entities.
 *  \param [in] nameArr - the array of the names.
 *  \throw If there are no mesh entities of \a meshDimRelToMaxExt dimension in \a this mesh.
 *  \throw If \a nameArr has an invalid size.
 */
void MEDFileUMesh::setNameFieldAtLevel(int meshDimRelToMaxExt, DataArrayAsciiChar *nameArr)
{
  if(meshDimRelToMaxExt==1)
    {
      if(!nameArr)
        {
          _name_coords=0;
          return ;
        }
      DataArrayDouble *coo(_coords);
      if(!coo)
        throw INTERP_KERNEL::Exception("MEDFileUMesh::setNameFieldAtLevel : the coordinates have not been set !");
      nameArr->checkNbOfTuplesAndComp(coo->getNumberOfTuples(),MED_SNAME_SIZE,"MEDFileUMesh::setNameFieldAtLevel : Problem in size of node numbering arr ! ");
      nameArr->incrRef();
      _name_coords=nameArr;
      return ;
    }
  if(meshDimRelToMaxExt>1)
    throw INTERP_KERNEL::Exception("MEDFileUMesh::setNameFieldAtLevel : Dimension request is invalid (>1) !");
  int traducedRk=-meshDimRelToMaxExt;
  if(traducedRk>=(int)_ms.size())
    throw INTERP_KERNEL::Exception("Invalid mesh dim relative to max given ! Too low !");
  if((MEDFileUMeshSplitL1 *)_ms[traducedRk]==0)
    throw INTERP_KERNEL::Exception("On specified lev (or entity) no cells exists !");
  return _ms[traducedRk]->setNameArr(nameArr);
}

void MEDFileUMesh::synchronizeTinyInfoOnLeaves() const
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::const_iterator it=_ms.begin();it!=_ms.end();it++)
    if((const MEDFileUMeshSplitL1 *)(*it))
      (*it)->synchronizeTinyInfo(*this);
}

/*!
 * This method is called by MEDFileMesh::changeFamilyId. It performs only one part of the family id modification.
 */
void MEDFileUMesh::changeFamilyIdArr(int oldId, int newId)
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

std::list< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > MEDFileUMesh::getAllNonNullFamilyIds() const
{
  std::list< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > ret;
  const DataArrayInt *da(_fam_coords);
  if(da)
    { da->incrRef(); ret.push_back(MEDCouplingAutoRefCountObjectPtr<DataArrayInt>(const_cast<DataArrayInt *>(da))); }
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> >::const_iterator it=_ms.begin();it!=_ms.end();it++)
    {
      const MEDFileUMeshSplitL1 *elt(*it);
      if(elt)
        {
          da=elt->getFamilyField();
          if(da)
            { da->incrRef(); ret.push_back(MEDCouplingAutoRefCountObjectPtr<DataArrayInt>(const_cast<DataArrayInt *>(da))); }
        }
    }
  return ret;
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

std::size_t MEDFileStructuredMesh::getHeapMemorySizeWithoutChildren() const
{
  return MEDFileMesh::getHeapMemorySizeWithoutChildren();
}

std::vector<const BigMemoryObject *> MEDFileStructuredMesh::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret(MEDFileMesh::getDirectChildrenWithNull());
  ret.push_back((const DataArrayInt *)_fam_nodes);
  ret.push_back((const DataArrayInt *)_num_nodes);
  ret.push_back((const DataArrayAsciiChar *)_names_nodes);
  ret.push_back((const DataArrayInt *)_fam_cells);
  ret.push_back((const DataArrayInt *)_num_cells);
  ret.push_back((const DataArrayAsciiChar *)_names_cells);
  ret.push_back((const DataArrayInt *)_fam_faces);
  ret.push_back((const DataArrayInt *)_num_faces);
  ret.push_back((const DataArrayInt *)_rev_num_nodes);
  ret.push_back((const DataArrayAsciiChar *)_names_faces);
  ret.push_back((const DataArrayInt *)_rev_num_cells);
  ret.push_back((const MEDCoupling1SGTUMesh*)_faces_if_necessary);
  return ret;
}

int MEDFileStructuredMesh::getMaxAbsFamilyIdInArrays() const
{
  int ret=-std::numeric_limits<int>::max(),tmp=-1;
  if((const DataArrayInt *)_fam_nodes)
    {
      int val=_fam_nodes->getMaxValue(tmp);
      ret=std::max(ret,std::abs(val));
    }
  if((const DataArrayInt *)_fam_cells)
    {
      int val=_fam_cells->getMaxValue(tmp);
      ret=std::max(ret,std::abs(val));
    }
  if((const DataArrayInt *)_fam_faces)
    {
      int val=_fam_faces->getMaxValue(tmp);
      ret=std::max(ret,std::abs(val));
    }
  return ret;
}

int MEDFileStructuredMesh::getMaxFamilyIdInArrays() const
{
  int ret=-std::numeric_limits<int>::max(),tmp=-1;
  if((const DataArrayInt *)_fam_nodes)
    {
      int val=_fam_nodes->getMaxValue(tmp);
      ret=std::max(ret,val);
    }
  if((const DataArrayInt *)_fam_cells)
    {
      int val=_fam_cells->getMaxValue(tmp);
      ret=std::max(ret,val);
    }
  if((const DataArrayInt *)_fam_faces)
    {
      int val=_fam_faces->getMaxValue(tmp);
      ret=std::max(ret,val);
    }
  return ret;
}

int MEDFileStructuredMesh::getMinFamilyIdInArrays() const
{
  int ret=std::numeric_limits<int>::max(),tmp=-1;
  if((const DataArrayInt *)_fam_nodes)
    {
      int val=_fam_nodes->getMinValue(tmp);
      ret=std::min(ret,val);
    }
  if((const DataArrayInt *)_fam_cells)
    {
      int val=_fam_cells->getMinValue(tmp);
      ret=std::min(ret,val);
    }
  if((const DataArrayInt *)_fam_faces)
    {
      int val=_fam_faces->getMinValue(tmp);
      ret=std::min(ret,val);
    }
  return ret;
}

bool MEDFileStructuredMesh::isEqual(const MEDFileMesh *other, double eps, std::string& what) const
{
  if(!MEDFileMesh::isEqual(other,eps,what))
    return false;
  const MEDFileStructuredMesh *otherC=dynamic_cast<const  MEDFileStructuredMesh *>(other);
  if(!otherC)
    {
      what="Mesh types differ ! This is structured and other is NOT !";
      return false;
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
  famc1=_fam_faces;
  famc2=otherC->_fam_faces;
  if((famc1==0 && famc2!=0) || (famc1!=0 && famc2==0))
    {
      what="Mismatch of families arr on faces ! One is defined and not other !";
      return false;
    }
  if(famc1)
    {
      bool ret=famc1->isEqual(*famc2);
      if(!ret)
        {
          what="Families arr on faces differ !";
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
  famc1=_num_faces;
  famc2=otherC->_num_faces;
  if((famc1==0 && famc2!=0) || (famc1!=0 && famc2==0))
    {
      what="Mismatch of numbering arr on faces ! One is defined and not other !";
      return false;
    }
  if(famc1)
    {
      bool ret=famc1->isEqual(*famc2);
      if(!ret)
        {
          what="Numbering arr on faces differ !";
          return false;
        }
    }
  const DataArrayAsciiChar *d1=_names_cells;
  const DataArrayAsciiChar *d2=otherC->_names_cells;
  if((d1==0 && d2!=0) || (d1!=0 && d2==0))
    {
      what="Mismatch of naming arr on cells ! One is defined and not other !";
      return false;
    }
  if(d1)
    {
      bool ret=d1->isEqual(*d2);
      if(!ret)
        {
          what="Naming arr on cells differ !";
          return false;
        }
    }
  d1=_names_faces;
  d2=otherC->_names_faces;
  if((d1==0 && d2!=0) || (d1!=0 && d2==0))
    {
      what="Mismatch of naming arr on faces ! One is defined and not other !";
      return false;
    }
  if(d1)
    {
      bool ret=d1->isEqual(*d2);
      if(!ret)
        {
          what="Naming arr on faces differ !";
          return false;
        }
    }
  d1=_names_nodes;
  d2=otherC->_names_nodes;
  if((d1==0 && d2!=0) || (d1!=0 && d2==0))
    {
      what="Mismatch of naming arr on nodes ! One is defined and not other !";
      return false;
    }
  if(d1)
    {
      bool ret=d1->isEqual(*d2);
      if(!ret)
        {
          what="Naming arr on nodes differ !";
          return false;
        }
    }
  return true;
}

void MEDFileStructuredMesh::clearNonDiscrAttributes() const
{
  MEDFileMesh::clearNonDiscrAttributes();
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
  tmp=_fam_faces;
  if(tmp)
    (const_cast<DataArrayInt *>(tmp))->setName("");
  tmp=_num_faces;
  if(tmp)
    (const_cast<DataArrayInt *>(tmp))->setName("");
}

/*!
 * Returns ids of mesh entities contained in given families of a given dimension.
 *  \param [in] meshDimRelToMaxExt - a relative dimension of the mesh entities whose ids
 *          are required.
 *  \param [in] fams - the names of the families of interest.
 *  \param [in] renum - if \c true, the optional numbers of entities, if available, are
 *          returned instead of ids.
 *  \return DataArrayInt * - a new instance of DataArrayInt holding either ids or
 *          numbers, if available and required, of mesh entities of the families. The caller
 *          is to delete this array using decrRef() as it is no more needed. 
 *  \throw If the family field is missing for \a meshDimRelToMaxExt.
 */
DataArrayInt *MEDFileStructuredMesh::getFamiliesArr(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum) const
{
  std::vector<int> famIds(getFamiliesIds(fams));
  switch(meshDimRelToMaxExt)
  {
    case 1:
      {
        if((const DataArrayInt *)_fam_nodes)
          {
            MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da;
            if(!famIds.empty())
              da=_fam_nodes->getIdsEqualList(&famIds[0],&famIds[0]+famIds.size());
            else
              da=_fam_nodes->getIdsEqualList(0,0);
            if(renum)
              return MEDFileUMeshSplitL1::Renumber(_num_nodes,da);
            else
              return da.retn();
          }
        else
          throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::getFamiliesArr : no family array specified on nodes !");
        break;
      }
    case 0:
      {
        if((const DataArrayInt *)_fam_cells)
          {
            MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da;
            if(!famIds.empty())
              da=_fam_cells->getIdsEqualList(&famIds[0],&famIds[0]+famIds.size());
            else
              da=_fam_cells->getIdsEqualList(0,0);
            if(renum)
              return MEDFileUMeshSplitL1::Renumber(_num_cells,da);
            else
              return da.retn();
          }
        else
          throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::getFamiliesArr : no family array specified on cells !");
        break;
      }
    case -1:
      {
        if((const DataArrayInt *)_fam_faces)
          {
            MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da;
            if(!famIds.empty())
              da=_fam_faces->getIdsEqualList(&famIds[0],&famIds[0]+famIds.size());
            else
              da=_fam_faces->getIdsEqualList(0,0);
            if(renum)
              return MEDFileUMeshSplitL1::Renumber(_num_faces,da);
            else
              return da.retn();
          }
        else
          throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::getFamiliesArr : no family array specified on faces !");
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::getFamiliesArr : input meshDimRelative must be in [0,1,-1] !");
  }
  throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::getFamiliesArr : unmanaged case !");
}

/*!
 * Sets the family field of a given relative dimension.
 *  \param [in] meshDimRelToMaxExt - the relative dimension of entities for which
 *          the family field is set.
 *  \param [in] famArr - the array of the family field.
 *  \throw If there are no mesh entities of \a meshDimRelToMaxExt dimension in \a this mesh.
 *  \throw If \a famArr has an invalid size.
 *  \throw If \a meshDimRelToMaxExt != 0 and \a meshDimRelToMaxExt != 1 and \a meshDimRelToMaxExt != -1.
 */
void MEDFileStructuredMesh::setFamilyFieldArr(int meshDimRelToMaxExt, DataArrayInt *famArr)
{
  const MEDCouplingStructuredMesh *mesh(getStructuredMesh());
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::setFamilyFieldArr : no structured mesh specified ! Impossible to set family array !");
  switch(meshDimRelToMaxExt)
  {
    case 0:
      {
        int nbCells=mesh->getNumberOfCells();
        famArr->checkNbOfTuplesAndComp(nbCells,1,"MEDFileStructuredMesh::setFamilyFieldArr : Problem in size of Family arr ! Mismatch with number of cells of mesh !");
        _fam_cells=famArr;
        break;
      }
    case 1:
      {
        int nbNodes=mesh->getNumberOfNodes();
        famArr->checkNbOfTuplesAndComp(nbNodes,1,"MEDFileStructuredMesh::setFamilyFieldArr : Problem in size of Family arr ! Mismatch with number of nodes of mesh !");
        _fam_nodes=famArr;
        break;
      }
    case -1:
      {
        int nbCells=mesh->getNumberOfCellsOfSubLevelMesh();
        famArr->checkNbOfTuplesAndComp(nbCells,1,"MEDFileStructuredMesh::setFamilyFieldArr : Problem in size of Family arr ! Mismatch with number of faces of mesh !");
        _fam_faces=famArr;
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::setFamilyFieldArr : Only available for levels 0 or 1 or -1 !");
  }
  if(famArr)
    famArr->incrRef();
}

/*!
 * Sets the optional numbers of mesh entities of a given dimension.
 *  \param [in] meshDimRelToMaxExt - the relative dimension of mesh entities.
 *  \param [in] renumArr - the array of the numbers.
 *  \throw If there are no mesh entities of \a meshDimRelToMaxExt dimension in \a this mesh.
 *  \throw If \a renumArr has an invalid size.
 *  \throw If \a meshDimRelToMaxExt != 0 and \a meshDimRelToMaxExt != 1.
 */
void MEDFileStructuredMesh::setRenumFieldArr(int meshDimRelToMaxExt, DataArrayInt *renumArr)
{
  const MEDCouplingStructuredMesh *mesh=getStructuredMesh();
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::setRenumFieldArr : no structured mesh specified ! Impossible to set number array !");
  switch(meshDimRelToMaxExt)
  {
    case 0:
      {
        int nbCells=mesh->getNumberOfCells();
        renumArr->checkNbOfTuplesAndComp(nbCells,1,"MEDFileStructuredMesh::setRenumFieldArr : Problem in size of Renum arr ! Mismatch with number of cells of mesh !");
        _num_cells=renumArr;
        break;
      }
    case 1:
      {
        int nbNodes=mesh->getNumberOfNodes();
        renumArr->checkNbOfTuplesAndComp(nbNodes,1,"MEDFileStructuredMesh::setRenumFieldArr : Problem in size of Family arr ! Mismatch with number of nodes of mesh !");
        _num_nodes=renumArr;
        break;
      }
    case -1:
      {
        int nbCells=mesh->getNumberOfCellsOfSubLevelMesh();
        renumArr->checkNbOfTuplesAndComp(nbCells,1,"MEDFileStructuredMesh::setRenumFieldArr : Problem in size of Renum arr ! Mismatch with number of faces of mesh !");
        _num_faces=renumArr;
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::setRenumFieldArr : Only available for levels 0 or 1 or -1 !");
  }
  if(renumArr)
    renumArr->incrRef();
}

/*!
 * Sets the optional names of mesh entities of a given dimension.
 *  \param [in] meshDimRelToMaxExt - the relative dimension of mesh entities.
 *  \param [in] nameArr - the array of the names.
 *  \throw If there are no mesh entities of \a meshDimRelToMaxExt dimension in \a this mesh.
 *  \throw If \a nameArr has an invalid size.
 */
void MEDFileStructuredMesh::setNameFieldAtLevel(int meshDimRelToMaxExt, DataArrayAsciiChar *nameArr)
{
  const MEDCouplingStructuredMesh *mesh(getStructuredMesh());
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::setNameFieldAtLevel : no structured mesh specified ! Impossible to set names array !");
  switch(meshDimRelToMaxExt)
  {
    case 0:
      {
        int nbCells=mesh->getNumberOfCells();
        nameArr->checkNbOfTuplesAndComp(nbCells,MED_SNAME_SIZE,"MEDFileStructuredMesh::setNameFieldAtLevel : Problem in size of names arr ! Mismatch with number of cells of mesh !");
        _names_cells=nameArr;
        break;
      }
    case 1:
      {
        int nbNodes=mesh->getNumberOfNodes();
        nameArr->checkNbOfTuplesAndComp(nbNodes,MED_SNAME_SIZE,"MEDFileStructuredMesh::setNameFieldAtLevel : Problem in size of names arr ! Mismatch with number of nodes of mesh !");
        _names_nodes=nameArr;
        break;
      }
    case -1:
      {
        int nbCells=mesh->getNumberOfCellsOfSubLevelMesh();
        nameArr->checkNbOfTuplesAndComp(nbCells,MED_SNAME_SIZE,"MEDFileStructuredMesh::setNameFieldAtLevel : Problem in size of names arr ! Mismatch with number of faces of mesh !");
        _names_cells=nameArr;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::setNameFieldAtLevel : Only available for levels 0 or 1 or -1 !");
  }
  if(nameArr)
    nameArr->incrRef();
}

/*!
 * Adds a group of nodes to \a this mesh.
 *  \param [in] ids - a DataArrayInt providing ids and a name of the group to add.
 *          The ids should be sorted and different each other (MED file norm).
 *
 *  \warning this method can alter default "FAMILLE_ZERO" family.
 *  For users sensitive to this a call to MEDFileMesh::rearrangeFamilies will be necessary after addGroup session.
 *
 *  \throw If the node coordinates array is not set.
 *  \throw If \a ids == \c NULL.
 *  \throw If \a ids->getName() == "".
 *  \throw If \a ids does not respect the MED file norm.
 *  \throw If a group with name \a ids->getName() already exists.
 */
void MEDFileStructuredMesh::addNodeGroup(const DataArrayInt *ids)
{
  addGroup(1,ids);
}

/*!
 * Adds a group of nodes/cells/faces/edges to \a this mesh.
 *
 *  \param [in] ids - a DataArrayInt providing ids and a name of the group to add.
 *          The ids should be sorted and different each other (MED file norm).
 *
 * \warning this method can alter default "FAMILLE_ZERO" family.
 * For users sensitive to this a call to MEDFileMesh::rearrangeFamilies will be necessary after addGroup session.
 *
 *  \throw If the node coordinates array is not set.
 *  \throw If \a ids == \c NULL.
 *  \throw If \a ids->getName() == "".
 *  \throw If \a ids does not respect the MED file norm.
 *  \throw If a group with name \a ids->getName() already exists.
 */
void MEDFileStructuredMesh::addGroup(int meshDimRelToMaxExt, const DataArrayInt *ids)
{
  DataArrayInt *fam(getOrCreateAndGetFamilyFieldAtLevel(meshDimRelToMaxExt));
  addGroupUnderground(false,ids,fam);
  return ;
}

/*!
 * Returns the family field for mesh entities of a given dimension.
 *  \param [in] meshDimRelToMaxExt - the relative dimension of mesh entities.
 *  \return const DataArrayInt * - the family field. It is an array of ids of families
 *          each mesh entity belongs to. It can be \c NULL.
 *  \throw If \a meshDimRelToMaxExt != 0 and \a meshDimRelToMaxExt != 1.
 */
const DataArrayInt *MEDFileStructuredMesh::getFamilyFieldAtLevel(int meshDimRelToMaxExt) const
{
  switch(meshDimRelToMaxExt)
  {
    case 0:
      return _fam_cells;
    case 1:
      return _fam_nodes;
    case -1:
      return _fam_faces;
    default:
      throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::getFamilyFieldAtLevel : Only available for levels 0 or 1 or -1 !");
  }
}

/*!
 * Returns the family field for mesh entities of a given dimension.
 *  \param [in] meshDimRelToMaxExt - the relative dimension of mesh entities.
 *  \return const DataArrayInt * - the family field. It is an array of ids of families
 *          each mesh entity belongs to. It can be \c NULL.
 *  \throw If \a meshDimRelToMaxExt != 0 and \a meshDimRelToMaxExt != 1.
 */
DataArrayInt *MEDFileStructuredMesh::getFamilyFieldAtLevel(int meshDimRelToMaxExt)
{
  switch(meshDimRelToMaxExt)
  {
    case 0:
      return _fam_cells;
    case 1:
      return _fam_nodes;
    case -1:
      return _fam_faces;
    default:
      throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::getFamilyFieldAtLevel : Only available for levels 0 or 1 or -1 !");
  }
}

/*!
 * Returns the optional numbers of mesh entities of a given dimension.
 *  \param [in] meshDimRelToMaxExt - the relative dimension of mesh entities.
 *  \return const DataArrayInt * - the array of the entity numbers.
 *  \throw If there are no mesh entities of \a meshDimRelToMaxExt dimension in \a this mesh.
 *  \throw If \a meshDimRelToMaxExt != 0 and \a meshDimRelToMaxExt != 1.
 */
const DataArrayInt *MEDFileStructuredMesh::getNumberFieldAtLevel(int meshDimRelToMaxExt) const
{
  switch(meshDimRelToMaxExt)
  {
    case 0:
      return _num_cells;
    case 1:
      return _num_nodes;
    case -1:
      return _num_faces;
    default:
      throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::getNumberFieldAtLevel : Only available for levels 0 or 1 or -1 !");
  }
}

/*!
 * Returns the optional numbers of mesh entities of a given dimension transformed using
 * DataArrayInt::invertArrayN2O2O2N().
 *  \param [in] meshDimRelToMaxExt - the relative dimension of mesh entities.
 *  \return const DataArrayInt * - the array of the entity numbers transformed using
 *          DataArrayInt::invertArrayN2O2O2N().
 *  \throw If \a meshDimRelToMaxExt != 0 and \a meshDimRelToMaxExt != 1.
 *  \throw If there are no mesh entities of \a meshDimRelToMaxExt dimension in \a this mesh.
 */
const DataArrayInt *MEDFileStructuredMesh::getRevNumberFieldAtLevel(int meshDimRelToMaxExt) const
{
  if(meshDimRelToMaxExt!=0 && meshDimRelToMaxExt!=1)
    throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::getRevNumberFieldAtLevel : Only available for levels 0 or 1 !");
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

const DataArrayAsciiChar *MEDFileStructuredMesh::getNameFieldAtLevel(int meshDimRelToMaxExt) const
{
  switch(meshDimRelToMaxExt)
  {
    case 0:
      return _names_cells;
    case 1:
      return _names_nodes;
    case -1:
      return _names_faces;
    default:
      throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::getNameFieldAtLevel : Only available for levels 0 or 1 or -1 !");
  }
}

/*!
 * Returns relative dimensions of mesh entities (excluding nodes) present in \a this mesh.
 *  \return std::vector<int> - a sequence of the relative dimensions: [0].
 */
std::vector<int> MEDFileStructuredMesh::getNonEmptyLevels() const
{
  std::vector<int> ret(1);
  return ret;
}

/*!
 * Returns relative dimensions of mesh entities (including nodes) present in \a this mesh.
 *  \return std::vector<int> - a sequence of the relative dimensions: [1,0].
 */
std::vector<int> MEDFileStructuredMesh::getNonEmptyLevelsExt() const
{
  std::vector<int> ret(2);
  ret[0]=1;
  return ret;
}

/*!
 * Returns the set of extensive levels (nodes included) where not NULL family arr are defined.
 */
std::vector<int> MEDFileStructuredMesh::getFamArrNonEmptyLevelsExt() const
{
  std::vector<int> ret;
  const DataArrayInt *famNodes(_fam_nodes),*famCells(_fam_cells),*famFaces(_fam_faces);
  if(famNodes)
    ret.push_back(1);
  if(famCells)
    ret.push_back(0);
  if(famFaces)
    ret.push_back(-1);
  return ret;
}

/*!
 * Returns the set of extensive levels (nodes included) where not NULL numbering arr are defined.
 */
std::vector<int> MEDFileStructuredMesh::getNumArrNonEmptyLevelsExt() const
{
  std::vector<int> ret;
  const DataArrayInt *numNodes(_num_nodes),*numCells(_num_cells),*numFaces(_num_faces);
  if(numNodes)
    ret.push_back(1);
  if(numCells)
    ret.push_back(0);
  if(numFaces)
    ret.push_back(-1);
  return ret;
}

/*!
 * Returns the set of extensive levels (nodes included) where not NULL naming arr are defined.
 */
std::vector<int> MEDFileStructuredMesh::getNameArrNonEmptyLevelsExt() const
{
  std::vector<int> ret;
  const DataArrayAsciiChar *namesNodes(_names_nodes),*namesCells(_names_cells),*namesFaces(_names_faces);
  if(namesNodes)
    ret.push_back(1);
  if(namesCells)
    ret.push_back(0);
  if(namesFaces)
    ret.push_back(-1);
  return ret;
}

/*!
 * no implementation here, it is not a bug, but intresically no polyhedra in \a this.
 */
bool MEDFileStructuredMesh::unPolyze(std::vector<int>& oldCode, std::vector<int>& newCode, DataArrayInt *& o2nRenumCell)
{
  oldCode.clear(); newCode.clear(); o2nRenumCell=0;
  return false;
}

void MEDFileStructuredMesh::changeFamilyIdArr(int oldId, int newId)
{
  DataArrayInt *arr=_fam_nodes;
  if(arr)
    arr->changeValue(oldId,newId);
  arr=_fam_cells;
  if(arr)
    arr->changeValue(oldId,newId);
  arr=_fam_faces;
  if(arr)
    arr->changeValue(oldId,newId);
}

std::list< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > MEDFileStructuredMesh::getAllNonNullFamilyIds() const
{
  std::list< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > ret;
  const DataArrayInt *da(_fam_nodes);
  if(da)
    { da->incrRef(); ret.push_back(MEDCouplingAutoRefCountObjectPtr<DataArrayInt>(const_cast<DataArrayInt *>(da))); }
  da=_fam_cells;
  if(da)
    { da->incrRef(); ret.push_back(MEDCouplingAutoRefCountObjectPtr<DataArrayInt>(const_cast<DataArrayInt *>(da))); }
  da=_fam_faces;
  if(da)
    { da->incrRef(); ret.push_back(MEDCouplingAutoRefCountObjectPtr<DataArrayInt>(const_cast<DataArrayInt *>(da))); }
  return ret;
}

void MEDFileStructuredMesh::deepCpyAttributes()
{
  if((const DataArrayInt*)_fam_nodes)
    _fam_nodes=_fam_nodes->deepCpy();
  if((const DataArrayInt*)_num_nodes)
    _num_nodes=_num_nodes->deepCpy();
  if((const DataArrayAsciiChar*)_names_nodes)
    _names_nodes=_names_nodes->deepCpy();
  if((const DataArrayInt*)_fam_cells)
    _fam_cells=_fam_cells->deepCpy();
  if((const DataArrayInt*)_num_cells)
    _num_cells=_num_cells->deepCpy();
  if((const DataArrayAsciiChar*)_names_cells)
    _names_cells=_names_cells->deepCpy();
  if((const DataArrayInt*)_fam_faces)
    _fam_faces=_fam_faces->deepCpy();
  if((const DataArrayInt*)_num_faces)
    _num_faces=_num_faces->deepCpy();
  if((const DataArrayAsciiChar*)_names_faces)
    _names_faces=_names_faces->deepCpy();
  if((const DataArrayInt*)_rev_num_nodes)
    _rev_num_nodes=_rev_num_nodes->deepCpy();
  if((const DataArrayInt*)_rev_num_cells)
    _rev_num_cells=_rev_num_cells->deepCpy();
}

/*!
 * Returns a pointer to mesh at the specified level (here 0 is compulsary for cartesian mesh).
 * 
 * \return a pointer to cartesian mesh that need to be managed by the caller.
 * \warning the returned pointer has to be managed by the caller.
 */

/*!
 * Returns a pointer to MEDCouplingStructuredMesh held by \a this. 
 *  \param [in] meshDimRelToMax - it must be \c 0 or \c -1.
 *  \param [in] renum - it must be \c false.
 *  \return MEDCouplingMesh * - a pointer to MEDCouplingMesh that the caller is to
 *          delete using decrRef() as it is no more needed. 
 */
MEDCouplingMesh *MEDFileStructuredMesh::getGenMeshAtLevel(int meshDimRelToMax, bool renum) const
{
  if(renum)
    throw INTERP_KERNEL::Exception("MEDFileCurveLinearMesh does not support renumbering ! To do it perform request of renum array directly !");
  const MEDCouplingStructuredMesh *m(getStructuredMesh());
  switch(meshDimRelToMax)
  {
    case 0:
      {
        if(m)
          m->incrRef();
        return const_cast<MEDCouplingStructuredMesh *>(m);
      }
    case -1:
      {
        if(!m)
          throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::getGenMeshAtLevel : level -1 requested must be non empty to be able to compute unstructured sub mesh !");
        buildMinusOneImplicitPartIfNeeded();
        MEDCouplingMesh *ret(_faces_if_necessary);
        if(ret)
          ret->incrRef();
        return ret;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDFileCurveLinearMesh does not support multi level for mesh 0 expected as input !");
  }
}

/*!
 * Returns number of mesh entities of a given relative dimension in \a this mesh.
 *  \param [in] meshDimRelToMaxExt - the relative dimension of interest.
 *  \return int - the number of entities.
 *  \throw If no mesh entities of dimension \a meshDimRelToMaxExt are available in \a this mesh.
 */
int MEDFileStructuredMesh::getSizeAtLevel(int meshDimRelToMaxExt) const
{
  const MEDCouplingStructuredMesh *cmesh(getStructuredMesh());
  if(!cmesh)
    throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::getSizeAtLevel : No structured mesh set !");
  switch(meshDimRelToMaxExt)
  {
    case 0:
      return cmesh->getNumberOfCells();
    case 1:
      return cmesh->getNumberOfNodes();
    case -1:
      return cmesh->getNumberOfCellsOfSubLevelMesh();
    default:
      throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::getSizeAtLevel : Only available for levels 0 or 1 or -1 !");
  }
}

int MEDFileStructuredMesh::getNumberOfNodes() const
{
  const MEDCouplingStructuredMesh *cmesh(getStructuredMesh());
  if(!cmesh)
    throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::getNumberOfNodes : no cartesian mesh set !");
  return cmesh->getNumberOfNodes();
}

int MEDFileStructuredMesh::getNumberOfCellsAtLevel(int meshDimRelToMaxExt) const
{
  const MEDCouplingStructuredMesh *cmesh(getStructuredMesh());
  if(!cmesh)
    throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::getNumberOfNodes : no cartesian mesh set !");
  switch(meshDimRelToMaxExt)
  {
    case 0:
      return cmesh->getNumberOfCells();
    case -1:
      return cmesh->getNumberOfCellsOfSubLevelMesh();
    default:
      throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::getNumberOfNodes : only meshDimRelToMax=0 and meshDimRelToMax=-1 supported !");
    }
}

bool MEDFileStructuredMesh::hasImplicitPart() const
{
  return true;
}

/*!
 * \sa MEDFileStructuredMesh::getImplicitFaceMesh
 */
int MEDFileStructuredMesh::buildImplicitPartIfAny(INTERP_KERNEL::NormalizedCellType gt) const
{
  static const char MSG[]="MEDFileStructuredMesh::buildImplicitPartIfAny : the given geo type is not manageable by a structured mesh !";
  const MEDCoupling1SGTUMesh *zeFaceMesh(_faces_if_necessary);
  if(!zeFaceMesh)
    {
      const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel(MEDCouplingStructuredMesh::GetGeoTypeGivenMeshDimension(getMeshDimension())));
      if(cm.getReverseExtrudedType()!=gt)
        throw INTERP_KERNEL::Exception(MSG);
      buildImplicitPart();
      return getStructuredMesh()->getNumberOfCellsOfSubLevelMesh();
    }
  else
    {
      if(gt!=zeFaceMesh->getCellModelEnum())
        throw INTERP_KERNEL::Exception(MSG);
      return zeFaceMesh->getNumberOfCells();
    }
}

void MEDFileStructuredMesh::buildMinusOneImplicitPartIfNeeded() const
{
  const MEDCoupling1SGTUMesh *zeFaceMesh(_faces_if_necessary);
  if(!zeFaceMesh)
    buildImplicitPart();
}

void MEDFileStructuredMesh::buildImplicitPart() const
{
  const MEDCouplingStructuredMesh *mcmesh(getStructuredMesh());
  if(!mcmesh)
    throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::buildImplicitPart : Unable to build the implicit part of structured mesh because no structured mesh at level 0 defined !");
  _faces_if_necessary=mcmesh->build1SGTSubLevelMesh();
}

void MEDFileStructuredMesh::releaseImplicitPartIfAny() const
{
  _faces_if_necessary=0;
}

/*!
 * Retrieves the internal pointer (no decrRef requested) of the implicit face mesh if any.
 * To force to build it you can invoke MEDFileStructuredMesh::buildImplicitPartIfAny method.
 * 
 * \sa MEDFileStructuredMesh::buildImplicitPartIfAny
 */
MEDCoupling1SGTUMesh *MEDFileStructuredMesh::getImplicitFaceMesh() const
{
  return _faces_if_necessary;
}

std::vector<INTERP_KERNEL::NormalizedCellType> MEDFileStructuredMesh::getGeoTypesAtLevel(int meshDimRelToMax) const
{
  const MEDCouplingStructuredMesh *cmesh(getStructuredMesh());
  if(!cmesh)
    throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::getGeoTypesAtLevel : No structured mesh set !");
  switch(meshDimRelToMax)
  {
    case 0:
      {
        std::vector<INTERP_KERNEL::NormalizedCellType> ret(1,cmesh->getTypeOfCell(0));
        return ret;
      }
    case -1:
      {
        int mdim(cmesh->getMeshDimension());
        if(mdim<1)
          throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::getGeoTypesAtLevel : only one level available for structured meshes ! Input 0 is mandatory or 0D mesh !");
        std::vector<INTERP_KERNEL::NormalizedCellType> ret(1,MEDCouplingStructuredMesh::GetGeoTypeGivenMeshDimension(mdim-1));
        return ret;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::getGeoTypesAtLevel : only 2 levels available at most : 0 and -1 !");
  }
}

void MEDFileStructuredMesh::whichAreNodesFetched(const MEDFileField1TSStructItem& st, const MEDFileFieldGlobsReal *globs, std::vector<bool>& nodesFetched) const
{
  if(st.getNumberOfItems()!=1)
    throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::whichAreNodesFetched : The sturture of field is not lying on single geo type ! it is not managed yet for structured mesh !");
  if(st[0].getGeo()!=MEDCouplingStructuredMesh::GetGeoTypeGivenMeshDimension(getMeshDimension()))
    throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::whichAreNodesFetched : The sturture of field is not lying on expected geo type !");
  if(getNumberOfNodes()!=(int)nodesFetched.size())
    throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::whichAreNodesFetched : invalid size of array !");
  if(st[0].getPflName().empty())
    {
      std::fill(nodesFetched.begin(),nodesFetched.end(),true);
      return ;
    }
  const DataArrayInt *arr(globs->getProfile(st[0].getPflName()));
  const MEDCouplingStructuredMesh *cmesh=getStructuredMesh();//cmesh not null because getNumberOfNodes called before
  int sz(nodesFetched.size());
  for(const int *work=arr->begin();work!=arr->end();work++)
    {
      std::vector<int> conn;
      cmesh->getNodeIdsOfCell(*work,conn);
      for(std::vector<int>::const_iterator it=conn.begin();it!=conn.end();it++)
        if(*it>=0 && *it<sz)
          nodesFetched[*it]=true;
        else
          throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::whichAreNodesFetched : internal error !");
    }
}

med_geometry_type MEDFileStructuredMesh::GetGeoTypeFromMeshDim(int meshDim)
{
  INTERP_KERNEL::NormalizedCellType ct(MEDCouplingStructuredMesh::GetGeoTypeGivenMeshDimension(meshDim));
  return typmai3[ct];
}

void MEDFileStructuredMesh::LoadStrMeshDAFromFile(med_idt fid, int meshDim, int dt, int it, const std::string& mName, MEDFileMeshReadSelector *mrs,
                                                  MEDCouplingAutoRefCountObjectPtr<DataArrayInt>& famCells, MEDCouplingAutoRefCountObjectPtr<DataArrayInt>& numCells, MEDCouplingAutoRefCountObjectPtr<DataArrayAsciiChar>& namesCells)
{
  med_bool chgt=MED_FALSE,trsf=MED_FALSE;
  med_geometry_type geoTypeReq=MEDFileStructuredMesh::GetGeoTypeFromMeshDim(meshDim);
  int nbOfElt(0);
  nbOfElt=MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_CELL,geoTypeReq,MED_FAMILY_NUMBER,MED_NODAL,&chgt,&trsf);
  if(nbOfElt>0)
    {
      if(!mrs || mrs->isCellFamilyFieldReading())
        {
          famCells=DataArrayInt::New();
          famCells->alloc(nbOfElt,1);
          MEDFILESAFECALLERRD0(MEDmeshEntityFamilyNumberRd,(fid,mName.c_str(),dt,it,MED_CELL,geoTypeReq,famCells->getPointer()));
        }
    }
  nbOfElt=MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_CELL,geoTypeReq,MED_NUMBER,MED_NODAL,&chgt,&trsf);
  if(nbOfElt>0)
    {
      if(!mrs || mrs->isCellNumFieldReading())
        {
          numCells=DataArrayInt::New();
          numCells->alloc(nbOfElt,1);
          MEDFILESAFECALLERRD0(MEDmeshEntityNumberRd,(fid,mName.c_str(),dt,it,MED_CELL,geoTypeReq,numCells->getPointer()));
        }
    }
  nbOfElt=MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_CELL,geoTypeReq,MED_NAME,MED_NODAL,&chgt,&trsf);
  if(nbOfElt>0)
    {
      if(!mrs || mrs->isCellNameFieldReading())
        {
          namesCells=DataArrayAsciiChar::New();
          namesCells->alloc(nbOfElt+1,MED_SNAME_SIZE);//not a bug to avoid the memory corruption due to last \0 at the end
          MEDFILESAFECALLERRD0(MEDmeshEntityNameRd,(fid,mName.c_str(),dt,it,MED_CELL,geoTypeReq,namesCells->getPointer()));
          namesCells->reAlloc(nbOfElt);//not a bug to avoid the memory corruption due to last \0 at the end
        }
    }
}

void MEDFileStructuredMesh::loadStrMeshFromFile(MEDFileStrMeshL2 *strm, med_idt fid, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs)
{
  setName(strm->getName());
  setDescription(strm->getDescription());
  setUnivName(strm->getUnivName());
  setIteration(strm->getIteration());
  setOrder(strm->getOrder());
  setTimeValue(strm->getTime());
  setTimeUnit(strm->getTimeUnit());
  MEDFileMeshL2::ReadFamiliesAndGrps(fid,mName,_families,_groups,mrs);
  med_bool chgt=MED_FALSE,trsf=MED_FALSE;
  int nbOfElt(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NONE,MED_FAMILY_NUMBER,MED_NODAL,&chgt,&trsf));
  if(nbOfElt>0)
    {
      if(!mrs || mrs->isNodeFamilyFieldReading())
        {
          int nbNodes(getNumberOfNodes());
          if(nbOfElt>nbNodes)
            throw INTERP_KERNEL::Exception("MEDFileStructuredMesh::loadStrMeshFromFile : invalid size of family node array regarding number of nodes in this ! File seems to be corrupted !");
          _fam_nodes=DataArrayInt::New();
          _fam_nodes->alloc(nbNodes,1);//yes nbNodes and not nbOfElt see next line.
          if(nbNodes>nbOfElt)//yes it appends some times... It explains surely the mdump implementation. Bug revealed by PARAVIS EDF #2475 on structured.med file where only 12 first nodes are !=0 so nbOfElt=12 and nbOfNodes=378...
            _fam_nodes->fillWithZero();
          MEDFILESAFECALLERRD0(MEDmeshEntityFamilyNumberRd,(fid,mName.c_str(),dt,it,MED_NODE,MED_NONE,_fam_nodes->getPointer()));
        }
    }
  nbOfElt=MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NONE,MED_NUMBER,MED_NODAL,&chgt,&trsf);
  if(nbOfElt>0)
    {
      if(!mrs || mrs->isNodeNumFieldReading())
        {
          _num_nodes=DataArrayInt::New();
          _num_nodes->alloc(nbOfElt,1);
          MEDFILESAFECALLERRD0(MEDmeshEntityNumberRd,(fid,mName.c_str(),dt,it,MED_NODE,MED_NONE,_num_nodes->getPointer()));
        }
    }
  nbOfElt=MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NONE,MED_NAME,MED_NODAL,&chgt,&trsf);
  if(nbOfElt>0)
    {
      if(!mrs || mrs->isNodeNameFieldReading())
        {
          _names_nodes=DataArrayAsciiChar::New();
          _names_nodes->alloc(nbOfElt+1,MED_SNAME_SIZE);//not a bug to avoid the memory corruption due to last \0 at the end
          MEDFILESAFECALLERRD0(MEDmeshEntityNameRd,(fid,mName.c_str(),dt,it,MED_NODE,MED_NONE,_names_nodes->getPointer()));
          _names_nodes->reAlloc(nbOfElt);//not a bug to avoid the memory corruption due to last \0 at the end
        }
    }
  int meshDim(getStructuredMesh()->getMeshDimension());
  LoadStrMeshDAFromFile(fid,meshDim,dt,it,mName,mrs,_fam_cells,_num_cells,_names_cells);
  if(meshDim>=1)
    LoadStrMeshDAFromFile(fid,meshDim-1,dt,it,mName,mrs,_fam_faces,_num_faces,_names_faces);
}

void MEDFileStructuredMesh::writeStructuredLL(med_idt fid, const std::string& maa) const
{
  int meshDim(getStructuredMesh()->getMeshDimension());
  med_geometry_type geoTypeReq(GetGeoTypeFromMeshDim(meshDim)),geoTypeReq2(GetGeoTypeFromMeshDim(meshDim-1));
  //
  if((const DataArrayInt *)_fam_cells)
    MEDFILESAFECALLERWR0(MEDmeshEntityFamilyNumberWr,(fid,maa.c_str(),_iteration,_order,MED_CELL,geoTypeReq,_fam_cells->getNumberOfTuples(),_fam_cells->getConstPointer()));
  if((const DataArrayInt *)_fam_faces)
    MEDFILESAFECALLERWR0(MEDmeshEntityFamilyNumberWr,(fid,maa.c_str(),_iteration,_order,MED_CELL,geoTypeReq2,_fam_faces->getNumberOfTuples(),_fam_faces->getConstPointer()));
  if((const DataArrayInt *)_fam_nodes)
    MEDFILESAFECALLERWR0(MEDmeshEntityFamilyNumberWr,(fid,maa.c_str(),_iteration,_order,MED_NODE,MED_NONE,_fam_nodes->getNumberOfTuples(),_fam_nodes->getConstPointer()));
  if((const DataArrayInt *)_num_cells)
    MEDFILESAFECALLERWR0(MEDmeshEntityNumberWr,(fid,maa.c_str(),_iteration,_order,MED_CELL,geoTypeReq,_num_cells->getNumberOfTuples(),_num_cells->getConstPointer()));
  if((const DataArrayInt *)_num_faces)
    MEDFILESAFECALLERWR0(MEDmeshEntityNumberWr,(fid,maa.c_str(),_iteration,_order,MED_CELL,geoTypeReq2,_num_faces->getNumberOfTuples(),_num_faces->getConstPointer()));
  if((const DataArrayInt *)_num_nodes)
    MEDFILESAFECALLERWR0(MEDmeshEntityNumberWr,(fid,maa.c_str(),_iteration,_order,MED_NODE,MED_NONE,_num_nodes->getNumberOfTuples(),_num_nodes->getConstPointer()));
  if((const DataArrayAsciiChar *)_names_cells)
    {
      if(_names_cells->getNumberOfComponents()!=MED_SNAME_SIZE)
        {
          std::ostringstream oss; oss << "MEDFileStructuredMesh::writeStructuredLL : expected a name field on cells with number of components set to " << MED_SNAME_SIZE;
          oss << " ! The array has " << _names_cells->getNumberOfComponents() << " components !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      MEDFILESAFECALLERWR0(MEDmeshEntityNameWr,(fid,maa.c_str(),_iteration,_order,MED_CELL,geoTypeReq,_names_cells->getNumberOfTuples(),_names_cells->getConstPointer()));
    }
  if((const DataArrayAsciiChar *)_names_faces)
    {
      if(_names_faces->getNumberOfComponents()!=MED_SNAME_SIZE)
        {
          std::ostringstream oss; oss << "MEDFileStructuredMesh::writeStructuredLL : expected a name field on faces with number of components set to " << MED_SNAME_SIZE;
          oss << " ! The array has " << _names_faces->getNumberOfComponents() << " components !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      MEDFILESAFECALLERWR0(MEDmeshEntityNameWr,(fid,maa.c_str(),_iteration,_order,MED_CELL,geoTypeReq2,_names_faces->getNumberOfTuples(),_names_faces->getConstPointer()));
    }
  if((const DataArrayAsciiChar *)_names_nodes)
    {
      if(_names_nodes->getNumberOfComponents()!=MED_SNAME_SIZE)
        {
          std::ostringstream oss; oss << "MEDFileStructuredMesh::writeStructuredLL : expected a name field on nodes with number of components set to " << MED_SNAME_SIZE;
          oss << " ! The array has " << _names_cells->getNumberOfComponents() << " components !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      MEDFILESAFECALLERWR0(MEDmeshEntityNameWr,(fid,maa.c_str(),_iteration,_order,MED_NODE,MED_NONE,_names_nodes->getNumberOfTuples(),_names_nodes->getConstPointer()));
    }
  //
  MEDFileUMeshL2::WriteFamiliesAndGrps(fid,maa.c_str(),_families,_groups,_too_long_str);
}

/*!
 * Returns an empty instance of MEDFileCMesh.
 *  \return MEDFileCMesh * - a new instance of MEDFileCMesh. The caller is to delete this
 *          mesh using decrRef() as it is no more needed. 
 */
MEDFileCMesh *MEDFileCMesh::New()
{
  return new MEDFileCMesh;
}

/*!
 * Returns a new MEDFileCMesh holding the mesh data that has been read from a given MED
 * file. The first mesh in the file is loaded.
 *  \param [in] fileName - the name of MED file to read.
 *  \return MEDFileCMesh * - a new instance of MEDFileCMesh. The caller is to delete this
 *          mesh using decrRef() as it is no more needed. 
 *  \throw If the file is not readable.
 *  \throw If there is no meshes in the file.
 *  \throw If the mesh in the file is not a Cartesian one.
 */
MEDFileCMesh *MEDFileCMesh::New(const std::string& fileName, MEDFileMeshReadSelector *mrs)
{
  std::vector<std::string> ms=MEDLoader::GetMeshNames(fileName);
  if(ms.empty())
    {
      std::ostringstream oss; oss << "MEDFileUMesh::New : no meshes in file \"" << fileName << "\" !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName.c_str(),MED_ACC_RDONLY);
  int dt,it;
  ParaMEDMEM::MEDCouplingMeshType meshType;
  std::string dummy2;
  MEDFileMeshL2::GetMeshIdFromName(fid,ms.front(),meshType,dt,it,dummy2);
  return new MEDFileCMesh(fid,ms.front(),dt,it,mrs);
}

/*!
 * Returns a new MEDFileCMesh holding the mesh data that has been read from a given MED
 * file. The mesh to load is specified by its name and numbers of a time step and an
 * iteration.
 *  \param [in] fileName - the name of MED file to read.
 *  \param [in] mName - the name of the mesh to read.
 *  \param [in] dt - the number of a time step.
 *  \param [in] it - the number of an iteration.
 *  \return MEDFileCMesh * - a new instance of MEDFileCMesh. The caller is to delete this
 *          mesh using decrRef() as it is no more needed. 
 *  \throw If the file is not readable.
 *  \throw If there is no mesh with given attributes in the file.
 *  \throw If the mesh in the file is not a Cartesian one.
 */
MEDFileCMesh *MEDFileCMesh::New(const std::string& fileName, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName.c_str(),MED_ACC_RDONLY);
  return new MEDFileCMesh(fid,mName,dt,it,mrs);
}

std::size_t MEDFileCMesh::getHeapMemorySizeWithoutChildren() const
{
  return MEDFileStructuredMesh::getHeapMemorySizeWithoutChildren();
}

std::vector<const BigMemoryObject *> MEDFileCMesh::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret(MEDFileStructuredMesh::getDirectChildrenWithNull());
  ret.push_back((const MEDCouplingCMesh *)_cmesh);
  return ret;
}

/*!
 * Returns the dimension on cells in \a this mesh.
 *  \return int - the mesh dimension.
 *  \throw If there are no cells in this mesh.
 */
int MEDFileCMesh::getMeshDimension() const
{
  if(!((const MEDCouplingCMesh*)_cmesh))
    throw INTERP_KERNEL::Exception("MEDFileCMesh::getMeshDimension : unable to get meshdimension because no mesh set !");
  return _cmesh->getMeshDimension();
}

/*!
 * Returns the dimension on nodes in \a this mesh.
 *  \return int - the space dimension.
 *  \throw If there are no cells in this mesh.
 */
int MEDFileCMesh::getSpaceDimension() const
{
  if(!((const MEDCouplingCMesh*)_cmesh))
    throw INTERP_KERNEL::Exception("MEDFileCMesh::getSpaceDimension : unable to get spacedimension because no mesh set !");
  return _cmesh->getSpaceDimension();
}

/*!
 * Returns a string describing \a this mesh.
 *  \return std::string - the mesh information string.
 */
std::string MEDFileCMesh::simpleRepr() const
{
  return MEDFileStructuredMesh::simpleRepr();
}

/*!
 * Returns a full textual description of \a this mesh.
 *  \return std::string - the string holding the mesh description.
 */
std::string MEDFileCMesh::advancedRepr() const
{
  return simpleRepr();
}

MEDFileMesh *MEDFileCMesh::shallowCpy() const
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileCMesh> ret=new MEDFileCMesh(*this);
  return ret.retn();
}

MEDFileMesh *MEDFileCMesh::createNewEmpty() const
{
  return new MEDFileCMesh;
}

MEDFileMesh *MEDFileCMesh::deepCpy() const
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileCMesh> ret=new MEDFileCMesh(*this);
  if((const MEDCouplingCMesh*)_cmesh)
    ret->_cmesh=static_cast<MEDCouplingCMesh*>(_cmesh->deepCpy());
  ret->deepCpyAttributes();
  return ret.retn();
}

/*!
 * Checks if \a this and another mesh are equal.
 *  \param [in] other - the mesh to compare with.
 *  \param [in] eps - a precision used to compare real values.
 *  \param [in,out] what - the string returning description of unequal data.
 *  \return bool - \c true if the meshes are equal, \c false, else.
 */
bool MEDFileCMesh::isEqual(const MEDFileMesh *other, double eps, std::string& what) const
{
  if(!MEDFileStructuredMesh::isEqual(other,eps,what))
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
  return true;
}

/*!
 * Clears redundant attributes of incorporated data arrays.
 */
void MEDFileCMesh::clearNonDiscrAttributes() const
{
  MEDFileStructuredMesh::clearNonDiscrAttributes();
  MEDFileUMeshSplitL1::ClearNonDiscrAttributes(_cmesh);//to it is not a bug umeshsplit have already the method implemented
}

MEDFileCMesh::MEDFileCMesh()
{
}

MEDFileCMesh::MEDFileCMesh(med_idt fid, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs)
try
{
    loadCMeshFromFile(fid,mName,dt,it,mrs);
    loadJointsFromFile(fid);
}
catch(INTERP_KERNEL::Exception& e)
{
    throw e;
}

void MEDFileCMesh::loadCMeshFromFile(med_idt fid, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs)
{
  ParaMEDMEM::MEDCouplingMeshType meshType;
  int dummy0,dummy1;
  std::string dtunit;
  int mid=MEDFileMeshL2::GetMeshIdFromName(fid,mName,meshType,dummy0,dummy1,dtunit);
  if(meshType!=CARTESIAN)
    {
      std::ostringstream oss; oss << "Trying to load as cartesian an existing mesh with name '" << mName << "' that is NOT cartesian !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  MEDFileCMeshL2 loaderl2;
  loaderl2.loadAll(fid,mid,mName,dt,it);
  MEDCouplingCMesh *mesh=loaderl2.getMesh();
  mesh->incrRef();
  _cmesh=mesh;
  loadStrMeshFromFile(&loaderl2,fid,mName,dt,it,mrs);
}

/*!
 * Returns a const pointer to MEDCouplingCMesh held by \a this mesh.
 *  \return const MEDCouplingCMesh * - a pointer to the held MEDCouplingCMesh.
 */
const MEDCouplingCMesh *MEDFileCMesh::getMesh() const
{
  synchronizeTinyInfoOnLeaves();
  return _cmesh;
}

const MEDCouplingStructuredMesh *MEDFileCMesh::getStructuredMesh() const
{
  synchronizeTinyInfoOnLeaves();
  return _cmesh;
}

/*!
 * Sets the MEDCouplingCMesh holding the data of \a this mesh.
 *  \param [in] m - the new MEDCouplingCMesh to refer to.
 *  \throw If the name or the description of \a this mesh and \a m are not empty and are
 *         different. 
 */
void MEDFileCMesh::setMesh(MEDCouplingCMesh *m)
{
  dealWithTinyInfo(m);
  if(m)
    m->incrRef();
  _cmesh=m;
}

void MEDFileCMesh::writeLL(med_idt fid) const
{
  INTERP_KERNEL::AutoPtr<char> maa=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> desc=MEDLoaderBase::buildEmptyString(MED_COMMENT_SIZE);
  INTERP_KERNEL::AutoPtr<char> dtunit=MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE);
  MEDLoaderBase::safeStrCpy(_name.c_str(),MED_NAME_SIZE,maa,_too_long_str);
  MEDLoaderBase::safeStrCpy(_desc_name.c_str(),MED_COMMENT_SIZE,desc,_too_long_str);
  MEDLoaderBase::safeStrCpy(_dt_unit.c_str(),MED_LNAME_SIZE,dtunit,_too_long_str);
  int spaceDim(_cmesh->getSpaceDimension());
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
  MEDFILESAFECALLERWR0(MEDmeshCr,(fid,maa,spaceDim,spaceDim,MED_STRUCTURED_MESH,desc,dtunit,MED_SORT_DTIT,MED_CARTESIAN,comp,unit));
  if(_univ_wr_status)
    MEDFILESAFECALLERWR0(MEDmeshUniversalNameWr,(fid,maa));
  MEDFILESAFECALLERWR0(MEDmeshGridTypeWr,(fid,maa,MED_CARTESIAN_GRID));
  for(int i=0;i<spaceDim;i++)
    {
      const DataArrayDouble *da=_cmesh->getCoordsAt(i);
      MEDFILESAFECALLERWR0(MEDmeshGridIndexCoordinateWr,(fid,maa,_iteration,_order,_time,i+1,da->getNumberOfTuples(),da->getConstPointer()));
    }
  //
  std::string meshName(MEDLoaderBase::buildStringFromFortran(maa,MED_NAME_SIZE));
  MEDFileStructuredMesh::writeStructuredLL(fid,meshName);

  writeJoints(fid);
}

void MEDFileCMesh::synchronizeTinyInfoOnLeaves() const
{
  const MEDCouplingCMesh *cmesh=_cmesh;
  if(!cmesh)
    return;
  (const_cast<MEDCouplingCMesh *>(cmesh))->setName(_name);
  (const_cast<MEDCouplingCMesh *>(cmesh))->setDescription(_desc_name);
  (const_cast<MEDCouplingCMesh *>(cmesh))->setTime(_time,_iteration,_order);
  (const_cast<MEDCouplingCMesh *>(cmesh))->setTimeUnit(_dt_unit);
}

MEDFileCurveLinearMesh *MEDFileCurveLinearMesh::New()
{
  return new MEDFileCurveLinearMesh;
}

MEDFileCurveLinearMesh *MEDFileCurveLinearMesh::New(const std::string& fileName, MEDFileMeshReadSelector *mrs)
{
  std::vector<std::string> ms=MEDLoader::GetMeshNames(fileName);
  if(ms.empty())
    {
      std::ostringstream oss; oss << "MEDFileUMesh::New : no meshes in file \"" << fileName << "\" !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName.c_str(),MED_ACC_RDONLY);
  int dt,it;
  ParaMEDMEM::MEDCouplingMeshType meshType;
  std::string dummy2;
  MEDFileMeshL2::GetMeshIdFromName(fid,ms.front(),meshType,dt,it,dummy2);
  return new MEDFileCurveLinearMesh(fid,ms.front(),dt,it,mrs);
}

MEDFileCurveLinearMesh *MEDFileCurveLinearMesh::New(const std::string& fileName, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName.c_str(),MED_ACC_RDONLY);
  return new MEDFileCurveLinearMesh(fid,mName,dt,it,mrs);
}

std::size_t MEDFileCurveLinearMesh::getHeapMemorySizeWithoutChildren() const
{
  return MEDFileStructuredMesh::getHeapMemorySizeWithoutChildren();
}

std::vector<const BigMemoryObject *> MEDFileCurveLinearMesh::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret(MEDFileStructuredMesh::getDirectChildrenWithNull());
  ret.push_back((const MEDCouplingCurveLinearMesh *)_clmesh);
  return ret;
}

MEDFileMesh *MEDFileCurveLinearMesh::shallowCpy() const
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileCurveLinearMesh> ret=new MEDFileCurveLinearMesh(*this);
  return ret.retn();
}

MEDFileMesh *MEDFileCurveLinearMesh::createNewEmpty() const
{
  return new MEDFileCurveLinearMesh;
}

MEDFileMesh *MEDFileCurveLinearMesh::deepCpy() const
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileCurveLinearMesh> ret=new MEDFileCurveLinearMesh(*this);
  if((const MEDCouplingCurveLinearMesh*)_clmesh)
    ret->_clmesh=static_cast<MEDCouplingCurveLinearMesh*>(_clmesh->deepCpy());
  ret->deepCpyAttributes();
  return ret.retn();
}

int MEDFileCurveLinearMesh::getMeshDimension() const
{
  if(!((const MEDCouplingCurveLinearMesh*)_clmesh))
    throw INTERP_KERNEL::Exception("MEDFileCurveLinearMesh::getMeshDimension : unable to get meshdimension because no mesh set !");
  return _clmesh->getMeshDimension();
}

std::string MEDFileCurveLinearMesh::simpleRepr() const
{
  return MEDFileStructuredMesh::simpleRepr();
}

std::string MEDFileCurveLinearMesh::advancedRepr() const
{
  return simpleRepr();
}

bool MEDFileCurveLinearMesh::isEqual(const MEDFileMesh *other, double eps, std::string& what) const
{
  if(!MEDFileStructuredMesh::isEqual(other,eps,what))
    return false;
  const MEDFileCurveLinearMesh *otherC=dynamic_cast<const MEDFileCurveLinearMesh *>(other);
  if(!otherC)
    {
      what="Mesh types differ ! This is curve linear and other is NOT !";
      return false;
    }
  clearNonDiscrAttributes();
  otherC->clearNonDiscrAttributes();
  const MEDCouplingCurveLinearMesh *coo1=_clmesh;
  const MEDCouplingCurveLinearMesh *coo2=otherC->_clmesh;
  if((coo1==0 && coo2!=0) || (coo1!=0 && coo2==0))
    {
      what="Mismatch of curve linear meshes ! One is defined and not other !";
      return false;
    }
  if(coo1)
    {
      bool ret=coo1->isEqual(coo2,eps);
      if(!ret)
        {
          what="curve linear meshes differ !";
          return false;
        }
    }
  return true;
}

void MEDFileCurveLinearMesh::clearNonDiscrAttributes() const
{
  MEDFileStructuredMesh::clearNonDiscrAttributes();
  MEDFileUMeshSplitL1::ClearNonDiscrAttributes(_clmesh);//to it is not a bug umeshsplit have already the method implemented
}

void MEDFileCurveLinearMesh::synchronizeTinyInfoOnLeaves() const
{
  const MEDCouplingCurveLinearMesh *clmesh=_clmesh;
  if(!clmesh)
    return;
  (const_cast<MEDCouplingCurveLinearMesh *>(clmesh))->setName(_name);
  (const_cast<MEDCouplingCurveLinearMesh *>(clmesh))->setDescription(_desc_name);
  (const_cast<MEDCouplingCurveLinearMesh *>(clmesh))->setTime(_time,_iteration,_order);
  (const_cast<MEDCouplingCurveLinearMesh *>(clmesh))->setTimeUnit(_dt_unit);
}

const MEDCouplingCurveLinearMesh *MEDFileCurveLinearMesh::getMesh() const
{
  synchronizeTinyInfoOnLeaves();
  return _clmesh;
}

void MEDFileCurveLinearMesh::setMesh(MEDCouplingCurveLinearMesh *m)
{
  dealWithTinyInfo(m);
  if(m)
    m->incrRef();
  _clmesh=m;
}

const MEDCouplingStructuredMesh *MEDFileCurveLinearMesh::getStructuredMesh() const
{
  synchronizeTinyInfoOnLeaves();
  return _clmesh;
}

MEDFileCurveLinearMesh::MEDFileCurveLinearMesh()
{
}

MEDFileCurveLinearMesh::MEDFileCurveLinearMesh(med_idt fid, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs)
try
{
    loadCLMeshFromFile(fid,mName,dt,it,mrs);
    loadJointsFromFile(fid);
}
catch(INTERP_KERNEL::Exception& e)
{
    throw e;
}

void MEDFileCurveLinearMesh::writeLL(med_idt fid) const
{
  INTERP_KERNEL::AutoPtr<char> maa=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> desc=MEDLoaderBase::buildEmptyString(MED_COMMENT_SIZE);
  INTERP_KERNEL::AutoPtr<char> dtunit=MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE);
  MEDLoaderBase::safeStrCpy(_name.c_str(),MED_NAME_SIZE,maa,_too_long_str);
  MEDLoaderBase::safeStrCpy(_desc_name.c_str(),MED_COMMENT_SIZE,desc,_too_long_str);
  MEDLoaderBase::safeStrCpy(_dt_unit.c_str(),MED_LNAME_SIZE,dtunit,_too_long_str);
  int spaceDim=_clmesh->getSpaceDimension();
  int meshDim=_clmesh->getMeshDimension();
  INTERP_KERNEL::AutoPtr<char> comp=MEDLoaderBase::buildEmptyString(spaceDim*MED_SNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> unit=MEDLoaderBase::buildEmptyString(spaceDim*MED_SNAME_SIZE);
  const DataArrayDouble *coords=_clmesh->getCoords();
  if(!coords)
    throw INTERP_KERNEL::Exception("MEDFileCurveLinearMesh::writeLL : no coordinates set !");
  for(int i=0;i<spaceDim;i++)
    {
      std::string info(_clmesh->getCoords()->getInfoOnComponent(i));
      std::string c,u;
      MEDLoaderBase::splitIntoNameAndUnit(info,c,u);
      MEDLoaderBase::safeStrCpy2(c.c_str(),MED_SNAME_SIZE-1,comp+i*MED_SNAME_SIZE,_too_long_str);//MED_TAILLE_PNOM-1 to avoid to write '\0' on next compo
      MEDLoaderBase::safeStrCpy2(u.c_str(),MED_SNAME_SIZE-1,unit+i*MED_SNAME_SIZE,_too_long_str);//MED_TAILLE_PNOM-1 to avoid to write '\0' on next compo
    }
  MEDFILESAFECALLERWR0(MEDmeshCr,(fid,maa,spaceDim,meshDim,MED_STRUCTURED_MESH,desc,dtunit,MED_SORT_DTIT,MED_CARTESIAN,comp,unit));
  if(_univ_wr_status)
    MEDFILESAFECALLERWR0(MEDmeshUniversalNameWr,(fid,maa));
  MEDFILESAFECALLERWR0(MEDmeshGridTypeWr,(fid,maa,MED_CURVILINEAR_GRID));
  std::vector<int> nodeGridSt=_clmesh->getNodeGridStructure();
  MEDFILESAFECALLERWR0(MEDmeshGridStructWr,(fid,maa,_iteration,_order,_time,&nodeGridSt[0]));

  MEDFILESAFECALLERWR0(MEDmeshNodeCoordinateWr,(fid,maa,_iteration,_order,_time,MED_FULL_INTERLACE,coords->getNumberOfTuples(),coords->begin()));
  //
  std::string meshName(MEDLoaderBase::buildStringFromFortran(maa,MED_NAME_SIZE));
  MEDFileStructuredMesh::writeStructuredLL(fid,meshName);

  writeJoints(fid);
}

void MEDFileCurveLinearMesh::loadCLMeshFromFile(med_idt fid, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs)
{
  ParaMEDMEM::MEDCouplingMeshType meshType;
  int dummy0,dummy1;
  std::string dtunit;
  int mid=MEDFileMeshL2::GetMeshIdFromName(fid,mName,meshType,dummy0,dummy1,dtunit);
  if(meshType!=CURVE_LINEAR)
    {
      std::ostringstream oss; oss << "Trying to load as curve linear an existing mesh with name '" << mName << "' that is NOT curve linear !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  MEDFileCLMeshL2 loaderl2;
  loaderl2.loadAll(fid,mid,mName,dt,it);
  MEDCouplingCurveLinearMesh *mesh=loaderl2.getMesh();
  mesh->incrRef();
  _clmesh=mesh;
  loadStrMeshFromFile(&loaderl2,fid,mName,dt,it,mrs);
}

MEDFileMeshMultiTS *MEDFileMeshMultiTS::New()
{
  return new MEDFileMeshMultiTS;
}

MEDFileMeshMultiTS *MEDFileMeshMultiTS::New(const std::string& fileName)
{
  return new MEDFileMeshMultiTS(fileName);
}

MEDFileMeshMultiTS *MEDFileMeshMultiTS::New(const std::string& fileName, const std::string& mName)
{
  return new MEDFileMeshMultiTS(fileName,mName);
}

MEDFileMeshMultiTS *MEDFileMeshMultiTS::deepCpy() const
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileMeshMultiTS> ret=MEDFileMeshMultiTS::New();
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> > meshOneTs(_mesh_one_ts.size());
  std::size_t i=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> >::const_iterator it=_mesh_one_ts.begin();it!=_mesh_one_ts.end();it++,i++)
    if((const MEDFileMesh *)*it)
      meshOneTs[i]=(*it)->deepCpy();
  ret->_mesh_one_ts=meshOneTs;
  return ret.retn();
}

std::size_t MEDFileMeshMultiTS::getHeapMemorySizeWithoutChildren() const
{
  return _mesh_one_ts.capacity()*sizeof(MEDCouplingAutoRefCountObjectPtr<MEDFileMesh>);
}

std::vector<const BigMemoryObject *> MEDFileMeshMultiTS::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> >::const_iterator it=_mesh_one_ts.begin();it!=_mesh_one_ts.end();it++)
    ret.push_back((const MEDFileMesh *)*it);
  return ret;
}

std::string MEDFileMeshMultiTS::getName() const
{
  if(_mesh_one_ts.empty())
    throw INTERP_KERNEL::Exception("MEDFileMeshMultiTS::getName : no time steps set !");
  return _mesh_one_ts[0]->getName();
}

void MEDFileMeshMultiTS::setName(const std::string& newMeshName)
{
  std::string oldName(getName());
  std::vector< std::pair<std::string,std::string> > v(1);
  v[0].first=oldName; v[0].second=newMeshName;
  changeNames(v);
}

bool MEDFileMeshMultiTS::changeNames(const std::vector< std::pair<std::string,std::string> >& modifTab)
{
  bool ret=false;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> >::iterator it=_mesh_one_ts.begin();it!=_mesh_one_ts.end();it++)
    {
      MEDFileMesh *cur(*it);
      if(cur)
        ret=cur->changeNames(modifTab) || ret;
    }
  return ret;
}

MEDFileMesh *MEDFileMeshMultiTS::getOneTimeStep() const
{
  if(_mesh_one_ts.empty())
    throw INTERP_KERNEL::Exception("MEDFileMeshMultiTS::getOneTimeStep : empty time step set !");
  return const_cast<MEDFileMesh *>(static_cast<const MEDFileMesh *>(_mesh_one_ts[0]));
}

void MEDFileMeshMultiTS::setOneTimeStep(MEDFileMesh *mesh1TimeStep)
{
  if(!mesh1TimeStep)
    throw INTERP_KERNEL::Exception("MEDFileMeshMultiTS::setOneTimeStep : input pointer should be different from 0 !");
  _mesh_one_ts.resize(1);
  mesh1TimeStep->incrRef();
  //MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> toto=mesh1TimeStep;
  _mesh_one_ts[0]=mesh1TimeStep;
}

MEDFileJoints * MEDFileMeshMultiTS::getJoints() const
{
  if ( MEDFileMesh* m = getOneTimeStep() )
    return m->getJoints();
  return 0;
}

/*!
 * \brief Set Joints that are common to all time-stamps
 */
void MEDFileMeshMultiTS::setJoints( MEDFileJoints* joints )
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> >::iterator it=_mesh_one_ts.begin();it!=_mesh_one_ts.end();it++)
    {
      (*it)->setJoints( joints );
    }
}

void MEDFileMeshMultiTS::write(med_idt fid) const
{
  MEDFileJoints *joints(getJoints());
  bool jointsWritten(false);

  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> >::const_iterator it=_mesh_one_ts.begin();it!=_mesh_one_ts.end();it++)
    {
      if ( jointsWritten )
        const_cast<MEDFileMesh&>(**it).setJoints( 0 );
      else
        jointsWritten = true;

      (*it)->copyOptionsFrom(*this);
      (*it)->write(fid);
    }

  (const_cast<MEDFileMeshMultiTS*>(this))->setJoints( joints ); // restore joints
}

void MEDFileMeshMultiTS::write(const std::string& fileName, int mode) const
{
  med_access_mode medmod=MEDFileUtilities::TraduceWriteMode(mode);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName.c_str(),medmod);
  std::ostringstream oss; oss << "MEDFileMesh : error on attempt to write in file : \"" << fileName << "\""; 
  MEDFileUtilities::CheckMEDCode(fid,fid,oss.str());
  write(fid);
}

void MEDFileMeshMultiTS::loadFromFile(const std::string& fileName, const std::string& mName)
{
  MEDFileJoints* joints = 0;
  if ( !_mesh_one_ts.empty() && getOneTimeStep() )
    {
      // joints of mName already read, pass them to MEDFileMesh::New() to prevent repeated reading
      joints = getOneTimeStep()->getJoints();
    }

  _mesh_one_ts.clear();  //for the moment to be improved
  _mesh_one_ts.push_back( MEDFileMesh::New(fileName,mName,-1,-1,0, joints ));
}

MEDFileMeshMultiTS::MEDFileMeshMultiTS()
{
}

MEDFileMeshMultiTS::MEDFileMeshMultiTS(const std::string& fileName)
try
{
    std::vector<std::string> ms=MEDLoader::GetMeshNames(fileName);
    if(ms.empty())
      {
        std::ostringstream oss; oss << "MEDFileUMesh::New : no meshes in file \"" << fileName << "\" !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    MEDFileUtilities::CheckFileForRead(fileName);
    MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName.c_str(),MED_ACC_RDONLY);
    int dt,it;
    ParaMEDMEM::MEDCouplingMeshType meshType;
    std::string dummy2;
    MEDFileMeshL2::GetMeshIdFromName(fid,ms.front(),meshType,dt,it,dummy2);
    loadFromFile(fileName,ms.front());
}
catch(INTERP_KERNEL::Exception& e)
{
    throw e;
}

MEDFileMeshMultiTS::MEDFileMeshMultiTS(const std::string& fileName, const std::string& mName)
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

MEDFileMeshes *MEDFileMeshes::New(const std::string& fileName)
{
  return new MEDFileMeshes(fileName);
}

void MEDFileMeshes::write(med_idt fid) const
{
  checkCoherency();
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileMeshMultiTS> >::const_iterator it=_meshes.begin();it!=_meshes.end();it++)
    {
      (*it)->copyOptionsFrom(*this);
      (*it)->write(fid);
    }
}

void MEDFileMeshes::write(const std::string& fileName, int mode) const
{
  med_access_mode medmod=MEDFileUtilities::TraduceWriteMode(mode);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName.c_str(),medmod);
  std::ostringstream oss; oss << "MEDFileMesh : error on attempt to write in file : \"" << fileName << "\""; 
  MEDFileUtilities::CheckMEDCode(fid,fid,oss.str());
  checkCoherency();
  write(fid);
}

int MEDFileMeshes::getNumberOfMeshes() const
{
  return _meshes.size();
}

MEDFileMeshesIterator *MEDFileMeshes::iterator()
{
  return new MEDFileMeshesIterator(this);
}

/** Return a borrowed reference (caller is not responsible) */
MEDFileMesh *MEDFileMeshes::getMeshAtPos(int i) const
{
  if(i<0 || i>=(int)_meshes.size())
    {
      std::ostringstream oss; oss << "MEDFileMeshes::getMeshAtPos : invalid mesh id given in parameter ! Should be in [0;" << _meshes.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return _meshes[i]->getOneTimeStep();
}

/** Return a borrowed reference (caller is not responsible) */
MEDFileMesh *MEDFileMeshes::getMeshWithName(const std::string& mname) const
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

std::vector<std::string> MEDFileMeshes::getMeshesNames() const
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
/*const MEDFileJoints* MEDFileMeshes::getJoints() const
{
  const MEDFileJoints *ret=_joints;
  if(!ret)
  {
    std::ostringstream oss; oss << "MEDFileMeshes::getJoints : joints is not defined !";
    throw INTERP_KERNEL::Exception(oss.str().c_str());
  }
  return ret;
}*/

bool MEDFileMeshes::changeNames(const std::vector< std::pair<std::string,std::string> >& modifTab)
{
  bool ret=false;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileMeshMultiTS> >::iterator it=_meshes.begin();it!=_meshes.end();it++)
    {
      MEDFileMeshMultiTS *cur(*it);
      if(cur)
        ret=cur->changeNames(modifTab) || ret;
    }
  return ret;
}

void MEDFileMeshes::resize(int newSize)
{
  _meshes.resize(newSize);
}

void MEDFileMeshes::pushMesh(MEDFileMesh *mesh)
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDFileMeshes::pushMesh : invalid input pointer ! should be different from 0 !");
  MEDFileMeshMultiTS *elt=MEDFileMeshMultiTS::New();
  elt->setOneTimeStep(mesh);
  _meshes.push_back(elt);
}

void MEDFileMeshes::setMeshAtPos(int i, MEDFileMesh *mesh)
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDFileMeshes::setMeshAtPos : invalid input pointer ! should be different from 0 !");
  if(i>=(int)_meshes.size())
    _meshes.resize(i+1);
  MEDFileMeshMultiTS *elt=MEDFileMeshMultiTS::New();
  elt->setOneTimeStep(mesh);
  _meshes[i]=elt;
}

void MEDFileMeshes::destroyMeshAtPos(int i)
{
  if(i<0 || i>=(int)_meshes.size())
    {
      std::ostringstream oss; oss << "MEDFileMeshes::destroyMeshAtPos : Invalid given id in input (" << i << ") should be in [0," << _meshes.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  _meshes.erase(_meshes.begin()+i);
}

void MEDFileMeshes::loadFromFile(const std::string& fileName)
{
  std::vector<std::string> ms=MEDLoader::GetMeshNames(fileName);
  int i=0;
  _meshes.resize(ms.size());
  for(std::vector<std::string>::const_iterator it=ms.begin();it!=ms.end();it++,i++)
    _meshes[i]=MEDFileMeshMultiTS::New(fileName,(*it));
}

MEDFileMeshes::MEDFileMeshes()
{
}

MEDFileMeshes::MEDFileMeshes(const std::string& fileName)
try
{
    loadFromFile(fileName);
}
catch(INTERP_KERNEL::Exception& /*e*/)
{
}

MEDFileMeshes *MEDFileMeshes::deepCpy() const
{
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileMeshMultiTS> > meshes(_meshes.size());
  std::size_t i=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileMeshMultiTS> >::const_iterator it=_meshes.begin();it!=_meshes.end();it++,i++)
    if((const MEDFileMeshMultiTS *)*it)
      meshes[i]=(*it)->deepCpy();
  MEDCouplingAutoRefCountObjectPtr<MEDFileMeshes> ret=MEDFileMeshes::New();
  ret->_meshes=meshes;
  return ret.retn();
}

std::size_t MEDFileMeshes::getHeapMemorySizeWithoutChildren() const
{
  return _meshes.capacity()*(sizeof(MEDCouplingAutoRefCountObjectPtr<MEDFileMeshMultiTS>));
}

std::vector<const BigMemoryObject *> MEDFileMeshes::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileMeshMultiTS> >::const_iterator it=_meshes.begin();it!=_meshes.end();it++)
    ret.push_back((const MEDFileMeshMultiTS *)*it);
  return ret;
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

void MEDFileMeshes::checkCoherency() const
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

MEDFileMeshesIterator::MEDFileMeshesIterator(MEDFileMeshes *ms):_ms(ms),_iter_id(0),_nb_iter(0)
{
  if(ms)
    {
      ms->incrRef();
      _nb_iter=ms->getNumberOfMeshes();
    }
}

MEDFileMeshesIterator::~MEDFileMeshesIterator()
{
}

MEDFileMesh *MEDFileMeshesIterator::nextt()
{
  if(_iter_id<_nb_iter)
    {
      MEDFileMeshes *ms(_ms);
      if(ms)
        return ms->getMeshAtPos(_iter_id++);
      else
        return 0;
    }
  else
    return 0;
}
