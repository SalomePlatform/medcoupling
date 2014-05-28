// Copyright (C) 2007-2014  CEA/DEN, EDF R&D
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
// Author : Anthony Geay

#include "MEDCouplingAMRAttribute.hxx"
#include "MEDCouplingMemArray.hxx"

using namespace ParaMEDMEM;

DataArrayDoubleCollection *DataArrayDoubleCollection::New(const std::vector< std::pair<std::string,int> >& fieldNames)
{
  return new DataArrayDoubleCollection(fieldNames);
}

void DataArrayDoubleCollection::allocTuples(int nbOfTuples)
{
  std::size_t sz(_arrs.size());
  for(std::size_t i=0;i<sz;i++)
    _arrs[i]->reAlloc(nbOfTuples);
}

void DataArrayDoubleCollection::dellocTuples()
{
  std::size_t sz(_arrs.size());
  for(std::size_t i=0;i<sz;i++)
    _arrs[i]->reAlloc(0);
}

void DataArrayDoubleCollection::spillInfoOnComponents(const std::vector< std::vector<std::string> >& compNames)
{
  std::size_t sz(_arrs.size());
  if(sz!=compNames.size())
    throw INTERP_KERNEL::Exception("DataArrayDoubleCollection::spillInfoOnComponents : first size of compNames has to be equal to the number of fields defined !");
  for(std::size_t i=0;i<sz;i++)
    {
      const std::vector<std::string>& names(compNames[i]);
      _arrs[i]->setInfoOnComponents(names);
    }
}

DataArrayDoubleCollection::DataArrayDoubleCollection(const std::vector< std::pair<std::string,int> >& fieldNames):_arrs(fieldNames.size())
{
  std::size_t sz(fieldNames.size());
  std::vector<std::string> names(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      const std::pair<std::string,int>& info(fieldNames[i]);
      _arrs[i]=DataArrayDouble::New();
      _arrs[i]->alloc(0,info.second);
      _arrs[i]->setName(info.first);
      names[i]=info.second;
    }
  CheckDiscriminantNames(names);
}

std::size_t DataArrayDoubleCollection::getHeapMemorySizeWithoutChildren() const
{
  std::size_t ret(sizeof(DataArrayDoubleCollection));
  ret+=_arrs.capacity()*sizeof(MEDCouplingAutoRefCountObjectPtr<DataArrayDouble>);
  return ret;
}

std::vector<const BigMemoryObject *> DataArrayDoubleCollection::getDirectChildren() const
{
  std::vector<const BigMemoryObject *> ret;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> >::const_iterator it=_arrs.begin();it!=_arrs.end();it++)
    {
      const DataArrayDouble *pt(*it);
      if(pt)
        ret.push_back(pt);
    }
  return ret;
}

void DataArrayDoubleCollection::updateTime() const
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> >::const_iterator it=_arrs.begin();it!=_arrs.end();it++)
    {
      const DataArrayDouble *pt(*it);
      if(pt)
        updateTimeWith(*pt);
    }
}

void DataArrayDoubleCollection::CheckDiscriminantNames(const std::vector<std::string>& names)
{
  std::set<std::string> s(names.begin(),names.end());
  if(s.size()!=names.size())
    throw INTERP_KERNEL::Exception("DataArrayDoubleCollection::CheckDiscriminantNames : The names of fields must be different each other ! It is not the case !");
}

MEDCouplingGridCollection *MEDCouplingGridCollection::New(const std::vector<const MEDCouplingCartesianAMRMeshGen *>& ms, const std::vector< std::pair<std::string,int> >& fieldNames)
{
  return new MEDCouplingGridCollection(ms,fieldNames);
}

void MEDCouplingGridCollection::alloc(int ghostLev)
{
  for(std::vector< std::pair<const MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> > >::iterator it=_map_of_dadc.begin();it!=_map_of_dadc.end();it++)
    {
      int nbTuples((*it).first->getNumberOfCellsAtCurrentLevelGhost(ghostLev));
      DataArrayDoubleCollection *dadc((*it).second);
      if(dadc)
        dadc->allocTuples(nbTuples);
      else
        throw INTERP_KERNEL::Exception("MEDCouplingGridCollection::alloc : internal error !");
    }
}

void MEDCouplingGridCollection::dealloc()
{
  for(std::vector< std::pair<const MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> > >::iterator it=_map_of_dadc.begin();it!=_map_of_dadc.end();it++)
    {
      DataArrayDoubleCollection *dadc((*it).second);
      if(dadc)
        dadc->dellocTuples();
      else
        throw INTERP_KERNEL::Exception("MEDCouplingGridCollection::dealloc : internal error !");
    }
}

void MEDCouplingGridCollection::spillInfoOnComponents(const std::vector< std::vector<std::string> >& compNames)
{
  for(std::vector< std::pair<const MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> > >::iterator it=_map_of_dadc.begin();it!=_map_of_dadc.end();it++)
    (*it).second->spillInfoOnComponents(compNames);
}

MEDCouplingGridCollection::MEDCouplingGridCollection(const std::vector<const MEDCouplingCartesianAMRMeshGen *>& ms, const std::vector< std::pair<std::string,int> >& fieldNames):_map_of_dadc(ms.size())
{
  std::size_t sz(ms.size());
  for(std::size_t i=0;i<sz;i++)
    {
      if(!ms[i])
        throw INTERP_KERNEL::Exception("MEDCouplingGridCollection constructor : presence of NULL MEDCouplingCartesianAMRMeshGen instance !");
      _map_of_dadc[i].first=ms[i];
      _map_of_dadc[i].second=DataArrayDoubleCollection::New(fieldNames);
    }
}

std::size_t MEDCouplingGridCollection::getHeapMemorySizeWithoutChildren() const
{
  std::size_t ret(sizeof(MEDCouplingGridCollection));
  ret+=_map_of_dadc.capacity()*sizeof(std::pair<const MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> >);
  return ret;
}

std::vector<const BigMemoryObject *> MEDCouplingGridCollection::getDirectChildren() const
{
  std::vector<const BigMemoryObject *> ret;
  for(std::vector< std::pair<const MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> > >::const_iterator it=_map_of_dadc.begin();it!=_map_of_dadc.end();it++)
    {
      const DataArrayDoubleCollection *col((*it).second);
      if(col)
        ret.push_back(col);
    }
  return ret;
}

void MEDCouplingGridCollection::updateTime() const
{
  for(std::vector< std::pair<const MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> > >::const_iterator it=_map_of_dadc.begin();it!=_map_of_dadc.end();it++)
    {
      const MEDCouplingCartesianAMRMeshGen *a((*it).first);
      if(a)
        updateTimeWith(*a);
      const DataArrayDoubleCollection *b((*it).second);
      if(b)
        updateTimeWith(*b);
    }
}

/*!
 * This method creates, attach to a main AMR mesh \a gf ( called god father :-) ) and returns a data linked to \a gf ready for the computation.
 */
MEDCouplingAMRAttribute *MEDCouplingAMRAttribute::New(MEDCouplingCartesianAMRMesh *gf, const std::vector< std::pair<std::string,int> >& fieldNames)
{
  return new MEDCouplingAMRAttribute(gf,fieldNames);
}

/*!
 * Assign the info on components for all DataArrayDouble instance recursively stored in \a this.
 * The first dim of input \a compNames is the field id in the same order than those implicitely specified in \a fieldNames parameter of MEDCouplingAMRAttribute::New.
 * The second dim of \a compNames represent the component names component per component corresponding to the field. The size of this 2nd dimension has
 * to perfectly fit with those specified in MEDCouplingAMRAttribute::New.
 */
void MEDCouplingAMRAttribute::spillInfoOnComponents(const std::vector< std::vector<std::string> >& compNames)
{
  _tlc.checkConst();
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingGridCollection> >::iterator it=_levs.begin();it!=_levs.end();it++)
    (*it)->spillInfoOnComponents(compNames);
}

/*!
 * This method allocates all DataArrayDouble instances stored recursively in \a this.
 *
 * \param [in] ghostLev - The size of ghost zone.
 *
 * \sa dealloc
 */
void MEDCouplingAMRAttribute::alloc(int ghostLev)
{
  _tlc.resetState();
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingGridCollection> >::iterator it=_levs.begin();it!=_levs.end();it++)
    {
      MEDCouplingGridCollection *elt(*it);
      if(elt)
        elt->alloc(ghostLev);
      else
        throw INTERP_KERNEL::Exception("MEDCouplingAMRAttribute::alloc : internal error !");
    }
}

/*!
 * This method deallocates all DataArrayDouble instances stored recursively in \a this.
 * \sa alloc
 */
void MEDCouplingAMRAttribute::dealloc()
{
  _tlc.checkConst();
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingGridCollection> >::iterator it=_levs.begin();it!=_levs.end();it++)
    {
      MEDCouplingGridCollection *elt(*it);
      if(elt)
        elt->dealloc();
      else
        throw INTERP_KERNEL::Exception("MEDCouplingAMRAttribute::dealloc : internal error !");
    }
}

bool MEDCouplingAMRAttribute::changeGodFather(MEDCouplingCartesianAMRMesh *gf)
{
  bool ret(MEDCouplingDataForGodFather::changeGodFather(gf));
  return ret;
}

std::size_t MEDCouplingAMRAttribute::getHeapMemorySizeWithoutChildren() const
{
  std::size_t ret(sizeof(MEDCouplingAMRAttribute));
  ret+=_levs.capacity()*sizeof(MEDCouplingAutoRefCountObjectPtr<MEDCouplingGridCollection>);
  return ret;
}

std::vector<const BigMemoryObject *> MEDCouplingAMRAttribute::getDirectChildren() const
{
  std::vector<const BigMemoryObject *> ret;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingGridCollection> >::const_iterator it=_levs.begin();it!=_levs.end();it++)
    {
      const MEDCouplingGridCollection *elt(*it);
      if(elt)
        ret.push_back(elt);
    }
  return ret;
}

void MEDCouplingAMRAttribute::updateTime() const
{//tony
}

MEDCouplingAMRAttribute::MEDCouplingAMRAttribute(MEDCouplingCartesianAMRMesh *gf, const std::vector< std::pair<std::string,int> >& fieldNames):MEDCouplingDataForGodFather(gf)
{
  //gf non empty, checked by constructor
  int maxLev(gf->getMaxNumberOfLevelsRelativeToThis()+1);
  _levs.resize(maxLev+1);
  for(int i=0;i<maxLev;i++)
    {
      std::vector<MEDCouplingCartesianAMRPatchGen *> patches(gf->retrieveGridsAt(i));
      std::size_t sz(patches.size());
      std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatchGen> > patchesSafe(patches.size());
      for(std::size_t j=0;j<sz;j++)
        patchesSafe[j]=patches[j];
      std::vector<const MEDCouplingCartesianAMRMeshGen *> ms(sz);
      for(std::size_t j=0;j<sz;j++)
        {
          ms[j]=patches[j]->getMesh();
        }
      _levs[i]=MEDCouplingGridCollection::New(ms,fieldNames);
    }
}
