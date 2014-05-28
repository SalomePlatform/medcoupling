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

void DataArrayDoubleCollection::CheckDiscriminantNames(const std::vector<std::string>& names)
{
  std::set<std::string> s(names.begin(),names.end());
  if(s.size()!=names.size())
    throw INTERP_KERNEL::Exception("DataArrayDoubleCollection::CheckDiscriminantNames : The names of fields must be different each other ! It is not the case !");
}

std::size_t MEDCouplingGridCollection::getHeapMemorySizeWithoutChildren() const
{
  std::size_t ret(sizeof(MEDCouplingGridCollection));
  ret+=_map_of_dadc.capacity()*sizeof(std::pair<MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> >);
  return ret;
}

std::vector<const BigMemoryObject *> MEDCouplingGridCollection::getDirectChildren() const
{
  std::vector<const BigMemoryObject *> ret;
  for(std::vector< std::pair<MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> > >::const_iterator it=_map_of_dadc.begin();it!=_map_of_dadc.end();it++)
    {
      const DataArrayDoubleCollection *col((*it).second);
      if(col)
        ret.push_back(col);
    }
  return ret;
}

MEDCouplingAMRAttribute *MEDCouplingAMRAttribute::New(MEDCouplingCartesianAMRMesh *gf, const std::vector< std::pair<std::string,int> >& fieldNames)
{
  return new MEDCouplingAMRAttribute(gf,fieldNames);
}

void MEDCouplingAMRAttribute::alloc()
{
  _tlc.resetState();
}

void MEDCouplingAMRAttribute::dealloc()
{//tony
}

bool MEDCouplingAMRAttribute::changeGodFather(MEDCouplingCartesianAMRMesh *gf)
{
  bool ret(MEDCouplingDataForGodFather::changeGodFather(gf));
  return ret;
}

std::size_t MEDCouplingAMRAttribute::getHeapMemorySizeWithoutChildren() const
{
  std::size_t ret(sizeof(MEDCouplingAMRAttribute));
  return ret;
}

std::vector<const BigMemoryObject *> MEDCouplingAMRAttribute::getDirectChildren() const
{//tony
  return std::vector<const BigMemoryObject *>();
}

void MEDCouplingAMRAttribute::updateTime() const
{//tony
}

MEDCouplingAMRAttribute::MEDCouplingAMRAttribute(MEDCouplingCartesianAMRMesh *gf, const std::vector< std::pair<std::string,int> >& fieldNames):MEDCouplingDataForGodFather(gf)
{
}
