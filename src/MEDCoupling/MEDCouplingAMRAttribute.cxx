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
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingIMesh.hxx"

#include <sstream>

using namespace ParaMEDMEM;

DataArrayDoubleCollection *DataArrayDoubleCollection::New(const std::vector< std::pair<std::string,int> >& fieldNames)
{
  return new DataArrayDoubleCollection(fieldNames);
}

void DataArrayDoubleCollection::allocTuples(int nbOfTuples)
{
  std::size_t sz(_arrs.size());
  for(std::size_t i=0;i<sz;i++)
    _arrs[i].first->reAlloc(nbOfTuples);
}

void DataArrayDoubleCollection::dellocTuples()
{
  std::size_t sz(_arrs.size());
  for(std::size_t i=0;i<sz;i++)
    _arrs[i].first->reAlloc(0);
}

void DataArrayDoubleCollection::spillInfoOnComponents(const std::vector< std::vector<std::string> >& compNames)
{
  std::size_t sz(_arrs.size());
  if(sz!=compNames.size())
    throw INTERP_KERNEL::Exception("DataArrayDoubleCollection::spillInfoOnComponents : first size of compNames has to be equal to the number of fields defined !");
  for(std::size_t i=0;i<sz;i++)
    {
      const std::vector<std::string>& names(compNames[i]);
      _arrs[i].first->setInfoOnComponents(names);
    }
}

void DataArrayDoubleCollection::spillNatures(const std::vector<NatureOfField>& nfs)
{
  std::size_t sz(_arrs.size());
  if(sz!=nfs.size())
    throw INTERP_KERNEL::Exception("DataArrayDoubleCollection::spillNatures : first size of vector of NatureOfField has to be equal to the number of fields defined !");
  for(std::size_t i=0;i<sz;i++)
    _arrs[i].second=nfs[i];
}

std::vector<DataArrayDouble *> DataArrayDoubleCollection::retrieveFields() const
{
  std::size_t sz(_arrs.size());
  std::vector<DataArrayDouble *> ret(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      const DataArrayDouble *tmp(_arrs[i].first);
      ret[i]=const_cast<DataArrayDouble *>(tmp);
      if(ret[i])
        ret[i]->incrRef();
    }
  return ret;
}

const DataArrayDouble *DataArrayDoubleCollection::getFieldWithName(const std::string& name) const
{
  std::vector<std::string> vec;
  for(std::vector< std::pair< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble>, NatureOfField > >::const_iterator it=_arrs.begin();it!=_arrs.end();it++)
    {
      const DataArrayDouble *obj((*it).first);
      if(obj)
        {
          if(obj->getName()==name)
            return obj;
          else
            vec.push_back(obj->getName());
        }
    }
  std::ostringstream oss; oss << "DataArrayDoubleCollection::getFieldWithName : fieldName \"" << name << "\" does not exist in this ! Possibilities are :";
  std::copy(vec.begin(),vec.end(),std::ostream_iterator<std::string>(oss," "));
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

void DataArrayDoubleCollection::SynchronizeFineToCoarse(int ghostLev, const MEDCouplingCartesianAMRMeshGen *fatherOfFineMesh, int patchId, const DataArrayDoubleCollection *fine, DataArrayDoubleCollection *coarse)
{
  if(!fine || !coarse)
    throw INTERP_KERNEL::Exception("DataArrayDoubleCollection::SynchronizeFineToCoarse : the input DataArrayDouble collections must be non NULL !");
  std::size_t sz(coarse->_arrs.size());
  if(fine->_arrs.size()!=sz)
    throw INTERP_KERNEL::Exception("DataArrayDoubleCollection::SynchronizeFineToCoarse : the input DataArrayDouble collection must have the same size !");
  for(std::size_t i=0;i<sz;i++)
    {
      CheckSameNatures(fine->_arrs[i].second,coarse->_arrs[i].second);
      fatherOfFineMesh->fillCellFieldComingFromPatchGhost(patchId,fine->_arrs[i].first,coarse->_arrs[i].first,ghostLev,IsConservativeNature(coarse->_arrs[i].second));
    }
}

void DataArrayDoubleCollection::SynchronizeCoarseToFine(int ghostLev, const MEDCouplingCartesianAMRMeshGen *fatherOfFineMesh, int patchId, const DataArrayDoubleCollection *coarse, DataArrayDoubleCollection *fine)
{
  if(!fine || !coarse)
    throw INTERP_KERNEL::Exception("DataArrayDoubleCollection::SynchronizeCoarseToFine : the input DataArrayDouble collections must be non NULL !");
  std::size_t sz(coarse->_arrs.size());
  if(fine->_arrs.size()!=sz)
    throw INTERP_KERNEL::Exception("DataArrayDoubleCollection::SynchronizeCoarseToFine : the input DataArrayDouble collection must have the same size !");
  for(std::size_t i=0;i<sz;i++)
    {
      CheckSameNatures(fine->_arrs[i].second,coarse->_arrs[i].second);
      fatherOfFineMesh->fillCellFieldOnPatchGhost(patchId,coarse->_arrs[i].first,fine->_arrs[i].first,ghostLev,IsConservativeNature(coarse->_arrs[i].second));
    }
}

void DataArrayDoubleCollection::SynchronizeFineEachOther(int patchId, int ghostLev, const MEDCouplingCartesianAMRMeshGen *fatherOfFineMesh, const std::vector<const MEDCouplingCartesianAMRMeshGen *>& children, const std::vector<DataArrayDoubleCollection *>& fieldsOnFine)
{
  if(!fatherOfFineMesh)
    throw INTERP_KERNEL::Exception("DataArrayDoubleCollection::SynchronizeFineEachOther : father is NULL !");
  std::size_t sz(children.size());
  if(fieldsOnFine.size()!=sz)
    throw INTERP_KERNEL::Exception("DataArrayDoubleCollection::SynchronizeFineEachOther : sizes of vectors mismatch !");
  if(sz<=1)
    return ;
  std::size_t nbOfCall(fieldsOnFine[0]->_arrs.size());
  for(std::size_t i=0;i<sz;i++)
    if(fatherOfFineMesh->getPatchIdFromChildMesh(children[i])!=(int)i)
      throw INTERP_KERNEL::Exception("DataArrayDoubleCollection::SynchronizeFineEachOther : internal error !");
  for(std::size_t i=1;i<sz;i++)
    if(nbOfCall!=fieldsOnFine[i]->_arrs.size())
      throw INTERP_KERNEL::Exception("DataArrayDoubleCollection::SynchronizeFineEachOther : the collection of DataArrayDouble must have all the same size !");
  for(std::size_t i=0;i<nbOfCall;i++)
    {
      std::vector<const DataArrayDouble *> arrs(sz);
      for(std::size_t j=0;j<sz;j++)
        arrs[j]=fieldsOnFine[j]->_arrs[i].first;
      fatherOfFineMesh->fillCellFieldOnPatchOnlyGhostAdv(patchId,ghostLev,arrs);
    }
}

/*!
 * This method updates \a p1dac ghost zone parts using \a p2dac (which is really const). \a p2 is in the neighborhood of \a p1 (which size is defined by \a ghostLev).
 */
void DataArrayDoubleCollection::SynchronizeGhostZoneOfOneUsingTwo(int ghostLev, const MEDCouplingCartesianAMRPatch *p1, const DataArrayDoubleCollection *p1dac, const MEDCouplingCartesianAMRPatch *p2, const DataArrayDoubleCollection *p2dac)
{
  if(!p1 || !p1dac || !p2 || !p2dac)
    throw INTERP_KERNEL::Exception("DataArrayDoubleCollection::SynchronizeGhostZoneOfOneUsingTwo : input pointer must be not NULL !");
  std::size_t sz(p1dac->_arrs.size());
  if(p2dac->_arrs.size()!=sz)
    throw INTERP_KERNEL::Exception("DataArrayDoubleCollection::SynchronizeGhostZoneOfOneUsingTwo : size of DataArrayDouble Collection must be the same !");
  for(std::size_t i=0;i<sz;i++)
    {
      const DataArrayDouble *zeArrWhichGhostsWillBeUpdated(p1dac->_arrs[i].first);
      MEDCouplingCartesianAMRPatch::UpdateNeighborsOfOneWithTwoMixedLev(ghostLev,p1,p2,const_cast<DataArrayDouble *>(zeArrWhichGhostsWillBeUpdated),p2dac->_arrs[i].first);
    }
}

void DataArrayDoubleCollection::SynchronizeCoarseToFineOnlyInGhostZone(int ghostLev, const MEDCouplingCartesianAMRMeshGen *fatherOfFineMesh, int patchId, const DataArrayDoubleCollection *coarse, DataArrayDoubleCollection *fine)
{
  if(!fine || !coarse)
    throw INTERP_KERNEL::Exception("DataArrayDoubleCollection::SynchronizeCoarseToFineOnlyInGhostZone : the input DataArrayDouble collections must be non NULL !");
  std::size_t sz(coarse->_arrs.size());
  if(fine->_arrs.size()!=sz)
    throw INTERP_KERNEL::Exception("DataArrayDoubleCollection::SynchronizeCoarseToFineOnlyInGhostZone : the input DataArrayDouble collection must have the same size !");
  for(std::size_t i=0;i<sz;i++)
    fatherOfFineMesh->fillCellFieldOnPatchOnlyOnGhostZone(patchId,coarse->_arrs[i].first,fine->_arrs[i].first,ghostLev);
}

void DataArrayDoubleCollection::synchronizeMyGhostZoneUsing(int ghostLev, const DataArrayDoubleCollection& other, const MEDCouplingCartesianAMRPatch *thisp, const MEDCouplingCartesianAMRPatch *otherp, const MEDCouplingCartesianAMRMeshGen *father) const
{
  DataArrayDoubleCollection *thisNC(const_cast<DataArrayDoubleCollection *>(this));
  std::size_t sz(_arrs.size());
  if(other._arrs.size()!=sz)
    throw INTERP_KERNEL::Exception("DataArrayDoubleCollection::synchronizeMyGhostZoneUsing : sizes of collections must match !");
  for(std::size_t i=0;i<sz;i++)
    father->fillCellFieldOnPatchOnlyOnGhostZoneWith(ghostLev,thisp,otherp,thisNC->_arrs[i].first,other._arrs[i].first);
}

void DataArrayDoubleCollection::synchronizeMyGhostZoneUsingExt(int ghostLev, const DataArrayDoubleCollection& other, const MEDCouplingCartesianAMRPatch *thisp, const MEDCouplingCartesianAMRPatch *otherp) const
{
  DataArrayDoubleCollection *thisNC(const_cast<DataArrayDoubleCollection *>(this));
  std::size_t sz(_arrs.size());
  if(other._arrs.size()!=sz)
    throw INTERP_KERNEL::Exception("DataArrayDoubleCollection::synchronizeMyGhostZoneUsingExt : sizes of collections must match !");
  for(std::size_t i=0;i<sz;i++)
    MEDCouplingCartesianAMRPatch::UpdateNeighborsOfOneWithTwoExt(ghostLev,thisp,otherp,thisNC->_arrs[i].first,other._arrs[i].first);
}

DataArrayDoubleCollection::DataArrayDoubleCollection(const std::vector< std::pair<std::string,int> >& fieldNames):_arrs(fieldNames.size())
{
  std::size_t sz(fieldNames.size());
  std::vector<std::string> names(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      const std::pair<std::string,int>& info(fieldNames[i]);
      _arrs[i].first=DataArrayDouble::New();
      _arrs[i].first->alloc(0,info.second);
      _arrs[i].first->setName(info.first);
      names[i]=info.second;
      _arrs[i].second=ConservativeVolumic;
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
  for(std::vector< std::pair< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble>, NatureOfField > >::const_iterator it=_arrs.begin();it!=_arrs.end();it++)
    {
      const DataArrayDouble *pt((*it).first);
      if(pt)
        ret.push_back(pt);
    }
  return ret;
}

void DataArrayDoubleCollection::updateTime() const
{
  for(std::vector< std::pair< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble>, NatureOfField > >::const_iterator it=_arrs.begin();it!=_arrs.end();it++)
    {
      const DataArrayDouble *pt((*it).first);
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

bool DataArrayDoubleCollection::IsConservativeNature(NatureOfField n)
{
  CheckValidNature(n);
  return n==RevIntegral || n==IntegralGlobConstraint;
}

void DataArrayDoubleCollection::CheckSameNatures(NatureOfField n1, NatureOfField n2)
{
  CheckValidNature(n1);
  CheckValidNature(n2);
  if(n1!=n2)
    throw INTERP_KERNEL::Exception("DataArrayDoubleCollection::CheckSameNatures : natures are not the same !");
}

void DataArrayDoubleCollection::CheckValidNature(NatureOfField n)
{
  if(n!=ConservativeVolumic && n!=Integral && n!=IntegralGlobConstraint && n!=RevIntegral)
    throw INTERP_KERNEL::Exception("DataArrayDoubleCollection::CheckValidNature : unrecognized nature !");
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

void MEDCouplingGridCollection::spillNatures(const std::vector<NatureOfField>& nfs)
{
  for(std::vector< std::pair<const MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> > >::iterator it=_map_of_dadc.begin();it!=_map_of_dadc.end();it++)
    (*it).second->spillNatures(nfs);
}

bool MEDCouplingGridCollection::presenceOf(const MEDCouplingCartesianAMRMeshGen *m, int& pos) const
{
  int ret(0);
  for(std::vector< std::pair<const MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> > >::const_iterator it=_map_of_dadc.begin();it!=_map_of_dadc.end();it++,ret++)
    {
      if((*it).first==m)
        {
          pos=ret;
          return true;
        }
    }
  return false;
}

const DataArrayDoubleCollection& MEDCouplingGridCollection::getFieldsAt(int pos) const
{
  if(pos<0 || pos>(int)_map_of_dadc.size())
    throw INTERP_KERNEL::Exception("MEDCouplingGridCollection::getFieldsAt : invalid pos given in input ! Must be in [0,size) !");
  return *_map_of_dadc[pos].second;
}

void MEDCouplingGridCollection::SynchronizeFineToCoarse(int ghostLev, const MEDCouplingGridCollection *fine, const MEDCouplingGridCollection *coarse)
{
  if(!fine || !coarse)
    throw INTERP_KERNEL::Exception("MEDCouplingGridCollection::SynchronizeFineToCoarse : one or more input pointer is NULL !");
  const std::vector< std::pair<const MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> > >& mf(fine->_map_of_dadc);
  const std::vector< std::pair<const MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> > >& mc(coarse->_map_of_dadc);
  for(std::vector< std::pair<const MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> > >::const_iterator it=mf.begin();it!=mf.end();it++)
    {
      const MEDCouplingCartesianAMRMeshGen *fineMesh((*it).first);
      const MEDCouplingCartesianAMRMeshGen *fatherOfFineMesh(fineMesh->getFather());
      bool found(false);
      for(std::vector< std::pair<const MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> > >::const_iterator it0=mc.begin();it0!=mc.end() && !found;it0++)
        {
          if((*it0).first==fatherOfFineMesh)
            {
              found=true;
              int patchId(fatherOfFineMesh->getPatchIdFromChildMesh(fineMesh));
              const DataArrayDoubleCollection *coarseDaCol((*it0).second);
              DataArrayDoubleCollection *coarseModified(const_cast<DataArrayDoubleCollection *>(coarseDaCol));//coarse values in DataArrayDouble will be altered
              DataArrayDoubleCollection::SynchronizeFineToCoarse(ghostLev,fatherOfFineMesh,patchId,(*it).second,coarseModified);
            }
        }
      if(!found)
        throw INTERP_KERNEL::Exception("MEDCouplingGridCollection::SynchronizeFineToCoarse : a fine mesh is orphan regarding given coarse meshes !");
    }
}

void MEDCouplingGridCollection::SynchronizeCoarseToFine(int ghostLev, const MEDCouplingGridCollection *coarse, const MEDCouplingGridCollection *fine)
{
  if(!fine || !coarse)
    throw INTERP_KERNEL::Exception("MEDCouplingGridCollection::SynchronizeCoarseToFine : one or more input pointer is NULL !");
  const std::vector< std::pair<const MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> > >& mf(fine->_map_of_dadc);
  const std::vector< std::pair<const MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> > >& mc(coarse->_map_of_dadc);
  for(std::vector< std::pair<const MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> > >::const_iterator it=mf.begin();it!=mf.end();it++)
    {
      const MEDCouplingCartesianAMRMeshGen *fineMesh((*it).first);
      const MEDCouplingCartesianAMRMeshGen *fatherOfFineMesh(fineMesh->getFather());
      bool found(false);
      for(std::vector< std::pair<const MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> > >::const_iterator it0=mc.begin();it0!=mc.end() && !found;it0++)
        {
          if((*it0).first==fatherOfFineMesh)
            {
              found=true;
              int patchId(fatherOfFineMesh->getPatchIdFromChildMesh(fineMesh));
              const DataArrayDoubleCollection *fineDaCol((*it).second);
              DataArrayDoubleCollection *fineModified(const_cast<DataArrayDoubleCollection *>(fineDaCol));//fine values in DataArrayDouble will be altered
              DataArrayDoubleCollection::SynchronizeCoarseToFine(ghostLev,fatherOfFineMesh,patchId,(*it0).second,fineModified);
            }
        }
      if(!found)
        throw INTERP_KERNEL::Exception("MEDCouplingGridCollection::SynchronizeCoarseToFine : a fine mesh is orphan regarding given coarse meshes !");
    }
}

/*!
 * All the pairs in \a ps must share the same father. If not call synchronizeFineEachOtherExt instead.
 *
 * \sa synchronizeFineEachOtherExt
 */
void MEDCouplingGridCollection::synchronizeFineEachOther(int ghostLev, const std::vector< std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *> >& ps) const
{
  for(std::vector< std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *> >::const_iterator it=ps.begin();it!=ps.end();it++)
    {
      int p1,p2;
      if(!presenceOf((*it).first->getMesh(),p1))
        throw INTERP_KERNEL::Exception("MEDCouplingGridCollection::synchronizeFineEachOther : internal error #1 !");
      if(!presenceOf((*it).second->getMesh(),p2))
        throw INTERP_KERNEL::Exception("MEDCouplingGridCollection::synchronizeFineEachOther : internal error #2 !");
      const DataArrayDoubleCollection& col1(getFieldsAt(p1));
      const DataArrayDoubleCollection& col2(getFieldsAt(p2));
      col1.synchronizeMyGhostZoneUsing(ghostLev,col2,(*it).first,(*it).second,(*it).first->getMesh()->getFather());
    }
}

/*!
 * This method is a generalization of synchronizeFineEachOther because elements in pairs are \b not sharing the same father but are neighbors nevertheless.
 *
 * \sa synchronizeFineEachOther
 */
void MEDCouplingGridCollection::synchronizeFineEachOtherExt(int ghostLev, const std::vector< std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *> >& ps) const
{
  for(std::vector< std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *> >::const_iterator it=ps.begin();it!=ps.end();it++)
    {
      int p1,p2;
      if(!presenceOf((*it).first->getMesh(),p1))
        throw INTERP_KERNEL::Exception("MEDCouplingGridCollection::synchronizeFineEachOtherExt : internal error #1 !");
      if(!presenceOf((*it).second->getMesh(),p2))
        throw INTERP_KERNEL::Exception("MEDCouplingGridCollection::synchronizeFineEachOtherExt : internal error #2 !");
      const DataArrayDoubleCollection& col1(getFieldsAt(p1));
      const DataArrayDoubleCollection& col2(getFieldsAt(p2));
      col1.synchronizeMyGhostZoneUsingExt(ghostLev,col2,(*it).first,(*it).second);
    }
}

/*!
 * The pairs returned share the same direct father. The number of returned elements must be even.
 */
std::vector< std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *> > MEDCouplingGridCollection::findNeighbors(int ghostLev) const
{
  std::vector< std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *> > ret;
  std::map<const MEDCouplingCartesianAMRMeshGen *,std::vector< const MEDCouplingCartesianAMRMeshGen * > > m;
  for(std::vector< std::pair<const MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> > >::const_iterator it=_map_of_dadc.begin();it!=_map_of_dadc.end();it++)
    {
      const MEDCouplingCartesianAMRMeshGen *fineMesh((*it).first);
      const MEDCouplingCartesianAMRMeshGen *fatherOfFineMesh(fineMesh->getFather());
      m[fatherOfFineMesh].push_back(fineMesh);
    }
  for(std::map<const MEDCouplingCartesianAMRMeshGen *,std::vector< const MEDCouplingCartesianAMRMeshGen * > >::const_iterator it0=m.begin();it0!=m.end();it0++)
    {
      for(std::vector<const MEDCouplingCartesianAMRMeshGen *>::const_iterator it1=(*it0).second.begin();it1!=(*it0).second.end();it1++)
        {
          int patchId((*it0).first->getPatchIdFromChildMesh(*it1));
          std::vector<int> neighs((*it0).first->getPatchIdsInTheNeighborhoodOf(patchId,ghostLev));
          const MEDCouplingCartesianAMRPatch *pRef((*it0).first->getPatch(patchId));
          for(std::vector<int>::const_iterator it2=neighs.begin();it2!=neighs.end();it2++)
            {
              const MEDCouplingCartesianAMRPatch *pLoc((*it0).first->getPatch(*it2));
              ret.push_back(std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *>(pRef,pLoc));
            }
        }
    }
  if(ret.size()%2!=0)
    throw INTERP_KERNEL::Exception("MEDCouplingGridCollection::findNeighbors : something is wrong ! The number of neighbor pairs must be %2 ==0 !");
  return ret;
}

void MEDCouplingGridCollection::SynchronizeCoarseToFineOnlyInGhostZone(int ghostLev, const MEDCouplingGridCollection *coarse, const MEDCouplingGridCollection *fine)
{
  if(!fine || !coarse)
    throw INTERP_KERNEL::Exception("MEDCouplingGridCollection::SynchronizeCoarseToFineOnlyInGhostZone : one or more input pointer is NULL !");
  const std::vector< std::pair<const MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> > >& mf(fine->_map_of_dadc);
  const std::vector< std::pair<const MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> > >& mc(coarse->_map_of_dadc);
  for(std::vector< std::pair<const MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> > >::const_iterator it=mf.begin();it!=mf.end();it++)
    {
      const MEDCouplingCartesianAMRMeshGen *fineMesh((*it).first);
      const MEDCouplingCartesianAMRMeshGen *fatherOfFineMesh(fineMesh->getFather());
      bool found(false);
      for(std::vector< std::pair<const MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> > >::const_iterator it0=mc.begin();it0!=mc.end() && !found;it0++)
        {
          if((*it0).first==fatherOfFineMesh)
            {
              found=true;
              int patchId(fatherOfFineMesh->getPatchIdFromChildMesh(fineMesh));
              const DataArrayDoubleCollection *fineDaCol((*it).second);
              DataArrayDoubleCollection *fineModified(const_cast<DataArrayDoubleCollection *>(fineDaCol));//fine values in DataArrayDouble will be altered
              DataArrayDoubleCollection::SynchronizeCoarseToFineOnlyInGhostZone(ghostLev,fatherOfFineMesh,patchId,(*it0).second,fineModified);
            }
        }
      if(!found)
        throw INTERP_KERNEL::Exception("MEDCouplingGridCollection::SynchronizeCoarseToFineOnlyInGhostZone : a fine mesh is orphan regarding given coarse meshes !");
    }
}

void MEDCouplingGridCollection::fillIfInTheProgenyOf(const std::string& fieldName, const MEDCouplingCartesianAMRMeshGen *head, std::vector<const DataArrayDouble *>& recurseArrs) const
{
  for(std::vector< std::pair<const MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> > >::const_iterator it=_map_of_dadc.begin();it!=_map_of_dadc.end();it++)
    {
      const MEDCouplingCartesianAMRMeshGen *a((*it).first);
      if(head==a || head->isObjectInTheProgeny(a))
        {
          const DataArrayDoubleCollection *gc((*it).second);
          recurseArrs.push_back(gc->getFieldWithName(fieldName));
        }
    }
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

MEDCouplingCartesianAMRPatchGF::MEDCouplingCartesianAMRPatchGF(MEDCouplingCartesianAMRMesh *mesh):MEDCouplingCartesianAMRPatchGen(mesh)
{
}

std::size_t MEDCouplingCartesianAMRPatchGF::getHeapMemorySizeWithoutChildren() const
{
  return sizeof(MEDCouplingCartesianAMRPatchGF);
}

MEDCouplingDataForGodFather::MEDCouplingDataForGodFather(MEDCouplingCartesianAMRMeshGen *gf):_gf(gf),_tlc(gf)
{
  if(!gf)
    throw INTERP_KERNEL::Exception("MEDCouplingDataForGodFather constructor : A data has to be attached to a AMR Mesh instance !");
  gf->incrRef();
}

void MEDCouplingDataForGodFather::checkGodFatherFrozen() const
{
  _tlc.checkConst();
}

bool MEDCouplingDataForGodFather::changeGodFather(MEDCouplingCartesianAMRMeshGen *gf)
{
  bool ret(_tlc.keepTrackOfNewTL(gf));
  if(ret)
    {
      _gf=gf;
      if(gf)
        gf->incrRef();
    }
  return ret;
}

/*!
 * This method creates, attach to a main AMR mesh \a gf ( called god father :-) ) and returns a data linked to \a gf ready for the computation.
 */
MEDCouplingAMRAttribute *MEDCouplingAMRAttribute::New(MEDCouplingCartesianAMRMeshGen *gf, const std::vector< std::pair<std::string,int> >& fieldNames, int ghostLev)
{
  return new MEDCouplingAMRAttribute(gf,fieldNames,ghostLev);
}

MEDCouplingAMRAttribute *MEDCouplingAMRAttribute::New(MEDCouplingCartesianAMRMeshGen *gf, const std::vector< std::pair<std::string, std::vector<std::string> > >& fieldNames, int ghostLev)
{
  std::size_t sz(fieldNames.size());
  std::vector< std::pair<std::string,int> > fieldNames2(sz);
  std::vector< std::vector<std::string> > compNames(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      fieldNames2[i].first=fieldNames[i].first;
      fieldNames2[i].second=(int)fieldNames[i].second.size();
      compNames[i]=fieldNames[i].second;
    }
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingAMRAttribute> ret(New(gf,fieldNames2,ghostLev));
  ret->spillInfoOnComponents(compNames);
  return ret.retn();
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
 * Assign nature for each fields in \a this.
 * \param [in] nfs
 */
void MEDCouplingAMRAttribute::spillNatures(const std::vector<NatureOfField>& nfs)
{
  _tlc.checkConst();
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingGridCollection> >::iterator it=_levs.begin();it!=_levs.end();it++)
    (*it)->spillNatures(nfs);
}

/*!
 * This method returns all DataArrayDouble instances lying on the specified mesh \a mesh.
 * If \a mesh is not part of the progeny of god father object given at construction of \a this an exception will be thrown.
 *
 * \return std::vector<DataArrayDouble *> - DataArrayDouble instances to be deallocated by the caller (using decrRef).
 * \sa retrieveFieldOn
 */
std::vector<DataArrayDouble *> MEDCouplingAMRAttribute::retrieveFieldsOn(MEDCouplingCartesianAMRMeshGen *mesh) const
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingGridCollection> >::const_iterator it=_levs.begin();it!=_levs.end();it++)
    {
      int tmp(-1);
      if((*it)->presenceOf(mesh,tmp))
        {
          const DataArrayDoubleCollection& ddc((*it)->getFieldsAt(tmp));
          return ddc.retrieveFields();
        }
    }
  throw INTERP_KERNEL::Exception("MEDCouplingAMRAttribute::retrieveFieldsOn : the mesh specified is not in the progeny of this !");
}

/*!
 * \sa retrieveFieldsOn
 */
const DataArrayDouble *MEDCouplingAMRAttribute::getFieldOn(MEDCouplingCartesianAMRMeshGen *mesh, const std::string& fieldName) const
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingGridCollection> >::const_iterator it=_levs.begin();it!=_levs.end();it++)
    {
      int tmp(-1);
      if((*it)->presenceOf(mesh,tmp))
        {
          const DataArrayDoubleCollection& ddc((*it)->getFieldsAt(tmp));
          return ddc.getFieldWithName(fieldName);
        }
    }
  throw INTERP_KERNEL::Exception("MEDCouplingAMRAttribute::retrieveFieldOn : the mesh specified is not in the progeny of this !");
}

/*!
 * This method returns a field on an unstructured mesh the most refined as possible without overlap.
 * Ghost part are not visible here.
 *
 * \return MEDCouplingFieldDouble * - a field on cells that the caller has to deal with (deallocate it).
 */
MEDCouplingFieldDouble *MEDCouplingAMRAttribute::buildCellFieldOnRecurseWithoutOverlapWithoutGhost(MEDCouplingCartesianAMRMeshGen *mesh, const std::string& fieldName) const
{
  std::vector<const DataArrayDouble *> recurseArrs;
  std::size_t lev(0);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingGridCollection> >::const_iterator it=_levs.begin();it!=_levs.end();it++,lev++)
    {
      int tmp(-1);
      if((*it)->presenceOf(mesh,tmp))
        {
          const DataArrayDoubleCollection& ddc((*it)->getFieldsAt(tmp));
          recurseArrs.push_back(ddc.getFieldWithName(fieldName));
          break;
        }
    }
  lev++;
  for(std::size_t i=lev;i<_levs.size();i++)
    {
      const MEDCouplingGridCollection *gc(_levs[i]);
      gc->fillIfInTheProgenyOf(fieldName,mesh,recurseArrs);
    }
  return mesh->buildCellFieldOnRecurseWithoutOverlapWithoutGhost(_ghost_lev,recurseArrs);
}

/*!
 * This method builds a newly created field on cell just lying on mesh \a mesh without its eventual refinement.
 * The output field also displays ghost cells.
 *
 * \return MEDCouplingFieldDouble * - a field on cells that the caller has to deal with (deallocate it).
 *
 * \sa buildCellFieldOnWithoutGhost
 */
MEDCouplingFieldDouble *MEDCouplingAMRAttribute::buildCellFieldOnWithGhost(MEDCouplingCartesianAMRMeshGen *mesh, const std::string& fieldName) const
{
  const DataArrayDouble *arr(0);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingGridCollection> >::const_iterator it=_levs.begin();it!=_levs.end();it++)
    {
      int tmp(-1);
      if((*it)->presenceOf(mesh,tmp))
        {
          const DataArrayDoubleCollection& ddc((*it)->getFieldsAt(tmp));
          arr=ddc.getFieldWithName(fieldName);
        }
    }
  if(!arr)
    throw INTERP_KERNEL::Exception("MEDCouplingAMRAttribute::buildCellFieldOnWithGhost : the mesh specified is not in the progeny of this !");
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingIMesh> im(mesh->getImageMesh()->buildWithGhost(_ghost_lev));
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret(MEDCouplingFieldDouble::New(ON_CELLS));
  ret->setMesh(im);
  ret->setArray(const_cast<DataArrayDouble *>(arr));
  ret->setName(arr->getName());
  return ret.retn();
}

/*!
 * This method builds a newly created field on cell just lying on mesh \a mesh without its eventual refinement.
 * The output field does not display ghost cells.
 *
 * \return MEDCouplingFieldDouble * - a field on cells that the caller has to deal with (deallocate it).
 *
 * \sa buildCellFieldOnWithGhost
 */
MEDCouplingFieldDouble *MEDCouplingAMRAttribute::buildCellFieldOnWithoutGhost(MEDCouplingCartesianAMRMeshGen *mesh, const std::string& fieldName) const
{
  //tony
  const DataArrayDouble *arr(0);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingGridCollection> >::const_iterator it=_levs.begin();it!=_levs.end();it++)
    {
      int tmp(-1);
      if((*it)->presenceOf(mesh,tmp))
        {
          const DataArrayDoubleCollection& ddc((*it)->getFieldsAt(tmp));
          arr=ddc.getFieldWithName(fieldName);
        }
    }
  if(!arr)
    throw INTERP_KERNEL::Exception("MEDCouplingAMRAttribute::buildCellFieldOnWithoutGhost : the mesh specified is not in the progeny of this !");
  //
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingIMesh> im(mesh->getImageMesh()->buildWithGhost(_ghost_lev));
  std::vector<int> cgs(mesh->getImageMesh()->getCellGridStructure()),cgsWG(im->getCellGridStructure());
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr2(DataArrayDouble::New());
  arr2->alloc(mesh->getImageMesh()->getNumberOfCells(),arr->getNumberOfComponents());
  std::vector< std::pair<int,int> > cgs2(MEDCouplingStructuredMesh::GetCompactFrmtFromDimensions(cgs));
  MEDCouplingStructuredMesh::ApplyGhostOnCompactFrmt(cgs2,_ghost_lev);
  std::vector<int> fakeFactors(mesh->getImageMesh()->getSpaceDimension(),1);
  MEDCouplingIMesh::SpreadCoarseToFine(arr,cgsWG,arr2,cgs2,fakeFactors);
  arr2->copyStringInfoFrom(*arr);
  //
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret(MEDCouplingFieldDouble::New(ON_CELLS));
  ret->setMesh(mesh->getImageMesh());
  ret->setArray(arr2);
  ret->setName(arr->getName());
  return ret.retn();
}

/*!
 * This method synchronizes from fine to coarse direction arrays. This method makes the hypothesis that \a this has been allocated before using
 * MEDCouplingAMRAttribute::alloc method.
 */
void MEDCouplingAMRAttribute::synchronizeFineToCoarse()
{
  if(_levs.empty())
    throw INTERP_KERNEL::Exception("MEDCouplingAMRAttribute::synchronizeFineToCoarse : not any levels in this !");
  std::size_t sz(_levs.size());
  //
  while(sz>1)
    {
      sz--;
      const MEDCouplingGridCollection *fine(_levs[sz]),*coarse(_levs[sz-1]);
      MEDCouplingGridCollection::SynchronizeFineToCoarse(_ghost_lev,fine,coarse);
    }
}

/*!
 * This method synchronizes from coarse to fine arrays and fine to fine each other (if _ghost_lev is >0). This method makes the hypothesis that \a this has been allocated before using
 * MEDCouplingAMRAttribute::alloc method.
 */
void MEDCouplingAMRAttribute::synchronizeCoarseToFine()
{
  if(_levs.empty())
    throw INTERP_KERNEL::Exception("MEDCouplingAMRAttribute::synchronizeCoarseToFine : not any levels in this !");
  std::size_t sz(_levs.size());
  //
  for(std::size_t i=1;i<sz;i++)
    {
      const MEDCouplingGridCollection *fine(_levs[i]),*coarse(_levs[i-1]);
      MEDCouplingGridCollection::SynchronizeCoarseToFine(_ghost_lev,coarse,fine);
    }
}

/*!
 * This method synchronizes the ghost zone of all patches (excepted the god father one).
 * This method operates in 4 steps. Larger is the number of steps more accurate is the information in the ghost zone.
 *
 * - firstly coarse to fine with no interactions between brother patches.
 * - secondly connected brother patches in a same master patch are updated.
 * - thirdly connected nephew patches are updated each other.
 * - forthly nth generation cousin patches are updated each other.
 *
 * This method makes the hypothesis that \a this has been allocated before using MEDCouplingAMRAttribute::alloc method.
 * So if \a _ghost_lev == 0 this method has no effect.
 */
void MEDCouplingAMRAttribute::synchronizeAllGhostZones()
{
  if(_levs.empty())
    throw INTERP_KERNEL::Exception("MEDCouplingAMRAttribute::synchronizeFineEachOther : not any levels in this !");
  // 1st - synchronize from coarse to the finest all the patches (excepted the god father one)
  std::size_t sz(_levs.size());
  for(std::size_t i=1;i<sz;i++)
    {
      const MEDCouplingGridCollection *fine(_levs[i]),*coarse(_levs[i-1]);
      MEDCouplingGridCollection::SynchronizeCoarseToFineOnlyInGhostZone(_ghost_lev,coarse,fine);
    }
  // 2nd - classical direct sublevel inside common patch
  for(std::size_t i=1;i<sz;i++)
    {
      const MEDCouplingGridCollection *fine(_levs[i]);
      if(!fine)
        throw INTERP_KERNEL::Exception("MEDCouplingAMRAttribute::synchronizeFineEachOtherInGhostZone : presence of a NULL element !");
      fine->synchronizeFineEachOther(_ghost_lev,_neighbors[i]);
    }
  // 3rd - mixed level
  for(std::vector< std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *> >::const_iterator it=_mixed_lev_neighbors.begin();it!=_mixed_lev_neighbors.end();it++)
    {
      const DataArrayDoubleCollection *firstDAC(&findCollectionAttachedTo((*it).first->getMesh())),*secondDAC(&findCollectionAttachedTo((*it).second->getMesh()));
      DataArrayDoubleCollection::SynchronizeGhostZoneOfOneUsingTwo(_ghost_lev,(*it).first,firstDAC,(*it).second,secondDAC);
    }
  // 4th - same level but with far ancestor.
  for(std::size_t i=1;i<sz;i++)
    {
      const MEDCouplingGridCollection *fine(_levs[i]);
      fine->synchronizeFineEachOtherExt(_ghost_lev,_cross_lev_neighbors[i]);
    }
}

/*!
 * This method allocates all DataArrayDouble instances stored recursively in \a this.
 *
 * \sa dealloc
 */
void MEDCouplingAMRAttribute::alloc()
{
  _tlc.resetState();
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingGridCollection> >::iterator it=_levs.begin();it!=_levs.end();it++)
    {
      MEDCouplingGridCollection *elt(*it);
      if(elt)
        elt->alloc(_ghost_lev);
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

bool MEDCouplingAMRAttribute::changeGodFather(MEDCouplingCartesianAMRMeshGen *gf)
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

MEDCouplingAMRAttribute::MEDCouplingAMRAttribute(MEDCouplingCartesianAMRMeshGen *gf, const std::vector< std::pair<std::string,int> >& fieldNames, int ghostLev):MEDCouplingDataForGodFather(gf),_ghost_lev(ghostLev)
{
  //gf non empty, checked by constructor
  int maxLev(gf->getMaxNumberOfLevelsRelativeToThis());
  _levs.resize(maxLev);
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
  // updates cross levels neighbors
  _neighbors.resize(_levs.size());
  _cross_lev_neighbors.resize(_levs.size());
  if(_levs.empty())
    throw INTERP_KERNEL::Exception("constructor of MEDCouplingAMRAttribute : not any levels in this !");
  std::size_t sz(_levs.size());
  for(std::size_t i=1;i<sz;i++)
    {
      const MEDCouplingGridCollection *fine(_levs[i]);
      if(!fine)
        throw INTERP_KERNEL::Exception("constructor of MEDCouplingAMRAttribute : presence of a NULL element !");
      _neighbors[i]=fine->findNeighbors(_ghost_lev);
      if(i!=sz-1)
        {
          for(std::vector< std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *> >::const_iterator it=_neighbors[i].begin();it!=_neighbors[i].end();it++)
            {
              MEDCouplingCartesianAMRPatch::FindNeighborsOfSubPatchesOf(_ghost_lev,(*it).first,(*it).second,_mixed_lev_neighbors);
              std::vector< std::vector < std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *> > > neighs2(MEDCouplingCartesianAMRPatch::FindNeighborsOfSubPatchesOfSameLev(_ghost_lev,(*it).first,(*it).second));
              std::size_t fullLev(i+neighs2.size());
              if(fullLev>=sz)
                throw INTERP_KERNEL::Exception("constructor of MEDCouplingAMRAttribute : internal error ! something is wrong in computation of cross level neighbors !");
              std::size_t ii(i+1);
              for(std::vector< std::vector< std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *> > >::const_iterator it0=neighs2.begin();it0!=neighs2.end();it0++,ii++)
                _cross_lev_neighbors[ii].insert(_cross_lev_neighbors[ii].end(),(*it0).begin(),(*it0).end());
            }
        }
    }
}

const DataArrayDoubleCollection& MEDCouplingAMRAttribute::findCollectionAttachedTo(const MEDCouplingCartesianAMRMeshGen *m) const
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingGridCollection> >::const_iterator it=_levs.begin();it!=_levs.end();it++)
    {
      const MEDCouplingGridCollection *elt(*it);
      if(elt)
        {
          int tmp(-1);
          if(elt->presenceOf(m,tmp))
            {
              return elt->getFieldsAt(tmp);
            }
        }
    }
  throw INTERP_KERNEL::Exception("MEDCouplingAMRAttribute::findCollectionAttachedTo : unable to find such part of mesh in this !");
}
