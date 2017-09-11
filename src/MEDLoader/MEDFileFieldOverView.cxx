// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
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

#include "MEDFileFieldOverView.hxx"
#include "MEDFileField.hxx"
#include "MEDFileMesh.hxx"

#include "MEDCouplingFieldDiscretization.hxx"
#include "CellModel.hxx"

using namespace MEDCoupling;

const unsigned char MEDMeshMultiLev::PARAMEDMEM_2_VTKTYPE[MEDMeshMultiLev::PARAMEDMEM_2_VTKTYPE_LGTH]=
{1,3,21,5,9,7,22,34,23,28,35,255,255,255,10,14,13,255,12,255,24,255,16,27,255,26,255,29,32,255,25,42,36,4};

const unsigned char MEDMeshMultiLev::HEXA27_PERM_ARRAY[27]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,24,22,21,23,20,25,26};

const char MEDFileField1TSStructItem2::NEWLY_CREATED_PFL_NAME[]="???";

MEDFileMeshStruct *MEDFileMeshStruct::New(const MEDFileMesh *mesh)
{
  return new MEDFileMeshStruct(mesh);
}

std::size_t MEDFileMeshStruct::getHeapMemorySizeWithoutChildren() const
{
  std::size_t ret(0);
  for(std::vector< std::vector<int> >::const_iterator it0=_geo_types_distrib.begin();it0!=_geo_types_distrib.end();it0++)
    ret+=(*it0).capacity()*sizeof(int);
  ret+=_geo_types_distrib.capacity()*sizeof(std::vector<int>);
  return ret;
}

std::vector<const BigMemoryObject *> MEDFileMeshStruct::getDirectChildrenWithNull() const
{
  return std::vector<const BigMemoryObject *>();
}

MEDFileMeshStruct::MEDFileMeshStruct(const MEDFileMesh *mesh):_mesh(mesh)
{
  std::vector<int> levs(mesh->getNonEmptyLevels());
  _name=mesh->getName();
  _nb_nodes=mesh->getNumberOfNodes();
  if(!levs.empty())
    {
      _geo_types_distrib.resize(-(*std::min_element(levs.begin(),levs.end()))+1);
      for(std::vector<int>::const_iterator lev=levs.begin();lev!=levs.end();lev++)
        _geo_types_distrib[-(*lev)]=mesh->getDistributionOfTypes(*lev);
    }
}

int MEDFileMeshStruct::getLevelOfGeoType(INTERP_KERNEL::NormalizedCellType t) const
{
  int j=0;
  for(std::vector< std::vector<int> >::const_iterator it1=_geo_types_distrib.begin();it1!=_geo_types_distrib.end();it1++,j--)
    {
      std::size_t sz=(*it1).size();
      if(sz%3!=0)
        throw INTERP_KERNEL::Exception("MEDFileMeshStruct::getLevelOfGeoType : internal error in code !");
      std::size_t nbGeo=sz/3;
      for(std::size_t i=0;i<nbGeo;i++)
        if((*it1)[3*i]==(int)t)
          return j;
    }
  throw INTERP_KERNEL::Exception("MEDFileMeshStruct::getLevelOfGeoType : The specified geometric type is not present in the mesh structure !");
}

/*!
 * \sa MEDFileMeshStruct::doesManageGeoType
 */
int MEDFileMeshStruct::getNumberOfElemsOfGeoType(INTERP_KERNEL::NormalizedCellType t) const
{
  for(std::vector< std::vector<int> >::const_iterator it1=_geo_types_distrib.begin();it1!=_geo_types_distrib.end();it1++)
    {
      std::size_t sz=(*it1).size();
      if(sz%3!=0)
        throw INTERP_KERNEL::Exception("MEDFileMeshStruct::getNumberOfElemsOfGeoType : internal error in code !");
      std::size_t nbGeo=sz/3;
      for(std::size_t i=0;i<nbGeo;i++)
        if((*it1)[3*i]==(int)t)
          return (*it1)[3*i+1];
    }
  throw INTERP_KERNEL::Exception("The specified geometric type is not present in the mesh structure !");
}

/*!
 * \sa MEDFileMeshStruct::getNumberOfElemsOfGeoType
 */
bool MEDFileMeshStruct::doesManageGeoType(INTERP_KERNEL::NormalizedCellType t) const
{
  for(std::vector< std::vector<int> >::const_iterator it1=_geo_types_distrib.begin();it1!=_geo_types_distrib.end();it1++)
    {
      std::size_t sz=(*it1).size();
      if(sz%3!=0)
        throw INTERP_KERNEL::Exception("MEDFileMeshStruct::doesManageGeoType : internal error in code !");
      std::size_t nbGeo=sz/3;
      for(std::size_t i=0;i<nbGeo;i++)
        if((*it1)[3*i]==(int)t)
          return true;
    }
  return false;
}

void MEDFileMeshStruct::appendIfImplicitType(INTERP_KERNEL::NormalizedCellType t)
{
  if(!_mesh->hasImplicitPart())
    throw INTERP_KERNEL::Exception("MEDFileMeshStruct::appendIfImplicitType : by default no implicit geo type can be appended !");
  static const char MSG[]="MEDFileMeshStruct::appendIfImplicitType : the distribution does not looks like structured standard !";
  if(_geo_types_distrib.size()!=1)
    throw INTERP_KERNEL::Exception(MSG);
  std::size_t sz(_geo_types_distrib[0].size());
  if(sz%3!=0)
    throw INTERP_KERNEL::Exception("MEDFileMeshStruct::appendIfImplicitType : internal error in code !");
  std::size_t nbGeo(sz/3);
  if(nbGeo!=1)
    throw INTERP_KERNEL::Exception(MSG);
  std::vector<int> arr(3); arr[0]=(int)t; arr[1]=_mesh->buildImplicitPartIfAny(t); arr[2]=-1;
  _geo_types_distrib.push_back(arr);
}


int MEDFileMeshStruct::getNumberOfLevs() const
{
  return (int)_geo_types_distrib.size();
}

int MEDFileMeshStruct::getNumberOfGeoTypesInLev(int relativeLev) const
{
  int pos(-relativeLev);
  if(pos<0 || pos>=(int)_geo_types_distrib.size())
    throw INTERP_KERNEL::Exception("MEDFileMeshStruct::getNumberOfGeoTypesInLev : invalid level specified !");
  std::size_t sz=_geo_types_distrib[pos].size();
  if(sz%3!=0)
    throw INTERP_KERNEL::Exception("MEDFileMeshStruct::getNumberOfGeoTypesInLev : internal error in code !");
  return (int)(sz/3);
}

//=

std::size_t MEDMeshMultiLev::getHeapMemorySizeWithoutChildren() const
{
  return 0;
}

std::vector<const BigMemoryObject *> MEDMeshMultiLev::getDirectChildrenWithNull() const
{
  return std::vector<const BigMemoryObject *>();
}

MEDMeshMultiLev *MEDMeshMultiLev::New(const MEDFileMesh *m, const std::vector<int>& levs)
{
  if(!m)
    throw INTERP_KERNEL::Exception("MEDMeshMultiLev::New : null input pointer !");
  const MEDFileUMesh *um(dynamic_cast<const MEDFileUMesh *>(m));
  if(um)
    return MEDUMeshMultiLev::New(um,levs);
  const MEDFileCMesh *cm(dynamic_cast<const MEDFileCMesh *>(m));
  if(cm)
    return MEDCMeshMultiLev::New(cm,levs);
  const MEDFileCurveLinearMesh *clm(dynamic_cast<const MEDFileCurveLinearMesh *>(m));
  if(clm)
    return MEDCurveLinearMeshMultiLev::New(clm,levs);
  throw INTERP_KERNEL::Exception("MEDMeshMultiLev::New : unrecognized type of mesh ! Must be in [MEDFileUMesh,MEDFileCMesh,MEDFileCurveLinearMesh] !");
}

MEDMeshMultiLev *MEDMeshMultiLev::New(const MEDFileMesh *m, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayInt *>& pfls, const std::vector<int>& nbEntities)
{
  if(!m)
    throw INTERP_KERNEL::Exception("MEDMeshMultiLev::New 2 : null input pointer !");
  const MEDFileUMesh *um(dynamic_cast<const MEDFileUMesh *>(m));
  if(um)
    return MEDUMeshMultiLev::New(um,gts,pfls,nbEntities);
  const MEDFileCMesh *cm(dynamic_cast<const MEDFileCMesh *>(m));
  if(cm)
    return MEDCMeshMultiLev::New(cm,gts,pfls,nbEntities);
  const MEDFileCurveLinearMesh *clm(dynamic_cast<const MEDFileCurveLinearMesh *>(m));
  if(clm)
    return MEDCurveLinearMeshMultiLev::New(clm,gts,pfls,nbEntities);
  throw INTERP_KERNEL::Exception("MEDMeshMultiLev::New 2 : unrecognized type of mesh ! Must be in [MEDFileUMesh,MEDFileCMesh,MEDFileCurveLinearMesh] !");
}

MEDMeshMultiLev *MEDMeshMultiLev::NewOnlyOnNode(const MEDFileMesh *m, const DataArrayInt *pflOnNode)
{
  MCAuto<MEDMeshMultiLev> ret(MEDMeshMultiLev::New(m,m->getNonEmptyLevels()));
  ret->selectPartOfNodes(pflOnNode);
  return ret.retn();
}

void MEDMeshMultiLev::setNodeReduction(const DataArrayInt *nr)
{
  if(nr)
    nr->incrRef();
  _node_reduction=const_cast<DataArrayInt*>(nr);
}

void MEDMeshMultiLev::setCellReduction(const DataArrayInt *cr)
{
  if(_pfls.size()!=1)
    throw INTERP_KERNEL::Exception("MEDMeshMultiLev::setCellReduction : can be used only for single geo type mesh !");
  _pfls[0]=const_cast<DataArrayInt*>(cr);
  if(cr)
    cr->incrRef();
}

bool MEDMeshMultiLev::isFastlyTheSameStruct(const MEDFileField1TSStructItem& fst, const MEDFileFieldGlobsReal *globs) const
{
  if(fst.getType()==ON_NODES)
    {
      if(fst.getNumberOfItems()!=1)
        throw INTERP_KERNEL::Exception("MEDMeshMultiLev::isFastlyTheSameStruct : unexpected situation for nodes !");
      const MEDFileField1TSStructItem2& p(fst[0]);
      std::string pflName(p.getPflName());
      const DataArrayInt *nr(_node_reduction);
      if(pflName.empty() && !nr)
        return true;
      if(!pflName.empty() && !nr)
        return false;
      if(pflName==nr->getName())
        return true;
      return false;
    }
  else
    {
      std::size_t sz(fst.getNumberOfItems());
      if(sz!=_geo_types.size())
        return false;
      int strt(0);
      for(std::size_t i=0;i<sz;i++)
        {
          const MEDFileField1TSStructItem2& p(fst[i]);
          if(!p.isFastlyEqual(strt,_geo_types[i],getPflNameOfId(i).c_str()))
            return false;
        }
      return true;
    }
}

DataArray *MEDMeshMultiLev::buildDataArray(const MEDFileField1TSStructItem& fst, const MEDFileFieldGlobsReal *globs, const DataArray *vals) const
{
  MCAuto<DataArray> ret(const_cast<DataArray *>(vals)); ret->incrRef();
  if(isFastlyTheSameStruct(fst,globs))
    return ret.retn();
  else
    return constructDataArray(fst,globs,vals);
}

/*!
 * \param [out] famIds - Can be null. If not null the instance has to be dealt by the caller (decrRef).
 * \param [out] isWithoutCopy - When true the returned instance \a famIds if not null is directly those in the data structure.
 */
void MEDMeshMultiLev::retrieveFamilyIdsOnCells(DataArrayInt *& famIds, bool& isWithoutCopy) const
{
  const DataArrayInt *fids(_cell_fam_ids);
  if(!fids)
    { famIds=0; isWithoutCopy=true; return ; }
  std::size_t sz(_geo_types.size());
  bool presenceOfPfls(false);
  for(std::size_t i=0;i<sz && !presenceOfPfls;i++)
    {
      const DataArrayInt *pfl(_pfls[i]);
      if(pfl)
        presenceOfPfls=true;
    }
  if(!presenceOfPfls)
    { famIds=const_cast<DataArrayInt *>(fids); famIds->incrRef(); isWithoutCopy=_mesh->isObjectInTheProgeny(famIds); return ; }
  //bad luck the slowest part
  isWithoutCopy=false;
  std::vector< MCAuto<DataArrayInt> > retSafe(sz);
  std::vector< const DataArrayInt *> ret(sz);
  int start(0);
  for(std::size_t i=0;i<sz;i++)
    {
      const DataArrayInt *pfl(_pfls[i]);
      int lgth(_nb_entities[i]);
      if(pfl)
        {
          MCAuto<DataArrayInt> tmp(fids->selectByTupleIdSafeSlice(start,start+lgth,1));
          retSafe[i]=tmp->selectByTupleIdSafe(pfl->begin(),pfl->end());
        }
      else
        {
          retSafe[i]=fids->selectByTupleIdSafeSlice(start,start+lgth,1);
        }
      ret[i]=retSafe[i];
      start+=lgth;
    }
  famIds=DataArrayInt::Aggregate(ret);
}

/*!
 * \param [out] numIds - Can be null. If not null the instance has to be dealt by the caller (decrRef).
 * \param [out] isWithoutCopy - When true the returned instance \a numIds if not null is directly those in the data structure.
 */
void MEDMeshMultiLev::retrieveNumberIdsOnCells(DataArrayInt *& numIds, bool& isWithoutCopy) const
{
  const DataArrayInt *nids(_cell_num_ids);
  if(!nids)
    { numIds=0; isWithoutCopy=true; return ; }
  std::size_t sz(_geo_types.size());
  bool presenceOfPfls(false);
  for(std::size_t i=0;i<sz && !presenceOfPfls;i++)
    {
      const DataArrayInt *pfl(_pfls[i]);
      if(pfl)
        presenceOfPfls=true;
    }
  if(!presenceOfPfls)
    { numIds=const_cast<DataArrayInt *>(nids); numIds->incrRef(); isWithoutCopy=_mesh->isObjectInTheProgeny(numIds); return ; }
  //bad luck the slowest part
  isWithoutCopy=false;
  std::vector< MCAuto<DataArrayInt> > retSafe(sz);
  std::vector< const DataArrayInt *> ret(sz);
  int start(0);
  for(std::size_t i=0;i<sz;i++)
    {
      const DataArrayInt *pfl(_pfls[i]);
      int lgth(_nb_entities[i]);
      if(pfl)
        {
          MCAuto<DataArrayInt> tmp(nids->selectByTupleIdSafeSlice(start,start+lgth,1));
          retSafe[i]=tmp->selectByTupleIdSafe(pfl->begin(),pfl->end());
        }
      else
        {
          retSafe[i]=nids->selectByTupleIdSafeSlice(start,start+lgth,1);
        }
      ret[i]=retSafe[i];
      start+=lgth;
    }
  numIds=DataArrayInt::Aggregate(ret);
}

/*!
 * \param [out] famIds - Can be null. If not null the instance has to be dealt by the caller (decrRef).
 * \param [out] isWithoutCopy - When true the returned instance \a famIds if not null is directly those in the data structure.
 */
void MEDMeshMultiLev::retrieveFamilyIdsOnNodes(DataArrayInt *& famIds, bool& isWithoutCopy) const
{
  const DataArrayInt *fids(_node_fam_ids);
  if(!fids)
    { famIds=0; isWithoutCopy=true; return ; }
  const DataArrayInt *nr(_node_reduction);
  if(nr)
    {
      isWithoutCopy=false;
      famIds=fids->selectByTupleIdSafe(nr->begin(),nr->end());
    }
  else
    {
      famIds=const_cast<DataArrayInt *>(fids); famIds->incrRef();
      isWithoutCopy=_mesh->isObjectInTheProgeny(famIds);
    }
}

/*!
 * \param [out] numIds - Can be null. If not null the instance has to be dealt by the caller (decrRef).
 * \param [out] isWithoutCopy - When true the returned instance \a numIds if not null is directly those in the data structure.
 */
void MEDMeshMultiLev::retrieveNumberIdsOnNodes(DataArrayInt *& numIds, bool& isWithoutCopy) const
{
  const DataArrayInt *fids(_node_num_ids);
  if(!fids)
    { numIds=0; isWithoutCopy=true; return ; }
  const DataArrayInt *nr(_node_reduction);
  if(nr)
    {
      isWithoutCopy=false;
      numIds=fids->selectByTupleIdSafe(nr->begin(),nr->end());
    }
  else
    {
      numIds=const_cast<DataArrayInt *>(fids); numIds->incrRef();
      isWithoutCopy=_mesh->isObjectInTheProgeny(numIds);
    }
}

/*!
 * This method returns, if any, a new object containing the global node ids **BUT CONTRARY TO OTHER RETRIEVE METHODS** the returned object is always a NON AGGREGATED object. So the returned object if not null
 * can be used as this safely.
 */
DataArrayInt *MEDMeshMultiLev::retrieveGlobalNodeIdsIfAny() const
{
  const MEDFileUMesh *umesh(dynamic_cast<const MEDFileUMesh *>(_mesh));
  if(!umesh)
    return 0;
  const PartDefinition *pd(umesh->getPartDefAtLevel(1));
  if(!pd)
    return 0;
  MCAuto<DataArrayInt> tmp(pd->toDAI());
  const DataArrayInt *tmpCpp(tmp);
  if(!tmpCpp)
    return 0;
  //
  const DataArrayInt *nr(_node_reduction);
  if(nr)
    return tmp->selectByTupleIdSafe(nr->begin(),nr->end());
  else
    return tmp->deepCopy();// Yes a deep copy is needed because this method has to return a non aggregated object !
}

std::vector< INTERP_KERNEL::NormalizedCellType > MEDMeshMultiLev::getGeoTypes() const
{
  return _geo_types;
}

void MEDMeshMultiLev::setFamilyIdsOnCells(DataArrayInt *famIds)
{
  _cell_fam_ids=famIds;
  if(famIds)
    famIds->incrRef();
}

void MEDMeshMultiLev::setNumberIdsOnCells(DataArrayInt *numIds)
{
  _cell_num_ids=numIds;
  if(numIds)
    numIds->incrRef();
}

void MEDMeshMultiLev::setFamilyIdsOnNodes(DataArrayInt *famIds)
{
  _node_fam_ids=famIds;
  if(famIds)
    famIds->incrRef();
}

void MEDMeshMultiLev::setNumberIdsOnNodes(DataArrayInt *numIds)
{
  _node_num_ids=numIds;
  if(numIds)
    numIds->incrRef();
}

std::string MEDMeshMultiLev::getPflNameOfId(int id) const
{
  std::size_t sz(_pfls.size());
  if(id<0 || id>=(int)sz)
    throw INTERP_KERNEL::Exception("MEDMeshMultiLev::getPflNameOfId : invalid input id !");
  const DataArrayInt *pfl(_pfls[id]);
  if(!pfl)
    return std::string("");
  return pfl->getName();
}

/*!
 * Returns the number of cells having geometric type \a t.
 * The profiles are **NOT** taken into account here.
 */
int MEDMeshMultiLev::getNumberOfCells(INTERP_KERNEL::NormalizedCellType t) const
{
  std::size_t sz(_nb_entities.size());
  for(std::size_t i=0;i<sz;i++)
    if(_geo_types[i]==t)
      return _nb_entities[i];
  throw INTERP_KERNEL::Exception("MEDMeshMultiLev::getNumberOfCells : not existing geometric type in this !");
}

int MEDMeshMultiLev::getNumberOfNodes() const
{
  return _nb_nodes;
}

DataArray *MEDMeshMultiLev::constructDataArray(const MEDFileField1TSStructItem& fst, const MEDFileFieldGlobsReal *globs, const DataArray *vals) const
{
  if(fst.getType()==ON_NODES)
    {
      if(fst.getNumberOfItems()!=1)
        throw INTERP_KERNEL::Exception("MEDMeshMultiLev::constructDataArray : unexpected situation for nodes !");
      const MEDFileField1TSStructItem2& p(fst[0]);
      std::string pflName(p.getPflName());
      const DataArrayInt *nr(_node_reduction);
      if(pflName.empty() && !nr)
        return vals->deepCopy();
      if(pflName.empty() && nr)
        throw INTERP_KERNEL::Exception("MEDMeshMultiLev::constructDataArray : unexpected situation for nodes 2 !");
      if(!pflName.empty() && nr)
        {
          MCAuto<DataArrayInt> p1(globs->getProfile(pflName.c_str())->deepCopy());
          MCAuto<DataArrayInt> p2(nr->deepCopy());
          p1->sort(true); p2->sort(true);
          if(!p1->isEqualWithoutConsideringStr(*p2))
            throw INTERP_KERNEL::Exception("MEDMeshMultiLev::constructDataArray : it appears that a profile on nodes does not cover the cells correctly !");
          p1=DataArrayInt::FindPermutationFromFirstToSecond(globs->getProfile(pflName.c_str()),nr);
          MCAuto<DataArray> ret(vals->deepCopy());
          ret->renumberInPlace(p1->begin());
          return ret.retn();
        }
      if(!pflName.empty() && !nr)
        {
          MCAuto<DataArrayInt> p1(globs->getProfile(pflName.c_str())->deepCopy());
          p1->sort(true);
          if(!p1->isIota(getNumberOfNodes()))
            throw INTERP_KERNEL::Exception("MEDMeshMultiLev::constructDataArray : unexpected situation for nodes 4 !");
          MCAuto<DataArray> ret(vals->deepCopy());
          ret->renumberInPlace(globs->getProfile(pflName.c_str())->begin());
          return ret.retn();
        }
      throw INTERP_KERNEL::Exception("MEDMeshMultiLev::constructDataArray : unexpected situation for nodes 5 !");
    }
  else
    {
      std::size_t sz(fst.getNumberOfItems());
      std::set<INTERP_KERNEL::NormalizedCellType> s(_geo_types.begin(),_geo_types.end());
      if(s.size()!=_geo_types.size())
        throw INTERP_KERNEL::Exception("MEDMeshMultiLev::constructDataArray : unexpected situation for cells 2 !");
      std::vector< const DataArray *> arr(s.size());
      std::vector< MCAuto<DataArray> > arrSafe(s.size());
      int iii(0);
      int nc(vals->getNumberOfComponents());
      std::vector<std::string> compInfo(vals->getInfoOnComponents());
      for(std::vector< INTERP_KERNEL::NormalizedCellType >::const_iterator it=_geo_types.begin();it!=_geo_types.end();it++,iii++)
        {
          const DataArrayInt *thisP(_pfls[iii]);
          std::vector<const MEDFileField1TSStructItem2 *> ps;
          for(std::size_t i=0;i<sz;i++)
            {
              const MEDFileField1TSStructItem2& p(fst[i]);
              if(p.getGeo()==*it)
                ps.push_back(&p);
            }
          if(ps.empty())
            throw INTERP_KERNEL::Exception("MEDMeshMultiLev::constructDataArray : unexpected situation for cells 1 !");
          if(ps.size()==1)
            {
              int nbi(ps[0]->getNbOfIntegrationPts(globs));
              const DataArrayInt *otherP(ps[0]->getPfl(globs));
              const std::pair<int,int>& strtStop(ps[0]->getStartStop());
              MCAuto<DataArray> ret(vals->selectByTupleIdSafeSlice(strtStop.first,strtStop.second,1));
              if(!thisP && !otherP)
                {
                  arrSafe[iii]=ret; arr[iii]=ret;
                  continue;
                }
              if(thisP && otherP)
                {
                  MCAuto<DataArrayInt> p1(otherP->invertArrayN2O2O2N(getNumberOfCells(ps[0]->getGeo())));
                  MCAuto<DataArrayInt> p2(thisP->deepCopy());
                  p2->transformWithIndArr(p1->begin(),p1->end());
                  //p1=p2->findIdsNotEqual(-1);
                  //p1=p2->selectByTupleIdSafe(p1->begin(),p1->end());
                  ret->rearrange(nbi*nc); ret=ret->selectByTupleIdSafe(p2->begin(),p2->end()); ret->rearrange(nc); ret->setInfoOnComponents(compInfo);
                  arrSafe[iii]=ret; arr[iii]=ret;
                  continue;
                }
              if(!thisP && otherP)
                {
                  MCAuto<DataArrayInt> p1(otherP->deepCopy());
                  p1->sort(true);
                  p1->checkAllIdsInRange(0,getNumberOfCells(ps[0]->getGeo()));
                  p1=DataArrayInt::FindPermutationFromFirstToSecond(otherP,p1);
                  ret->rearrange(nbi*nc); ret->renumberInPlace(p1->begin()); ret->rearrange(nc); ret->setInfoOnComponents(compInfo);
                  arrSafe[iii]=ret; arr[iii]=ret;
                  continue;
                }
              throw INTERP_KERNEL::Exception("MEDMeshMultiLev::constructDataArray : unexpected situation for cells 3 !");
            }
          else
            {
              std::vector< const DataArrayInt * >otherPS(ps.size());
              std::vector< const DataArray * > arr2(ps.size());
              std::vector< MCAuto<DataArray> > arr2Safe(ps.size());
              std::vector< const DataArrayInt * > nbis(ps.size());
              std::vector< MCAuto<DataArrayInt> > nbisSafe(ps.size());
              int jj(0);
              for(std::vector<const MEDFileField1TSStructItem2 *>::const_iterator it2=ps.begin();it2!=ps.end();it2++,jj++)
                {
                  int nbi((*it2)->getNbOfIntegrationPts(globs));
                  const DataArrayInt *otherPfl((*it2)->getPfl(globs));
                  const std::pair<int,int>& strtStop((*it2)->getStartStop());
                  MCAuto<DataArray> ret2(vals->selectByTupleIdSafeSlice(strtStop.first,strtStop.second,1));
                  if(!otherPfl)
                    throw INTERP_KERNEL::Exception("MEDMeshMultiLev::constructDataArray : unexpected situation for cells 4 !");
                  arr2[jj]=ret2; arr2Safe[jj]=ret2; otherPS[jj]=otherPfl;
                  nbisSafe[jj]=DataArrayInt::New(); nbisSafe[jj]->alloc(otherPfl->getNumberOfTuples(),1); nbisSafe[jj]->fillWithValue(nbi);
                  nbis[jj]=nbisSafe[jj];
                }
              MCAuto<DataArray> arr3(DataArray::Aggregate(arr2));
              MCAuto<DataArrayInt> otherP(DataArrayInt::Aggregate(otherPS));
              MCAuto<DataArrayInt> zenbis(DataArrayInt::Aggregate(nbis));
              MCAuto<DataArrayInt> otherPN(otherP->invertArrayN2O2O2N(getNumberOfCells(*it)));
              MCAuto<DataArrayInt> p1;
              if(thisP)
                p1=DataArrayInt::FindPermutationFromFirstToSecond(otherP,thisP);
              else
                p1=otherP->deepCopy();
              MCAuto<DataArrayInt> zenbisN(zenbis->renumber(p1->begin()));
              zenbisN->computeOffsetsFull();
              jj=0;
              for(std::vector<const MEDFileField1TSStructItem2 *>::const_iterator it2=ps.begin();it2!=ps.end();it2++,jj++)
                {
                  //int nbi((*it2)->getNbOfIntegrationPts(globs));
                  const DataArrayInt *otherPfl((*it2)->getPfl(globs));
                  const std::pair<int,int>& strtStop((*it2)->getStartStop());
                  MCAuto<DataArray> ret2(vals->selectByTupleIdSafeSlice(strtStop.first,strtStop.second,1));
                  //
                  MCAuto<DataArrayInt> p2(otherPfl->deepCopy());
                  p2->transformWithIndArr(otherPN->begin(),otherPN->end());
                  p2->transformWithIndArr(p1->begin(),p1->end());
                  MCAuto<DataArrayInt> idsN(p2->buildExplicitArrByRanges(zenbisN));
                  arr3->setPartOfValuesBase3(ret2,idsN->begin(),idsN->end(),0,nc,1);
                }
              arrSafe[iii]=arr3; arr[iii]=arr3;
              continue;
            }
        }
      return DataArray::Aggregate(arr);
    }
}

/*!
 * This method is called to add NORM_POINT1 cells in \a this so that orphan nodes in \a verticesToAdd will be fetched.
 */
void MEDMeshMultiLev::appendVertices(const DataArrayInt *verticesToAdd, DataArrayInt *nr)
{
  int nbOfVertices(verticesToAdd->getNumberOfTuples());
  std::size_t sz(_pfls.size());
  _pfls.resize(sz+1);
  _geo_types.resize(sz+1,INTERP_KERNEL::NORM_POINT1);
  _nb_entities.resize(sz+1,nbOfVertices);
  _node_reduction=nr; nr->incrRef();
  _nb_nodes+=nbOfVertices;
  const DataArrayInt *cf(_cell_fam_ids),*cn(_cell_num_ids),*nf(_node_fam_ids),*nn(_node_num_ids);
  if(cf)
    {
      MCAuto<DataArrayInt> tmp;
      std::vector<const DataArrayInt *> a(2);
      a[0]=cf;
      if(nf)
        tmp=nf->selectByTupleIdSafe(verticesToAdd->begin(),verticesToAdd->end());
      else
        {
          tmp=DataArrayInt::New(); tmp->alloc(nbOfVertices,1); tmp->fillWithZero();
        }
      a[1]=tmp;
      _cell_fam_ids=DataArrayInt::Aggregate(a);
    }
  if(cn)
    {
      MCAuto<DataArrayInt> tmp;
      std::vector<const DataArrayInt *> a(2);
      a[0]=cn;
      if(nn)
        tmp=nn->selectByTupleIdSafe(verticesToAdd->begin(),verticesToAdd->end());
      else
        {
          tmp=DataArrayInt::New(); tmp->alloc(nbOfVertices,1); tmp->fillWithZero();
        }
      a[1]=tmp;
      _cell_num_ids=DataArrayInt::Aggregate(a);
    }
}

MEDMeshMultiLev::MEDMeshMultiLev(const MEDFileMesh *mesh):_mesh(mesh),_nb_nodes(0)
{
}

MEDMeshMultiLev::MEDMeshMultiLev(const MEDFileMesh *mesh, int nbNodes, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayInt *>& pfls, const std::vector<int>& nbEntities):_mesh(mesh),_geo_types(gts),_nb_entities(nbEntities),_nb_nodes(nbNodes)
{
  std::size_t sz(_geo_types.size());
  if(sz!=pfls.size() || sz!=nbEntities.size())
    throw INTERP_KERNEL::Exception("MEDMeshMultiLev::MEDMeshMultiLev : input vector must have the same size !");
  _pfls.resize(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      if(pfls[i])
        pfls[i]->incrRef();
      _pfls[i]=const_cast<DataArrayInt *>(pfls[i]);
    }
}

MEDMeshMultiLev::MEDMeshMultiLev(const MEDMeshMultiLev& other):RefCountObject(other),_mesh(other._mesh),_pfls(other._pfls),_geo_types(other._geo_types),_nb_entities(other._nb_entities),_node_reduction(other._node_reduction),_nb_nodes(other._nb_nodes),_cell_fam_ids(other._cell_fam_ids),_cell_num_ids(other._cell_num_ids),_node_fam_ids(other._node_fam_ids),_node_num_ids(other._node_num_ids)
{
}

//=

MEDUMeshMultiLev *MEDUMeshMultiLev::New(const MEDFileUMesh *m, const std::vector<int>& levs)
{
  return new MEDUMeshMultiLev(m,levs);
}

MEDUMeshMultiLev::MEDUMeshMultiLev(const MEDFileUMesh *m, const std::vector<int>& levs):MEDMeshMultiLev(m)
{
  if(!m)
    throw INTERP_KERNEL::Exception("MEDUMeshMultiLev constructor : null input pointer !");
  std::vector<MEDCoupling1GTUMesh *> v;
  for(std::vector<int>::const_iterator it=levs.begin();it!=levs.end();it++)
    {
      std::vector<MEDCoupling1GTUMesh *> vTmp(m->getDirectUndergroundSingleGeoTypeMeshes(*it));
      v.insert(v.end(),vTmp.begin(),vTmp.end());
    }
  std::size_t sz(v.size());
  if(v.empty())
    {
      _coords=m->getCoords(); _coords->incrRef();
    }
  _parts.resize(sz);
  _pfls.resize(sz);
  _geo_types.resize(sz);
  _nb_entities.resize(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      MEDCoupling1GTUMesh *obj(v[i]);
      if(obj)
        obj->incrRef();
      else
        throw INTERP_KERNEL::Exception("MEDUMeshMultiLev constructor : presence of a null pointer !");
      _parts[i]=obj;
      _geo_types[i]=obj->getCellModelEnum();
      _nb_entities[i]=obj->getNumberOfCells();
    }
  // ids fields management
  bool cellFamIdsNoCpy(levs.size()==1);
  if(cellFamIdsNoCpy)
    {
      const DataArrayInt *tmp(m->getFamilyFieldAtLevel(levs[0]));
      if(tmp)
        {
          tmp->incrRef();
          _cell_fam_ids=(const_cast<DataArrayInt *>(tmp));
        }
    }
  else
    {
      std::vector<const DataArrayInt *> tmps(levs.size());
      bool f(true);
      for(std::size_t i=0;i<levs.size();i++)
        {
          tmps[i]=m->getFamilyFieldAtLevel(levs[i]);
          if(!tmps[i])
            f=false;
        }
      if(f && !tmps.empty())
        _cell_fam_ids=DataArrayInt::Aggregate(tmps);
    }
  bool cellNumIdsNoCpy(levs.size()==1);
  if(cellNumIdsNoCpy)
    {
      const DataArrayInt *tmp(m->getNumberFieldAtLevel(levs[0]));
      if(tmp)
        {
          tmp->incrRef();
          _cell_num_ids=(const_cast<DataArrayInt *>(tmp));
        }
    }
  else
    {
      std::vector<const DataArrayInt *> tmps(levs.size());
      bool n(true);
      for(std::size_t i=0;i<levs.size();i++)
        {
          tmps[i]=m->getNumberFieldAtLevel(levs[i]);
          if(!tmps[i])
            n=false;
        }
      if(n && !tmps.empty())
        _cell_num_ids=DataArrayInt::Aggregate(tmps);
    }
  // node part
  {
    const DataArrayInt *tmp(m->getFamilyFieldAtLevel(1));
    if(tmp)
      {
        tmp->incrRef();
        _node_fam_ids=(const_cast<DataArrayInt *>(tmp));
      }
  }
  {
    const DataArrayInt *tmp(m->getNumberFieldAtLevel(1));
    if(tmp)
      {
        tmp->incrRef();
        _node_num_ids=(const_cast<DataArrayInt *>(tmp));
      }
  }
}

MEDUMeshMultiLev *MEDUMeshMultiLev::New(const MEDFileUMesh *m, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayInt *>& pfls, const std::vector<int>& nbEntities)
{
  return new MEDUMeshMultiLev(m,gts,pfls,nbEntities);
}

MEDUMeshMultiLev::MEDUMeshMultiLev(const MEDFileUMesh *m, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayInt *>& pfls, const std::vector<int>& nbEntities):MEDMeshMultiLev(m,m->getNumberOfNodes(),gts,pfls,nbEntities)
{
  std::size_t sz(gts.size());
  if(sz<1)
    throw INTERP_KERNEL::Exception("constructor of MEDUMeshMultiLev : number of different geo type must be >= 1 !");
  unsigned dim(INTERP_KERNEL::CellModel::GetCellModel(gts[0]).getDimension());
  _parts.resize(sz);
  bool isSameDim(true),isNoPfl(true);
  for(std::size_t i=0;i<sz;i++)
    {
      MEDCoupling1GTUMesh *elt(m->getDirectUndergroundSingleGeoTypeMesh(gts[i]));
      if(INTERP_KERNEL::CellModel::GetCellModel(gts[i]).getDimension()!=dim)
        isSameDim=false;
      if(pfls[i])
        isNoPfl=false;
      if(elt)
        elt->incrRef();
      _parts[i]=elt;
    }
  // ids fields management
  int lev((int)dim-m->getMeshDimension());
  if(isSameDim && isNoPfl && m->getGeoTypesAtLevel(lev)==gts)//optimized part
    {
      const DataArrayInt *famIds(m->getFamilyFieldAtLevel(lev));
      if(famIds)
        { _cell_fam_ids=const_cast<DataArrayInt*>(famIds); famIds->incrRef(); }
      const DataArrayInt *numIds(m->getNumberFieldAtLevel(lev));
      if(numIds)
        { _cell_num_ids=const_cast<DataArrayInt*>(numIds); numIds->incrRef(); }
      famIds=m->getFamilyFieldAtLevel(1);
      if(famIds)
        { _node_fam_ids=const_cast<DataArrayInt*>(famIds); famIds->incrRef(); }
      numIds=m->getNumberFieldAtLevel(1);
      if(numIds)
        { _node_num_ids=const_cast<DataArrayInt*>(numIds); numIds->incrRef(); }
      return ;
    }
  //
  std::vector< MCAuto<DataArrayInt> > famIdsSafe(sz);
  std::vector<const DataArrayInt *> famIds(sz);
  bool f(true);
  for(std::size_t i=0;i<sz;i++)
    {
      famIdsSafe[i]=m->extractFamilyFieldOnGeoType(gts[i]);
      famIds[i]=famIdsSafe[i];
      if(!famIds[i])
        f=false;
    }
  if(f)
    _cell_fam_ids=DataArrayInt::Aggregate(famIds);
  std::vector< MCAuto<DataArrayInt> > numIdsSafe(sz);
  std::vector<const DataArrayInt *> numIds(sz);
  bool n(true);
  for(std::size_t i=0;i<sz;i++)
    {
      numIdsSafe[i]=m->extractNumberFieldOnGeoType(gts[i]);
      numIds[i]=numIdsSafe[i];
      if(!numIds[i])
        n=false;
    }
  if(n)
    _cell_num_ids=DataArrayInt::Aggregate(numIds);
  // node ids management
  const DataArrayInt *nodeFamIds(m->getFamilyFieldAtLevel(1));
  if(nodeFamIds)
    { _node_fam_ids=const_cast<DataArrayInt*>(nodeFamIds); nodeFamIds->incrRef(); }
  const DataArrayInt *nodeNumIds(m->getNumberFieldAtLevel(1));
  if(nodeNumIds)
    { _node_num_ids=const_cast<DataArrayInt*>(nodeNumIds); nodeNumIds->incrRef(); }
}

void MEDUMeshMultiLev::selectPartOfNodes(const DataArrayInt *pflNodes)
{
  if(!pflNodes || !pflNodes->isAllocated())
    return ;
  std::size_t sz(_parts.size());
  std::vector< MCAuto<DataArrayInt> > a(sz);
  std::vector< const DataArrayInt *> aa(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      const DataArrayInt *pfl(_pfls[i]);
      MCAuto<MEDCoupling1GTUMesh> m(_parts[i]);
      if(pfl)
        m=dynamic_cast<MEDCoupling1GTUMesh *>(_parts[i]->buildPartOfMySelfKeepCoords(pfl->begin(),pfl->end()));
      DataArrayInt *cellIds=0;
      m->fillCellIdsToKeepFromNodeIds(pflNodes->begin(),pflNodes->end(),true,cellIds);
      MCAuto<DataArrayInt> cellIdsSafe(cellIds);
      MCAuto<MEDCouplingPointSet> m2(m->buildPartOfMySelfKeepCoords(cellIds->begin(),cellIds->end()));
      int tmp=-1;
      MCAuto<DataArrayInt> o2n(m2->getNodeIdsInUse(tmp));
      a[i]=o2n->invertArrayO2N2N2O(tmp); aa[i]=a[i];
      if(pfl)
        _pfls[i]=pfl->selectByTupleIdSafe(cellIds->begin(),cellIds->end());
      else
        _pfls[i]=cellIdsSafe;
    }
  if(!aa.empty())
    _node_reduction=DataArrayInt::Aggregate(aa);//general case
  else
    _node_reduction=pflNodes->deepCopy();//case where no cells in read mesh.
  _node_reduction->sort(true);
  _node_reduction=_node_reduction->buildUnique();
  if(_node_reduction->getNumberOfTuples()==pflNodes->getNumberOfTuples())
    return ;//This is the classical case where the input node profile corresponds perfectly to a subset of cells in _parts
  if(_node_reduction->getNumberOfTuples()>pflNodes->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("MEDUMeshMultiLev::selectPartOfNodes : internal error in MEDCoupling during cell select from a list of nodes !");
  // Here the cells available in _parts is not enough to cover all the nodes in pflNodes. So adding vertices cells in _parts...
  MCAuto<DataArrayInt> pflNodes2(pflNodes->deepCopy());
  pflNodes2->sort(true);
  MCAuto<DataArrayInt> diff(pflNodes2->buildSubstractionOptimized(_node_reduction));
  appendVertices(diff,pflNodes2);
}

MEDMeshMultiLev *MEDUMeshMultiLev::prepare() const
{
  return new MEDUMeshMultiLev(*this);
}

MEDUMeshMultiLev::MEDUMeshMultiLev(const MEDUMeshMultiLev& other):MEDMeshMultiLev(other),_parts(other._parts),_coords(other._coords)
{
}

MEDUMeshMultiLev::MEDUMeshMultiLev(const MEDStructuredMeshMultiLev& other, const MCAuto<MEDCoupling1GTUMesh>& part):MEDMeshMultiLev(other)
{
  _parts.resize(1);
  _parts[0]=part;
  _geo_types.resize(1); _geo_types[0]=part->getCellModelEnum();
  _nb_entities.resize(1); _nb_entities[0]=part->getNumberOfCells();
  _pfls.resize(1); _pfls[0]=0;
}

/*! 
 * To be called only once ! Because due to some optimizations (sometimes aggressive) the internal state can be changed...
 * If returned value is false output pointer \a coords is not the internal pointer. If returned value is true output pointer \a coords is directly the internal pointer.
 * If true is returned, the \a coords output parameter should be used with care (non const method call) to avoid to change the internal state of MEDFileUMesh instance.
 */
bool MEDUMeshMultiLev::buildVTUArrays(DataArrayDouble *& coords, DataArrayByte *&types, DataArrayInt *&cellLocations, DataArrayInt *& cells, DataArrayInt *&faceLocations, DataArrayInt *&faces) const
{
  const DataArrayDouble *tmp(0);
  if(_parts.empty())
    tmp=_coords;
  else
    tmp=_parts[0]->getCoords();
  if(!tmp)
    throw INTERP_KERNEL::Exception("MEDUMeshMultiLev::getVTUArrays : the coordinates are null !");
  MCAuto<DataArrayDouble> a(const_cast<DataArrayDouble *>(tmp)); tmp->incrRef();
  int szBCE(0),szD(0),szF(0);
  bool isPolyh(false);
  int iii(0);
  for(std::vector< MCAuto<MEDCoupling1GTUMesh> >::const_iterator it=_parts.begin();it!=_parts.end();it++,iii++)
    {
      const MEDCoupling1GTUMesh *cur(*it);
      if(!cur)
        throw INTERP_KERNEL::Exception("MEDUMeshMultiLev::getVTUArrays : a part is null !");
      //
      const DataArrayInt *pfl(_pfls[iii]);
      MCAuto<MEDCoupling1GTUMesh> cur2;
      if(!pfl)
        { cur2=const_cast<MEDCoupling1GTUMesh *>(cur); cur2->incrRef(); }
      else
        { cur2=dynamic_cast<MEDCoupling1GTUMesh *>(cur->buildPartOfMySelfKeepCoords(pfl->begin(),pfl->end())); cur=cur2; }
      //
      int curNbCells(cur->getNumberOfCells());
      szBCE+=curNbCells;
      if((*it)->getCellModelEnum()!=INTERP_KERNEL::NORM_POLYHED)
        szD+=cur->getNodalConnectivity()->getNumberOfTuples()+curNbCells;
      else
        {
          isPolyh=true;
          MCAuto<DataArrayInt> tmp2(cur->computeEffectiveNbOfNodesPerCell());
          szD+=tmp2->accumulate(0)+curNbCells;
          szF+=2*curNbCells+cur->getNodalConnectivity()->getNumberOfTuples();
        }
    }
  MCAuto<DataArrayByte> b(DataArrayByte::New()); b->alloc(szBCE,1); char *bPtr(b->getPointer());
  MCAuto<DataArrayInt> c(DataArrayInt::New()); c->alloc(szBCE,1); int *cPtr(c->getPointer());
  MCAuto<DataArrayInt> d(DataArrayInt::New()); d->alloc(szD,1); int *dPtr(d->getPointer());
  MCAuto<DataArrayInt> e(DataArrayInt::New()),f(DataArrayInt::New()); int *ePtr(0),*fPtr(0);
  if(isPolyh)
    { e->alloc(szBCE,1); ePtr=e->getPointer(); f->alloc(szF,1); fPtr=f->getPointer(); }
  int k(0);
  iii=0;
  for(std::vector< MCAuto<MEDCoupling1GTUMesh> >::const_iterator it=_parts.begin();it!=_parts.end();it++,iii++)
    {
      const MEDCoupling1GTUMesh *cur(*it);
      //
      const DataArrayInt *pfl(_pfls[iii]);
      MCAuto<MEDCoupling1GTUMesh> cur2;
      if(!pfl)
        { cur2=const_cast<MEDCoupling1GTUMesh *>(cur); cur2->incrRef(); }
      else
        { cur2=dynamic_cast<MEDCoupling1GTUMesh *>(cur->buildPartOfMySelfKeepCoords(pfl->begin(),pfl->end())); cur=cur2; }
      //
      int curNbCells(cur->getNumberOfCells());
      int gt((int)cur->getCellModelEnum());
      if(gt<0 || gt>=PARAMEDMEM_2_VTKTYPE_LGTH)
        throw INTERP_KERNEL::Exception("MEDUMeshMultiLev::getVTUArrays : invalid geometric type !");
      unsigned char gtvtk(PARAMEDMEM_2_VTKTYPE[gt]);
      if(gtvtk==255)
        throw INTERP_KERNEL::Exception("MEDUMeshMultiLev::getVTUArrays : no VTK type for the requested INTERP_KERNEL geometric type !");
      std::fill(bPtr,bPtr+curNbCells,gtvtk); bPtr+=curNbCells;
      const MEDCoupling1SGTUMesh *scur(dynamic_cast<const MEDCoupling1SGTUMesh *>(cur));
      const MEDCoupling1DGTUMesh *dcur(dynamic_cast<const MEDCoupling1DGTUMesh *>(cur));
      const int *connPtr(cur->getNodalConnectivity()->begin());
      if(!scur && !dcur)
        throw INTERP_KERNEL::Exception("MEDUMeshMultiLev::getVTUArrays : internal error !");
      if(scur)
        {
          if(cur->getCellModelEnum()!=INTERP_KERNEL::NORM_HEXA27)
            {
              int nnpc(scur->getNumberOfNodesPerCell());
              for(int i=0;i<curNbCells;i++,connPtr+=nnpc)
                {
                  *dPtr++=nnpc;
                  dPtr=std::copy(connPtr,connPtr+nnpc,dPtr);
                  *cPtr++=k; k+=nnpc+1;
                }
            }
          else
            {
              for(int i=0;i<curNbCells;i++,connPtr+=27)
                {
                  *dPtr++=27;
                  for(int j=0;j<27;j++,dPtr++)
                    *dPtr=connPtr[HEXA27_PERM_ARRAY[j]];
                  *cPtr++=k; k+=28;
                }
            }
          if(isPolyh)
            { std::fill(ePtr,ePtr+curNbCells,-1); ePtr+=curNbCells; }
        }
      else
        {
          const int *connIPtr(dcur->getNodalConnectivityIndex()->begin());
          if(cur->getCellModelEnum()!=INTERP_KERNEL::NORM_POLYHED)
            {
              for(int i=0;i<curNbCells;i++,connIPtr++)
                {
                  *dPtr++=connIPtr[1]-connIPtr[0];
                  dPtr=std::copy(connPtr+connIPtr[0],connPtr+connIPtr[1],dPtr);
                  *cPtr++=k; k+=connIPtr[1]-connIPtr[0]+1;
                }
            }
          else
            {
              for(int i=0;i<curNbCells;i++,connIPtr++)
                {
                  std::set<int> s(connPtr+connIPtr[0],connPtr+connIPtr[1]); s.erase(-1);
                  *dPtr++=(int)s.size();
                  dPtr=std::copy(s.begin(),s.end(),dPtr);
                  *cPtr++=k; k+=(int)s.size()+1;
                }
            }
          if(isPolyh)
            {
              connIPtr=dcur->getNodalConnectivityIndex()->begin();
              if(cur->getCellModelEnum()!=INTERP_KERNEL::NORM_POLYHED)
                { std::fill(ePtr,ePtr+curNbCells,-1); ePtr+=curNbCells; }
              else
                {
                  int kk(0);
                  for(int i=0;i<curNbCells;i++,connIPtr++)
                    {
                      int nbFace(std::count(connPtr+connIPtr[0],connPtr+connIPtr[1],-1)+1);
                      *fPtr++=nbFace;
                      const int *work(connPtr+connIPtr[0]);
                      for(int j=0;j<nbFace;j++)
                        {
                          const int *work2=std::find(work,connPtr+connIPtr[1],-1);
                          *fPtr++=std::distance(work,work2);
                          fPtr=std::copy(work,work2,fPtr);
                          work=work2+1;
                        }
                      *ePtr++=kk; kk+=connIPtr[1]-connIPtr[0]+2;
                    }
                }
            }
        }
    }
  if(!isPolyh)
    reorderNodesIfNecessary(a,d,0);
  else
    reorderNodesIfNecessary(a,d,f);
  if(a->getNumberOfComponents()!=3)
    a=a->changeNbOfComponents(3,0.);
  coords=a.retn(); types=b.retn(); cellLocations=c.retn(); cells=d.retn();
  if(!isPolyh)
    { faceLocations=0; faces=0; }
  else
    { faceLocations=e.retn(); faces=f.retn(); }
  return _mesh->isObjectInTheProgeny(coords);
}

void MEDUMeshMultiLev::reorderNodesIfNecessary(MCAuto<DataArrayDouble>& coords, DataArrayInt *nodalConnVTK, DataArrayInt *polyhedNodalConnVTK) const
{
  const DataArrayInt *nr(_node_reduction);
  if(!nr)
    return ;
  if(nodalConnVTK->empty() && !polyhedNodalConnVTK)
    {
      coords=(coords->selectByTupleIdSafe(nr->begin(),nr->end()));
      return ;
    }
  int sz(coords->getNumberOfTuples());
  std::vector<bool> b(sz,false);
  const int *work(nodalConnVTK->begin()),*endW(nodalConnVTK->end());
  while(work!=endW)
    {
      int nb(*work++);
      for(int i=0;i<nb && work!=endW;i++,work++)
        {
          if(*work>=0 && *work<sz)
            b[*work]=true;
          else
            throw INTERP_KERNEL::Exception("MEDUMeshMultiLev::reorderNodesIfNecessary : internal error !");
        }
    }
  if(polyhedNodalConnVTK)
    {
      work=polyhedNodalConnVTK->begin(); endW=polyhedNodalConnVTK->end();
      while(work!=endW)
        {
          int nb(*work++);
          for(int i=0;i<nb && work!=endW;i++)
            {
              int nb2(*work++);
              for(int j=0;j<nb2 && work!=endW;j++,work++)
                {
                  if(*work>=0 && *work<sz)
                    b[*work]=true;
                  else
                    throw INTERP_KERNEL::Exception("MEDUMeshMultiLev::reorderNodesIfNecessary : internal error #2 !");
                }
            }
        }
    }
  int szExp(std::count(b.begin(),b.end(),true));
  if(szExp!=nr->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("MEDUMeshMultiLev::reorderNodesIfNecessary : internal error #3 !");
  // Go renumbering !
  MCAuto<DataArrayInt> o2n(DataArrayInt::New()); o2n->alloc(sz,1);
  int *o2nPtr(o2n->getPointer());
  int newId(0);
  for(int i=0;i<sz;i++,o2nPtr++)
    if(b[i]) *o2nPtr=newId++; else *o2nPtr=-1;
  const int *o2nPtrc(o2n->begin());
  MCAuto<DataArrayInt> n2o(o2n->invertArrayO2N2N2O(nr->getNumberOfTuples()));
  MCAuto<DataArrayInt> perm(DataArrayInt::FindPermutationFromFirstToSecond(n2o,nr));
  const int *permPtr(perm->begin());
  int *work2(nodalConnVTK->getPointer()),*endW2(nodalConnVTK->getPointer()+nodalConnVTK->getNumberOfTuples());
  while(work2!=endW2)
    {
      int nb(*work2++);
      for(int i=0;i<nb && work2!=endW2;i++,work2++)
        *work2=permPtr[o2nPtrc[*work2]];
    }
  if(polyhedNodalConnVTK)
    {
      work2=polyhedNodalConnVTK->getPointer(); endW2=polyhedNodalConnVTK->getPointer()+polyhedNodalConnVTK->getNumberOfTuples();
      while(work2!=endW2)
        {
          int nb(*work2++);
          for(int i=0;i<nb && work2!=endW2;i++)
            {
              int nb2(*work2++);
              for(int j=0;j<nb2 && work2!=endW2;j++,work2++)
                *work2=permPtr[o2nPtrc[*work2]];
            }
        }
    }
  coords=(coords->selectByTupleIdSafe(nr->begin(),nr->end()));
}


void MEDUMeshMultiLev::appendVertices(const DataArrayInt *verticesToAdd, DataArrayInt *nr)
{
  int nbOfCells(verticesToAdd->getNumberOfTuples());//it is not a bug cells are NORM_POINT1
  MEDMeshMultiLev::appendVertices(verticesToAdd,nr);
  MCAuto<MEDCoupling1SGTUMesh> elt(MEDCoupling1SGTUMesh::New("",INTERP_KERNEL::NORM_POINT1));
  elt->allocateCells(nbOfCells);
  for(int i=0;i<nbOfCells;i++)
    {
      int pt(verticesToAdd->getIJ(i,0));
      elt->insertNextCell(&pt,&pt+1);
    }
  if(_parts.empty())
    throw INTERP_KERNEL::Exception("MEDUMeshMultiLev::appendVertices : parts are empty !");
  elt->setCoords(_parts[0]->getCoords());
  MCAuto<MEDCoupling1GTUMesh> elt2((MEDCoupling1SGTUMesh *)elt); elt2->incrRef();
  _parts.push_back(elt2);
}

//=

MEDStructuredMeshMultiLev::MEDStructuredMeshMultiLev(const MEDFileStructuredMesh *m, const std::vector<int>& lev):MEDMeshMultiLev(m),_is_internal(true)
{
  initStdFieldOfIntegers(m);
}

MEDStructuredMeshMultiLev::MEDStructuredMeshMultiLev(const MEDFileStructuredMesh *m, int nbOfNodes, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayInt *>& pfls, const std::vector<int>& nbEntities):MEDMeshMultiLev(m,nbOfNodes,gts,pfls,nbEntities),_is_internal(true)
{
  initStdFieldOfIntegers(m);
}

MEDStructuredMeshMultiLev::MEDStructuredMeshMultiLev(const MEDStructuredMeshMultiLev& other):MEDMeshMultiLev(other),_is_internal(true),_face_fam_ids(other._face_fam_ids),_face_num_ids(other._face_num_ids)
{
}

void MEDStructuredMeshMultiLev::initStdFieldOfIntegers(const MEDFileStructuredMesh *m)
{
  // ids fields management
  const DataArrayInt *tmp(0);
  tmp=m->getFamilyFieldAtLevel(0);
  if(tmp)
    {
      tmp->incrRef();
      _cell_fam_ids=const_cast<DataArrayInt *>(tmp);
    }
  tmp=m->getNumberFieldAtLevel(0);
  if(tmp)
    {
      tmp->incrRef();
      _cell_num_ids=const_cast<DataArrayInt *>(tmp);
    }
  //
  tmp=0;
  tmp=m->getFamilyFieldAtLevel(1);
  if(tmp)
    {
      tmp->incrRef();
      _node_fam_ids=const_cast<DataArrayInt *>(tmp);
    }
  tmp=m->getNumberFieldAtLevel(1);
  if(tmp)
    {
      tmp->incrRef();
      _node_num_ids=const_cast<DataArrayInt *>(tmp);
    }
  // faces (if any)
  tmp=m->getFamilyFieldAtLevel(-1);
  if(tmp)
    {
      tmp->incrRef();
      _face_fam_ids=const_cast<DataArrayInt *>(tmp);
    }
  tmp=m->getNumberFieldAtLevel(-1);
  if(tmp)
    {
      tmp->incrRef();
      _face_num_ids=const_cast<DataArrayInt *>(tmp);
    }
}

void MEDStructuredMeshMultiLev::moveFaceToCell() const
{
  const_cast<MEDStructuredMeshMultiLev *>(this)->_cell_fam_ids=_face_fam_ids; const_cast<MEDStructuredMeshMultiLev *>(this)->_face_fam_ids=0;
  const_cast<MEDStructuredMeshMultiLev *>(this)->_cell_num_ids=_face_num_ids; const_cast<MEDStructuredMeshMultiLev *>(this)->_face_num_ids=0;
}

bool MEDStructuredMeshMultiLev::prepareForImplicitUnstructuredMeshCase(MEDMeshMultiLev *&ret) const
{
  ret=0;
  if(_geo_types.empty())
    return false;
  if(_geo_types.size()!=1)
    throw INTERP_KERNEL::Exception("MEDStructuredMeshMultiLev::prepareForImplicitUnstructuredMeshCase only one geo types supported at most supported for the moment !");
  INTERP_KERNEL::NormalizedCellType gt(MEDCouplingStructuredMesh::GetGeoTypeGivenMeshDimension(_mesh->getMeshDimension()));
  if(_geo_types[0]==gt)
    return false;
  MEDCoupling1GTUMesh *facesIfPresent((static_cast<const MEDFileStructuredMesh *>(_mesh))->getImplicitFaceMesh());
  if(!facesIfPresent)
    return false;
  const DataArrayInt *pfl(0),*nr(_node_reduction);
  if(!_pfls.empty())
    pfl=_pfls[0];
  MCAuto<MEDCoupling1GTUMesh> facesIfPresent2(facesIfPresent); facesIfPresent->incrRef();
  moveFaceToCell();
  MCAuto<MEDUMeshMultiLev> ret2(new MEDUMeshMultiLev(*this,facesIfPresent2));
  if(pfl)
    ret2->setCellReduction(pfl);
  if(nr)
    throw INTERP_KERNEL::Exception("MEDStructuredMeshMultiLev::prepareForImplicitUnstructuredMeshCase : case is not treated yet for node reduction on implicit unstructured mesh.");
  ret=ret2.retn();
  return true;
}

void MEDStructuredMeshMultiLev::dealWithImplicitUnstructuredMesh(const MEDFileMesh *m)
{
  const DataArrayInt *tmp(0);
  tmp=m->getFamilyFieldAtLevel(-1);
  if(tmp)
    {
      tmp->incrRef();
      _cell_fam_ids=const_cast<DataArrayInt *>(tmp);
    }
  tmp=m->getNumberFieldAtLevel(-1);
  if(tmp)
    {
      tmp->incrRef();
      _cell_num_ids=const_cast<DataArrayInt *>(tmp);
    }
}

void MEDStructuredMeshMultiLev::selectPartOfNodes(const DataArrayInt *pflNodes)
{
  if(!pflNodes || !pflNodes->isAllocated())
    return ;
  std::vector<int> ngs(getNodeGridStructure());
  MCAuto<DataArrayInt> conn(MEDCouplingStructuredMesh::Build1GTNodalConnectivity(&ngs[0],&ngs[0]+ngs.size()));
  MCAuto<MEDCoupling1SGTUMesh> m(MEDCoupling1SGTUMesh::New("",MEDCouplingStructuredMesh::GetGeoTypeGivenMeshDimension(ngs.size())));
  m->setNodalConnectivity(conn);
  const DataArrayInt *pfl(_pfls[0]);
  if(pfl)
    {
      m=dynamic_cast<MEDCoupling1SGTUMesh *>(m->buildPartOfMySelfKeepCoords(pfl->begin(),pfl->end()));
    }
  DataArrayInt *cellIds=0;
  m->fillCellIdsToKeepFromNodeIds(pflNodes->begin(),pflNodes->end(),true,cellIds);
  MCAuto<DataArrayInt> cellIdsSafe(cellIds);
  MCAuto<MEDCouplingPointSet> m2(m->buildPartOfMySelfKeepCoords(cellIds->begin(),cellIds->end()));
  int tmp=-1;
  _node_reduction=m2->getNodeIdsInUse(tmp);
  if(pfl)
    _pfls[0]=pfl->selectByTupleIdSafe(cellIds->begin(),cellIds->end());
  else
    _pfls[0]=cellIdsSafe;
}

//=

MEDCMeshMultiLev *MEDCMeshMultiLev::New(const MEDFileCMesh *m, const std::vector<int>& levs)
{
  return new MEDCMeshMultiLev(m,levs);
}

MEDCMeshMultiLev *MEDCMeshMultiLev::New(const MEDFileCMesh *m, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayInt *>& pfls, const std::vector<int>& nbEntities)
{
  return new MEDCMeshMultiLev(m,gts,pfls,nbEntities);
}

MEDCMeshMultiLev::MEDCMeshMultiLev(const MEDFileCMesh *m, const std::vector<int>& levs):MEDStructuredMeshMultiLev(m,levs)
{
  if(!m)
    throw INTERP_KERNEL::Exception("MEDCMeshMultiLev constructor : null input pointer !");
  if(levs.size()!=1 || levs[0]!=0)
    throw INTERP_KERNEL::Exception("MEDCMeshMultiLev constructor : levels supported is 0 only !");
  int sdim(m->getSpaceDimension());
  _coords.resize(sdim);
  for(int i=0;i<sdim;i++)
    {
      DataArrayDouble *elt(const_cast<DataArrayDouble *>(m->getMesh()->getCoordsAt(i)));
      if(!elt)
        throw INTERP_KERNEL::Exception("MEDCMeshMultiLev constructor 2 : presence of null pointer for an vector of double along an axis !");
      elt->incrRef();
      _coords[i]=elt;
    }
}

MEDCMeshMultiLev::MEDCMeshMultiLev(const MEDFileCMesh *m, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayInt *>& pfls, const std::vector<int>& nbEntities):MEDStructuredMeshMultiLev(m,m->getNumberOfNodes(),gts,pfls,nbEntities)
{
  if(!m)
    throw INTERP_KERNEL::Exception("MEDCMeshMultiLev constructor 2 : null input pointer !");
  if(gts.size()!=1 || pfls.size()!=1)
    throw INTERP_KERNEL::Exception("MEDCMeshMultiLev constructor 2 : lengthes of gts and pfls must be equal to one !");
  int mdim(m->getMeshDimension());
  INTERP_KERNEL::NormalizedCellType gt(MEDCouplingStructuredMesh::GetGeoTypeGivenMeshDimension(mdim));
  if(gt==gts[0])
    {
      _coords.resize(mdim);
      for(int i=0;i<mdim;i++)
        {
          DataArrayDouble *elt(const_cast<DataArrayDouble *>(m->getMesh()->getCoordsAt(i)));
          if(!elt)
            throw INTERP_KERNEL::Exception("MEDCMeshMultiLev constructor 2 : presence of null pointer for an vector of double along an axis !");
          _coords[i]=elt; _coords[i]->incrRef();
        }
    }
  else
    dealWithImplicitUnstructuredMesh(m);
}

MEDCMeshMultiLev::MEDCMeshMultiLev(const MEDCMeshMultiLev& other):MEDStructuredMeshMultiLev(other),_coords(other._coords)
{
}

std::vector<int> MEDCMeshMultiLev::getNodeGridStructure() const
{
  std::vector<int> ret(_coords.size());
  for(std::size_t i=0;i<_coords.size();i++)
    ret[i]=_coords[i]->getNumberOfTuples();
  return ret;
}

MEDMeshMultiLev *MEDCMeshMultiLev::prepare() const
{
  MEDMeshMultiLev *retSpecific(0);
  if(prepareForImplicitUnstructuredMeshCase(retSpecific))
    return retSpecific;
  const DataArrayInt *pfl(0),*nr(_node_reduction);
  if(!_pfls.empty())
    pfl=_pfls[0];
  MCAuto<DataArrayInt> nnr;
  std::vector<int> cgs,ngs(getNodeGridStructure());
  cgs.resize(ngs.size());
  std::transform(ngs.begin(),ngs.end(),cgs.begin(),std::bind2nd(std::plus<int>(),-1));
  if(pfl)
    {
      std::vector< std::pair<int,int> > cellParts;
      MCAuto<MEDMeshMultiLev> ret2;
      if(MEDCouplingStructuredMesh::IsPartStructured(pfl->begin(),pfl->end(),cgs,cellParts))
        {
          MCAuto<MEDCMeshMultiLev> ret(new MEDCMeshMultiLev(*this));
          ret->_is_internal=false;
          if(nr)
            { nnr=nr->deepCopy(); nnr->sort(true); ret->setNodeReduction(nnr); }
          ret->_nb_entities[0]=pfl->getNumberOfTuples();
          ret->_pfls[0]=0;
          std::vector< MCAuto<DataArrayDouble> > coords(_coords.size());
          for(std::size_t i=0;i<_coords.size();i++)
            coords[i]=_coords[i]->selectByTupleIdSafeSlice(cellParts[i].first,cellParts[i].second+1,1);
          ret->_coords=coords;
          ret2=(MEDCMeshMultiLev *)ret; ret2->incrRef();
        }
      else
        {
          MCAuto<MEDCouplingCMesh> m(MEDCouplingCMesh::New());
          for(std::size_t i=0;i<ngs.size();i++)
            m->setCoordsAt(i,_coords[i]);
          MCAuto<MEDCoupling1SGTUMesh> m2(m->build1SGTUnstructured());
          MCAuto<MEDCoupling1GTUMesh> m3=dynamic_cast<MEDCoupling1GTUMesh *>(m2->buildPartOfMySelfKeepCoords(pfl->begin(),pfl->end()));
          MCAuto<MEDUMeshMultiLev> ret(new MEDUMeshMultiLev(*this,m3));
          if(nr)
            { m3->zipCoords(); nnr=nr->deepCopy(); nnr->sort(true); ret->setNodeReduction(nnr); }
          ret2=(MEDUMeshMultiLev *)ret; ret2->incrRef();
        }
      const DataArrayInt *famIds(_cell_fam_ids),*numIds(_cell_num_ids);
      if(famIds)
        {
          MCAuto<DataArrayInt> tmp(famIds->selectByTupleIdSafe(pfl->begin(),pfl->end()));
          ret2->setFamilyIdsOnCells(tmp);
        }
      if(numIds)
        {
          MCAuto<DataArrayInt> tmp(numIds->selectByTupleIdSafe(pfl->begin(),pfl->end()));
          ret2->setNumberIdsOnCells(tmp);
        }
      return ret2.retn();

    }
  else
    {
      MCAuto<MEDCMeshMultiLev> ret(new MEDCMeshMultiLev(*this));
      if(nr)
        { nnr=nr->deepCopy(); nnr->sort(true); ret->setNodeReduction(nnr); }
      return ret.retn();
    }
}

/*!
 * \a param [out] isInternal if true the returned pointers are those in main data structure. If false those pointers have been built espacially for that method.
 */
std::vector< DataArrayDouble * > MEDCMeshMultiLev::buildVTUArrays(bool& isInternal) const
{
  isInternal=_is_internal;
  std::size_t sz(_coords.size());
  std::vector< DataArrayDouble * > ret(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      ret[i]=const_cast<DataArrayDouble *>((const DataArrayDouble *)_coords[i]);
      ret[i]->incrRef();
    }
  return ret;
}

//=

MEDCurveLinearMeshMultiLev *MEDCurveLinearMeshMultiLev::New(const MEDFileCurveLinearMesh *m, const std::vector<int>& levs)
{
  return new MEDCurveLinearMeshMultiLev(m,levs);
}

MEDCurveLinearMeshMultiLev *MEDCurveLinearMeshMultiLev::New(const MEDFileCurveLinearMesh *m, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayInt *>& pfls, const std::vector<int>& nbEntities)
{
  return new MEDCurveLinearMeshMultiLev(m,gts,pfls,nbEntities);
}

MEDCurveLinearMeshMultiLev::MEDCurveLinearMeshMultiLev(const MEDFileCurveLinearMesh *m, const std::vector<int>& levs):MEDStructuredMeshMultiLev(m,levs)
{
  if(!m)
    throw INTERP_KERNEL::Exception("MEDCurveLinearMeshMultiLev constructor : null input pointer !");
  if(levs.size()!=1 || levs[0]!=0)
    throw INTERP_KERNEL::Exception("MEDCurveLinearMeshMultiLev constructor : levels supported is 0 only !");
  DataArrayDouble *coords(const_cast<DataArrayDouble *>(m->getMesh()->getCoords()));
  if(!coords)
    throw INTERP_KERNEL::Exception("MEDCurveLinearMeshMultiLev constructor 2 : no coords set !");
  coords->incrRef();
  _coords=coords;
  _structure=m->getMesh()->getNodeGridStructure();
}

MEDCurveLinearMeshMultiLev::MEDCurveLinearMeshMultiLev(const MEDFileCurveLinearMesh *m, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayInt *>& pfls, const std::vector<int>& nbEntities):MEDStructuredMeshMultiLev(m,m->getNumberOfNodes(),gts,pfls,nbEntities)
{
  if(!m)
    throw INTERP_KERNEL::Exception("MEDCurveLinearMeshMultiLev constructor 2 : null input pointer !");
  if(gts.size()!=1 || pfls.size()!=1)
    throw INTERP_KERNEL::Exception("MEDCurveLinearMeshMultiLev constructor 2 : lengthes of gts and pfls must be equal to one !");
  INTERP_KERNEL::NormalizedCellType gt(MEDCouplingStructuredMesh::GetGeoTypeGivenMeshDimension(m->getMeshDimension()));
  if(gt==gts[0])
    {
      DataArrayDouble *coords(const_cast<DataArrayDouble *>(m->getMesh()->getCoords()));
      if(!coords)
        throw INTERP_KERNEL::Exception("MEDCurveLinearMeshMultiLev constructor 2 : no coords set !");
      coords->incrRef();
      _coords=coords;
      _structure=m->getMesh()->getNodeGridStructure();
    }
  else
    dealWithImplicitUnstructuredMesh(m);
}

MEDCurveLinearMeshMultiLev::MEDCurveLinearMeshMultiLev(const MEDCurveLinearMeshMultiLev& other):MEDStructuredMeshMultiLev(other),_coords(other._coords),_structure(other._structure)
{
}

std::vector<int> MEDCurveLinearMeshMultiLev::getNodeGridStructure() const
{
  return _structure;
}

MEDMeshMultiLev *MEDCurveLinearMeshMultiLev::prepare() const
{
  MEDMeshMultiLev *retSpecific(0);
  if(prepareForImplicitUnstructuredMeshCase(retSpecific))
    return retSpecific;
  const DataArrayInt *pfl(0),*nr(_node_reduction);
  if(!_pfls.empty())
    pfl=_pfls[0];
  MCAuto<DataArrayInt> nnr;
  std::vector<int> cgs,ngs(getNodeGridStructure());
  cgs.resize(ngs.size());
  std::transform(ngs.begin(),ngs.end(),cgs.begin(),std::bind2nd(std::plus<int>(),-1));
  if(pfl)
    {
      std::vector< std::pair<int,int> > cellParts,nodeParts;
      MCAuto<MEDMeshMultiLev> ret2;
      if(MEDCouplingStructuredMesh::IsPartStructured(pfl->begin(),pfl->end(),cgs,cellParts))
        {
          nodeParts=cellParts;
          std::vector<int> st(ngs.size());
          for(std::size_t i=0;i<ngs.size();i++)
            {
              nodeParts[i].second++;
              st[i]=nodeParts[i].second-nodeParts[i].first;
            }
          MCAuto<DataArrayInt> p(MEDCouplingStructuredMesh::BuildExplicitIdsFrom(ngs,nodeParts));
          MCAuto<MEDCurveLinearMeshMultiLev> ret(new MEDCurveLinearMeshMultiLev(*this));
          ret->_is_internal=false;
          if(nr)
            { nnr=nr->deepCopy(); nnr->sort(true); ret->setNodeReduction(nnr); }
          ret->_nb_entities[0]=pfl->getNumberOfTuples();
          ret->_pfls[0]=0;
          ret->_coords=_coords->selectByTupleIdSafe(p->begin(),p->end());
          ret->_structure=st;
          ret2=(MEDCurveLinearMeshMultiLev *)ret; ret2->incrRef();
        }
      else
        {
          MCAuto<MEDCouplingCurveLinearMesh> m(MEDCouplingCurveLinearMesh::New());
          m->setCoords(_coords); m->setNodeGridStructure(&_structure[0],&_structure[0]+_structure.size());
          MCAuto<MEDCoupling1SGTUMesh> m2(m->build1SGTUnstructured());
          MCAuto<MEDCoupling1GTUMesh> m3=dynamic_cast<MEDCoupling1GTUMesh *>(m2->buildPartOfMySelfKeepCoords(pfl->begin(),pfl->end()));
          MCAuto<MEDUMeshMultiLev> ret(new MEDUMeshMultiLev(*this,m3));
          if(nr)
            { m3->zipCoords(); nnr=nr->deepCopy(); nnr->sort(true); ret->setNodeReduction(nnr); }
          ret2=(MEDUMeshMultiLev *)ret; ret2->incrRef();
        }
      const DataArrayInt *famIds(_cell_fam_ids),*numIds(_cell_num_ids);
      if(famIds)
        {
          MCAuto<DataArrayInt> tmp(famIds->selectByTupleIdSafe(pfl->begin(),pfl->end()));
          ret2->setFamilyIdsOnCells(tmp);
        }
      if(numIds)
        {
          MCAuto<DataArrayInt> tmp(numIds->selectByTupleIdSafe(pfl->begin(),pfl->end()));
          ret2->setNumberIdsOnCells(tmp);
        }
      return ret2.retn();
    }
  else
    {
      MCAuto<MEDCurveLinearMeshMultiLev> ret(new MEDCurveLinearMeshMultiLev(*this));
      if(nr)
        { nnr=nr->deepCopy(); nnr->sort(true); ret->setNodeReduction(nnr); }
      return ret.retn();
    }
}

void MEDCurveLinearMeshMultiLev::buildVTUArrays(DataArrayDouble *&coords, std::vector<int>& nodeStrct, bool& isInternal) const
{
  isInternal=_is_internal;
  nodeStrct=_structure;
  const DataArrayDouble *coo(_coords);
  if(!coo)
    throw INTERP_KERNEL::Exception("MEDCurveLinearMeshMultiLev::buildVTUArrays : null pointer on coordinates !");
  coords=const_cast<DataArrayDouble *>(coo); coords->incrRef();
}

//=

MEDFileField1TSStructItem2::MEDFileField1TSStructItem2()
{
}

MEDFileField1TSStructItem2::MEDFileField1TSStructItem2(INTERP_KERNEL::NormalizedCellType a, const std::pair<int,int>& b, const std::string& c, const std::string& d):_geo_type(a),_start_end(b),_pfl(DataArrayInt::New()),_loc(d),_nb_of_entity(-1)
{
  _pfl->setName(c.c_str());
}

void MEDFileField1TSStructItem2::checkWithMeshStructForCells(const MEDFileMeshStruct *mst, const MEDFileFieldGlobsReal *globs)
{
  if(!mst->doesManageGeoType(_geo_type))
    {
      MEDFileMeshStruct *mstUnConstCasted(const_cast<MEDFileMeshStruct *>(mst));
      mstUnConstCasted->appendIfImplicitType(_geo_type);
    }
  int nbOfEnt=mst->getNumberOfElemsOfGeoType(_geo_type);
  checkInRange(nbOfEnt,1,globs);
}

void MEDFileField1TSStructItem2::checkWithMeshStructForGaussNE(const MEDFileMeshStruct *mst, const MEDFileFieldGlobsReal *globs)
{
  int nbOfEnt=mst->getNumberOfElemsOfGeoType(_geo_type);
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(_geo_type);
  checkInRange(nbOfEnt,(int)cm.getNumberOfNodes(),globs);
}

void MEDFileField1TSStructItem2::checkWithMeshStructForGaussPT(const MEDFileMeshStruct *mst, const MEDFileFieldGlobsReal *globs)
{
  if(!globs)
    throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem2::checkWithMeshStructForGaussPT : no globals specified !");
  if(_loc.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem2::checkWithMeshStructForGaussPT : no localization specified !");
  const MEDFileFieldLoc& loc=globs->getLocalization(_loc.c_str());
  int nbOfEnt=mst->getNumberOfElemsOfGeoType(_geo_type);
  checkInRange(nbOfEnt,loc.getNumberOfGaussPoints(),globs);
}

int MEDFileField1TSStructItem2::getNbOfIntegrationPts(const MEDFileFieldGlobsReal *globs) const
{
  if(_loc.empty())
    {
      if(getPflName().empty())
        return (_start_end.second-_start_end.first)/_nb_of_entity;
      else
        return (_start_end.second-_start_end.first)/getPfl(globs)->getNumberOfTuples();
    }
  else
    {
      const MEDFileFieldLoc& loc(globs->getLocalization(_loc.c_str()));
      return loc.getNumberOfGaussPoints();
    }
}

std::string MEDFileField1TSStructItem2::getPflName() const
{
  return _pfl->getName();
}

const DataArrayInt *MEDFileField1TSStructItem2::getPfl(const MEDFileFieldGlobsReal *globs) const
{
  if(!_pfl->isAllocated())
    {
      if(_pfl->getName().empty())
        return 0;
      else
        return globs->getProfile(_pfl->getName().c_str());
    }
  else
    return _pfl;
}

/*!
 * \param [in] nbOfEntity - number of entity that can be either cells or nodes. Not other possiblity.
 * \param [in] nip - number of integration points. 1 for ON_CELLS and NO_NODES
 */
void MEDFileField1TSStructItem2::checkInRange(int nbOfEntity, int nip, const MEDFileFieldGlobsReal *globs)
{
  _nb_of_entity=nbOfEntity;
  if(_pfl->getName().empty())
    {
      if(nbOfEntity!=(_start_end.second-_start_end.first)/nip)
        throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem2::checkInRange : Mismatch between number of entities and size of field !");
      return ;
    }
  else
    {
      if(!globs)
        throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem2::checkInRange : Presence of a profile on field whereas no globals found in file !");
      const DataArrayInt *pfl=globs->getProfile(_pfl->getName().c_str());
      if(!pfl)
        throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem2::checkInRange : Presence of a profile on field whereas no such profile found in file !");
      pfl->checkAllIdsInRange(0,nbOfEntity);
    }
}

bool MEDFileField1TSStructItem2::isFastlyEqual(int& startExp, INTERP_KERNEL::NormalizedCellType gt, const std::string& pflName) const
{
  if(startExp!=_start_end.first)
    return false;
  if(gt!=_geo_type)
    return false;
  if(getPflName()!=pflName)
    return false;
  startExp=_start_end.second;
  return true;
}

bool MEDFileField1TSStructItem2::operator==(const MEDFileField1TSStructItem2& other) const
{
  //_nb_of_entity is not taken into account here. It is not a bug, because no mesh consideration needed here to perform fast compare.
  //idem for _loc. It is not an effective attribute for support comparison.
  return _geo_type==other._geo_type && _start_end==other._start_end && _pfl->getName()==other._pfl->getName();
}

bool MEDFileField1TSStructItem2::isCellSupportEqual(const MEDFileField1TSStructItem2& other, const MEDFileFieldGlobsReal *globs) const
{
  if(_geo_type!=other._geo_type)
    return false;
  if(_nb_of_entity!=other._nb_of_entity)
    return false;
  if((_pfl->getName().empty() && !other._pfl->getName().empty()) || (!_pfl->getName().empty() && other._pfl->getName().empty()))
    return false;
  if(_pfl->getName().empty() && other._pfl->getName().empty())
    return true;
  const DataArrayInt *pfl1(getPfl(globs)),*pfl2(other.getPfl(globs));
  return pfl1->isEqualWithoutConsideringStr(*pfl2);
}

bool MEDFileField1TSStructItem2::isNodeSupportEqual(const MEDFileField1TSStructItem2& other, const MEDFileFieldGlobsReal *globs) const
{
  return isCellSupportEqual(other,globs);
}

/*!
 * \a objs must be non empty. \a objs should contain items having same geometric type.
 */
MEDFileField1TSStructItem2 MEDFileField1TSStructItem2::BuildAggregationOf(const std::vector<const MEDFileField1TSStructItem2 *>& objs, const MEDFileFieldGlobsReal *globs)
{
  if(objs.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem2::BuildAggregationOf : empty input !");
  if(objs.size()==1)
    return MEDFileField1TSStructItem2(*objs[0]);
  INTERP_KERNEL::NormalizedCellType gt(objs[0]->_geo_type);
  int nbEntityRef(objs[0]->_nb_of_entity);
  std::size_t sz(objs.size());
  std::vector<const DataArrayInt *> arrs(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      const MEDFileField1TSStructItem2 *obj(objs[i]);
      if(gt!=obj->_geo_type)
        throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem2::BuildAggregationOf : invalid situation ! All input must have the same geo type !");
      if(nbEntityRef!=obj->_nb_of_entity)
        throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem2::BuildAggregationOf : invalid situation ! All input must have the global nb of entity !");
      if(obj->_pfl->getName().empty())
        throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem2::BuildAggregationOf : invalid situation ! Several same geo type chunk must all lie on profiles !");
      arrs[i]=globs->getProfile(obj->_pfl->getName().c_str());
    }
  MCAuto<DataArrayInt> arr(DataArrayInt::Aggregate(arrs));
  arr->sort();
  int oldNbTuples(arr->getNumberOfTuples());
  arr=arr->buildUnique();
  if(oldNbTuples!=arr->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem2::BuildAggregationOf : some entities are present several times !");
  if(arr->isIota(nbEntityRef))
    {
      std::pair<int,int> p(0,nbEntityRef);
      std::string a,b;
      MEDFileField1TSStructItem2 ret(gt,p,a,b);
      ret._nb_of_entity=nbEntityRef;
      return ret;
    }
  else
    {
      arr->setName(NEWLY_CREATED_PFL_NAME);
      std::pair<int,int> p(0,oldNbTuples);
      std::string a,b;
      MEDFileField1TSStructItem2 ret(gt,p,a,b);
      ret._nb_of_entity=nbEntityRef;
      ret._pfl=arr;
      return ret;
    }
}

std::size_t MEDFileField1TSStructItem2::getHeapMemorySizeWithoutChildren() const
{
  std::size_t ret(_loc.capacity());
  return ret;
}

std::vector<const BigMemoryObject *> MEDFileField1TSStructItem2::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  ret.push_back((const DataArrayInt *)_pfl);
  return ret;
}

//=

MEDFileField1TSStructItem::MEDFileField1TSStructItem(TypeOfField a, const std::vector< MEDFileField1TSStructItem2 >& b):_computed(false),_type(a),_items(b)
{
}

void MEDFileField1TSStructItem::checkWithMeshStruct(const MEDFileMeshStruct *mst, const MEDFileFieldGlobsReal *globs)
{
  switch(_type)
  {
    case ON_NODES:
      {
        int nbOfEnt=mst->getNumberOfNodes();
        if(_items.size()!=1)
          throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem::checkWithMeshStruct : for nodes field only one subdivision supported !");
        _items[0].checkInRange(nbOfEnt,1,globs);
        break ;
      }
    case ON_CELLS:
      {
        for(std::vector< MEDFileField1TSStructItem2 >::iterator it=_items.begin();it!=_items.end();it++)
          (*it).checkWithMeshStructForCells(mst,globs);
        break;
      }
    case ON_GAUSS_NE:
      {
        for(std::vector< MEDFileField1TSStructItem2 >::iterator it=_items.begin();it!=_items.end();it++)
          (*it).checkWithMeshStructForGaussNE(mst,globs);
        break;
      }
    case ON_GAUSS_PT:
      {
        for(std::vector< MEDFileField1TSStructItem2 >::iterator it=_items.begin();it!=_items.end();it++)
          (*it).checkWithMeshStructForGaussPT(mst,globs);
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem::checkWithMeshStruct : not managed field type !");
  }
}

bool MEDFileField1TSStructItem::operator==(const MEDFileField1TSStructItem& other) const
{
  if(_type!=other._type)
    return false;
  if(_items.size()!=other._items.size())
    return false;
  for(std::size_t i=0;i<_items.size();i++)
    if(!(_items[i]==other._items[i]))
      return false;
  return true;
}

bool MEDFileField1TSStructItem::isCellSupportEqual(const MEDFileField1TSStructItem& other, const MEDFileFieldGlobsReal *globs) const
{
  if(_type!=other._type)
    return false;
  if(_items.size()!=other._items.size())
    return false;
  for(std::size_t i=0;i<_items.size();i++)
    if(!(_items[i].isCellSupportEqual(other._items[i],globs)))
      return false;
  return true;
}

bool MEDFileField1TSStructItem::isNodeSupportEqual(const MEDFileField1TSStructItem& other, const MEDFileFieldGlobsReal *globs) const
{
  if(_type!=other._type)
    return false;
  if(_items.size()!=other._items.size())
    return false;
  for(std::size_t i=0;i<_items.size();i++)
    if(!(_items[i].isNodeSupportEqual(other._items[i],globs)))
      return false;
  return true;
}

bool MEDFileField1TSStructItem::isEntityCell() const
{
  if(_type==ON_NODES)
    return false;
  else
    return true;
}

class CmpGeo
{
public:
  CmpGeo(INTERP_KERNEL::NormalizedCellType geoTyp):_geo_type(geoTyp) { }
  bool operator()(const std::pair< INTERP_KERNEL::NormalizedCellType, std::vector<std::size_t> > & v) const { return _geo_type==v.first; }
private:
  INTERP_KERNEL::NormalizedCellType _geo_type;
};

MEDFileField1TSStructItem MEDFileField1TSStructItem::simplifyMeOnCellEntity(const MEDFileFieldGlobsReal *globs) const
{
  if(!isEntityCell())
    throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem::simplifyMeOnCellEntity : must be on ON_CELLS, ON_GAUSS_NE or ON_GAUSS_PT !");
  std::vector< std::pair< INTERP_KERNEL::NormalizedCellType, std::vector<std::size_t> > > m;
  std::size_t i=0;
  for(std::vector< MEDFileField1TSStructItem2 >::const_iterator it=_items.begin();it!=_items.end();it++,i++)
    {
      std::vector< std::pair< INTERP_KERNEL::NormalizedCellType, std::vector<std::size_t> > >::iterator it0(std::find_if(m.begin(),m.end(),CmpGeo((*it).getGeo())));
      if(it0==m.end())
        m.push_back(std::pair< INTERP_KERNEL::NormalizedCellType, std::vector<std::size_t> >((*it).getGeo(),std::vector<std::size_t>(1,i)));
      else
        (*it0).second.push_back(i);
    }
  if(m.size()==_items.size())
    {
      MEDFileField1TSStructItem ret(*this);
      ret._type=ON_CELLS;
      return ret;
    }
  std::size_t sz(m.size());
  std::vector< MEDFileField1TSStructItem2 > items(sz);
  for(i=0;i<sz;i++)
    {
      const std::vector<std::size_t>& ids=m[i].second;
      std::vector<const MEDFileField1TSStructItem2 *>objs(ids.size());
      for(std::size_t j=0;j<ids.size();j++)
        objs[j]=&_items[ids[j]];
      items[i]=MEDFileField1TSStructItem2::BuildAggregationOf(objs,globs);
    }
  MEDFileField1TSStructItem ret(ON_CELLS,items);
  ret._computed=true;
  return ret;
}

/*!
 * \a this is expected to be ON_CELLS and simplified.
 */
bool MEDFileField1TSStructItem::isCompatibleWithNodesDiscr(const MEDFileField1TSStructItem& other, const MEDFileMeshStruct *meshSt, const MEDFileFieldGlobsReal *globs) const
{
  if(other._type!=ON_NODES)
    throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem::isCompatibleWithNodesDiscr : other must be on nodes !");
  if(other._items.size()!=1)
    throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem::isCompatibleWithNodesDiscr : other is on nodes but number of subparts !");
  int theFirstLevFull;
  bool ret0=isFullyOnOneLev(meshSt,theFirstLevFull);
  const MEDFileField1TSStructItem2& otherNodeIt(other._items[0]);
  int nbOfNodes(meshSt->getNumberOfNodes());
  if(otherNodeIt.getPflName().empty())
    {//on all nodes
      if(!ret0)
        return false;
      std::vector<bool> nodesFetched(nbOfNodes,false);
      meshSt->getTheMesh()->whichAreNodesFetched(*this,globs,nodesFetched);
      if(std::find(nodesFetched.begin(),nodesFetched.end(),false)==nodesFetched.end())
        return theFirstLevFull==0;
      else
        return false;
    }
  else
    {
      const DataArrayInt *pfl=globs->getProfile(otherNodeIt.getPflName().c_str());
      MCAuto<DataArrayInt> cpyPfl(pfl->deepCopy());
      cpyPfl->sort();
      if(cpyPfl->isIota(nbOfNodes))
        {//on all nodes also !
          if(!ret0)
            return false;
          return theFirstLevFull==0;
        }
      std::vector<bool> nodesFetched(nbOfNodes,false);
      meshSt->getTheMesh()->whichAreNodesFetched(*this,globs,nodesFetched);
      return cpyPfl->isFittingWith(nodesFetched);
    }
}

bool MEDFileField1TSStructItem::isFullyOnOneLev(const MEDFileMeshStruct *meshSt, int& theFirstLevFull) const
{
  if(_type!=ON_CELLS)
    throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem::isFullyOnOneLev : works only for ON_CELLS discretization !");
  if(_items.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem::isFullyOnOneLev : items vector is empty !");
  int nbOfLevs(meshSt->getNumberOfLevs());
  if(nbOfLevs==0)
    throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem::isFullyOnOneLev : no levels in input mesh structure !");
  std::vector<int> levs(nbOfLevs);
  theFirstLevFull=1;
  std::set<INTERP_KERNEL::NormalizedCellType> gts;
  for(std::vector< MEDFileField1TSStructItem2 >::const_iterator it=_items.begin();it!=_items.end();it++)
    {
      if(!(*it).getPflName().empty())
        return false;
      INTERP_KERNEL::NormalizedCellType gt((*it).getGeo());
      if(gts.find(gt)!=gts.end())
        throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem::isFullyOnOneLev : internal error !");
      gts.insert(gt);
      int pos(meshSt->getLevelOfGeoType((*it).getGeo()));
      levs[-pos]++;
    }
  for(int i=0;i<nbOfLevs;i++)
    if(meshSt->getNumberOfGeoTypesInLev(-i)==levs[i])
      { theFirstLevFull=-i; return true; }
  return false;
}

const MEDFileField1TSStructItem2& MEDFileField1TSStructItem::operator[](std::size_t i) const
{
  if(i>=_items.size())
    throw INTERP_KERNEL::Exception("MEDFileField1TSStructItem::operator[] : input is not in valid range !");
  return _items[i];
}

std::size_t MEDFileField1TSStructItem::getHeapMemorySizeWithoutChildren() const
{
  std::size_t ret(_items.size()*sizeof(MEDFileField1TSStructItem2));
  return ret;
}

std::vector<const BigMemoryObject *> MEDFileField1TSStructItem::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  for(std::vector< MEDFileField1TSStructItem2 >::const_iterator it=_items.begin();it!=_items.end();it++)
    ret.push_back(&(*it));
  return ret;
}

MEDMeshMultiLev *MEDFileField1TSStructItem::buildFromScratchDataSetSupportOnCells(const MEDFileMeshStruct *mst, const MEDFileFieldGlobsReal *globs) const
{
  std::size_t sz(_items.size());
  std::vector<INTERP_KERNEL::NormalizedCellType> a0(sz);
  std::vector<const DataArrayInt *> a1(sz);
  std::vector<int> a2(sz);
  std::size_t i(0);
  for(std::vector< MEDFileField1TSStructItem2 >::const_iterator it=_items.begin();it!=_items.end();it++,i++)
    {
      a0[i]=(*it).getGeo();
      a1[i]=(*it).getPfl(globs);
      a2[i]=mst->getNumberOfElemsOfGeoType((*it).getGeo());
    }
  return MEDMeshMultiLev::New(mst->getTheMesh(),a0,a1,a2);
}

std::vector<INTERP_KERNEL::NormalizedCellType> MEDFileField1TSStructItem::getGeoTypes(const MEDFileMesh *m) const
{
  std::vector<INTERP_KERNEL::NormalizedCellType> ret;
  if(_type==ON_NODES)
    {
      if(!_items.empty() && _items[0].getPflName().empty())
        {
          if(m)
            return m->getAllGeoTypes();
          else
            return ret;
        }
      else
        return ret;
    }
  for(std::vector< MEDFileField1TSStructItem2 >::const_iterator it=_items.begin();it!=_items.end();it++)
    {
      INTERP_KERNEL::NormalizedCellType elt((*it).getGeo());
      std::vector<INTERP_KERNEL::NormalizedCellType>::iterator it2(std::find(ret.begin(),ret.end(),elt));
      if(it2==ret.end())
        ret.push_back(elt);
    }
  return ret;
}

MEDFileField1TSStructItem MEDFileField1TSStructItem::BuildItemFrom(const MEDFileAnyTypeField1TS *ref, const MEDFileMeshStruct *meshSt)
{
  std::vector< MEDFileField1TSStructItem2 > anItems;
  //
  std::vector< std::vector<std::string> > pfls,locs;
  std::vector< std::vector<TypeOfField> > typesF;
  std::vector<INTERP_KERNEL::NormalizedCellType> geoTypes;
  std::vector< std::vector<std::pair<int,int> > > strtEnds=ref->getFieldSplitedByType(std::string(),geoTypes,typesF,pfls,locs);
  std::size_t nbOfGeoTypes(geoTypes.size());
  if(nbOfGeoTypes==0)
    throw INTERP_KERNEL::Exception("MEDFileField1TSStruct : not null by empty ref  !");
  if(typesF[0].empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSStruct : internal error #1 bis !");
  TypeOfField atype(typesF[0][0]);
  for(std::size_t i=0;i<nbOfGeoTypes;i++)
    {
      std::size_t sz=typesF[i].size();
      if(strtEnds[i].size()<1 || sz<1 || pfls[i].size()<1)
        throw INTERP_KERNEL::Exception("MEDFileField1TSStruct : internal error #1 !");
      //
      for(std::size_t j=0;j<sz;j++)
        {
          if(atype==typesF[i][j])
            anItems.push_back(MEDFileField1TSStructItem2(geoTypes[i],strtEnds[i][j],pfls[i][j],locs[i][j]));
          else
            throw INTERP_KERNEL::Exception("MEDFileField1TSStruct : can be applied only on single spatial discretization fields ! Call SplitPerDiscretization method !");
        }
    }
  MEDFileField1TSStructItem ret(atype,anItems);
  try
    {
      ret.checkWithMeshStruct(meshSt,ref);
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      std::ostringstream oss; oss << e.what() << " (" << MEDCouplingFieldDiscretization::GetTypeOfFieldRepr(ret.getType()) << ")";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return ret;
}

//=

MEDFileField1TSStruct *MEDFileField1TSStruct::New(const MEDFileAnyTypeField1TS *ref, MEDFileMeshStruct *mst)
{
  return new MEDFileField1TSStruct(ref,mst);
}

MEDFileField1TSStruct::MEDFileField1TSStruct(const MEDFileAnyTypeField1TS *ref, MEDFileMeshStruct *mst)
{
  _already_checked.push_back(MEDFileField1TSStructItem::BuildItemFrom(ref,mst));
}

void MEDFileField1TSStruct::checkWithMeshStruct(MEDFileMeshStruct *mst, const MEDFileFieldGlobsReal *globs)
{
  if(_already_checked.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSStruct::checkWithMeshStruct : not correctly initialized !");
  _already_checked.back().checkWithMeshStruct(mst,globs);
}

bool MEDFileField1TSStruct::isEqualConsideringThePast(const MEDFileAnyTypeField1TS *other, const MEDFileMeshStruct *mst) const
{
  MEDFileField1TSStructItem b(MEDFileField1TSStructItem::BuildItemFrom(other,mst));
  for(std::vector<MEDFileField1TSStructItem>::const_iterator it=_already_checked.begin();it!=_already_checked.end();it++)
    {
      if((*it)==b)
        return true;
    }
  return false;
}

/*!
 * Not const because \a other structure will be added to the \c _already_checked attribute in case of success.
 */
bool MEDFileField1TSStruct::isSupportSameAs(const MEDFileAnyTypeField1TS *other, const MEDFileMeshStruct *meshSt)
{
  if(_already_checked.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSStruct::isSupportSameAs : no ref !");
  MEDFileField1TSStructItem b(MEDFileField1TSStructItem::BuildItemFrom(other,meshSt));
  if(!_already_checked[0].isEntityCell() || !b.isEntityCell())
    throw INTERP_KERNEL::Exception("MEDFileField1TSStruct::isSupportSameAs : only available on cell entities !");
  MEDFileField1TSStructItem other1(b.simplifyMeOnCellEntity(other));
  int found=-1,i=0;
  for(std::vector<MEDFileField1TSStructItem>::const_iterator it=_already_checked.begin();it!=_already_checked.end();it++,i++)
    if((*it).isComputed())
      { found=i; break; }
  bool ret(false);
  if(found==-1)
    {
      MEDFileField1TSStructItem this1(_already_checked[0].simplifyMeOnCellEntity(other));
      ret=this1.isCellSupportEqual(other1,other);
      if(ret)
        _already_checked.push_back(this1);
    }
  else
    ret=_already_checked[found].isCellSupportEqual(other1,other);
  if(ret)
    _already_checked.push_back(b);
  return ret;
}

/*!
 * \param [in] other - a field with only one spatial discretization : ON_NODES.
 */
bool MEDFileField1TSStruct::isCompatibleWithNodesDiscr(const MEDFileAnyTypeField1TS *other, const MEDFileMeshStruct *meshSt)
{
  if(_already_checked.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSStruct::isCompatibleWithNodesDiscr : no ref !");
  MEDFileField1TSStructItem other1(MEDFileField1TSStructItem::BuildItemFrom(other,meshSt));
  if(_already_checked[0].isEntityCell())
    {
      int found=-1,i=0;
      for(std::vector<MEDFileField1TSStructItem>::const_iterator it=_already_checked.begin();it!=_already_checked.end();it++,i++)
        if((*it).isComputed())
          { found=i; break; }
      bool ret(false);
      if(found==-1)
        {
          MEDFileField1TSStructItem this1(_already_checked[0].simplifyMeOnCellEntity(other));
          ret=this1.isCompatibleWithNodesDiscr(other1,meshSt,other);
          if(ret)
            _already_checked.push_back(this1);
        }
      else
        ret=_already_checked[found].isCompatibleWithNodesDiscr(other1,meshSt,other);
      if(ret)
        _already_checked.push_back(other1);
      return ret;
    }
  else
    return _already_checked[0].isNodeSupportEqual(other1,other);
}

std::size_t MEDFileField1TSStruct::getHeapMemorySizeWithoutChildren() const
{
  std::size_t ret(_already_checked.capacity()*sizeof(MEDFileField1TSStructItem));
  return ret;
}

std::vector<const BigMemoryObject *> MEDFileField1TSStruct::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  for(std::vector<MEDFileField1TSStructItem>::const_iterator it=_already_checked.begin();it!=_already_checked.end();it++)
    ret.push_back(&(*it));
  return ret;
}

MEDMeshMultiLev *MEDFileField1TSStruct::buildFromScratchDataSetSupport(const MEDFileMeshStruct *mst, const MEDFileFieldGlobsReal *globs) const
{
  if(_already_checked.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSStruct::buildFromScratchDataSetSupport : No outline structure in this !");
  int pos0(-1),pos1(-1);
  if(presenceOfCellDiscr(pos0))
    {
      MCAuto<MEDMeshMultiLev> ret(_already_checked[pos0].buildFromScratchDataSetSupportOnCells(mst,globs));
      if(presenceOfPartialNodeDiscr(pos1))
        ret->setNodeReduction(_already_checked[pos1][0].getPfl(globs));
      return ret.retn();
    }
  else
    {
      if(!presenceOfPartialNodeDiscr(pos1))
        {//we have only all nodes, no cell definition info -> all existing levels !;
          return MEDMeshMultiLev::New(mst->getTheMesh(),mst->getTheMesh()->getNonEmptyLevels());
        }
      else
        return MEDMeshMultiLev::NewOnlyOnNode(mst->getTheMesh(),_already_checked[pos1][0].getPfl(globs));
    }
}

bool MEDFileField1TSStruct::isDataSetSupportFastlyEqualTo(const MEDFileField1TSStruct& other, const MEDFileFieldGlobsReal *globs) const
{
  int b0,b1;
  bool a0(presenceOfCellDiscr(b0)),a1(presenceOfPartialNodeDiscr(b1));
  int d0,d1;
  bool c0(other.presenceOfCellDiscr(d0)),c1(other.presenceOfPartialNodeDiscr(d1)); 
  if(a0!=c0 || a1!=c1)
    return false;
  if(a0)
    if(!_already_checked[b0].isCellSupportEqual(other._already_checked[d0],globs))
      return false;
  if(a1)
    if(!_already_checked[b1].isNodeSupportEqual(other._already_checked[d1],globs))
      return false;
  return true;
}

std::vector<INTERP_KERNEL::NormalizedCellType> MEDFileField1TSStruct::getGeoTypes(const MEDFileMesh *m) const
{
  std::vector<INTERP_KERNEL::NormalizedCellType> ret;
  for(std::vector<MEDFileField1TSStructItem>::const_iterator it=_already_checked.begin();it!=_already_checked.end();it++)
    {
      std::vector<INTERP_KERNEL::NormalizedCellType> ret2((*it).getGeoTypes(m));
      for(std::vector<INTERP_KERNEL::NormalizedCellType>::const_iterator it2=ret2.begin();it2!=ret2.end();it2++)
        {
          if(*it2==INTERP_KERNEL::NORM_ERROR)
            continue;
          std::vector<INTERP_KERNEL::NormalizedCellType>::iterator it3(std::find(ret.begin(),ret.end(),*it2));
          if(it3==ret.end())
            ret.push_back(*it2);
        }
    }
  return ret;
}

/*!
 * Returns true if presence in \a this of discretization ON_CELLS, ON_GAUSS_PT, ON_GAUSS_NE.
 * If true is returned the pos of the easiest is returned. The easiest is the first element in \a this having the less splitted subparts.
 */
bool MEDFileField1TSStruct::presenceOfCellDiscr(int& pos) const
{
  std::size_t refSz(std::numeric_limits<std::size_t>::max());
  bool ret(false);
  int i(0);
  for(std::vector<MEDFileField1TSStructItem>::const_iterator it=_already_checked.begin();it!=_already_checked.end();it++,i++)
    {
      if((*it).getType()!=ON_NODES)
        {
          ret=true;
          std::size_t sz((*it).getNumberOfItems());
          if(refSz>sz)
            { pos=i; refSz=sz; }
        }
    }
  if(refSz==0)
    throw INTERP_KERNEL::Exception("MEDFileField1TSStruct::presenceOfCellDiscr : an element in this on entity CELL is empty !");
  return ret;
}

/*!
 * Returns true if presence in \a this of discretization ON_NODES.
 * If true is returned the pos of the first element containing the single subpart.
 */
bool MEDFileField1TSStruct::presenceOfPartialNodeDiscr(int& pos) const
{
  int i(0);
  for(std::vector<MEDFileField1TSStructItem>::const_iterator it=_already_checked.begin();it!=_already_checked.end();it++,i++)
    {
      if((*it).getType()==ON_NODES)
        {
          std::size_t sz((*it).getNumberOfItems());
          if(sz==1)
            {
              if(!(*it)[0].getPflName().empty())
                { pos=i; return true; }
            }
          else
            throw INTERP_KERNEL::Exception("MEDFileField1TSStruct::presenceOfPartialNodeDiscr : an element in this on entity NODE is split into several parts !");
        }
    }
  return false;
}

//=

MEDFileFastCellSupportComparator *MEDFileFastCellSupportComparator::New(const MEDFileMeshStruct *m, const MEDFileAnyTypeFieldMultiTS *ref)
{
  return new MEDFileFastCellSupportComparator(m,ref);
}

MEDFileFastCellSupportComparator::MEDFileFastCellSupportComparator(const MEDFileMeshStruct *m, const MEDFileAnyTypeFieldMultiTS *ref)
{
  if(!m)
    throw INTERP_KERNEL::Exception("MEDFileFastCellSupportComparator constructor : null input mesh struct !");
  _mesh_comp=const_cast<MEDFileMeshStruct *>(m); _mesh_comp->incrRef();
  int nbPts=ref->getNumberOfTS();
  _f1ts_cmps.resize(nbPts);
  for(int i=0;i<nbPts;i++)
    {
      MCAuto<MEDFileAnyTypeField1TS> elt=ref->getTimeStepAtPos(i);
      try
        {
          _f1ts_cmps[i]=MEDFileField1TSStruct::New(elt,_mesh_comp);
          _f1ts_cmps[i]->checkWithMeshStruct(_mesh_comp,elt);
        }
      catch(INTERP_KERNEL::Exception& e)
        {
          std::ostringstream oss; oss << "Problem in field with name \"" << ref->getName() << "\"" << std::endl;
          oss << "More Details : " << e.what();
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
}

std::size_t MEDFileFastCellSupportComparator::getHeapMemorySizeWithoutChildren() const
{
  std::size_t ret(_f1ts_cmps.capacity()*sizeof(MCAuto<MEDFileField1TSStruct>));
  return ret;
}

std::vector<const BigMemoryObject *> MEDFileFastCellSupportComparator::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  const MEDFileMeshStruct *mst(_mesh_comp);
  if(mst)
    ret.push_back(mst);
  for(std::vector< MCAuto<MEDFileField1TSStruct> >::const_iterator it=_f1ts_cmps.begin();it!=_f1ts_cmps.end();it++)
    ret.push_back((const MEDFileField1TSStruct *)*it);
  return ret;
}

bool MEDFileFastCellSupportComparator::isEqual(const MEDFileAnyTypeFieldMultiTS *other)
{
  int nbPts=other->getNumberOfTS();
  if(nbPts!=(int)_f1ts_cmps.size())
    {
      std::ostringstream oss; oss << "MEDFileFastCellSupportComparator::isEqual : unexpected nb of time steps in  input ! Should be " << _f1ts_cmps.size() << " it is in reality " << nbPts << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  for(int i=0;i<nbPts;i++)
    {
      MCAuto<MEDFileAnyTypeField1TS> elt=other->getTimeStepAtPos(i);
      if(!_f1ts_cmps[i]->isEqualConsideringThePast(elt,_mesh_comp))
        if(!_f1ts_cmps[i]->isSupportSameAs(elt,_mesh_comp))
          return false;
    }
  return true;
}

bool MEDFileFastCellSupportComparator::isCompatibleWithNodesDiscr(const MEDFileAnyTypeFieldMultiTS *other)
{
  int nbPts=other->getNumberOfTS();
  if(nbPts!=(int)_f1ts_cmps.size())
    {
      std::ostringstream oss; oss << "MEDFileFastCellSupportComparator::isCompatibleWithNodesDiscr : unexpected nb of time steps in  input ! Should be " << _f1ts_cmps.size() << " it is in reality " << nbPts << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  for(int i=0;i<nbPts;i++)
    {
      MCAuto<MEDFileAnyTypeField1TS> elt=other->getTimeStepAtPos(i);
      if(!_f1ts_cmps[i]->isCompatibleWithNodesDiscr(elt,_mesh_comp))
        return false;
    }
  return true;
}

MEDMeshMultiLev *MEDFileFastCellSupportComparator::buildFromScratchDataSetSupport(int timeStepId, const MEDFileFieldGlobsReal *globs) const
{
  if(timeStepId<0 || timeStepId>=(int)_f1ts_cmps.size())
    {
      std::ostringstream oss; oss << "MEDFileFastCellSupportComparator::buildFromScratchDataSetSupport : requested time step id #" << timeStepId << " is not in [0," << _f1ts_cmps.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  const MEDFileField1TSStruct *obj(_f1ts_cmps[timeStepId]);
  if(!obj)
    {
      std::ostringstream oss; oss << "MEDFileFastCellSupportComparator::buildFromScratchDataSetSupport : at time step id #" << timeStepId << " no field structure overview defined !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return obj->buildFromScratchDataSetSupport(_mesh_comp,globs);
}

bool MEDFileFastCellSupportComparator::isDataSetSupportEqualToThePreviousOne(int timeStepId, const MEDFileFieldGlobsReal *globs) const
{
  if(timeStepId<=0 || timeStepId>=(int)_f1ts_cmps.size())
    {
      std::ostringstream oss; oss << "MEDFileFastCellSupportComparator::isDataSetSupportEqualToThePreviousOne : requested time step id #" << timeStepId << " is not in [1," << _f1ts_cmps.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  const MEDFileField1TSStruct *obj(_f1ts_cmps[timeStepId]);
  const MEDFileField1TSStruct *objRef(_f1ts_cmps[timeStepId-1]);
  return objRef->isDataSetSupportFastlyEqualTo(*obj,globs);
}

int MEDFileFastCellSupportComparator::getNumberOfTS() const
{
  return _f1ts_cmps.size();
}

std::vector<INTERP_KERNEL::NormalizedCellType> MEDFileFastCellSupportComparator::getGeoTypesAt(int timeStepId, const MEDFileMesh *m) const
{
  if(timeStepId<0 || timeStepId>=(int)_f1ts_cmps.size())
    {
      std::ostringstream oss; oss << "MEDFileFastCellSupportComparator::getGeoTypesAt : requested time step id #" << timeStepId << " is not in [0," << _f1ts_cmps.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  const MEDFileField1TSStruct *elt(_f1ts_cmps[timeStepId]);
  if(!elt)
    {
      std::ostringstream oss; oss << "MEDFileFastCellSupportComparator::getGeoTypesAt : requested time step id #" << timeStepId << " points to a NULL pointer !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return elt->getGeoTypes(m);
}
