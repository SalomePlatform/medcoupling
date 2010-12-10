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

#include "MEDFileMeshLL.hxx"
#include "MEDLoaderBase.hxx"

#include "MEDCouplingUMesh.hxx"

#include "InterpKernelAutoPtr.hxx"

#include <set>

extern med_geometrie_element typmai[MED_NBR_GEOMETRIE_MAILLE+2];
extern INTERP_KERNEL::NormalizedCellType typmai2[MED_NBR_GEOMETRIE_MAILLE+2];
extern med_geometrie_element typmainoeud[1];

using namespace ParaMEDMEM;

MEDFileMeshL2::MEDFileMeshL2():_name(MED_TAILLE_NOM),_description(MED_TAILLE_DESC)
{
}

MEDFileUMeshL2::MEDFileUMeshL2()
{
}

void MEDFileUMeshL2::loadAll(med_idt fid, int mId, const char *mName)
{
  _name.set(mName);
  med_maillage type_maillage;
  med_int Mdim;
  if(MEDmaaInfo(fid,mId,(char *)mName,&Mdim,&type_maillage,_description.getPointer())!=0)
    throw INTERP_KERNEL::Exception("A problem has been detected when trying to get info on mesh !");
  if(type_maillage!=MED_NON_STRUCTURE)
    throw INTERP_KERNEL::Exception("Invalid mesh type ! You are expected an unstructured one whereas in file it is not an unstructured !");
  loadConnectivity(fid,Mdim,mName);
  loadCoords(fid,mId,Mdim,mName);
}

void MEDFileUMeshL2::loadConnectivity(med_idt fid, int mdim, const char *mName)
{
  _per_type_mesh.resize(1);
  _per_type_mesh[0].clear();
  for(int j=0;j<MED_NBR_GEOMETRIE_MAILLE+2;j++)
    {
      MEDFileUMeshPerType *tmp=MEDFileUMeshPerType::New(fid,mName,mdim,typmai[j],typmai2[j]);
      if(tmp)
        _per_type_mesh[0].push_back(tmp);
    }
  sortTypes();
}

void MEDFileUMeshL2::loadCoords(med_idt fid, int mId, int mdim, const char *mName) throw(INTERP_KERNEL::Exception)
{
  med_int edim=MEDdimEspaceLire(fid,(char *)mName);
  int spaceDim=std::max((int)mdim,(int)edim);
  int nCoords=MEDnEntMaa(fid,(char *)mName,MED_COOR,MED_NOEUD,(med_geometrie_element)0,(med_connectivite)0);
  _coords=DataArrayDouble::New();
  _coords->alloc(nCoords,spaceDim);
  double *coordsPtr=_coords->getPointer();
  med_repere repere;
  char *comp=MEDLoaderBase::buildEmptyString(spaceDim*MED_TAILLE_PNOM);
  char *unit=MEDLoaderBase::buildEmptyString(spaceDim*MED_TAILLE_PNOM);
  MEDcoordLire(fid,(char *)mName,spaceDim,coordsPtr,MED_FULL_INTERLACE,MED_ALL,NULL,0,&repere,comp,unit);
  _fam_coords=DataArrayInt::New();
  _fam_coords->alloc(nCoords,1);
  _num_coords=DataArrayInt::New();
  _num_coords->alloc(nCoords,1);
  MEDfamLire(fid,(char *)mName,_fam_coords->getPointer(),nCoords,MED_NOEUD,MED_NONE);
  if(MEDnumLire(fid,(char *)mName,_num_coords->getPointer(),nCoords,MED_NOEUD,MED_NONE)!=0)
    _num_coords=0;
  for(int i=0;i<spaceDim;i++)
    {
      std::string n,u;
      std::string info=MEDLoaderBase::buildUnionUnit(comp+i*MED_TAILLE_PNOM,MED_TAILLE_PNOM,unit+i*MED_TAILLE_PNOM,MED_TAILLE_PNOM);
      _coords->setInfoOnComponent(i,info.c_str());
    }
  delete [] comp;
  delete [] unit;
}

void MEDFileUMeshL2::sortTypes()
{
  std::set<int> mdims;
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshPerType> > tmp(_per_type_mesh[0]);
  _per_type_mesh.clear();
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshPerType> >::const_iterator it=tmp.begin();it!=tmp.end();it++)
    mdims.insert((*it)->getDim());
  if(mdims.empty())
    return;
  int mdim=*mdims.rbegin();
  _per_type_mesh.resize(mdim+1);
  for(int dim=mdim+1;dim!=0;dim--)
    {
      std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshPerType> >& elt=_per_type_mesh[mdim+1-dim];
      for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshPerType> >::const_iterator it=tmp.begin();it!=tmp.end();it++)
        if((*it)->getDim()==dim-1)
          elt.push_back(*it);
    }
}

int MEDFileUMeshL2::getMeshIdFromName(med_idt fid, const char *mname) throw(INTERP_KERNEL::Exception)
{
  med_maillage type_maillage;
  char maillage_description[MED_TAILLE_DESC+1];
  med_int dim;
  char nommaa[MED_TAILLE_NOM+1];
  med_int n=MEDnMaa(fid);
  bool found=false;
  int ret=-1;
  std::vector<std::string> ms;
  for(int i=0;i<n;i++)
    {
      MEDmaaInfo(fid,i+1,nommaa,&dim,&type_maillage,maillage_description);
      std::string cur=MEDLoaderBase::buildStringFromFortran(nommaa,sizeof(nommaa));
      ms.push_back(cur);
      if(cur==mname)
        {
          found=true;
          ret=i+1;
        }
    }
  if(!found)
    {
      std::ostringstream oss;
      oss << "No such meshname (" << mname <<  ") in file ! Must be in :";
      std::copy(ms.begin(),ms.end(),std::ostream_iterator<std::string>(oss,", "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return ret;
}

void MEDFileUMeshL2::readFamiliesAndGrps(med_idt fid, const char *meshName, std::map<std::string,int>& fams, std::map<std::string, std::vector<std::string> >& grps)
{
  char nomfam[MED_TAILLE_NOM+1];
  med_int numfam;
  int nfam=MEDnFam(fid,(char *)meshName);
  for(int i=0;i<nfam;i++)
    {
      int ngro=MEDnGroupe(fid,(char *)meshName,i+1);
      med_int natt=MEDnAttribut(fid,(char *)meshName,i+1);
      med_int *attide=new int[natt];
      med_int *attval=new int[natt];
      char *attdes=new char[MED_TAILLE_DESC*natt+1];
      char *gro=new char[MED_TAILLE_LNOM*ngro+1];
      MEDfamInfo(fid,(char *)meshName,i+1,nomfam,&numfam,attide,attval,attdes,&natt,gro,&ngro);
      std::string famName=MEDLoaderBase::buildStringFromFortran(nomfam,MED_TAILLE_LNOM);
      fams[famName]=numfam;
      for(int j=0;j<ngro;j++)
        {
          std::string groupname=MEDLoaderBase::buildStringFromFortran(gro+j*MED_TAILLE_LNOM,MED_TAILLE_LNOM);
          grps[groupname].push_back(famName);
        }
      delete [] attdes;
      delete [] gro;
      delete [] attide;
      delete [] attval;
    }
}

void MEDFileUMeshL2::writeFamiliesAndGrps(med_idt fid, const char *mname, const std::map<std::string,int>& fams, const std::map<std::string, std::vector<std::string> >& grps)
{
  for(std::map<std::string,int>::const_iterator it=fams.begin();it!=fams.end();it++)
    {
      std::vector<std::string> grpsOfFam;
      for(std::map<std::string, std::vector<std::string> >::const_iterator it1=grps.begin();it1!=grps.end();it1++)
        {
          if(std::find((*it1).second.begin(),(*it1).second.end(),(*it).first)!=(*it1).second.end())
            grpsOfFam.push_back((*it1).first);
        }
      int ngro=grpsOfFam.size();
      INTERP_KERNEL::AutoPtr<char> groName=MEDLoaderBase::buildEmptyString(MED_TAILLE_LNOM*ngro);
      int i=0;
      for(std::vector<std::string>::const_iterator it2=grpsOfFam.begin();it2!=grpsOfFam.end();it2++,i++)
        MEDLoaderBase::safeStrCpy((*it2).c_str(),MED_TAILLE_LNOM-1,groName+i*MED_TAILLE_LNOM,0);//tony too long
      INTERP_KERNEL::AutoPtr<char> famName=MEDLoaderBase::buildEmptyString(MED_TAILLE_NOM);
      MEDLoaderBase::safeStrCpy((*it).first.c_str(),MED_TAILLE_NOM,famName,0);//tony too long
      MEDfamCr(fid,(char *)mname,famName,(*it).second,0,0,0,0,groName,ngro);
    }
}

void MEDFileUMeshL2::writeCoords(med_idt fid, const char *mname, const DataArrayDouble *coords, const DataArrayInt *famCoords, const DataArrayInt *numCoords)
{
  if(!coords)
    return ;
  int spaceDim=coords->getNumberOfComponents();
  INTERP_KERNEL::AutoPtr<char> comp=MEDLoaderBase::buildEmptyString(spaceDim*MED_TAILLE_PNOM);
  INTERP_KERNEL::AutoPtr<char> unit=MEDLoaderBase::buildEmptyString(spaceDim*MED_TAILLE_PNOM);
  for(int i=0;i<spaceDim;i++)
    {
      std::string info=coords->getInfoOnComponent(i);
      std::string c,u;
      MEDLoaderBase::splitIntoNameAndUnit(info,c,u);
      MEDLoaderBase::safeStrCpy(c.c_str(),MED_TAILLE_PNOM-1,comp+i*MED_TAILLE_PNOM,0);//MED_TAILLE_PNOM-1 to avoid to write '\0' on next compo
      MEDLoaderBase::safeStrCpy(u.c_str(),MED_TAILLE_PNOM-1,unit+i*MED_TAILLE_PNOM,0);//MED_TAILLE_PNOM-1 to avoid to write '\0' on next compo
    }
  MEDcoordEcr(fid,(char *)mname,spaceDim,coords->getPointer(),MED_FULL_INTERLACE,coords->getNumberOfTuples(),MED_CART,comp,unit);
  MEDfamEcr(fid,(char *)mname,famCoords->getPointer(),famCoords->getNumberOfTuples(),MED_NOEUD,MED_NONE);
  if(numCoords)
    MEDnumEcr(fid,(char *)mname,numCoords->getPointer(),numCoords->getNumberOfTuples(),MED_NOEUD,MED_NONE);
}

bool MEDFileUMeshL2::isFamDefinedOnLev(int levId) const
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshPerType> >::const_iterator it=_per_type_mesh[levId].begin();it!=_per_type_mesh[levId].end();it++)
    if((*it)->getFam()==0)
      return false;
  return true;
}

bool MEDFileUMeshL2::isNumDefinedOnLev(int levId) const
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshPerType> >::const_iterator it=_per_type_mesh[levId].begin();it!=_per_type_mesh[levId].end();it++)
    if((*it)->getNum()==0)
      return false;
  return true;
}

MEDFileUMeshSplitL1::MEDFileUMeshSplitL1(const MEDFileUMeshL2& l2, const char *mName, int id)
{
  const std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshPerType> >& v=l2.getLev(id);
  if(v.empty())
    return;
  int sz=v.size();
  std::vector<MEDCouplingUMesh *> ms(sz);
  for(int i=0;i<sz;i++)
    {
      MEDCouplingUMesh *tmp=MEDCouplingUMesh::New("",v[i]->getDim());
      MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> tmp2=l2.getCoords();
      tmp->setCoords(tmp2);
      tmp->setConnectivity(const_cast<DataArrayInt *>(v[i]->getNodal()),const_cast<DataArrayInt *>(v[i]->getNodalIndex()));
      ms[i]=tmp;
    }
  _m_by_types=MEDCouplingUMesh::mergeUMeshesOnSameCoords(ms);
  _m_by_types->setName(mName);
  if(l2.isFamDefinedOnLev(id))
    {
      int nbOfCells=_m_by_types->getNumberOfCells();
      _fam=DataArrayInt::New();
      _fam->alloc(nbOfCells,1);
      int *w=_fam->getPointer();
      for(int i=0;i<sz;i++)
        w=std::copy(v[i]->getFam()->getConstPointer(),v[i]->getFam()->getConstPointer()+v[i]->getFam()->getNumberOfTuples(),w);
    }
  if(l2.isNumDefinedOnLev(id))
    {
      int nbOfCells=_m_by_types->getNumberOfCells();
      _num=DataArrayInt::New();
      _num->alloc(nbOfCells,1);
      int *w=_num->getPointer();
      for(int i=0;i<sz;i++)
        w=std::copy(v[i]->getNum()->getConstPointer(),v[i]->getNum()->getConstPointer()+v[i]->getNum()->getNumberOfTuples(),w);
      _m=(MEDCouplingUMesh *)_m_by_types->deepCpy();
      _m->renumberCells(_num->getConstPointer(),true);
    }
  else
    _m=_m_by_types;
  for(int i=0;i<sz;i++)
    ms[i]->decrRef();
}

MEDFileUMeshSplitL1::MEDFileUMeshSplitL1(MEDCouplingUMesh *m)
{
  m->incrRef();
  _m=m;
  _m_by_types=(MEDCouplingUMesh *)_m->deepCpy();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da=_m_by_types->getRenumArrForConsecutiveCellTypesSpec(typmai2,typmai2+MED_NBR_GEOMETRIE_MAILLE+2);
  _num=da->invertArrayO2N2N2O(m->getNumberOfCells());
  _fam=DataArrayInt::New();
  _fam->alloc(m->getNumberOfCells(),1);
  _fam->fillWithValue(0);
  _m_by_types->renumberCells(da->getConstPointer(),false);
}

bool MEDFileUMeshSplitL1::empty() const
{
  return ((const MEDCouplingUMesh *)_m_by_types)==0;
}

int MEDFileUMeshSplitL1::getMeshDimension() const
{
  return _m_by_types->getMeshDimension();
}

MEDCouplingUMesh *MEDFileUMeshSplitL1::getFamilyPart(const std::vector<int>& ids, bool renum) const
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> eltsToKeep=_fam->getIdsEqualList(ids);
  MEDCouplingUMesh *m=(MEDCouplingUMesh *)_m_by_types->buildPartOfMySelf(eltsToKeep->getConstPointer(),eltsToKeep->getConstPointer()+eltsToKeep->getNumberOfTuples(),true);
  if(renum)
    return renumIfNeeded(m,eltsToKeep->getConstPointer());
  return m;
}

DataArrayInt *MEDFileUMeshSplitL1::getFamilyPartArr(const std::vector<int>& ids, bool renum) const
{
  DataArrayInt *da=_fam->getIdsEqualList(ids);
  if(renum)
    return renumIfNeededArr(da);
  return da;
}

MEDCouplingUMesh *MEDFileUMeshSplitL1::getWholeMesh(bool renum) const
{
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> tmp;
  if(renum)
    tmp=_m;
  else
    tmp=_m_by_types;
  tmp->incrRef();
  return tmp;
}

/*!
 * This method ignores _m and _m_by_types.
 */
void MEDFileUMeshSplitL1::setGroupsFromScratch(const std::vector<MEDCouplingUMesh *>& ms, std::map<std::string,int>& familyIds,
                                               std::map<std::string, std::vector<std::string> >& groups) throw(INTERP_KERNEL::Exception)
{
  int sz=ms.size();
  std::vector< DataArrayInt * > corr;
  _m=MEDCouplingUMesh::fuseUMeshesOnSameCoords(ms,0,corr);
  std::vector< std::vector<int> > fidsOfGroups;
  _fam=DataArrayInt::makePartition(corr,_m->getNumberOfCells(),fidsOfGroups);
  int nbOfCells=_m->getNumberOfCells();
  std::map<int,std::string> newfams;
  std::map<int,int> famIdTrad;
  traduceFamilyNumber(fidsOfGroups,familyIds,famIdTrad,newfams);
  for(int i=0;i<sz;i++)
    corr[i]->decrRef();
  int *w=_fam->getPointer();
  for(int i=0;i<nbOfCells;i++,w++)
    *w=famIdTrad[*w];
}

void MEDFileUMeshSplitL1::write(med_idt fid, const char *mName, int mdim) const
{
  std::vector<MEDCouplingUMesh *> ms=_m_by_types->splitByType();
  int start=0;
  for(std::vector<MEDCouplingUMesh *>::const_iterator it=ms.begin();it!=ms.end();it++)
    {
      int nbCells=(*it)->getNumberOfCells();
      int end=start+nbCells;
      DataArrayInt *fam=0,*num=0;
      if((const DataArrayInt *)_fam)
        fam=_fam->substr(start,end);
      if((const DataArrayInt *)_num)
        num=_num->substr(start,end);
      MEDFileUMeshPerType::write(fid,mName,mdim,(*it),fam,num);
      if(fam)
        fam->decrRef();
      if(num)
        num->decrRef();
      (*it)->decrRef();
      start=end;
    }
}

MEDCouplingUMesh *MEDFileUMeshSplitL1::renumIfNeeded(MEDCouplingUMesh *m, const int *cellIds) const
{
  if((const DataArrayInt *)_num==0)
    return m;
  if(cellIds==0)
    m->renumberCells(_num->getConstPointer(),true);
  else
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> locnum=_num->selectByTupleId(cellIds,cellIds+m->getNumberOfCells());
      m->renumberCells(locnum->getConstPointer(),true);
    }
  return m;
}

DataArrayInt *MEDFileUMeshSplitL1::renumber(const DataArrayInt *renum, DataArrayInt *da)
{
  if((const DataArrayInt *)renum==0)
    return da;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> locnum=renum->selectByTupleId(da->getConstPointer(),da->getConstPointer()+da->getNumberOfTuples());
  da->decrRef();
  locnum->incrRef();
  return locnum;
}

DataArrayInt *MEDFileUMeshSplitL1::renumIfNeededArr(DataArrayInt *da) const
{
  return renumber(_num,da);
}

std::vector<int> MEDFileUMeshSplitL1::getNewFamiliesNumber(int nb, const std::map<std::string,int>& families)
{
  int id=-1;
  for(std::map<std::string,int>::const_iterator it=families.begin();it!=families.end();it++)
    id=std::max(id,(*it).second);
  if(id==-1)
    id=0;
  std::vector<int> ret(nb);
  for(int i=1;i<=nb;i++)
    ret[i]=id+i;
  return ret;
}

void MEDFileUMeshSplitL1::traduceFamilyNumber(const std::vector< std::vector<int> >& fidsGrps, std::map<std::string,int>& familyIds,
                                              std::map<int,int>& famIdTrad, std::map<int,std::string>& newfams)
{
  std::set<int> allfids;
  
}
