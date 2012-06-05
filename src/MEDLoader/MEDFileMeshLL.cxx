// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
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

#include "MEDFileMeshLL.hxx"
#include "MEDFileMesh.hxx"
#include "MEDLoaderBase.hxx"

#include "MEDCouplingUMesh.hxx"

#include "InterpKernelAutoPtr.hxx"
#include "CellModel.hxx"

#include <set>

extern med_geometry_type typmai[MED_N_CELL_FIXED_GEO];
extern INTERP_KERNEL::NormalizedCellType typmai2[MED_N_CELL_FIXED_GEO];
extern med_geometry_type typmainoeud[1];

using namespace ParaMEDMEM;

MEDFileMeshL2::MEDFileMeshL2():_name(MED_NAME_SIZE),_description(MED_COMMENT_SIZE),_dt_unit(MED_LNAME_SIZE)
{
}

int MEDFileMeshL2::GetMeshIdFromName(med_idt fid, const char *mname, ParaMEDMEM::MEDCouplingMeshType& meshType, int& dt, int& it, std::string& dtunit1) throw(INTERP_KERNEL::Exception)
{
  med_mesh_type type_maillage;
  char maillage_description[MED_COMMENT_SIZE+1];
  char dtunit[MED_LNAME_SIZE+1];
  med_int spaceDim,dim;
  char nommaa[MED_NAME_SIZE+1];
  med_int n=MEDnMesh(fid);
  bool found=false;
  int ret=-1;
  med_sorting_type stype;
  std::vector<std::string> ms;
  int nstep;
  med_axis_type axistype;
  for(int i=0;i<n && !found;i++)
    {
      int naxis=MEDmeshnAxis(fid,i+1);
      INTERP_KERNEL::AutoPtr<char> axisname=MEDLoaderBase::buildEmptyString(naxis*MED_SNAME_SIZE);
      INTERP_KERNEL::AutoPtr<char> axisunit=MEDLoaderBase::buildEmptyString(naxis*MED_SNAME_SIZE);
      MEDmeshInfo(fid,i+1,nommaa,&spaceDim,&dim,&type_maillage,maillage_description,dtunit,&stype,&nstep,&axistype,axisname,axisunit);
      dtunit1=MEDLoaderBase::buildStringFromFortran(dtunit,sizeof(dtunit));
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
  switch(type_maillage)
    {
    case MED_UNSTRUCTURED_MESH:
      meshType=UNSTRUCTURED;
      break;
    case MED_STRUCTURED_MESH:
      meshType=CARTESIAN;
      break;
    default:
      throw INTERP_KERNEL::Exception("MEDFileUMeshL2::getMeshIdFromName : unrecognized mesh type !");
    }
  med_int numdt,numit;
  med_float dtt;
  MEDmeshComputationStepInfo(fid,mname,1,&numdt,&numit,&dtt);
  dt=numdt; it=numit;
  return ret;
}

double MEDFileMeshL2::CheckMeshTimeStep(med_idt fid, const char *mName, int nstep, int dt, int it) throw(INTERP_KERNEL::Exception)
{
  bool found=false;
  med_int numdt,numit;
  med_float dtt;
  std::vector< std::pair<int,int> > p(nstep);
  for(int i=0;i<nstep;i++)
    {
      MEDmeshComputationStepInfo(fid,mName,i+1,&numdt,&numit,&dtt);
      p[i]=std::make_pair<int,int>(numdt,numit);
      found=(numdt==dt) && (numit==numit);
    }
  if(!found)
    {
      std::ostringstream oss; oss << "No such iteration=" << dt << ",order=" << it << " numbers found for mesh '" << mName << "' ! ";
      oss << "Possibilities are : ";
      for(int i=0;i<nstep;i++)
        oss << "(" << p[i].first << "," << p[i].second << "), ";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return dtt;
}

std::vector<std::string> MEDFileMeshL2::getAxisInfoOnMesh(med_idt fid, int mId, const char *mName, ParaMEDMEM::MEDCouplingMeshType& meshType, int& nstep, int& Mdim) throw(INTERP_KERNEL::Exception)
{
  med_mesh_type type_maillage;
  med_int spaceDim;
  med_sorting_type stype;
  med_axis_type axistype;
  int naxis=MEDmeshnAxis(fid,mId);
  INTERP_KERNEL::AutoPtr<char> nameTmp=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> axisname=MEDLoaderBase::buildEmptyString(naxis*MED_SNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> axisunit=MEDLoaderBase::buildEmptyString(naxis*MED_SNAME_SIZE);
  if(MEDmeshInfo(fid,mId,nameTmp,&spaceDim,&Mdim,&type_maillage,_description.getPointer(),_dt_unit.getPointer(),
                 &stype,&nstep,&axistype,axisname,axisunit)!=0)
    throw INTERP_KERNEL::Exception("A problem has been detected when trying to get info on mesh !");
  switch(type_maillage)
    {
    case MED_UNSTRUCTURED_MESH:
      meshType=UNSTRUCTURED;
      break;
    case MED_STRUCTURED_MESH:
      meshType=CARTESIAN;
      break;
    default:
      throw INTERP_KERNEL::Exception("MEDFileUMeshL2::getMeshIdFromName : unrecognized mesh type !");
    }
  //
  std::vector<std::string> infosOnComp(naxis);
  for(int i=0;i<naxis;i++)
    {
      std::string info=MEDLoaderBase::buildUnionUnit(((char *)axisname)+i*MED_SNAME_SIZE,MED_SNAME_SIZE,((char *)axisunit)+i*MED_SNAME_SIZE,MED_SNAME_SIZE);
      infosOnComp[i]=info;
    }
  return infosOnComp;
}

void MEDFileMeshL2::ReadFamiliesAndGrps(med_idt fid, const char *meshName, std::map<std::string,int>& fams, std::map<std::string, std::vector<std::string> >& grps)
{
  char nomfam[MED_NAME_SIZE+1];
  med_int numfam;
  int nfam=MEDnFamily(fid,meshName);
  for(int i=0;i<nfam;i++)
    {
      int ngro=MEDnFamilyGroup(fid,meshName,i+1);
      med_int natt=MEDnFamily23Attribute(fid,meshName,i+1);
      INTERP_KERNEL::AutoPtr<med_int> attide=new med_int[natt];
      INTERP_KERNEL::AutoPtr<med_int> attval=new med_int[natt];
      INTERP_KERNEL::AutoPtr<char> attdes=new char[MED_COMMENT_SIZE*natt+1];
      INTERP_KERNEL::AutoPtr<char> gro=new char[MED_LNAME_SIZE*ngro+1];
      MEDfamily23Info(fid,meshName,i+1,nomfam,attide,attval,attdes,&numfam,gro);
      std::string famName=MEDLoaderBase::buildStringFromFortran(nomfam,MED_NAME_SIZE);
      fams[famName]=numfam;
      for(int j=0;j<ngro;j++)
        {
          std::string groupname=MEDLoaderBase::buildStringFromFortran(gro+j*MED_LNAME_SIZE,MED_LNAME_SIZE);
          grps[groupname].push_back(famName);
        }
    }
}

void MEDFileMeshL2::WriteFamiliesAndGrps(med_idt fid, const char *mname, const std::map<std::string,int>& fams, const std::map<std::string, std::vector<std::string> >& grps, int tooLongStrPol)
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
      INTERP_KERNEL::AutoPtr<char> groName=MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE*ngro);
      int i=0;
      for(std::vector<std::string>::const_iterator it2=grpsOfFam.begin();it2!=grpsOfFam.end();it2++,i++)
        MEDLoaderBase::safeStrCpy2((*it2).c_str(),MED_LNAME_SIZE-1,groName+i*MED_LNAME_SIZE,tooLongStrPol);
      INTERP_KERNEL::AutoPtr<char> famName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
      MEDLoaderBase::safeStrCpy((*it).first.c_str(),MED_NAME_SIZE,famName,tooLongStrPol);
      int ret=MEDfamilyCr(fid,mname,famName,(*it).second,ngro,groName);
      ret++;
    }
}

MEDFileUMeshL2::MEDFileUMeshL2()
{
}

void MEDFileUMeshL2::loadAll(med_idt fid, int mId, const char *mName, int dt, int it)
{
  _name.set(mName);
  int nstep;
  int Mdim;
  ParaMEDMEM::MEDCouplingMeshType meshType;
  std::vector<std::string> infosOnComp=getAxisInfoOnMesh(fid,mId,mName,meshType,nstep,Mdim);
  if(meshType!=UNSTRUCTURED)
    throw INTERP_KERNEL::Exception("Invalid mesh type ! You are expected an unstructured one whereas in file it is not an unstructured !");
  _time=CheckMeshTimeStep(fid,mName,nstep,dt,it);
  _iteration=dt;
  _order=it;
  loadConnectivity(fid,Mdim,mName,dt,it);//to improve check (dt,it) coherency
  loadCoords(fid,mId,infosOnComp,mName,dt,it);
}

void MEDFileUMeshL2::loadConnectivity(med_idt fid, int mdim, const char *mName, int dt, int it)
{
  _per_type_mesh.resize(1);
  _per_type_mesh[0].clear();
  for(int j=0;j<MED_N_CELL_FIXED_GEO;j++)
    {
      MEDFileUMeshPerType *tmp=MEDFileUMeshPerType::New(fid,mName,dt,it,mdim,typmai[j],typmai2[j]);
      if(tmp)
        _per_type_mesh[0].push_back(tmp);
    }
  sortTypes();
}

void MEDFileUMeshL2::loadCoords(med_idt fid, int mId, const std::vector<std::string>& infosOnComp, const char *mName, int dt, int it) throw(INTERP_KERNEL::Exception)
{
  int spaceDim=infosOnComp.size();
  med_bool changement,transformation;
  int nCoords=MEDmeshnEntity(fid,mName,dt,it,MED_NODE,MED_NONE,MED_COORDINATE,MED_NO_CMODE,&changement,&transformation);
  _coords=DataArrayDouble::New();
  _coords->alloc(nCoords,spaceDim);
  double *coordsPtr=_coords->getPointer();
  MEDmeshNodeCoordinateRd(fid,mName,dt,it,MED_FULL_INTERLACE,coordsPtr);
  _fam_coords=DataArrayInt::New();
  _fam_coords->alloc(nCoords,1);
  _num_coords=DataArrayInt::New();
  _num_coords->alloc(nCoords,1);
  if(MEDmeshnEntity(fid,mName,dt,it,MED_NODE,MED_NO_GEOTYPE,MED_FAMILY_NUMBER,MED_NODAL,&changement,&transformation)>0)
    MEDmeshEntityFamilyNumberRd(fid,mName,dt,it,MED_NODE,MED_NO_GEOTYPE,_fam_coords->getPointer());
  else
    _fam_coords=0;
  if(MEDmeshnEntity(fid,mName,dt,it,MED_NODE,MED_NO_GEOTYPE,MED_NUMBER,MED_NODAL,&changement,&transformation)>0)
    MEDmeshEntityNumberRd(fid,mName,dt,it,MED_NODE,MED_NO_GEOTYPE,_num_coords->getPointer());
  else
    _num_coords=0;
  for(int i=0;i<spaceDim;i++)
    _coords->setInfoOnComponent(i,infosOnComp[i].c_str());
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
  // suppression of contiguous empty levels at the end of _per_type_mesh.
  int nbOfUselessLev=0;
  bool isFirst=true;
  for(std::vector< std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshPerType> > >::reverse_iterator it2=_per_type_mesh.rbegin();it2!=_per_type_mesh.rend();it2++)
    {
      if((*it2).empty() && isFirst)
        {
          nbOfUselessLev++;
        }
      else
        isFirst=false;
    }
  _per_type_mesh.resize(_per_type_mesh.size()-nbOfUselessLev);
}

void MEDFileUMeshL2::WriteCoords(med_idt fid, const char *mname, int dt, int it, double time, const DataArrayDouble *coords, const DataArrayInt *famCoords, const DataArrayInt *numCoords)
{
  if(!coords)
    return ;
  MEDmeshNodeCoordinateWr(fid,mname,dt,it,time,MED_FULL_INTERLACE,coords->getNumberOfTuples(),coords->getConstPointer());
  if(famCoords)
    MEDmeshEntityFamilyNumberWr(fid,mname,dt,it,MED_NODE,MED_NO_GEOTYPE,famCoords->getNumberOfTuples(),famCoords->getConstPointer());
  if(numCoords)
    MEDmeshEntityNumberWr(fid,mname,dt,it,MED_NODE,MED_NO_GEOTYPE,numCoords->getNumberOfTuples(),numCoords->getConstPointer());
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

MEDFileCMeshL2::MEDFileCMeshL2()
{
}

void MEDFileCMeshL2::loadAll(med_idt fid, int mId, const char *mName, int dt, int it) throw(INTERP_KERNEL::Exception)
{
  _name.set(mName);
  int nstep;
  int Mdim;
  ParaMEDMEM::MEDCouplingMeshType meshType;
  std::vector<std::string> infosOnComp=getAxisInfoOnMesh(fid,mId,mName,meshType,nstep,Mdim);
  if(meshType!=CARTESIAN)
    throw INTERP_KERNEL::Exception("Invalid mesh type ! You are expected a structured one whereas in file it is not a structured !");
  _time=CheckMeshTimeStep(fid,mName,nstep,dt,it);
  _iteration=dt;
  _order=it;
  //
  med_grid_type gridtype;
  MEDmeshGridTypeRd(fid,mName,&gridtype);
  if(gridtype!=MED_CARTESIAN_GRID)
    throw INTERP_KERNEL::Exception("Invalid cartesion mesh type ! Only Cartesian Grid supported ! Curvilinear grid will come soon !");
  _cmesh=MEDCouplingCMesh::New();
  for(int i=0;i<Mdim;i++)
    {
      med_data_type dataTypeReq=GetDataTypeCorrespondingToSpaceId(i);
      med_bool chgt=MED_FALSE,trsf=MED_FALSE;
      int nbOfElt=MEDmeshnEntity(fid,mName,dt,it,MED_NODE,MED_NONE,dataTypeReq,MED_NO_CMODE,&chgt,&trsf);
      MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> da=DataArrayDouble::New();
      da->alloc(nbOfElt,1);
      da->setInfoOnComponent(0,infosOnComp[i].c_str());
      MEDmeshGridIndexCoordinateRd(fid,mName,dt,it,i+1,da->getPointer());
      _cmesh->setCoordsAt(i,da);
    }
}

med_data_type MEDFileCMeshL2::GetDataTypeCorrespondingToSpaceId(int id) throw(INTERP_KERNEL::Exception)
{
  switch(id)
    {
    case 0:
      return MED_COORDINATE_AXIS1;
    case 1:
      return MED_COORDINATE_AXIS2;
    case 2:
      return MED_COORDINATE_AXIS3;
    default:
      throw INTERP_KERNEL::Exception("Invalid meshdim detected in Cartesian Grid !");
    }
}

MEDFileUMeshPermCompute::MEDFileUMeshPermCompute(const MEDFileUMeshSplitL1* st):_st(st),_mpt_time(0),_num_time(0)
{
}

/*!
 * Warning it returns an instance to deallocate !!!!
 */
MEDFileUMeshPermCompute::operator MEDCouplingUMesh *() const
{
  _st->_m_by_types->updateTime();
  _st->_num->updateTime();
  if((MEDCouplingUMesh *)_m==0)
    {
      updateTime();
      MEDCouplingUMesh *ret=(MEDCouplingUMesh *)_st->_m_by_types->deepCpy();
      _m=ret;
      _m->renumberCells(_st->_num->getConstPointer(),true);
      ret->incrRef();
      return ret;
    }
  else
    {
      if(_mpt_time==_st->_m_by_types->getTimeOfThis() && _num_time==_st->_num->getTimeOfThis())
        {
          _m->incrRef();
          return _m;
        }
      else
        {
          updateTime();
          MEDCouplingUMesh *ret=(MEDCouplingUMesh *)_st->_m_by_types->deepCpy();
          _m=ret;
          _m->renumberCells(_st->_num->getConstPointer(),true);
          ret->incrRef();
          return ret;
        }
    }
}

void MEDFileUMeshPermCompute::operator=(MEDCouplingUMesh *m)
{
  _m=m;
}

void MEDFileUMeshPermCompute::updateTime() const
{
  _mpt_time=_st->_m_by_types->getTimeOfThis();
  _num_time=_st->_num->getTimeOfThis();
}

MEDFileUMeshSplitL1::MEDFileUMeshSplitL1(const MEDFileUMeshL2& l2, const char *mName, int id):_m(this)
{
  const std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshPerType> >& v=l2.getLev(id);
  if(v.empty())
    return;
  int sz=v.size();
  std::vector<const MEDCouplingUMesh *> ms(sz);
  for(int i=0;i<sz;i++)
    {
      MEDCouplingUMesh *tmp=MEDCouplingUMesh::New("",v[i]->getDim());
      MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> tmp2=l2.getCoords();
      tmp->setCoords(tmp2);
      tmp->setConnectivity(const_cast<DataArrayInt *>(v[i]->getNodal()),const_cast<DataArrayInt *>(v[i]->getNodalIndex()));
      ms[i]=tmp;
    }
  _m_by_types=MEDCouplingUMesh::MergeUMeshesOnSameCoords(ms);
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
      computeRevNum();
    }
  for(int i=0;i<sz;i++)
    (const_cast<MEDCouplingUMesh *>(ms[i]))->decrRef();//const cast under control to avoid a copy of array
}

MEDFileUMeshSplitL1::MEDFileUMeshSplitL1(MEDCouplingUMesh *m):_m(this)
{
  assignMesh(m,true);
}

MEDFileUMeshSplitL1::MEDFileUMeshSplitL1(MEDCouplingUMesh *m, bool newOrOld):_m(this)
{
  assignMesh(m,newOrOld);
}

bool MEDFileUMeshSplitL1::isEqual(const MEDFileUMeshSplitL1 *other, double eps, std::string& what) const
{
  const MEDCouplingUMesh *m1=_m_by_types;
  const MEDCouplingUMesh *m2=other->_m_by_types;
  if((m1==0 && m2!=0) || (m1!=0 && m2==0))
    {
      what="Presence of mesh in one sublevel and not in other!";
      return false;
    }
  if(m1)
    if(!m1->isEqual(m2,eps))
      {
        what="meshes at a sublevel are not deeply equal !";
        return false;
      }
  const DataArrayInt *d1=_fam;
  const DataArrayInt *d2=other->_fam;
  if((d1==0 && d2!=0) || (d1!=0 && d2==0))
    {
      what="Presence of family arr in one sublevel and not in other!";
      return false;
    }
  if(d1)
    if(!d1->isEqual(*d2))
      {
        what="family arr at a sublevel are not deeply equal !";
        return false;
      }
  d1=_num;
  d2=other->_num;
  if((d1==0 && d2!=0) || (d1!=0 && d2==0))
    {
      what="Presence of cell numbering arr in one sublevel and not in other!";
      return false;
    }
  if(d1)
    if(!d1->isEqual(*d2))
      {
        what="Numbering cell arr at a sublevel are not deeply equal !";
        return false;
      }
  return true;
}

void MEDFileUMeshSplitL1::synchronizeTinyInfo(const MEDFileMesh& master) const
{
  const MEDCouplingUMesh *tmp=_m_by_types;
  if(!tmp)
    return ;
  (const_cast<MEDCouplingUMesh *>(tmp))->setName(master.getName());
  (const_cast<MEDCouplingUMesh *>(tmp))->setDescription(master.getDescription());
  (const_cast<MEDCouplingUMesh *>(tmp))->setTime(master.getTimeValue(),master.getIteration(),master.getOrder());
  (const_cast<MEDCouplingUMesh *>(tmp))->setTimeUnit(master.getTimeUnit());
}

void MEDFileUMeshSplitL1::clearNonDiscrAttributes() const
{
  ClearNonDiscrAttributes(_m_by_types);
}

void MEDFileUMeshSplitL1::ClearNonDiscrAttributes(const MEDCouplingMesh *tmp)
{
  if(!tmp)
    return ;
  (const_cast<MEDCouplingMesh *>(tmp))->setName("");
  (const_cast<MEDCouplingMesh *>(tmp))->setDescription("");
  (const_cast<MEDCouplingMesh *>(tmp))->setTime(0.,-1,-1);
  (const_cast<MEDCouplingMesh *>(tmp))->setTimeUnit("");
}

void MEDFileUMeshSplitL1::assignMesh(MEDCouplingUMesh *m, bool newOrOld) throw(INTERP_KERNEL::Exception)
{
  if(newOrOld)
    {
      m->incrRef();
      _m=m;
      _m_by_types=(MEDCouplingUMesh *)m->deepCpy();
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da=_m_by_types->getRenumArrForConsecutiveCellTypesSpec(typmai2,typmai2+MED_N_CELL_FIXED_GEO);
      if(!da->isIdentity())
        {
          _num=da->invertArrayO2N2N2O(m->getNumberOfCells());
          _m.updateTime();
          computeRevNum();
          _m_by_types->renumberCells(da->getConstPointer(),false);
        }
    }
  else
    {
      if(!m->checkConsecutiveCellTypesAndOrder(typmai2,typmai2+MED_N_CELL_FIXED_GEO))
        throw INTERP_KERNEL::Exception("MEDFileUMeshSplitL1::assignMesh : the mode of mesh setting expects to follow the MED file numbering convention ! it is not the case !");
      m->incrRef();
      _m_by_types=m;
    }
  _fam=DataArrayInt::New();
  _fam->alloc(m->getNumberOfCells(),1);
  _fam->fillWithValue(0);
}

bool MEDFileUMeshSplitL1::empty() const
{
  return ((const MEDCouplingUMesh *)_m_by_types)==0;
}

bool MEDFileUMeshSplitL1::presenceOfOneFams(const std::vector<int>& ids) const
{
  const DataArrayInt *fam=_fam;
  if(!fam)
    return false;
  return fam->presenceOfValue(ids);
}

int MEDFileUMeshSplitL1::getMeshDimension() const
{
  return _m_by_types->getMeshDimension();
}

void MEDFileUMeshSplitL1::simpleRepr(std::ostream& oss) const
{
  std::vector<int> code=_m_by_types->getDistributionOfTypes();
  int nbOfTypes=code.size()/3;
  for(int i=0;i<nbOfTypes;i++)
    {
      INTERP_KERNEL::NormalizedCellType typ=(INTERP_KERNEL::NormalizedCellType) code[3*i];
      oss << "    - Number of cells with type " << INTERP_KERNEL::CellModel::GetCellModel(typ).getRepr() << " : " << code[3*i+1] << std::endl;
    }
}

int MEDFileUMeshSplitL1::getSize() const throw(INTERP_KERNEL::Exception)
{
  if((const MEDCouplingUMesh *)_m_by_types==0)
    throw INTERP_KERNEL::Exception("MEDFileUMeshSplitL1::getSize : no mesh specified at level !");
  return _m_by_types->getNumberOfCells();
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
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da=_fam->getIdsEqualList(ids);
  if(renum)
    return renumIfNeededArr(da);
  da->incrRef();
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

const DataArrayInt *MEDFileUMeshSplitL1::getFamilyField() const
{
  return _fam;
}

const DataArrayInt *MEDFileUMeshSplitL1::getNumberField() const
{
  return _num;
}

const DataArrayInt *MEDFileUMeshSplitL1::getRevNumberField() const
{
  return _rev_num;
}

void MEDFileUMeshSplitL1::eraseFamilyField()
{
  _fam->fillWithZero();
}

/*!
 * This method ignores _m and _m_by_types.
 */
void MEDFileUMeshSplitL1::setGroupsFromScratch(const std::vector<const MEDCouplingUMesh *>& ms, std::map<std::string,int>& familyIds,
                                               std::map<std::string, std::vector<std::string> >& groups) throw(INTERP_KERNEL::Exception)
{
  int sz=ms.size();
  std::vector< DataArrayInt * > corr;
  _m=MEDCouplingUMesh::FuseUMeshesOnSameCoords(ms,0,corr);
  std::vector< std::vector<int> > fidsOfGroups;
  std::vector< const DataArrayInt * > corr2(corr.begin(),corr.end());
  _fam=DataArrayInt::MakePartition(corr2,((MEDCouplingUMesh *)_m)->getNumberOfCells(),fidsOfGroups);
  int nbOfCells=((MEDCouplingUMesh *)_m)->getNumberOfCells();
  std::map<int,std::string> newfams;
  std::map<int,int> famIdTrad;
  TraduceFamilyNumber(fidsOfGroups,familyIds,famIdTrad,newfams);
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

void MEDFileUMeshSplitL1::changeFamilyIdArr(int oldId, int newId) throw(INTERP_KERNEL::Exception)
{
  DataArrayInt *arr=_fam;
  if(arr)
    arr->changeValue(oldId,newId);
}

void MEDFileUMeshSplitL1::setFamilyArr(DataArrayInt *famArr)
{
  famArr->incrRef();
  _fam=famArr;
}

void MEDFileUMeshSplitL1::setRenumArr(DataArrayInt *renumArr)
{
  renumArr->incrRef();
  _num=renumArr;
  computeRevNum();
}

MEDCouplingUMesh *MEDFileUMeshSplitL1::Renumber2(const DataArrayInt *renum, MEDCouplingUMesh *m, const int *cellIds)
{
  if(renum==0)
    return m;
  if(cellIds==0)
    m->renumberCells(renum->getConstPointer(),true);
  else
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> locnum=renum->selectByTupleId(cellIds,cellIds+m->getNumberOfCells());
      m->renumberCells(locnum->getConstPointer(),true);
    }
  return m;
}

MEDCouplingUMesh *MEDFileUMeshSplitL1::renumIfNeeded(MEDCouplingUMesh *m, const int *cellIds) const
{
  return Renumber2(_num,m,cellIds);
}

DataArrayInt *MEDFileUMeshSplitL1::Renumber(const DataArrayInt *renum, const DataArrayInt *da)
{
  if((const DataArrayInt *)renum==0)
    {
      da->incrRef();
      return const_cast<DataArrayInt *>(da);
    }
  return renum->selectByTupleId(da->getConstPointer(),da->getConstPointer()+da->getNumberOfTuples());
}

DataArrayInt *MEDFileUMeshSplitL1::renumIfNeededArr(const DataArrayInt *da) const
{
  return Renumber(_num,da);
}

std::vector<int> MEDFileUMeshSplitL1::GetNewFamiliesNumber(int nb, const std::map<std::string,int>& families)
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

void MEDFileUMeshSplitL1::TraduceFamilyNumber(const std::vector< std::vector<int> >& fidsGrps, std::map<std::string,int>& familyIds,
                                              std::map<int,int>& famIdTrad, std::map<int,std::string>& newfams)
{
  std::set<int> allfids;
  
}

void MEDFileUMeshSplitL1::computeRevNum() const
{
  int pos;
  int maxValue=_num->getMaxValue(pos);
  _rev_num=_num->invertArrayN2O2O2N(maxValue+1);
}
