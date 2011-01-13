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

#include "MEDFileField.hxx"
#include "MEDLoaderBase.hxx"
#include "MEDFileUtilities.hxx"

#include "InterpKernelAutoPtr.hxx"

#include <algorithm>

extern med_geometrie_element typmai[MED_NBR_GEOMETRIE_MAILLE+2];
extern INTERP_KERNEL::NormalizedCellType typmai2[MED_NBR_GEOMETRIE_MAILLE+2];
extern med_geometrie_element typmainoeud[1];

using namespace ParaMEDMEM;

MEDFileFieldPerMeshPerTypePerDisc *MEDFileFieldPerMeshPerTypePerDisc::New(MEDFileFieldPerMeshPerType *fath, med_idt fid) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldPerMeshPerTypePerDisc(fath,fid);
}

MEDFileFieldPerMeshPerTypePerDisc::MEDFileFieldPerMeshPerTypePerDisc(MEDFileFieldPerMeshPerType *fath, med_idt fid) throw(INTERP_KERNEL::Exception)
try:_father(fath)
{
  //to implement
}
catch(INTERP_KERNEL::Exception& e)
{
  throw e;
}

const MEDFileFieldPerMeshPerType *MEDFileFieldPerMeshPerTypePerDisc::getFather() const
{
  return _father;
}

MEDFileFieldPerMeshPerType *MEDFileFieldPerMeshPerType::New(MEDFileFieldPerMesh *fath, TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldPerMeshPerType(fath,type,geoType);
}

const MEDFileFieldPerMesh *MEDFileFieldPerMeshPerType::getFather() const
{
  return _father;
}

MEDFileFieldPerMeshPerType::MEDFileFieldPerMeshPerType(MEDFileFieldPerMesh *fath, TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception):_father(fath),_type(type),_geo_type(geoType)
{
}

void MEDFileFieldPerMeshPerType::finishLoading(med_idt fid) throw(INTERP_KERNEL::Exception)
{
  //for MED3 porting this limitation will be killed !
  _field_pm_pt_pd.resize(1);
  _field_pm_pt_pd[0]=MEDFileFieldPerMeshPerTypePerDisc::New(this,fid);
}

MEDFileFieldPerMesh *MEDFileFieldPerMesh::New(MEDFileField1TSWithoutDAS *fath, const char *meshName, double time)
{
  return new MEDFileFieldPerMesh(fath,meshName,time);
}

void MEDFileFieldPerMesh::pushBack(TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType)
{
  _types.push_back(type);
  _geo_types.push_back(geoType);
}

void MEDFileFieldPerMesh::finishLoading(med_idt fid) throw(INTERP_KERNEL::Exception)
{
  int sz=_types.size();
  _field_pm_pt.resize(sz);
  for(int i=0;i<sz;i++)
    {
      _field_pm_pt[i]=MEDFileFieldPerMeshPerType::New(this,_types[i],_geo_types[i]);
      _field_pm_pt[i]->finishLoading(fid);
    }
}

std::string MEDFileFieldPerMesh::getMeshName() const
{
  return _mesh_name;
}

double MEDFileFieldPerMesh::getTime() const
{
  return _time;
}

MEDFileFieldPerMesh::MEDFileFieldPerMesh(MEDFileField1TSWithoutDAS *fath, const char *meshName, double time):_father(fath),_mesh_name(meshName),_time(time)
{
}

MEDFileField1TSWithoutDAS *MEDFileField1TSWithoutDAS::New(const char *fieldName, int iteration, int order)
{
  return new MEDFileField1TSWithoutDAS(fieldName,iteration,order);
}

void MEDFileField1TSWithoutDAS::pushBack(TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType, double time, const char *meshName) const
{
  _types.push_back(type);
  _geo_types.push_back(geoType);
  _times.push_back(time);
  _meshes.push_back(meshName);
}

void MEDFileField1TSWithoutDAS::finishLoading(med_idt fid) throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> meshesOnlyOnce;
  int nbOfTurn=_meshes.size();
  for(int i=0;i<nbOfTurn;i++)
    {
      std::vector<std::string>::iterator it=std::find(meshesOnlyOnce.begin(),meshesOnlyOnce.end(),_meshes[i]);
      if(it==meshesOnlyOnce.end())
        {
          _field_per_mesh.push_back(MEDFileFieldPerMesh::New(this,_meshes[i].c_str(),_times[i]));
          meshesOnlyOnce.push_back(_meshes[i]);
          _field_per_mesh.back()->pushBack(_types[i],_geo_types[i]);
        }
      else
        {
          int w=std::distance(meshesOnlyOnce.begin(),it);
          _field_per_mesh[w]->pushBack(_types[i],_geo_types[i]);
        }
    }
  int nbOfMeshes=_field_per_mesh.size();
  for(int i=0;i<nbOfMeshes;i++)
    _field_per_mesh[i]->finishLoading(fid);
}

MEDFileField1TSWithoutDAS::MEDFileField1TSWithoutDAS(const char *fieldName, int iteration, int order):_name(fieldName),_iteration(iteration),_order(order)
{
}

MEDFileFieldMultiTSWithoutDAS *MEDFileFieldMultiTSWithoutDAS::New(med_idt fid, const char *fieldName, int id, const std::vector<std::string>& infos) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldMultiTSWithoutDAS(fid,fieldName,id,infos);
}

MEDFileFieldMultiTSWithoutDAS::MEDFileFieldMultiTSWithoutDAS(med_idt fid, const char *fieldName, int id,  const std::vector<std::string>& infos) throw(INTERP_KERNEL::Exception)
try:_name(fieldName)
{
  for(int i=0;i<MED_NBR_GEOMETRIE_MAILLE+2;i++)
    {
      int nbPdt=MEDnPasdetemps(fid,(char *)fieldName,MED_MAILLE,typmai[i]);
      for(int j=0;j<nbPdt;j++)
        appendTimeStepEntry(fid,fieldName,MED_MAILLE,i,j);
      nbPdt=MEDnPasdetemps(fid,(char *)fieldName,MED_NOEUD,MED_NONE);
      for(int j=0;j<nbPdt;j++)
        appendTimeStepEntry(fid,fieldName,MED_NOEUD,i,j);
    }
  int nbOfTimeSteps=_time_steps.size();
  for(int i=0;i<nbOfTimeSteps;i++)
    _time_steps[i]->finishLoading(fid);
}
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

void MEDFileFieldMultiTSWithoutDAS::appendTimeStepEntry(med_idt fid, const char *fieldName, med_entite_maillage entity, int i, int j) throw(INTERP_KERNEL::Exception)
{
  std::vector< std::pair<int,int> > ts;
  med_int ngauss=0;
  med_int numdt=0,numo=0,nbrefmaa;
  med_float dt=0.0;
  med_booleen local;
  INTERP_KERNEL::AutoPtr<char> maa_ass=MEDLoaderBase::buildEmptyString(MED_TAILLE_NOM);
  INTERP_KERNEL::AutoPtr<char> dt_unit=MEDLoaderBase::buildEmptyString(MED_TAILLE_PNOM);
  MEDpasdetempsInfo(fid,(char *)fieldName,MED_MAILLE,typmai[i],j+1,&ngauss,&numdt,&numo,dt_unit,&dt,maa_ass,&local,&nbrefmaa);
  std::pair<int,int> p(numdt,numo);
  std::vector< std::pair<int,int> >::iterator where=std::find(ts.begin(),ts.end(),p);
  if(where==ts.end())
    {
      ts.push_back(p);
      _time_steps.push_back(MEDFileField1TSWithoutDAS::New(fieldName,numdt,numo));
      _time_steps.back()->pushBack(entity==MED_MAILLE?ON_CELLS:ON_NODES,entity==MED_MAILLE?typmai2[i]:INTERP_KERNEL::NORM_ERROR,dt,maa_ass);
    }
  else
    {
      int w=std::distance(ts.begin(),where);
      _time_steps[w]->pushBack(entity==MED_MAILLE?ON_CELLS:ON_NODES,entity==MED_MAILLE?typmai2[i]:INTERP_KERNEL::NORM_ERROR,dt,maa_ass);
    }
}

MEDFileFields *MEDFileFields::New(const char *fileName) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFields(fileName);
}

int MEDFileFields::getNumberOfFields() const
{
  return _fields.size();
}

MEDFileFields::MEDFileFields(const char *fileName) throw(INTERP_KERNEL::Exception)
try
  {
    MEDFileUtilities::CheckFileForRead(fileName);
    med_idt fid=MEDouvrir((char *)fileName,MED_LECTURE);
    int nbFields=MEDnChamp(fid,0);
    _fields.resize(nbFields);
    INTERP_KERNEL::AutoPtr<char> nomcha=MEDLoaderBase::buildEmptyString(MED_TAILLE_NOM);
    med_type_champ typcha;
    for(int i=0;i<nbFields;i++)
      {
        int ncomp=MEDnChamp(fid,i+1);
        INTERP_KERNEL::AutoPtr<char> comp=new char[ncomp*MED_TAILLE_PNOM+1];
        INTERP_KERNEL::AutoPtr<char> unit=new char[ncomp*MED_TAILLE_PNOM+1];
        MEDchampInfo(fid,i+1,nomcha,&typcha,comp,unit,ncomp);
        std::vector<std::string> infos(ncomp);
        for(int j=0;j<ncomp;j++)
          infos[j]=MEDLoaderBase::buildUnionUnit((char *)comp+i*MED_TAILLE_PNOM,MED_TAILLE_PNOM,(char *)unit+i*MED_TAILLE_PNOM,MED_TAILLE_PNOM);
        _fields[i]=MEDFileFieldMultiTSWithoutDAS::New(fid,nomcha,i,infos);
      }
    MEDfermer(fid);
  }
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }
