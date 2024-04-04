// Copyright (C) 2007-2024  CEA, EDF
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

#include "MEDFileMeshLL.hxx"
#include "MEDCouplingCurveLinearMesh.hxx"
#include "MEDCouplingMesh.hxx"
#include "MEDCouplingMemArray.txx"
#include "MEDFileBasis.hxx"
#include "MCType.hxx"
#include "MCAuto.hxx"
#include "MEDCouplingMap.txx"
#include "MEDCouplingPartDefinition.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "MCIdType.hxx"
#include "MEDCoupling1GTUMesh.hxx"
#include "MEDFileMesh.hxx"
#include "MEDFileMeshElt.hxx"
#include "MEDLoaderBase.hxx"
#include "MEDFileSafeCaller.txx"
#include "MEDFileMeshReadSelector.hxx"
#include "MEDFileStructureElement.hxx"
#include "MEDFileMeshSupport.hxx"

#include "MEDCouplingUMesh.hxx"

#include "InterpKernelAutoPtr.hxx"
#include "CellModel.hxx"

#include "MEDFilterEntity.hxx"
#include "med.h"
#include "NormalizedGeometricTypes"
#include "medmesh.h"
#include <algorithm>
#include <cstddef>
#include "medfamily.h"
#include "medstructelement.h"
#include <istream>
#include <iterator>
#include <map>
#include <ostream>
#include <set>
#include <string>
#include <sstream>
#include <vector>
#include <utility>
#include <iomanip>

// From MEDLOader.cxx TU
extern med_geometry_type typmai[MED_N_CELL_FIXED_GEO];
extern INTERP_KERNEL::NormalizedCellType typmai2[MED_N_CELL_FIXED_GEO];
extern med_geometry_type typmainoeud[1];

using namespace MEDCoupling;

const char MEDFileMeshL2::ZE_SEP_FOR_FAMILY_KILLERS[]="!/__\\!";//important start by - because ord('!')==33 the smallest (!=' ') to preserve orders at most.

int MEDFileMeshL2::ZE_SEP2_FOR_FAMILY_KILLERS=4;

std::vector<std::string> MeshCls::getAxisInfoOnMesh(med_idt fid, const std::string& mName, MEDCoupling::MEDCouplingMeshType& meshType, MEDCoupling::MEDCouplingAxisType& axType, int& nstep, int& Mdim, MEDFileString& description, MEDFileString& dtunit, MEDFileString& univName) const
{
  med_mesh_type type_maillage;
  med_int spaceDim, meshDim, nbSteps;
  med_sorting_type stype;
  med_axis_type axistype;
  med_int const naxis(MEDmeshnAxis(fid,getID()));
  INTERP_KERNEL::AutoPtr<char> nameTmp(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> axisname(MEDLoaderBase::buildEmptyString(naxis*MED_SNAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> axisunit(MEDLoaderBase::buildEmptyString(naxis*MED_SNAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> const univTmp(MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE));
  if(MEDmeshInfo(fid,getID(),nameTmp,&spaceDim,&meshDim,&type_maillage,description.getPointer(),dtunit.getPointer(),
      &stype,&nbSteps,&axistype,axisname,axisunit)!=0)
    throw INTERP_KERNEL::Exception("A problem has been detected when trying to get info on mesh !");
  Mdim=FromMedInt<int>(meshDim);
  nstep=FromMedInt<int>(nbSteps);
  MEDmeshUniversalNameRd(fid,nameTmp,univName.getPointer());// do not protect  MEDFILESAFECALLERRD0 call : Thanks to fra.med.
  axType=MEDFileMeshL2::TraduceAxisType(axistype);
  switch(type_maillage)
  {
    case MED_UNSTRUCTURED_MESH:
      meshType=UNSTRUCTURED;
      break;
    case MED_STRUCTURED_MESH:
      {
        med_grid_type gt;
        MEDFILESAFECALLERRD0(MEDmeshGridTypeRd,(fid,mName.c_str(),&gt));
        switch(gt)
        {
          case MED_CARTESIAN_GRID:
            meshType=CARTESIAN;
            break;
          case MED_CURVILINEAR_GRID:
            meshType=CURVE_LINEAR;
            break;
        case MED_POLAR_GRID:// this is not a bug. A MED file POLAR_GRID is deal by CARTESIAN MEDLoader
            meshType=CARTESIAN;
            break;
          default:
            throw INTERP_KERNEL::Exception("MEDFileMeshL2::getAxisInfoOnMesh : unrecognized structured mesh type ! Supported are :\n - cartesian\n - curve linear\n");
        }
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDFileMeshL2::getMeshIdFromName : unrecognized mesh type !");
  }
  //
  std::vector<std::string> infosOnComp(naxis);
  for(int i=0;i<naxis;i++)
    {
      std::string const info(MEDLoaderBase::buildUnionUnit(((char *)axisname)+i*MED_SNAME_SIZE,MED_SNAME_SIZE,((char *)axisunit)+i*MED_SNAME_SIZE,MED_SNAME_SIZE));
      infosOnComp[i]=info;
    }
  return infosOnComp;
}

double MeshCls::checkMeshTimeStep(med_idt fid, const std::string& mName, int nstep, int dt, int it) const
{
  bool found=false;
  med_int numdt,numit;
  med_float dtt=0.0;
  std::vector< std::pair<int,int> > p(nstep);
  for(int i=0;i<nstep;i++)
    {
      MEDFILESAFECALLERRD0(MEDmeshComputationStepInfo,(fid,mName.c_str(),i+1,&numdt,&numit,&dtt));
      p[i]=std::make_pair((int)numdt,(int)numit);
      found=(numdt==dt) && (numit==it);
      if (found) break;
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

std::vector<std::string> StructMeshCls::getAxisInfoOnMesh(med_idt fid, const std::string&  /*mName*/, MEDCoupling::MEDCouplingMeshType& meshType, MEDCoupling::MEDCouplingAxisType& axType, int& nstep, int&  /*Mdim*/, MEDFileString& description, MEDFileString& dtunit, MEDFileString& univName) const
{
  INTERP_KERNEL::AutoPtr<char> msn(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> zeDescription(MEDLoaderBase::buildEmptyString(MED_COMMENT_SIZE));
  med_axis_type medAxType;
  med_int const nAxis(MEDsupportMeshnAxis(fid,getID()));
  INTERP_KERNEL::AutoPtr<char> axisName(new char[MED_SNAME_SIZE*nAxis+1]),axisUnit(new char[MED_SNAME_SIZE*nAxis+1]);
  med_int spaceDim(0),meshDim(0);
  MEDFILESAFECALLERRD0(MEDsupportMeshInfo,(fid,getID(),msn,&spaceDim,&meshDim,zeDescription,&medAxType,axisName,axisUnit));
  std::string const descriptionCpp(MEDLoaderBase::buildStringFromFortran(zeDescription,MED_COMMENT_SIZE));
  description.set(descriptionCpp.c_str());
  dtunit.clear(); univName.clear(); meshType=UNSTRUCTURED; nstep=1;
  axType=MEDFileMeshL2::TraduceAxisType(medAxType);
  //int nmodels(0);
  //med_bool chgt=MED_FALSE,trsf=MED_FALSE;
  //nmodels=MEDmeshnEntity(fid,_name.c_str(),MED_NO_DT,MED_NO_IT,MED_STRUCT_ELEMENT,MED_GEO_ALL,MED_CONNECTIVITY,MED_NODAL,&chgt,&trsf);
  std::vector<std::string> ret;
  for(int i=0;i<nAxis;i++)
    {
      std::string const info(DataArray::BuildInfoFromVarAndUnit(MEDLoaderBase::buildStringFromFortran(axisName+i*MED_SNAME_SIZE,MED_SNAME_SIZE),
                                                          MEDLoaderBase::buildStringFromFortran(axisUnit+i*MED_SNAME_SIZE,MED_SNAME_SIZE)));
      ret.push_back(info);
    }
  return ret;
}

double StructMeshCls::checkMeshTimeStep(med_idt  /*fid*/, const std::string&  /*mName*/, int  /*nstep*/, int  /*dt*/, int  /*it*/) const
{
  return 0.;
}

MEDFileMeshL2::MEDFileMeshL2():_name(MED_NAME_SIZE),_description(MED_COMMENT_SIZE),_univ_name(MED_LNAME_SIZE),_dt_unit(MED_LNAME_SIZE)
{
}

std::size_t MEDFileMeshL2::getHeapMemorySizeWithoutChildren() const
{
  return 0;
}

std::vector<const BigMemoryObject *> MEDFileMeshL2::getDirectChildrenWithNull() const
{
  return std::vector<const BigMemoryObject *>();
}

INTERP_KERNEL::AutoCppPtr<MeshOrStructMeshCls> MEDFileMeshL2::GetMeshIdFromName(med_idt fid, const std::string& mName, MEDCoupling::MEDCouplingMeshType& meshType, MEDCoupling::MEDCouplingAxisType& axType, int& dt, int& it, std::string& dtunit1)
{
  med_mesh_type type_maillage=MED_UNDEF_MESH_TYPE;
  char maillage_description[MED_COMMENT_SIZE+1];
  char dtunit[MED_LNAME_SIZE+1];
  med_int spaceDim,dim;
  char nommaa[MED_NAME_SIZE+1];
  med_int const n=MEDnMesh(fid);
  char found(0);
  int ret=-1;
  med_sorting_type stype;
  std::vector<std::string> ms;
  med_int nstep;
  med_axis_type axistype=MED_UNDEF_AXIS_TYPE;
  for(int i=0;i<n && found==0;i++)
    {
      med_int const naxis(MEDmeshnAxis(fid,i+1));
      INTERP_KERNEL::AutoPtr<char> axisname(MEDLoaderBase::buildEmptyString(naxis*MED_SNAME_SIZE)),axisunit(MEDLoaderBase::buildEmptyString(naxis*MED_SNAME_SIZE));
      MEDFILESAFECALLERRD0(MEDmeshInfo,(fid,i+1,nommaa,&spaceDim,&dim,&type_maillage,maillage_description,dtunit,&stype,&nstep,&axistype,axisname,axisunit));      
      dtunit1=MEDLoaderBase::buildStringFromFortran(dtunit,sizeof(dtunit));
      std::string const cur(MEDLoaderBase::buildStringFromFortran(nommaa,sizeof(nommaa)));
      ms.push_back(cur);
      if(cur==mName)
        {
          found=1;
          ret=i+1;
        }
    }
  if(found==0)
    {//last chance ! Is it a support mesh ?
      med_int const nbSM(MEDnSupportMesh(fid));
      for(int i=0;i<nbSM && found==0;i++)
        {
          med_int const naxis(MEDsupportMeshnAxis(fid,i+1));
          INTERP_KERNEL::AutoPtr<char> axisname(MEDLoaderBase::buildEmptyString(naxis*MED_SNAME_SIZE)),axisunit(MEDLoaderBase::buildEmptyString(naxis*MED_SNAME_SIZE));
          MEDFILESAFECALLERRD0(MEDsupportMeshInfo,(fid,i+1,nommaa,&spaceDim,&dim,maillage_description,&axistype,axisname,axisunit));
          std::string const cur(MEDLoaderBase::buildStringFromFortran(nommaa,sizeof(nommaa)));
          ms.push_back(cur);
          if(cur==mName)
            {
              found=2;
              ret=i+1;
            }
        }
    }
  ////////////////////////
  switch(found)
    {
    case 1:
      {
        axType=TraduceAxisType(axistype);
        switch(type_maillage)
          {
          case MED_UNSTRUCTURED_MESH:
            meshType=UNSTRUCTURED;
            break;
          case MED_STRUCTURED_MESH:
            {
              med_grid_type gt;
              MEDFILESAFECALLERRD0(MEDmeshGridTypeRd,(fid,mName.c_str(),&gt));
              switch(gt)
                {
                case MED_CARTESIAN_GRID:
                  meshType=CARTESIAN;
                  break;
                case MED_CURVILINEAR_GRID:
                  meshType=CURVE_LINEAR;
                  break;
                case MED_POLAR_GRID:// this is not a bug. A MED file POLAR_GRID is deal by CARTESIAN MEDLoader
                  meshType=CARTESIAN;
                  break;
                default:
                  throw INTERP_KERNEL::Exception("MEDFileMeshL2::getMeshIdFromName : unrecognized structured mesh type ! Supported are :\n - cartesian\n - curve linear\n");
                }
              break;
            }
          default:
            throw INTERP_KERNEL::Exception("MEDFileMeshL2::getMeshIdFromName : unrecognized mesh type !");
          }
        med_int numdt,numit;
        med_float dtt;
        MEDFILESAFECALLERRD0(MEDmeshComputationStepInfo,(fid,mName.c_str(),1,&numdt,&numit,&dtt));
        dt=FromMedInt<int>(numdt); it=FromMedInt<int>(numit);
        return new MeshCls(ret);
      }
    case 2:
      {
        meshType=UNSTRUCTURED;
        dt=MED_NO_DT; it=MED_NO_IT; dtunit1.clear();
        axType=TraduceAxisType(axistype);
        return new StructMeshCls(ret);
      }
    default:
      {
        std::ostringstream oss;
        oss << "No such meshname (" << mName <<  ") in file ! Must be in : ";
        std::copy(ms.begin(),ms.end(),std::ostream_iterator<std::string>(oss,", "));
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    }
  
}

/*!
 * non static and non const method because _description, _dt_unit... are set in this method.
 */
std::vector<std::string> MEDFileMeshL2::getAxisInfoOnMesh(med_idt fid, const MeshOrStructMeshCls *mId, const std::string& mName, MEDCoupling::MEDCouplingMeshType& meshType, MEDCoupling::MEDCouplingAxisType& axType, int& nstep, int& Mdim)
{
  return mId->getAxisInfoOnMesh(fid,mName,meshType,axType,nstep,Mdim,_description,_dt_unit,_univ_name);
}

void MEDFileMeshL2::ReadFamiliesAndGrps(med_idt fid, const std::string& meshName, std::map<std::string,mcIdType>& fams, std::map<std::string, std::vector<std::string> >& grps, MEDFileMeshReadSelector *mrs)
{
  if(mrs && !(mrs->isCellFamilyFieldReading() || mrs->isNodeFamilyFieldReading()))
    return ;
  char nomfam[MED_NAME_SIZE+1];
  med_int numfam;
  med_int const nfam=MEDnFamily(fid,meshName.c_str());
  std::vector< std::pair<std::string,std::pair<mcIdType,std::vector<std::string> > > > crudeFams(nfam);
  for(int i=0;i<nfam;i++)
    {
      med_int const ngro=MEDnFamilyGroup(fid,meshName.c_str(),i+1);
      med_int const natt=MEDnFamily23Attribute(fid,meshName.c_str(),i+1);
      INTERP_KERNEL::AutoPtr<med_int> attide=new med_int[natt];
      INTERP_KERNEL::AutoPtr<med_int> attval=new med_int[natt];
      INTERP_KERNEL::AutoPtr<char> attdes=new char[MED_COMMENT_SIZE*natt+1];
      INTERP_KERNEL::AutoPtr<char> gro=new char[MED_LNAME_SIZE*ngro+1];
      MEDfamily23Info(fid,meshName.c_str(),i+1,nomfam,attide,attval,attdes,&numfam,gro);
      std::string const famName(MEDLoaderBase::buildStringFromFortran(nomfam,MED_NAME_SIZE));
      std::vector<std::string> grps2(ngro);
      for(int j=0;j<ngro;j++)
        grps2[j]=MEDLoaderBase::buildStringFromFortran(gro+j*MED_LNAME_SIZE,MED_LNAME_SIZE);
      crudeFams[i]=std::pair<std::string,std::pair<mcIdType,std::vector<std::string> > >(famName,std::pair<mcIdType,std::vector<std::string> >(numfam,grps2));
    }
  RenameFamiliesFromFileToMemInternal(crudeFams);
  for(const auto & crudeFam : crudeFams)
    {
      fams[crudeFam.first]=crudeFam.second.first;
      for(auto it1=crudeFam.second.second.begin();it1!=crudeFam.second.second.end();it1++)
        grps[*it1].push_back(crudeFam.first);
    }
}

void MEDFileMeshL2::WriteFamiliesAndGrps(med_idt fid, const std::string& mname, const std::map<std::string,mcIdType>& fams, const std::map<std::string, std::vector<std::string> >& grps, int tooLongStrPol)
{
  std::vector< std::pair<std::string,std::pair<mcIdType,std::vector<std::string> > > > crudeFams(fams.size());
  std::size_t ii(0);
  for(auto it=fams.begin();it!=fams.end();it++,ii++)
    {
      std::vector<std::string> grpsOfFam;
      for(const auto & grp : grps)
        {
          if(std::find(grp.second.begin(),grp.second.end(),(*it).first)!=grp.second.end())
            grpsOfFam.push_back(grp.first);
        }
      crudeFams[ii]=std::pair<std::string,std::pair<mcIdType,std::vector<std::string> > >((*it).first,std::pair<mcIdType,std::vector<std::string> >((*it).second,grpsOfFam));
    }
  RenameFamiliesFromMemToFileInternal(crudeFams);
  for(const auto & crudeFam : crudeFams)
    {
      std::size_t const ngro(crudeFam.second.second.size());
      INTERP_KERNEL::AutoPtr<char> groName=MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE*ngro);
      int i=0;
      for(auto it2=crudeFam.second.second.begin();it2!=crudeFam.second.second.end();it2++,i++)
        MEDLoaderBase::safeStrCpy2((*it2).c_str(),MED_LNAME_SIZE,groName+i*MED_LNAME_SIZE,tooLongStrPol);
      INTERP_KERNEL::AutoPtr<char> famName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
      MEDLoaderBase::safeStrCpy(crudeFam.first.c_str(),MED_NAME_SIZE,famName,tooLongStrPol);
      med_int ret=MEDfamilyCr(fid,mname.c_str(),famName,ToMedInt(crudeFam.second.first),ToMedInt(ngro),groName);
      ret++;
    }
}

void MEDFileMeshL2::RenameFamiliesPatternInternal(std::vector< std::pair<std::string,std::pair<mcIdType,std::vector<std::string> > > >& crudeFams, RenameFamiliesPatternFunc func)
{
  std::size_t ii(0);
  std::vector<std::string> fams(crudeFams.size());
  for(std::vector< std::pair<std::string,std::pair<mcIdType,std::vector<std::string> > > >::const_iterator it=crudeFams.begin();it!=crudeFams.end();it++,ii++)
    fams[ii]=(*it).first;
  if(!func(fams))
    return ;
  ii=0;
  for(auto it=crudeFams.begin();it!=crudeFams.end();it++,ii++)
    (*it).first=fams[ii];
}

/*!
 * This method is dedicated to the killers that use a same family name to store different family ids. MED file API authorizes it.
 * So this method renames families (if needed generally not !) in order to have a discriminant name for families.
 */
void MEDFileMeshL2::RenameFamiliesFromFileToMemInternal(std::vector< std::pair<std::string,std::pair<mcIdType,std::vector<std::string> > > >& crudeFams)
{
  RenameFamiliesPatternInternal(crudeFams,RenameFamiliesFromFileToMem);
}

bool MEDFileMeshL2::RenameFamiliesFromFileToMem(std::vector< std::string >& famNames)
{
  std::map<std::string,mcIdType> m;
  std::set<std::string> s;
  for(const auto & famName : famNames)
    {
      if(s.find(famName)!=s.end())
        m[famName]=0;
      s.insert(famName);
    }
  if(m.empty())
    return false;// the general case !
  for(auto & famName : famNames)
    {
      auto const it2(m.find(famName));
      if(it2!=m.end())
        {
          std::ostringstream oss; oss << famName << ZE_SEP_FOR_FAMILY_KILLERS << std::setfill('0') << std::setw(ZE_SEP2_FOR_FAMILY_KILLERS) << (*it2).second++;
          famName=oss.str();
        }
    }
  return true;
}

/*!
 * This method is dedicated to the killers that use a same family name to store different family ids. MED file API authorizes it.
 * So this method renames families (if needed generally not !) in order to have a discriminant name for families.
 */
void MEDFileMeshL2::RenameFamiliesFromMemToFileInternal(std::vector< std::pair<std::string,std::pair<mcIdType,std::vector<std::string> > > >& crudeFams)
{
  RenameFamiliesPatternInternal(crudeFams,RenameFamiliesFromMemToFile);
}

bool MEDFileMeshL2::RenameFamiliesFromMemToFile(std::vector< std::string >& famNames)
{
  bool isSmthingStrange(false);
  for(const auto & famName : famNames)
    {
      std::size_t const found(famName.find(ZE_SEP_FOR_FAMILY_KILLERS));
      if(found!=std::string::npos)
        isSmthingStrange=true;
    }
  if(!isSmthingStrange)
    return false;
  // pattern matching
  std::map< std::string, std::vector<std::string> > m;
  for(const auto & famName : famNames)
    {
      std::size_t const found(famName.find(ZE_SEP_FOR_FAMILY_KILLERS));
      if(found!=std::string::npos && found>=1)
        {
          std::string const s1(famName.substr(found+sizeof(ZE_SEP_FOR_FAMILY_KILLERS)-1));
          if((int)s1.size()!=ZE_SEP2_FOR_FAMILY_KILLERS)
            continue;
          int k(-1);
          std::istringstream iss(s1);
          iss >> k;
          bool const isOK((iss.rdstate() & ( std::istream::failbit | std::istream::eofbit)) == std::istream::eofbit);
          if(isOK && k>=0)
            {
              std::string const s0(famName.substr(0,found));
              m[s0].push_back(famName);
            }
        }
    }
  if(m.empty())
    return false;
  // filtering
  std::map<std::string,std::string> zeMap;
  for(std::map< std::string, std::vector<std::string> >::const_iterator it=m.begin();it!=m.end();it++)
    {
      if((*it).second.size()==1)
        continue;
      for(auto it1=(*it).second.begin();it1!=(*it).second.end();it1++)
        zeMap[*it1]=(*it).first;
    }
  if(zeMap.empty())
    return false;
  // traduce
  for(auto & famName : famNames)
    {
      auto const it1(zeMap.find(famName));
      if(it1!=zeMap.end())
        famName=(*it1).second;
    }    
  return true;
}

MEDCoupling::MEDCouplingAxisType MEDFileMeshL2::TraduceAxisType(med_axis_type at)
{
  switch(at)
    {
    case MED_CARTESIAN:
      return AX_CART;
    case MED_CYLINDRICAL:
      return AX_CYL;
    case MED_SPHERICAL:
      return AX_SPHER;
    case MED_UNDEF_AXIS_TYPE:
      return AX_CART;
    default:
      throw INTERP_KERNEL::Exception("MEDFileMeshL2::TraduceAxisType : unrecognized axis type !");
    }
}

MEDCoupling::MEDCouplingAxisType MEDFileMeshL2::TraduceAxisTypeStruct(med_grid_type gt)
{
  switch(gt)
    {
    case MED_CARTESIAN_GRID:
      return AX_CART;
    case MED_POLAR_GRID:
      return AX_CYL;
    default:
      throw INTERP_KERNEL::Exception("MEDFileMeshL2::TraduceAxisTypeStruct : only Cartesian and Cylindrical supported by MED file !");
    }
}

med_axis_type MEDFileMeshL2::TraduceAxisTypeRev(MEDCoupling::MEDCouplingAxisType at)
{
  switch(at)
    {
    case AX_CART:
      return MED_CARTESIAN;
    case AX_CYL:
      return MED_CYLINDRICAL;
    case AX_SPHER:
      return MED_SPHERICAL;
    default:
      throw INTERP_KERNEL::Exception("MEDFileMeshL2::TraduceAxisTypeRev : unrecognized axis type !");
    }
}

med_grid_type MEDFileMeshL2::TraduceAxisTypeRevStruct(MEDCoupling::MEDCouplingAxisType at)
{
  switch(at)
    {
    case AX_CART:
      return MED_CARTESIAN_GRID;
    case AX_CYL:
      return MED_POLAR_GRID;
    case AX_SPHER:
      return MED_POLAR_GRID;
    default:
      throw INTERP_KERNEL::Exception("MEDFileMeshL2::TraduceAxisTypeRevStruct : unrecognized axis type !");
    }
}

MEDFileUMeshL2::MEDFileUMeshL2()
= default;

std::vector<std::string> MEDFileUMeshL2::loadCommonPart(med_idt fid, const MeshOrStructMeshCls *mId, const std::string& mName, int dt, int it, int& Mdim)
{
  Mdim=-3;
  _name.set(mName.c_str());
  int nstep;
  MEDCoupling::MEDCouplingMeshType meshType;
  MEDCoupling::MEDCouplingAxisType dummy3;
  std::vector<std::string> ret(getAxisInfoOnMesh(fid,mId,mName.c_str(),meshType,dummy3,nstep,Mdim));
  if(nstep==0)
    {
      Mdim=-4;
      return std::vector<std::string>();
    }
  if(meshType!=UNSTRUCTURED)
    throw INTERP_KERNEL::Exception("Invalid mesh type ! You are expected an unstructured one whereas in file it is not an unstructured !");
  _time=mId->checkMeshTimeStep(fid,mName,nstep,dt,it);
  _iteration=dt;
  _order=it;
  return ret;
}

void MEDFileUMeshL2::loadAll(med_idt fid, const MeshOrStructMeshCls *mId, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs)
{
  int Mdim;
  std::vector<std::string> const infosOnComp(loadCommonPart(fid,mId,mName,dt,it,Mdim));
  if(Mdim==-4)
    return ;
  loadConnectivity(fid,Mdim,mName,dt,it,mrs);//to improve check (dt,it) coherency
  loadCoords(fid,infosOnComp,mName,dt,it);
}

/*!
 * This method is expected to be invoked after the load of connectivity.
 * This method is in charge of :
 *  - dealing with optimized load of coordinates (loading only points fetched by the already loaded cells)
 *  - update the connectivity in \a this to fit the coordinates loaded just above
 */
void MEDFileUMeshL2::dealWithCoordsInLoadPart(med_idt fid, const MeshOrStructMeshCls * /*mId*/, const std::string& mName, const std::vector<std::string>& infosOnComp, const std::vector<INTERP_KERNEL::NormalizedCellType>&  /*types*/, const std::vector<mcIdType>&  /*slicPerTyp*/, int dt, int it, MEDFileMeshReadSelector *mrs)
{
  med_bool changement,transformation;
  mcIdType const nCoords(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NONE,MED_COORDINATE,MED_NO_CMODE,&changement,&transformation));
  std::vector<bool> fetchedNodeIds(nCoords,false);
  for(const auto & it0 : _per_type_mesh)
    for(const auto & it1 : it0)
      it1->getMesh()->computeNodeIdsAlg(fetchedNodeIds);
  if(!mrs || mrs->getNumberOfCoordsLoadSessions()==1)
  {
    mcIdType const nMin(ToIdType(std::distance(fetchedNodeIds.begin(),std::find(fetchedNodeIds.begin(),fetchedNodeIds.end(),true))));
    mcIdType nMax(ToIdType(std::distance(fetchedNodeIds.rbegin(),std::find(fetchedNodeIds.rbegin(),fetchedNodeIds.rend(),true))));
    nMax=nCoords-nMax;
    for(const auto & it0 : _per_type_mesh)
      for(const auto & it1 : it0)
        it1->getMesh()->renumberNodesWithOffsetInConn(-nMin);
    this->loadPartCoords(fid,infosOnComp,mName,dt,it,nMin,nMax);
  }
  else
  {
    mcIdType const nbOfCooLS(mrs->getNumberOfCoordsLoadSessions());
    MCAuto<DataArrayIdType> fni(DataArrayIdType::BuildListOfSwitchedOn(fetchedNodeIds));
    MCAuto< MapKeyVal<mcIdType, mcIdType> > o2n(fni->invertArrayN2O2O2NOptimized());
    for(const auto & it0 : _per_type_mesh)
      for(const auto & it1 : it0)
        it1->getMesh()->renumberNodesInConn(o2n->data());
    this->loadPartCoordsSlice(fid,infosOnComp,mName,dt,it,fni,nbOfCooLS);
  }
}

std::vector<std::string> MEDFileUMeshL2::loadPartConnectivityOnly(med_idt fid, const MeshOrStructMeshCls *mId, const std::string& mName, const std::vector<INTERP_KERNEL::NormalizedCellType>& types, const std::vector<mcIdType>& slicPerTyp, int dt, int it, MEDFileMeshReadSelector *mrs, int& Mdim)
{
  std::vector<std::string> infosOnComp(loadCommonPart(fid,mId,mName,dt,it,Mdim));
  if(Mdim==-4)
    return infosOnComp;
  loadPartOfConnectivity(fid,Mdim,mName,types,slicPerTyp,dt,it,mrs);
  return infosOnComp;
}

void MEDFileUMeshL2::loadPart(med_idt fid, const MeshOrStructMeshCls *mId, const std::string& mName, const std::vector<INTERP_KERNEL::NormalizedCellType>& types, const std::vector<mcIdType>& slicPerTyp, int dt, int it, MEDFileMeshReadSelector *mrs)
{
  int Mdim;
  std::vector<std::string> const infosOnComp(loadPartConnectivityOnly(fid,mId,mName,types,slicPerTyp,dt,it,mrs,Mdim));
  if(Mdim==-4)
    return ;
  loadPartOfConnectivity(fid,Mdim,mName,types,slicPerTyp,dt,it,mrs);
  dealWithCoordsInLoadPart(fid,mId,mName,infosOnComp,types,slicPerTyp,dt,it,mrs);
}

/*!
 * This method loads from file \a fid a part of the mesh (made of same geometrical type cells \a type) called \a mName. The loading is done in 2 steps:
 * First, we load the connectivity of nodes.
 * Second, we load coordinates of nodes lying in the specified cells (same as MEDFileUMeshL2::dealWithCoordsInLoadPart, except in this case, we're not limited to slice of nodes)
 * \throw exception if multiple load sessions are requested
 */
void MEDFileUMeshL2::loadPartFromUserDistrib(med_idt fid, const MeshOrStructMeshCls *mId, const std::string& mName, const std::map<INTERP_KERNEL::NormalizedCellType,std::vector<mcIdType>>& distrib, int dt, int it, MEDFileMeshReadSelector *mrs)
{
  int Mdim;
  std::vector<std::string> const infosOnComp(loadCommonPart(fid,mId,mName,dt,it,Mdim));
  if(Mdim==-4)
    return ;

  /* First step : loading connectivity of nodes, ie building a new mesh of one geometrical type with only the specified cells in distrib */
  loadPartOfConnectivityFromUserDistrib(fid,Mdim,mName,distrib,dt,it,mrs);

  /* Second step : loading nodes */
  med_bool changement,transformation;
  mcIdType const nCoords(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NONE,MED_COORDINATE,MED_NO_CMODE,&changement,&transformation));
  std::vector<bool> fetchedNodeIds(nCoords,false);
  for(const auto & it0 : _per_type_mesh)
    for(const auto & it1 : it0)
      it1->getMesh()->computeNodeIdsAlg(fetchedNodeIds);    // for each node in the original mesh, which ones are laying on the current single geometrical type partial mesh

  if(!mrs || mrs->getNumberOfCoordsLoadSessions()==1)
    {
      // renumbering nodes inside the connectivity of the partial mesh:
      // until now, the numbering used in the connectivity of the cells was that of the integral mesh,
      // so it might be sparsed as some original nodes are missing in the partial mesh,
      // thus we want each node to be renumbered so that the sequence of their numbers form a range
      MCAuto<DataArrayIdType> fni(DataArrayIdType::BuildListOfSwitchedOn(fetchedNodeIds));
      MCAuto< MapKeyVal<mcIdType, mcIdType> > o2n(fni->invertArrayN2O2O2NOptimized());
      for(const auto & it0 : _per_type_mesh)
        for(const auto & it1 : it0)
          it1->getMesh()->renumberNodesInConn(o2n->data());

      // loading coordinates of fetched nodes
      std::vector<mcIdType> distribNodes;
      for(std::map<mcIdType,mcIdType>::const_iterator mapIter = o2n->data().begin(); mapIter != o2n->data().end(); ++mapIter)
          distribNodes.push_back(mapIter->first);
      this->loadPartCoords(fid,infosOnComp,mName,dt,it,distribNodes);
    }
  else
    throw INTERP_KERNEL::Exception("MEDFileUMeshL2::loadPartFromUserDistrib: multiple load sessions not handled!");
}

void MEDFileUMeshL2::loadConnectivity(med_idt fid, int mdim, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs)
{
  _per_type_mesh.resize(1);
  _per_type_mesh[0].clear();
  for(int j=0;j<MED_N_CELL_FIXED_GEO;j++)
    {
      MEDFileUMeshPerType *tmp(MEDFileUMeshPerType::New(fid,mName.c_str(),dt,it,mdim,typmai[j],typmai2[j],mrs));
      if(tmp)
        _per_type_mesh[0].push_back(tmp);
    }
  sortTypes();
}

void MEDFileUMeshL2::loadPartOfConnectivity(med_idt fid, int mdim, const std::string& mName, const std::vector<INTERP_KERNEL::NormalizedCellType>& types, const std::vector<mcIdType>& slicPerTyp, int dt, int it, MEDFileMeshReadSelector *mrs)
{
  std::size_t const nbOfTypes(types.size());
  if(slicPerTyp.size()!=3*nbOfTypes)
    throw INTERP_KERNEL::Exception("MEDFileUMeshL2::loadPartOfConnectivity : The size of slicPerTyp array is expected to be equal to 3 times size of array types !");
  std::set<INTERP_KERNEL::NormalizedCellType> const types2(types.begin(),types.end());
  if(types2.size()!=nbOfTypes)
    throw INTERP_KERNEL::Exception("MEDFileUMeshL2::loadPartOfConnectivity : the geometric types in types array must appear once !");
  _per_type_mesh.resize(1);
  _per_type_mesh[0].clear();
  for(std::size_t ii=0;ii<nbOfTypes;ii++)
    {
      mcIdType strt(slicPerTyp[3*ii+0]),stp(slicPerTyp[3*ii+1]),step(slicPerTyp[3*ii+2]);
      MCAuto<MEDFileUMeshPerType> const tmp(MEDFileUMeshPerType::NewPart(fid,mName.c_str(),dt,it,mdim,types[ii],strt,stp,step,mrs));
      _per_type_mesh[0].push_back(tmp);
    }
  sortTypes();
}

/*!
 * This method builds a new mesh of single geometrical type based on the partition of cells \a distrib, from mesh \a mName in file \a fid.
 * This distribution is not necessarily a slice.
 */
void MEDFileUMeshL2::loadPartOfConnectivityFromUserDistrib(med_idt fid, int mdim, const std::string& mName, const std::map<INTERP_KERNEL::NormalizedCellType,std::vector<mcIdType>>& distrib, int dt, int it, MEDFileMeshReadSelector *mrs)
{
  _per_type_mesh.resize(1);
  _per_type_mesh[0].clear();
  std::map<INTERP_KERNEL::NormalizedCellType,std::vector<mcIdType>>::const_iterator iter;
  for (iter = distrib.begin(); iter != distrib.end(); iter++)
    {
        MCAuto<MEDFileUMeshPerType> const tmp(MEDFileUMeshPerType::NewPart(fid,mName.c_str(),dt,it,mdim,iter->first/*type*/,iter->second/*distrib over the current type*/,mrs));
        _per_type_mesh[0].push_back(tmp);
    }
  sortTypes();
}

void MEDFileUMeshL2::loadCoords(med_idt fid, const std::vector<std::string>& infosOnComp, const std::string& mName, int dt, int it)
{
  int const spaceDim((int)infosOnComp.size());
  med_bool changement,transformation;
  med_int const nCoords(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NONE,MED_COORDINATE,MED_NO_CMODE,&changement,&transformation));
  _coords=DataArrayDouble::New();
  _coords->alloc(nCoords,spaceDim);
  double *coordsPtr(_coords->getPointer());
  if (nCoords)
    MEDFILESAFECALLERRD0(MEDmeshNodeCoordinateRd,(fid,mName.c_str(),dt,it,MED_FULL_INTERLACE,coordsPtr));
  if(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,MED_FAMILY_NUMBER,MED_NODAL,&changement,&transformation)>0)
    {
      MCAuto<DataArrayMedInt> miFamCoord=DataArrayMedInt::New();
      miFamCoord->alloc(nCoords,1);
      MEDFILESAFECALLERRD0(MEDmeshEntityFamilyNumberRd,(fid,mName.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,miFamCoord->getPointer()));
      _fam_coords=FromMedIntArray<mcIdType>(miFamCoord);
    }
  else
    _fam_coords=nullptr;
  if(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,MED_NUMBER,MED_NODAL,&changement,&transformation)>0)
    {
      MCAuto<DataArrayMedInt> miNumCoord=DataArrayMedInt::New();
      miNumCoord->alloc(nCoords,1);
      MEDFILESAFECALLERRD0(MEDmeshEntityNumberRd,(fid,mName.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,miNumCoord->getPointer()));
      _num_coords=FromMedIntArray<mcIdType>(miNumCoord);
    }
  else
    _num_coords=nullptr;
  if(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,MED_NAME,MED_NODAL,&changement,&transformation)>0)
    {
      _name_coords=DataArrayAsciiChar::New();
      _name_coords->alloc(nCoords+1,MED_SNAME_SIZE);//not a bug to avoid the memory corruption due to last \0 at the end
      MEDFILESAFECALLERRD0(MEDmeshEntityNameRd,(fid,mName.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,_name_coords->getPointer()));
      _name_coords->reAlloc(nCoords);//not a bug to avoid the memory corruption due to last \0 at the end
    }
  else
    _name_coords=nullptr;
  if(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,MED_GLOBAL_NUMBER,MED_NODAL,&changement,&transformation)>0)
    {
      MCAuto<DataArrayMedInt> miNumCoord=DataArrayMedInt::New();
      miNumCoord->alloc(nCoords,1);
      MEDFILESAFECALLERRD0(MEDmeshGlobalNumberRd,(fid,mName.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,miNumCoord->getPointer()));
      _global_num_coords=FromMedIntArray<mcIdType>(miNumCoord);
    }
  for(int i=0;i<spaceDim;i++)
    _coords->setInfoOnComponent(i,infosOnComp[i]);
}


void MEDFileUMeshL2::LoadPartCoords(med_idt fid, const std::vector<std::string>& infosOnComp, const std::string& mName, int dt, int it, const std::vector<mcIdType>& distribNodes,
MCAuto<DataArrayDouble>& _coords, MCAuto<PartDefinition>& _part_coords, MCAuto<DataArrayIdType>& _fam_coords, MCAuto<DataArrayIdType>& _num_coords, MCAuto<DataArrayAsciiChar>& _name_coords)
{
  med_int const spaceDim((int)infosOnComp.size());
  allocCoordsPartCoords(spaceDim,distribNodes,_coords,_part_coords);
  _coords->setInfoOnComponents(infosOnComp);
  fillPartCoords(fid,spaceDim,mName,dt,it,_part_coords,_coords,_fam_coords,_num_coords,_name_coords);
}

void MEDFileUMeshL2::LoadPartCoords(med_idt fid, const std::vector<std::string>& infosOnComp, const std::string& mName, int dt, int it, mcIdType nMin, mcIdType nMax,
MCAuto<DataArrayDouble>& _coords, MCAuto<PartDefinition>& _part_coords, MCAuto<DataArrayIdType>& _fam_coords, MCAuto<DataArrayIdType>& _num_coords, MCAuto<DataArrayAsciiChar>& _name_coords)
{
  med_int const spaceDim((int)infosOnComp.size());
  allocCoordsPartCoords(spaceDim,nMin,nMax,_coords,_part_coords);
  _coords->setInfoOnComponents(infosOnComp);
  fillPartCoords(fid,spaceDim,mName,dt,it,_part_coords,_coords,_fam_coords,_num_coords,_name_coords);
}

/*!
 * This method allocates the space needed to load coordinates of nodes specified in the vector \a nodeIds and creates a PartDefinition object to store the ids in \a nodeIds
 */
void MEDFileUMeshL2::allocCoordsPartCoords(mcIdType spaceDim, const std::vector<mcIdType>& nodeIds, MCAuto<DataArrayDouble>& _coords, MCAuto<PartDefinition>& _part_coords)
{
  mcIdType const nbNodesToLoad(nodeIds.size());
  _coords=DataArrayDouble::New();
  _coords->alloc(nbNodesToLoad,spaceDim);

  MCAuto<DataArrayIdType> nodeIdsArray=DataArrayIdType::New();
  nodeIdsArray->useArray(nodeIds.data(),false,DeallocType::C_DEALLOC,nbNodesToLoad,1);
  _part_coords=PartDefinition::New(nodeIdsArray);
}

/*!
 * This method allocates the space needed to load coordinates of all nodes between \a nMin and \a nMax and creates a PartDefinition object to store them
 */
void MEDFileUMeshL2::allocCoordsPartCoords(mcIdType spaceDim, mcIdType nMin, mcIdType nMax, MCAuto<DataArrayDouble>& _coords, MCAuto<PartDefinition>& _part_coords)
{
  _coords=DataArrayDouble::New();
  mcIdType const nbNodesToLoad(nMax-nMin);
  _coords->alloc(nbNodesToLoad,spaceDim);

  _part_coords=PartDefinition::New(nMin,nMax,1);
}

/*!
 * This method loads coordinates of every node in \a partCoords and additionnal low-level information
 */
void MEDFileUMeshL2::fillPartCoords(med_idt fid, mcIdType spaceDim, const std::string& mName, int dt, int it, const PartDefinition *partCoords,
                                    MCAuto<DataArrayDouble>& _coords, MCAuto<DataArrayIdType>& _fam_coords, MCAuto<DataArrayIdType>& _num_coords, MCAuto<DataArrayAsciiChar>& _name_coords)
{
  med_bool changement,transformation;
  med_int const nCoords(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NONE,MED_COORDINATE,MED_NO_CMODE,&changement,&transformation));
  mcIdType const nbNodesToLoad = partCoords->getNumberOfElems();

  // Based on the ids in \a partCoords, defining the appropriate med_filter (filter of block if the ids form a slice, a generic filter otherwise)
  {
    MEDFilterEntity filter1;
    filter1.fill(fid,/*nentity*/nCoords,/*nvaluesperentity*/1,/*nconstituentpervalue*/spaceDim,
        MED_ALL_CONSTITUENT,MED_FULL_INTERLACE,MED_COMPACT_STMODE,MED_NO_PROFILE,
        partCoords);
    // With the filter defined above, retrieve coordinates of nodes
    MEDFILESAFECALLERRD0(MEDmeshNodeCoordinateAdvancedRd,(fid,mName.c_str(),dt,it,filter1.getPtr(),_coords->getPointer()));
  }

  {
    MEDFilterEntity filter2;
    filter2.fill(fid,nCoords,1,1,
        MED_ALL_CONSTITUENT,MED_FULL_INTERLACE,MED_COMPACT_STMODE,MED_NO_PROFILE, partCoords);

    // Retrieve additional information regarding nodes
    if(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,MED_FAMILY_NUMBER,MED_NODAL,&changement,&transformation)>0)
      {
        MCAuto<DataArrayMedInt> miFamCoord=DataArrayMedInt::New();
        miFamCoord->alloc(nbNodesToLoad,1);
        MEDFILESAFECALLERRD0(MEDmeshEntityAttributeAdvancedRd,(fid,mName.c_str(),MED_FAMILY_NUMBER,dt,it,MED_NODE,MED_NO_GEOTYPE,filter2.getPtr(),miFamCoord->getPointer()));
        _fam_coords=FromMedIntArray<mcIdType>(miFamCoord);
      }
    else
      _fam_coords=nullptr;
    if(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,MED_NUMBER,MED_NODAL,&changement,&transformation)>0)
      {
        MCAuto<DataArrayMedInt> miNumCoord=DataArrayMedInt::New();
        miNumCoord->alloc(nbNodesToLoad,1);
        MEDFILESAFECALLERRD0(MEDmeshEntityAttributeAdvancedRd,(fid,mName.c_str(),MED_NUMBER,dt,it,MED_NODE,MED_NO_GEOTYPE,filter2.getPtr(),miNumCoord->getPointer()));
        _num_coords=FromMedIntArray<mcIdType>(miNumCoord);
      }
    else
      _num_coords=nullptr;
    if(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,MED_NAME,MED_NODAL,&changement,&transformation)>0)
      {
        _name_coords=DataArrayAsciiChar::New();
        _name_coords->alloc(nbNodesToLoad+1,MED_SNAME_SIZE);//not a bug to avoid the memory corruption due to last \0 at the end
        MEDFILESAFECALLERRD0(MEDmeshEntityAttributeAdvancedRd,(fid,mName.c_str(),MED_NAME,dt,it,MED_NODE,MED_NO_GEOTYPE,filter2.getPtr(),_name_coords->getPointer()));
        _name_coords->reAlloc(nbNodesToLoad);//not a bug to avoid the memory corruption due to last \0 at the end
      }
    else
      _name_coords=nullptr;
  }  // filter2
}


/*!
 * For performance reasons LoadPartCoordsArray method calls LoadPartCoords
 */
void MEDFileUMeshL2::LoadPartCoordsArray(med_idt fid, const std::vector<std::string>& infosOnComp, const std::string& mName, int dt, int it, const DataArrayIdType *nodeIds,
MCAuto<DataArrayDouble>& _coords, MCAuto<DataArrayIdType>& _fam_coords, MCAuto<DataArrayIdType>& _num_coords, MCAuto<DataArrayAsciiChar>& _name_coords)
{
  MCAuto<PartDefinition> useless;
  nodeIds->checkAllocated();
  nodeIds->checkNbOfComps(1,"loadPartCoordsSlice : Only one component expected !");
  mcIdType nMin(0),nMax(0);
  if(!nodeIds->empty())
  { nMin = nodeIds->front(); nMax = nodeIds->back()+1; }
  LoadPartCoords(fid,infosOnComp,mName,dt,it,nMin,nMax,_coords,useless,_fam_coords,_num_coords,_name_coords);
  if(nodeIds->empty())
    return ;
  MCAuto<DataArrayIdType> nodeIds2(nodeIds->deepCopy());
  nodeIds2->applyLin(1,-nMin);
  _coords = _coords->selectByTupleIdSafe(nodeIds2->begin(),nodeIds2->end());
  if(_fam_coords.isNotNull())
    _fam_coords = _fam_coords->selectByTupleIdSafe(nodeIds2->begin(),nodeIds2->end());
  if(_num_coords.isNotNull())
    _num_coords = _num_coords->selectByTupleIdSafe(nodeIds2->begin(),nodeIds2->end());
  if(_name_coords.isNotNull())
    {
      MCAuto<DataArrayChar> tmp(_name_coords->selectByTupleIdSafe(nodeIds2->begin(),nodeIds2->end()));
      _name_coords = DynamicCastSafe<DataArrayChar,DataArrayAsciiChar>( tmp );
    }
}

void MEDFileUMeshL2::loadPartCoords(med_idt fid, const std::vector<std::string>& infosOnComp, const std::string& mName, int dt, int it, mcIdType nMin, mcIdType nMax)
{
  LoadPartCoords(fid,infosOnComp,mName,dt,it,nMin,nMax,_coords,_part_coords,_fam_coords,_num_coords,_name_coords);
}

void MEDFileUMeshL2::loadPartCoords(med_idt fid, const std::vector<std::string>& infosOnComp, const std::string& mName, int dt, int it, const std::vector<mcIdType>& distribNodes)
{
  LoadPartCoords(fid,infosOnComp,mName,dt,it,distribNodes,_coords,_part_coords,_fam_coords,_num_coords,_name_coords);
}


void MEDFileUMeshL2::loadPartCoordsSlice(med_idt fid, const std::vector<std::string>& infosOnComp, const std::string& mName, int dt, int it, const DataArrayIdType *nodeIds, mcIdType nbOfCoordLS)
{
  nodeIds->checkAllocated();
  nodeIds->checkNbOfComps(1,"loadPartCoordsSlice : Only one component expected !");
  if(nodeIds->empty())
    return ;
  if( nbOfCoordLS<1 )
    throw INTERP_KERNEL::Exception("MEDFileUMeshL2::loadPartCoordsSlice : nb of coords load session must be >=1 !");
  mcIdType nMin(nodeIds->front()),nMax(nodeIds->back()+1);
  std::vector< MCAuto<DataArrayDouble> > coords(nbOfCoordLS);
  std::vector< MCAuto<DataArrayIdType> > famCoords(nbOfCoordLS);
  std::vector< MCAuto<DataArrayIdType> > numCoords(nbOfCoordLS);
  std::vector< MCAuto<DataArrayAsciiChar> > nameCoords(nbOfCoordLS);
  for(mcIdType ipart = 0 ; ipart < nbOfCoordLS ; ++ipart)
    {
      mcIdType partStart,partStop;
      DataArray::GetSlice(nMin,nMax,1,ipart,nbOfCoordLS,partStart,partStop);
      MCAuto<DataArrayIdType> idsNodeIdsToKeep(nodeIds->findIdsInRange(partStart,partStop));
      MCAuto<DataArrayIdType> const nodeIdsToKeep( nodeIds->selectByTupleIdSafe(idsNodeIdsToKeep->begin(),idsNodeIdsToKeep->end()) );
      LoadPartCoordsArray(fid,infosOnComp,mName,dt,it,nodeIdsToKeep,coords[ipart],famCoords[ipart],numCoords[ipart],nameCoords[ipart]);
    }
  _coords = DataArrayDouble::Aggregate(ToConstVect<DataArrayDouble>(coords));
  if(famCoords[0].isNotNull())
    _fam_coords = DataArrayIdType::Aggregate(ToConstVect<DataArrayIdType>(famCoords));
  if(numCoords[0].isNotNull())
    _num_coords = DataArrayIdType::Aggregate(ToConstVect<DataArrayIdType>(numCoords));
  if(nameCoords[0].isNotNull())
  {
    std::vector< MCAuto<DataArrayChar> > nameCoords2(nameCoords.begin(),nameCoords.end());
    std::for_each(nameCoords2.begin(),nameCoords2.end(),[](MCAuto<DataArrayChar>& elt){ elt->incrRef(); });
    MCAuto<DataArrayChar> tmp( DataArrayChar::Aggregate(ToConstVect<DataArrayChar>(nameCoords2)) );
    _name_coords = DynamicCastSafe<DataArrayChar,DataArrayAsciiChar>( tmp );
  }
  _part_coords = DataArrayPartDefinition::New( const_cast<DataArrayIdType *>(nodeIds) );
}

void MEDFileUMeshL2::sortTypes()
{
  std::set<int> mdims;
  std::vector< MCAuto<MEDFileUMeshPerType> > const tmp(_per_type_mesh[0]);
  _per_type_mesh.clear();
  for(const auto & it : tmp)
    mdims.insert(it->getDim());
  if(mdims.empty())
    return;
  int const mdim=*mdims.rbegin();
  _per_type_mesh.resize(mdim+1);
  for(int dim=mdim+1;dim!=0;dim--)
    {
      std::vector< MCAuto<MEDFileUMeshPerType> >& elt=_per_type_mesh[mdim+1-dim];
      for(const auto & it : tmp)
        if(it->getDim()==dim-1)
          elt.push_back(it);
    }
  // suppression of contiguous empty levels at the end of _per_type_mesh.
  int nbOfUselessLev=0;
  bool isFirst=true;
  for(auto it2=_per_type_mesh.rbegin();it2!=_per_type_mesh.rend();it2++)
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

void MEDFileUMeshL2::WriteCoords(med_idt fid, const std::string& mname, int dt, int it, double time, const DataArrayDouble *coords, const DataArrayIdType *famCoords, const DataArrayIdType *numCoords, const DataArrayAsciiChar *nameCoords, const DataArrayIdType *globalNumCoords)
{
  if(!coords)
    return ;
  MEDFILESAFECALLERWR0(MEDmeshNodeCoordinateWr,(fid,mname.c_str(),dt,it,time,MED_FULL_INTERLACE,ToMedInt(coords->getNumberOfTuples()),coords->begin()));
  if(famCoords)
    MEDFILESAFECALLERWR0(MEDmeshEntityFamilyNumberWr,(fid,mname.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,ToMedInt(famCoords->getNumberOfTuples()),ToMedIntArray<mcIdType>(famCoords)->begin()));
  if(numCoords)
    MEDFILESAFECALLERWR0(MEDmeshEntityNumberWr,(fid,mname.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,ToMedInt(numCoords->getNumberOfTuples()),ToMedIntArray<mcIdType>(numCoords)->begin()));
  if(nameCoords)
    {
      if(nameCoords->getNumberOfComponents()!=MED_SNAME_SIZE)
        {
          std::ostringstream oss; oss << " MEDFileUMeshL2::WriteCoords : expected a name field on nodes with number of components set to " << MED_SNAME_SIZE;
          oss << " ! The array has " << nameCoords->getNumberOfComponents() << " components !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      MEDFILESAFECALLERWR0(MEDmeshEntityNameWr,(fid,mname.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,ToMedInt(nameCoords->getNumberOfTuples()),nameCoords->begin()));
    }
  if(globalNumCoords)
    MEDFILESAFECALLERWR0(MEDmeshGlobalNumberWr,(fid,mname.c_str(),dt,it,MED_NODE,MED_NONE,ToMedInt(globalNumCoords->getNumberOfTuples()),ToMedIntArray<mcIdType>(globalNumCoords)->begin()));
}

bool MEDFileUMeshL2::isFamDefinedOnLev(int levId) const
{
  for(const auto & it : _per_type_mesh[levId])
    if(it->getFam()==nullptr)
      return false;
  return true;
}

bool MEDFileUMeshL2::isNumDefinedOnLev(int levId) const
{
  for(const auto & it : _per_type_mesh[levId])
    if(it->getNum()==nullptr)
      return false;
  return true;
}

bool MEDFileUMeshL2::isNamesDefinedOnLev(int levId) const
{
  for(const auto & it : _per_type_mesh[levId])
    if(it->getNames()==nullptr)
      return false;
  return true;
}

MEDFileCMeshL2::MEDFileCMeshL2():_ax_type(AX_CART)
{
}

void MEDFileCMeshL2::loadAll(med_idt fid, const MeshOrStructMeshCls *mId, const std::string& mName, int dt, int it)
{
  _name.set(mName.c_str());
  int nstep;
  int Mdim;
  MEDCoupling::MEDCouplingMeshType meshType;
  MEDCoupling::MEDCouplingAxisType dummy3;
  std::vector<std::string> infosOnComp(getAxisInfoOnMesh(fid,mId,mName.c_str(),meshType,dummy3,nstep,Mdim));
  if(meshType!=CARTESIAN)
    throw INTERP_KERNEL::Exception("Invalid mesh type ! You are expected a structured one whereas in file it is not a structured !");
  _time=mId->checkMeshTimeStep(fid,mName,nstep,dt,it);
  _iteration=dt;
  _order=it;
  //
  med_grid_type gridtype;
  MEDFILESAFECALLERRD0(MEDmeshGridTypeRd,(fid,mName.c_str(),&gridtype));
  if(gridtype!=MED_CARTESIAN_GRID && gridtype!=MED_POLAR_GRID)
    throw INTERP_KERNEL::Exception("Invalid rectilinear mesh ! Only cartesian and polar are supported !");
  _ax_type=TraduceAxisTypeStruct(gridtype);
  _cmesh=MEDCouplingCMesh::New();
  for(int i=0;i<Mdim;i++)
    {
      med_data_type const dataTypeReq=GetDataTypeCorrespondingToSpaceId(i);
      med_bool chgt=MED_FALSE,trsf=MED_FALSE;
      med_int const nbOfElt(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NONE,dataTypeReq,MED_NO_CMODE,&chgt,&trsf));
      MCAuto<DataArrayDouble> da=DataArrayDouble::New();
      da->alloc(nbOfElt,1);
      da->setInfoOnComponent(0,infosOnComp[i]);
      MEDFILESAFECALLERRD0(MEDmeshGridIndexCoordinateRd,(fid,mName.c_str(),dt,it,i+1,da->getPointer()));
      _cmesh->setCoordsAt(i,da);
    }
}

med_data_type MEDFileCMeshL2::GetDataTypeCorrespondingToSpaceId(int id)
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

MEDFileCLMeshL2::MEDFileCLMeshL2()
= default;

void MEDFileCLMeshL2::loadAll(med_idt fid, const MeshOrStructMeshCls *mId, const std::string& mName, int dt, int it)
{
  _name.set(mName.c_str());
  int nstep;
  int Mdim;
  MEDCoupling::MEDCouplingMeshType meshType;
  MEDCoupling::MEDCouplingAxisType dummy3;
  std::vector<std::string> const infosOnComp(getAxisInfoOnMesh(fid,mId,mName,meshType,dummy3,nstep,Mdim));
  if(meshType!=CURVE_LINEAR)
    throw INTERP_KERNEL::Exception("Invalid mesh type ! You are expected a structured one whereas in file it is not a structured !");
  _time=mId->checkMeshTimeStep(fid,mName,nstep,dt,it);
  _iteration=dt;
  _order=it;
  //
  _clmesh=MEDCouplingCurveLinearMesh::New();
  MCAuto<DataArrayMedInt> miStGrid=DataArrayMedInt::New();
  miStGrid->alloc(Mdim,1);
  MEDFILESAFECALLERRD0(MEDmeshGridStructRd,(fid,mName.c_str(),dt,it,miStGrid->getPointer()));
  MCAuto<DataArrayIdType> stGrid=FromMedIntArray<mcIdType>(miStGrid);
  _clmesh->setNodeGridStructure(stGrid->begin(),stGrid->end());
  med_bool chgt=MED_FALSE,trsf=MED_FALSE;
  med_int const nbNodes(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NONE,MED_COORDINATE,MED_NO_CMODE,&chgt,&trsf));
  MCAuto<DataArrayDouble> da=DataArrayDouble::New();
  da->alloc(nbNodes,infosOnComp.size());
  da->setInfoOnComponents(infosOnComp);
  MEDFILESAFECALLERRD0(MEDmeshNodeCoordinateRd,(fid,mName.c_str(),dt,it,MED_FULL_INTERLACE,da->getPointer()));
  _clmesh->setCoords(da);
}

MEDFileUMeshPermCompute::MEDFileUMeshPermCompute(const MEDFileUMeshSplitL1* st):_st(st),_mpt_time(0),_num_time(0)
{
}

/*!
 * Warning it returns an instance to deallocate !!!!
 */
MEDFileUMeshPermCompute::operator MEDCouplingUMesh *() const
{
  _st->_num->updateTime();
  if((MEDCouplingUMesh *)_m==nullptr)
    {
      updateTime();
      _m=static_cast<MEDCouplingUMesh *>(_st->_m_by_types.getUmesh()->deepCopy());
      _m->renumberCells(_st->_num->begin(),true);
      return _m.retn();
    }
  else
    {
      if(_mpt_time==_st->_m_by_types.getTimeOfThis() && _num_time==_st->_num->getTimeOfThis())
        return _m.retn();
      else
        {
          updateTime();
          _m=static_cast<MEDCouplingUMesh *>(_st->_m_by_types.getUmesh()->deepCopy());
          _m->renumberCells(_st->_num->begin(),true);
          return _m.retn();
        }
    }
}

void MEDFileUMeshPermCompute::operator=(MEDCouplingUMesh *m)
{
  _m=m;
}

void MEDFileUMeshPermCompute::updateTime() const
{
  _mpt_time=_st->_m_by_types.getTimeOfThis();
  _num_time=_st->_num->getTimeOfThis();
}

std::vector<const BigMemoryObject *> MEDFileUMeshPermCompute::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  ret.push_back((const MEDCouplingUMesh *)_m);
  return ret;
}

std::size_t MEDFileUMeshPermCompute::getHeapMemorySizeWithoutChildren() const
{
  return sizeof(MEDFileUMeshPermCompute);
}

MEDFileUMeshSplitL1::MEDFileUMeshSplitL1(const MEDFileUMeshSplitL1& other):RefCountObject(other),_m_by_types(other._m_by_types),_fam(other._fam),_num(other._num),_names(other._names),_rev_num(other._rev_num),_m(this)
{
}

MEDFileUMeshSplitL1::MEDFileUMeshSplitL1(const MEDFileUMeshL2& l2, const std::string&  /*mName*/, int id):_m(this)
{
  const std::vector< MCAuto<MEDFileUMeshPerType> >& v=l2.getLev(id);
  if(v.empty())
    return;
  std::size_t const sz=v.size();
  std::vector<const MEDCoupling1GTUMesh *> ms(sz);
  std::vector<const DataArrayIdType *> fams(sz),nums(sz);
  std::vector<const DataArrayChar *> names(sz);
  std::vector<const PartDefinition *> pds(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      MEDCoupling1GTUMesh *elt(v[i]->getMesh());
      MCAuto<DataArrayDouble> tmp2=l2.getCoords();
      elt->setCoords(tmp2);
      ms[i]=elt;
      pds[i]=v[i]->getPartDef();
    }
  _m_by_types.assignParts(ms);
  _m_by_types.assignDefParts(pds);
  if(l2.isFamDefinedOnLev(id))
    {
      for(std::size_t i=0;i<sz;i++)
        fams[i]=v[i]->getFam();
      if(sz!=1)
        _fam=DataArrayIdType::Aggregate(fams);
      else
        {
          fams[0]->incrRef();
          _fam=const_cast<DataArrayIdType *>(fams[0]);
        }
    }
  if(l2.isNumDefinedOnLev(id))
    {
      for(std::size_t i=0;i<sz;i++)
        nums[i]=v[i]->getNum();
      if(sz!=1)
        _num=DataArrayIdType::Aggregate(nums);
      else
        {
          nums[0]->incrRef();
          _num=const_cast<DataArrayIdType *>(nums[0]);
        }
      computeRevNum();
    }
  if(l2.isNamesDefinedOnLev(id))
    {
      for(std::size_t i=0;i<sz;i++)
        names[i]=v[i]->getNames();
      _names=dynamic_cast<DataArrayAsciiChar *>(DataArrayChar::Aggregate(names));
    }
}

MEDFileUMeshSplitL1::MEDFileUMeshSplitL1(MEDCoupling1GTUMesh *m):_m(this)
{
  std::vector< const MEDCoupling1GTUMesh * > v(1);
  v[0]=m;
  assignParts(v);
}

MEDFileUMeshSplitL1::MEDFileUMeshSplitL1(MEDCouplingUMesh *m):_m(this)
{
  assignMesh(m,true);
}

MEDFileUMeshSplitL1::MEDFileUMeshSplitL1(MEDCouplingUMesh *m, bool newOrOld):_m(this)
{
  assignMesh(m,newOrOld);
}

void MEDFileUMeshSplitL1::setName(const std::string& name)
{
  _m_by_types.setName(name);
}

std::size_t MEDFileUMeshSplitL1::getHeapMemorySizeWithoutChildren() const
{
  return 0;
}

std::vector<const BigMemoryObject *> MEDFileUMeshSplitL1::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  ret.push_back(&_m_by_types);
  ret.push_back(&_m);
  ret.push_back((const DataArrayIdType*)_fam);
  ret.push_back((const DataArrayIdType*)_num);
  ret.push_back((const DataArrayIdType*)_rev_num);
  ret.push_back((const DataArrayAsciiChar*)_names);
  return ret;
}

MEDFileUMeshSplitL1 *MEDFileUMeshSplitL1::shallowCpyUsingCoords(DataArrayDouble *coords) const
{
  MCAuto<MEDFileUMeshSplitL1> ret(new MEDFileUMeshSplitL1(*this));
  ret->_m_by_types.shallowCpyMeshes();
  ret->_m_by_types.setCoords(coords);
  return ret.retn();
}

MEDFileUMeshSplitL1 *MEDFileUMeshSplitL1::deepCopy(DataArrayDouble *coords) const
{
  MCAuto<MEDFileUMeshSplitL1> ret(new MEDFileUMeshSplitL1(*this));
  ret->_m_by_types=_m_by_types.deepCopy(coords);
  if((const DataArrayIdType *)_fam)
    ret->_fam=_fam->deepCopy();
  if((const DataArrayIdType *)_num)
    ret->_num=_num->deepCopy();
  if((const DataArrayIdType *)_rev_num)
    ret->_rev_num=_rev_num->deepCopy();
  if((const DataArrayAsciiChar *)_names)
    ret->_names=_names->deepCopy();
  return ret.retn();
}

void MEDFileUMeshSplitL1::checkConsistency() const
{
  if (!_fam || _fam->getNumberOfTuples() != getSize())
    throw INTERP_KERNEL::Exception("MEDFileUMeshSplitL1::checkConsistency(): internal family array has an invalid size!");
  mcIdType const nbCells = getSize();
  if (_num)
    {
      _num->checkNbOfTuplesAndComp(nbCells,1,"MEDFileUMeshSplitL1::checkConsistency(): inconsistent internal node numbering array!");
      mcIdType pos;
      mcIdType const maxValue=_num->getMaxValue(pos);
      if (!_rev_num || _rev_num->getNumberOfTuples() != (maxValue+1))
        throw INTERP_KERNEL::Exception("MEDFileUMeshSplitL1::checkConsistency(): inconsistent internal revert node numbering array!");
    }
  if ((_num && !_rev_num) || (!_num && _rev_num))
    throw INTERP_KERNEL::Exception("MEDFileUMeshSplitL1::checkConsistency(): inconsistent internal numbering arrays (one is null)!");
  if (_num && !_num->hasUniqueValues())
    throw INTERP_KERNEL::Exception("MEDFileUMeshSplitL1::checkConsistency(): inconsistent internal node numbering array: duplicates found!");
  if (_names)
    _names->checkNbOfTuplesAndComp(nbCells,1,"MEDFileUMeshSplitL1::checkConsistency(): internal cell naming array has an invalid size!");

  _m_by_types.checkConsistency();
}

bool MEDFileUMeshSplitL1::isEqual(const MEDFileUMeshSplitL1 *other, double eps, std::string& what) const
{
  if(!_m_by_types.isEqual(other->_m_by_types,eps,what))
    return false;
  const DataArrayIdType *d1=_fam;
  const DataArrayIdType *d2=other->_fam;
  if((d1==nullptr && d2!=nullptr) || (d1!=nullptr && d2==nullptr))
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
  if((d1==nullptr && d2!=nullptr) || (d1!=nullptr && d2==nullptr))
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
  const DataArrayAsciiChar *e1=_names;
  const DataArrayAsciiChar *e2=other->_names;
  if((e1==nullptr && e2!=nullptr) || (e1!=nullptr && e2==nullptr))
    {
      what="Presence of cell names arr in one sublevel and not in other!";
      return false;
    }
  if(e1)
    if(!e1->isEqual(*e2))
      {
        what="Name cell arr at a sublevel are not deeply equal !";
        return false;
      }
  return true;
}

void MEDFileUMeshSplitL1::synchronizeTinyInfo(const MEDFileMesh& master) const
{
  _m_by_types.synchronizeTinyInfo(master);
}

void MEDFileUMeshSplitL1::clearNonDiscrAttributes() const
{
  _m_by_types.clearNonDiscrAttributes();
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

void MEDFileUMeshSplitL1::setCoords(DataArrayDouble *coords)
{
  _m_by_types.setCoords(coords);
}

void MEDFileUMeshSplitL1::assignMesh(MEDCouplingUMesh *m, bool newOrOld)
{
  if(newOrOld)
    {
      m->incrRef();
      _m=m;
      _m_by_types.assignUMesh(dynamic_cast<MEDCouplingUMesh *>(m->deepCopy()));
      MCAuto<DataArrayIdType> da=_m_by_types.getUmesh()->getRenumArrForConsecutiveCellTypesSpec(typmai2,typmai2+MED_N_CELL_FIXED_GEO);
      if(!da->isIota(m->getNumberOfCells()))
        {
          _num=da->invertArrayO2N2N2O(m->getNumberOfCells());
          _m.updateTime();
          computeRevNum();
          _m_by_types.getUmesh()->renumberCells(da->begin(),false);
        }
    }
  else
    {
      if(!m->checkConsecutiveCellTypesAndOrder(typmai2,typmai2+MED_N_CELL_FIXED_GEO))
        throw INTERP_KERNEL::Exception("MEDFileUMeshSplitL1::assignMesh(): the mesh does not follow the MED file numbering convention! Invoke sortCellsInMEDFileFrmt() first!");
      m->incrRef();
      _m_by_types.assignUMesh(m);
    }
  assignCommonPart();
}

void MEDFileUMeshSplitL1::forceComputationOfParts() const
{
  _m_by_types.forceComputationOfPartsFromUMesh();
}

void MEDFileUMeshSplitL1::declarePartsUpdated() const
{
  _m_by_types.declarePartsUpdated();
}

void MEDFileUMeshSplitL1::assignParts(const std::vector< const MEDCoupling1GTUMesh * >& mParts)
{
  _m_by_types.assignParts(mParts);
  assignCommonPart();
}

MEDFileUMeshSplitL1::MEDFileUMeshSplitL1():_m(this)
{
}

void MEDFileUMeshSplitL1::assignCommonPart()
{
  _fam=DataArrayIdType::New();
  _fam->alloc(_m_by_types.getSize(),1);
  _fam->fillWithValue(0);
}

bool MEDFileUMeshSplitL1::empty() const
{
  return _m_by_types.empty();
}

bool MEDFileUMeshSplitL1::presenceOfOneFams(const std::vector<mcIdType>& ids) const
{
  const DataArrayIdType *fam=_fam;
  if(!fam)
    return false;
  return fam->presenceOfValue(ids);
}

int MEDFileUMeshSplitL1::getMeshDimension() const
{
  return _m_by_types.getMeshDimension();
}

void MEDFileUMeshSplitL1::simpleRepr(std::ostream& oss) const
{
  std::vector<mcIdType> code=_m_by_types.getDistributionOfTypes();
  std::size_t const nbOfTypes=code.size()/3;
  for(std::size_t i=0;i<nbOfTypes;i++)
    {
      auto typ=(INTERP_KERNEL::NormalizedCellType) code[3*i];
      oss << "    - Number of cells with type " << INTERP_KERNEL::CellModel::GetCellModel(typ).getRepr() << " : " << code[3*i+1] << std::endl;
    }
}

mcIdType MEDFileUMeshSplitL1::getSize() const
{
  return _m_by_types.getSize();
}

MEDCouplingUMesh *MEDFileUMeshSplitL1::getFamilyPart(const mcIdType *idsBg, const mcIdType *idsEnd, bool renum) const
{
  MCAuto<DataArrayIdType> eltsToKeep=_fam->findIdsEqualList(idsBg,idsEnd);
  auto *m=(MEDCouplingUMesh *)_m_by_types.getUmesh()->buildPartOfMySelf(eltsToKeep->begin(),eltsToKeep->end(),true);
  if(renum)
    return renumIfNeeded(m,eltsToKeep->begin());
  return m;
}

DataArrayIdType *MEDFileUMeshSplitL1::getFamilyPartArr(const mcIdType *idsBg, const mcIdType *idsEnd, bool renum) const
{
  MCAuto<DataArrayIdType> da=_fam->findIdsEqualList(idsBg,idsEnd);
  if(renum)
    return renumIfNeededArr(da);
  return da.retn();
}

std::vector<INTERP_KERNEL::NormalizedCellType> MEDFileUMeshSplitL1::getGeoTypes() const
{
  return _m_by_types.getGeoTypes();
}

mcIdType MEDFileUMeshSplitL1::getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType ct) const
{
  return _m_by_types.getNumberOfCellsWithType(ct);
}

MEDCouplingUMesh *MEDFileUMeshSplitL1::getWholeMesh(bool renum) const
{
  MCAuto<MEDCouplingUMesh> tmp;
  if(renum && ((const DataArrayIdType *)_num))
    tmp=_m;
  else
    { tmp=_m_by_types.getUmesh(); if(tmp) tmp->incrRef(); }
  return tmp.retn();
}

mcIdType MEDFileUMeshSplitL1::getNumberOfCells() const
{
  return _m_by_types.getNumberOfCells();
}

DataArrayIdType *MEDFileUMeshSplitL1::extractFamilyFieldOnGeoType(INTERP_KERNEL::NormalizedCellType gt) const
{
  const DataArrayIdType *fam(_fam);
  if(!fam)
    return nullptr;
  mcIdType start(0),stop(0);
  _m_by_types.getStartStopOfGeoTypeWithoutComputation(gt,start,stop);
  return fam->selectByTupleIdSafeSlice(start,stop,1);
}

DataArrayIdType *MEDFileUMeshSplitL1::extractNumberFieldOnGeoType(INTERP_KERNEL::NormalizedCellType gt) const
{
  const DataArrayIdType *num(_num);
  if(!num)
    return nullptr;
  mcIdType start(0),stop(0);
  _m_by_types.getStartStopOfGeoTypeWithoutComputation(gt,start,stop);
  return num->selectByTupleIdSafeSlice(start,stop,1);
}

DataArrayIdType *MEDFileUMeshSplitL1::getOrCreateAndGetFamilyField()
{
  if((DataArrayIdType *)_fam)
    return _fam;
  mcIdType const nbOfTuples=_m_by_types.getSize();
  _fam=DataArrayIdType::New(); _fam->alloc(nbOfTuples,1); _fam->fillWithZero();
  return _fam;
}

const DataArrayIdType *MEDFileUMeshSplitL1::getFamilyField() const
{
  return _fam;
}

const DataArrayIdType *MEDFileUMeshSplitL1::getNumberField() const
{
  return _num;
}

const DataArrayIdType *MEDFileUMeshSplitL1::getRevNumberField() const
{
  return _rev_num;
}

const DataArrayAsciiChar *MEDFileUMeshSplitL1::getNameField() const
{
  return _names;
}

const PartDefinition *MEDFileUMeshSplitL1::getPartDef(INTERP_KERNEL::NormalizedCellType gt) const
{
  return _m_by_types.getPartDefOfWithoutComputation(gt);
}

void MEDFileUMeshSplitL1::eraseFamilyField()
{
  _fam->fillWithZero();
}

/*!
 * This method ignores _m and _m_by_types.
 */
void MEDFileUMeshSplitL1::setGroupsFromScratch(const std::vector<const MEDCouplingUMesh *>& ms, std::map<std::string,mcIdType>& familyIds,
                                               std::map<std::string, std::vector<std::string> >&  /*groups*/)
{
  std::vector< DataArrayIdType * > corr;
  _m=MEDCouplingUMesh::FuseUMeshesOnSameCoords(ms,0,corr);
  std::vector< MCAuto<DataArrayIdType> > const corrMSafe(corr.begin(),corr.end());
  std::vector< std::vector<mcIdType> > fidsOfGroups;
  std::vector< const DataArrayIdType * > const corr2(corr.begin(),corr.end());
  _fam=DataArrayIdType::MakePartition(corr2,((MEDCouplingUMesh *)_m)->getNumberOfCells(),fidsOfGroups);
  mcIdType const nbOfCells=((MEDCouplingUMesh *)_m)->getNumberOfCells();
  std::map<mcIdType,std::string> newfams;
  std::map<mcIdType,mcIdType> famIdTrad;
  TraduceFamilyNumber(fidsOfGroups,familyIds,famIdTrad,newfams);
  mcIdType *w=_fam->getPointer();
  for(mcIdType i=0;i<nbOfCells;i++,w++)
    *w=famIdTrad[*w];
}

void MEDFileUMeshSplitL1::checkCoordsConsistency(const DataArrayDouble *coords) const
{
  std::vector<MEDCoupling1GTUMesh *> const ms(_m_by_types.getParts());
  for(auto mesh : ms)
  {
    if(mesh)
      if(mesh->getCoords() != coords)
        mesh->getCoords()->checkNbOfTuplesAndComp(*coords,"MEDFileUMeshSplitL1::checkCoordsConsistency : mismatch between coordinates instance in MEDFileUMesh and instance in subparts");
  }
}

void MEDFileUMeshSplitL1::write(med_idt fid, const std::string& mName, int mdim) const
{
  std::vector<MEDCoupling1GTUMesh *> const ms(_m_by_types.getParts());
  mcIdType start=0;
  for(auto m : ms)
    {
      mcIdType const nbCells=m->getNumberOfCells();
      mcIdType const end=start+nbCells;
      MCAuto<DataArrayIdType> fam,num;
      MCAuto<DataArrayAsciiChar> names;
      if((const DataArrayIdType *)_fam)
        fam=_fam->subArray(start,end);
      if((const DataArrayIdType *)_num)
        num=_num->subArray(start,end);
      if((const DataArrayAsciiChar *)_names)
        names=static_cast<DataArrayAsciiChar *>(_names->subArray(start,end));
      MEDFileUMeshPerType::Write(fid,mName,mdim,m,fam,num,names);
      start=end;
    }
}

void MEDFileUMeshSplitL1::renumberNodesInConn(const mcIdType *newNodeNumbersO2N)
{
  _m_by_types.renumberNodesInConnWithoutComputation(newNodeNumbersO2N);
}

void MEDFileUMeshSplitL1::serialize(std::vector<mcIdType>& tinyInt, std::vector< MCAuto<DataArrayIdType> >& bigArraysI) const
{
  bigArraysI.push_back(_fam);
  bigArraysI.push_back(_num);
  _m_by_types.serialize(tinyInt,bigArraysI);
}

void MEDFileUMeshSplitL1::unserialize(const std::string& name, DataArrayDouble *coo, std::vector<mcIdType>& tinyInt, std::vector< MCAuto<DataArrayIdType> >& bigArraysI)
{
  _fam=bigArraysI.back(); bigArraysI.pop_back();
  _num=bigArraysI.back(); bigArraysI.pop_back();
  _m_by_types.unserialize(name,coo,tinyInt,bigArraysI);
}

void MEDFileUMeshSplitL1::changeFamilyIdArr(mcIdType oldId, mcIdType newId)
{
  DataArrayIdType *arr=_fam;
  if(arr)
    arr->changeValue(oldId,newId);
}

void MEDFileUMeshSplitL1::setFamilyArr(DataArrayIdType *famArr)
{
  if(!famArr)
    {
      _fam=nullptr;
      return ;
    }
  mcIdType const sz(_m_by_types.getSize());
  famArr->checkNbOfTuplesAndComp(sz,1,"MEDFileUMeshSplitL1::setFamilyArr : Problem in size of Family arr ! ");
  famArr->incrRef();
  _fam=famArr;
}

DataArrayIdType *MEDFileUMeshSplitL1::getFamilyField()
{
  return _fam;
}

void MEDFileUMeshSplitL1::setRenumArr(DataArrayIdType *renumArr)
{
  if(!renumArr)
    {
      _num=nullptr;
      _rev_num=nullptr;
      return ;
    }
  mcIdType const sz(_m_by_types.getSize());
  renumArr->checkNbOfTuplesAndComp(sz,1,"MEDFileUMeshSplitL1::setRenumArr : Problem in size of numbering arr ! ");
  renumArr->incrRef();
  _num=renumArr;
  computeRevNum();
}

void MEDFileUMeshSplitL1::setNameArr(DataArrayAsciiChar *nameArr)
{
  if(!nameArr)
    {
      _names=nullptr;
      return ;
    }
  mcIdType const sz(_m_by_types.getSize());
  nameArr->checkNbOfTuplesAndComp(sz,MED_SNAME_SIZE,"MEDFileUMeshSplitL1::setNameArr : Problem in size of name arr ! ");
  nameArr->incrRef();
  _names=nameArr;
}

MEDCouplingUMesh *MEDFileUMeshSplitL1::Renumber2(const DataArrayIdType *renum, MEDCouplingUMesh *m, const mcIdType *cellIds)
{
  if(renum==nullptr)
    return m;
  if(cellIds==nullptr)
    m->renumberCells(renum->begin(),true);
  else
    {
      MCAuto<DataArrayIdType> locnum=renum->selectByTupleId(cellIds,cellIds+m->getNumberOfCells());
      m->renumberCells(locnum->begin(),true);
    }
  return m;
}

MEDFileUMeshSplitL1 *MEDFileUMeshSplitL1::Unserialize(const std::string& name, DataArrayDouble *coo, std::vector<mcIdType>& tinyInt, std::vector< MCAuto<DataArrayIdType> >& bigArraysI)
{
  MCAuto<MEDFileUMeshSplitL1> ret(new MEDFileUMeshSplitL1);
  ret->unserialize(name,coo,tinyInt,bigArraysI);
  return ret.retn();
}

MEDCouplingUMesh *MEDFileUMeshSplitL1::renumIfNeeded(MEDCouplingUMesh *m, const mcIdType *cellIds) const
{
  return Renumber2(_num,m,cellIds);
}

DataArrayIdType *MEDFileUMeshSplitL1::Renumber(const DataArrayIdType *renum, const DataArrayIdType *da)
{
  if((const DataArrayIdType *)renum==nullptr)
    {
      da->incrRef();
      return const_cast<DataArrayIdType *>(da);
    }
  return renum->selectByTupleId(da->begin(),da->end());
}

DataArrayIdType *MEDFileUMeshSplitL1::renumIfNeededArr(const DataArrayIdType *da) const
{
  return Renumber(_num,da);
}

std::vector<mcIdType> MEDFileUMeshSplitL1::GetNewFamiliesNumber(mcIdType nb, const std::map<std::string,mcIdType>& families)
{
  mcIdType id=-1;
  for(const auto & familie : families)
    id=std::max(id,familie.second);
  if(id==-1)
    id=0;
  std::vector<mcIdType> ret(nb);
  for(mcIdType i=1;i<=nb;i++)
    ret[i]=id+i;
  return ret;
}

void MEDFileUMeshSplitL1::TraduceFamilyNumber(const std::vector< std::vector<mcIdType> >&  /*fidsGrps*/, std::map<std::string,mcIdType>&  /*familyIds*/,
                                              std::map<mcIdType,mcIdType>&  /*famIdTrad*/, std::map<mcIdType,std::string>&  /*newfams*/)
{
  std::set<mcIdType> const allfids;
  //tony
}

void MEDFileUMeshSplitL1::computeRevNum() const
{
  mcIdType pos;
  if(!_num->empty())
  {
    mcIdType const maxValue=_num->getMaxValue(pos);
    _rev_num=_num->invertArrayN2O2O2N(maxValue+1);
  }
  else
  {
    _rev_num = DataArrayIdType::New();
    _rev_num->alloc(0,1);
  }
  
}

//=

MEDFileUMeshAggregateCompute::MEDFileUMeshAggregateCompute():_mp_time(0),_m_time(0)
{
}

void MEDFileUMeshAggregateCompute::setName(const std::string& name)
{
  if(_m_time>=_mp_time)
    {
      MEDCouplingUMesh *um(_m);
      if(um)
        um->setName(name);
    }
  if(_mp_time>=_m_time)
    {
      for(auto & _m_part : _m_parts)
        {
          MEDCoupling1GTUMesh *tmp(_m_part);
          if(tmp)
            tmp->setName(name);
        }
    }
}

void MEDFileUMeshAggregateCompute::assignParts(const std::vector< const MEDCoupling1GTUMesh * >& mParts)
{
  std::size_t const sz(mParts.size());
  std::vector< MCAuto<MEDCoupling1GTUMesh> > ret(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      const MEDCoupling1GTUMesh *elt(mParts[i]);
      if(!elt)
        throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::assignParts : presence of null pointer !");
      ret[i]=const_cast<MEDCoupling1GTUMesh *>(elt); elt->incrRef();
    }
  _m_parts=ret;
  _part_def.clear(); _part_def.resize(sz);
  _mp_time=std::max(_mp_time,_m_time)+1;
  _m=nullptr;
}

void MEDFileUMeshAggregateCompute::assignDefParts(const std::vector<const PartDefinition *>& partDefs)
{
  if(_mp_time<_m_time)
    throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::assignDefParts : the parts require a computation !");
  std::size_t const sz(partDefs.size());
  if(_part_def.size()!=partDefs.size() || _part_def.size()!=_m_parts.size())
    throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::assignDefParts : sizes of vectors of part definition mismatch !");
  for(std::size_t i=0;i<sz;i++)
    {
      const PartDefinition *elt(partDefs[i]);
      if(elt)
        elt->incrRef();
      _part_def[i]=const_cast<PartDefinition*>(elt);
    }
}

void MEDFileUMeshAggregateCompute::assignUMesh(MEDCouplingUMesh *m)
{
  _m=m;
  _m_parts.clear();
  _m_time=std::max(_mp_time,_m_time)+1;
}

MEDCouplingUMesh *MEDFileUMeshAggregateCompute::getUmesh() const
{
  if(_mp_time<=_m_time)
    return _m;
  std::vector< const MEDCoupling1GTUMesh *> mp(_m_parts.size());
  std::copy(_m_parts.begin(),_m_parts.end(),mp.begin());
  _m=MEDCoupling1GTUMesh::AggregateOnSameCoordsToUMesh(mp);
  _m_parts.clear();//to avoid memory peak !
  _m_time=_mp_time+1;//+1 is important ! That is to say that only _m is OK not _m_parts because cleared !
  return _m;
}

mcIdType MEDFileUMeshAggregateCompute::getNumberOfCells() const
{
  if(_mp_time<=_m_time)
    return _m->getNumberOfCells();
  mcIdType ret(0);
  for(const auto & _m_part : _m_parts)
    ret+=_m_part->getNumberOfCells();
  return ret;
}

std::vector<INTERP_KERNEL::NormalizedCellType> MEDFileUMeshAggregateCompute::getGeoTypes() const
{
  if(_mp_time>=_m_time)
    {
      std::size_t const sz(_m_parts.size());
      std::vector<INTERP_KERNEL::NormalizedCellType> ret(sz);
      for(std::size_t i=0;i<sz;i++)
        ret[i]=_m_parts[i]->getCellModelEnum();
      return ret;
    }
  else
    return _m->getAllGeoTypesSorted();
}

mcIdType MEDFileUMeshAggregateCompute::getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType ct) const
{
  if(_mp_time>=_m_time)
    {
      for(const auto & _m_part : _m_parts)
        {
          const MEDCoupling1GTUMesh *elt(_m_part);
          if(elt && elt->getCellModelEnum()==ct)
            return elt->getNumberOfCells();
        }
      return 0;
    }
  else
    return _m->getNumberOfCellsWithType(ct);
}

std::vector<MEDCoupling1GTUMesh *> MEDFileUMeshAggregateCompute::retrievePartsWithoutComputation() const
{
  if(_mp_time<_m_time)
    throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::getPartsWithoutComputation : the parts require a computation !");
  //
  std::vector<MEDCoupling1GTUMesh *> ret(_m_parts.size());
  std::size_t i(0);
  for(std::vector< MCAuto<MEDCoupling1GTUMesh> >::const_iterator it=_m_parts.begin();it!=_m_parts.end();it++,i++)
    {
      const MEDCoupling1GTUMesh *elt(*it);
      ret[i]=const_cast<MEDCoupling1GTUMesh *>(elt);
    }
  return ret;
}

std::vector<MEDCoupling1GTUMesh *> MEDFileUMeshAggregateCompute::getParts() const
{
  if(_mp_time<_m_time)
    forceComputationOfPartsFromUMesh();
  return retrievePartsWithoutComputation();
}

void MEDFileUMeshAggregateCompute::highlightUsedNodes(std::vector<bool>& nodesToBeHighlighted) const
{
  if(_mp_time<_m_time)
    forceComputationOfPartsFromUMesh();
  for(auto part : this->_m_parts)
  {
    part->computeNodeIdsAlg(nodesToBeHighlighted);
  }
}

MEDCoupling1GTUMesh *MEDFileUMeshAggregateCompute::retrievePartWithoutComputation(INTERP_KERNEL::NormalizedCellType gt) const
{
  std::vector<MEDCoupling1GTUMesh *> v(retrievePartsWithoutComputation());
  std::size_t const sz(v.size());
  for(std::size_t i=0;i<sz;i++)
    {
      if(v[i])
        if(v[i]->getCellModelEnum()==gt)
          return v[i];
    }
  throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::getPartWithoutComputation : the geometric type is not existing !");
}

void MEDFileUMeshAggregateCompute::getStartStopOfGeoTypeWithoutComputation(INTERP_KERNEL::NormalizedCellType gt, mcIdType& start, mcIdType& stop) const
{
  start=0; stop=0;
  std::vector<MEDCoupling1GTUMesh *> v(retrievePartsWithoutComputation());
  std::size_t const sz(v.size());
  for(std::size_t i=0;i<sz;i++)
    {
      if(v[i])
        {
          if(v[i]->getCellModelEnum()==gt)
            {
              stop=start+v[i]->getNumberOfCells();
              return;
            }
          else
            start+=v[i]->getNumberOfCells();
        }
    }
  throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::getStartStopOfGeoTypeWithoutComputation : the geometric type is not existing !");
}

void MEDFileUMeshAggregateCompute::renumberNodesInConnWithoutComputation(const mcIdType *newNodeNumbersO2N)
{
  if(_mp_time>_m_time)
    {
      for(auto & _m_part : _m_parts)
        {
          MEDCoupling1GTUMesh *m(_m_part);
          if(m)
            m->renumberNodesInConn(newNodeNumbersO2N);
        }
    }
  else
    {
      MEDCouplingUMesh *m(getUmesh());
      if(!m)
        return;
      m->renumberNodesInConn(newNodeNumbersO2N);
      // if _mp_time == _m_time notify for future clients that _m_parts is obsolete
      _m_parts.clear();
      _m_time = std::max(_m_time,_mp_time+1);
    }
}

void MEDFileUMeshAggregateCompute::forceComputationOfPartsFromUMesh() const
{
  const MEDCouplingUMesh *m(_m);
  if(!m)
    {
      if(_m_parts.empty())
        throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::forceComputationOfPartsFromUMesh : null UMesh !");
      else
        return ;// no needs to compte parts they are already here !
    }
  std::vector<MEDCouplingUMesh *> ms(m->splitByType());
  std::vector< MCAuto<MEDCouplingUMesh> > const msMSafe(ms.begin(),ms.end());
  std::size_t const sz(msMSafe.size());
  _m_parts.resize(sz);
  for(std::size_t i=0;i<sz;i++)
    _m_parts[i]=MEDCoupling1GTUMesh::New(ms[i]);
  _part_def.clear();
  _part_def.resize(_m_parts.size());
  _mp_time=std::max(_mp_time,_m_time);
}

void MEDFileUMeshAggregateCompute::declarePartsUpdated() const
{
  _mp_time=std::max(_mp_time,_m_time) + 1;
  _m.nullify();
}

const PartDefinition *MEDFileUMeshAggregateCompute::getPartDefOfWithoutComputation(INTERP_KERNEL::NormalizedCellType gt) const
{
  if(_mp_time<_m_time)
    throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::getPartDefOfWithoutComputation : the parts require a computation !");
  if(_m_parts.size()!=_part_def.size())
    throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::getPartDefOfWithoutComputation : size of arrays are expected to be the same !");
  std::size_t const sz(_m_parts.size());
  for(std::size_t i=0;i<sz;i++)
    {
      const MEDCoupling1GTUMesh *mesh(_m_parts[i]);
      if(mesh)
        if(mesh->getCellModelEnum()==gt)
          return _part_def[i];
    }
  throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::getPartDefOfWithoutComputation : The input geo type is not existing in this !");
}

void MEDFileUMeshAggregateCompute::serialize(std::vector<mcIdType>& tinyInt, std::vector< MCAuto<DataArrayIdType> >& bigArraysI) const
{
  if(_mp_time<_m_time)
    throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::serialize : the parts require a computation !");
  std::size_t const sz(_m_parts.size());
  tinyInt.push_back((mcIdType)sz);
  for(std::size_t i=0;i<sz;i++)
    {
      const MEDCoupling1GTUMesh *mesh(_m_parts[i]);
      if(!mesh)
        throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::serialize : one part is empty !");
      tinyInt.push_back(mesh->getCellModelEnum());
      const auto *mesh1(dynamic_cast<const MEDCoupling1SGTUMesh *>(mesh));
      const auto *mesh2(dynamic_cast<const MEDCoupling1DGTUMesh *>(mesh));
      if(mesh1)
        {
          DataArrayIdType *elt(mesh1->getNodalConnectivity());
          if(elt)
            elt->incrRef();
          MCAuto<DataArrayIdType> const elt1(elt);
          bigArraysI.push_back(elt1);
        }
      else if(mesh2)
        {
          DataArrayIdType *elt1(mesh2->getNodalConnectivity()),*elt2(mesh2->getNodalConnectivityIndex());
          if(elt1)
            elt1->incrRef();
          if(elt2)
            elt2->incrRef();
          MCAuto<DataArrayIdType> elt11(elt1),elt22(elt2);
          bigArraysI.push_back(elt11); bigArraysI.push_back(elt22);
        }
      else
        throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::serialize : unrecognized single geo type mesh !");
      const PartDefinition *pd(_part_def[i]);
      if(!pd)
        tinyInt.push_back(-1);
      else
        {
          std::vector<mcIdType> tinyTmp;
          pd->serialize(tinyTmp,bigArraysI);
          tinyInt.push_back((mcIdType)tinyTmp.size());
          tinyInt.insert(tinyInt.end(),tinyTmp.begin(),tinyTmp.end());
        }
    }
}

void MEDFileUMeshAggregateCompute::unserialize(const std::string& name, DataArrayDouble *coo, std::vector<mcIdType>& tinyInt, std::vector< MCAuto<DataArrayIdType> >& bigArraysI)
{
  mcIdType const nbParts(tinyInt.back()); tinyInt.pop_back();
  _part_def.clear(); _part_def.resize(nbParts);
  _m_parts.clear(); _m_parts.resize(nbParts);
  for(mcIdType i=0;i<nbParts;i++)
    {
      auto tp((INTERP_KERNEL::NormalizedCellType) tinyInt.back()); tinyInt.pop_back();
      MCAuto<MEDCoupling1GTUMesh> mesh(MEDCoupling1GTUMesh::New(name,tp));
      mesh->setCoords(coo);
      auto *mesh1(dynamic_cast<MEDCoupling1SGTUMesh *>((MEDCoupling1GTUMesh *) mesh));
      auto *mesh2(dynamic_cast<MEDCoupling1DGTUMesh *>((MEDCoupling1GTUMesh *) mesh));
      if(mesh1)
        {
          mesh1->setNodalConnectivity(bigArraysI.back()); bigArraysI.pop_back();
        }
      else if(mesh2)
        {
          MCAuto<DataArrayIdType> elt0,elt1;
          elt0=bigArraysI.back(); bigArraysI.pop_back();
          elt1=bigArraysI.back(); bigArraysI.pop_back();
          mesh2->setNodalConnectivity(elt0,elt1);
        }
      else
        throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::unserialize : unrecognized single geo type mesh !");
      _m_parts[i]=mesh;
      mcIdType const pdid(tinyInt.back()); tinyInt.pop_back();
      if(pdid!=-1)
        _part_def[i]=PartDefinition::Unserialize(tinyInt,bigArraysI);
      _mp_time=std::max(_mp_time,_m_time)+1;
    }
}

/*!
 * This method returns true if \a this is stored split by type false if stored in a merged unstructured mesh.
 */
bool MEDFileUMeshAggregateCompute::isStoredSplitByType() const
{
  return _mp_time>=_m_time;
}

std::size_t MEDFileUMeshAggregateCompute::getTimeOfThis() const
{
  if(_mp_time>_m_time)
    return getTimeOfParts();
  if(_m_time>_mp_time)
    return getTimeOfUMesh();
  return std::max(getTimeOfParts(),getTimeOfUMesh());
}

std::size_t MEDFileUMeshAggregateCompute::getTimeOfParts() const
{
  std::size_t ret(0);
  for(const auto & _m_part : _m_parts)
    {
      const MEDCoupling1GTUMesh *elt(_m_part);
      if(!elt)
        throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::getTimeOfParts : null obj in parts !");
      ret=std::max(ret,elt->getTimeOfThis());
    }
  if(ret==0)
    throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::getTimeOfParts : parts is empty !");
  return ret;
}

std::size_t MEDFileUMeshAggregateCompute::getTimeOfUMesh() const
{
  const MEDCouplingUMesh *m(_m);
  if(!m)
    throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::getTimeOfUMesh : unmesh is null !");
  return m->getTimeOfThis();
}

std::size_t MEDFileUMeshAggregateCompute::getHeapMemorySizeWithoutChildren() const
{
  std::size_t const ret(_m_parts.size()*sizeof(MCAuto<MEDCoupling1GTUMesh>));
  return ret;
}

std::vector<const BigMemoryObject *> MEDFileUMeshAggregateCompute::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  for(const auto & _m_part : _m_parts)
    ret.push_back((const MEDCoupling1GTUMesh *)_m_part);
  ret.push_back((const MEDCouplingUMesh *)_m);
  return ret;
}

MEDFileUMeshAggregateCompute MEDFileUMeshAggregateCompute::deepCopy(DataArrayDouble *coords) const
{
  MEDFileUMeshAggregateCompute ret;
  ret._m_parts.resize(_m_parts.size());
  for(std::size_t i=0;i<_m_parts.size();i++)
    {
      const MEDCoupling1GTUMesh *elt(_m_parts[i]);
      if(elt)
        {
          ret._m_parts[i]=static_cast<MEDCoupling::MEDCoupling1GTUMesh*>(elt->deepCopy());
          ret._m_parts[i]->setCoords(coords);
        }
    }
  ret._mp_time=_mp_time; ret._m_time=_m_time;
  if((const MEDCouplingUMesh *)_m)
    {
      ret._m=static_cast<MEDCoupling::MEDCouplingUMesh*>(_m->deepCopy());
      ret._m->setCoords(coords);
    }
  std::size_t const sz(_part_def.size());
  ret._part_def.clear(); ret._part_def.resize(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      const PartDefinition *elt(_part_def[i]);
      if(elt)
        ret._part_def[i]=elt->deepCopy();
    }
  return ret;
}

void MEDFileUMeshAggregateCompute::shallowCpyMeshes()
{
  for(auto & _m_part : _m_parts)
    {
      const MEDCoupling1GTUMesh *elt(_m_part);
      if(elt)
        {
          MCAuto<MEDCouplingMesh> elt2(elt->clone(false));
          _m_part=DynamicCastSafe<MEDCouplingMesh,MEDCoupling1GTUMesh>(elt2);
        }
    }
  const MEDCouplingUMesh *m(_m);
  if(m)
    _m=m->clone(false);
}

bool MEDFileUMeshAggregateCompute::isEqual(const MEDFileUMeshAggregateCompute& other, double eps, std::string& what) const
{
  const MEDCouplingUMesh *m1(getUmesh());
  const MEDCouplingUMesh *m2(other.getUmesh());
  if((m1==nullptr && m2!=nullptr) || (m1!=nullptr && m2==nullptr))
    {
      what="Presence of mesh in one sublevel and not in other!";
      return false;
    }
  if(m1)
    {
      std::string what2;
      if(!m1->isEqualIfNotWhy(m2,eps,what2))
        {
          what=std::string("meshes at a sublevel are not deeply equal (")+what2+std::string(")!");
          return false;
        }
    }
  std::size_t const sz(_part_def.size());
  if(sz!=other._part_def.size())
    {
      what=std::string("number of subdivision per geo type for part definition is not the same !");
      return false;
    }
  for(std::size_t i=0;i<sz;i++)
    {
      const PartDefinition *pd0(_part_def[i]),*pd1(other._part_def[i]);
      if(!pd0 && !pd1)
        continue;
      if((!pd0 && pd1) || (pd0 && !pd1))
        {
          what=std::string("a cell part def is defined only for one among this or other !");
          return false;
        }
      bool const ret(pd0->isEqual(pd1,what));
      if(!ret)
        return false;
    }
  return true;
}

void MEDFileUMeshAggregateCompute::checkConsistency() const
{
  if(_mp_time >= _m_time)
    for(const auto & _m_part : _m_parts)
      _m_part->checkConsistency();
  else
    _m->checkConsistency();
}

void MEDFileUMeshAggregateCompute::clearNonDiscrAttributes() const
{
  for(const auto & _m_part : _m_parts)
    MEDFileUMeshSplitL1::ClearNonDiscrAttributes(_m_part);
  MEDFileUMeshSplitL1::ClearNonDiscrAttributes(_m);
}

void MEDFileUMeshAggregateCompute::synchronizeTinyInfo(const MEDFileMesh& master) const
{
  for(const auto & _m_part : _m_parts)
    {
      const MEDCoupling1GTUMesh *tmp(_m_part);
      if(tmp)
        {
          (const_cast<MEDCoupling1GTUMesh *>(tmp))->setName(master.getName().c_str());
          (const_cast<MEDCoupling1GTUMesh *>(tmp))->setDescription(master.getDescription().c_str());
          (const_cast<MEDCoupling1GTUMesh *>(tmp))->setTime(master.getTimeValue(),master.getIteration(),master.getOrder());
          (const_cast<MEDCoupling1GTUMesh *>(tmp))->setTimeUnit(master.getTimeUnit());
        }
    }
  const MEDCouplingUMesh *m(_m);
  if(m)
    {
      (const_cast<MEDCouplingUMesh *>(m))->setName(master.getName().c_str());
      (const_cast<MEDCouplingUMesh *>(m))->setDescription(master.getDescription().c_str());
      (const_cast<MEDCouplingUMesh *>(m))->setTime(master.getTimeValue(),master.getIteration(),master.getOrder());
      (const_cast<MEDCouplingUMesh *>(m))->setTimeUnit(master.getTimeUnit());
    }
}

bool MEDFileUMeshAggregateCompute::empty() const
{
  if(_mp_time<_m_time)
    return ((const MEDCouplingUMesh *)_m)==nullptr;
  //else _mp_time>=_m_time)
  return _m_parts.empty();
}

int MEDFileUMeshAggregateCompute::getMeshDimension() const
{
  if(_mp_time<_m_time)
    {
      const MEDCouplingUMesh *m(_m);
      if(!m)
        throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::getMeshDimension : no umesh in this !");
      return m->getMeshDimension();
    }
  else
    {
      if(_m_parts.empty())
        throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::getMeshDimension : part mesh is empty !");
      const MEDCoupling1GTUMesh *m(_m_parts[0]);
      if(!m)
        throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::getMeshDimension : part mesh contains null instance !");
      return m->getMeshDimension();
    }
}

std::vector<mcIdType> MEDFileUMeshAggregateCompute::getDistributionOfTypes() const
{
  if(_mp_time<_m_time)
    {
      const MEDCouplingUMesh *m(_m);
      if(!m)
        throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::getDistributionOfTypes : no umesh in this !");
      return m->getDistributionOfTypes();
    }
  else
    {
      std::vector<mcIdType> ret;
      for(const auto & _m_part : _m_parts)
        {
          const MEDCoupling1GTUMesh *tmp(_m_part);
          if(!tmp)
            throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::getDistributionOfTypes : part mesh contains null instance !");
          std::vector<mcIdType> ret0(tmp->getDistributionOfTypes());
          ret.insert(ret.end(),ret0.begin(),ret0.end());
        }
      return ret;
    }
}

mcIdType MEDFileUMeshAggregateCompute::getSize() const
{
  if(_mp_time<_m_time)
    {
      const MEDCouplingUMesh *m(_m);
      if(!m)
        throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::getSize : no umesh in this !");
      return m->getNumberOfCells();
    }
  else
    {
      mcIdType ret=0;
      for(const auto & _m_part : _m_parts)
        {
          const MEDCoupling1GTUMesh *m(_m_part);
          if(!m)
            throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::getSize : part mesh contains null instance !");
          ret+=m->getNumberOfCells();
        }
      return ret;
    }
}

void MEDFileUMeshAggregateCompute::setCoords(DataArrayDouble *coords)
{
  for(auto & _m_part : _m_parts)
    {
      MEDCoupling1GTUMesh *tmp(_m_part);
      if(tmp)
        _m_part->setCoords(coords);
    }
  MEDCouplingUMesh *m(_m);
  if(m)
    m->setCoords(coords);
}

MEDFileEltStruct4Mesh *MEDFileEltStruct4Mesh::New(med_idt fid, const std::string& mName, int dt, int it, int iterOnStEltOfMesh, MEDFileMeshReadSelector *mrs)
{
  return new MEDFileEltStruct4Mesh(fid,mName,dt,it,iterOnStEltOfMesh,mrs);
}

std::size_t MEDFileEltStruct4Mesh::getHeapMemorySizeWithoutChildren() const
{
  return _geo_type_name.capacity()+_vars.capacity()*sizeof(MCAuto<DataArray>);
}

std::vector<const MEDCoupling::BigMemoryObject*> MEDFileEltStruct4Mesh::getDirectChildrenWithNull() const
{
  std::vector<const MEDCoupling::BigMemoryObject*> ret;
  ret.push_back(_conn);
  ret.push_back(_common);
  for(const auto & _var : _vars)
    ret.push_back(_var);
  return ret;
}

MEDFileEltStruct4Mesh::MEDFileEltStruct4Mesh(med_idt fid, const std::string& mName, int dt, int it, int iterOnStEltOfMesh, MEDFileMeshReadSelector *mrs)
{
  med_geometry_type geoType;
  INTERP_KERNEL::AutoPtr<char> geoTypeName(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  MEDFILESAFECALLERRD0(MEDmeshEntityInfo,(fid,mName.c_str(),dt,it,MED_STRUCT_ELEMENT,iterOnStEltOfMesh+1,geoTypeName,&geoType));
  _geo_type=geoType;
  _geo_type_name=MEDLoaderBase::buildStringFromFortran(geoTypeName,MED_NAME_SIZE);
  mcIdType nCells(0);
  {
    med_bool chgt=MED_FALSE,trsf=MED_FALSE;
    nCells=MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_STRUCT_ELEMENT,geoType,MED_CONNECTIVITY,MED_NODAL,&chgt,&trsf);
  }
  MCAuto<MEDFileMeshSupports> mss(MEDFileMeshSupports::New(fid));
  MCAuto<MEDFileStructureElements> mse(MEDFileStructureElements::New(fid,mss));
  mcIdType const nbEntities(mse->getNumberOfNodesPerSE(_geo_type_name));
  MCAuto<DataArrayMedInt> miConn=DataArrayMedInt::New(); miConn->alloc(nCells*nbEntities);
  MEDFILESAFECALLERRD0(MEDmeshElementConnectivityRd,(fid,mName.c_str(),dt,it,MED_STRUCT_ELEMENT,_geo_type,MED_NODAL,MED_FULL_INTERLACE,miConn->getPointer()));
  _conn=FromMedIntArray<mcIdType>(miConn);
  _conn->applyLin(1,-1);
  _conn->rearrange(nbEntities);
  _common=MEDFileUMeshPerTypeCommon::New();
  _common->loadCommonPart(fid,mName.c_str(),dt,it,nCells,geoType,MED_STRUCT_ELEMENT,mrs);
  std::vector<std::string> vns(mse->getVarAttsOf(_geo_type_name));
  std::size_t const sz(vns.size());
  _vars.resize(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      const MEDFileSEVarAtt *var(mse->getVarAttOf(_geo_type_name,vns[i]));
      MCAuto<DataArray> gen(var->getGenerator());
      MCAuto<DataArray> arr(gen->buildNewEmptyInstance());
      arr->alloc(nCells,var->getNbOfComponents());
      arr->setName(vns[i]);
      MEDFILESAFECALLERRD0(MEDmeshStructElementVarAttRd,(fid,mName.c_str(),dt,it,_geo_type,vns[i].c_str(),arr->getVoidStarPointer()));
      _vars[i]=arr;
    }
}
