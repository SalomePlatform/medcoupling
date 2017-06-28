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

#include "MEDFileMeshLL.hxx"
#include "MEDFileMesh.hxx"
#include "MEDLoaderBase.hxx"
#include "MEDFileSafeCaller.txx"
#include "MEDFileMeshReadSelector.hxx"
#include "MEDFileStructureElement.hxx"
#include "MEDFileMeshSupport.hxx"

#include "MEDCouplingUMesh.hxx"

#include "InterpKernelAutoPtr.hxx"
#include "CellModel.hxx"

#include <set>
#include <iomanip>

extern med_geometry_type typmai[MED_N_CELL_FIXED_GEO];
extern INTERP_KERNEL::NormalizedCellType typmai2[MED_N_CELL_FIXED_GEO];
extern med_geometry_type typmainoeud[1];

using namespace MEDCoupling;

const char MEDFileMeshL2::ZE_SEP_FOR_FAMILY_KILLERS[]="!/__\\!";//important start by - because ord('!')==33 the smallest (!=' ') to preserve orders at most.

int MEDFileMeshL2::ZE_SEP2_FOR_FAMILY_KILLERS=4;

std::vector<std::string> MeshCls::getAxisInfoOnMesh(med_idt fid, const std::string& mName, MEDCoupling::MEDCouplingMeshType& meshType, MEDCoupling::MEDCouplingAxisType& axType, int& nstep, int& Mdim, MEDFileString& description, MEDFileString& dtunit, MEDFileString& univName) const
{
  med_mesh_type type_maillage;
  med_int spaceDim;
  med_sorting_type stype;
  med_axis_type axistype;
  int naxis(MEDmeshnAxis(fid,getID()));
  INTERP_KERNEL::AutoPtr<char> nameTmp(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> axisname(MEDLoaderBase::buildEmptyString(naxis*MED_SNAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> axisunit(MEDLoaderBase::buildEmptyString(naxis*MED_SNAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> univTmp(MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE));
  if(MEDmeshInfo(fid,getID(),nameTmp,&spaceDim,&Mdim,&type_maillage,description.getPointer(),dtunit.getPointer(),
      &stype,&nstep,&axistype,axisname,axisunit)!=0)
    throw INTERP_KERNEL::Exception("A problem has been detected when trying to get info on mesh !");
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
      std::string info(MEDLoaderBase::buildUnionUnit(((char *)axisname)+i*MED_SNAME_SIZE,MED_SNAME_SIZE,((char *)axisunit)+i*MED_SNAME_SIZE,MED_SNAME_SIZE));
      infosOnComp[i]=info;
    }
  return infosOnComp;
}

double MeshCls::checkMeshTimeStep(med_idt fid, const std::string& mName, int nstep, int dt, int it) const
{
  bool found=false;
  med_int numdt,numit;
  med_float dtt;
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

std::vector<std::string> StructMeshCls::getAxisInfoOnMesh(med_idt fid, const std::string& mName, MEDCoupling::MEDCouplingMeshType& meshType, MEDCoupling::MEDCouplingAxisType& axType, int& nstep, int& Mdim, MEDFileString& description, MEDFileString& dtunit, MEDFileString& univName) const
{
  INTERP_KERNEL::AutoPtr<char> msn(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> zeDescription(MEDLoaderBase::buildEmptyString(MED_COMMENT_SIZE));
  med_axis_type medAxType;
  int nAxis(MEDsupportMeshnAxis(fid,getID()));
  INTERP_KERNEL::AutoPtr<char> axisName(new char[MED_SNAME_SIZE*nAxis+1]),axisUnit(new char[MED_SNAME_SIZE*nAxis+1]);
  int spaceDim(0),meshDim(0);
  MEDFILESAFECALLERRD0(MEDsupportMeshInfo,(fid,getID(),msn,&spaceDim,&meshDim,zeDescription,&medAxType,axisName,axisUnit));
  std::string descriptionCpp(MEDLoaderBase::buildStringFromFortran(zeDescription,MED_COMMENT_SIZE));
  description.set(descriptionCpp.c_str());
  dtunit.clear(); univName.clear(); meshType=UNSTRUCTURED; nstep=1;
  axType=MEDFileMeshL2::TraduceAxisType(medAxType);
  //int nmodels(0);
  //med_bool chgt=MED_FALSE,trsf=MED_FALSE;
  //nmodels=MEDmeshnEntity(fid,_name.c_str(),MED_NO_DT,MED_NO_IT,MED_STRUCT_ELEMENT,MED_GEO_ALL,MED_CONNECTIVITY,MED_NODAL,&chgt,&trsf);
  std::vector<std::string> ret;
  for(int i=0;i<nAxis;i++)
    {
      std::string info(DataArray::BuildInfoFromVarAndUnit(MEDLoaderBase::buildStringFromFortran(axisName+i*MED_SNAME_SIZE,MED_SNAME_SIZE),
                                                          MEDLoaderBase::buildStringFromFortran(axisUnit+i*MED_SNAME_SIZE,MED_SNAME_SIZE)));
      ret.push_back(info);
    }
  return ret;
}

double StructMeshCls::checkMeshTimeStep(med_idt fid, const std::string& mName, int nstep, int dt, int it) const
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
  med_mesh_type type_maillage;
  char maillage_description[MED_COMMENT_SIZE+1];
  char dtunit[MED_LNAME_SIZE+1];
  med_int spaceDim,dim;
  char nommaa[MED_NAME_SIZE+1];
  med_int n=MEDnMesh(fid);
  char found(0);
  int ret=-1;
  med_sorting_type stype;
  std::vector<std::string> ms;
  int nstep;
  med_axis_type axistype;
  for(int i=0;i<n && found==0;i++)
    {
      int naxis(MEDmeshnAxis(fid,i+1));
      INTERP_KERNEL::AutoPtr<char> axisname(MEDLoaderBase::buildEmptyString(naxis*MED_SNAME_SIZE)),axisunit(MEDLoaderBase::buildEmptyString(naxis*MED_SNAME_SIZE));
      MEDFILESAFECALLERRD0(MEDmeshInfo,(fid,i+1,nommaa,&spaceDim,&dim,&type_maillage,maillage_description,dtunit,&stype,&nstep,&axistype,axisname,axisunit));      
      dtunit1=MEDLoaderBase::buildStringFromFortran(dtunit,sizeof(dtunit));
      std::string cur(MEDLoaderBase::buildStringFromFortran(nommaa,sizeof(nommaa)));
      ms.push_back(cur);
      if(cur==mName)
        {
          found=1;
          ret=i+1;
        }
    }
  if(found==0)
    {//last chance ! Is it a support mesh ?
      int nbSM(MEDnSupportMesh(fid));
      for(int i=0;i<nbSM && found==0;i++)
        {
          int naxis(MEDsupportMeshnAxis(fid,i+1));
          INTERP_KERNEL::AutoPtr<char> axisname(MEDLoaderBase::buildEmptyString(naxis*MED_SNAME_SIZE)),axisunit(MEDLoaderBase::buildEmptyString(naxis*MED_SNAME_SIZE));
          MEDFILESAFECALLERRD0(MEDsupportMeshInfo,(fid,i+1,nommaa,&spaceDim,&dim,maillage_description,&axistype,axisname,axisunit));
          std::string cur(MEDLoaderBase::buildStringFromFortran(nommaa,sizeof(nommaa)));
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
        dt=numdt; it=numit;
        return new MeshCls(ret);
      }
    case 2:
      {
        meshType=UNSTRUCTURED;
        dt=MED_NO_DT; it=MED_NO_IT; dtunit1.clear();
        axType=axType=TraduceAxisType(axistype);
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

void MEDFileMeshL2::ReadFamiliesAndGrps(med_idt fid, const std::string& meshName, std::map<std::string,int>& fams, std::map<std::string, std::vector<std::string> >& grps, MEDFileMeshReadSelector *mrs)
{
  if(mrs && !(mrs->isCellFamilyFieldReading() || mrs->isNodeFamilyFieldReading()))
    return ;
  char nomfam[MED_NAME_SIZE+1];
  med_int numfam;
  int nfam=MEDnFamily(fid,meshName.c_str());
  std::vector< std::pair<std::string,std::pair<int,std::vector<std::string> > > > crudeFams(nfam);
  for(int i=0;i<nfam;i++)
    {
      int ngro=MEDnFamilyGroup(fid,meshName.c_str(),i+1);
      med_int natt=MEDnFamily23Attribute(fid,meshName.c_str(),i+1);
      INTERP_KERNEL::AutoPtr<med_int> attide=new med_int[natt];
      INTERP_KERNEL::AutoPtr<med_int> attval=new med_int[natt];
      INTERP_KERNEL::AutoPtr<char> attdes=new char[MED_COMMENT_SIZE*natt+1];
      INTERP_KERNEL::AutoPtr<char> gro=new char[MED_LNAME_SIZE*ngro+1];
      MEDfamily23Info(fid,meshName.c_str(),i+1,nomfam,attide,attval,attdes,&numfam,gro);
      std::string famName(MEDLoaderBase::buildStringFromFortran(nomfam,MED_NAME_SIZE));
      std::vector<std::string> grps(ngro);
      for(int j=0;j<ngro;j++)
        grps[j]=MEDLoaderBase::buildStringFromFortran(gro+j*MED_LNAME_SIZE,MED_LNAME_SIZE);
      crudeFams[i]=std::pair<std::string,std::pair<int,std::vector<std::string> > >(famName,std::pair<int,std::vector<std::string> >(numfam,grps));
    }
  RenameFamiliesFromFileToMemInternal(crudeFams);
  for(std::vector< std::pair<std::string,std::pair<int,std::vector<std::string> > > >::const_iterator it0=crudeFams.begin();it0!=crudeFams.end();it0++)
    {
      fams[(*it0).first]=(*it0).second.first;
      for(std::vector<std::string>::const_iterator it1=(*it0).second.second.begin();it1!=(*it0).second.second.end();it1++)
        grps[*it1].push_back((*it0).first);
    }
}

void MEDFileMeshL2::WriteFamiliesAndGrps(med_idt fid, const std::string& mname, const std::map<std::string,int>& fams, const std::map<std::string, std::vector<std::string> >& grps, int tooLongStrPol)
{
  std::vector< std::pair<std::string,std::pair<int,std::vector<std::string> > > > crudeFams(fams.size());
  std::size_t ii(0);
  for(std::map<std::string,int>::const_iterator it=fams.begin();it!=fams.end();it++,ii++)
    {
      std::vector<std::string> grpsOfFam;
      for(std::map<std::string, std::vector<std::string> >::const_iterator it1=grps.begin();it1!=grps.end();it1++)
        {
          if(std::find((*it1).second.begin(),(*it1).second.end(),(*it).first)!=(*it1).second.end())
            grpsOfFam.push_back((*it1).first);
        }
      crudeFams[ii]=std::pair<std::string,std::pair<int,std::vector<std::string> > >((*it).first,std::pair<int,std::vector<std::string> >((*it).second,grpsOfFam));
    }
  RenameFamiliesFromMemToFileInternal(crudeFams);
  for(std::vector< std::pair<std::string,std::pair<int,std::vector<std::string> > > >::const_iterator it=crudeFams.begin();it!=crudeFams.end();it++)
    {
      int ngro((*it).second.second.size());
      INTERP_KERNEL::AutoPtr<char> groName=MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE*ngro);
      int i=0;
      for(std::vector<std::string>::const_iterator it2=(*it).second.second.begin();it2!=(*it).second.second.end();it2++,i++)
        MEDLoaderBase::safeStrCpy2((*it2).c_str(),MED_LNAME_SIZE,groName+i*MED_LNAME_SIZE,tooLongStrPol);
      INTERP_KERNEL::AutoPtr<char> famName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
      MEDLoaderBase::safeStrCpy((*it).first.c_str(),MED_NAME_SIZE,famName,tooLongStrPol);
      int ret=MEDfamilyCr(fid,mname.c_str(),famName,(*it).second.first,ngro,groName);
      ret++;
    }
}

void MEDFileMeshL2::RenameFamiliesPatternInternal(std::vector< std::pair<std::string,std::pair<int,std::vector<std::string> > > >& crudeFams, RenameFamiliesPatternFunc func)
{
  std::size_t ii(0);
  std::vector<std::string> fams(crudeFams.size());
  for(std::vector< std::pair<std::string,std::pair<int,std::vector<std::string> > > >::const_iterator it=crudeFams.begin();it!=crudeFams.end();it++,ii++)
    fams[ii]=(*it).first;
  if(!func(fams))
    return ;
  ii=0;
  for(std::vector< std::pair<std::string,std::pair<int,std::vector<std::string> > > >::iterator it=crudeFams.begin();it!=crudeFams.end();it++,ii++)
    (*it).first=fams[ii];
}

/*!
 * This method is dedicated to the killers that use a same family name to store different family ids. MED file API authorizes it.
 * So this method renames families (if needed generaly not !) in order to have a discriminant name for families.
 */
void MEDFileMeshL2::RenameFamiliesFromFileToMemInternal(std::vector< std::pair<std::string,std::pair<int,std::vector<std::string> > > >& crudeFams)
{
  RenameFamiliesPatternInternal(crudeFams,RenameFamiliesFromFileToMem);
}

bool MEDFileMeshL2::RenameFamiliesFromFileToMem(std::vector< std::string >& famNames)
{
  std::map<std::string,int> m;
  std::set<std::string> s;
  for(std::vector< std::string >::const_iterator it=famNames.begin();it!=famNames.end();it++)
    {
      if(s.find(*it)!=s.end())
        m[*it]=0;
      s.insert(*it);
    }
  if(m.empty())
    return false;// the general case !
  for(std::vector< std::string >::iterator it=famNames.begin();it!=famNames.end();it++)
    {
      std::map<std::string,int>::iterator it2(m.find(*it));
      if(it2!=m.end())
        {
          std::ostringstream oss; oss << *it << ZE_SEP_FOR_FAMILY_KILLERS << std::setfill('0') << std::setw(ZE_SEP2_FOR_FAMILY_KILLERS) << (*it2).second++;
          *it=oss.str();
        }
    }
  return true;
}

/*!
 * This method is dedicated to the killers that use a same family name to store different family ids. MED file API authorizes it.
 * So this method renames families (if needed generaly not !) in order to have a discriminant name for families.
 */
void MEDFileMeshL2::RenameFamiliesFromMemToFileInternal(std::vector< std::pair<std::string,std::pair<int,std::vector<std::string> > > >& crudeFams)
{
  RenameFamiliesPatternInternal(crudeFams,RenameFamiliesFromMemToFile);
}

bool MEDFileMeshL2::RenameFamiliesFromMemToFile(std::vector< std::string >& famNames)
{
  bool isSmthingStrange(false);
  for(std::vector< std::string >::const_iterator it=famNames.begin();it!=famNames.end();it++)
    {
      std::size_t found((*it).find(ZE_SEP_FOR_FAMILY_KILLERS));
      if(found!=std::string::npos)
        isSmthingStrange=true;
    }
  if(!isSmthingStrange)
    return false;
  // pattern matching
  std::map< std::string, std::vector<std::string> > m;
  for(std::vector< std::string >::const_iterator it=famNames.begin();it!=famNames.end();it++)
    {
      std::size_t found((*it).find(ZE_SEP_FOR_FAMILY_KILLERS));
      if(found!=std::string::npos && found>=1)
        {
          std::string s1((*it).substr(found+sizeof(ZE_SEP_FOR_FAMILY_KILLERS)-1));
          if(s1.size()!=ZE_SEP2_FOR_FAMILY_KILLERS)
            continue;
          int k(-1);
          std::istringstream iss(s1);
          iss >> k;
          bool isOK((iss.rdstate() & ( std::istream::failbit | std::istream::eofbit)) == std::istream::eofbit);
          if(isOK && k>=0)
            {
              std::string s0((*it).substr(0,found));
              m[s0].push_back(*it);
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
      for(std::vector<std::string>::const_iterator it1=(*it).second.begin();it1!=(*it).second.end();it1++)
        zeMap[*it1]=(*it).first;
    }
  if(zeMap.empty())
    return false;
  // traduce
  for(std::vector< std::string >::iterator it=famNames.begin();it!=famNames.end();it++)
    {
      std::map<std::string,std::string>::iterator it1(zeMap.find(*it));
      if(it1!=zeMap.end())
        *it=(*it1).second;
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
{
}

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
  std::vector<std::string> infosOnComp(loadCommonPart(fid,mId,mName,dt,it,Mdim));
  if(Mdim==-4)
    return ;
  loadConnectivity(fid,Mdim,mName,dt,it,mrs);//to improve check (dt,it) coherency
  loadCoords(fid,infosOnComp,mName,dt,it);
}

void MEDFileUMeshL2::loadPart(med_idt fid, const MeshOrStructMeshCls *mId, const std::string& mName, const std::vector<INTERP_KERNEL::NormalizedCellType>& types, const std::vector<int>& slicPerTyp, int dt, int it, MEDFileMeshReadSelector *mrs)
{
  int Mdim;
  std::vector<std::string> infosOnComp(loadCommonPart(fid,mId,mName,dt,it,Mdim));
  if(Mdim==-4)
    return ;
  loadPartOfConnectivity(fid,Mdim,mName,types,slicPerTyp,dt,it,mrs);
  med_bool changement,transformation;
  int nCoords(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NONE,MED_COORDINATE,MED_NO_CMODE,&changement,&transformation));
  std::vector<bool> fetchedNodeIds(nCoords,false);
  for(std::vector< std::vector< MCAuto<MEDFileUMeshPerType> > >::const_iterator it0=_per_type_mesh.begin();it0!=_per_type_mesh.end();it0++)
    for(std::vector< MCAuto<MEDFileUMeshPerType> >::const_iterator it1=(*it0).begin();it1!=(*it0).end();it1++)
      (*it1)->getMesh()->computeNodeIdsAlg(fetchedNodeIds);
  int nMin(std::distance(fetchedNodeIds.begin(),std::find(fetchedNodeIds.begin(),fetchedNodeIds.end(),true)));
  int nMax(std::distance(fetchedNodeIds.rbegin(),std::find(fetchedNodeIds.rbegin(),fetchedNodeIds.rend(),true)));
  nMax=nCoords-nMax;
  for(std::vector< std::vector< MCAuto<MEDFileUMeshPerType> > >::const_iterator it0=_per_type_mesh.begin();it0!=_per_type_mesh.end();it0++)
    for(std::vector< MCAuto<MEDFileUMeshPerType> >::const_iterator it1=(*it0).begin();it1!=(*it0).end();it1++)
      (*it1)->getMesh()->renumberNodesWithOffsetInConn(-nMin);
  loadPartCoords(fid,infosOnComp,mName,dt,it,nMin,nMax);
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

void MEDFileUMeshL2::loadPartOfConnectivity(med_idt fid, int mdim, const std::string& mName, const std::vector<INTERP_KERNEL::NormalizedCellType>& types, const std::vector<int>& slicPerTyp, int dt, int it, MEDFileMeshReadSelector *mrs)
{
  std::size_t nbOfTypes(types.size());
  if(slicPerTyp.size()!=3*nbOfTypes)
    throw INTERP_KERNEL::Exception("MEDFileUMeshL2::loadPartOfConnectivity : The size of slicPerTyp array is expected to be equal to 3 times size of array types !");
  std::set<INTERP_KERNEL::NormalizedCellType> types2(types.begin(),types.end());
  if(types2.size()!=nbOfTypes)
    throw INTERP_KERNEL::Exception("MEDFileUMeshL2::loadPartOfConnectivity : the geometric types in types array must appear once !");
  _per_type_mesh.resize(1);
  _per_type_mesh[0].clear();
  for(std::size_t ii=0;ii<nbOfTypes;ii++)
    {
      int strt(slicPerTyp[3*ii+0]),stp(slicPerTyp[3*ii+1]),step(slicPerTyp[3*ii+2]);
      MCAuto<MEDFileUMeshPerType> tmp(MEDFileUMeshPerType::NewPart(fid,mName.c_str(),dt,it,mdim,types[ii],strt,stp,step,mrs));
      _per_type_mesh[0].push_back(tmp);
    }
  sortTypes();
}

void MEDFileUMeshL2::loadCoords(med_idt fid, const std::vector<std::string>& infosOnComp, const std::string& mName, int dt, int it)
{
  int spaceDim((int)infosOnComp.size());
  med_bool changement,transformation;
  int nCoords(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NONE,MED_COORDINATE,MED_NO_CMODE,&changement,&transformation));
  _coords=DataArrayDouble::New();
  _coords->alloc(nCoords,spaceDim);
  double *coordsPtr(_coords->getPointer());
  if (nCoords)
    MEDFILESAFECALLERRD0(MEDmeshNodeCoordinateRd,(fid,mName.c_str(),dt,it,MED_FULL_INTERLACE,coordsPtr));
  if(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,MED_FAMILY_NUMBER,MED_NODAL,&changement,&transformation)>0)
    {
      _fam_coords=DataArrayInt::New();
      _fam_coords->alloc(nCoords,1);
      MEDFILESAFECALLERRD0(MEDmeshEntityFamilyNumberRd,(fid,mName.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,_fam_coords->getPointer()));
    }
  else
    _fam_coords=0;
  if(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,MED_NUMBER,MED_NODAL,&changement,&transformation)>0)
    {
      _num_coords=DataArrayInt::New();
      _num_coords->alloc(nCoords,1);
      MEDFILESAFECALLERRD0(MEDmeshEntityNumberRd,(fid,mName.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,_num_coords->getPointer()));
    }
  else
    _num_coords=0;
  if(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,MED_NAME,MED_NODAL,&changement,&transformation)>0)
    {
      _name_coords=DataArrayAsciiChar::New();
      _name_coords->alloc(nCoords+1,MED_SNAME_SIZE);//not a bug to avoid the memory corruption due to last \0 at the end
      MEDFILESAFECALLERRD0(MEDmeshEntityNameRd,(fid,mName.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,_name_coords->getPointer()));
      _name_coords->reAlloc(nCoords);//not a bug to avoid the memory corruption due to last \0 at the end
    }
  else
    _name_coords=0;
  if(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,MED_GLOBAL_NUMBER,MED_NODAL,&changement,&transformation)>0)
    {
      _global_num_coords=DataArrayInt::New();
      _global_num_coords->alloc(nCoords,1);
      MEDFILESAFECALLERRD0(MEDmeshGlobalNumberRd,(fid,mName.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,_global_num_coords->getPointer()));
    }
  for(int i=0;i<spaceDim;i++)
    _coords->setInfoOnComponent(i,infosOnComp[i]);
}

void MEDFileUMeshL2::loadPartCoords(med_idt fid, const std::vector<std::string>& infosOnComp, const std::string& mName, int dt, int it, int nMin, int nMax)
{
  med_bool changement,transformation;
  int spaceDim((int)infosOnComp.size()),nCoords(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NONE,MED_COORDINATE,MED_NO_CMODE,&changement,&transformation));
  _coords=DataArrayDouble::New();
  int nbNodesToLoad(nMax-nMin);
  _coords->alloc(nbNodesToLoad,spaceDim);
  med_filter filter=MED_FILTER_INIT,filter2=MED_FILTER_INIT;
  MEDfilterBlockOfEntityCr(fid,/*nentity*/nCoords,/*nvaluesperentity*/1,/*nconstituentpervalue*/spaceDim,
                           MED_ALL_CONSTITUENT,MED_FULL_INTERLACE,MED_COMPACT_STMODE,MED_NO_PROFILE,
                           /*start*/nMin+1,/*stride*/1,/*count*/1,/*blocksize*/nbNodesToLoad,
                           /*lastblocksize=useless because count=1*/0,&filter);
  MEDFILESAFECALLERRD0(MEDmeshNodeCoordinateAdvancedRd,(fid,mName.c_str(),dt,it,&filter,_coords->getPointer()));
  _part_coords=PartDefinition::New(nMin,nMax,1);
  MEDfilterClose(&filter);
  MEDfilterBlockOfEntityCr(fid,nCoords,1,1,MED_ALL_CONSTITUENT,MED_FULL_INTERLACE,MED_COMPACT_STMODE,
                           MED_NO_PROFILE,nMin+1,1,1,nbNodesToLoad,0,&filter2);
  if(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,MED_FAMILY_NUMBER,MED_NODAL,&changement,&transformation)>0)
    {
      _fam_coords=DataArrayInt::New();
      _fam_coords->alloc(nbNodesToLoad,1);
      MEDFILESAFECALLERRD0(MEDmeshEntityAttributeAdvancedRd,(fid,mName.c_str(),MED_FAMILY_NUMBER,dt,it,MED_NODE,MED_NO_GEOTYPE,&filter2,_fam_coords->getPointer()));
    }
  else
    _fam_coords=0;
  if(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,MED_NUMBER,MED_NODAL,&changement,&transformation)>0)
    {
      _num_coords=DataArrayInt::New();
      _num_coords->alloc(nbNodesToLoad,1);
      MEDFILESAFECALLERRD0(MEDmeshEntityAttributeAdvancedRd,(fid,mName.c_str(),MED_NUMBER,dt,it,MED_NODE,MED_NO_GEOTYPE,&filter2,_num_coords->getPointer()));
    }
  else
    _num_coords=0;
  if(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,MED_NAME,MED_NODAL,&changement,&transformation)>0)
    {
      _name_coords=DataArrayAsciiChar::New();
      _name_coords->alloc(nbNodesToLoad+1,MED_SNAME_SIZE);//not a bug to avoid the memory corruption due to last \0 at the end
      MEDFILESAFECALLERRD0(MEDmeshEntityAttributeAdvancedRd,(fid,mName.c_str(),MED_NAME,dt,it,MED_NODE,MED_NO_GEOTYPE,&filter2,_name_coords->getPointer()));
      _name_coords->reAlloc(nbNodesToLoad);//not a bug to avoid the memory corruption due to last \0 at the end
    }
  else
    _name_coords=0;
  MEDfilterClose(&filter2);
  _coords->setInfoOnComponents(infosOnComp);
}

void MEDFileUMeshL2::sortTypes()
{
  std::set<int> mdims;
  std::vector< MCAuto<MEDFileUMeshPerType> > tmp(_per_type_mesh[0]);
  _per_type_mesh.clear();
  for(std::vector< MCAuto<MEDFileUMeshPerType> >::const_iterator it=tmp.begin();it!=tmp.end();it++)
    mdims.insert((*it)->getDim());
  if(mdims.empty())
    return;
  int mdim=*mdims.rbegin();
  _per_type_mesh.resize(mdim+1);
  for(int dim=mdim+1;dim!=0;dim--)
    {
      std::vector< MCAuto<MEDFileUMeshPerType> >& elt=_per_type_mesh[mdim+1-dim];
      for(std::vector< MCAuto<MEDFileUMeshPerType> >::const_iterator it=tmp.begin();it!=tmp.end();it++)
        if((*it)->getDim()==dim-1)
          elt.push_back(*it);
    }
  // suppression of contiguous empty levels at the end of _per_type_mesh.
  int nbOfUselessLev=0;
  bool isFirst=true;
  for(std::vector< std::vector< MCAuto<MEDFileUMeshPerType> > >::reverse_iterator it2=_per_type_mesh.rbegin();it2!=_per_type_mesh.rend();it2++)
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

void MEDFileUMeshL2::WriteCoords(med_idt fid, const std::string& mname, int dt, int it, double time, const DataArrayDouble *coords, const DataArrayInt *famCoords, const DataArrayInt *numCoords, const DataArrayAsciiChar *nameCoords, const DataArrayInt *globalNumCoords)
{
  if(!coords)
    return ;
  MEDFILESAFECALLERWR0(MEDmeshNodeCoordinateWr,(fid,mname.c_str(),dt,it,time,MED_FULL_INTERLACE,coords->getNumberOfTuples(),coords->begin()));
  if(famCoords)
    MEDFILESAFECALLERWR0(MEDmeshEntityFamilyNumberWr,(fid,mname.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,famCoords->getNumberOfTuples(),famCoords->begin()));
  if(numCoords)
    MEDFILESAFECALLERWR0(MEDmeshEntityNumberWr,(fid,mname.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,numCoords->getNumberOfTuples(),numCoords->begin()));
  if(nameCoords)
    {
      if(nameCoords->getNumberOfComponents()!=MED_SNAME_SIZE)
        {
          std::ostringstream oss; oss << " MEDFileUMeshL2::WriteCoords : expected a name field on nodes with number of components set to " << MED_SNAME_SIZE;
          oss << " ! The array has " << nameCoords->getNumberOfComponents() << " components !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      MEDFILESAFECALLERWR0(MEDmeshEntityNameWr,(fid,mname.c_str(),dt,it,MED_NODE,MED_NO_GEOTYPE,nameCoords->getNumberOfTuples(),nameCoords->begin()));
    }
  if(globalNumCoords)
    MEDFILESAFECALLERWR0(MEDmeshGlobalNumberWr,(fid,mname.c_str(),dt,it,MED_NODE,MED_NONE,globalNumCoords->getNumberOfTuples(),globalNumCoords->begin()));
}

bool MEDFileUMeshL2::isFamDefinedOnLev(int levId) const
{
  for(std::vector< MCAuto<MEDFileUMeshPerType> >::const_iterator it=_per_type_mesh[levId].begin();it!=_per_type_mesh[levId].end();it++)
    if((*it)->getFam()==0)
      return false;
  return true;
}

bool MEDFileUMeshL2::isNumDefinedOnLev(int levId) const
{
  for(std::vector< MCAuto<MEDFileUMeshPerType> >::const_iterator it=_per_type_mesh[levId].begin();it!=_per_type_mesh[levId].end();it++)
    if((*it)->getNum()==0)
      return false;
  return true;
}

bool MEDFileUMeshL2::isNamesDefinedOnLev(int levId) const
{
  for(std::vector< MCAuto<MEDFileUMeshPerType> >::const_iterator it=_per_type_mesh[levId].begin();it!=_per_type_mesh[levId].end();it++)
    if((*it)->getNames()==0)
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
      med_data_type dataTypeReq=GetDataTypeCorrespondingToSpaceId(i);
      med_bool chgt=MED_FALSE,trsf=MED_FALSE;
      int nbOfElt(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NONE,dataTypeReq,MED_NO_CMODE,&chgt,&trsf));
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
{
}

void MEDFileCLMeshL2::loadAll(med_idt fid, const MeshOrStructMeshCls *mId, const std::string& mName, int dt, int it)
{
  _name.set(mName.c_str());
  int nstep;
  int Mdim;
  MEDCoupling::MEDCouplingMeshType meshType;
  MEDCoupling::MEDCouplingAxisType dummy3;
  std::vector<std::string> infosOnComp(getAxisInfoOnMesh(fid,mId,mName,meshType,dummy3,nstep,Mdim));
  if(meshType!=CURVE_LINEAR)
    throw INTERP_KERNEL::Exception("Invalid mesh type ! You are expected a structured one whereas in file it is not a structured !");
  _time=mId->checkMeshTimeStep(fid,mName,nstep,dt,it);
  _iteration=dt;
  _order=it;
  //
  _clmesh=MEDCouplingCurveLinearMesh::New();
  INTERP_KERNEL::AutoPtr<int> stGrid=new int[Mdim];
  MEDFILESAFECALLERRD0(MEDmeshGridStructRd,(fid,mName.c_str(),dt,it,stGrid));
  _clmesh->setNodeGridStructure(stGrid,((int *)stGrid)+Mdim);
  med_bool chgt=MED_FALSE,trsf=MED_FALSE;
  int nbNodes(MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_NODE,MED_NONE,MED_COORDINATE,MED_NO_CMODE,&chgt,&trsf));
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
  if((MEDCouplingUMesh *)_m==0)
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

MEDFileUMeshSplitL1::MEDFileUMeshSplitL1(const MEDFileUMeshL2& l2, const std::string& mName, int id):_m(this)
{
  const std::vector< MCAuto<MEDFileUMeshPerType> >& v=l2.getLev(id);
  if(v.empty())
    return;
  int sz=v.size();
  std::vector<const MEDCoupling1GTUMesh *> ms(sz);
  std::vector<const DataArrayInt *> fams(sz),nums(sz);
  std::vector<const DataArrayChar *> names(sz);
  std::vector<const PartDefinition *> pds(sz);
  for(int i=0;i<sz;i++)
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
      for(int i=0;i<sz;i++)
        fams[i]=v[i]->getFam();
      if(sz!=1)
        _fam=DataArrayInt::Aggregate(fams);
      else
        {
          fams[0]->incrRef();
          _fam=const_cast<DataArrayInt *>(fams[0]);
        }
    }
  if(l2.isNumDefinedOnLev(id))
    {
      for(int i=0;i<sz;i++)
        nums[i]=v[i]->getNum();
      if(sz!=1)
        _num=DataArrayInt::Aggregate(nums);
      else
        {
          nums[0]->incrRef();
          _num=const_cast<DataArrayInt *>(nums[0]);
        }
      computeRevNum();
    }
  if(l2.isNamesDefinedOnLev(id))
    {
      for(int i=0;i<sz;i++)
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
  ret.push_back((const DataArrayInt*)_fam);
  ret.push_back((const DataArrayInt*)_num);
  ret.push_back((const DataArrayInt*)_rev_num);
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
  if((const DataArrayInt *)_fam)
    ret->_fam=_fam->deepCopy();
  if((const DataArrayInt *)_num)
    ret->_num=_num->deepCopy();
  if((const DataArrayInt *)_rev_num)
    ret->_rev_num=_rev_num->deepCopy();
  if((const DataArrayAsciiChar *)_names)
    ret->_names=_names->deepCopy();
  return ret.retn();
}

void MEDFileUMeshSplitL1::checkConsistency() const
{
  if (!_fam || _fam->getNumberOfTuples() != getSize())
    throw INTERP_KERNEL::Exception("MEDFileUMeshSplitL1::checkConsistency(): internal family array has an invalid size!");
  int nbCells = getSize();
  if (_num)
    {
      _num->checkNbOfTuplesAndComp(nbCells,1,"MEDFileUMeshSplitL1::checkConsistency(): inconsistent internal node numbering array!");
      int pos;
      int maxValue=_num->getMaxValue(pos);
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
  const DataArrayAsciiChar *e1=_names;
  const DataArrayAsciiChar *e2=other->_names;
  if((e1==0 && e2!=0) || (e1!=0 && e2==0))
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
      MCAuto<DataArrayInt> da=_m_by_types.getUmesh()->getRenumArrForConsecutiveCellTypesSpec(typmai2,typmai2+MED_N_CELL_FIXED_GEO);
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
  _fam=DataArrayInt::New();
  _fam->alloc(_m_by_types.getSize(),1);
  _fam->fillWithValue(0);
}

bool MEDFileUMeshSplitL1::empty() const
{
  return _m_by_types.empty();
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
  return _m_by_types.getMeshDimension();
}

void MEDFileUMeshSplitL1::simpleRepr(std::ostream& oss) const
{
  std::vector<int> code=_m_by_types.getDistributionOfTypes();
  int nbOfTypes=code.size()/3;
  for(int i=0;i<nbOfTypes;i++)
    {
      INTERP_KERNEL::NormalizedCellType typ=(INTERP_KERNEL::NormalizedCellType) code[3*i];
      oss << "    - Number of cells with type " << INTERP_KERNEL::CellModel::GetCellModel(typ).getRepr() << " : " << code[3*i+1] << std::endl;
    }
}

int MEDFileUMeshSplitL1::getSize() const
{
  return _m_by_types.getSize();
}

MEDCouplingUMesh *MEDFileUMeshSplitL1::getFamilyPart(const int *idsBg, const int *idsEnd, bool renum) const
{
  MCAuto<DataArrayInt> eltsToKeep=_fam->findIdsEqualList(idsBg,idsEnd);
  MEDCouplingUMesh *m=(MEDCouplingUMesh *)_m_by_types.getUmesh()->buildPartOfMySelf(eltsToKeep->begin(),eltsToKeep->end(),true);
  if(renum)
    return renumIfNeeded(m,eltsToKeep->begin());
  return m;
}

DataArrayInt *MEDFileUMeshSplitL1::getFamilyPartArr(const int *idsBg, const int *idsEnd, bool renum) const
{
  MCAuto<DataArrayInt> da=_fam->findIdsEqualList(idsBg,idsEnd);
  if(renum)
    return renumIfNeededArr(da);
  return da.retn();
}

std::vector<INTERP_KERNEL::NormalizedCellType> MEDFileUMeshSplitL1::getGeoTypes() const
{
  return _m_by_types.getGeoTypes();
}

int MEDFileUMeshSplitL1::getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType ct) const
{
  return _m_by_types.getNumberOfCellsWithType(ct);
}

MEDCouplingUMesh *MEDFileUMeshSplitL1::getWholeMesh(bool renum) const
{
  MCAuto<MEDCouplingUMesh> tmp;
  if(renum && ((const DataArrayInt *)_num))
    tmp=_m;
  else
    { tmp=_m_by_types.getUmesh(); if(tmp) tmp->incrRef(); }
  return tmp.retn();
}

int MEDFileUMeshSplitL1::getNumberOfCells() const
{
  return _m_by_types.getNumberOfCells();
}

DataArrayInt *MEDFileUMeshSplitL1::extractFamilyFieldOnGeoType(INTERP_KERNEL::NormalizedCellType gt) const
{
  const DataArrayInt *fam(_fam);
  if(!fam)
    return 0;
  int start(0),stop(0);
  _m_by_types.getStartStopOfGeoTypeWithoutComputation(gt,start,stop);
  return fam->selectByTupleIdSafeSlice(start,stop,1);
}

DataArrayInt *MEDFileUMeshSplitL1::extractNumberFieldOnGeoType(INTERP_KERNEL::NormalizedCellType gt) const
{
  const DataArrayInt *num(_num);
  if(!num)
    return 0;
  int start(0),stop(0);
  _m_by_types.getStartStopOfGeoTypeWithoutComputation(gt,start,stop);
  return num->selectByTupleIdSafeSlice(start,stop,1);
}

DataArrayInt *MEDFileUMeshSplitL1::getOrCreateAndGetFamilyField()
{
  if((DataArrayInt *)_fam)
    return _fam;
  int nbOfTuples=_m_by_types.getSize();
  _fam=DataArrayInt::New(); _fam->alloc(nbOfTuples,1); _fam->fillWithZero();
  return _fam;
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
void MEDFileUMeshSplitL1::setGroupsFromScratch(const std::vector<const MEDCouplingUMesh *>& ms, std::map<std::string,int>& familyIds,
                                               std::map<std::string, std::vector<std::string> >& groups)
{
  std::vector< DataArrayInt * > corr;
  _m=MEDCouplingUMesh::FuseUMeshesOnSameCoords(ms,0,corr);
  std::vector< MCAuto<DataArrayInt> > corrMSafe(corr.begin(),corr.end());
  std::vector< std::vector<int> > fidsOfGroups;
  std::vector< const DataArrayInt * > corr2(corr.begin(),corr.end());
  _fam=DataArrayInt::MakePartition(corr2,((MEDCouplingUMesh *)_m)->getNumberOfCells(),fidsOfGroups);
  int nbOfCells=((MEDCouplingUMesh *)_m)->getNumberOfCells();
  std::map<int,std::string> newfams;
  std::map<int,int> famIdTrad;
  TraduceFamilyNumber(fidsOfGroups,familyIds,famIdTrad,newfams);
  int *w=_fam->getPointer();
  for(int i=0;i<nbOfCells;i++,w++)
    *w=famIdTrad[*w];
}

void MEDFileUMeshSplitL1::write(med_idt fid, const std::string& mName, int mdim) const
{
  std::vector<MEDCoupling1GTUMesh *> ms(_m_by_types.getParts());
  int start=0;
  for(std::vector<MEDCoupling1GTUMesh *>::const_iterator it=ms.begin();it!=ms.end();it++)
    {
      int nbCells=(*it)->getNumberOfCells();
      int end=start+nbCells;
      MCAuto<DataArrayInt> fam,num;
      MCAuto<DataArrayAsciiChar> names;
      if((const DataArrayInt *)_fam)
        fam=_fam->subArray(start,end);
      if((const DataArrayInt *)_num)
        num=_num->subArray(start,end);
      if((const DataArrayAsciiChar *)_names)
        names=static_cast<DataArrayAsciiChar *>(_names->subArray(start,end));
      MEDFileUMeshPerType::Write(fid,mName,mdim,(*it),fam,num,names);
      start=end;
    }
}

void MEDFileUMeshSplitL1::renumberNodesInConn(const int *newNodeNumbersO2N)
{
  _m_by_types.renumberNodesInConnWithoutComputation(newNodeNumbersO2N);
}

void MEDFileUMeshSplitL1::serialize(std::vector<int>& tinyInt, std::vector< MCAuto<DataArrayInt> >& bigArraysI) const
{
  bigArraysI.push_back(_fam);
  bigArraysI.push_back(_num);
  _m_by_types.serialize(tinyInt,bigArraysI);
}

void MEDFileUMeshSplitL1::unserialize(const std::string& name, DataArrayDouble *coo, std::vector<int>& tinyInt, std::vector< MCAuto<DataArrayInt> >& bigArraysI)
{
  _fam=bigArraysI.back(); bigArraysI.pop_back();
  _num=bigArraysI.back(); bigArraysI.pop_back();
  _m_by_types.unserialize(name,coo,tinyInt,bigArraysI);
}

void MEDFileUMeshSplitL1::changeFamilyIdArr(int oldId, int newId)
{
  DataArrayInt *arr=_fam;
  if(arr)
    arr->changeValue(oldId,newId);
}

void MEDFileUMeshSplitL1::setFamilyArr(DataArrayInt *famArr)
{
  if(!famArr)
    {
      _fam=0;
      return ;
    }
  int sz(_m_by_types.getSize());
  famArr->checkNbOfTuplesAndComp(sz,1,"MEDFileUMeshSplitL1::setFamilyArr : Problem in size of Family arr ! ");
  famArr->incrRef();
  _fam=famArr;
}

DataArrayInt *MEDFileUMeshSplitL1::getFamilyField()
{
  return _fam;
}

void MEDFileUMeshSplitL1::setRenumArr(DataArrayInt *renumArr)
{
  if(!renumArr)
    {
      _num=0;
      _rev_num=0;
      return ;
    }
  int sz(_m_by_types.getSize());
  renumArr->checkNbOfTuplesAndComp(sz,1,"MEDFileUMeshSplitL1::setRenumArr : Problem in size of numbering arr ! ");
  renumArr->incrRef();
  _num=renumArr;
  computeRevNum();
}

void MEDFileUMeshSplitL1::setNameArr(DataArrayAsciiChar *nameArr)
{
  if(!nameArr)
    {
      _names=0;
      return ;
    }
  int sz(_m_by_types.getSize());
  nameArr->checkNbOfTuplesAndComp(sz,MED_SNAME_SIZE,"MEDFileUMeshSplitL1::setNameArr : Problem in size of name arr ! ");
  nameArr->incrRef();
  _names=nameArr;
}

MEDCouplingUMesh *MEDFileUMeshSplitL1::Renumber2(const DataArrayInt *renum, MEDCouplingUMesh *m, const int *cellIds)
{
  if(renum==0)
    return m;
  if(cellIds==0)
    m->renumberCells(renum->begin(),true);
  else
    {
      MCAuto<DataArrayInt> locnum=renum->selectByTupleId(cellIds,cellIds+m->getNumberOfCells());
      m->renumberCells(locnum->begin(),true);
    }
  return m;
}

MEDFileUMeshSplitL1 *MEDFileUMeshSplitL1::Unserialize(const std::string& name, DataArrayDouble *coo, std::vector<int>& tinyInt, std::vector< MCAuto<DataArrayInt> >& bigArraysI)
{
  MCAuto<MEDFileUMeshSplitL1> ret(new MEDFileUMeshSplitL1);
  ret->unserialize(name,coo,tinyInt,bigArraysI);
  return ret.retn();
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
  return renum->selectByTupleId(da->begin(),da->end());
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
  //tony
}

void MEDFileUMeshSplitL1::computeRevNum() const
{
  int pos;
  int maxValue=_num->getMaxValue(pos);
  _rev_num=_num->invertArrayN2O2O2N(maxValue+1);
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
      for(std::vector< MCAuto<MEDCoupling1GTUMesh> >::iterator it=_m_parts.begin();it!=_m_parts.end();it++)
        {
          MEDCoupling1GTUMesh *tmp(*it);
          if(tmp)
            tmp->setName(name);
        }
    }
}

void MEDFileUMeshAggregateCompute::assignParts(const std::vector< const MEDCoupling1GTUMesh * >& mParts)
{
  std::size_t sz(mParts.size());
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
  _m=0;
}

void MEDFileUMeshAggregateCompute::assignDefParts(const std::vector<const PartDefinition *>& partDefs)
{
  if(_mp_time<_m_time)
    throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::assignDefParts : the parts require a computation !");
  std::size_t sz(partDefs.size());
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

int MEDFileUMeshAggregateCompute::getNumberOfCells() const
{
  if(_mp_time<=_m_time)
    return _m->getNumberOfCells();
  int ret(0);
  for(std::vector< MCAuto<MEDCoupling1GTUMesh> >::const_iterator it=_m_parts.begin();it!=_m_parts.end();it++)
    ret+=(*it)->getNumberOfCells();
  return ret;
}

std::vector<INTERP_KERNEL::NormalizedCellType> MEDFileUMeshAggregateCompute::getGeoTypes() const
{
  if(_mp_time>=_m_time)
    {
      std::size_t sz(_m_parts.size());
      std::vector<INTERP_KERNEL::NormalizedCellType> ret(sz);
      for(std::size_t i=0;i<sz;i++)
        ret[i]=_m_parts[i]->getCellModelEnum();
      return ret;
    }
  else
    return _m->getAllGeoTypesSorted();
}

int MEDFileUMeshAggregateCompute::getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType ct) const
{
  if(_mp_time>=_m_time)
    {
      for(std::vector< MCAuto<MEDCoupling1GTUMesh> >::const_iterator it=_m_parts.begin();it!=_m_parts.end();it++)
        {
          const MEDCoupling1GTUMesh *elt(*it);
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

MEDCoupling1GTUMesh *MEDFileUMeshAggregateCompute::retrievePartWithoutComputation(INTERP_KERNEL::NormalizedCellType gt) const
{
  std::vector<MEDCoupling1GTUMesh *> v(retrievePartsWithoutComputation());
  std::size_t sz(v.size());
  for(std::size_t i=0;i<sz;i++)
    {
      if(v[i])
        if(v[i]->getCellModelEnum()==gt)
          return v[i];
    }
  throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::getPartWithoutComputation : the geometric type is not existing !");
}

void MEDFileUMeshAggregateCompute::getStartStopOfGeoTypeWithoutComputation(INTERP_KERNEL::NormalizedCellType gt, int& start, int& stop) const
{
  start=0; stop=0;
  std::vector<MEDCoupling1GTUMesh *> v(retrievePartsWithoutComputation());
  std::size_t sz(v.size());
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

void MEDFileUMeshAggregateCompute::renumberNodesInConnWithoutComputation(const int *newNodeNumbersO2N)
{
  if(_mp_time>_m_time)
    {
      for(std::vector< MCAuto<MEDCoupling1GTUMesh> >::iterator it=_m_parts.begin();it!=_m_parts.end();it++)
        {
          MEDCoupling1GTUMesh *m(*it);
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
  std::vector< MCAuto<MEDCouplingUMesh> > msMSafe(ms.begin(),ms.end());
  std::size_t sz(msMSafe.size());
  _m_parts.resize(sz);
  for(std::size_t i=0;i<sz;i++)
    _m_parts[i]=MEDCoupling1GTUMesh::New(ms[i]);
  _part_def.clear();
  _part_def.resize(_m_parts.size());
  _mp_time=std::max(_mp_time,_m_time);
}

const PartDefinition *MEDFileUMeshAggregateCompute::getPartDefOfWithoutComputation(INTERP_KERNEL::NormalizedCellType gt) const
{
  if(_mp_time<_m_time)
    throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::getPartDefOfWithoutComputation : the parts require a computation !");
  if(_m_parts.size()!=_part_def.size())
    throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::getPartDefOfWithoutComputation : size of arrays are expected to be the same !");
  std::size_t sz(_m_parts.size());
  for(std::size_t i=0;i<sz;i++)
    {
      const MEDCoupling1GTUMesh *mesh(_m_parts[i]);
      if(mesh)
        if(mesh->getCellModelEnum()==gt)
          return _part_def[i];
    }
  throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::getPartDefOfWithoutComputation : The input geo type is not existing in this !");
}

void MEDFileUMeshAggregateCompute::serialize(std::vector<int>& tinyInt, std::vector< MCAuto<DataArrayInt> >& bigArraysI) const
{
  if(_mp_time<_m_time)
    throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::serialize : the parts require a computation !");
  std::size_t sz(_m_parts.size());
  tinyInt.push_back((int)sz);
  for(std::size_t i=0;i<sz;i++)
    {
      const MEDCoupling1GTUMesh *mesh(_m_parts[i]);
      if(!mesh)
        throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::serialize : one part is empty !");
      tinyInt.push_back(mesh->getCellModelEnum());
      const MEDCoupling1SGTUMesh *mesh1(dynamic_cast<const MEDCoupling1SGTUMesh *>(mesh));
      const MEDCoupling1DGTUMesh *mesh2(dynamic_cast<const MEDCoupling1DGTUMesh *>(mesh));
      if(mesh1)
        {
          DataArrayInt *elt(mesh1->getNodalConnectivity());
          if(elt)
            elt->incrRef();
          MCAuto<DataArrayInt> elt1(elt);
          bigArraysI.push_back(elt1);
        }
      else if(mesh2)
        {
          DataArrayInt *elt1(mesh2->getNodalConnectivity()),*elt2(mesh2->getNodalConnectivityIndex());
          if(elt1)
            elt1->incrRef();
          if(elt2)
            elt2->incrRef();
          MCAuto<DataArrayInt> elt11(elt1),elt22(elt2);
          bigArraysI.push_back(elt11); bigArraysI.push_back(elt22);
        }
      else
        throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::serialize : unrecognized single geo type mesh !");
      const PartDefinition *pd(_part_def[i]);
      if(!pd)
        tinyInt.push_back(-1);
      else
        {
          std::vector<int> tinyTmp;
          pd->serialize(tinyTmp,bigArraysI);
          tinyInt.push_back((int)tinyTmp.size());
          tinyInt.insert(tinyInt.end(),tinyTmp.begin(),tinyTmp.end());
        }
    }
}

void MEDFileUMeshAggregateCompute::unserialize(const std::string& name, DataArrayDouble *coo, std::vector<int>& tinyInt, std::vector< MCAuto<DataArrayInt> >& bigArraysI)
{
  int nbParts(tinyInt.back()); tinyInt.pop_back();
  _part_def.clear(); _part_def.resize(nbParts);
  _m_parts.clear(); _m_parts.resize(nbParts);
  for(int i=0;i<nbParts;i++)
    {
      INTERP_KERNEL::NormalizedCellType tp((INTERP_KERNEL::NormalizedCellType) tinyInt.back()); tinyInt.pop_back();
      MCAuto<MEDCoupling1GTUMesh> mesh(MEDCoupling1GTUMesh::New(name,tp));
      mesh->setCoords(coo);
      MEDCoupling1SGTUMesh *mesh1(dynamic_cast<MEDCoupling1SGTUMesh *>((MEDCoupling1GTUMesh *) mesh));
      MEDCoupling1DGTUMesh *mesh2(dynamic_cast<MEDCoupling1DGTUMesh *>((MEDCoupling1GTUMesh *) mesh));
      if(mesh1)
        {
          mesh1->setNodalConnectivity(bigArraysI.back()); bigArraysI.pop_back();
        }
      else if(mesh2)
        {
          MCAuto<DataArrayInt> elt0,elt1;
          elt0=bigArraysI.back(); bigArraysI.pop_back();
          elt1=bigArraysI.back(); bigArraysI.pop_back();
          mesh2->setNodalConnectivity(elt0,elt1);
        }
      else
        throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::unserialize : unrecognized single geo type mesh !");
      _m_parts[i]=mesh;
      int pdid(tinyInt.back()); tinyInt.pop_back();
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
  for(std::vector< MCAuto<MEDCoupling1GTUMesh> >::const_iterator it=_m_parts.begin();it!=_m_parts.end();it++)
    {
      const MEDCoupling1GTUMesh *elt(*it);
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
  std::size_t ret(_m_parts.size()*sizeof(MCAuto<MEDCoupling1GTUMesh>));
  return ret;
}

std::vector<const BigMemoryObject *> MEDFileUMeshAggregateCompute::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  for(std::vector< MCAuto<MEDCoupling1GTUMesh> >::const_iterator it=_m_parts.begin();it!=_m_parts.end();it++)
    ret.push_back((const MEDCoupling1GTUMesh *)*it);
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
  std::size_t sz(_part_def.size());
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
  for(std::vector< MCAuto<MEDCoupling1GTUMesh> >::iterator it=_m_parts.begin();it!=_m_parts.end();it++)
    {
      const MEDCoupling1GTUMesh *elt(*it);
      if(elt)
        {
          MCAuto<MEDCouplingMesh> elt2(elt->clone(false));
          *it=DynamicCastSafe<MEDCouplingMesh,MEDCoupling1GTUMesh>(elt2);
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
  if((m1==0 && m2!=0) || (m1!=0 && m2==0))
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
  std::size_t sz(_part_def.size());
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
      bool ret(pd0->isEqual(pd1,what));
      if(!ret)
        return false;
    }
  return true;
}

void MEDFileUMeshAggregateCompute::checkConsistency() const
{
  if(_mp_time >= _m_time)
    for(std::vector< MCAuto<MEDCoupling1GTUMesh> >::const_iterator it=_m_parts.begin();
        it!=_m_parts.end(); it++)
      (*it)->checkConsistency();
  else
    _m->checkConsistency();
}

void MEDFileUMeshAggregateCompute::clearNonDiscrAttributes() const
{
  for(std::vector< MCAuto<MEDCoupling1GTUMesh> >::const_iterator it=_m_parts.begin();it!=_m_parts.end();it++)
    MEDFileUMeshSplitL1::ClearNonDiscrAttributes(*it);
  MEDFileUMeshSplitL1::ClearNonDiscrAttributes(_m);
}

void MEDFileUMeshAggregateCompute::synchronizeTinyInfo(const MEDFileMesh& master) const
{
  for(std::vector< MCAuto<MEDCoupling1GTUMesh> >::const_iterator it=_m_parts.begin();it!=_m_parts.end();it++)
    {
      const MEDCoupling1GTUMesh *tmp(*it);
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
    return ((const MEDCouplingUMesh *)_m)==0;
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

std::vector<int> MEDFileUMeshAggregateCompute::getDistributionOfTypes() const
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
      std::vector<int> ret;
      for(std::vector< MCAuto<MEDCoupling1GTUMesh> >::const_iterator it=_m_parts.begin();it!=_m_parts.end();it++)
        {
          const MEDCoupling1GTUMesh *tmp(*it);
          if(!tmp)
            throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::getDistributionOfTypes : part mesh contains null instance !");
          std::vector<int> ret0(tmp->getDistributionOfTypes());
          ret.insert(ret.end(),ret0.begin(),ret0.end());
        }
      return ret;
    }
}

int MEDFileUMeshAggregateCompute::getSize() const
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
      int ret=0;
      for(std::vector< MCAuto<MEDCoupling1GTUMesh> >::const_iterator it=_m_parts.begin();it!=_m_parts.end();it++)
        {
          const MEDCoupling1GTUMesh *m(*it);
          if(!m)
            throw INTERP_KERNEL::Exception("MEDFileUMeshAggregateCompute::getSize : part mesh contains null instance !");
          ret+=m->getNumberOfCells();
        }
      return ret;
    }
}

void MEDFileUMeshAggregateCompute::setCoords(DataArrayDouble *coords)
{
  for(std::vector< MCAuto<MEDCoupling1GTUMesh> >::iterator it=_m_parts.begin();it!=_m_parts.end();it++)
    {
      MEDCoupling1GTUMesh *tmp(*it);
      if(tmp)
        (*it)->setCoords(coords);
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
  for(std::vector< MCAuto<DataArray> >::const_iterator it=_vars.begin();it!=_vars.end();it++)
    ret.push_back(*it);
  return ret;
}

MEDFileEltStruct4Mesh::MEDFileEltStruct4Mesh(med_idt fid, const std::string& mName, int dt, int it, int iterOnStEltOfMesh, MEDFileMeshReadSelector *mrs)
{
  med_geometry_type geoType;
  INTERP_KERNEL::AutoPtr<char> geoTypeName(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  MEDFILESAFECALLERRD0(MEDmeshEntityInfo,(fid,mName.c_str(),dt,it,MED_STRUCT_ELEMENT,iterOnStEltOfMesh+1,geoTypeName,&geoType));
  _geo_type=geoType;
  _geo_type_name=MEDLoaderBase::buildStringFromFortran(geoTypeName,MED_NAME_SIZE);
  int nCells(0);
  {
    med_bool chgt=MED_FALSE,trsf=MED_FALSE;
    nCells=MEDmeshnEntity(fid,mName.c_str(),dt,it,MED_STRUCT_ELEMENT,geoType,MED_CONNECTIVITY,MED_NODAL,&chgt,&trsf);
  }
  MCAuto<MEDFileMeshSupports> mss(MEDFileMeshSupports::New(fid));
  MCAuto<MEDFileStructureElements> mse(MEDFileStructureElements::New(fid,mss));
  int nbEntities(mse->getNumberOfNodesPerSE(_geo_type_name));
  _conn=DataArrayInt::New(); _conn->alloc(nCells*nbEntities);
  MEDFILESAFECALLERRD0(MEDmeshElementConnectivityRd,(fid,mName.c_str(),dt,it,MED_STRUCT_ELEMENT,_geo_type,MED_NODAL,MED_FULL_INTERLACE,_conn->getPointer()));
  _conn->applyLin(1,-1);
  _conn->rearrange(nbEntities);
  _common=MEDFileUMeshPerTypeCommon::New();
  _common->loadCommonPart(fid,mName.c_str(),dt,it,nCells,geoType,MED_STRUCT_ELEMENT,mrs);
  std::vector<std::string> vns(mse->getVarAttsOf(_geo_type_name));
  std::size_t sz(vns.size());
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
