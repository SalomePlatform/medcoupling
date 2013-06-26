// Copyright (C) 2007-2013  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (CEA/DEN)

#include "MEDLoader.hxx"
#include "MEDLoaderBase.hxx"
#include "MEDFileUtilities.hxx"
#include "MEDFileMesh.hxx"
#include "MEDFileField.hxx"
#include "CellModel.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingGaussLocalization.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"

#include "InterpKernelAutoPtr.hxx"

#include "med.h"

#include <string>
#include <limits>
#include <cstring>
#include <sstream>
#include <fstream>
#include <numeric>
#include <iterator>
#include <algorithm>

med_geometry_type typmai[MED_N_CELL_FIXED_GEO] = { MED_POINT1,
                                                   MED_SEG2,
                                                   MED_SEG3,
                                                   MED_SEG4,
                                                   MED_TRIA3,
                                                   MED_QUAD4,
                                                   MED_TRIA6,
                                                   MED_TRIA7,
                                                   MED_QUAD8,
                                                   MED_QUAD9,
                                                   MED_TETRA4,
                                                   MED_PYRA5,
                                                   MED_PENTA6,
                                                   MED_HEXA8,
                                                   MED_OCTA12,
                                                   MED_TETRA10,
                                                   MED_PYRA13,
                                                   MED_PENTA15,
                                                   MED_HEXA20,
                                                   MED_HEXA27,
                                                   MED_POLYGON,
                                                   MED_POLYHEDRON };

med_geometry_type typmainoeud[1] = { MED_NONE };

INTERP_KERNEL::NormalizedCellType typmai2[MED_N_CELL_FIXED_GEO] = { INTERP_KERNEL::NORM_POINT1,
                                                                    INTERP_KERNEL::NORM_SEG2,
                                                                    INTERP_KERNEL::NORM_SEG3,
                                                                    INTERP_KERNEL::NORM_SEG4,
                                                                    INTERP_KERNEL::NORM_TRI3,
                                                                    INTERP_KERNEL::NORM_QUAD4,
                                                                    INTERP_KERNEL::NORM_TRI6,
                                                                    INTERP_KERNEL::NORM_TRI7,
                                                                    INTERP_KERNEL::NORM_QUAD8,
                                                                    INTERP_KERNEL::NORM_QUAD9,
                                                                    INTERP_KERNEL::NORM_TETRA4,
                                                                    INTERP_KERNEL::NORM_PYRA5,
                                                                    INTERP_KERNEL::NORM_PENTA6,
                                                                    INTERP_KERNEL::NORM_HEXA8,
                                                                    INTERP_KERNEL::NORM_HEXGP12,
                                                                    INTERP_KERNEL::NORM_TETRA10,
                                                                    INTERP_KERNEL::NORM_PYRA13,
                                                                    INTERP_KERNEL::NORM_PENTA15,
                                                                    INTERP_KERNEL::NORM_HEXA20,
                                                                    INTERP_KERNEL::NORM_HEXA27,
                                                                    INTERP_KERNEL::NORM_POLYGON,
                                                                    INTERP_KERNEL::NORM_POLYHED };

med_geometry_type typmai3[34] = { MED_POINT1,//0
                                  MED_SEG2,//1
                                  MED_SEG3,//2
                                  MED_TRIA3,//3
                                  MED_QUAD4,//4
                                  MED_POLYGON,//5
                                  MED_TRIA6,//6
                                  MED_TRIA7,//7
                                  MED_QUAD8,//8
                                  MED_QUAD9,//9
                                  MED_SEG4,//10
                                  MED_NONE,//11
                                  MED_NONE,//12
                                  MED_NONE,//13
                                  MED_TETRA4,//14
                                  MED_PYRA5,//15
                                  MED_PENTA6,//16
                                  MED_NONE,//17
                                  MED_HEXA8,//18
                                  MED_NONE,//19
                                  MED_TETRA10,//20
                                  MED_NONE,//21
                                  MED_OCTA12,//22
                                  MED_PYRA13,//23
                                  MED_NONE,//24
                                  MED_PENTA15,//25
                                  MED_NONE,//26
                                  MED_HEXA27,//27
                                  MED_NONE,//28
                                  MED_NONE,//29
                                  MED_HEXA20,//30
                                  MED_POLYHEDRON,//31
                                  MED_NONE,//32
                                  MED_NONE//33
};

double MEDLoader::_EPS_FOR_NODE_COMP=1.e-12;

int MEDLoader::_COMP_FOR_CELL=0;

int MEDLoader::_TOO_LONG_STR=0;

using namespace ParaMEDMEM;

/// @cond INTERNAL

namespace MEDLoaderNS
{
  int readUMeshDimFromFile(const char *fileName, const char *meshName, std::vector<int>& possibilities);
  void dispatchElems(int nbOfElemCell, int nbOfElemFace, int& nbOfElem, med_entity_type& whichEntity);
  void writeFieldWithoutReadingAndMappingOfMeshInFile(const char *fileName, const ParaMEDMEM::MEDCouplingFieldDouble *f, bool writeFromScratch);
  med_int getIdFromMeshName(med_idt fid, const char *meshName, std::string& trueMeshName) throw(INTERP_KERNEL::Exception);
  std::vector<std::string> getMeshNamesFid(med_idt fid);
}

/// @endcond


/// @cond INTERNAL

/*!
 * This method returns a first quick overview of mesh with name 'meshName' into the file 'fileName'.
 * @param possibilities the relativeToMeshDim authorized to returned maxdim. This vector is systematically cleared at the begin of this method.
 * @return the maximal mesh dimension of specified mesh. If nothing found -1 is returned.
 */
int MEDLoaderNS::readUMeshDimFromFile(const char *fileName, const char *meshName, std::vector<int>& possibilities)
{
  possibilities.clear();
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  int ret;
  std::set<int> poss;
  char nommaa[MED_NAME_SIZE+1];
  char maillage_description[MED_COMMENT_SIZE+1];
  med_mesh_type type_maillage;
  med_int Sdim,Mdim;
  std::string trueMeshName;
  med_int meshId=getIdFromMeshName(fid,meshName,trueMeshName);
  INTERP_KERNEL::AutoPtr<char> dt_unit=MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE);
  med_sorting_type sortingType;
  med_int nstep;
  med_axis_type axisType;
  int naxis=MEDmeshnAxis(fid,meshId);
  INTERP_KERNEL::AutoPtr<char> axisname=MEDLoaderBase::buildEmptyString(naxis*MED_SNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> axisunit=MEDLoaderBase::buildEmptyString(naxis*MED_SNAME_SIZE);
  MEDmeshInfo(fid,meshId,nommaa,&Sdim,&Mdim,&type_maillage,maillage_description,dt_unit,&sortingType,&nstep,&axisType,axisname,axisunit);
  // limitation
  if(nstep!=1)
    {
      throw INTERP_KERNEL::Exception("multisteps on mesh not managed yet !");
    } 
  med_int numdt,numit;
  med_float dt;
  MEDmeshComputationStepInfo(fid,nommaa,1,&numdt,&numit,&dt);
  // endlimitation
  for(int i=0;i<MED_N_CELL_GEO_FIXED_CON;i++)
    {
      med_geometry_type curMedType=typmai[i];
      med_bool changement,transformation;
      int curNbOfElemM=MEDmeshnEntity(fid,nommaa,numdt,numit,MED_CELL,curMedType,MED_CONNECTIVITY,MED_NODAL,&changement,&transformation);
      int curNbOfElemF=MEDmeshnEntity(fid,nommaa,numdt,numit,MED_CELL,curMedType,MED_CONNECTIVITY,MED_NODAL,&changement,&transformation);//limitation
      int curNbOfElem;
      med_entity_type whichEntity;
      MEDLoaderNS::dispatchElems(curNbOfElemM,curNbOfElemF,curNbOfElem,whichEntity);
      if(curNbOfElem>0)
        {
          INTERP_KERNEL::NormalizedCellType type=typmai2[i];
          int curDim=(int)INTERP_KERNEL::CellModel::GetCellModel(type).getDimension();
          poss.insert(curDim);
        }
    }
  if(!poss.empty())
    {
      ret=*poss.rbegin();
      for(std::set<int>::const_reverse_iterator it=poss.rbegin();it!=poss.rend();it++)
        possibilities.push_back(*it-ret);
    }
  else
    ret=-2;
  return ret;
}

med_int MEDLoaderNS::getIdFromMeshName(med_idt fid, const char *meshName, std::string& trueMeshName) throw(INTERP_KERNEL::Exception)
{
  if(meshName==0)
    {
      std::vector<std::string> meshes=getMeshNamesFid(fid);
      if(meshes.empty())
        throw INTERP_KERNEL::Exception("No mesh in file");
      trueMeshName=meshes[0];
      return 1;
    }
  std::string meshNameStr(meshName);
  std::vector<std::string> meshes=getMeshNamesFid(fid);
  if(meshes.empty())
    throw INTERP_KERNEL::Exception("No mesh in file");
  std::vector<std::string>::iterator iter=std::find(meshes.begin(),meshes.end(),meshNameStr);
  if(iter==meshes.end())
    {
      std::ostringstream os2;
      os2 << "MeshName '" << meshName << "' not in file : meshes available : ";
      std::copy(meshes.begin(),meshes.end(),std::ostream_iterator<std::string>(os2," "));
      throw INTERP_KERNEL::Exception(os2.str().c_str());
    }
  trueMeshName=meshName;
  return iter-meshes.begin()+1;
}

std::vector<std::string> MEDLoaderNS::getMeshNamesFid(med_idt fid)
{
  med_mesh_type type_maillage;
  char maillage_description[MED_COMMENT_SIZE+1];
  char dtunit[MED_COMMENT_SIZE+1];
  med_int space_dim;
  med_int mesh_dim;
  char nommaa[MED_NAME_SIZE+1];
  med_axis_type axistype;
  med_sorting_type stype;
  med_int n=MEDnMesh(fid);
  std::vector<std::string> ret(n);
  for(int i=0;i<n;i++)
    {
      int naxis=MEDmeshnAxis(fid,i+1);
      INTERP_KERNEL::AutoPtr<char> axisname=MEDLoaderBase::buildEmptyString(naxis*MED_SNAME_SIZE);
      INTERP_KERNEL::AutoPtr<char> axisunit=MEDLoaderBase::buildEmptyString(naxis*MED_SNAME_SIZE);
      int nstep;
      MEDmeshInfo(fid,i+1,nommaa,&space_dim,&mesh_dim,&type_maillage,maillage_description,dtunit,&stype,&nstep,&axistype,axisname,axisunit);
      std::string cur=MEDLoaderBase::buildStringFromFortran(nommaa,sizeof(nommaa));
      ret[i]=cur;
    }
  return ret;
}

/*!
 * This methods allows to merger all entities and to considerate only cell types.
 */
void MEDLoaderNS::dispatchElems(int nbOfElemCell, int nbOfElemFace, int& nbOfElem, med_entity_type& whichEntity)
{
  if(nbOfElemCell>=nbOfElemFace)
    {
      whichEntity=MED_CELL;
      nbOfElem=nbOfElemCell;
    }
  else
    {
      whichEntity=MED_CELL;
      nbOfElem=nbOfElemFace;
    }
}

/// @endcond

/*!
 * This method sets the epsilon value used for node comparison when trying to buid a profile for a field on node/cell on an already written mesh.
 */
void MEDLoader::SetEpsilonForNodeComp(double val) throw(INTERP_KERNEL::Exception)
{
  _EPS_FOR_NODE_COMP=val;
}

/*!
 * This method sets the policy comparison when trying to fit the already written mesh on a field. The semantic of the policy is specified in MEDCouplingUMesh::zipConnectivityTraducer.
 */
void MEDLoader::SetCompPolicyForCell(int val) throw(INTERP_KERNEL::Exception)
{
  _COMP_FOR_CELL=val;
}

/*!
 * This method set the behaviour of MEDLoader when a too long string is seen in datastructure before copy it in MED file.
 * By default (0) an exception is thrown. If equal to 1 a warning is emitted in std_err but no exception is thrown.
 */
void MEDLoader::SetTooLongStrPolicy(int val) throw(INTERP_KERNEL::Exception)
{
  _TOO_LONG_STR=val;
}

/*!
 * Given a 'fileName' and a 'meshName' this method returns global information concerning this mesh.
 * It returns, in this order :
 * - number of cells sorted by dimension and by geometry type. The first entry in the vector is the maximal dimension, the 2nd in the vector is the maximal dimension-1...
 * - the mesh dimension
 * - the space dimension
 * - the number of nodes
 */
std::vector< std::vector< std::pair<INTERP_KERNEL::NormalizedCellType,int> > > MEDLoader::GetUMeshGlobalInfo(const char *fileName, const char *meshName, int &meshDim, int& spaceDim, int& numberOfNodes) throw(INTERP_KERNEL::Exception)
{
  CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  std::set<int> poss;
  char nommaa[MED_NAME_SIZE+1];
  char maillage_description[MED_COMMENT_SIZE+1];
  med_mesh_type type_maillage;
  std::string trueMeshName;
  med_int meshId=MEDLoaderNS::getIdFromMeshName(fid,meshName,trueMeshName);
  INTERP_KERNEL::AutoPtr<char> dt_unit=MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE);
  med_sorting_type sortingType;
  med_int nstep;
  med_axis_type axisType;
  int naxis=MEDmeshnAxis(fid,meshId);
  INTERP_KERNEL::AutoPtr<char> axisname=MEDLoaderBase::buildEmptyString(naxis*MED_SNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> axisunit=MEDLoaderBase::buildEmptyString(naxis*MED_SNAME_SIZE);
  MEDmeshInfo(fid,meshId,nommaa,&spaceDim,&meshDim,&type_maillage,maillage_description,dt_unit,&sortingType,&nstep,&axisType,axisname,axisunit); 
  if(type_maillage!=MED_UNSTRUCTURED_MESH)
    {
      std::ostringstream oss; oss << "MEDLoader::GetUMeshGlobalInfo : Mesh \""<< meshName << "\" in file \"" << fileName;
      oss << "\" exists but it is not an unstructured mesh ! This method is not relevant for mesh types that are not unstructured !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  // limitation
  if(nstep!=1)
    throw INTERP_KERNEL::Exception("MEDLoader::GetUMeshGlobalInfo : multisteps on mesh not managed !");
  med_int numdt,numit;
  med_float dt;
  MEDmeshComputationStepInfo(fid,nommaa,1,&numdt,&numit,&dt);
  // endlimitation
  std::vector<int> dims;
  std::vector< std::pair<INTERP_KERNEL::NormalizedCellType,int> > geoTypes;
  med_bool changement,transformation;
  for(int i=0;i<MED_N_CELL_FIXED_GEO;i++)
    {
      med_geometry_type curMedType=typmai[i];
      int curNbOfElemM=MEDmeshnEntity(fid,nommaa,numdt,numit,MED_CELL,curMedType,MED_CONNECTIVITY,MED_NODAL,&changement,&transformation);
      if(curNbOfElemM>0)
        {
          INTERP_KERNEL::NormalizedCellType typp=typmai2[i];
          int mdimCell=INTERP_KERNEL::CellModel::GetCellModel(typp).getDimension();
          dims.push_back(mdimCell);
          geoTypes.push_back(std::pair<INTERP_KERNEL::NormalizedCellType,int>(typp,curNbOfElemM));
        }
    }
  int maxLev=*std::max_element(dims.begin(),dims.end());
  int lowLev=*std::min_element(dims.begin(),dims.end());
  int nbOfLevels=maxLev-lowLev+1;
  std::vector< std::vector< std::pair<INTERP_KERNEL::NormalizedCellType,int> > > ret(nbOfLevels);
  for(std::size_t i=0;i<dims.size();i++)
    {
      ret[maxLev-dims[i]].push_back(geoTypes[i]);
    }
  numberOfNodes=MEDmeshnEntity(fid,nommaa,numdt,numit,MED_NODE,MED_NONE,MED_COORDINATE,MED_NO_CMODE,&changement,&transformation);
  return ret;
}

void MEDLoader::CheckFileForRead(const char *fileName) throw(INTERP_KERNEL::Exception)
{
  MEDFileUtilities::CheckFileForRead(fileName);
}

std::vector<std::string> MEDLoader::GetMeshNames(const char *fileName) throw(INTERP_KERNEL::Exception)
{
  CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  std::vector<std::string> ret=MEDLoaderNS::getMeshNamesFid(fid);
  return ret;
}

std::vector< std::pair<std::string,std::string> > MEDLoader::GetComponentsNamesOfField(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
{
  CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  med_int nbFields=MEDnField(fid);
  std::vector<std::string> fields(nbFields);
  med_field_type typcha;
  for(int i=0;i<nbFields;i++)
    {
      med_int ncomp=MEDfieldnComponent(fid,i+1);
      INTERP_KERNEL::AutoPtr<char> comp=new char[ncomp*MED_SNAME_SIZE+1];
      INTERP_KERNEL::AutoPtr<char> unit=new char[ncomp*MED_SNAME_SIZE+1];
      INTERP_KERNEL::AutoPtr<char> dt_unit=MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE);
      med_int nbPdt;
      med_bool localmesh;
      INTERP_KERNEL::AutoPtr<char> maa_ass=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
      INTERP_KERNEL::AutoPtr<char> nomcha=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
      MEDfieldInfo(fid,i+1,nomcha,maa_ass,&localmesh,&typcha,comp,unit,dt_unit,&nbPdt);
      std::string meshName=MEDLoaderBase::buildStringFromFortran(maa_ass,MED_NAME_SIZE);
      std::string curFieldName=MEDLoaderBase::buildStringFromFortran(nomcha,MED_NAME_SIZE+1);
      if(curFieldName==fieldName)
        {
          std::vector< std::pair<std::string,std::string> > ret(ncomp);
          for(int j=0;j<ncomp;j++)
            ret[j]=std::pair<std::string,std::string>(MEDLoaderBase::buildStringFromFortran(((char *)comp)+j*MED_SNAME_SIZE,MED_SNAME_SIZE),
                                                      MEDLoaderBase::buildStringFromFortran(((char *)unit)+j*MED_SNAME_SIZE,MED_SNAME_SIZE));
          return ret;
        }
      fields[i]=curFieldName;
    }
  std::ostringstream oss; oss << "MEDLoader::GetComponentsNamesOfField : no such field \"" << fieldName << "\" in file \"" << fileName << "\" !" << std::endl;
  oss << "Possible field names are : " << std::endl;
  std::copy(fields.begin(),fields.end(),std::ostream_iterator<std::string>(oss," "));
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

std::vector<std::string> MEDLoader::GetMeshNamesOnField(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
{
  CheckFileForRead(fileName);
  std::vector<std::string> ret;
  //
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  med_int nbFields=MEDnField(fid);
  //
  med_field_type typcha;
  INTERP_KERNEL::AutoPtr<char> dt_unit=MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> nomcha=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  med_bool localmesh;
  //
  for(int i=0;i<nbFields;i++)
    {
      med_int ncomp=MEDfieldnComponent(fid,i+1);
      INTERP_KERNEL::AutoPtr<char> comp=new char[ncomp*MED_SNAME_SIZE+1];
      INTERP_KERNEL::AutoPtr<char> unit=new char[ncomp*MED_SNAME_SIZE+1];
      med_int nbPdt;
      INTERP_KERNEL::AutoPtr<char> maa_ass=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
      MEDfieldInfo(fid,i+1,nomcha,maa_ass,&localmesh,&typcha,comp,unit,dt_unit,&nbPdt);
      std::string meshName=MEDLoaderBase::buildStringFromFortran(maa_ass,MED_NAME_SIZE);
      std::string curFieldName=MEDLoaderBase::buildStringFromFortran(nomcha,MED_NAME_SIZE+1);
      if(curFieldName==fieldName)
        ret.push_back(meshName);
    }
  return ret;
}

std::vector<std::string> MEDLoader::GetMeshFamiliesNames(const char *fileName, const char *meshName) throw(INTERP_KERNEL::Exception)
{
  CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  med_int nfam=MEDnFamily(fid,meshName);
  std::vector<std::string> ret(nfam);
  char nomfam[MED_NAME_SIZE+1];
  med_int numfam;
  for(int i=0;i<nfam;i++)
    {
      int ngro=MEDnFamilyGroup(fid,meshName,i+1);
      med_int natt=MEDnFamily23Attribute(fid,meshName,i+1);
      INTERP_KERNEL::AutoPtr<med_int> attide=new med_int[natt];
      INTERP_KERNEL::AutoPtr<med_int> attval=new med_int[natt];
      INTERP_KERNEL::AutoPtr<char> attdes=new char[MED_COMMENT_SIZE*natt+1];
      INTERP_KERNEL::AutoPtr<char> gro=new char[MED_LNAME_SIZE*ngro+1];
      MEDfamily23Info(fid,meshName,i+1,nomfam,attide,attval,attdes,&numfam,gro);
      std::string cur=MEDLoaderBase::buildStringFromFortran(nomfam,sizeof(nomfam));
      ret[i]=cur;
    }
  return ret;
}


std::vector<std::string> MEDLoader::GetMeshFamiliesNamesOnGroup(const char *fileName, const char *meshName, const char *grpName) throw(INTERP_KERNEL::Exception)
{
  CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  med_int nfam=MEDnFamily(fid,meshName);
  std::vector<std::string> ret;
  char nomfam[MED_NAME_SIZE+1];
  med_int numfam;
  for(int i=0;i<nfam;i++)
    {
      int ngro=MEDnFamilyGroup(fid,meshName,i+1);
      med_int natt=MEDnFamily23Attribute(fid,meshName,i+1);
      INTERP_KERNEL::AutoPtr<med_int> attide=new med_int[natt];
      INTERP_KERNEL::AutoPtr<med_int> attval=new med_int[natt];
      INTERP_KERNEL::AutoPtr<char> attdes=new char[MED_COMMENT_SIZE*natt+1];
      INTERP_KERNEL::AutoPtr<char> gro=new char[MED_LNAME_SIZE*ngro+1];
      MEDfamily23Info(fid,meshName,i+1,nomfam,attide,attval,attdes,&numfam,gro);
      std::string cur=MEDLoaderBase::buildStringFromFortran(nomfam,sizeof(nomfam));
      for(int j=0;j<ngro;j++)
        {
          std::string cur2=MEDLoaderBase::buildStringFromFortran(gro+j*MED_LNAME_SIZE,MED_LNAME_SIZE);
          if(cur2==grpName)
            ret.push_back(cur);
        }
    }
  return ret;
}

std::vector<std::string> MEDLoader::GetMeshGroupsNamesOnFamily(const char *fileName, const char *meshName, const char *famName) throw(INTERP_KERNEL::Exception)
{
  CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  med_int nfam=MEDnFamily(fid,meshName);
  std::vector<std::string> ret;
  char nomfam[MED_NAME_SIZE+1];
  med_int numfam;
  bool found=false;
  for(int i=0;i<nfam && !found;i++)
    {
      int ngro=MEDnFamilyGroup(fid,meshName,i+1);
      med_int natt=MEDnFamily23Attribute(fid,meshName,i+1);
      INTERP_KERNEL::AutoPtr<med_int> attide=new med_int[natt];
      INTERP_KERNEL::AutoPtr<med_int> attval=new med_int[natt];
      INTERP_KERNEL::AutoPtr<char> attdes=new char[MED_COMMENT_SIZE*natt+1];
      INTERP_KERNEL::AutoPtr<char> gro=new char[MED_LNAME_SIZE*ngro+1];
      MEDfamily23Info(fid,meshName,i+1,nomfam,attide,attval,attdes,&numfam,gro);
      std::string cur=MEDLoaderBase::buildStringFromFortran(nomfam,sizeof(nomfam));
      found=(cur==famName);
      if(found)
        for(int j=0;j<ngro;j++)
          {
            std::string cur2=MEDLoaderBase::buildStringFromFortran(gro+j*MED_LNAME_SIZE,MED_LNAME_SIZE);
            ret.push_back(cur2);
          }
    }
  if(!found)
    {
      std::ostringstream oss;
      oss << "MEDLoader::GetMeshGroupsNamesOnFamily : no such family \"" << famName << "\" in file \"" << fileName << "\" in mesh \"" << meshName << "\" !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return ret;
}

  
std::vector<std::string> MEDLoader::GetMeshGroupsNames(const char *fileName, const char *meshName) throw(INTERP_KERNEL::Exception)
{
  CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  med_int nfam=MEDnFamily(fid,meshName);
  std::vector<std::string> ret;
  char nomfam[MED_NAME_SIZE+1];
  med_int numfam;
  for(int i=0;i<nfam;i++)
    {
      int ngro=MEDnFamilyGroup(fid,meshName,i+1);
      med_int natt=MEDnFamily23Attribute(fid,meshName,i+1);
      INTERP_KERNEL::AutoPtr<med_int> attide=new med_int[natt];
      INTERP_KERNEL::AutoPtr<med_int> attval=new med_int[natt];
      INTERP_KERNEL::AutoPtr<char> attdes=new char[MED_COMMENT_SIZE*natt+1];
      INTERP_KERNEL::AutoPtr<char> gro=new char[MED_LNAME_SIZE*ngro+1];
      MEDfamily23Info(fid,meshName,i+1,nomfam,attide,attval,attdes,&numfam,gro);
      for(int j=0;j<ngro;j++)
        {
          std::string cur=MEDLoaderBase::buildStringFromFortran(gro+j*MED_LNAME_SIZE,MED_LNAME_SIZE);
          if(std::find(ret.begin(),ret.end(),cur)==ret.end())
            ret.push_back(cur);
        }
    }
  return ret;
}
std::vector<ParaMEDMEM::TypeOfField> MEDLoader::GetTypesOfField(const char *fileName, const char *meshName, const char *fieldName) throw(INTERP_KERNEL::Exception)
{
  CheckFileForRead(fileName);
  std::vector<ParaMEDMEM::TypeOfField> ret;
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  med_int nbFields=MEDnField(fid);
  //
  med_field_type typcha;
  //med_int nbpdtnor=0,pflsize,*pflval,lnsize;
  med_int numdt=0,numo=0;
  med_float dt=0.0;
  char pflname[MED_NAME_SIZE+1]="";
  char locname[MED_NAME_SIZE+1]="";
  char *maa_ass=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  char *nomcha=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  med_bool localmesh;
  //
  for(int i=0;i<nbFields;i++)
    {
      med_int ncomp=MEDfieldnComponent(fid,i+1);
      INTERP_KERNEL::AutoPtr<char> comp=new char[ncomp*MED_SNAME_SIZE+1];
      INTERP_KERNEL::AutoPtr<char> unit=new char[ncomp*MED_SNAME_SIZE+1];
      INTERP_KERNEL::AutoPtr<char> dt_unit=new char[MED_LNAME_SIZE+1];
      med_int nbPdt;
      MEDfieldInfo(fid,i+1,nomcha,maa_ass,&localmesh,&typcha,comp,unit,dt_unit,&nbPdt);
      std::string curFieldName=MEDLoaderBase::buildStringFromFortran(nomcha,MED_NAME_SIZE+1);
      std::string curMeshName=MEDLoaderBase::buildStringFromFortran(maa_ass,MED_NAME_SIZE+1);
      if(curMeshName==meshName)
        {
          if(curFieldName==fieldName)
            {
              int profilesize,nbi;
              if(nbPdt>0)
                {
                  bool found=false;
                  for(int ii=0;ii<nbPdt && !found;ii++)
                    {
                      MEDfieldComputingStepInfo(fid,nomcha,1,&numdt,&numo,&dt);
                      med_int nbOfVal=MEDfieldnValueWithProfile(fid,nomcha,numdt,numo,MED_NODE,MED_NONE,1,MED_COMPACT_PFLMODE,
                                                                pflname,&profilesize,locname,&nbi);
                      if(nbOfVal>0)
                        {
                          ret.push_back(ON_NODES);
                          found=true;
                        }
                    }
                }
              bool found=false;
              for(int j=0;j<MED_N_CELL_FIXED_GEO && !found;j++)
                {
                  if(nbPdt>0)
                    {
                      MEDfieldComputingStepInfo(fid,nomcha,1,&numdt,&numo,&dt);
                      med_int nbOfVal=MEDfieldnValueWithProfile(fid,nomcha,numdt,numo,MED_CELL,typmai[j],1,MED_COMPACT_PFLMODE,
                                                                pflname,&profilesize,locname,&nbi);
                      if(nbOfVal>0)
                        {
                          found=true;
                          ret.push_back(ON_CELLS);
                        }
                    }
                }
            }
        }
    }
  delete [] maa_ass;
  delete [] nomcha;
  return ret;
}

std::vector<std::string> MEDLoader::GetAllFieldNames(const char *fileName) throw(INTERP_KERNEL::Exception)
{
  CheckFileForRead(fileName);
  std::vector<std::string> ret;
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  med_int nbFields=MEDnField(fid);
  med_field_type typcha;
  for(int i=0;i<nbFields;i++)
    {
      med_int ncomp=MEDfieldnComponent(fid,i+1);
      INTERP_KERNEL::AutoPtr<char> comp=new char[ncomp*MED_SNAME_SIZE+1];
      INTERP_KERNEL::AutoPtr<char> unit=new char[ncomp*MED_SNAME_SIZE+1];
      INTERP_KERNEL::AutoPtr<char> nomcha=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
      INTERP_KERNEL::AutoPtr<char> maa_ass=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
      INTERP_KERNEL::AutoPtr<char> dt_unit=new char[MED_LNAME_SIZE+1];
      med_int nbPdt;
      med_bool localmesh;
      MEDfieldInfo(fid,i+1,nomcha,maa_ass,&localmesh,&typcha,comp,unit,dt_unit,&nbPdt);
      ret.push_back(std::string(nomcha));
    }
  return ret;
}

std::vector<std::string> MEDLoader::GetAllFieldNamesOnMesh(const char *fileName, const char *meshName) throw(INTERP_KERNEL::Exception)
{
  CheckFileForRead(fileName);
  std::vector<std::string> ret;
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  med_int nbFields=MEDnField(fid);
  //
  med_field_type typcha;
  char *maa_ass=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  char *nomcha=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  //
  for(int i=0;i<nbFields;i++)
    {
      med_int ncomp=MEDfieldnComponent(fid,i+1);
      INTERP_KERNEL::AutoPtr<char> comp=new char[ncomp*MED_SNAME_SIZE+1];
      INTERP_KERNEL::AutoPtr<char> unit=new char[ncomp*MED_SNAME_SIZE+1];
      INTERP_KERNEL::AutoPtr<char> dt_unit=new char[MED_LNAME_SIZE+1];
      med_int nbPdt;
      med_bool localmesh;
      MEDfieldInfo(fid,i+1,nomcha,maa_ass,&localmesh,&typcha,comp,unit,dt_unit,&nbPdt);
      std::string curFieldName=MEDLoaderBase::buildStringFromFortran(nomcha,MED_NAME_SIZE+1);
      std::string curMeshName=MEDLoaderBase::buildStringFromFortran(maa_ass,MED_NAME_SIZE+1);
      //
      if(curMeshName==meshName)
        ret.push_back(curFieldName);
    }
  delete [] maa_ass;
  delete [] nomcha;
  return ret;
}

std::vector<std::string> MEDLoader::GetFieldNamesOnMesh(ParaMEDMEM::TypeOfField type, const char *fileName, const char *meshName) throw(INTERP_KERNEL::Exception)
{
  CheckFileForRead(fileName);
  switch(type)
    {
    case ON_CELLS:
      return GetCellFieldNamesOnMesh(fileName,meshName);
    case ON_NODES:
      return GetNodeFieldNamesOnMesh(fileName,meshName);
    default:
      throw INTERP_KERNEL::Exception("Type of field specified not managed ! manages are ON_NODES or ON_CELLS !");
    } 
}

std::vector<std::string> MEDLoader::GetCellFieldNamesOnMesh(const char *fileName, const char *meshName) throw(INTERP_KERNEL::Exception)
{
  CheckFileForRead(fileName);
  std::vector<std::string> ret;
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  med_int nbFields=MEDnField(fid);
  //
  med_field_type typcha;
  //med_int nbpdtnor=0,pflsize,*pflval,lnsize;
  med_int numdt=0,numo=0;
  med_float dt=0.0;
  char pflname[MED_NAME_SIZE+1]="";
  char locname[MED_NAME_SIZE+1]="";
  INTERP_KERNEL::AutoPtr<char> maa_ass=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> dt_unit=MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> nomcha=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  med_bool localmesh;
  med_int nbPdt;
  //
  for(int i=0;i<nbFields;i++)
    {
      med_int ncomp=MEDfieldnComponent(fid,i+1);
      INTERP_KERNEL::AutoPtr<char> comp=new char[ncomp*MED_SNAME_SIZE+1];
      INTERP_KERNEL::AutoPtr<char> unit=new char[ncomp*MED_SNAME_SIZE+1];
      MEDfieldInfo(fid,i+1,nomcha,maa_ass,&localmesh,&typcha,comp,unit,dt_unit,&nbPdt);
      std::string curFieldName=MEDLoaderBase::buildStringFromFortran(nomcha,MED_NAME_SIZE+1);
      std::string curMeshName=MEDLoaderBase::buildStringFromFortran(maa_ass,MED_NAME_SIZE+1);
      int profilesize,nbi;
      if(curMeshName==meshName)
        {
          bool found=false;
          for(int j=0;j<MED_N_CELL_FIXED_GEO && !found;j++)
            {
              if(nbPdt>0)
                {
                  MEDfieldComputingStepInfo(fid,nomcha,1,&numdt,&numo,&dt);
                  med_int nbOfVal=MEDfieldnValueWithProfile(fid,nomcha,numdt,numo,MED_CELL,typmai[j],1,MED_COMPACT_PFLMODE,
                                                            pflname,&profilesize,locname,&nbi);
                  if(nbOfVal>0)
                    {
                      found=true;
                      ret.push_back(curFieldName);
                    }
                }
            }
        }
    }
  return ret;
}

std::vector<std::string> MEDLoader::GetNodeFieldNamesOnMesh(const char *fileName, const char *meshName) throw(INTERP_KERNEL::Exception)
{
  CheckFileForRead(fileName);
  std::vector<std::string> ret;
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  med_int nbFields=MEDnField(fid);
  char pflname[MED_NAME_SIZE+1]="";
  char locname[MED_NAME_SIZE+1]="";
  //
  med_field_type typcha;
  med_int numdt=0,numo=0;
  med_float dt=0.0;
  INTERP_KERNEL::AutoPtr<char> maa_ass=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> dt_unit=MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> nomcha=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  med_bool localmesh;
  //
  for(int i=0;i<nbFields;i++)
    {
      med_int ncomp=MEDfieldnComponent(fid,i+1);
      INTERP_KERNEL::AutoPtr<char> comp=new char[ncomp*MED_SNAME_SIZE+1];
      INTERP_KERNEL::AutoPtr<char> unit=new char[ncomp*MED_SNAME_SIZE+1];
      med_int nbPdt;
      MEDfieldInfo(fid,i+1,nomcha,maa_ass,&localmesh,&typcha,comp,unit,dt_unit,&nbPdt);
      std::string curFieldName=MEDLoaderBase::buildStringFromFortran(nomcha,MED_NAME_SIZE+1);
      std::string curMeshName=MEDLoaderBase::buildStringFromFortran(maa_ass,MED_NAME_SIZE+1);
      if(nbPdt>0)
        {
          int profilesize,nbi;
          MEDfieldComputingStepInfo(fid,nomcha,1,&numdt,&numo,&dt);
          med_int nbOfVal=MEDfieldnValueWithProfile(fid,nomcha,numdt,numo,MED_NODE,MED_NONE,1,MED_COMPACT_PFLMODE,
                                                    pflname,&profilesize,locname,&nbi);
          if(curMeshName==meshName && nbOfVal>0)
            {
              ret.push_back(curFieldName);
            }
        }
    }
  return ret;
}

std::vector< std::pair< std::pair<int,int>, double> > MEDLoader::GetAllFieldIterations(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
{
  CheckFileForRead(fileName);
  std::vector< std::pair< std::pair<int,int>, double > > ret;
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  med_int nbFields=MEDnField(fid);
  //
  med_field_type typcha;
  med_int numdt=0,numo=0;
  med_float dt=0.0;
  INTERP_KERNEL::AutoPtr<char> maa_ass=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> dt_unit=MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> nomcha=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  med_bool localmesh;
  //
  std::ostringstream oss; oss << "MEDLoader::GetAllFieldIterations : No field with name \"" << fieldName<< "\" in file \"" << fileName << "\" ! Possible fields are : ";
  for(int i=0;i<nbFields;i++)
    {
      med_int ncomp=MEDfieldnComponent(fid,i+1);
      INTERP_KERNEL::AutoPtr<char> comp=new char[ncomp*MED_SNAME_SIZE+1];
      INTERP_KERNEL::AutoPtr<char> unit=new char[ncomp*MED_SNAME_SIZE+1];
      med_int nbPdt;
      MEDfieldInfo(fid,i+1,nomcha,maa_ass,&localmesh,&typcha,comp,unit,dt_unit,&nbPdt);
      std::string curFieldName=MEDLoaderBase::buildStringFromFortran(nomcha,MED_NAME_SIZE+1);
      if(curFieldName==fieldName)
        {
          for(int k=0;k<nbPdt;k++)
            {
              MEDfieldComputingStepInfo(fid,nomcha,k+1,&numdt,&numo,&dt);
              ret.push_back(std::make_pair(std::make_pair(numdt,numo),dt));
            }
          return ret;
        }
      else
        {
          oss << "\"" << curFieldName << "\"";
          if(i!=nbFields-1) oss << ", ";
        }
    }
  oss << " !";
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

double MEDLoader::GetTimeAttachedOnFieldIteration(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
{
  CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  med_int nbFields=MEDnField(fid);
  //
  med_field_type typcha;
  med_int numdt=0,numo=0;
  med_float dt=0.0;
  med_bool local;
  INTERP_KERNEL::AutoPtr<char> maa_ass=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> dt_unit=MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> nomcha=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  //
  bool found=false;
  bool found2=false;
  double ret=std::numeric_limits<double>::max();
  for(int i=0;i<nbFields && !found;i++)
    {
      med_int ncomp=MEDfieldnComponent(fid,i+1);
      INTERP_KERNEL::AutoPtr<char> comp=new char[ncomp*MED_SNAME_SIZE+1];
      INTERP_KERNEL::AutoPtr<char> unit=new char[ncomp*MED_SNAME_SIZE+1];
      med_int nbPdt;
      MEDfieldInfo(fid,i+1,nomcha,maa_ass,&local,&typcha,comp,unit,dt_unit,&nbPdt);
      std::string curFieldName=MEDLoaderBase::buildStringFromFortran(nomcha,MED_NAME_SIZE+1);
      if(curFieldName==fieldName)
        {
          found=true;
          for(int k=0;k<nbPdt;k++)
            {
              MEDfieldComputingStepInfo(fid,nomcha,k+1,&numdt,&numo,&dt);
              if(numdt==iteration && numo==order)
                {
                  found2=true;
                  ret=dt;
                }
            }
        }
    }
  if(!found || !found2)
    {
      std::ostringstream oss;
      oss << "No such field with name \"" << fieldName << "\" and iteration,order=(" << iteration << "," << order << ") exists in file \"" << fileName << "\" !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return ret;
}

std::vector< std::pair<int,int> > MEDLoader::GetFieldIterations(ParaMEDMEM::TypeOfField type, const char *fileName, const char *meshName, const char *fieldName) throw(INTERP_KERNEL::Exception)
{
  CheckFileForRead(fileName);
  switch(type)
    {
    case ON_CELLS:
      return GetCellFieldIterations(fileName,meshName,fieldName);
    case ON_NODES:
      return GetNodeFieldIterations(fileName,meshName,fieldName);
    default:
      throw INTERP_KERNEL::Exception("Type of field specified not managed ! manages are ON_NODES or ON_CELLS !");
    }
}

std::vector< std::pair<int,int> > MEDLoader::GetCellFieldIterations(const char *fileName, const char *meshName, const char *fieldName) throw(INTERP_KERNEL::Exception)
{
  CheckFileForRead(fileName);
  std::string meshNameCpp(meshName);
  std::vector< std::pair<int,int> > ret;
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  med_int nbFields=MEDnField(fid);
  //
  med_field_type typcha;
  med_int numdt=0,numo=0;
  med_float dt=0.0;
  char pflname[MED_NAME_SIZE+1]="";
  char locname[MED_NAME_SIZE+1]="";
  INTERP_KERNEL::AutoPtr<char> maa_ass=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> dt_unit=MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> nomcha=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  med_bool localmesh;
  //
  std::ostringstream oss; oss << "MEDLoader::GetCellFieldIterations : No cell Field field with name \"" << fieldName<< "\" in file \"" << fileName << "\" ! Possible fields are : ";
  std::set<std::string> s2;
  for(int i=0;i<nbFields;i++)
    {
      med_int ncomp=MEDfieldnComponent(fid,i+1);
      INTERP_KERNEL::AutoPtr<char> comp=new char[ncomp*MED_SNAME_SIZE+1];
      INTERP_KERNEL::AutoPtr<char> unit=new char[ncomp*MED_SNAME_SIZE+1];
      med_int nbPdt;
      MEDfieldInfo(fid,i+1,nomcha,maa_ass,&localmesh,&typcha,comp,unit,dt_unit,&nbPdt);
      std::string curFieldName=MEDLoaderBase::buildStringFromFortran(nomcha,MED_NAME_SIZE+1);
      if(curFieldName==fieldName)
        {
          bool found=false;
          for(int j=0;j<MED_N_CELL_FIXED_GEO && !found;j++)
            {
              for(int k=0;k<nbPdt;k++)
                {
                  int profilesize,nbi;
                  MEDfieldComputingStepInfo(fid,nomcha,k+1,&numdt,&numo,&dt);
                  med_int nbOfVal=MEDfieldnValueWithProfile(fid,nomcha,numdt,numo,MED_CELL,typmai[j],1,MED_COMPACT_PFLMODE,
                                                            pflname,&profilesize,locname,&nbi);
                  std::string maa_ass_cpp(maa_ass);
                  if(nbOfVal>0)
                    {
                      if(meshNameCpp==maa_ass_cpp)
                        {
                          found=true;
                          ret.push_back(std::make_pair(numdt,numo));
                        }
                      else
                        s2.insert(maa_ass_cpp);
                    }
                }
            }
        }
      else
        {
          oss << "\"" << curFieldName << "\"";
          if(i!=nbFields-1) oss << ", ";
        }
    }
  if(ret.empty())
    {
      if(!s2.empty())
        {
          oss << ". Cell Field \"" << fieldName << "\" exists but lies on meshes with names : \"";
          std::copy(s2.begin(),s2.end(),std::ostream_iterator<std::string>(oss,"\", \""));
        }
      oss << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return ret;
}

std::vector< std::pair<int,int> > MEDLoader::GetNodeFieldIterations(const char *fileName, const char *meshName, const char *fieldName) throw(INTERP_KERNEL::Exception)
{
  CheckFileForRead(fileName);
  std::string meshNameCpp(meshName);
  std::vector< std::pair<int,int> > ret;
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  med_int nbFields=MEDnField(fid);
  //
  med_field_type typcha;
  med_int numdt=0,numo=0;
  med_float dt=0.0;
  char pflname[MED_NAME_SIZE+1]="";
  char locname[MED_NAME_SIZE+1]="";
  INTERP_KERNEL::AutoPtr<char> maa_ass=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> dt_unit=MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> nomcha=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  med_bool localmesh;
  //
  std::ostringstream oss; oss << "MEDLoader::GetNodeFieldIterations : No node Field field with name \"" << fieldName<< "\" in file \"" << fileName << "\" ! Possible fields are : ";
  std::set<std::string> s2;
  for(int i=0;i<nbFields;i++)
    {
      med_int ncomp=MEDfieldnComponent(fid,i+1);
      INTERP_KERNEL::AutoPtr<char> comp=new char[ncomp*MED_SNAME_SIZE+1];
      INTERP_KERNEL::AutoPtr<char> unit=new char[ncomp*MED_SNAME_SIZE+1];
      med_int nbPdt;
      MEDfieldInfo(fid,i+1,nomcha,maa_ass,&localmesh,&typcha,comp,unit,dt_unit,&nbPdt);
      std::string curFieldName=MEDLoaderBase::buildStringFromFortran(nomcha,MED_NAME_SIZE+1);
      if(curFieldName==fieldName)
        {
          for(int k=0;k<nbPdt;k++)
            {
              int profilesize,nbi;
              MEDfieldComputingStepInfo(fid,nomcha,k+1,&numdt,&numo,&dt);
              med_int nbOfVal=MEDfieldnValueWithProfile(fid,nomcha,numdt,numo,MED_NODE,MED_NONE,1,MED_COMPACT_PFLMODE,
                                                        pflname,&profilesize,locname,&nbi);
               std::string maa_ass_cpp(maa_ass);
               if(nbOfVal>0)
                 {
                   if(meshNameCpp==maa_ass_cpp)
                     { ret.push_back(std::make_pair(numdt,numo)); }
                   else
                     s2.insert(maa_ass_cpp);
                 }
            }
        }
      else
        {
          oss << "\"" << curFieldName << "\"";
          if(i!=nbFields-1) oss << ", ";
        }
    }
  if(ret.empty())
    {
      if(!s2.empty())
        {
          oss << ". Node Field \"" << fieldName << "\" exists but lies on meshes with names : \"";
          std::copy(s2.begin(),s2.end(),std::ostream_iterator<std::string>(oss,"\", \""));
        }
      oss << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return ret;
}

MEDCouplingUMesh *MEDLoader::ReadUMeshFromFile(const char *fileName, const char *meshName, int meshDimRelToMax) throw(INTERP_KERNEL::Exception)
{
  CheckFileForRead(fileName);
  MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> mm=MEDFileMesh::New(fileName,meshName);
  MEDFileMesh *mmPtr(mm);
  MEDFileUMesh *mmuPtr=dynamic_cast<MEDFileUMesh *>(mmPtr);
  if(!mmuPtr)
    {
      std::ostringstream oss; oss << "MEDLoader::ReadUMeshFromFile : With fileName=\""<< fileName << "\", meshName=\""<< meshName << "\" exists but it is not an unstructured mesh !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
 return  mmuPtr->getMeshAtLevel(meshDimRelToMax,true);
}

ParaMEDMEM::MEDCouplingUMesh *MEDLoader::ReadUMeshFromFile(const char *fileName, int meshDimRelToMax) throw(INTERP_KERNEL::Exception)
{
  CheckFileForRead(fileName);
  MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> mm=MEDFileMesh::New(fileName);
  MEDFileMesh *mmPtr(mm);
  MEDFileUMesh *mmuPtr=dynamic_cast<MEDFileUMesh *>(mmPtr);
  if(!mmuPtr)
    {
      std::ostringstream oss; oss << "MEDLoader::ReadUMeshFromFile : With fileName=\""<< fileName << "\", meshName (the first) =\""<< mm->getName() << "\" exists but it is not an unstructured mesh !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
 return  mmuPtr->getMeshAtLevel(meshDimRelToMax,true);
}

int MEDLoader::ReadUMeshDimFromFile(const char *fileName, const char *meshName) throw(INTERP_KERNEL::Exception)
{
  CheckFileForRead(fileName);
  std::vector<int> poss;
  return MEDLoaderNS::readUMeshDimFromFile(fileName,meshName,poss);
}

ParaMEDMEM::MEDCouplingUMesh *MEDLoader::ReadUMeshFromFamilies(const char *fileName, const char *meshName, int meshDimRelToMax, const std::vector<std::string>& fams) throw(INTERP_KERNEL::Exception)
{
  CheckFileForRead(fileName);
  MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> mm=MEDFileMesh::New(fileName,meshName);
  MEDFileMesh *mmPtr(mm);
  MEDFileUMesh *mmuPtr=dynamic_cast<MEDFileUMesh *>(mmPtr);
  if(!mmuPtr)
    {
      std::ostringstream oss; oss << "MEDLoader::ReadUMeshFromFamilies : With fileName=\""<< fileName << "\", meshName (the first) =\""<< mm->getName() << "\" exists but it is not an unstructured mesh !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
    return mmuPtr->getFamilies(meshDimRelToMax,fams,true);
}

ParaMEDMEM::MEDCouplingUMesh *MEDLoader::ReadUMeshFromGroups(const char *fileName, const char *meshName, int meshDimRelToMax, const std::vector<std::string>& grps) throw(INTERP_KERNEL::Exception)
{
  CheckFileForRead(fileName);
  MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> mm=MEDFileMesh::New(fileName,meshName);
  MEDFileMesh *mmPtr(mm);
  MEDFileUMesh *mmuPtr=dynamic_cast<MEDFileUMesh *>(mmPtr);
  if(!mmuPtr)
    {
      std::ostringstream oss; oss << "MEDLoader::ReadUMeshFromGroups : With fileName=\""<< fileName << "\", meshName (the first) =\""<< mm->getName() << "\" exists but it is not an unstructured mesh !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
    return mmuPtr->getGroups(meshDimRelToMax,grps,true);
}

ParaMEDMEM::MEDCouplingFieldDouble *MEDLoader::ReadField(ParaMEDMEM::TypeOfField type, const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
{
  CheckFileForRead(fileName);
  switch(type)
    {
    case ON_CELLS:
      return ReadFieldCell(fileName,meshName,meshDimRelToMax,fieldName,iteration,order);
    case ON_NODES:
      return ReadFieldNode(fileName,meshName,meshDimRelToMax,fieldName,iteration,order);
    case ON_GAUSS_PT:
      return ReadFieldGauss(fileName,meshName,meshDimRelToMax,fieldName,iteration,order);
    case ON_GAUSS_NE:
      return ReadFieldGaussNE(fileName,meshName,meshDimRelToMax,fieldName,iteration,order);
    default:
      throw INTERP_KERNEL::Exception("Type of field specified not managed ! manages are ON_NODES, ON_CELLS, ON_GAUSS_PT or ON_GAUSS_NE !");
    } 
}

std::vector<ParaMEDMEM::MEDCouplingFieldDouble *> MEDLoader::ReadFieldsOnSameMesh(ParaMEDMEM::TypeOfField type, const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName,
                                                                                  const std::vector<std::pair<int,int> >& its) throw(INTERP_KERNEL::Exception)
{
  if(its.empty())
    return std::vector<ParaMEDMEM::MEDCouplingFieldDouble *>();
  CheckFileForRead(fileName);
  std::vector<ParaMEDMEM::MEDCouplingFieldDouble *> ret(its.size());
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> > retSafe(its.size());
  if(its.empty())
    return ret;
  //Retrieving mesh of rank 0 and field on rank 0 too.
  MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> mm=MEDFileMesh::New(fileName,meshName);
  MEDFileMesh *mmPtr(mm);
  MEDFileUMesh *mmuPtr=dynamic_cast<MEDFileUMesh *>(mmPtr);
  if(!mmuPtr)
    throw INTERP_KERNEL::Exception("MEDLoader::ReadFieldsOnSameMesh : only unstructured mesh is managed !");
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m=mmuPtr->getMeshAtLevel(meshDimRelToMax);
  const DataArrayInt *o2n=mmuPtr->getNumberFieldAtLevel(meshDimRelToMax);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m2(m->clone(true));
  if(o2n)
    m2->renumberCells(o2n->begin(),true);
  int i=0;
  for(std::vector<std::pair<int,int> >::const_iterator it=its.begin();it!=its.end();it++,i++)
    {
      MEDCouplingAutoRefCountObjectPtr<MEDFileField1TS> ff=MEDFileField1TS::New(fileName,fieldName,(*it).first,(*it).second);
      MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> retElt=ff->getFieldOnMeshAtLevel(type,m);
      if(o2n)
        retElt->renumberCells(o2n->begin(),true);
      retElt->setMesh(m2);
      retSafe[i]=retElt;
    }
  i=0;
  for(std::vector<std::pair<int,int> >::const_iterator it=its.begin();it!=its.end();it++,i++)
    ret[i]=retSafe[i].retn();
  return ret;
}

std::vector<ParaMEDMEM::MEDCouplingFieldDouble *> MEDLoader::ReadFieldsCellOnSameMesh(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName,
                                                                                            const std::vector<std::pair<int,int> >& its) throw(INTERP_KERNEL::Exception)
{
  return ReadFieldsOnSameMesh(ON_CELLS,fileName,meshName,meshDimRelToMax,fieldName,its);
}

std::vector<ParaMEDMEM::MEDCouplingFieldDouble *> MEDLoader::ReadFieldsNodeOnSameMesh(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName,
                                                                                      const std::vector<std::pair<int,int> >& its) throw(INTERP_KERNEL::Exception)
{
  return ReadFieldsOnSameMesh(ON_NODES,fileName,meshName,meshDimRelToMax,fieldName,its);
}

std::vector<ParaMEDMEM::MEDCouplingFieldDouble *> MEDLoader::ReadFieldsGaussOnSameMesh(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName,
                                                                                       const std::vector<std::pair<int,int> >& its) throw(INTERP_KERNEL::Exception)
{
  return ReadFieldsOnSameMesh(ON_GAUSS_PT,fileName,meshName,meshDimRelToMax,fieldName,its);
}

std::vector<ParaMEDMEM::MEDCouplingFieldDouble *> MEDLoader::ReadFieldsGaussNEOnSameMesh(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName,
                                                                                         const std::vector<std::pair<int,int> >& its) throw(INTERP_KERNEL::Exception)
{
  return ReadFieldsOnSameMesh(ON_GAUSS_NE,fileName,meshName,meshDimRelToMax,fieldName,its);
}

ParaMEDMEM::MEDCouplingFieldDouble *MEDLoader::ReadFieldCell(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileField1TS> ff=MEDFileField1TS::New(fileName,fieldName,iteration,order);
  MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> mm=MEDFileMesh::New(fileName,meshName);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> m=mm->getGenMeshAtLevel(meshDimRelToMax,false);
  MEDFileMesh *mPtr(mm);
  MEDFileUMesh *muPtr=dynamic_cast<MEDFileUMesh *>(mPtr);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=ff->getFieldOnMeshAtLevel(ON_CELLS,m);
  if(muPtr)
    {
      const DataArrayInt *num=muPtr->getNumberFieldAtLevel(meshDimRelToMax);
      if(num)
        ret->renumberCells(num->begin());
    }
  return ret.retn();
}

ParaMEDMEM::MEDCouplingFieldDouble *MEDLoader::ReadFieldNode(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileField1TS> ff=MEDFileField1TS::New(fileName,fieldName,iteration,order);
  MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> mm=MEDFileMesh::New(fileName,meshName);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> m=mm->getGenMeshAtLevel(meshDimRelToMax,false);
  MEDFileMesh *mPtr(mm);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=ff->getFieldOnMeshAtLevel(ON_NODES,m);
  MEDFileUMesh *muPtr=dynamic_cast<MEDFileUMesh *>(mPtr);
  if(ff->getPflsReallyUsed().empty())
    {
      if(muPtr)
        {
          const DataArrayInt *num=muPtr->getNumberFieldAtLevel(meshDimRelToMax);
          if(num)
            ret->renumberCells(num->begin());
        }
    }
  else
    {
      DataArrayInt *pfl=0,*arr2=0;
      MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=ff->getFieldWithProfile(ON_NODES,meshDimRelToMax,mm,pfl);
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> pflSafe(pfl);
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> mp=m->getCellIdsFullyIncludedInNodeIds(pfl->begin(),pfl->end());
      MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> mzip=static_cast<MEDCouplingUMesh *>(m->buildPartAndReduceNodes(mp->begin(),mp->end(),arr2));
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr2Safe(arr2);
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr3=arr2->invertArrayO2N2N2O(mzip->getNumberOfNodes());
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> pflSorted(pflSafe->deepCpy()); pflSorted->sort(true);
      if(!arr3->isEqualWithoutConsideringStr(*pflSorted))
        throw INTERP_KERNEL::Exception("MEDLoader::ReadFieldNode : not implemented yet !");
      if(!arr3->isEqualWithoutConsideringStr(*pflSafe))
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> o2n2=pflSafe->checkAndPreparePermutation();
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> n2o2=o2n2->invertArrayO2N2N2O(o2n2->getNumberOfTuples());
          mzip->renumberNodes(n2o2->begin(),n2o2->getNumberOfTuples());
          arr->setName("");
          ret->setArray(arr);
        }
      ret->setMesh(mzip);
    }
  return ret.retn();
}

ParaMEDMEM::MEDCouplingFieldDouble *MEDLoader::ReadFieldGauss(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileField1TS> ff=MEDFileField1TS::New(fileName,fieldName,iteration,order);
  MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> mm=MEDFileMesh::New(fileName,meshName);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> m=mm->getGenMeshAtLevel(meshDimRelToMax,false);
  MEDFileMesh *mPtr(mm);
  MEDFileUMesh *muPtr=dynamic_cast<MEDFileUMesh *>(mPtr);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=ff->getFieldOnMeshAtLevel(ON_GAUSS_PT,m);
  if(muPtr)
    {
      const DataArrayInt *num=muPtr->getNumberFieldAtLevel(meshDimRelToMax);
      if(num)
        ret->renumberCells(num->begin());
    }
  return ret.retn();
}

ParaMEDMEM::MEDCouplingFieldDouble *MEDLoader::ReadFieldGaussNE(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileField1TS> ff=MEDFileField1TS::New(fileName,fieldName,iteration,order);
  MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> mm=MEDFileMesh::New(fileName,meshName);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> m=mm->getGenMeshAtLevel(meshDimRelToMax,false);
  MEDFileMesh *mPtr(mm);
  MEDFileUMesh *muPtr=dynamic_cast<MEDFileUMesh *>(mPtr);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=ff->getFieldOnMeshAtLevel(ON_GAUSS_NE,m);
  if(muPtr)
    {
      const DataArrayInt *num=muPtr->getNumberFieldAtLevel(meshDimRelToMax);
      if(num)
        ret->renumberCells(num->begin());
    }
  return ret.retn();
}

void MEDLoader::WriteUMesh(const char *fileName, const ParaMEDMEM::MEDCouplingUMesh *mesh, bool writeFromScratch) throw(INTERP_KERNEL::Exception)
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDLoader::WriteUMesh : input mesh is null !");
  int mod=writeFromScratch?2:0;
  MEDCouplingAutoRefCountObjectPtr<MEDFileUMesh> m=MEDFileUMesh::New();
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> mcpy(static_cast<MEDCouplingUMesh *>(mesh->deepCpy()));
  m->setMeshAtLevel(0,mcpy,true);
  m->write(fileName,mod);
}

void MEDLoader::WriteUMeshDep(const char *fileName, const ParaMEDMEM::MEDCouplingUMesh *mesh, bool writeFromScratch) throw(INTERP_KERNEL::Exception)
{
  MEDLoader::WriteUMesh(fileName,mesh,writeFromScratch);
}

void MEDLoader::WriteUMeshesPartition(const char *fileName, const char *meshNameC, const std::vector<const ParaMEDMEM::MEDCouplingUMesh *>& meshes, bool writeFromScratch) throw(INTERP_KERNEL::Exception)
{
  std::string meshName(meshNameC);
  if(meshName.empty())
    throw INTERP_KERNEL::Exception("Trying to write a unstructured mesh with no name ! MED file format needs a not empty mesh name : change 2nd parameter !");
  int status=MEDLoaderBase::getStatusOfFile(fileName);
  if(status!=MEDLoaderBase::EXIST_RW && status!=MEDLoaderBase::NOT_EXIST)
    {
      std::ostringstream oss; oss << "File with name \'" << fileName << "\' has not valid permissions !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  MEDCouplingAutoRefCountObjectPtr<MEDFileUMesh> m=MEDFileUMesh::New();
  m->setGroupsFromScratch(0,meshes,true);
  m->setName(meshNameC);
  int mod=writeFromScratch?2:0;
  m->write(fileName,mod);
}

void MEDLoader::WriteUMeshesPartitionDep(const char *fileName, const char *meshNameC, const std::vector<const ParaMEDMEM::MEDCouplingUMesh *>& meshes, bool writeFromScratch) throw(INTERP_KERNEL::Exception)
{
  WriteUMeshesPartition(fileName,meshNameC,meshes,writeFromScratch);
}

void MEDLoader::WriteUMeshes(const char *fileName, const std::vector<const ParaMEDMEM::MEDCouplingUMesh *>& meshes, bool writeFromScratch) throw(INTERP_KERNEL::Exception)
{
  int mod=writeFromScratch?2:0;
  MEDCouplingAutoRefCountObjectPtr<MEDFileUMesh> m=MEDFileUMesh::New();
  m->setMeshes(meshes,true);
  m->write(fileName,mod);
}

void MEDLoaderNS::writeFieldWithoutReadingAndMappingOfMeshInFile(const char *fileName, const ParaMEDMEM::MEDCouplingFieldDouble *f, bool writeFromScratch)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileField1TS> ff=MEDFileField1TS::New();
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> f2(f->deepCpy());
  const MEDCouplingMesh *m=f2->getMesh();
  const MEDCouplingUMesh *um=dynamic_cast<const MEDCouplingUMesh *>(m);
  MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> mm;
  int mod=writeFromScratch?2:0;
  if(um)
    {
      MEDCouplingAutoRefCountObjectPtr<MEDFileUMesh> mmu=MEDFileUMesh::New();
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> o2n=um->getRenumArrForMEDFileFrmt();
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> n2o=o2n->invertArrayO2N2N2O(o2n->getNumberOfTuples());
      f2->renumberCells(o2n->begin(),false);
      mmu->setMeshAtLevel(0,const_cast<MEDCouplingUMesh *>(static_cast<const MEDCouplingUMesh *>(f2->getMesh())));
      mmu->setRenumFieldArr(0,n2o);
      ff->setFieldNoProfileSBT(f2);
      mmu->write(fileName,mod);
    }
  else
    throw INTERP_KERNEL::Exception("MEDLoaderNS::writeFieldWithoutReadingAndMappingOfMeshInFile : only umeshes are dealed in this API for the moment !");
  ff->write(fileName,0);
}

void MEDLoader::WriteField(const char *fileName, const ParaMEDMEM::MEDCouplingFieldDouble *f, bool writeFromScratch) throw(INTERP_KERNEL::Exception)
{
  if(!f)
    throw INTERP_KERNEL::Exception("MEDLoader::WriteField : input field is NULL !");
  f->checkCoherency();
  int status=MEDLoaderBase::getStatusOfFile(fileName);
  if(status!=MEDLoaderBase::EXIST_RW && status!=MEDLoaderBase::NOT_EXIST)
    {
      std::ostringstream oss; oss << "File with name \'" << fileName << "\' has not valid permissions !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  if(writeFromScratch || (!writeFromScratch && status==MEDLoaderBase::NOT_EXIST))
    {
      MEDLoaderNS::writeFieldWithoutReadingAndMappingOfMeshInFile(fileName,f,true);
    }
  else
    {
      std::vector<std::string> meshNames=GetMeshNames(fileName);
      if(!f->getMesh())
        throw INTERP_KERNEL::Exception("MEDLoader::WriteField : trying to write a field with no mesh !");
      std::string fileNameCpp(f->getMesh()->getName());
      if(std::find(meshNames.begin(),meshNames.end(),fileNameCpp)==meshNames.end())
        MEDLoaderNS::writeFieldWithoutReadingAndMappingOfMeshInFile(fileName,f,false);
      else
        {
          MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> mm=MEDFileMesh::New(fileName,f->getMesh()->getName());
          const MEDFileMesh *mmPtr(mm);
          const MEDFileUMesh *mmuPtr=dynamic_cast<const MEDFileUMesh *>(mmPtr);
          if(!mmuPtr)
            throw INTERP_KERNEL::Exception("MEDLoader::WriteField : only umeshes are supported now !");
          MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> f2(f->deepCpy());
          MEDCouplingUMesh *m=dynamic_cast<MEDCouplingUMesh *>(const_cast<MEDCouplingMesh *>(f2->getMesh()));
          if(!m)
            throw INTERP_KERNEL::Exception("MEDLoader::WriteField : only umesh in input field supported !");
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> o2n=m->getRenumArrForMEDFileFrmt();
          f2->renumberCells(o2n->begin(),false);
          m=static_cast<MEDCouplingUMesh *>(const_cast<MEDCouplingMesh *>(f2->getMesh()));
          MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> mread=mmuPtr->getMeshAtLevel(m->getMeshDimension()-mm->getMeshDimension());
          if(f2->getTypeOfField()!=ON_NODES)
            {
              m->tryToShareSameCoordsPermute(*mread,_EPS_FOR_NODE_COMP);
              DataArrayInt *part=0;
              bool b=mread->areCellsIncludedIn(m,_COMP_FOR_CELL,part);
              MEDCouplingAutoRefCountObjectPtr<DataArrayInt> partSafe(part);
              if(!b)
                {
                  std::ostringstream oss; oss << "MEDLoader::WriteField : The file \""<< fileName << "\" already contains a mesh named \""<< f->getMesh()->getName() << "\" and this mesh in the file is not compatible (a subpart) with the mesh you intend to write ! This is maybe due to a too strict policy ! Try with to lease it by calling SetCompPolicyForCell !";
                  throw INTERP_KERNEL::Exception(oss.str().c_str());
                }
              MEDCouplingAutoRefCountObjectPtr<MEDFileField1TS> f1ts=MEDFileField1TS::New();
              if(part->isIdentity() && part->getNumberOfTuples()==mread->getNumberOfCells())
                f1ts->setFieldNoProfileSBT(f2);
              else
                {
                  part->setName(f1ts->createNewNameOfPfl().c_str());
                  f1ts->setFieldProfile(f2,mm,m->getMeshDimension()-mm->getMeshDimension(),part);
                }
              f1ts->write(fileName,0);
              return ;
            }
          else
            {
              DataArrayInt *part=0;
              bool b=mread->getCoords()->areIncludedInMe(m->getCoords(),_EPS_FOR_NODE_COMP,part);
              MEDCouplingAutoRefCountObjectPtr<DataArrayInt> partSafe(part);
              if(!b)
                {
                  std::ostringstream oss; oss << "MEDLoader::WriteField : The file \""<< fileName << "\" already contains a mesh named \""<< f->getMesh()->getName() << "\" and this mesh in the file is not compatible (a subpart regarding nodes) with the mesh you intend to write ! This is maybe due to a too strict epsilon ! Try with to lease it by calling SetEpsilonForNodeComp !";
                  throw INTERP_KERNEL::Exception(oss.str().c_str());
                }
              MEDCouplingAutoRefCountObjectPtr<MEDFileField1TS> f1ts=MEDFileField1TS::New();
              if(part->isIdentity() && part->getNumberOfTuples()==mread->getNumberOfNodes())
                f1ts->setFieldNoProfileSBT(f2);
              else
                {
                  part->setName(f1ts->createNewNameOfPfl().c_str());
                  f1ts->setFieldProfile(f2,mm,m->getMeshDimension()-mm->getMeshDimension(),part);
                }
              f1ts->write(fileName,0);
            }
        }
    }
}

void MEDLoader::WriteFieldDep(const char *fileName, const ParaMEDMEM::MEDCouplingFieldDouble *f, bool writeFromScratch) throw(INTERP_KERNEL::Exception)
{
  WriteField(fileName,f,writeFromScratch);
}

void MEDLoader::WriteFieldUsingAlreadyWrittenMesh(const char *fileName, const ParaMEDMEM::MEDCouplingFieldDouble *f) throw(INTERP_KERNEL::Exception)
{
  if(!f)
    throw INTERP_KERNEL::Exception("MEDLoader::WriteFieldUsingAlreadyWrittenMesh : input field is null !");
  f->checkCoherency();
  int status=MEDLoaderBase::getStatusOfFile(fileName);
  if(status!=MEDLoaderBase::EXIST_RW)
    {
      std::ostringstream oss; oss << "File with name \'" << fileName << "\' has not valid permissions or not exists !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  MEDCouplingAutoRefCountObjectPtr<MEDFileField1TS> f1ts=MEDFileField1TS::New();
  MEDCouplingUMesh *m=dynamic_cast<MEDCouplingUMesh *>(const_cast<MEDCouplingMesh *>(f->getMesh()));
  if(!m)
    throw INTERP_KERNEL::Exception("MEDLoader::WriteFieldUsingAlreadyWrittenMesh : only umesh in input field supported !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> o2n=m->getRenumArrForMEDFileFrmt();
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> f2(f->deepCpy());
  f2->renumberCells(o2n->begin(),false);
  f1ts->setFieldNoProfileSBT(f2);
  f1ts->write(fileName,0);
}
