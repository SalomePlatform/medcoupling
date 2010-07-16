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

#include "MEDLoader.hxx"
#include "MEDLoaderBase.hxx"
#include "CellModel.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingFieldDouble.hxx"

extern "C"
{
#include "med.h"
}

#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <numeric>

med_geometrie_element typmai[MED_NBR_GEOMETRIE_MAILLE+2] = { MED_POINT1,
                                                             MED_SEG2,
                                                             MED_SEG3,
                                                             MED_TRIA3,
                                                             MED_TRIA6,
                                                             MED_QUAD4,
                                                             MED_QUAD8,
                                                             MED_TETRA4,
                                                             MED_TETRA10,
                                                             MED_HEXA8,
                                                             MED_HEXA20,
                                                             MED_PENTA6,
                                                             MED_PENTA15,
                                                             MED_PYRA5,
                                                             MED_PYRA13,
                                                             MED_POLYGONE,
                                                             MED_POLYEDRE };

med_geometrie_element typmainoeud[1] = { MED_NONE };

INTERP_KERNEL::NormalizedCellType typmai2[MED_NBR_GEOMETRIE_MAILLE+2] = { INTERP_KERNEL::NORM_ERROR,
                                                                          INTERP_KERNEL::NORM_SEG2,
                                                                          INTERP_KERNEL::NORM_SEG3,
                                                                          INTERP_KERNEL::NORM_TRI3,
                                                                          INTERP_KERNEL::NORM_TRI6,
                                                                          INTERP_KERNEL::NORM_QUAD4,
                                                                          INTERP_KERNEL::NORM_QUAD8,
                                                                          INTERP_KERNEL::NORM_TETRA4,
                                                                          INTERP_KERNEL::NORM_TETRA10,
                                                                          INTERP_KERNEL::NORM_HEXA8,
                                                                          INTERP_KERNEL::NORM_HEXA20,
                                                                          INTERP_KERNEL::NORM_PENTA6,
                                                                          INTERP_KERNEL::NORM_PENTA15,
                                                                          INTERP_KERNEL::NORM_PYRA5,
                                                                          INTERP_KERNEL::NORM_PYRA13,
                                                                          INTERP_KERNEL::NORM_POLYGON,
                                                                          INTERP_KERNEL::NORM_POLYHED };

med_geometrie_element typmai3[32] = { MED_POINT1,//0
                                      MED_SEG2,//1
                                      MED_SEG3,//2
                                      MED_TRIA3,//3
                                      MED_QUAD4,//4
                                      MED_POLYGONE,//5
                                      MED_TRIA6,//6
                                      MED_NONE,//7
                                      MED_QUAD8,//8
                                      MED_NONE,//9
                                      MED_NONE,//10
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
                                      MED_NONE,//22
                                      MED_PYRA13,//23
                                      MED_NONE,//24
                                      MED_PENTA15,//25
                                      MED_NONE,//26
                                      MED_NONE,//27
                                      MED_NONE,//28
                                      MED_NONE,//29
                                      MED_HEXA20,//30
                                      MED_POLYEDRE//31
};

using namespace ParaMEDMEM;

namespace MEDLoaderNS
{
  class FieldPerTypeCopier
  {
  public:
    FieldPerTypeCopier(double *ptr):_ptr(ptr) { }
    void operator()(const MEDLoader::MEDFieldDoublePerCellType& elt) { _ptr=std::copy(elt.getArray(),elt.getArray()+elt.getNbOfValues(),_ptr); }
  private:
    double *_ptr;
  };
 
  class ConnReaderML
  {
  public:
    ConnReaderML(const int *c, int val):_conn(c),_val(val) { }
    bool operator() (const int& pos) { return _conn[pos]!=_val; }
  private:
    const int *_conn;
    int _val;
  };
  
  std::string buildStringFromFortran(const char *expr, int lgth);
  std::vector<std::string> getMeshNamesFid(med_idt fid);
  void readFieldDoubleDataInMedFile(const char *fileName, const char *meshName, const char *fieldName, std::list<MEDLoader::MEDFieldDoublePerCellType>& field,
                                    int iteration, int order, ParaMEDMEM::TypeOfField typeOfOutField, double& time);
  std::vector<int> getIdsFromFamilies(const char *fileName, const char *meshName, const std::vector<std::string>& fams);
  std::vector<int> getIdsFromGroups(const char *fileName, const char *meshName, const std::vector<std::string>& grps);
  med_int getIdFromMeshName(med_idt fid, const char *meshName, std::string& trueMeshName) throw(INTERP_KERNEL::Exception);
  void dispatchElems(int nbOfElemCell, int nbOfElemFace, int& nbOfElem, med_entite_maillage& whichEntity);
  void readUMeshDataInMedFile(med_idt fid, med_int meshId, double *&coords, int& nCoords, int& spaceDim, std::list<MEDLoader::MEDConnOfOneElemType>& conn);
  int buildMEDSubConnectivityOfOneType(const DataArrayInt *conn, const DataArrayInt *connIndex, const DataArrayInt *families, INTERP_KERNEL::NormalizedCellType type,
                                       std::vector<int>& conn4MEDFile, std::vector<int>& connIndex4MEDFile, std::vector<int>& connIndexRk24MEDFile,
                                       std::vector<int>& fam4MEDFile);
  MEDCouplingUMesh *readUMeshFromFileLev1(const char *fileName, const char *meshName, int meshDimRelToMax, const std::vector<int>& ids,
                                          const std::vector<INTERP_KERNEL::NormalizedCellType>& typesToKeep, unsigned& meshDimExtract) throw(INTERP_KERNEL::Exception);
  void tradMEDFileCoreFrmt2MEDCouplingUMesh(const std::list<MEDLoader::MEDConnOfOneElemType>& medConnFrmt,
                                            DataArrayInt* &conn,
                                            DataArrayInt* &connIndex,
                                            const std::vector<int>& familiesToKeep);
  ParaMEDMEM::DataArrayDouble *buildArrayFromRawData(const std::list<MEDLoader::MEDFieldDoublePerCellType>& fieldPerType);
  int buildMEDSubConnectivityOfOneTypesPolyg(const DataArrayInt *conn, const DataArrayInt *connIndex, const DataArrayInt *families,
                                             std::vector<int>& conn4MEDFile, std::vector<int>& connIndex4MEDFile, std::vector<int>& fam4MEDFile);
  int buildMEDSubConnectivityOfOneTypesPolyh(const DataArrayInt *conn, const DataArrayInt *connIndex, const DataArrayInt *families, std::vector<int>& conn4MEDFile,
                                             std::vector<int>& connIndex4MEDFile, std::vector<int>& connIndexRk24MEDFile,
                                             std::vector<int>& fam4MEDFile);
  int buildMEDSubConnectivityOfOneTypeStaticTypes(const DataArrayInt *conn, const DataArrayInt *connIndex, const DataArrayInt *families,
                                                  INTERP_KERNEL::NormalizedCellType type, std::vector<int>& conn4MEDFile, std::vector<int>& fam4MEDFile);
  ParaMEDMEM::MEDCouplingFieldDouble *readFieldDoubleLev1(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order,
                                                          ParaMEDMEM::TypeOfField typeOfOutField);
  void appendFieldDirectly(const char *fileName, ParaMEDMEM::MEDCouplingFieldDouble *f);
  void prepareCellFieldDoubleForWriting(const ParaMEDMEM::MEDCouplingFieldDouble *f, std::list<MEDLoader::MEDFieldDoublePerCellType>& split);
  void writeUMeshDirectly(const char *fileName, ParaMEDMEM::MEDCouplingUMesh *mesh, const DataArrayInt *families, bool forceFromScratch);
  void writeUMeshesDirectly(const char *fileName, const char *meshName, const std::vector<ParaMEDMEM::MEDCouplingUMesh *>& meshes, bool forceFromScratch);
  void writeFieldAndMeshDirectly(const char *fileName, ParaMEDMEM::MEDCouplingFieldDouble *f, bool forceFromScratch);
}

const char WHITE_SPACES[]=" \n";

/*!
 * @param lgth is the size of fam tab. For classical types conn is size of 'lgth'*number_of_nodes_in_type.
 * @param index is optionnal only for polys. Set it to 0 if it is not the case.
 * @param connLgth is the size of conn in the case of poly. Unsued if it is not the case.
 */
MEDLoader::MEDConnOfOneElemType::MEDConnOfOneElemType(INTERP_KERNEL::NormalizedCellType type, int *conn, int *index, int *fam, int lgth, int connLgth):_lgth(lgth),_fam(fam),
                                                                                                                                                       _conn(conn),_index(index),
                                                                                                                                                       _global(0),_type(type),
                                                                                                                                                       _conn_lgth(connLgth)
{
}

void MEDLoader::MEDConnOfOneElemType::setGlobal(int *global)
{
  if(_global!=global)
    {
      if(_global)
        delete [] _global;
      _global=global;
    }
}

void MEDLoader::MEDConnOfOneElemType::releaseArray()
{
  delete [] _fam;
  delete [] _conn;
  delete [] _index;
  delete [] _global;
}

MEDLoader::MEDFieldDoublePerCellType::MEDFieldDoublePerCellType(INTERP_KERNEL::NormalizedCellType type, double *values, int ncomp, int ntuple):_ntuple(ntuple),_ncomp(ncomp),_values(values),_type(type)
{
}

void MEDLoader::MEDFieldDoublePerCellType::releaseArray()
{
  delete [] _values;
}


std::string MEDLoaderNS::buildStringFromFortran(const char *expr, int lgth)
{
  std::string ret(expr,lgth);
  std::string whiteSpaces(WHITE_SPACES);
  std::size_t lgthReal=strlen(ret.c_str());
  std::string ret2=ret.substr(0,lgthReal);
  std::size_t found=ret2.find_last_not_of(whiteSpaces);
  if (found!=std::string::npos)
    ret2.erase(found+1);
  else
    ret2.clear();//ret is all whitespace
  return ret2;
}

std::vector<std::string> MEDLoaderNS::getMeshNamesFid(med_idt fid)
{
  med_maillage type_maillage;
  char maillage_description[MED_TAILLE_DESC+1];
  med_int dim;
  char nommaa[MED_TAILLE_NOM+1];
  med_int n=MEDnMaa(fid);
  std::vector<std::string> ret(n);
  for(int i=0;i<n;i++)
    {
      MEDmaaInfo(fid,i+1,nommaa,&dim,&type_maillage,maillage_description);
      std::string cur=buildStringFromFortran(nommaa,sizeof(nommaa));
      ret[i]=cur;
    }
  return ret;
}

std::vector<std::string> MEDLoader::GetMeshNames(const char *fileName)
{
  med_idt fid=MEDouvrir((char *)fileName,MED_LECTURE);
  std::vector<std::string> ret=MEDLoaderNS::getMeshNamesFid(fid);
  MEDfermer(fid);
  return ret;
}

std::vector<std::string> MEDLoader::GetMeshFamilyNames(const char *fileName, const char *meshName)
{
  med_idt fid=MEDouvrir((char *)fileName,MED_LECTURE);
  med_int nfam=MEDnFam(fid,(char *)meshName);
  std::vector<std::string> ret(nfam);
  char nomfam[MED_TAILLE_NOM+1];
  med_int numfam;
  for(int i=0;i<nfam;i++)
    {
      int ngro=MEDnGroupe(fid,(char *)meshName,i+1);
      med_int natt=MEDnAttribut(fid,(char *)meshName,i+1);
      med_int *attide=new int[natt];
      med_int *attval=new int[natt];
      char *attdes=new char[MED_TAILLE_DESC*natt+1];
      char *gro=new char[MED_TAILLE_LNOM*ngro+1];
      MEDfamInfo(fid,(char *)meshName,i+1,nomfam,&numfam,attide,attval,attdes,&natt,gro,&ngro);
      std::string cur=MEDLoaderNS::buildStringFromFortran(nomfam,sizeof(nomfam));
      ret[i]=cur;
      delete [] attdes;
      delete [] gro;
      delete [] attide;
      delete [] attval;
    }
  MEDfermer(fid);
  return ret;
}
  
std::vector<std::string> MEDLoader::GetMeshGroupsNames(const char *fileName, const char *meshName)
{
  med_idt fid=MEDouvrir((char *)fileName,MED_LECTURE);
  med_int nfam=MEDnFam(fid,(char *)meshName);
  std::vector<std::string> ret;
  char nomfam[MED_TAILLE_NOM+1];
  med_int numfam;
  for(int i=0;i<nfam;i++)
    {
      int ngro=MEDnGroupe(fid,(char *)meshName,i+1);
      med_int natt=MEDnAttribut(fid,(char *)meshName,i+1);
      med_int *attide=new int[natt];
      med_int *attval=new int[natt];
      char *attdes=new char[MED_TAILLE_DESC*natt+1];
      char *gro=new char[MED_TAILLE_LNOM*ngro+1];
      MEDfamInfo(fid,(char *)meshName,i+1,nomfam,&numfam,attide,attval,attdes,&natt,gro,&ngro);
      for(int j=0;j<ngro;j++)
        {
          std::string cur=MEDLoaderNS::buildStringFromFortran(gro+j*MED_TAILLE_LNOM,MED_TAILLE_LNOM);
          if(std::find(ret.begin(),ret.end(),cur)==ret.end())
            ret.push_back(cur);
        }
      delete [] attdes;
      delete [] gro;
      delete [] attide;
      delete [] attval;
    }
  MEDfermer(fid);
  return ret;
}

std::vector<std::string> MEDLoader::GetCellFieldNamesOnMesh(const char *fileName, const char *meshName)
{
  std::vector<std::string> ret;
  med_idt fid=MEDouvrir((char *)fileName,MED_LECTURE);
  med_int nbFields=MEDnChamp(fid,0);
  //
  med_type_champ typcha;
  //med_int nbpdtnor=0,pflsize,*pflval,lnsize;
  med_int ngauss=0;
  med_int numdt=0,numo=0,nbrefmaa;
  med_float dt=0.0;
  med_booleen local;
  //char pflname[MED_TAILLE_NOM+1]="";
  //char locname[MED_TAILLE_NOM+1]="";
  char *maa_ass=MEDLoaderBase::buildEmptyString(MED_TAILLE_NOM);
  char *dt_unit=MEDLoaderBase::buildEmptyString(MED_TAILLE_PNOM);
  char *nomcha=MEDLoaderBase::buildEmptyString(MED_TAILLE_NOM);
  //
  for(int i=0;i<nbFields;i++)
    {
      med_int ncomp=MEDnChamp(fid,i+1);
      char *comp=new char[ncomp*MED_TAILLE_PNOM+1];
      char *unit=new char[ncomp*MED_TAILLE_PNOM+1];
      MEDchampInfo(fid,i+1,nomcha,&typcha,comp,unit,ncomp);
      std::string curFieldName=MEDLoaderNS::buildStringFromFortran(nomcha,MED_TAILLE_NOM+1);
      delete [] comp;
      delete [] unit;
      bool found=false;
      for(int j=0;j<MED_NBR_GEOMETRIE_MAILLE+2 && !found;j++)
        {
          med_int nbPdt=MEDnPasdetemps(fid,nomcha,MED_MAILLE,typmai[j]);
          if(nbPdt>0)
            {
              MEDpasdetempsInfo(fid,nomcha,MED_MAILLE,typmai[j],1, &ngauss, &numdt, &numo, dt_unit,&dt, maa_ass, &local, &nbrefmaa);
              std::string curMeshName=MEDLoaderNS::buildStringFromFortran(maa_ass,MED_TAILLE_NOM+1);
              if(curMeshName==meshName)
                {
                  found=true;
                  ret.push_back(curFieldName);
                }
            }
        }
    }
  delete [] maa_ass;
  delete [] dt_unit;
  delete [] nomcha;
  MEDfermer(fid);
  return ret;
}

std::vector<std::string> MEDLoader::GetNodeFieldNamesOnMesh(const char *fileName, const char *meshName)
{
  std::vector<std::string> ret;
  med_idt fid=MEDouvrir((char *)fileName,MED_LECTURE);
  med_int nbFields=MEDnChamp(fid,0);
  //
  med_type_champ typcha;
  med_int ngauss=0;
  med_int numdt=0,numo=0,nbrefmaa;
  med_float dt=0.0;
  med_booleen local;
  char *maa_ass=MEDLoaderBase::buildEmptyString(MED_TAILLE_NOM);
  char *dt_unit=MEDLoaderBase::buildEmptyString(MED_TAILLE_PNOM);
  char *nomcha=MEDLoaderBase::buildEmptyString(MED_TAILLE_NOM);
  //
  for(int i=0;i<nbFields;i++)
    {
      med_int ncomp=MEDnChamp(fid,i+1);
      char *comp=new char[ncomp*MED_TAILLE_PNOM+1];
      char *unit=new char[ncomp*MED_TAILLE_PNOM+1];
      MEDchampInfo(fid,i+1,nomcha,&typcha,comp,unit,ncomp);
      std::string curFieldName=MEDLoaderNS::buildStringFromFortran(nomcha,MED_TAILLE_NOM+1);
      delete [] comp;
      delete [] unit;
      bool found=false;
      med_int nbPdt=MEDnPasdetemps(fid,nomcha,MED_NOEUD,MED_NONE);
      if(nbPdt>0)
        {
          MEDpasdetempsInfo(fid,nomcha,MED_NOEUD,MED_NONE,1, &ngauss, &numdt, &numo, dt_unit,&dt, maa_ass, &local, &nbrefmaa);
          std::string curMeshName=MEDLoaderNS::buildStringFromFortran(maa_ass,MED_TAILLE_NOM+1);
          if(curMeshName==meshName)
            {
              found=true;
              ret.push_back(curFieldName);
            }
        }
    }
  delete [] maa_ass;
  delete [] dt_unit;
  delete [] nomcha;
  MEDfermer(fid);
  return ret;
}

std::vector< std::pair<int,int> > MEDLoader::GetFieldIterations(ParaMEDMEM::TypeOfField type, const char *fileName, const char *meshName, const char *fieldName)
{
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

std::vector< std::pair<int,int> > MEDLoader::GetCellFieldIterations(const char *fileName, const char *meshName, const char *fieldName)
{
  std::string meshNameCpp(meshName);
  std::vector< std::pair<int,int> > ret;
  med_idt fid=MEDouvrir((char *)fileName,MED_LECTURE);
  med_int nbFields=MEDnChamp(fid,0);
  //
  med_type_champ typcha;
  med_int ngauss=0;
  med_int numdt=0,numo=0,nbrefmaa;
  med_float dt=0.0;
  med_booleen local;
  char *maa_ass=MEDLoaderBase::buildEmptyString(MED_TAILLE_NOM);
  char *dt_unit=MEDLoaderBase::buildEmptyString(MED_TAILLE_PNOM);
  char *nomcha=MEDLoaderBase::buildEmptyString(MED_TAILLE_NOM);
  //
  for(int i=0;i<nbFields;i++)
    {
      med_int ncomp=MEDnChamp(fid,i+1);
      char *comp=new char[ncomp*MED_TAILLE_PNOM+1];
      char *unit=new char[ncomp*MED_TAILLE_PNOM+1];
      MEDchampInfo(fid,i+1,nomcha,&typcha,comp,unit,ncomp);
      std::string curFieldName=MEDLoaderNS::buildStringFromFortran(nomcha,MED_TAILLE_NOM+1);
      delete [] comp;
      delete [] unit;
      if(curFieldName==fieldName)
        {
          bool found=false;
          for(int j=0;j<MED_NBR_GEOMETRIE_MAILLE+2 && !found;j++)
            {
              med_int nbPdt=MEDnPasdetemps(fid,nomcha,MED_MAILLE,typmai[j]);
              for(int k=0;k<nbPdt;k++)
                {
                  MEDpasdetempsInfo(fid,nomcha,MED_MAILLE,typmai[j],k+1, &ngauss, &numdt, &numo, dt_unit,&dt, maa_ass, &local, &nbrefmaa);
                  std::string maa_ass_cpp(maa_ass);
                  if(meshNameCpp==maa_ass_cpp)
                    {
                      found=true;
                      ret.push_back(std::make_pair(numdt,numo));
                    }
                }
            }
        }
    }
  delete [] maa_ass;
  delete [] dt_unit;
  delete [] nomcha;
  MEDfermer(fid);
  return ret;
}

std::vector< std::pair<int,int> > MEDLoader::GetNodeFieldIterations(const char *fileName, const char *meshName, const char *fieldName)
{
  std::string meshNameCpp(meshName);
  std::vector< std::pair<int,int> > ret;
  med_idt fid=MEDouvrir((char *)fileName,MED_LECTURE);
  med_int nbFields=MEDnChamp(fid,0);
  //
  med_type_champ typcha;
  med_int ngauss=0;
  med_int numdt=0,numo=0,nbrefmaa;
  med_float dt=0.0;
  med_booleen local;
  char *maa_ass=MEDLoaderBase::buildEmptyString(MED_TAILLE_NOM);
  char *dt_unit=MEDLoaderBase::buildEmptyString(MED_TAILLE_PNOM);
  char *nomcha=MEDLoaderBase::buildEmptyString(MED_TAILLE_NOM);
  //
  for(int i=0;i<nbFields;i++)
    {
      med_int ncomp=MEDnChamp(fid,i+1);
      char *comp=new char[ncomp*MED_TAILLE_PNOM+1];
      char *unit=new char[ncomp*MED_TAILLE_PNOM+1];
      MEDchampInfo(fid,i+1,nomcha,&typcha,comp,unit,ncomp);
      std::string curFieldName=MEDLoaderNS::buildStringFromFortran(nomcha,MED_TAILLE_NOM+1);
      delete [] comp;
      delete [] unit;
      if(curFieldName==fieldName)
        {
          med_int nbPdt=MEDnPasdetemps(fid,nomcha,MED_NOEUD,MED_NONE);
          for(int k=0;k<nbPdt;k++)
            {
              MEDpasdetempsInfo(fid,nomcha,MED_NOEUD,MED_NONE,k+1, &ngauss, &numdt, &numo, dt_unit,&dt, maa_ass, &local, &nbrefmaa);
               std::string maa_ass_cpp(maa_ass);
               if(meshNameCpp==maa_ass_cpp)
                 {
                   ret.push_back(std::make_pair(numdt,numo));
                 }
            }
        }
    }
  delete [] maa_ass;
  delete [] dt_unit;
  delete [] nomcha;
  MEDfermer(fid);
  return ret;
}

void MEDLoaderNS::readFieldDoubleDataInMedFile(const char *fileName, const char *meshName, const char *fieldName, std::list<MEDLoader::MEDFieldDoublePerCellType>& field,
                                               int iteration, int order, ParaMEDMEM::TypeOfField typeOfOutField, double& time)
{
  time=0.;
  med_idt fid=MEDouvrir((char *)fileName,MED_LECTURE);
  med_int nbFields=MEDnChamp(fid,0);
  //
  med_type_champ typcha;
  char nomcha[MED_TAILLE_NOM+1]="";
  char pflname [MED_TAILLE_NOM+1]="";
  char locname [MED_TAILLE_NOM+1]="";
  std::map<ParaMEDMEM::TypeOfField, med_entite_maillage> tabEnt;
  std::map<ParaMEDMEM::TypeOfField, med_geometrie_element *> tabType;
  std::map<ParaMEDMEM::TypeOfField, int> tabTypeLgth;
  tabEnt[ON_CELLS]=MED_MAILLE;
  tabType[ON_CELLS]=typmai;
  tabTypeLgth[ON_CELLS]=MED_NBR_GEOMETRIE_MAILLE+2;
  tabEnt[ON_NODES]=MED_NOEUD;
  tabType[ON_NODES]=typmainoeud;
  tabTypeLgth[ON_NODES]=1;
  //
  for(int i=0;i<nbFields;i++)
    {
      med_int ncomp=MEDnChamp(fid,i+1);
      char *comp=new char[ncomp*MED_TAILLE_PNOM+1];
      char *unit=new char[ncomp*MED_TAILLE_PNOM+1];
      MEDchampInfo(fid,i+1,nomcha,&typcha,comp,unit,ncomp);
      std::string curFieldName=buildStringFromFortran(nomcha,MED_TAILLE_NOM+1);
      delete [] comp;
      delete [] unit;
      if(curFieldName==fieldName)
        {
          bool found=false;
          for(int j=0;j<tabTypeLgth[typeOfOutField] && !found;j++)
            {
              med_int nbPdt=MEDnPasdetemps(fid,nomcha,tabEnt[typeOfOutField],typmai[j]);
              if(nbPdt>0)
                {
                  int nval=MEDnVal(fid,(char *)fieldName,tabEnt[typeOfOutField],tabType[typeOfOutField][j],iteration,order,(char *)meshName,MED_COMPACT);
                  double *valr=new double[ncomp*nval];
                  //
                  med_int ngauss=0;
                  med_int numdt=0,numo=0,nbrefmaa;
                  char *dt_unit=MEDLoaderBase::buildEmptyString(MED_TAILLE_PNOM);
                  char *maa_ass=MEDLoaderBase::buildEmptyString(MED_TAILLE_NOM);
                  med_float dt=0.0;
                  med_booleen local;
                  med_int nbPdt=MEDnPasdetemps(fid,(char *)fieldName,tabEnt[typeOfOutField],tabType[typeOfOutField][j]);
                  bool found2=false;
                  for(int k=0;k<nbPdt && !found2;k++)
                    {
                      MEDpasdetempsInfo(fid,(char *)fieldName,tabEnt[typeOfOutField],tabType[typeOfOutField][j],k+1,&ngauss,
                                        &numdt,&numo,dt_unit,&dt,maa_ass,&local,&nbrefmaa);
                      found2=(numdt==iteration && numo==order);
                      if(found2)
                        time=dt;
                    }
                  MEDchampLire(fid,(char *)meshName,(char *)fieldName,(unsigned char*)valr,MED_FULL_INTERLACE,MED_ALL,locname,
                               pflname,MED_COMPACT,tabEnt[typeOfOutField],tabType[typeOfOutField][j],iteration,order);
                  field.push_back(MEDLoader::MEDFieldDoublePerCellType(typmai2[j],valr,ncomp,nval));
                  delete [] dt_unit;
                  delete [] maa_ass;
                }
            }
        }
    }
  MEDfermer(fid);
}

std::vector<int> MEDLoaderNS::getIdsFromFamilies(const char *fileName, const char *meshName, const std::vector<std::string>& fams)
{
  std::vector<int> ret;
  med_idt fid=MEDouvrir((char *)fileName,MED_LECTURE);
  med_int nfam=MEDnFam(fid,(char *)meshName);
  char nomfam[MED_TAILLE_NOM+1];
  med_int numfam;
  for(int i=0;i<nfam;i++)
    {
      int ngro=MEDnGroupe(fid,(char *)meshName,i+1);
      med_int natt=MEDnAttribut(fid,(char *)meshName,i+1);
      med_int *attide=new int[natt];
      med_int *attval=new int[natt];
      char *attdes=new char[MED_TAILLE_DESC*natt+1];
      char *gro=new char[MED_TAILLE_LNOM*ngro+1];
      MEDfamInfo(fid,(char *)meshName,i+1,nomfam,&numfam,attide,attval,attdes,&natt,gro,&ngro);
      std::string cur=buildStringFromFortran(nomfam,sizeof(nomfam));
      if(std::find(fams.begin(),fams.end(),cur)!=fams.end())
        ret.push_back(numfam);
      delete [] attdes;
      delete [] gro;
      delete [] attide;
      delete [] attval;
    }
  MEDfermer(fid);
  return ret;
}

std::vector<int> MEDLoaderNS::getIdsFromGroups(const char *fileName, const char *meshName, const std::vector<std::string>& grps)
{
  std::vector<int> ret;
  med_idt fid=MEDouvrir((char *)fileName,MED_LECTURE);
  med_int nfam=MEDnFam(fid,(char *)meshName);
  char nomfam[MED_TAILLE_NOM+1];
  med_int numfam;
  for(int i=0;i<nfam;i++)
    {
      int ngro=MEDnGroupe(fid,(char *)meshName,i+1);
      med_int natt=MEDnAttribut(fid,(char *)meshName,i+1);
      med_int *attide=new int[natt];
      med_int *attval=new int[natt];
      char *attdes=new char[MED_TAILLE_DESC*natt+1];
      char *gro=new char[MED_TAILLE_LNOM*ngro+1];
      MEDfamInfo(fid,(char *)meshName,i+1,nomfam,&numfam,attide,attval,attdes,&natt,gro,&ngro);
      std::string cur=buildStringFromFortran(nomfam,sizeof(nomfam));
      for(int j=0;j<ngro;j++)
        {
          std::string cur=buildStringFromFortran(gro+j*MED_TAILLE_LNOM,MED_TAILLE_LNOM);
          if(std::find(grps.begin(),grps.end(),cur)!=grps.end())
            {
              ret.push_back(numfam);
              break;
            }
        }
      delete [] attdes;
      delete [] gro;
      delete [] attide;
      delete [] attval;
    }
  MEDfermer(fid);
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

/*!
 * This methods allows to merger all entities and to considerate only cell types.
 */
void MEDLoaderNS::dispatchElems(int nbOfElemCell, int nbOfElemFace, int& nbOfElem, med_entite_maillage& whichEntity)
{
  if(nbOfElemCell>=nbOfElemFace)
    {
      whichEntity=MED_MAILLE;
      nbOfElem=nbOfElemCell;
    }
  else
    {
      whichEntity=MED_FACE;
      nbOfElem=nbOfElemFace;
    }
}

void MEDLoaderNS::readUMeshDataInMedFile(med_idt fid, med_int meshId, double *&coords, int& nCoords, int& spaceDim, std::list<MEDLoader::MEDConnOfOneElemType>& conn)
{
  char nommaa[MED_TAILLE_NOM+1];
  char maillage_description[MED_TAILLE_DESC+1];
  char comp[3*MED_TAILLE_PNOM+1];
  char unit[3*MED_TAILLE_PNOM+1];
  med_maillage type_maillage;
  med_int Mdim;
  MEDmaaInfo(fid,meshId,nommaa,&Mdim,&type_maillage,maillage_description);
  spaceDim=(int)Mdim;
  nCoords=MEDnEntMaa(fid,nommaa,MED_COOR,MED_NOEUD,(med_geometrie_element)0,(med_connectivite)0);
  coords=new double[nCoords*spaceDim];
  med_repere repere;
  MEDcoordLire(fid,nommaa,Mdim,coords,MED_FULL_INTERLACE,MED_ALL,NULL,0,&repere,comp,unit);
  med_booleen inoele, inuele;
  for(int i=0;i<MED_NBR_GEOMETRIE_MAILLE;i++)
    {
      med_geometrie_element curMedType=typmai[i];
      med_entite_maillage whichEntity;
      int curNbOfElemM=MEDnEntMaa(fid,nommaa,MED_CONN,MED_MAILLE,curMedType,MED_NOD);
      int curNbOfElemF=MEDnEntMaa(fid,nommaa,MED_CONN,MED_FACE,curMedType,MED_NOD);
      int curNbOfElem;
      MEDLoaderNS::dispatchElems(curNbOfElemM,curNbOfElemF,curNbOfElem,whichEntity);
      if(curNbOfElem>0)
        {
          int *connTab=new int[(curMedType%100)*curNbOfElem];
          int *fam=new int[curNbOfElem];
          MEDLoader::MEDConnOfOneElemType elem(typmai2[i],connTab,0,fam,curNbOfElem,-1);
          int *tmp=new int[curNbOfElem];
          char *noms=new char[MED_TAILLE_PNOM*curNbOfElem+1];
          MEDelementsLire(fid,nommaa,Mdim,connTab,MED_FULL_INTERLACE,noms,&inoele,tmp,&inuele,fam,curNbOfElem,whichEntity,curMedType,MED_NOD);
          delete [] tmp;
          delete [] noms;
          //trying to read global numbering
          int *globArr=new int[curNbOfElem];
          if(MEDglobalNumLire(fid,nommaa,globArr,curNbOfElem,whichEntity,curMedType)==0)
            elem.setGlobal(globArr);
          else
            delete [] globArr;
          conn.push_back(elem);
        }
    }
  int curNbOfPolyElem;
  int curNbOfPolyElemM=MEDnEntMaa(fid,nommaa,MED_CONN,MED_MAILLE,MED_POLYGONE,MED_NOD);
  int curNbOfPolyElemF=MEDnEntMaa(fid,nommaa,MED_CONN,MED_FACE,MED_POLYGONE,MED_NOD);
  med_entite_maillage whichPolyEntity;
  MEDLoaderNS::dispatchElems(curNbOfPolyElemM,curNbOfPolyElemF,curNbOfPolyElem,whichPolyEntity);
  if(curNbOfPolyElem>0)
    {
      med_int arraySize;
      MEDpolygoneInfo(fid,nommaa,whichPolyEntity,MED_NOD,&arraySize);
      int *index=new int[curNbOfPolyElem+1];
      int *locConn=new int[arraySize];
      int *fam=new int[curNbOfPolyElem];
      MEDLoader::MEDConnOfOneElemType elem(INTERP_KERNEL::NORM_POLYGON,locConn,index,fam,curNbOfPolyElem,arraySize);
      MEDpolygoneConnLire(fid,nommaa,index,curNbOfPolyElem+1,locConn,whichPolyEntity,MED_NOD);
      MEDfamLire(fid,nommaa,fam,curNbOfPolyElem,MED_MAILLE,MED_POLYGONE);
      conn.push_back(elem);
    }
  curNbOfPolyElem=MEDnEntMaa(fid,nommaa,MED_CONN,MED_MAILLE,MED_POLYEDRE,MED_NOD);
  if(curNbOfPolyElem>0)
    {
      med_int indexFaceLgth,connFaceLgth;
      MEDpolyedreInfo(fid,nommaa,MED_NOD,&indexFaceLgth,&connFaceLgth);
      int *index=new int[curNbOfPolyElem+1];
      int *indexFace=new int[indexFaceLgth];
      int *locConn=new int[connFaceLgth];
      int *fam=new int[curNbOfPolyElem];
      MEDpolyedreConnLire(fid,nommaa,index,curNbOfPolyElem+1,indexFace,indexFaceLgth,locConn,MED_NOD);
      MEDfamLire(fid,nommaa,fam,curNbOfPolyElem,MED_MAILLE,MED_POLYEDRE);
      int arraySize=connFaceLgth;
      for(int i=0;i<curNbOfPolyElem;i++)
        arraySize+=index[i+1]-index[i]-1;
      int *finalConn=new int[arraySize];
      int *finalIndex=new int[curNbOfPolyElem+1];
      finalIndex[0]=1;
      int *wFinalConn=finalConn;
      for(int i=0;i<curNbOfPolyElem;i++)
        {
          finalIndex[i+1]=finalIndex[i]+index[i+1]-index[i]-1+indexFace[index[i+1]-1]-indexFace[index[i]-1];
          wFinalConn=std::copy(locConn+indexFace[index[i]-1]-1,locConn+indexFace[index[i]]-1,wFinalConn);
          for(int j=index[i];j<index[i+1]-1;j++)
            {
              *wFinalConn++=0;
              wFinalConn=std::copy(locConn+indexFace[j]-1,locConn+indexFace[j+1]-1,wFinalConn);
            }
        }
      delete [] index;
      delete [] locConn;
      delete [] indexFace;
      MEDLoader::MEDConnOfOneElemType elem(INTERP_KERNEL::NORM_POLYHED,finalConn,finalIndex,fam,curNbOfPolyElem,arraySize);
      conn.push_back(elem);
    }
}

namespace MEDLoaderNS
{
  template<class T>
  unsigned calculateHighestMeshDim(const std::list<T>& conn)
  {
    unsigned ret=0;
    for(typename std::list<T>::const_iterator iter=conn.begin();iter!=conn.end();iter++)
      {
        unsigned curDim=INTERP_KERNEL::CellModel::getCellModel((*iter).getType()).getDimension();
        if(ret<curDim)
          ret=curDim;
      }
    return ret;
  }
  
  template<class T>
  void keepSpecifiedMeshDim(typename std::list<T>& conn, unsigned meshDim)
  {
    for(typename std::list<T>::iterator iter=conn.begin();iter!=conn.end();)
      {
        unsigned curDim=INTERP_KERNEL::CellModel::getCellModel((*iter).getType()).getDimension();
        if(curDim!=meshDim)
          {
            (*iter).releaseArray();
            iter=conn.erase(iter);
          }
        else
          iter++;
      }
  }
  
  template<class T>
  void keepTypes(typename std::list<T>& conn, const std::vector<INTERP_KERNEL::NormalizedCellType>& typesToKeep)
  {
    if(!typesToKeep.empty())
      {
        for(typename std::list<T>::iterator iter=conn.begin();iter!=conn.end();)
          {
            INTERP_KERNEL::NormalizedCellType curType=(*iter).getType();
            if(std::find(typesToKeep.begin(),typesToKeep.end(),curType)==typesToKeep.end())
              {
                (*iter).releaseArray();
                iter=conn.erase(iter);
              }
            else
              iter++;
          }
      }
  }
}

class FieldPerTypeAccumulator
{
public:
  int operator()(int res, const MEDLoader::MEDFieldDoublePerCellType& elt) { return res+elt.getNbOfTuple(); }
};

ParaMEDMEM::DataArrayDouble *MEDLoaderNS::buildArrayFromRawData(const std::list<MEDLoader::MEDFieldDoublePerCellType>& fieldPerType)
{
  ParaMEDMEM::DataArrayDouble *ret=ParaMEDMEM::DataArrayDouble::New();
  int totalNbOfTuple=std::accumulate(fieldPerType.begin(),fieldPerType.end(),0,FieldPerTypeAccumulator());
  int nbOfComp=(*fieldPerType.begin()).getNbComp();
  double *ptr=new double[nbOfComp*totalNbOfTuple];
  ret->useArray(ptr,true,ParaMEDMEM::CPP_DEALLOC,totalNbOfTuple,nbOfComp);
  std::for_each(fieldPerType.begin(),fieldPerType.end(),FieldPerTypeCopier(ptr));
  return ret;
}

class PolyCounterForFams
{
public:
  PolyCounterForFams(int id, const int *index):_id(id),_index(index),_count(0),_sigma(0) { }
  void operator()(int val) { if(val==_id) _sigma+=_index[_count+1]-_index[_count]; _count++; }
  int getSigma() const { return _sigma; }
private:
  int _id;
  const int *_index;
  int _count;
  int _sigma;
};

void MEDLoaderNS::tradMEDFileCoreFrmt2MEDCouplingUMesh(const std::list<MEDLoader::MEDConnOfOneElemType>& medConnFrmt,
                                                       DataArrayInt* &conn,
                                                       DataArrayInt* &connIndex,
                                                       const std::vector<int>& familiesToKeep)
{
  bool keepAll=familiesToKeep.empty();
  if(medConnFrmt.empty())
    {
      conn=0;
      connIndex=0;
      return ;
    }
  std::list<MEDLoader::MEDConnOfOneElemType>::const_iterator iter=medConnFrmt.begin();
  int totalNbOfCells=0;
  int totalNbOfMedConn=0;
  for(;iter!=medConnFrmt.end();iter++)
    {
      const INTERP_KERNEL::CellModel& cellMod=INTERP_KERNEL::CellModel::getCellModel((*iter).getType());
      if(keepAll)
        totalNbOfCells+=(*iter).getLength();
      else
        for(std::vector<int>::const_iterator iter2=familiesToKeep.begin();iter2!=familiesToKeep.end();iter2++)
          totalNbOfCells+=std::count((*iter).getFam(),(*iter).getFam()+(*iter).getLength(),*iter2);
      if(!cellMod.isDynamic())
        if(keepAll)
          totalNbOfMedConn+=(*iter).getLength()*cellMod.getNumberOfNodes();
        else
          for(std::vector<int>::const_iterator iter2=familiesToKeep.begin();iter2!=familiesToKeep.end();iter2++)
            totalNbOfMedConn+=std::count((*iter).getFam(),(*iter).getFam()+(*iter).getLength(),*iter2)*cellMod.getNumberOfNodes();
      else
        if(keepAll)
          totalNbOfMedConn+=(*iter).getConnLength();
        else
          for(std::vector<int>::const_iterator iter2=familiesToKeep.begin();iter2!=familiesToKeep.end();iter2++)
            {
              PolyCounterForFams res=std::for_each((*iter).getFam(),(*iter).getFam()+(*iter).getLength(),PolyCounterForFams(*iter2,(*iter).getIndex()));
              totalNbOfMedConn+=res.getSigma();
            }
    }
  connIndex=DataArrayInt::New();
  conn=DataArrayInt::New();
  connIndex->alloc(totalNbOfCells+1,1);
  int *connIdxPtr=connIndex->getPointer();
  int connFillId=0;
  conn->alloc(totalNbOfMedConn+totalNbOfCells,1);
  int *connPtr=conn->getPointer();
  for(iter=medConnFrmt.begin();iter!=medConnFrmt.end();iter++)
    {
      INTERP_KERNEL::NormalizedCellType type=(*iter).getType();
      const int *sourceConn=(*iter).getArray();
      const int *sourceIndex=(*iter).getIndex();
      const INTERP_KERNEL::CellModel& cellMod=INTERP_KERNEL::CellModel::getCellModel(type);
      int nbOfCellsInCurType;
      int nbOfNodesIn1Cell=cellMod.getNumberOfNodes();
      nbOfCellsInCurType=(*iter).getLength();
      bool isDyn=cellMod.isDynamic();
      int *tmpConnPtr;
      for(int i=0;i<nbOfCellsInCurType;i++)
        {
          if(keepAll)
            {
              *connIdxPtr=connFillId;
              *connPtr++=type;
              if(!isDyn)
                tmpConnPtr=std::transform(sourceConn,sourceConn+nbOfNodesIn1Cell,connPtr,std::bind2nd(std::minus<int>(),1));
              else
                tmpConnPtr=std::transform(sourceConn,sourceConn+sourceIndex[i+1]-sourceIndex[i],connPtr,std::bind2nd(std::minus<int>(),1));
              connIdxPtr++;
              nbOfNodesIn1Cell=tmpConnPtr-connPtr;
              connFillId+=nbOfNodesIn1Cell+1;
              connPtr=tmpConnPtr;
            }
          else if(std::find(familiesToKeep.begin(),familiesToKeep.end(),(*iter).getFam()[i])!=familiesToKeep.end())
            {
              *connIdxPtr=connFillId;
              *connPtr++=type;
              if(!isDyn)
                tmpConnPtr=std::transform(sourceConn,sourceConn+nbOfNodesIn1Cell,connPtr,std::bind2nd(std::minus<int>(),1));
              else
                tmpConnPtr=std::transform(sourceConn,sourceConn+sourceIndex[i+1]-sourceIndex[i],connPtr,std::bind2nd(std::minus<int>(),1));
              connIdxPtr++;
              nbOfNodesIn1Cell=tmpConnPtr-connPtr;
              connFillId+=nbOfNodesIn1Cell+1;
              connPtr=tmpConnPtr;
            }
          sourceConn+=nbOfNodesIn1Cell;
        }
      *connIdxPtr=connFillId;
    }
}

namespace MEDLoaderNS
{
  template<class T>
  void releaseMEDFileCoreFrmt(typename std::list<T>& medConnFrmt)
  {
    for(typename std::list<T>::iterator iter=medConnFrmt.begin();iter!=medConnFrmt.end();iter++)
      (*iter).releaseArray();
    medConnFrmt.clear();
  }
}

/*!
 * This method builds a sub set of connectivity for a given type 'type'.
 * @param conn input containing connectivity with MEDCoupling format.
 * @param connIndex input containing connectivity index in MEDCoupling format.
 * @param families input that may be equal to 0. This specifies an array specifying cell family foreach cell.
 * @param type input specifying which cell types will be extracted in conn4MEDFile. 
 * @param conn4MEDFile output containing the connectivity directly understandable by MEDFile; conn4MEDFile has to be empty before this method called.
 * @param connIndex4MEDFile output containing index connectivity understandable by MEDFile; only used by polygons and polyhedrons (it is face nodal connec).
 * @param connIndexRk24MEDFile output containing index of rank 2 understandable by MEDFile; only used by polyhedrons.
 * @param fam4MEDFile output containing family number of cells whose type is 'type'. This output is updated only if 'families' is different than 0.
 * @return nb of elements extracted.
 */
int MEDLoaderNS::buildMEDSubConnectivityOfOneTypeStaticTypes(const DataArrayInt *conn, const DataArrayInt *connIndex, const DataArrayInt *families, INTERP_KERNEL::NormalizedCellType type, std::vector<int>& conn4MEDFile,
                                                             std::vector<int>& fam4MEDFile)
{
  int ret=0;
  int nbOfElem=connIndex->getNbOfElems()-1;
  const int *connPtr=conn->getConstPointer();
  const int *connIdxPtr=connIndex->getConstPointer();
  const int *famPtr=0;
  if(families)
    famPtr=families->getConstPointer();
  for(int i=0;i<nbOfElem;i++)
    {
      int delta=connIdxPtr[1]-connIdxPtr[0];
      if(*connPtr==type)
        {
          conn4MEDFile.insert(conn4MEDFile.end(),connPtr+1,connPtr+delta);
          if(families)
            fam4MEDFile.push_back(famPtr[i]);
          ret++;
        }
      connIdxPtr++;
      connPtr+=delta;
    }
  std::transform(conn4MEDFile.begin(),conn4MEDFile.end(),conn4MEDFile.begin(),std::bind2nd(std::plus<int>(),1));
  return ret;
}

int MEDLoaderNS::buildMEDSubConnectivityOfOneTypesPolyg(const DataArrayInt *conn, const DataArrayInt *connIndex, const DataArrayInt *families, std::vector<int>& conn4MEDFile, std::vector<int>& connIndex4MEDFile,
                                                        std::vector<int>& fam4MEDFile)
{
  int ret=0;
  int nbOfElem=connIndex->getNbOfElems()-1;
  const int *connPtr=conn->getConstPointer();
  const int *connIdxPtr=connIndex->getConstPointer();
  connIndex4MEDFile.push_back(1);
  const int *famPtr=0;
  if(families)
    famPtr=families->getConstPointer();
  for(int i=0;i<nbOfElem;i++)
    {
      int delta=connIdxPtr[1]-connIdxPtr[0];
      if(*connPtr==INTERP_KERNEL::NORM_POLYGON)
        {
          conn4MEDFile.insert(conn4MEDFile.end(),connPtr+1,connPtr+delta);
          connIndex4MEDFile.push_back(connIndex4MEDFile.back()+delta-1);
          if(families)
            fam4MEDFile.push_back(famPtr[i]);
          ret++;
        }
      connIdxPtr++;
      connPtr+=delta;
    }
  std::transform(conn4MEDFile.begin(),conn4MEDFile.end(),conn4MEDFile.begin(),std::bind2nd(std::plus<int>(),1));
  return ret;
}
  
int MEDLoaderNS::buildMEDSubConnectivityOfOneTypesPolyh(const DataArrayInt *conn, const DataArrayInt *connIndex, const DataArrayInt *families,
                                                        std::vector<int>& conn4MEDFile, std::vector<int>& connIndex4MEDFile, std::vector<int>& connIndexRk24MEDFile,
                                                        std::vector<int>& fam4MEDFile)
{
  int ret=0;
  int nbOfElem=connIndex->getNbOfElems()-1;
  const int *connPtr=conn->getConstPointer();
  const int *connIdxPtr=connIndex->getConstPointer();
  connIndexRk24MEDFile.push_back(1);
  connIndex4MEDFile.push_back(1);
  const int *famPtr=0;
  if(families)
    famPtr=families->getConstPointer();
  for(int i=0;i<nbOfElem;i++)
    {
      int delta=connIdxPtr[1]-connIdxPtr[0];
      if(*connPtr==INTERP_KERNEL::NORM_POLYHED)
        {
          int nbOfFacesOfPolyh=std::count(connPtr+1,connPtr+delta,-1)+1;
          const int *work=connPtr+1;
          while(work!=connPtr+delta)
            {
              const int *end=std::find(work,connPtr+delta,-1);
              conn4MEDFile.insert(conn4MEDFile.end(),work,end);
              connIndex4MEDFile.push_back(connIndex4MEDFile.back()+std::distance(work,end));
              if(end==connPtr+delta)
                work=connPtr+delta;
              else
                work=end+1;
            }
          connIndexRk24MEDFile.push_back(connIndexRk24MEDFile.back()+nbOfFacesOfPolyh);
          if(families)
            fam4MEDFile.push_back(famPtr[i]);
          ret++;
        }
      connIdxPtr++;
      connPtr+=delta;
    }
  std::transform(conn4MEDFile.begin(),conn4MEDFile.end(),conn4MEDFile.begin(),std::bind2nd(std::plus<int>(),1));
  return ret;
}
  
/*!
 * This method builds a sub set of connectivity for a given type 'type'.
 * @param conn input containing connectivity with MEDCoupling format.
 * @param connIndex input containing connectivity index in MEDCoupling format.
 * @param families input containing, if any, the family number of each cells
 * @param type input specifying which cell types will be extracted in conn4MEDFile. 
 * @param conn4MEDFile output containing the connectivity directly understandable by MEDFile; conn4MEDFile has to be empty before this method called.
 * @param connIndex4MEDFile output containing index connectivity understandable by MEDFile; only used by polygons and polyhedrons (it is face nodal connec).
 * @param connIndexRk24MEDFile output containing index of rank 2 understandable by MEDFile; only used by polyhedrons.
 * @param fam4MEDFile output containing families id of cells whose type is 'type'.
 * @return nb of elements extracted.
 */
int MEDLoaderNS::buildMEDSubConnectivityOfOneType(const DataArrayInt *conn, const DataArrayInt *connIndex, const DataArrayInt *families,
                                                  INTERP_KERNEL::NormalizedCellType type, std::vector<int>& conn4MEDFile,
                                                  std::vector<int>& connIndex4MEDFile, std::vector<int>& connIndexRk24MEDFile, std::vector<int>& fam4MEDFile)
{
    
  const INTERP_KERNEL::CellModel& cellMod=INTERP_KERNEL::CellModel::getCellModel(type);
  if(!cellMod.isDynamic())
    return buildMEDSubConnectivityOfOneTypeStaticTypes(conn,connIndex,families,type,conn4MEDFile,fam4MEDFile);
  else
    {
      if(type==INTERP_KERNEL::NORM_POLYGON)
        return buildMEDSubConnectivityOfOneTypesPolyg(conn,connIndex,families,conn4MEDFile,connIndex4MEDFile,fam4MEDFile);
      else
        return buildMEDSubConnectivityOfOneTypesPolyh(conn,connIndex,families,conn4MEDFile,connIndex4MEDFile,connIndexRk24MEDFile,fam4MEDFile);
    }
}
  
/*!
 * @param ids is a in vector containing families ids whose cells have to be kept. If empty all cells are kept.
 * @param typesToKeep is a in vector that indicates which types to keep after dimension filtering.
 * @param meshDimExtract out parameter that gives the mesh dimension.
 */
MEDCouplingUMesh *MEDLoaderNS::readUMeshFromFileLev1(const char *fileName, const char *meshName, int meshDimRelToMax, const std::vector<int>& ids,
                                                     const std::vector<INTERP_KERNEL::NormalizedCellType>& typesToKeep, unsigned& meshDimExtract) throw(INTERP_KERNEL::Exception)
{
  //Extraction data from MED file.
  med_idt fid=MEDouvrir((char *)fileName,MED_LECTURE);
  std::string trueMeshName;
  med_int mid=getIdFromMeshName(fid,meshName,trueMeshName);
  double *coords;
  int nCoords;
  int spaceDim;
  std::list<MEDLoader::MEDConnOfOneElemType> conn;
  readUMeshDataInMedFile(fid,mid,coords,nCoords,spaceDim,conn);
  meshDimExtract=MEDLoaderNS::calculateHighestMeshDim<MEDLoader::MEDConnOfOneElemType>(conn);
  meshDimExtract=meshDimExtract+meshDimRelToMax;
  MEDLoaderNS::keepSpecifiedMeshDim<MEDLoader::MEDConnOfOneElemType>(conn,meshDimExtract);
  MEDLoaderNS::keepTypes<MEDLoader::MEDConnOfOneElemType>(conn,typesToKeep);
  MEDfermer(fid);
  //Put data in returned data structure.
  MEDCouplingUMesh *ret=MEDCouplingUMesh::New();
  ret->setName(trueMeshName.c_str());
  ret->setMeshDimension(meshDimExtract);
  //
  DataArrayDouble *coordsArr=DataArrayDouble::New();
  coordsArr->useArray(coords,true,ParaMEDMEM::CPP_DEALLOC,nCoords,spaceDim);
  ret->setCoords(coordsArr);
  coordsArr->decrRef();
  //
  DataArrayInt *connArr,*connIndexArr;
  tradMEDFileCoreFrmt2MEDCouplingUMesh(conn,connArr,connIndexArr,ids);
  ret->setConnectivity(connArr,connIndexArr);
  //clean-up
  if(connArr)
    connArr->decrRef();
  if(connIndexArr)
    connIndexArr->decrRef();
  releaseMEDFileCoreFrmt<MEDLoader::MEDConnOfOneElemType>(conn);
  return ret;
}

ParaMEDMEM::MEDCouplingFieldDouble *MEDLoaderNS::readFieldDoubleLev1(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order,
                                                                     ParaMEDMEM::TypeOfField typeOfOutField)
{
  std::list<MEDLoader::MEDFieldDoublePerCellType> fieldPerCellType;
  double time;
  readFieldDoubleDataInMedFile(fileName,meshName,fieldName,fieldPerCellType,iteration,order,typeOfOutField,time);
  std::vector<int> familiesToKeep;
  std::vector<INTERP_KERNEL::NormalizedCellType> typesToKeep;
  if(typeOfOutField==ON_CELLS)
    for(std::list<MEDLoader::MEDFieldDoublePerCellType>::const_iterator iter=fieldPerCellType.begin();iter!=fieldPerCellType.end();iter++)
      typesToKeep.push_back((*iter).getType());
  unsigned meshDim;
  ParaMEDMEM::MEDCouplingUMesh *mesh=readUMeshFromFileLev1(fileName,meshName,meshDimRelToMax,familiesToKeep,typesToKeep,meshDim);
  if(typeOfOutField==ON_CELLS)
    MEDLoaderNS::keepSpecifiedMeshDim<MEDLoader::MEDFieldDoublePerCellType>(fieldPerCellType,meshDim);
  ParaMEDMEM::MEDCouplingFieldDouble *ret=ParaMEDMEM::MEDCouplingFieldDouble::New(typeOfOutField,ONE_TIME);
  ret->setName(fieldName);
  ret->setTime(time,iteration,order);
  ret->setMesh(mesh);
  mesh->decrRef();
  ParaMEDMEM::DataArrayDouble *arr=buildArrayFromRawData(fieldPerCellType);
  ret->setArray(arr);
  arr->decrRef();
  releaseMEDFileCoreFrmt<MEDLoader::MEDFieldDoublePerCellType>(fieldPerCellType);
  return ret;
}

MEDCouplingUMesh *MEDLoader::ReadUMeshFromFile(const char *fileName, const char *meshName, int meshDimRelToMax) throw(INTERP_KERNEL::Exception)
{
  std::vector<int> familiesToKeep;
  std::vector<INTERP_KERNEL::NormalizedCellType> typesToKeep;
  unsigned meshDim;
  return MEDLoaderNS::readUMeshFromFileLev1(fileName,meshName,meshDimRelToMax,familiesToKeep,typesToKeep,meshDim);
}

ParaMEDMEM::MEDCouplingUMesh *MEDLoader::ReadUMeshFromFile(const char *fileName, int meshDimRelToMax) throw(INTERP_KERNEL::Exception)
{
  std::vector<int> familiesToKeep;
  std::vector<INTERP_KERNEL::NormalizedCellType> typesToKeep;
  unsigned meshDim;
  return MEDLoaderNS::readUMeshFromFileLev1(fileName,0,meshDimRelToMax,familiesToKeep,typesToKeep,meshDim);
}

ParaMEDMEM::MEDCouplingUMesh *MEDLoader::ReadUMeshFromFamilies(const char *fileName, const char *meshName, int meshDimRelToMax, const std::vector<std::string>& fams)
{
  std::vector<int> familiesToKeep=MEDLoaderNS::getIdsFromFamilies(fileName,meshName,fams);
  std::vector<INTERP_KERNEL::NormalizedCellType> typesToKeep;
  unsigned meshDim;
  return MEDLoaderNS::readUMeshFromFileLev1(fileName,meshName,meshDimRelToMax,familiesToKeep,typesToKeep,meshDim);
}

ParaMEDMEM::MEDCouplingUMesh *MEDLoader::ReadUMeshFromGroups(const char *fileName, const char *meshName, int meshDimRelToMax, const std::vector<std::string>& grps)
{
  std::vector<int> familiesToKeep=MEDLoaderNS::getIdsFromGroups(fileName,meshName,grps);
  std::vector<INTERP_KERNEL::NormalizedCellType> typesToKeep;
  unsigned meshDim;
  return MEDLoaderNS::readUMeshFromFileLev1(fileName,meshName,meshDimRelToMax,familiesToKeep,typesToKeep,meshDim);
}

ParaMEDMEM::MEDCouplingFieldDouble *MEDLoader::ReadFieldDouble(ParaMEDMEM::TypeOfField type, const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order)
{
  switch(type)
    {
    case ON_CELLS:
      return ReadFieldDoubleCell(fileName,meshName,meshDimRelToMax,fieldName,iteration,order);
    case ON_NODES:
      return ReadFieldDoubleNode(fileName,meshName,meshDimRelToMax,fieldName,iteration,order);
    default:
      throw INTERP_KERNEL::Exception("Type of field specified not managed ! manages are ON_NODES or ON_CELLS !");
    } 
}

ParaMEDMEM::MEDCouplingFieldDouble *MEDLoader::ReadFieldDoubleCell(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order)
{
  return MEDLoaderNS::readFieldDoubleLev1(fileName,meshName,meshDimRelToMax,fieldName,iteration,order,ON_CELLS);
}

ParaMEDMEM::MEDCouplingFieldDouble *MEDLoader::ReadFieldDoubleNode(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order)
{
  return MEDLoaderNS::readFieldDoubleLev1(fileName,meshName,meshDimRelToMax,fieldName,iteration,order,ON_NODES);
}

void MEDLoaderNS::writeUMeshDirectly(const char *fileName, ParaMEDMEM::MEDCouplingUMesh *mesh, const DataArrayInt *families, bool forceFromScratch)
{
  med_idt fid=MEDouvrir((char *)fileName,forceFromScratch?MED_CREATION:MED_LECTURE_ECRITURE);
  std::string meshName(mesh->getName());
  if(meshName=="")
    {
      MEDfermer(fid);
      throw INTERP_KERNEL::Exception("MEDCouplingMesh must have a not null name !");
    }
  char *maa=MEDLoaderBase::buildEmptyString(MED_TAILLE_NOM);
  char *desc=MEDLoaderBase::buildEmptyString(MED_TAILLE_DESC);
  strcpy(maa,meshName.c_str());
  strcpy(desc,meshName.c_str());
  MEDmaaCr(fid,maa,mesh->getSpaceDimension(),MED_NON_STRUCTURE,desc);
  MEDdimEspaceCr(fid,maa,mesh->getSpaceDimension());
  std::set<INTERP_KERNEL::NormalizedCellType> allTypes(mesh->getAllTypes());
  DataArrayInt *conn=mesh->getNodalConnectivity();
  DataArrayInt *connIndex=mesh->getNodalConnectivityIndex();
  char familyName[MED_TAILLE_NOM+1];
  std::fill(familyName,familyName+MED_TAILLE_NOM+1,'\0');
  const char DftFamilyName[]="DftFamily";
  std::copy(DftFamilyName,DftFamilyName+sizeof(DftFamilyName),familyName);
  for(int i=0;i<MED_NBR_GEOMETRIE_MAILLE+2;i++)
    {
      med_geometrie_element curMedType=typmai[i];
      INTERP_KERNEL::NormalizedCellType curType=typmai2[i];
      if(allTypes.find(curType)!=allTypes.end())
        {
          std::vector<int> medConn;
          std::vector<int> medConnIndex;
          std::vector<int> medConnIndex2;
          std::vector<int> fam;
          int nbOfElt=MEDLoaderNS::buildMEDSubConnectivityOfOneType(conn,connIndex,families,curType,medConn,medConnIndex,medConnIndex2,fam);
          if(curMedType!=MED_POLYGONE && curMedType!=MED_POLYEDRE)
            MEDconnEcr(fid,maa,mesh->getMeshDimension(),&medConn[0],MED_FULL_INTERLACE,nbOfElt,MED_MAILLE,curMedType,MED_NOD);
          else
            {
              if(curMedType==MED_POLYGONE)
                MEDpolygoneConnEcr(fid,maa,&medConnIndex[0],medConnIndex.size(),&medConn[0],MED_MAILLE,MED_NOD);
              if(curMedType==MED_POLYEDRE)
                {
                  MEDpolyedreConnEcr(fid,maa,&medConnIndex2[0],medConnIndex2.size(),&medConnIndex[0],medConnIndex.size(),
                                     &medConn[0],MED_NOD);
                }
            }
          if(families)
            MEDfamEcr(fid,maa,&fam[0],nbOfElt,MED_MAILLE,curMedType);
        }
    }
  MEDfamCr(fid,maa,familyName,0,0,0,0,0,0,0);
  DataArrayDouble *arr=mesh->getCoords();
  char comp[2*MED_TAILLE_PNOM+1];
  char unit[2*MED_TAILLE_PNOM+1];
  std::fill(comp,comp+2*MED_TAILLE_PNOM,' ');
  comp[2*MED_TAILLE_PNOM]='\0';
  char *work=comp;
  for(int i=0;i<mesh->getSpaceDimension();i++,work+=3)
    *work='X'+i;
  std::fill(unit,unit+2*MED_TAILLE_PNOM+1,'\0');
  MEDcoordEcr(fid,maa,mesh->getSpaceDimension(),arr->getPointer(),MED_FULL_INTERLACE,mesh->getNumberOfNodes(),MED_CART,comp,unit);
  delete [] maa;
  delete [] desc;
  MEDfermer(fid);
}

/*!
 * In this method meshes are assumed to shared the same coords.
 * This method makes the assumption that 'meshes' is not empty, no check on that is done (responsability of the caller)
 */
void MEDLoaderNS::writeUMeshesDirectly(const char *fileName, const char *meshName, const std::vector<ParaMEDMEM::MEDCouplingUMesh *>& meshes, bool forceFromScratch)
{
  std::string meshNameCpp(meshName);
  char *maa=MEDLoaderBase::buildEmptyString(MED_TAILLE_NOM);
  strcpy(maa,meshName);
  if(meshName=="")
    throw INTERP_KERNEL::Exception("writeUMeshesDirectly : Invalid meshName : Must be different from \"\" !");
  //MEDnumEcr(fid,maa,num,nele,_type_ent,typ_geo);
  std::vector< DataArrayInt * > corr;
  MEDCouplingUMesh *m=ParaMEDMEM::MEDCouplingUMesh::fuseUMeshesOnSameCoords(meshes,0,corr);
  m->setName(meshName);
  std::vector< std::vector<int> > fidsOfGroups;
  DataArrayInt *arr2=DataArrayInt::makePartition(corr,m->getNumberOfCells(),fidsOfGroups);
  for(std::vector< DataArrayInt * >::iterator it=corr.begin();it!=corr.end();it++)
    (*it)->decrRef();
  writeUMeshDirectly(fileName,m,arr2,forceFromScratch);
  // families creation
  std::set<int> familyIds;
  for(std::vector< std::vector<int> >::const_iterator it1=fidsOfGroups.begin();it1!=fidsOfGroups.end();it1++)
    for(std::vector<int>::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
      familyIds.insert(*it2);
  std::vector< std::vector<int> > gidsOfFamilies(familyIds.size());
  int fid=0;
  for(std::set<int>::const_iterator it=familyIds.begin();it!=familyIds.end();it++,fid++)
    {
      int gid=0;
      for(std::vector< std::vector<int> >::const_iterator it1=fidsOfGroups.begin();it1!=fidsOfGroups.end();it1++,gid++)
        for(std::vector<int>::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
          if(*it2==*it)
            gidsOfFamilies[fid].push_back(gid);
    }
  fid=0;
  med_idt fid2=MEDouvrir((char *)fileName,MED_LECTURE_ECRITURE);
  for(std::set<int>::const_iterator it=familyIds.begin();it!=familyIds.end();it++,fid++)
    {
      int ngro=gidsOfFamilies[fid].size();
      char *groName=MEDLoaderBase::buildEmptyString(MED_TAILLE_LNOM*ngro);
      for(int i=0;i<ngro;i++)
        strcpy(groName+i*MED_TAILLE_LNOM,meshes[gidsOfFamilies[fid][i]]->getName());
      std::ostringstream oss; oss << "Family_" << *it;
      char *famName=MEDLoaderBase::buildEmptyString(MED_TAILLE_NOM);
      strcpy(famName,oss.str().c_str());
      MEDfamCr(fid2,maa,famName,*it,0,0,0,0,groName,ngro);
      delete [] famName;
      delete [] groName;
    }
  MEDfermer(fid2);
  // end families creation
  delete [] maa;
  arr2->decrRef();
  m->decrRef();
}

void MEDLoaderNS::appendFieldDirectly(const char *fileName, ParaMEDMEM::MEDCouplingFieldDouble *f)
{
  med_idt fid=MEDouvrir((char *)fileName,MED_LECTURE_ECRITURE);
  int nbComp=f->getNumberOfComponents();
  char *comp=MEDLoaderBase::buildEmptyString(nbComp*MED_TAILLE_PNOM);
  char *unit=MEDLoaderBase::buildEmptyString(nbComp*MED_TAILLE_PNOM);
  MEDchampCr(fid,(char *)f->getName(),MED_FLOAT64,comp,unit,nbComp);
  med_int numdt,numo;
  med_float dt;
  ParaMEDMEM::TypeOfTimeDiscretization td=f->getTimeDiscretization();
  if(td==ParaMEDMEM::NO_TIME)
    {
      numdt=MED_NOPDT; numo=MED_NONOR; dt=0.0;
    }
  else if(td==ParaMEDMEM::ONE_TIME)
    {
      int tmp1,tmp2;
      double tmp0=f->getTime(tmp1,tmp2);
      numdt=(med_int)tmp1; numo=(med_int)tmp2;
      dt=(med_float)tmp0;
    }
  const double *pt=f->getArray()->getConstPointer();
  switch(f->getTypeOfField())
    {
    case ParaMEDMEM::ON_CELLS:
      {
        std::list<MEDLoader::MEDFieldDoublePerCellType> split;
        prepareCellFieldDoubleForWriting(f,split);
        for(std::list<MEDLoader::MEDFieldDoublePerCellType>::const_iterator iter=split.begin();iter!=split.end();iter++)
          {
            char *nommaa=MEDLoaderBase::buildEmptyString(MED_TAILLE_NOM);
            strcpy(nommaa,f->getMesh()->getName());
            MEDchampEcr(fid,nommaa,(char *)f->getName(),(unsigned char*)pt,MED_FULL_INTERLACE,(*iter).getNbOfTuple(),
                        (char *)MED_NOGAUSS,MED_ALL,(char *)MED_NOPFL,MED_NO_PFLMOD,MED_MAILLE,
                        typmai3[(int)(*iter).getType()],numdt,(char *)"",dt,numo);
            delete [] nommaa;
            pt+=(*iter).getNbOfTuple()*nbComp;
          }
        break;
      }
    case ParaMEDMEM::ON_NODES:
      {
        int nbOfTuples=f->getArray()->getNumberOfTuples();
        char *nommaa=MEDLoaderBase::buildEmptyString(MED_TAILLE_NOM);
        strcpy(nommaa,f->getMesh()->getName());
        MEDchampEcr(fid,nommaa,(char *)f->getName(),(unsigned char*)pt,MED_FULL_INTERLACE,nbOfTuples,(char *)MED_NOGAUSS,
                    MED_ALL,(char *)MED_NOPFL,MED_NO_PFLMOD,MED_NOEUD,MED_NONE,numdt,(char *)"",dt,numo);
        delete [] nommaa;
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("Not managed this type of FIELD !");
    }
  delete [] comp;
  delete [] unit;
  MEDfermer(fid);
}

void MEDLoaderNS::prepareCellFieldDoubleForWriting(const ParaMEDMEM::MEDCouplingFieldDouble *f, std::list<MEDLoader::MEDFieldDoublePerCellType>& split)
{
  int nbComp=f->getNumberOfComponents();
  const MEDCouplingMesh *mesh=f->getMesh();
  const MEDCouplingUMesh *meshC=dynamic_cast<const MEDCouplingUMesh *>(mesh);
  if(!meshC)
    throw INTERP_KERNEL::Exception("Not implemented yet for not unstructured mesh !");
  if(!meshC->checkConsecutiveCellTypes())
    throw INTERP_KERNEL::Exception("Unstructuded mesh has not consecutive cell types !");
  const int *connI=meshC->getNodalConnectivityIndex()->getConstPointer();
  const int *conn=meshC->getNodalConnectivity()->getConstPointer();
  int nbOfCells=meshC->getNumberOfCells();
  INTERP_KERNEL::NormalizedCellType curType;
  for(const int *pt=connI;pt!=connI+nbOfCells;)
    {
      curType=(INTERP_KERNEL::NormalizedCellType)conn[*pt];
      const int *pt2=std::find_if(pt+1,connI+nbOfCells,ConnReaderML(conn,(int)curType));
      split.push_back(MEDLoader::MEDFieldDoublePerCellType(curType,0,nbComp,pt2-pt));
      pt=pt2;
    }
}

void MEDLoaderNS::writeFieldAndMeshDirectly(const char *fileName, ParaMEDMEM::MEDCouplingFieldDouble *f, bool forceFromScratch)
{
  std::string meshName(f->getMesh()->getName());
  if(meshName.empty())
    throw INTERP_KERNEL::Exception("Trying to write a mesh (f->getMesh()) with no name ! MED file format needs a not empty mesh name !");
  std::string fieldName(f->getName());
  if(fieldName.empty())
    throw INTERP_KERNEL::Exception("Trying to write a field with no name ! MED file format needs a not empty field name !");
  MEDCouplingUMesh *mesh=dynamic_cast<MEDCouplingUMesh *>((MEDCouplingMesh *)f->getMesh());
  writeUMeshDirectly(fileName,mesh,0,forceFromScratch);
  appendFieldDirectly(fileName,f);
}

void MEDLoader::WriteUMesh(const char *fileName, ParaMEDMEM::MEDCouplingUMesh *mesh, bool writeFromScratch)
{
  std::string meshName(mesh->getName());
  if(meshName.empty())
    throw INTERP_KERNEL::Exception("Trying to write a unstructured mesh with no name ! MED file format needs a not empty mesh name !");
  int status=MEDLoaderBase::getStatusOfFile(fileName);
  if(status!=MEDLoaderBase::EXIST_RW && status!=MEDLoaderBase::NOT_EXIST)
    {
      std::ostringstream oss; oss << "File with name \'" << fileName << "\' has not valid permissions !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  if(writeFromScratch)
    {
      MEDLoaderNS::writeUMeshDirectly(fileName,mesh,0,true);
      return ;
    }
  if(status==MEDLoaderBase::NOT_EXIST)
    {
      MEDLoaderNS::writeUMeshDirectly(fileName,mesh,0,true);
      return;
    }
  else
    {
      std::vector<std::string> meshNames=GetMeshNames(fileName);
      if(std::find(meshNames.begin(),meshNames.end(),meshName)==meshNames.end())
        MEDLoaderNS::writeUMeshDirectly(fileName,mesh,0,false);
      else
        {
          std::ostringstream oss; oss << "File \'" << fileName << "\' already exists and has already a mesh called \"";
          oss << meshName << "\" !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
}

void MEDLoader::WriteUMeshes(const char *fileName, const char *meshNameC, const std::vector<ParaMEDMEM::MEDCouplingUMesh *>& meshes, bool writeFromScratch)
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
  if(meshes.empty())
    throw INTERP_KERNEL::Exception("List of meshes must be not empty !");
  DataArrayDouble *coords=meshes.front()->getCoords();
  for(std::vector<ParaMEDMEM::MEDCouplingUMesh *>::const_iterator iter=meshes.begin();iter!=meshes.end();iter++)
    if(coords!=(*iter)->getCoords())
      throw INTERP_KERNEL::Exception("Meshes does not not share the same coordinates : try method MEDCouplingPointSet::tryToShareSameCoords !");
  std::set<std::string> tmp;
  for(std::vector<ParaMEDMEM::MEDCouplingUMesh *>::const_iterator iter=meshes.begin();iter!=meshes.end();iter++)
    {
      if(tmp.find((*iter)->getName())==tmp.end())
        tmp.insert((*iter)->getName());
      else
        throw INTERP_KERNEL::Exception("The names of meshes must be different each other !");
    }
  tmp.clear();
  if(writeFromScratch)
    {
      MEDLoaderNS::writeUMeshesDirectly(fileName,meshNameC,meshes,true);
      return ;
    }
  if(status==MEDLoaderBase::NOT_EXIST)
    {
      MEDLoaderNS::writeUMeshesDirectly(fileName,meshNameC,meshes,true);
      return;
    }
  else
    {
      std::vector<std::string> meshNames=GetMeshNames(fileName);
      if(std::find(meshNames.begin(),meshNames.end(),meshName)==meshNames.end())
        MEDLoaderNS::writeUMeshesDirectly(fileName,meshNameC,meshes,false);
      else
        {
          std::ostringstream oss; oss << "File \'" << fileName << "\' already exists and has already a mesh called \"";
          oss << meshName << "\" !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
}

void MEDLoader::WriteField(const char *fileName, ParaMEDMEM::MEDCouplingFieldDouble *f, bool writeFromScratch)
{
  int status=MEDLoaderBase::getStatusOfFile(fileName);
  if(status!=MEDLoaderBase::EXIST_RW && status!=MEDLoaderBase::NOT_EXIST)
    {
      std::ostringstream oss; oss << "File with name \'" << fileName << "\' has not valid permissions !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  if(writeFromScratch)
    {
      MEDLoaderNS::writeFieldAndMeshDirectly(fileName,f,true);
     return ;
    }
  if(status==MEDLoaderBase::NOT_EXIST)
    {
     MEDLoaderNS::writeFieldAndMeshDirectly(fileName,f,true);
     return ;
    }
  else
    {
      std::vector<std::string> meshNames=GetMeshNames(fileName);
      std::string fileNameCpp(f->getMesh()->getName());
      if(std::find(meshNames.begin(),meshNames.end(),fileNameCpp)==meshNames.end())
        MEDLoaderNS::writeFieldAndMeshDirectly(fileName,f,false);
      else
        {
          std::ostringstream oss; oss << "File \'" << fileName << "\' already exists and has already a mesh called \"";
          oss << fileNameCpp << "\" !";
          throw INTERP_KERNEL::Exception(oss.str().c_str()); 
        }
    }
}

void MEDLoader::WriteFieldUsingAlreadyWrittenMesh(const char *fileName, ParaMEDMEM::MEDCouplingFieldDouble *f)
{
  f->checkCoherency();
  int status=MEDLoaderBase::getStatusOfFile(fileName);
  if(status!=MEDLoaderBase::EXIST_RW)
    {
      std::ostringstream oss; oss << "File with name \'" << fileName << "\' has not valid permissions or not exists !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  MEDLoaderNS::appendFieldDirectly(fileName,f);
}
