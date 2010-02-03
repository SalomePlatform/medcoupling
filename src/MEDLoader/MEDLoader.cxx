//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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
  
  std::string buildStringFromFortran(const char *expr, int lgth);
  std::vector<std::string> getMeshNamesFid(med_idt fid);
  void readFieldDoubleDataInMedFile(const char *fileName, const char *meshName, const char *fieldName, std::list<MEDLoader::MEDFieldDoublePerCellType>& field,
                                    int iteration, int order, ParaMEDMEM::TypeOfField typeOfOutField, double& time);
  std::vector<int> getIdsFromFamilies(const char *fileName, const char *meshName, const std::vector<std::string>& fams);
  std::vector<int> getIdsFromGroups(const char *fileName, const char *meshName, const std::vector<std::string>& grps);
  med_int getIdFromMeshName(med_idt fid, const char *meshName, std::string& trueMeshName) throw(INTERP_KERNEL::Exception);
  void dispatchElems(int nbOfElemCell, int nbOfElemFace, int& nbOfElem, med_entite_maillage& whichEntity);
  void readUMeshDataInMedFile(med_idt fid, med_int meshId, double *&coords, int& nCoords, int& spaceDim, std::list<MEDLoader::MEDConnOfOneElemType>& conn);
  int buildMEDSubConnectivityOfOneType(DataArrayInt *conn, DataArrayInt *connIndex, INTERP_KERNEL::NormalizedCellType type, std::vector<int>& conn4MEDFile,
                                       std::vector<int>& connIndex4MEDFile, std::vector<int>& connIndexRk24MEDFile);
  MEDCouplingUMesh *readUMeshFromFileLev1(const char *fileName, const char *meshName, int meshDimRelToMax, const std::vector<int>& ids,
                                          const std::vector<INTERP_KERNEL::NormalizedCellType>& typesToKeep, unsigned& meshDimExtract) throw(INTERP_KERNEL::Exception);
  void tradMEDFileCoreFrmt2MEDCouplingUMesh(const std::list<MEDLoader::MEDConnOfOneElemType>& medConnFrmt,
                                            DataArrayInt* &conn,
                                            DataArrayInt* &connIndex,
                                            const std::vector<int>& familiesToKeep);
  ParaMEDMEM::DataArrayDouble *buildArrayFromRawData(const std::list<MEDLoader::MEDFieldDoublePerCellType>& fieldPerType);
  int buildMEDSubConnectivityOfOneTypesPolyg(DataArrayInt *conn, DataArrayInt *connIndex, std::vector<int>& conn4MEDFile, std::vector<int>& connIndex4MEDFile);
  int buildMEDSubConnectivityOfOneTypesPolyh(DataArrayInt *conn, DataArrayInt *connIndex, std::vector<int>& conn4MEDFile,
                                             std::vector<int>& connIndex4MEDFile, std::vector<int>& connIndexRk24MEDFile);
  int buildMEDSubConnectivityOfOneTypeStaticTypes(DataArrayInt *conn, DataArrayInt *connIndex, INTERP_KERNEL::NormalizedCellType type, std::vector<int>& conn4MEDFile);
  ParaMEDMEM::MEDCouplingFieldDouble *readFieldDoubleLev1(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order,
                                                          ParaMEDMEM::TypeOfField typeOfOutField);
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

MEDLoader::MEDFieldDoublePerCellType::MEDFieldDoublePerCellType(INTERP_KERNEL::NormalizedCellType type, double *values, int ncomp, int nval):_nval(nval),_ncomp(ncomp),_values(values),_type(type)
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
  char maa_ass[MED_TAILLE_NOM+1]="";
  char dt_unit[MED_TAILLE_PNOM+1]="";
  char nomcha[MED_TAILLE_NOM+1]="";
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
  char maa_ass[MED_TAILLE_NOM+1]="";
  char dt_unit[MED_TAILLE_PNOM+1]="";
  char nomcha[MED_TAILLE_NOM+1]="";
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
  MEDfermer(fid);
  return ret;
}

std::vector< std::pair<int,int> > MEDLoader::GetCellFieldIterations(const char *fileName, const char *fieldName)
{
  std::vector< std::pair<int,int> > ret;
  med_idt fid=MEDouvrir((char *)fileName,MED_LECTURE);
  med_int nbFields=MEDnChamp(fid,0);
  //
  med_type_champ typcha;
  med_int ngauss=0;
  med_int numdt=0,numo=0,nbrefmaa;
  med_float dt=0.0;
  med_booleen local;
  char maa_ass[MED_TAILLE_NOM+1]="";
  char dt_unit[MED_TAILLE_PNOM+1]="";
  char nomcha[MED_TAILLE_NOM+1]="";
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
                  found=true;
                  ret.push_back(std::make_pair(numdt,numo));
                }
            }
        }
    }
  MEDfermer(fid);
  return ret;
}

std::vector< std::pair<int,int> > MEDLoader::GetNodeFieldIterations(const char *fileName, const char *fieldName)
{
  std::vector< std::pair<int,int> > ret;
  med_idt fid=MEDouvrir((char *)fileName,MED_LECTURE);
  med_int nbFields=MEDnChamp(fid,0);
  //
  med_type_champ typcha;
  med_int ngauss=0;
  med_int numdt=0,numo=0,nbrefmaa;
  med_float dt=0.0;
  med_booleen local;
  char maa_ass[MED_TAILLE_NOM+1]="";
  char dt_unit[MED_TAILLE_PNOM+1]="";
  char nomcha[MED_TAILLE_NOM+1]="";
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
              ret.push_back(std::make_pair(numdt,numo));
            }
        }
    }
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
                  MEDchampLire(fid,(char *)meshName,(char *)fieldName,(unsigned char*)valr,MED_FULL_INTERLACE,MED_ALL,locname,
                               pflname,MED_COMPACT,tabEnt[typeOfOutField],tabType[typeOfOutField][j],iteration,order);
                  field.push_back(MEDLoader::MEDFieldDoublePerCellType(typmai2[j],valr,ncomp,nval));
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
 * @param type input specifying which cell types will be extracted in conn4MEDFile. 
 * @param conn4MEDFile output containing the connectivity directly understandable by MEDFile; conn4MEDFile has to be empty before this method called.
 * @param connIndex4MEDFile output containing index connectivity understandable by MEDFile; only used by polygons and polyhedrons (it is face nodal connec).
 * @param connIndexRk24MEDFile output containing index of rank 2 understandable by MEDFile; only used by polyhedrons.
 * @return nb of elements extracted.
 */
int MEDLoaderNS::buildMEDSubConnectivityOfOneTypeStaticTypes(DataArrayInt *conn, DataArrayInt *connIndex, INTERP_KERNEL::NormalizedCellType type, std::vector<int>& conn4MEDFile)
{
  int ret=0;
  int nbOfElem=connIndex->getNbOfElems()-1;
  const int *connPtr=conn->getPointer();
  const int *connIdxPtr=connIndex->getPointer();
  for(int i=0;i<nbOfElem;i++)
    {
      int delta=connIdxPtr[1]-connIdxPtr[0];
      if(*connPtr==type)
        {
          conn4MEDFile.insert(conn4MEDFile.end(),connPtr+1,connPtr+delta);
          ret++;
        }
      connIdxPtr++;
      connPtr+=delta;
    }
  std::transform(conn4MEDFile.begin(),conn4MEDFile.end(),conn4MEDFile.begin(),std::bind2nd(std::plus<int>(),1));
  return ret;
}

int MEDLoaderNS::buildMEDSubConnectivityOfOneTypesPolyg(DataArrayInt *conn, DataArrayInt *connIndex, std::vector<int>& conn4MEDFile, std::vector<int>& connIndex4MEDFile)
{
  int ret=0;
  int nbOfElem=connIndex->getNbOfElems()-1;
  const int *connPtr=conn->getPointer();
  const int *connIdxPtr=connIndex->getPointer();
  connIndex4MEDFile.push_back(1);
  for(int i=0;i<nbOfElem;i++)
    {
      int delta=connIdxPtr[1]-connIdxPtr[0];
      if(*connPtr==INTERP_KERNEL::NORM_POLYGON)
        {
          conn4MEDFile.insert(conn4MEDFile.end(),connPtr+1,connPtr+delta);
          connIndex4MEDFile.push_back(connIndex4MEDFile.back()+delta-1);
          ret++;
        }
      connIdxPtr++;
      connPtr+=delta;
    }
  std::transform(conn4MEDFile.begin(),conn4MEDFile.end(),conn4MEDFile.begin(),std::bind2nd(std::plus<int>(),1));
  return ret;
}
  
int MEDLoaderNS::buildMEDSubConnectivityOfOneTypesPolyh(DataArrayInt *conn, DataArrayInt *connIndex, std::vector<int>& conn4MEDFile, std::vector<int>& connIndex4MEDFile, std::vector<int>& connIndexRk24MEDFile)
{
  return 0;
}
  
/*!
 * This method builds a sub set of connectivity for a given type 'type'.
 * @param conn input containing connectivity with MEDCoupling format.
 * @param connIndex input containing connectivity index in MEDCoupling format.
 * @param type input specifying which cell types will be extracted in conn4MEDFile. 
 * @param conn4MEDFile output containing the connectivity directly understandable by MEDFile; conn4MEDFile has to be empty before this method called.
 * @param connIndex4MEDFile output containing index connectivity understandable by MEDFile; only used by polygons and polyhedrons (it is face nodal connec).
 * @param connIndexRk24MEDFile output containing index of rank 2 understandable by MEDFile; only used by polyhedrons.
 * @return nb of elements extracted.
 */
int MEDLoaderNS::buildMEDSubConnectivityOfOneType(DataArrayInt *conn, DataArrayInt *connIndex, INTERP_KERNEL::NormalizedCellType type, std::vector<int>& conn4MEDFile,
                                                  std::vector<int>& connIndex4MEDFile, std::vector<int>& connIndexRk24MEDFile)
{
    
  const INTERP_KERNEL::CellModel& cellMod=INTERP_KERNEL::CellModel::getCellModel(type);
  if(!cellMod.isDynamic())
    return buildMEDSubConnectivityOfOneTypeStaticTypes(conn,connIndex,type,conn4MEDFile);
  else
    {
      if(type==INTERP_KERNEL::NORM_POLYGON)
        return buildMEDSubConnectivityOfOneTypesPolyg(conn,connIndex,conn4MEDFile,connIndex4MEDFile);
      else
        return buildMEDSubConnectivityOfOneTypesPolyh(conn,connIndex,conn4MEDFile,connIndex4MEDFile,connIndexRk24MEDFile);
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

ParaMEDMEM::MEDCouplingFieldDouble *MEDLoader::ReadFieldDoubleCell(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order)
{
  return MEDLoaderNS::readFieldDoubleLev1(fileName,meshName,meshDimRelToMax,fieldName,iteration,order,ON_CELLS);
}

ParaMEDMEM::MEDCouplingFieldDouble *MEDLoader::ReadFieldDoubleNode(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order)
{
  return MEDLoaderNS::readFieldDoubleLev1(fileName,meshName,meshDimRelToMax,fieldName,iteration,order,ON_NODES);
}

void MEDLoader::WriteUMesh(const char *fileName, ParaMEDMEM::MEDCouplingUMesh *mesh)
{
  med_idt fid=MEDouvrir((char *)fileName,MED_CREATION);
  std::string meshName(mesh->getName());
  if(meshName=="")
    {
      MEDfermer(fid);
      throw INTERP_KERNEL::Exception("MEDCouplingMesh must have a not null name !");
    }
  char maa[MED_TAILLE_NOM+1];
  strcpy(maa,meshName.c_str());
  MEDmaaCr(fid,maa,mesh->getSpaceDimension(),MED_NON_STRUCTURE,maa);
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
          int nbOfElt=MEDLoaderNS::buildMEDSubConnectivityOfOneType(conn,connIndex,curType,medConn,medConnIndex,medConnIndex2);
          if(curMedType!=MED_POLYGONE && curMedType!=MED_POLYEDRE)
            MEDconnEcr(fid,maa,mesh->getMeshDimension(),&medConn[0],MED_FULL_INTERLACE,nbOfElt,MED_MAILLE,curMedType,MED_NOD);
          else
            {
              if(curMedType==MED_POLYGONE)
                MEDpolygoneConnEcr(fid,maa,&medConnIndex[0],medConnIndex.size(),&medConn[0],MED_MAILLE,MED_NOD);
            }
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
  MEDfermer(fid);
}

void MEDLoader::WriteField(const char *fileName, ParaMEDMEM::MEDCouplingFieldDouble *f)
{
}
