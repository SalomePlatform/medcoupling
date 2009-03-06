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
#include "ParaMESH.hxx"
#include "BlockTopology.hxx"
#include "MEDCouplingUMesh.hxx"

extern "C"
{
#include "med.h"
}

#include <string>
#include <cstring>
#include <sstream>
#include <fstream>

med_geometrie_element typmai[MED_NBR_GEOMETRIE_MAILLE] = { MED_POINT1,
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
                                                           MED_PYRA13 };

INTERP_KERNEL::NormalizedCellType typmai2[MED_NBR_GEOMETRIE_MAILLE] = { INTERP_KERNEL::NORM_ERROR,
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
                                                                        INTERP_KERNEL::NORM_PYRA13 };

using namespace ParaMEDMEM;

const char WHITE_SPACES[]=" \n";

MEDLoader::MEDConnOfOneElemType::MEDConnOfOneElemType(INTERP_KERNEL::NormalizedCellType type, int *conn, int lgth):_lgth(lgth),
                                                                                                                   _conn(conn),_global(0),
                                                                                                                   _type(type)
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
  delete [] _conn;
  delete [] _global;
}

std::string buildStringFromFortran(const char *expr, int lgth)
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

namespace MEDLoader
{
  med_int getIdFromMeshName(med_idt fid, const char *meshName) throw(INTERP_KERNEL::Exception)
  {
    if(meshName==0)
      return 1;
    med_int n=MEDnMaa(fid);
    if(n==0)
      throw INTERP_KERNEL::Exception("No mesh in file.");
    med_maillage type_maillage;
    char maillage_description[MED_TAILLE_DESC+1];
    med_int dim;
    char nommaa[MED_TAILLE_NOM+1];
    std::ostringstream os;
    for(med_int i=1;i<=n;i++)
      {
        MEDmaaInfo(fid,i,nommaa,&dim,&type_maillage,maillage_description);
        std::string cur=buildStringFromFortran(nommaa,sizeof(nommaa));
        if(cur==meshName)
          return i;
        os << "\'" << cur.c_str() << "\' "; 
      }
    std::ostringstream os2;
    os2 << "MeshName '" << meshName << "' not in file : meshes available : " << os.str();
    throw INTERP_KERNEL::Exception(os2.str().c_str());
  }

  void readUMeshDataInMedFile(med_idt fid, med_int meshId, double *&coords, int& nCoords, int& spaceDim, std::list<MEDLoader::MEDConnOfOneElemType>& conn)
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
        if(curNbOfElemF>curNbOfElemM)
          {
            curNbOfElem=curNbOfElemF;
            whichEntity=MED_FACE;
          }
        else
          {
            curNbOfElem=curNbOfElemM;
            whichEntity=MED_MAILLE;
          }
        if(curNbOfElem>0)
          {
            int *connTab=new int[(curMedType%100)*curNbOfElem];
            MEDLoader::MEDConnOfOneElemType elem(typmai2[i],connTab,curNbOfElem);
            int *tmp=new int[curNbOfElem];
            int *fam=new int[curNbOfElem];
            char *noms=new char[MED_TAILLE_PNOM*curNbOfElem+1];
            MEDelementsLire(fid,nommaa,Mdim,connTab,MED_FULL_INTERLACE,noms,&inoele,tmp,&inuele,fam,curNbOfElem,whichEntity,curMedType,MED_NOD);
            delete [] tmp;
            delete [] fam;
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
  }
}

unsigned MEDLoader::calculateHighestMeshDim(const std::list<MEDLoader::MEDConnOfOneElemType>& conn)
{
  unsigned ret=0;
  for(std::list<MEDLoader::MEDConnOfOneElemType>::const_iterator iter=conn.begin();iter!=conn.end();iter++)
    {
      unsigned curDim=INTERP_KERNEL::CellModel::getCellModel((*iter).getType()).getDimension();
      if(ret<curDim)
        ret=curDim;
    }
  return ret;
}

void MEDLoader::keepSpecifiedMeshDim(std::list<MEDLoader::MEDConnOfOneElemType>& conn, unsigned meshDim)
{
  for(std::list<MEDLoader::MEDConnOfOneElemType>::iterator iter=conn.begin();iter!=conn.end();)
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

void MEDLoader::tradMEDFileCoreFrmt2MEDCouplingUMesh(const std::list<MEDLoader::MEDConnOfOneElemType>& medConnFrmt,
                                                     DataArrayInt* &conn,
                                                     DataArrayInt* &connIndex)
{
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
      totalNbOfCells+=(*iter).getLength();
      if(!cellMod.isDynamic())
        totalNbOfMedConn+=(*iter).getLength()*cellMod.getNumberOfNodes();
      else
        throw INTERP_KERNEL::Exception("Polyg/polh not implemented yet !");
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
      int *sourceConn=(*iter).getArray();
      const INTERP_KERNEL::CellModel& cellMod=INTERP_KERNEL::CellModel::getCellModel(type);
      int nbOfCellsInCurType;
      int nbOfNodesIn1Cell=cellMod.getNumberOfNodes();
      if(!cellMod.isDynamic())
        nbOfCellsInCurType=(*iter).getLength();
      else
        throw INTERP_KERNEL::Exception("Polyg/polh not implemented yet !");
      if(!cellMod.isDynamic())
        {
          for(int i=0;i<nbOfCellsInCurType;i++,connIdxPtr++)
            {
              *connIdxPtr=connFillId;
              *connPtr++=type;
              connPtr=std::transform(sourceConn,sourceConn+nbOfNodesIn1Cell,connPtr,std::bind2nd(std::minus<int>(),1));
              connFillId+=nbOfNodesIn1Cell+1;
              sourceConn+=nbOfNodesIn1Cell;
            }
          *connIdxPtr=connFillId;
        }
    }
}

void MEDLoader::releaseMEDFileCoreFrmt(std::list<MEDLoader::MEDConnOfOneElemType>& medConnFrmt)
{
  for(std::list<MEDLoader::MEDConnOfOneElemType>::iterator iter=medConnFrmt.begin();iter!=medConnFrmt.end();iter++)
    (*iter).releaseArray();
  medConnFrmt.clear();
}

/*!
 * This method builds a sub set of connectivity for a given type 'type'.
 * @param conn input containing connectivity with MEDCoupling format.
 * @param connIndex input containing connectivity index in MEDCoupling format.
 * @param type input specifying which cell types will be extracted in conn4MEDFile. 
 * @param conn4MEDFile output containing the connectivity directly understandable by MEDFile; conn4MEDFile has to be empty before this method called.
 * @return nb of elements extracted.
 */
int MEDLoader::buildMEDSubConnectivityOfOneType(DataArrayInt *conn, DataArrayInt *connIndex, INTERP_KERNEL::NormalizedCellType type, std::vector<int>& conn4MEDFile)
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

MEDCouplingUMesh *MEDLoader::ReadUMeshFromFile(const char *fileName, const char *meshName, int meshDimRelToMax) throw(INTERP_KERNEL::Exception)
{
  //Extraction data from MED file.
  med_idt fid=MEDouvrir((char *)fileName,MED_LECTURE);
  med_int mid=getIdFromMeshName(fid,meshName);
  unsigned meshDimExtract;
  double *coords;
  int nCoords;
  int spaceDim;
  std::list<MEDLoader::MEDConnOfOneElemType> conn;
  readUMeshDataInMedFile(fid,mid,coords,nCoords,spaceDim,conn);
  meshDimExtract=calculateHighestMeshDim(conn);
  meshDimExtract=meshDimExtract+meshDimRelToMax;
  keepSpecifiedMeshDim(conn,meshDimExtract);
  MEDfermer(fid);
  //Put data in returned data structure.
  MEDCouplingUMesh *ret=MEDCouplingUMesh::New();
  ret->setName(meshName);
  ret->setMeshDimension(meshDimExtract);
  //
  DataArrayDouble *coordsArr=DataArrayDouble::New();
  coordsArr->useArray(coords,true,ParaMEDMEM::CPP_DEALLOC,nCoords,spaceDim);
  ret->setCoords(coordsArr);
  coordsArr->decrRef();
  //
  DataArrayInt *connArr,*connIndexArr;
  tradMEDFileCoreFrmt2MEDCouplingUMesh(conn,connArr,connIndexArr);
  ret->setConnectivity(connArr,connIndexArr);
  //clean-up
  if(connArr)
    connArr->decrRef();
  if(connIndexArr)
    connIndexArr->decrRef();
  releaseMEDFileCoreFrmt(conn);
  return ret;
}

void MEDLoader::writeUMesh(const char *fileName, ParaMEDMEM::MEDCouplingUMesh *mesh)
{
  med_idt fid=MEDouvrir((char *)fileName,MED_CREATION);
  char maa[MED_TAILLE_NOM+1];
  std::fill(maa,maa+MED_TAILLE_NOM+1,'\0');
  const char *meshName=mesh->getName();
  strcpy(maa,meshName);
  MEDmaaCr(fid,maa,mesh->getMeshDimension(),MED_NON_STRUCTURE,maa);
  std::set<INTERP_KERNEL::NormalizedCellType> allTypes(mesh->getAllTypes());
  DataArrayInt *conn=mesh->getNodalConnectivity();
  DataArrayInt *connIndex=mesh->getNodalConnectivityIndex();
  char familyName[MED_TAILLE_NOM+1];
  std::fill(familyName,familyName+MED_TAILLE_NOM+1,'\0');
  const char DftFamilyName[]="DftFamily";
  std::copy(DftFamilyName,DftFamilyName+sizeof(DftFamilyName),familyName);
  for(int i=0;i<MED_NBR_GEOMETRIE_MAILLE;i++)
    {
      med_geometrie_element curMedType=typmai[i];
      INTERP_KERNEL::NormalizedCellType curType=typmai2[i];
      if(allTypes.find(curType)!=allTypes.end())
        {
          std::vector<int> medConn;
          int nbOfElt=buildMEDSubConnectivityOfOneType(conn,connIndex,curType,medConn);
          MEDconnEcr(fid,maa,mesh->getMeshDimension(),&medConn[0],MED_FULL_INTERLACE,nbOfElt,MED_MAILLE,curMedType,MED_NOD);
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

/*!
 * This method builds the master file 'fileName' of a parallel MED file defined in 'fileNames'.
 */
void MEDLoader::writeMasterFile(const char *fileName, const std::vector<std::string>& fileNames, const char *meshName)
{
  int nbOfDom=fileNames.size();
  std::ofstream fs(fileName);
  fs << "#MED Fichier V 2.3" << " " << std::endl;
  fs << "#"<<" " << std::endl;
  fs << nbOfDom <<" " << std::endl;
  for(int i=0;i<nbOfDom;i++)
    fs << meshName << " " << i+1 << " " << meshName << "_" << i+1 << " localhost " << fileNames[i] << " " << std::endl;
}

void MEDLoader::writeParaMesh(const char *fileName, ParaMEDMEM::ParaMESH *mesh)
{
  if(!mesh->getBlockTopology()->getProcGroup()->containsMyRank())
    return ;
  int myRank=mesh->getBlockTopology()->getProcGroup()->myRank();
  int nbDomains=mesh->getBlockTopology()->getProcGroup()->size();
  std::vector<std::string> fileNames(nbDomains);
  for(int i=0;i<nbDomains;i++)
    {
      std::ostringstream sstr;
      sstr << fileName << i+1 << ".med";
      fileNames[i]=sstr.str();
    }
  if(myRank==0)
    writeMasterFile(fileName,fileNames,mesh->getCellMesh()->getName());
  writeUMesh(fileNames[myRank].c_str(),mesh->getCellMesh());
}

void MEDLoader::writeParaField(const char *fileName, const char *meshName, ParaMEDMEM::ParaFIELD *f)
{
}
