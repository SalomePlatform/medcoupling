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

#include "MEDFileMeshElt.hxx"

#include "CellModel.hxx"

using namespace ParaMEDMEM;

MEDFileUMeshPerType *MEDFileUMeshPerType::New(med_idt fid, const char *mName, int mdim, med_geometrie_element geoElt, INTERP_KERNEL::NormalizedCellType geoElt2)
{
  med_entite_maillage whichEntity;
  if(!isExisting(fid,mName,geoElt,whichEntity))
    return 0;
  return new MEDFileUMeshPerType(fid,mName,mdim,geoElt,geoElt2,whichEntity);
}

bool MEDFileUMeshPerType::isExisting(med_idt fid, const char *mName, med_geometrie_element geoElt, med_entite_maillage& whichEntity)
{
  static const med_entite_maillage entities[3]={MED_MAILLE,MED_FACE,MED_ARETE};
  int nbOfElt=0;
  for(int i=0;i<3;i++)
    {
      int tmp=MEDnEntMaa(fid,(char *)mName,MED_CONN,entities[i],geoElt,MED_NOD);
      if(tmp>nbOfElt)
        {
          nbOfElt=tmp;
          whichEntity=entities[i];
        }
    }
  return nbOfElt>0;
}

int MEDFileUMeshPerType::getDim() const
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::getCellModel(_type);
  return cm.getDimension();
}

void MEDFileUMeshPerType::write(med_idt fid) const
{
}

MEDFileUMeshPerType::MEDFileUMeshPerType(med_idt fid, const char *mName, int mdim, med_geometrie_element geoElt, INTERP_KERNEL::NormalizedCellType type,
                                         med_entite_maillage entity):_type(type),_entity(entity)
{
  int curNbOfElem=MEDnEntMaa(fid,(char *)mName,MED_CONN,entity,geoElt,MED_NOD);
  if(type!=INTERP_KERNEL::NORM_POLYGON && type!=INTERP_KERNEL::NORM_POLYHED)
    {
      loadFromStaticType(fid,mName,mdim,curNbOfElem,geoElt,type,entity);
      return;
    }
  if(type==INTERP_KERNEL::NORM_POLYGON)
    {
      loadPolyg(fid,mName,mdim,curNbOfElem,geoElt,entity);
      return;
    }
  //if(type==INTERP_KERNEL::NORM_POLYHED)
  loadPolyh(fid,mName,mdim,curNbOfElem,geoElt,entity);
}

void MEDFileUMeshPerType::loadFromStaticType(med_idt fid, const char *mName, int mdim, int curNbOfElem, med_geometrie_element geoElt, INTERP_KERNEL::NormalizedCellType type,
                                             med_entite_maillage entity)
{
  _conn=DataArrayInt::New();
  int nbOfNodesPerCell=(geoElt%100);
  _conn->alloc((nbOfNodesPerCell+1)*curNbOfElem,1);
  _conn_index=DataArrayInt::New();
  _conn_index->alloc(curNbOfElem+1,1);
  int *connTab=new int[(nbOfNodesPerCell)*curNbOfElem];
  _num=DataArrayInt::New();
  _num->alloc(curNbOfElem,1);
  char *noms=new char[MED_TAILLE_PNOM*curNbOfElem+1];
  med_booleen inoele, inuele;
  _fam=DataArrayInt::New();
  _fam->alloc(curNbOfElem,1);
  MEDelementsLire(fid,(char *)mName,mdim,connTab,MED_FULL_INTERLACE,noms,&inoele,_num->getPointer(),&inuele,_fam->getPointer(),curNbOfElem,entity,geoElt,MED_NOD);
  delete [] noms;
  int *w1=_conn->getPointer();
  int *w2=_conn_index->getPointer();
  *w2++=0;
  const int *wi=connTab;
  for(int i=0;i<curNbOfElem;i++,wi+=nbOfNodesPerCell,w2++)
    {
      *w1++=(int)type;
      w1=std::transform(wi,wi+nbOfNodesPerCell,w1,std::bind2nd(std::plus<int>(),-1));
      *w2=w2[-1]+nbOfNodesPerCell+1;
    }
  delete [] connTab;
  if(!inuele)
    _num=0;
}

void MEDFileUMeshPerType::loadPolyg(med_idt fid, const char *mName, int mdim, int curNbOfElem, med_geometrie_element geoElt,
                                    med_entite_maillage entity)
{
  med_int arraySize;
  MEDpolygoneInfo(fid,(char *)mName,entity,MED_NOD,&arraySize);
  _conn_index=DataArrayInt::New();
  _conn_index->alloc(curNbOfElem+1,1);
  _conn=DataArrayInt::New();
  _conn->alloc(arraySize+curNbOfElem,1);
  _num=DataArrayInt::New();
  _num->alloc(curNbOfElem,1);
  _fam=DataArrayInt::New();
  _fam->alloc(curNbOfElem,1);
  int *locConn=new int[arraySize];
  MEDpolygoneConnLire(fid,(char *)mName,_conn_index->getPointer(),curNbOfElem+1,locConn,entity,MED_NOD);
  int *w1=_conn->getPointer();
  int *w2=_conn_index->getPointer();
  const int *wi=locConn;
  for(int i=0;i<curNbOfElem;i++,w2++)
    {
      *w1++=(int)INTERP_KERNEL::NORM_POLYGON;
      w1=std::transform(wi,wi+(w2[1]-w2[0]),w1,std::bind2nd(std::plus<int>(),-1));
      *w2=*w2-1+i;
    }
  *w2=*w2-1+curNbOfElem;
  delete [] locConn;
  MEDfamLire(fid,(char *)mName,_fam->getPointer(),curNbOfElem,entity,MED_POLYGONE);
  if(MEDnumLire(fid,(char *)mName,_num->getPointer(),curNbOfElem,entity,MED_POLYGONE)!=0)
    _num=0;
}

void MEDFileUMeshPerType::loadPolyh(med_idt fid, const char *mName, int mdim, int curNbOfElem, med_geometrie_element geoElt,
                                    med_entite_maillage entity)
{
  med_int indexFaceLgth,connFaceLgth;
  MEDpolyedreInfo(fid,(char*)mName,MED_NOD,&indexFaceLgth,&connFaceLgth);
  int *index=new int[curNbOfElem+1];
  int *indexFace=new int[indexFaceLgth];
  int *locConn=new int[connFaceLgth];
  _fam=DataArrayInt::New();
  _fam->alloc(curNbOfElem,1);
  MEDpolyedreConnLire(fid,(char *)mName,index,curNbOfElem+1,indexFace,indexFaceLgth,locConn,MED_NOD);
  MEDfamLire(fid,(char *)mName,_fam->getPointer(),curNbOfElem,MED_MAILLE,MED_POLYEDRE);
  int arraySize=connFaceLgth;
  for(int i=0;i<curNbOfElem;i++)
    arraySize+=index[i+1]-index[i]-1;
  _conn=DataArrayInt::New();
  _conn->alloc(arraySize+curNbOfElem,1);
  int *wFinalConn=_conn->getPointer();
  _conn_index=DataArrayInt::New();
  _conn_index->alloc(curNbOfElem+1,1);
  int *finalIndex=_conn_index->getPointer();
  finalIndex[0]=0;
  for(int i=0;i<curNbOfElem;i++)
    {
      *wFinalConn++=(int)INTERP_KERNEL::NORM_POLYHED;
      finalIndex[i+1]=finalIndex[i]+index[i+1]-index[i]-1+indexFace[index[i+1]-1]-indexFace[index[i]-1]+1;
      wFinalConn=std::transform(locConn+indexFace[index[i]-1]-1,locConn+indexFace[index[i]]-1,wFinalConn,std::bind2nd(std::plus<int>(),-1));
      for(int j=index[i];j<index[i+1]-1;j++)
        {
          *wFinalConn++=-1;
          wFinalConn=std::transform(locConn+indexFace[j]-1,locConn+indexFace[j+1]-1,wFinalConn,std::bind2nd(std::plus<int>(),-1));
        }
    }
  delete [] index;
  delete [] locConn;
  delete [] indexFace;
  _num=DataArrayInt::New();
  _num->alloc(curNbOfElem,1);
  if(MEDnumLire(fid,(char *)mName,_num->getPointer(),curNbOfElem,MED_MAILLE,MED_POLYEDRE)!=0)
    _num=0;
}
