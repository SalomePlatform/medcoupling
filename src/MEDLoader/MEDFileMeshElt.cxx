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

#include "MEDFileMeshElt.hxx"

#include "MEDCouplingUMesh.hxx"

#include "InterpKernelAutoPtr.hxx"
#include "CellModel.hxx"

#include <iostream>

extern med_geometry_type typmai3[32];

using namespace ParaMEDMEM;

MEDFileUMeshPerType *MEDFileUMeshPerType::New(med_idt fid, const char *mName, int dt, int it, int mdim, med_geometry_type geoElt, INTERP_KERNEL::NormalizedCellType geoElt2)
{
  med_entity_type whichEntity;
  if(!isExisting(fid,mName,dt,it,geoElt,whichEntity))
    return 0;
  return new MEDFileUMeshPerType(fid,mName,dt,it,mdim,geoElt,geoElt2,whichEntity);
}

bool MEDFileUMeshPerType::isExisting(med_idt fid, const char *mName, int dt, int it, med_geometry_type geoElt, med_entity_type& whichEntity)
{
  static const med_entity_type entities[3]={MED_CELL,MED_DESCENDING_FACE,MED_DESCENDING_EDGE};
  int nbOfElt=0;
  for(int i=0;i<3;i++)
    {
      med_bool changement,transformation;
      int tmp=MEDmeshnEntity(fid,mName,dt,it,entities[i],geoElt,MED_CONNECTIVITY,MED_NODAL,
                             &changement,&transformation);
      if(tmp>nbOfElt)
        {
          nbOfElt=tmp;
          whichEntity=entities[i];
          if(i>0)
            std::cerr << "WARNING : MEDFile has been detected to be no compilant with MED 3 : Please change entity in MEDFile for geotype " <<  geoElt << std::endl;
        }
    }
  return nbOfElt>0;
}

int MEDFileUMeshPerType::getDim() const
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(_type);
  return cm.getDimension();
}

MEDFileUMeshPerType::MEDFileUMeshPerType(med_idt fid, const char *mName, int dt, int it, int mdim, med_geometry_type geoElt, INTERP_KERNEL::NormalizedCellType type,
                                         med_entity_type entity):_type(type),_entity(entity)
{
  med_bool changement,transformation;
  int curNbOfElem=MEDmeshnEntity(fid,mName,dt,it,entity,geoElt,MED_CONNECTIVITY,MED_NODAL,
                                 &changement,&transformation);
  if(type!=INTERP_KERNEL::NORM_POLYGON && type!=INTERP_KERNEL::NORM_POLYHED)
    {
      loadFromStaticType(fid,mName,dt,it,mdim,curNbOfElem,geoElt,type,entity);
      return;
    }
  if(type==INTERP_KERNEL::NORM_POLYGON)
    {
      loadPolyg(fid,mName,dt,it,mdim,curNbOfElem,geoElt,entity);
      return;
    }
  //if(type==INTERP_KERNEL::NORM_POLYHED)
  loadPolyh(fid,mName,dt,it,mdim,curNbOfElem,geoElt,entity);
}

void MEDFileUMeshPerType::loadFromStaticType(med_idt fid, const char *mName, int dt, int it, int mdim, int curNbOfElem, med_geometry_type geoElt, INTERP_KERNEL::NormalizedCellType type,
                                             med_entity_type entity)
{
  _conn=DataArrayInt::New();
  int nbOfNodesPerCell=(geoElt%100);
  _conn->alloc((nbOfNodesPerCell+1)*curNbOfElem,1);
  _conn_index=DataArrayInt::New();
  _conn_index->alloc(curNbOfElem+1,1);
  INTERP_KERNEL::AutoPtr<int> connTab=new int[(nbOfNodesPerCell)*curNbOfElem];
  MEDmeshElementConnectivityRd(fid,mName,dt,it,entity,geoElt,MED_NODAL,MED_FULL_INTERLACE,connTab);
  loadCommonPart(fid,mName,dt,it,mdim,curNbOfElem,geoElt,entity);
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
}

void MEDFileUMeshPerType::loadCommonPart(med_idt fid, const char *mName, int dt, int it, int mdim, int curNbOfElem, med_geometry_type geoElt,
                                         med_entity_type entity)
{
  med_bool changement,transformation;
  _fam=DataArrayInt::New();
  _fam->alloc(curNbOfElem,1);
  if(MEDmeshnEntity(fid,mName,dt,it,entity,geoElt,MED_FAMILY_NUMBER,MED_NODAL,&changement,&transformation)>0)
    {
      if(MEDmeshEntityFamilyNumberRd(fid,mName,dt,it,entity,geoElt,_fam->getPointer())!=0)
        std::fill(_fam->getPointer(),_fam->getPointer()+curNbOfElem,0);
    }
  else
    std::fill(_fam->getPointer(),_fam->getPointer()+curNbOfElem,0);
  if(MEDmeshnEntity(fid,mName,dt,it,entity,geoElt,MED_NUMBER,MED_NODAL,&changement,&transformation)>0)
    {
      _num=DataArrayInt::New();
      _num->alloc(curNbOfElem,1);
      if(MEDmeshEntityNumberRd(fid,mName,dt,it,entity,geoElt,_num->getPointer())!=0)
        _num=0;
    }
  else
    _num=0;
  if(MEDmeshnEntity(fid,mName,dt,it,entity,geoElt,MED_NAME,MED_NODAL,&changement,&transformation)>0)
    {
      _names=DataArrayAsciiChar::New();
      _names->alloc(curNbOfElem+1,MED_SNAME_SIZE);//not a bug to avoid the memory corruption due to last \0 at the end
      if(MEDmeshEntityNameRd(fid,mName,dt,it,entity,geoElt,_names->getPointer())!=0)
        _names=0;
      else
        _names->reAlloc(curNbOfElem);//not a bug to avoid the memory corruption due to last \0 at the end
    }
}

void MEDFileUMeshPerType::loadPolyg(med_idt fid, const char *mName, int dt, int it, int mdim, int arraySize, med_geometry_type geoElt,
                                    med_entity_type entity)
{
  med_bool changement,transformation;
  med_int curNbOfElem=MEDmeshnEntity(fid,mName,dt,it,entity,MED_POLYGON,MED_INDEX_NODE,MED_NODAL,&changement,&transformation)-1;
  _conn_index=DataArrayInt::New();
  _conn_index->alloc(curNbOfElem+1,1);
  _conn=DataArrayInt::New();
  _conn->alloc(arraySize+curNbOfElem,1);
  INTERP_KERNEL::AutoPtr<int> locConn=new int[arraySize];
  MEDmeshPolygonRd(fid,mName,dt,it,MED_CELL,MED_NODAL,_conn_index->getPointer(),locConn);
  int *w1=_conn->getPointer();
  int *w2=_conn_index->getPointer();
  const int *wi=locConn;
  for(int i=0;i<curNbOfElem;i++,w2++)
    {
      *w1++=(int)INTERP_KERNEL::NORM_POLYGON;
      const int *wi2=wi+(w2[1]-w2[0]);
      w1=std::transform(wi,wi2,w1,std::bind2nd(std::plus<int>(),-1));
      wi=wi2;
      *w2=*w2-1+i;
    }
  *w2=*w2-1+curNbOfElem;
  loadCommonPart(fid,mName,dt,it,mdim,curNbOfElem,MED_POLYGON,entity);
}

void MEDFileUMeshPerType::loadPolyh(med_idt fid, const char *mName, int dt, int it, int mdim, int connFaceLgth, med_geometry_type geoElt,
                                    med_entity_type entity)
{
  med_bool changement,transformation;
  med_int indexFaceLgth=MEDmeshnEntity(fid,mName,dt,it,MED_CELL,MED_POLYHEDRON,MED_INDEX_NODE,MED_NODAL,&changement,&transformation);
  int curNbOfElem=MEDmeshnEntity(fid,mName,dt,it,MED_CELL,MED_POLYHEDRON,MED_INDEX_FACE,MED_NODAL,&changement,&transformation)-1;
  INTERP_KERNEL::AutoPtr<int> index=new int[curNbOfElem+1];
  INTERP_KERNEL::AutoPtr<int> indexFace=new int[indexFaceLgth];
  INTERP_KERNEL::AutoPtr<int> locConn=new int[connFaceLgth];
  MEDmeshPolyhedronRd(fid,mName,dt,it,MED_CELL,MED_NODAL,index,indexFace,locConn);
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
  loadCommonPart(fid,mName,dt,it,mdim,curNbOfElem,MED_POLYHEDRON,entity);
}

void MEDFileUMeshPerType::write(med_idt fid, const char *mname, int mdim, const MEDCouplingUMesh *m, const DataArrayInt *fam, const DataArrayInt *num, const DataArrayAsciiChar *names)
{
  int nbOfCells=m->getNumberOfCells();
  if(nbOfCells<1)
    return ;
  int dt,it;
  double timm=m->getTime(dt,it);
  INTERP_KERNEL::NormalizedCellType ikt=m->getTypeOfCell(0);
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(ikt);
  med_geometry_type curMedType=typmai3[(int)ikt];
  const int *conn=m->getNodalConnectivity()->getConstPointer();
  const int *connI=m->getNodalConnectivityIndex()->getConstPointer();
  if(ikt!=INTERP_KERNEL::NORM_POLYGON && ikt!=INTERP_KERNEL::NORM_POLYHED)
    {
      int nbNodesPerCell=cm.getNumberOfNodes();
      INTERP_KERNEL::AutoPtr<int> tab=new int[nbNodesPerCell*nbOfCells];
      int *w=tab;
      for(int i=0;i<nbOfCells;i++)
        w=std::transform(conn+connI[i]+1,conn+connI[i+1],w,std::bind2nd(std::plus<int>(),1));
      MEDmeshElementConnectivityWr(fid,mname,dt,it,timm,MED_CELL,curMedType,MED_NODAL,MED_FULL_INTERLACE,nbOfCells,tab);
    }
  else
    {
      if(ikt==INTERP_KERNEL::NORM_POLYGON)
        {
          INTERP_KERNEL::AutoPtr<int> tab1=new int[nbOfCells+1];
          INTERP_KERNEL::AutoPtr<int> tab2=new int[m->getMeshLength()];
          int *wI=tab1; *wI=1;
          int *w=tab2;
          for(int i=0;i<nbOfCells;i++,wI++)
            {
              wI[1]=wI[0]+connI[i+1]-connI[i]-1;
              w=std::transform(conn+connI[i]+1,conn+connI[i+1],w,std::bind2nd(std::plus<int>(),1));
            }
          MEDmeshPolygonWr(fid,mname,dt,it,timm,MED_CELL,MED_NODAL,nbOfCells+1,tab1,tab2);
        }
      else
        {
          int meshLgth=m->getMeshLength();
          int nbOfFaces=std::count(conn,conn+meshLgth,-1)+nbOfCells;
          INTERP_KERNEL::AutoPtr<int> tab1=new int[nbOfCells+1];
          int *w1=tab1; *w1=1;
          INTERP_KERNEL::AutoPtr<int> tab2=new int[nbOfFaces+1];
          int *w2=tab2; *w2=1;
          INTERP_KERNEL::AutoPtr<int> bigtab=new int[meshLgth-nbOfCells];
          int *bt=bigtab;
          for(int i=0;i<nbOfCells;i++,w1++)
            {
              int nbOfFaces2=0;
              for(const int *w=conn+connI[i]+1;w!=conn+connI[i+1];w2++)
                {
                  const int *wend=std::find(w,conn+connI[i+1],-1);
                  bt=std::transform(w,wend,bt,std::bind2nd(std::plus<int>(),1));
                  int nbOfNode=std::distance(w,wend);
                  w2[1]=w2[0]+nbOfNode;
                  if(wend!=conn+connI[i+1])
                    w=wend+1;
                  else
                    w=wend;
                  nbOfFaces2++;
                }
              w1[1]=w1[0]+nbOfFaces2;
            }
          MEDmeshPolyhedronWr(fid,mname,dt,it,timm,MED_CELL,MED_NODAL,nbOfCells+1,tab1,nbOfFaces+1,tab2,bigtab);
        }
    }
  if(fam)
    MEDmeshEntityFamilyNumberWr(fid,mname,dt,it,MED_CELL,curMedType,nbOfCells,fam->getConstPointer());
  if(num)
    MEDmeshEntityNumberWr(fid,mname,dt,it,MED_CELL,curMedType,nbOfCells,num->getConstPointer());
  if(names)
    {
      if(names->getNumberOfComponents()!=MED_SNAME_SIZE)
        {
          std::ostringstream oss; oss << "MEDFileUMeshPerType::write : expected a name field on cells with number of components set to " << MED_SNAME_SIZE;
          oss << " ! The array has " << names->getNumberOfComponents() << " components !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      MEDmeshEntityNameWr(fid,mname,dt,it,MED_CELL,curMedType,nbOfCells,names->getConstPointer());
    }
}
