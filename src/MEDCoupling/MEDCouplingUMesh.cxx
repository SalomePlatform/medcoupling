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
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "CellModel.hxx"
#include "VolSurfUser.txx"

#include <sstream>
#include <limits>
#include <list>

using namespace ParaMEDMEM;

const char MEDCouplingUMesh::PART_OF_NAME[]="PartOf_";

MEDCouplingUMesh *MEDCouplingUMesh::New()
{
  return new MEDCouplingUMesh;
}

MEDCouplingUMesh *MEDCouplingUMesh::New(const char *meshName, int meshDim)
{
  MEDCouplingUMesh *ret=new MEDCouplingUMesh;
  ret->setName(meshName);
  ret->setMeshDimension(meshDim);
  return ret;
}

MEDCouplingUMesh *MEDCouplingUMesh::clone(bool recDeepCpy) const
{
  return new MEDCouplingUMesh(*this,recDeepCpy);
}

void MEDCouplingUMesh::updateTime()
{
  MEDCouplingPointSet::updateTime();
  if(_nodal_connec)
    {
      updateTimeWith(*_nodal_connec);
    }
  if(_nodal_connec_index)
    {
      updateTimeWith(*_nodal_connec_index);
    }
}

MEDCouplingUMesh::MEDCouplingUMesh():_iterator(-1),_mesh_dim(-2),
                                     _nodal_connec(0),_nodal_connec_index(0)
{
}

void MEDCouplingUMesh::checkCoherency() const throw(INTERP_KERNEL::Exception)
{
  if(_mesh_dim<-1)
    throw INTERP_KERNEL::Exception("No mesh dimension specified !");
  for(std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator iter=_types.begin();iter!=_types.end();iter++)
    {
      if((int)INTERP_KERNEL::CellModel::getCellModel(*iter).getDimension()!=_mesh_dim)
        {
          std::ostringstream message;
          message << "Mesh invalid because dimension is " << _mesh_dim << " and there is presence of cell(s) with type " << (*iter);
          throw INTERP_KERNEL::Exception(message.str().c_str());
        }
    }
}

void MEDCouplingUMesh::setMeshDimension(int meshDim)
{
  if(meshDim<-1)
    throw INTERP_KERNEL::Exception("Invalid meshDim specified ! Must be greater or equal to -1 !");
  _mesh_dim=meshDim;
  declareAsNew();
}

void MEDCouplingUMesh::allocateCells(int nbOfCells)
{
  if(_nodal_connec_index)
    {
      _nodal_connec_index->decrRef();
    }
  if(_nodal_connec)
    {
      _nodal_connec->decrRef();
    }

  _nodal_connec_index=DataArrayInt::New();
  _nodal_connec_index->alloc(nbOfCells+1,1);
  int *pt=_nodal_connec_index->getPointer();
  pt[0]=0;
  _nodal_connec=DataArrayInt::New();
  _nodal_connec->alloc(2*nbOfCells,1);
  _iterator=0;
  _types.clear();
  declareAsNew();
}

void MEDCouplingUMesh::insertNextCell(INTERP_KERNEL::NormalizedCellType type, int size, const int *nodalConnOfCell)
{
  int *pt=_nodal_connec_index->getPointer();
  int idx=pt[_iterator];

  _nodal_connec->writeOnPlace(idx,type,nodalConnOfCell,size);
  _types.insert(type);
  pt[++_iterator]=idx+size+1;
}

void MEDCouplingUMesh::finishInsertingCells()
{
  const int *pt=_nodal_connec_index->getConstPointer();
  int idx=pt[_iterator];

  _nodal_connec->reAlloc(idx);
  _nodal_connec_index->reAlloc(_iterator+1);
  _iterator=-1;
}

bool MEDCouplingUMesh::isEqual(const MEDCouplingMesh *other, double prec) const
{
  //checkFullyDefined();
  const MEDCouplingUMesh *otherC=dynamic_cast<const MEDCouplingUMesh *>(other);
  if(!otherC)
    return false;
  //otherC->checkFullyDefined();
  if(!MEDCouplingMesh::isEqual(other,prec))
    return false;
  if(_mesh_dim!=otherC->_mesh_dim)
    return false;
  if(_types!=otherC->_types)
    return false;
  if(!areCoordsEqual(*otherC,prec))
    return false;
  if(_nodal_connec!=0 || otherC->_nodal_connec!=0)
    if(_nodal_connec==0 || otherC->_nodal_connec==0)
      return false;
  if(_nodal_connec!=otherC->_nodal_connec)
    if(!_nodal_connec->isEqual(*otherC->_nodal_connec))
      return false;
  if(_nodal_connec_index!=0 || otherC->_nodal_connec_index!=0)
    if(_nodal_connec_index==0 || otherC->_nodal_connec_index==0)
      return false;
  if(_nodal_connec_index!=otherC->_nodal_connec_index)
    if(!_nodal_connec_index->isEqual(*otherC->_nodal_connec_index))
      return false;
  return true;
}

/*!
 * \b WARNING this method do the assumption that connectivity lies on the coordinates set.
 * For speed reasons no check of this will be done.
 */
void MEDCouplingUMesh::getReverseNodalConnectivity(DataArrayInt *revNodal, DataArrayInt *revNodalIndx) const
{
  checkFullyDefined();
  int nbOfNodes=getNumberOfNodes();
  int *revNodalIndxPtr=new int[nbOfNodes+1];
  revNodalIndx->useArray(revNodalIndxPtr,true,CPP_DEALLOC,nbOfNodes+1,1);
  std::fill(revNodalIndxPtr,revNodalIndxPtr+nbOfNodes+1,0);
  const int *conn=_nodal_connec->getConstPointer();
  const int *connIndex=_nodal_connec_index->getConstPointer();
  int nbOfCells=getNumberOfCells();
  int nbOfEltsInRevNodal=0;
  for(int eltId=0;eltId<nbOfCells;eltId++)
    {
      const int *strtNdlConnOfCurCell=conn+connIndex[eltId]+1;
      const int *endNdlConnOfCurCell=conn+connIndex[eltId+1];
      for(const int *iter=strtNdlConnOfCurCell;iter!=endNdlConnOfCurCell;iter++)
        if(*iter>=0)//for polyhedrons
          {
            nbOfEltsInRevNodal++;
            revNodalIndxPtr[(*iter)+1]++;
          }
    }
  std::transform(revNodalIndxPtr+1,revNodalIndxPtr+nbOfNodes+1,revNodalIndxPtr,revNodalIndxPtr+1,std::plus<int>());
  int *revNodalPtr=new int[nbOfEltsInRevNodal];
  revNodal->useArray(revNodalPtr,true,CPP_DEALLOC,nbOfEltsInRevNodal,1);
  std::fill(revNodalPtr,revNodalPtr+nbOfEltsInRevNodal,-1);
  for(int eltId=0;eltId<nbOfCells;eltId++)
    {
      const int *strtNdlConnOfCurCell=conn+connIndex[eltId]+1;
      const int *endNdlConnOfCurCell=conn+connIndex[eltId+1];
      for(const int *iter=strtNdlConnOfCurCell;iter!=endNdlConnOfCurCell;iter++)
        if(*iter>=0)//for polyhedrons
          *std::find_if(revNodalPtr+revNodalIndxPtr[*iter],revNodalPtr+revNodalIndxPtr[*iter+1],std::bind2nd(std::equal_to<int>(),-1))=eltId;
    }
}

/*!
 * \b WARNING this method do the assumption that connectivity lies on the coordinates set.
 * For speed reasons no check of this will be done.
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildDescendingConnectivity(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx) const
{
  checkFullyDefined();
  int nbOfCells=getNumberOfCells();
  int nbOfNodes=getNumberOfNodes();
  const int *conn=_nodal_connec->getConstPointer();
  const int *connIndex=_nodal_connec_index->getConstPointer();
  std::vector< std::vector<int> > descMeshConnB(nbOfCells);
  std::vector< std::vector<int> > revDescMeshConnB;
  std::vector< std::vector<int> > revNodalB(nbOfNodes);
  std::vector<int> meshDM1Conn;
  std::vector<int> meshDM1ConnIndex(1); meshDM1ConnIndex[0]=0;
  std::vector<int> meshDM1Type;
  for(int eltId=0;eltId<nbOfCells;eltId++)
    {
      int pos=connIndex[eltId];
      int posP1=connIndex[eltId+1];
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::getCellModel((INTERP_KERNEL::NormalizedCellType)conn[pos]);
      unsigned nbOfSons=cm.getNumberOfSons2(conn+pos+1,posP1-pos-1);
      int *tmp=new int[posP1-pos];
      for(unsigned i=0;i<nbOfSons;i++)
        {
          INTERP_KERNEL::NormalizedCellType cmsId;
          unsigned nbOfNodesSon=cm.fillSonCellNodalConnectivity2(i,conn+pos+1,posP1-pos-1,tmp,cmsId);
          const INTERP_KERNEL::CellModel& cms=INTERP_KERNEL::CellModel::getCellModel(cmsId);
          std::set<int> shareableCells(revNodalB[tmp[0]].begin(),revNodalB[tmp[0]].end());
          for(unsigned j=1;j<nbOfNodesSon && !shareableCells.empty();j++)
            {
              std::set<int> tmp2(revNodalB[tmp[j]].begin(),revNodalB[tmp[j]].end());
              std::set<int> tmp3;
              std::set_intersection(tmp2.begin(),tmp2.end(),shareableCells.begin(),shareableCells.end(),inserter(tmp3,tmp3.begin()));
              shareableCells=tmp3;
            }
          std::list<int> shareableCellsL(shareableCells.begin(),shareableCells.end());
          std::set<int> ref(tmp,tmp+nbOfNodesSon);
          for(std::list<int>::iterator iter=shareableCellsL.begin();iter!=shareableCellsL.end();)
            {
              if(cms.isCompatibleWith((INTERP_KERNEL::NormalizedCellType)meshDM1Type[*iter]))
                {
                  std::set<int> ref2(meshDM1Conn.begin()+meshDM1ConnIndex[*iter],meshDM1Conn.begin()+meshDM1ConnIndex[(*iter)+1]);
                  if(ref==ref2)
                    break;
                  else
                    iter=shareableCellsL.erase(iter);
                }
              else
                iter=shareableCellsL.erase(iter);
            }
          if(shareableCellsL.empty())
            {
              meshDM1Conn.insert(meshDM1Conn.end(),tmp,tmp+nbOfNodesSon);
              meshDM1ConnIndex.push_back(meshDM1ConnIndex.back()+nbOfNodesSon);
              int cellDM1Id=meshDM1Type.size();
              meshDM1Type.push_back((int)cmsId);
              for(unsigned k=0;k<nbOfNodesSon;k++)
                revNodalB[tmp[k]].push_back(cellDM1Id);
              revDescMeshConnB.resize(cellDM1Id+1);
              revDescMeshConnB.back().push_back(eltId);
              descMeshConnB[eltId].push_back(cellDM1Id);
            }
          else
            {
              int DM1cellId=shareableCellsL.front();
              revDescMeshConnB[DM1cellId].push_back(eltId);
              descMeshConnB[eltId].push_back(DM1cellId);
            }
        }
      delete [] tmp;
    }
  revNodalB.clear();
  //
  std::string name="Mesh constituent of "; name+=getName();
  MEDCouplingUMesh *ret=MEDCouplingUMesh::New(name.c_str(),getMeshDimension()-1);
  ret->setCoords(getCoords());
  int nbOfCellsInConstituent=meshDM1Type.size();
  ret->allocateCells(nbOfCellsInConstituent);
  revDescIndx->alloc(nbOfCellsInConstituent+1,1);
  int *tmp3=revDescIndx->getPointer(); tmp3[0]=0;
  for(int ii=0;ii<nbOfCellsInConstituent;ii++)
    {
      ret->insertNextCell((INTERP_KERNEL::NormalizedCellType)meshDM1Type[ii],meshDM1ConnIndex[ii+1]-meshDM1ConnIndex[ii],&meshDM1Conn[meshDM1ConnIndex[ii]]);
      tmp3[ii+1]=tmp3[ii]+revDescMeshConnB[ii].size();
    }
  ret->finishInsertingCells();
  revDesc->alloc(tmp3[nbOfCellsInConstituent],1);
  tmp3=revDesc->getPointer();
  for(std::vector< std::vector<int> >::const_iterator iter2=revDescMeshConnB.begin();iter2!=revDescMeshConnB.end();iter2++)
    tmp3=std::copy((*iter2).begin(),(*iter2).end(),tmp3);
  meshDM1Type.clear(); meshDM1ConnIndex.clear(); meshDM1Conn.clear();
  descIndx->alloc(nbOfCells+1,1);
  tmp3=descIndx->getPointer(); tmp3[0]=0;
  for(int jj=0;jj<nbOfCells;jj++)
    tmp3[jj+1]=tmp3[jj]+descMeshConnB[jj].size();
  desc->alloc(tmp3[nbOfCells],1);
  tmp3=desc->getPointer();
  for(std::vector< std::vector<int> >::const_iterator iter3=descMeshConnB.begin();iter3!=descMeshConnB.end();iter3++)
    tmp3=std::copy((*iter3).begin(),(*iter3).end(),tmp3);
  //
  return ret;
}

struct MEDCouplingAccVisit
{
  MEDCouplingAccVisit():_new_nb_of_nodes(0) { }
  int operator()(int val) { if(val!=-1) return _new_nb_of_nodes++; else return -1; }
  int _new_nb_of_nodes;
};

/*!
 * This method convert this into dynamic types without changing geometry.
 * That is to say if 'this' is a 2D, mesh after the invocation of this method it will contain only polygons.
 * If 'this' is a 3D mesh after the invocation of this method it will contain only polyhedra.
 * If mesh dimension is not in [2,3] an exception is thrown.
 * Of course pay attention that the resulting mesh is slower than previous one.
 * This method is above all designed to test more extensively algorithms able to deal with polygons/polyhedra.
 */
void MEDCouplingUMesh::convertToPolyTypes(const std::vector<int>& cellIdsToConvert)
{
  checkFullyDefined();
  int dim=getMeshDimension();
  if(dim<2 || dim>3)
    throw INTERP_KERNEL::Exception("Invalid mesh dimension : must be 2 or 3 !");
  if(dim==2)
    {
      const int *connIndex=_nodal_connec_index->getConstPointer();
      int *conn=_nodal_connec->getPointer();
      for(std::vector<int>::const_iterator iter=cellIdsToConvert.begin();iter!=cellIdsToConvert.end();iter++)
        conn[connIndex[*iter]]=INTERP_KERNEL::NORM_POLYGON;
    }
  else
    {
      int *connIndex=_nodal_connec_index->getPointer();
      int connIndexLgth=_nodal_connec_index->getNbOfElems();
      const int *connOld=_nodal_connec->getConstPointer();
      int connOldLgth=_nodal_connec->getNbOfElems();
      std::vector<int> connNew(connOld,connOld+connOldLgth);
      for(std::vector<int>::const_iterator iter=cellIdsToConvert.begin();iter!=cellIdsToConvert.end();iter++)
        {
          int pos=connIndex[*iter];
          int posP1=connIndex[(*iter)+1];
          int lgthOld=posP1-pos-1;
          const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::getCellModel((INTERP_KERNEL::NormalizedCellType)connNew[pos]);
          connNew[pos]=INTERP_KERNEL::NORM_POLYHED;
          unsigned nbOfFaces=cm.getNumberOfSons2(&connNew[pos+1],lgthOld);
          int *tmp=new int[nbOfFaces*lgthOld];
          int *work=tmp;
          for(int j=0;j<(int)nbOfFaces;j++)
            {
              INTERP_KERNEL::NormalizedCellType type;
              unsigned offset=cm.fillSonCellNodalConnectivity2(j,&connNew[pos+1],lgthOld,work,type);
              work+=offset;
              *work++=-1;
            }
          unsigned newLgth=work-tmp-1;
          unsigned delta=newLgth-lgthOld;
          std::transform(connIndex+(*iter)+1,connIndex+connIndexLgth,connIndex+(*iter)+1,std::bind2nd(std::plus<int>(),delta));
          connNew.insert(connNew.begin()+posP1,tmp+lgthOld,tmp+newLgth);
          std::copy(tmp,tmp+lgthOld,connNew.begin()+pos+1);
          delete [] tmp;
        }
      _nodal_connec->alloc(connNew.size(),1);
      int *newConnPtr=_nodal_connec->getPointer();
      std::copy(connNew.begin(),connNew.end(),newConnPtr);
    }
}

/*!
 * Array returned is the correspondance old to new.
 * The maximum value stored in returned array is the number of nodes of 'this' minus 1 after call of this method.
 * The size of returned array is the number of nodes of the old (previous to the call of this method) number of nodes.
 * -1 values in returned array means that the corresponding old node is no more used.
 */
DataArrayInt *MEDCouplingUMesh::zipCoordsTraducer()
{
  DataArrayInt *ret=DataArrayInt::New();
  int nbOfNodes=getNumberOfNodes();
  int spaceDim=getSpaceDimension();
  int *traducer=new int[nbOfNodes];
  std::fill(traducer,traducer+nbOfNodes,-1);
  ret->useArray(traducer,true,CPP_DEALLOC,nbOfNodes,1);
  int nbOfCells=getNumberOfCells();
  const int *connIndex=_nodal_connec_index->getConstPointer();
  int *conn=_nodal_connec->getPointer();
  for(int i=0;i<nbOfCells;i++)
    for(int j=connIndex[i]+1;j<connIndex[i+1];j++)
      if(conn[j]>=0)
        traducer[conn[j]]=1;
  int newNbOfNodes=std::count(traducer,traducer+nbOfNodes,1);
  std::transform(traducer,traducer+nbOfNodes,traducer,MEDCouplingAccVisit());
  for(int i=0;i<nbOfCells;i++)
    for(int j=connIndex[i]+1;j<connIndex[i+1];j++)
      if(conn[j]>=0)
        conn[j]=traducer[conn[j]];
  DataArrayDouble *newCoords=DataArrayDouble::New();
  double *newCoordsPtr=new double[newNbOfNodes*spaceDim];
  const double *oldCoordsPtr=_coords->getConstPointer();
  newCoords->useArray(newCoordsPtr,true,CPP_DEALLOC,newNbOfNodes,spaceDim);
  int *work=std::find_if(traducer,traducer+nbOfNodes,std::bind2nd(std::not_equal_to<int>(),-1));
  for(;work!=traducer+nbOfNodes;work=std::find_if(work,traducer+nbOfNodes,std::bind2nd(std::not_equal_to<int>(),-1)))
    {
      newCoordsPtr=std::copy(oldCoordsPtr+spaceDim*(work-traducer),oldCoordsPtr+spaceDim*(work-traducer+1),newCoordsPtr);
      work++;
    }
  setCoords(newCoords);
  newCoords->decrRef();
  return ret;
}

MEDCouplingPointSet *MEDCouplingUMesh::buildPartOfMySelf(const int *start, const int *end, bool keepCoords) const
{
  MEDCouplingUMesh *ret=buildPartOfMySelfKeepCoords(start,end);
  if(!keepCoords)
    ret->zipCoords();
  return ret;
}

/*!
 * Keeps from 'this' only cells which constituing point id are in the ids specified by ['start','end').
 * The return newly allocated mesh will share the same coordinates as 'this'.
 * Parameter 'fullyIn' specifies if a cell that has part of its nodes in ids array is kept or not.
 * If 'fullyIn' is true only cells whose ids are \b fully contained in ['start','end') tab will be kept.
 */
MEDCouplingPointSet *MEDCouplingUMesh::buildPartOfMySelfNode(const int *start, const int *end, bool fullyIn) const
{
  std::set<int> fastFinder(start,end);
  const int *conn=getNodalConnectivity()->getConstPointer();
  const int *connIndex=getNodalConnectivityIndex()->getConstPointer();
  int nbOfCells=getNumberOfCells();
  std::vector<int> cellIdsKept;
  for(int i=0;i<nbOfCells;i++)
    {
      std::set<int> connOfCell(conn+connIndex[i]+1,conn+connIndex[i+1]);
      connOfCell.erase(-1);//polyhedron separator
      int refLgth=std::min(connOfCell.size(),fastFinder.size());
      std::set<int> locMerge;
      std::insert_iterator< std::set<int> > it(locMerge,locMerge.begin());
      std::set_intersection(connOfCell.begin(),connOfCell.end(),fastFinder.begin(),fastFinder.end(),it);
      if((int)locMerge.size()==refLgth && fullyIn || locMerge.size()!=0 && !fullyIn)
        cellIdsKept.push_back(i);
    }
  return buildPartOfMySelf(&cellIdsKept[0],&cellIdsKept[0]+cellIdsKept.size(),true);
}

MEDCouplingPointSet *MEDCouplingUMesh::buildBoundaryMesh(bool keepCoords) const
{
  DataArrayInt *desc=DataArrayInt::New();
  DataArrayInt *descIndx=DataArrayInt::New();
  DataArrayInt *revDesc=DataArrayInt::New();
  DataArrayInt *revDescIndx=DataArrayInt::New();
  //
  MEDCouplingUMesh *meshDM1=buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
  revDesc->decrRef();
  desc->decrRef();
  descIndx->decrRef();
  int nbOfCells=meshDM1->getNumberOfCells();
  const int *revDescIndxC=revDescIndx->getConstPointer();
  std::vector<int> boundaryCells;
  for(int i=0;i<nbOfCells;i++)
    if(revDescIndxC[i+1]-revDescIndxC[i]==1)
      boundaryCells.push_back(i);
  revDescIndx->decrRef();
  MEDCouplingPointSet *ret=meshDM1->buildPartOfMySelf(&boundaryCells[0],&boundaryCells[0]+boundaryCells.size(),keepCoords);
  meshDM1->decrRef();
  return ret;
}

void MEDCouplingUMesh::findBoundaryNodes(std::vector<int>& nodes) const
{
  DataArrayInt *desc=DataArrayInt::New();
  DataArrayInt *descIndx=DataArrayInt::New();
  DataArrayInt *revDesc=DataArrayInt::New();
  DataArrayInt *revDescIndx=DataArrayInt::New();
  //
  MEDCouplingUMesh *meshDM1=buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
  revDesc->decrRef();
  desc->decrRef();
  descIndx->decrRef();
  std::set<int> ret;
  int nbOfCells=meshDM1->getNumberOfCells();
  const int *revDescIndxC=revDescIndx->getConstPointer();
  std::vector<int> boundaryCells;
  for(int i=0;i<nbOfCells;i++)
    if(revDescIndxC[i+1]-revDescIndxC[i]==1)
      boundaryCells.push_back(i);
  revDescIndx->decrRef();
  const int *conn=meshDM1->getNodalConnectivity()->getConstPointer();
  const int *connIndx=meshDM1->getNodalConnectivityIndex()->getConstPointer();
  for(std::vector<int>::const_iterator iter=boundaryCells.begin();iter!=boundaryCells.end();iter++)
    for(int k=connIndx[*iter]+1;k<connIndx[(*iter)+1];k++)
      ret.insert(conn[k]);
  nodes.resize(ret.size());
  std::copy(ret.begin(),ret.end(),nodes.begin());
  //
  meshDM1->decrRef();
}

void MEDCouplingUMesh::renumberConnectivity(const int *newNodeNumbers)
{
  int *conn=getNodalConnectivity()->getPointer();
  const int *connIndex=getNodalConnectivityIndex()->getConstPointer();
  int nbOfCells=getNumberOfCells();
  for(int i=0;i<nbOfCells;i++)
    for(int iconn=connIndex[i]+1;iconn!=connIndex[i+1];iconn++)
      {
        int& node=conn[iconn];
        if(node>=0)//avoid polyhedron separator
          {
            node=newNodeNumbers[node];
          }
      }
}

/*!
 * Given a boundary box 'bbox' returns elements 'elems' contained in this 'bbox'.
 * Warning 'elems' is incremented during the call so if elems is not empty before call returned elements will be
 * added in 'elems' parameter.
 */
void MEDCouplingUMesh::giveElemsInBoundingBox(const double *bbox, double eps, std::vector<int>& elems)
{
  int dim=getSpaceDimension();
  double* elem_bb=new double[2*dim];
  const int* conn      = getNodalConnectivity()->getConstPointer();
  const int* conn_index= getNodalConnectivityIndex()->getConstPointer();
  const double* coords = getCoords()->getConstPointer();
  int nbOfCells=getNumberOfCells();
  for ( int ielem=0; ielem<nbOfCells;ielem++ )
    {
      for (int i=0; i<dim; i++)
        {
          elem_bb[i*2]=std::numeric_limits<double>::max();
          elem_bb[i*2+1]=-std::numeric_limits<double>::max();
        }

      for (int inode=conn_index[ielem]+1; inode<conn_index[ielem+1]; inode++)//+1 due to offset of cell type.
        {
          int node= conn[inode];
          if(node>=0)//avoid polyhedron separator
            {
              for (int idim=0; idim<dim; idim++)
                {
                  if ( coords[node*dim+idim] < elem_bb[idim*2] )
                    {
                      elem_bb[idim*2] = coords[node*dim+idim] ;
                    }
                  if ( coords[node*dim+idim] > elem_bb[idim*2+1] )
                    {
                      elem_bb[idim*2+1] = coords[node*dim+idim] ;
                    }
                }
            }
        }
      if (intersectsBoundingBox(elem_bb, bbox, dim, eps))
        {
          elems.push_back(ielem);
        }
    }
  delete [] elem_bb;
}

INTERP_KERNEL::NormalizedCellType MEDCouplingUMesh::getTypeOfCell(int cellId) const
{
  const int *ptI=_nodal_connec_index->getConstPointer();
  const int *pt=_nodal_connec->getConstPointer();
  return (INTERP_KERNEL::NormalizedCellType) pt[ptI[cellId]];
}

int MEDCouplingUMesh::getNumberOfNodesInCell(int cellId) const
{
  const int *ptI=_nodal_connec_index->getConstPointer();
  const int *pt=_nodal_connec->getConstPointer();
  if(pt[ptI[cellId]]!=INTERP_KERNEL::NORM_POLYHED)
    return ptI[cellId+1]-ptI[cellId]-1;
  else
    return std::count_if(pt+ptI[cellId]+1,pt+ptI[cellId+1],std::bind2nd(std::not_equal_to<int>(),-1));
}

void MEDCouplingUMesh::setConnectivity(DataArrayInt *conn, DataArrayInt *connIndex, bool isComputingTypes)
{
  DataArrayInt::setArrayIn(conn,_nodal_connec);
  DataArrayInt::setArrayIn(connIndex,_nodal_connec_index);
  if(isComputingTypes)
    computeTypes();
}

MEDCouplingUMesh::MEDCouplingUMesh(const MEDCouplingUMesh& other, bool deepCpy):MEDCouplingPointSet(other,deepCpy),_iterator(-1),_mesh_dim(other._mesh_dim),
                                                                                _nodal_connec(0),_nodal_connec_index(0),
                                                                                _types(other._types)
{
  if(other._nodal_connec)
    _nodal_connec=other._nodal_connec->performCpy(deepCpy);
  if(other._nodal_connec_index)
    _nodal_connec_index=other._nodal_connec_index->performCpy(deepCpy);
}

MEDCouplingUMesh::~MEDCouplingUMesh()
{
  if(_nodal_connec)
    _nodal_connec->decrRef();
  if(_nodal_connec_index)
    _nodal_connec_index->decrRef();
}

void MEDCouplingUMesh::computeTypes()
{
  if(_nodal_connec && _nodal_connec_index)
    {
      _types.clear();
      const int *conn=_nodal_connec->getConstPointer();
      const int *connIndex=_nodal_connec_index->getConstPointer();
      int nbOfElem=_nodal_connec_index->getNbOfElems()-1;
      for(const int *pt=connIndex;pt!=connIndex+nbOfElem;pt++)
        _types.insert((INTERP_KERNEL::NormalizedCellType)conn[*pt]);
    }
}

/*!
 * This method checks that all arrays are set. If yes nothing done if no an exception is thrown.
 */
void MEDCouplingUMesh::checkFullyDefined() const throw(INTERP_KERNEL::Exception)
{
  if(!_nodal_connec_index || !_nodal_connec || !_coords)
    throw INTERP_KERNEL::Exception("Reverse nodal connectivity computation requires full connectivity and coordinates set in unstructured mesh.");
}

int MEDCouplingUMesh::getNumberOfCells() const
{ 
  if(_nodal_connec_index)
    if(_iterator==-1)
      return _nodal_connec_index->getNumberOfTuples()-1;
    else
      return _iterator;
  else
    if(_mesh_dim==-1)
      return 1;
    else
      throw INTERP_KERNEL::Exception("Unable to get number of cells because no connectivity specified !");
}

int MEDCouplingUMesh::getMeshDimension() const
{
  if(_mesh_dim<-1)
    throw INTERP_KERNEL::Exception("No mesh dimension specified !");
  return _mesh_dim;
}

int MEDCouplingUMesh::getMeshLength() const
{
  return _nodal_connec->getNbOfElems();
}

void MEDCouplingUMesh::getTinySerializationInformation(std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const
{
  MEDCouplingPointSet::getTinySerializationInformation(tinyInfo,littleStrings);
  tinyInfo.push_back(getMeshDimension());
  tinyInfo.push_back(getNumberOfCells());
  if(_nodal_connec)
    tinyInfo.push_back(getMeshLength());
  else
    tinyInfo.push_back(-1);
}

/*!
 * @param tinyInfo must be equal to the result given by getTinySerializationInformation method.
 */
void MEDCouplingUMesh::resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings)
{
  MEDCouplingPointSet::resizeForUnserialization(tinyInfo,a1,a2,littleStrings);
  if(tinyInfo[5]!=-1)
    a1->alloc(tinyInfo[5]+tinyInfo[4]+1,1);
}

void MEDCouplingUMesh::serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const
{
  MEDCouplingPointSet::serialize(a1,a2);
  if(getMeshDimension()>-1)
    {
      a1=DataArrayInt::New();
      a1->alloc(getMeshLength()+getNumberOfCells()+1,1);
      int *ptA1=a1->getPointer();
      const int *conn=getNodalConnectivity()->getConstPointer();
      const int *index=getNodalConnectivityIndex()->getConstPointer();
      ptA1=std::copy(index,index+getNumberOfCells()+1,ptA1);
      std::copy(conn,conn+getMeshLength(),ptA1);
    }
  else
    a1=0;
}

/*!
 * @param tinyInfo must be equal to the result given by getTinySerializationInformation method.
 */
void MEDCouplingUMesh::unserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, const std::vector<std::string>& littleStrings)
{
  MEDCouplingPointSet::unserialization(tinyInfo,a1,a2,littleStrings);
  setMeshDimension(tinyInfo[3]);
  if(tinyInfo[5]!=-1)
    {
      // Connectivity
      const int *recvBuffer=a1->getConstPointer();
      DataArrayInt* myConnecIndex=DataArrayInt::New();
      myConnecIndex->alloc(tinyInfo[4]+1,1);
      std::copy(recvBuffer,recvBuffer+tinyInfo[4]+1,myConnecIndex->getPointer());
      DataArrayInt* myConnec=DataArrayInt::New();
      myConnec->alloc(tinyInfo[5],1);
      std::copy(recvBuffer+tinyInfo[4]+1,recvBuffer+tinyInfo[4]+1+tinyInfo[5],myConnec->getPointer());
      setConnectivity(myConnec, myConnecIndex) ;
      myConnec->decrRef();
      myConnecIndex->decrRef();
    }
}

MEDCouplingUMesh *MEDCouplingUMesh::buildPartOfMySelfKeepCoords(const int *start, const int *end) const
{
  checkFullyDefined();
  MEDCouplingUMesh *ret=MEDCouplingUMesh::New();
  std::ostringstream stream; stream << PART_OF_NAME << getName();
  ret->setName(stream.str().c_str());
  ret->_mesh_dim=_mesh_dim;
  ret->setCoords(_coords);
  int nbOfElemsRet=end-start;
  int *connIndexRet=new int[nbOfElemsRet+1];
  connIndexRet[0]=0;
  const int *conn=_nodal_connec->getConstPointer();
  const int *connIndex=_nodal_connec_index->getConstPointer();
  int newNbring=0;
  for(const int *work=start;work!=end;work++,newNbring++)
    connIndexRet[newNbring+1]=connIndexRet[newNbring]+connIndex[*work+1]-connIndex[*work];
  int *connRet=new int[connIndexRet[nbOfElemsRet]];
  int *connRetWork=connRet;
  std::set<INTERP_KERNEL::NormalizedCellType> types;
  for(const int *work=start;work!=end;work++)
    {
      types.insert((INTERP_KERNEL::NormalizedCellType)conn[connIndex[*work]]);
      connRetWork=std::copy(conn+connIndex[*work],conn+connIndex[*work+1],connRetWork);
    }
  DataArrayInt *connRetArr=DataArrayInt::New();
  connRetArr->useArray(connRet,true,CPP_DEALLOC,connIndexRet[nbOfElemsRet],1);
  DataArrayInt *connIndexRetArr=DataArrayInt::New();
  connIndexRetArr->useArray(connIndexRet,true,CPP_DEALLOC,nbOfElemsRet+1,1);
  ret->setConnectivity(connRetArr,connIndexRetArr,false);
  ret->_types=types;
  connRetArr->decrRef();
  connIndexRetArr->decrRef();
  return ret;
}

/*!
 * brief returns the volumes of the cells underlying the field \a field
 *
 * For 2D geometries, the returned field contains the areas.
 * For 3D geometries, the returned field contains the volumes.
 *
 * param field field on which cells the volumes are required
 * return field containing the volumes, area or length depending the meshdimension.
 */
MEDCouplingFieldDouble *MEDCouplingUMesh::getMeasureField(bool isAbs) const
{
  int ipt;
  INTERP_KERNEL::NormalizedCellType type;
  int nbelem=getNumberOfCells();
  int dim_space=getSpaceDimension();
  const double *coords=getCoords()->getConstPointer();
  const int *connec=getNodalConnectivity()->getConstPointer();
  const int *connec_index=getNodalConnectivityIndex()->getConstPointer();

  MEDCouplingFieldDouble* field=MEDCouplingFieldDouble::New(ON_CELLS);
  DataArrayDouble* array=DataArrayDouble::New();
  array->alloc(nbelem,1);
  double *area_vol = array->getPointer();
  for(int iel=0;iel<nbelem;iel++)
    {
      ipt=connec_index[iel];
      type=(INTERP_KERNEL::NormalizedCellType)connec[ipt];
      area_vol[iel]=INTERP_KERNEL::computeVolSurfOfCell2<int,INTERP_KERNEL::ALL_C_MODE>(type,connec+ipt+1,connec_index[iel+1]-ipt-1,coords,dim_space);
    }
  if(isAbs)
    for(int iel=0;iel<nbelem;iel++)
      area_vol[iel]=fabs(area_vol[iel]);
  field->setArray(array) ;
  array->decrRef();
  field->setMesh(const_cast<MEDCouplingUMesh *>(this));  
  return field ;
}
