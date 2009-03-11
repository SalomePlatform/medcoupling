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
#include "CellModel.hxx"

#include <sstream>

using namespace ParaMEDMEM;

const char MEDCouplingUMesh::PART_OF_NAME[]="PartOf_";

MEDCouplingUMesh *MEDCouplingUMesh::New()
{
  return new MEDCouplingUMesh;
}

MEDCouplingUMesh *MEDCouplingUMesh::clone(bool recDeepCpy) const
{
  return new MEDCouplingUMesh(*this,recDeepCpy);
}

void MEDCouplingUMesh::updateTime()
{
  if(_nodal_connec)
    {
      updateTimeWith(*_nodal_connec);
    }
  if(_nodal_connec_index)
    {
      updateTimeWith(*_nodal_connec_index);
    }
  if(_coords)
    {
      updateTimeWith(*_coords);
    }
}

MEDCouplingUMesh::MEDCouplingUMesh():_iterator(-1),_mesh_dim(-1),
                                     _nodal_connec(0),_nodal_connec_index(0),_coords(0)
{
}

void MEDCouplingUMesh::checkCoherency() const throw(INTERP_KERNEL::Exception)
{
  for(std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator iter=_types.begin();iter!=_types.end();iter++)
    {
      if(INTERP_KERNEL::CellModel::getCellModel(*iter).getDimension()!=_mesh_dim)
        {
          std::ostringstream message;
          message << "Mesh invalid because dimension is " << _mesh_dim << " and there is presence of cell(s) with type " << (*iter);
          throw INTERP_KERNEL::Exception(message.str().c_str());
        }
    }
}

void MEDCouplingUMesh::setMeshDimension(unsigned meshDim)
{
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

void MEDCouplingUMesh::setCoords(DataArrayDouble *coords)
{
  if( coords != _coords )
    {
      if (_coords)
        _coords->decrRef();
      _coords=coords;
      if(_coords)
        _coords->incrRef();
      declareAsNew();
    }
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
  int *pt=_nodal_connec_index->getPointer();
  int idx=pt[_iterator];

  _nodal_connec->reAlloc(idx);
  _nodal_connec_index->reAlloc(_iterator+1);
  _iterator=-1;
}

bool MEDCouplingUMesh::isEqual(const MEDCouplingMesh *other, double prec) const
{
  checkFullyDefined();
  const MEDCouplingUMesh *otherC=dynamic_cast<const MEDCouplingUMesh *>(other);
  if(!otherC)
    return false;
  otherC->checkFullyDefined();
  if(!MEDCouplingMesh::isEqual(other,prec))
    return false;
  if(_mesh_dim!=otherC->_mesh_dim)
    return false;
  if(_types!=otherC->_types)
    return false;
  if(!_coords->isEqual(otherC->_coords,prec))
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
  const int *conn=_nodal_connec->getPointer();
  const int *connIndex=_nodal_connec_index->getPointer();
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

void MEDCouplingUMesh::zipCoords()
{
  checkFullyDefined();
  DataArrayInt *traducer=zipCoordsTraducer();
  traducer->decrRef();
}

struct MEDCouplingAccVisit
{
  MEDCouplingAccVisit():_new_nb_of_nodes(0) { }
  int operator()(int val) { if(val!=-1) return _new_nb_of_nodes++; else return -1; }
  int _new_nb_of_nodes;
};

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
  const int *connIndex=_nodal_connec_index->getPointer();
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
  const double *oldCoordsPtr=_coords->getPointer();
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

MEDCouplingUMesh *MEDCouplingUMesh::buildPartOfMySelf(const int *start, const int *end, bool keepCoords) const
{
  MEDCouplingUMesh *ret=buildPartOfMySelfKeepCoords(start,end);
  if(!keepCoords)
    ret->zipCoords();
  return ret;
}

INTERP_KERNEL::NormalizedCellType MEDCouplingUMesh::getTypeOfCell(int cellId) const
{
  int *ptI=_nodal_connec_index->getPointer();
  int *pt=_nodal_connec->getPointer();
  return (INTERP_KERNEL::NormalizedCellType) pt[ptI[cellId]];
}

int MEDCouplingUMesh::getNumberOfNodesInCell(int cellId) const
{
  int *ptI=_nodal_connec_index->getPointer();
  return ptI[cellId+1]-ptI[cellId]-1;
}

void MEDCouplingUMesh::setConnectivity(DataArrayInt *conn, DataArrayInt *connIndex, bool isComputingTypes)
{
  if(_nodal_connec!=conn)
    {
      if(_nodal_connec)
        _nodal_connec->decrRef();
      _nodal_connec=conn;
      if(_nodal_connec)
        _nodal_connec->incrRef();
    }
  if(_nodal_connec_index!=connIndex)
    {
      if(_nodal_connec_index)
        _nodal_connec_index->decrRef();
      _nodal_connec_index=connIndex;
      if(_nodal_connec_index)
        _nodal_connec_index->incrRef();
    }
  if(isComputingTypes)
    computeTypes();
}

MEDCouplingUMesh::MEDCouplingUMesh(const MEDCouplingUMesh& other, bool deepCpy):MEDCouplingMesh(other),_iterator(-1),_mesh_dim(other._mesh_dim),
                                                                                _nodal_connec(0),_nodal_connec_index(0),_coords(0),
                                                                                _types(other._types)
{
  if(other._nodal_connec)
    _nodal_connec=other._nodal_connec->performCpy(deepCpy);
  if(other._nodal_connec_index)
    _nodal_connec_index=other._nodal_connec_index->performCpy(deepCpy);
  if(other._coords)
    _coords=other._coords->performCpy(deepCpy);
}

MEDCouplingUMesh::~MEDCouplingUMesh()
{
  if(_nodal_connec)
    _nodal_connec->decrRef();
  if(_nodal_connec_index)
    _nodal_connec_index->decrRef();
  if(_coords)
    _coords->decrRef();
}

void MEDCouplingUMesh::computeTypes()
{
  if(_nodal_connec && _nodal_connec_index)
    {
      _types.clear();
      const int *conn=_nodal_connec->getPointer();
      const int *connIndex=_nodal_connec_index->getPointer();
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

bool MEDCouplingUMesh::isStructured() const
{
  return false;
}

int MEDCouplingUMesh::getNumberOfCells() const
{ 
  if(_nodal_connec_index)
    if(_iterator==-1)
      return _nodal_connec_index->getNumberOfTuples()-1;
    else
      return _iterator;
  else
    throw INTERP_KERNEL::Exception("Unable to get number of cells because no coordinates specified !");
}

int MEDCouplingUMesh::getNumberOfNodes() const
{
  if(_coords)
    return _coords->getNumberOfTuples();
  else
    throw INTERP_KERNEL::Exception("Unable to get number of nodes because no coordinates specified !");
}

int MEDCouplingUMesh::getSpaceDimension() const
{
  if(_coords)
    return _coords->getNumberOfComponents();
  else
    throw INTERP_KERNEL::Exception("Unable to get space dimension because no coordinates specified !");
}

int MEDCouplingUMesh::getMeshLength() const
{
  return _nodal_connec->getNbOfElems();
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
  const int *conn=_nodal_connec->getPointer();
  const int *connIndex=_nodal_connec_index->getPointer();
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

