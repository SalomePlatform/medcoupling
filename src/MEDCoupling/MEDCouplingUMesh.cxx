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

MEDCouplingUMesh *MEDCouplingUMesh::New()
{
 return new MEDCouplingUMesh;
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
