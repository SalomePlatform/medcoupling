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
#ifndef __MEDCOUPLINGNORMALIZEDUNSTRUCTUREDMESH_TXX__
#define __MEDCOUPLINGNORMALIZEDUNSTRUCTUREDMESH_TXX__

#include "MEDCouplingNormalizedUnstructuredMesh.hxx"

#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingMemArray.hxx"

#include <limits>

template<int SPACEDIM,int MESHDIM>
MEDCouplingNormalizedUnstructuredMesh<SPACEDIM,MESHDIM>::MEDCouplingNormalizedUnstructuredMesh(const ParaMEDMEM::MEDCouplingUMesh *mesh):_mesh(mesh)
{
  if(_mesh)
    _mesh->incrRef();
  prepare();
}

template<int SPACEDIM,int MESHDIM>
void MEDCouplingNormalizedUnstructuredMesh<SPACEDIM,MESHDIM>::getBoundingBox(double *boundingBox) const
{
  for(int i=0;i<SPACEDIM;i++)
    {
      boundingBox[i]=std::numeric_limits<double>::max();
      boundingBox[SPACEDIM+i]=-std::numeric_limits<double>::max();
    }
  const ParaMEDMEM::DataArrayDouble *array=_mesh->getCoords();
  const double *ptr=array->getConstPointer();
  int nbOfPts=array->getNbOfElems()/SPACEDIM;
  for(int j=0;j<SPACEDIM;j++)
    {
      const double *work=ptr+j;
      for(int i=0;i<nbOfPts;i++,work+=SPACEDIM)
        {
          if(boundingBox[j]>*work)
            boundingBox[j]=*work;
          if(boundingBox[j+SPACEDIM]<*work)
            boundingBox[j+SPACEDIM]=*work;
        }
    }
}

template<int SPACEDIM,int MESHDIM>
INTERP_KERNEL::NormalizedCellType MEDCouplingNormalizedUnstructuredMesh<SPACEDIM,MESHDIM>::getTypeOfElement(int eltId) const
{
  return _mesh->getTypeOfCell(eltId);
}

template<int SPACEDIM,int MESHDIM>
unsigned char MEDCouplingNormalizedUnstructuredMesh<SPACEDIM,MESHDIM>::getNumberOfNodesOfElement(int eltId) const
{
  return _mesh->getNumberOfNodesInCell(eltId);
}

template<int SPACEDIM,int MESHDIM>
unsigned long MEDCouplingNormalizedUnstructuredMesh<SPACEDIM,MESHDIM>::getNumberOfElements() const
{
  return _mesh->getNumberOfCells();
}

template<int SPACEDIM,int MESHDIM>
unsigned long MEDCouplingNormalizedUnstructuredMesh<SPACEDIM,MESHDIM>::getNumberOfNodes() const
{
  return _mesh->getNumberOfNodes();
}

template<int SPACEDIM,int MESHDIM>
const int *MEDCouplingNormalizedUnstructuredMesh<SPACEDIM,MESHDIM>::getConnectivityPtr() const
{
  return _conn_for_interp;
}

template<int SPACEDIM,int MESHDIM>
const double *MEDCouplingNormalizedUnstructuredMesh<SPACEDIM,MESHDIM>::getCoordinatesPtr() const
{
  const ParaMEDMEM::DataArrayDouble *array=_mesh->getCoords();
  return array->getConstPointer();
}

template<int SPACEDIM,int MESHDIM>
const int *MEDCouplingNormalizedUnstructuredMesh<SPACEDIM,MESHDIM>::getConnectivityIndexPtr() const
{
  return _conn_index_for_interp;
}

template<int SPACEDIM,int MESHDIM>
void MEDCouplingNormalizedUnstructuredMesh<SPACEDIM,MESHDIM>::releaseTempArrays()
{
  delete [] _conn_for_interp;
  delete [] _conn_index_for_interp;
  _conn_for_interp=0;
  _conn_index_for_interp=0;
}

template<int SPACEDIM,int MESHDIM>
MEDCouplingNormalizedUnstructuredMesh<SPACEDIM,MESHDIM>::~MEDCouplingNormalizedUnstructuredMesh()
{
  if(_mesh)
    ((ParaMEDMEM::MEDCouplingUMesh *)_mesh)->decrRef();
  releaseTempArrays();
}

template<int SPACEDIM,int MESHDIM>
void MEDCouplingNormalizedUnstructuredMesh<SPACEDIM,MESHDIM>::prepare()
{
  int nbOfCell=_mesh->getNumberOfCells();
  int initialConnSize=_mesh->getNodalConnectivity()->getNbOfElems();
  _conn_for_interp=new int[initialConnSize-nbOfCell];
  _conn_index_for_interp=new int[nbOfCell+1];
  _conn_index_for_interp[0]=0;
  const int *work_conn=_mesh->getNodalConnectivity()->getConstPointer()+1;
  const int *work_conn_index=_mesh->getNodalConnectivityIndex()->getConstPointer();
  int *work_conn_for_interp=_conn_for_interp;
  int *work_conn_index_for_interp=_conn_index_for_interp;
  for(int i=0;i<nbOfCell;i++)
    {
      int nbOfValsToCopy=work_conn_index[1]-work_conn_index[0]-1;
      work_conn_for_interp=std::copy(work_conn,work_conn+nbOfValsToCopy,work_conn_for_interp);
      work_conn_index_for_interp[1]=work_conn_index_for_interp[0]+nbOfValsToCopy;
      work_conn_index++;
      work_conn+=nbOfValsToCopy+1;
      work_conn_index_for_interp++;
    }
}

#endif
