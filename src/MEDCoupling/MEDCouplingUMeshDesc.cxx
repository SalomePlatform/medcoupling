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
#include "MEDCouplingUMeshDesc.hxx"
#include "CellModel.hxx"
#include "MemArray.hxx"

using namespace ParaMEDMEM;

MEDCouplingUMeshDesc::MEDCouplingUMeshDesc():_mesh_dim(-1),_desc_connec(0),_desc_connec_index(0),
                                             _nodal_connec_face(0),_nodal_connec_face_index(0)
{
}

MEDCouplingUMeshDesc::~MEDCouplingUMeshDesc()
{
  if(_desc_connec)
    _desc_connec->decrRef();
  if(_desc_connec_index)
    _desc_connec_index->decrRef();
  if(_nodal_connec_face)
    _nodal_connec_face->decrRef();
  if(_nodal_connec_face_index)
    _nodal_connec_face_index->decrRef();
}

MEDCouplingUMeshDesc *MEDCouplingUMeshDesc::New()
{
  return new MEDCouplingUMeshDesc;
}

void MEDCouplingUMeshDesc::checkCoherency() const throw(INTERP_KERNEL::Exception)
{
  for(std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator iter=_types.begin();iter!=_types.end();iter++)
    {
      if(INTERP_KERNEL::CellModel::getCellModel(*iter).getDimension()!=_mesh_dim)
        {
          std::ostringstream message;
          message << "MeshDesc invalid because dimension is " << _mesh_dim << " and there is presence of cell(s) with type " << (*iter);
          throw INTERP_KERNEL::Exception(message.str().c_str());
        }
    }
}

void MEDCouplingUMeshDesc::setMeshDimension(unsigned meshDim)
{
  _mesh_dim=meshDim;
  declareAsNew();
}

int MEDCouplingUMeshDesc::getNumberOfCells() const
{
  if(_desc_connec_index)
    return _desc_connec_index->getNumberOfTuples()-1;
  else
    throw INTERP_KERNEL::Exception("Unable to get number of cells because no connectivity specified !");
}

int MEDCouplingUMeshDesc::getNumberOfFaces() const
{
  if(_nodal_connec_face_index)
    return _nodal_connec_face_index->getNumberOfTuples()-1;
  else
    throw INTERP_KERNEL::Exception("Unable to get number of faces because no connectivity specified !");
}

int MEDCouplingUMeshDesc::getCellMeshLength() const
{
  return _desc_connec->getNbOfElems();
}

int MEDCouplingUMeshDesc::getFaceMeshLength() const
{
  return _nodal_connec_face->getNbOfElems();
}

void MEDCouplingUMeshDesc::setConnectivity(DataArrayInt *descConn, DataArrayInt *descConnIndex, DataArrayInt *nodalFaceConn, DataArrayInt *nodalFaceConnIndx)
{
  DataArrayInt::setArrayIn(descConn,_desc_connec);
  DataArrayInt::setArrayIn(descConnIndex,_desc_connec_index);
  DataArrayInt::setArrayIn(nodalFaceConn,_nodal_connec_face);
  DataArrayInt::setArrayIn(nodalFaceConnIndx,_nodal_connec_face_index);
  computeTypes();
}

void MEDCouplingUMeshDesc::getTinySerializationInformation(std::vector<int>& tinyInfo) const
{
  tinyInfo.resize(7);
  tinyInfo[0]=getSpaceDimension();
  tinyInfo[1]=getMeshDimension();
  tinyInfo[2]=getNumberOfNodes();
  tinyInfo[3]=getNumberOfCells();
  tinyInfo[4]=getCellMeshLength();
  tinyInfo[5]=getNumberOfFaces();
  tinyInfo[6]=getFaceMeshLength();
}

void MEDCouplingUMeshDesc::resizeForSerialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2)
{
  a1->alloc(tinyInfo[4]+tinyInfo[3]+1+tinyInfo[6]+tinyInfo[5]+1,1);
  a2->alloc(tinyInfo[2],tinyInfo[0]);
}

void MEDCouplingUMeshDesc::serialize(DataArrayInt *&a1, DataArrayDouble *&a2)
{
  a1=DataArrayInt::New();
  a1->alloc(getCellMeshLength()+getNumberOfCells()+1+getFaceMeshLength()+getNumberOfFaces()+1,1);
  int *ptA1=a1->getPointer();
  const int *descConn=_desc_connec->getConstPointer();
  const int *descConnIndex=_desc_connec_index->getConstPointer();
  const int *faceConn=_nodal_connec_face->getConstPointer();
  const int *faceConnIndex=_nodal_connec_face_index->getConstPointer();
  ptA1=std::copy(descConn,descConn+getCellMeshLength(),ptA1);
  ptA1=std::copy(descConnIndex,descConnIndex+getNumberOfCells()+1,ptA1);
  ptA1=std::copy(faceConn,faceConn+getFaceMeshLength(),ptA1);
  std::copy(faceConnIndex,faceConnIndex+getNumberOfFaces()+1,ptA1);
  a2=getCoords();
  a2->incrRef();
}

MEDCouplingPointSet *MEDCouplingUMeshDesc::buildObjectFromUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2)
{
  MEDCouplingUMeshDesc *meshing=MEDCouplingUMeshDesc::New();
  meshing->setCoords(a2);
  const int *recvBuffer=a1->getConstPointer();
  DataArrayInt *descConn=DataArrayInt::New();
  descConn->alloc(tinyInfo[4],1);
  std::copy(recvBuffer,recvBuffer+tinyInfo[4],descConn->getPointer());
  DataArrayInt *descConnIndex=DataArrayInt::New();
  descConnIndex->alloc(tinyInfo[3]+1,1);
  std::copy(recvBuffer+tinyInfo[4],recvBuffer+tinyInfo[4]+tinyInfo[3]+1,descConnIndex->getPointer());
  DataArrayInt *faceConn=DataArrayInt::New();
  faceConn->alloc(tinyInfo[6],1);
  std::copy(recvBuffer+tinyInfo[4]+tinyInfo[3]+1,recvBuffer+tinyInfo[4]+tinyInfo[3]+1+tinyInfo[6],faceConn->getPointer());
  DataArrayInt *faceConnIndex=DataArrayInt::New();
  faceConnIndex->alloc(tinyInfo[5]+1,1);
  std::copy(recvBuffer+tinyInfo[4]+tinyInfo[3]+1+tinyInfo[6],
            recvBuffer+tinyInfo[4]+tinyInfo[3]+1+tinyInfo[6]+tinyInfo[5]+1,faceConnIndex->getPointer());
  meshing->setConnectivity(descConn,descConnIndex,faceConn,faceConnIndex);
  descConn->decrRef();
  descConnIndex->decrRef();
  faceConn->decrRef();
  faceConnIndex->decrRef();
  meshing->setMeshDimension(tinyInfo[1]);
  return meshing;
}

void MEDCouplingUMeshDesc::giveElemsInBoundingBox(const double *bbox, double eps, std::vector<int>& elems)
{
  int dim=getSpaceDimension();
  double* elem_bb=new double[2*dim];
  const int* conn      = _desc_connec->getConstPointer();
  const int* conn_index= _desc_connec_index->getConstPointer();
  const int* face      = _nodal_connec_face->getConstPointer();
  const int* face_index= _nodal_connec_face_index->getConstPointer();
  const double* coords = getCoords()->getConstPointer();
  int nbOfCells=getNumberOfCells();
  for ( int ielem=0; ielem<nbOfCells;ielem++ )
    {
      for (int i=0; i<dim; i++)
        {
          elem_bb[i*2]=std::numeric_limits<double>::max();
          elem_bb[i*2+1]=-std::numeric_limits<double>::max();
        }

      for (int iface=conn_index[ielem]+1; iface<conn_index[ielem+1]; iface++)//+1 due to offset of cell type.
        {
          for(int inode=face_index[iface]+1;inode<face_index[iface+1];inode++)
            {
              int node=face[inode];
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

MEDCouplingPointSet *MEDCouplingUMeshDesc::buildPartOfMySelf(const int *start, const int *end, bool keepCoords) const
{
  //not implemented yet.
  return 0;
}

MEDCouplingFieldDouble *MEDCouplingUMeshDesc::getMeasureField() const
{
  //not implemented yet.
  return 0;
}

void MEDCouplingUMeshDesc::computeTypes()
{
  if(_desc_connec && _desc_connec_index)
    {
      _types.clear();
      const int *conn=_desc_connec->getConstPointer();
      const int *connIndex=_desc_connec_index->getConstPointer();
      int nbOfElem=_desc_connec_index->getNbOfElems()-1;
      for(const int *pt=connIndex;pt!=connIndex+nbOfElem;pt++)
        _types.insert((INTERP_KERNEL::NormalizedCellType)conn[*pt]);
    }
}
