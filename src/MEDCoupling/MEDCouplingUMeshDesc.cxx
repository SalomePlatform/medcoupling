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

#include "MEDCouplingUMeshDesc.hxx"
#include "CellModel.hxx"
#include "MEDCouplingMemArray.hxx"

#include <limits>

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

MEDCouplingUMeshDesc *MEDCouplingUMeshDesc::New(const char *meshName, int meshDim)
{
  MEDCouplingUMeshDesc *ret=new MEDCouplingUMeshDesc;
  ret->setName(meshName);
  ret->setMeshDimension(meshDim);
  return ret;
}

/*!
 * not implemented
 */
MEDCouplingMesh *MEDCouplingUMeshDesc::deepCpy() const
{
  return 0;
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

void MEDCouplingUMeshDesc::checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                DataArrayInt *&cellCor, DataArrayInt *&nodeCor) const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("MEDCouplingUMeshDesc::checkDeepEquivalWith : not implemented yet !");
}

void MEDCouplingUMeshDesc::checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                           DataArrayInt *&cellCor) const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("MEDCouplingUMeshDesc::checkDeepEquivalOnSameNodesWith : not implemented yet !");
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

INTERP_KERNEL::NormalizedCellType MEDCouplingUMeshDesc::getTypeOfCell(int cellId) const
{
  const int *desc_connec=_desc_connec->getConstPointer();
  const int *desc_connec_index=_desc_connec_index->getConstPointer();
  return (INTERP_KERNEL::NormalizedCellType)desc_connec[desc_connec_index[cellId]+1];
}

int MEDCouplingUMeshDesc::getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const
{
  const int *desc_connec=_desc_connec->getConstPointer();
  const int *desc_connec_index=_desc_connec_index->getConstPointer();
  int nbOfCells=getNumberOfCells();
  int ret=0;
  for(int i=0;i<nbOfCells;i++)
    if((INTERP_KERNEL::NormalizedCellType) desc_connec[desc_connec_index[i]]==type)
      ret++;
  return ret;
}

void MEDCouplingUMeshDesc::getNodeIdsOfCell(int cellId, std::vector<int>& conn) const
{
  //not implemented yet.
}

void MEDCouplingUMeshDesc::getCoordinatesOfNode(int nodeId, std::vector<double>& coo) const
{
  //not implemented yet.
}

std::string MEDCouplingUMeshDesc::simpleRepr() const
{
  std::string ret("Unstructured mesh with descending connectivity : ");
  ret+=getName();
  ret+="\n";
  return ret;
}

std::string MEDCouplingUMeshDesc::advancedRepr() const
{
  std::string ret("Unstructured mesh with descending connectivity : ");
  ret+=getName();
  ret+="\n";
  return ret;
}

void MEDCouplingUMeshDesc::setConnectivity(DataArrayInt *descConn, DataArrayInt *descConnIndex, DataArrayInt *nodalFaceConn, DataArrayInt *nodalFaceConnIndx)
{
  DataArrayInt::SetArrayIn(descConn,_desc_connec);
  DataArrayInt::SetArrayIn(descConnIndex,_desc_connec_index);
  DataArrayInt::SetArrayIn(nodalFaceConn,_nodal_connec_face);
  DataArrayInt::SetArrayIn(nodalFaceConnIndx,_nodal_connec_face_index);
  computeTypes();
}

void MEDCouplingUMeshDesc::getTinySerializationInformation(std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const
{
  MEDCouplingPointSet::getTinySerializationInformation(tinyInfo,littleStrings);
  tinyInfo.push_back(getMeshDimension());
  tinyInfo.push_back(getNumberOfNodes());
  tinyInfo.push_back(getNumberOfCells());
  tinyInfo.push_back(getCellMeshLength());
  tinyInfo.push_back(getNumberOfFaces());
  tinyInfo.push_back(getFaceMeshLength());
}

bool MEDCouplingUMeshDesc::isEmptyMesh(const std::vector<int>& tinyInfo) const
{
  return tinyInfo[5]<=0;
}

void MEDCouplingUMeshDesc::resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings)
{
  std::vector<int> tinyInfoTmp(tinyInfo.begin()+1,tinyInfo.end());
  MEDCouplingPointSet::resizeForUnserialization(tinyInfoTmp,a1,a2,littleStrings);
  a1->alloc(tinyInfo[5]+tinyInfo[4]+1+tinyInfo[7]+tinyInfo[6]+1,1);
}

void MEDCouplingUMeshDesc::serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const
{
  MEDCouplingPointSet::serialize(a1,a2);
  //
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
}

void MEDCouplingUMeshDesc::unserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, const std::vector<std::string>& littleStrings)
{
  std::vector<int> tinyInfoTmp(tinyInfo.begin()+1,tinyInfo.end());
  MEDCouplingPointSet::unserialization(tinyInfoTmp,a1,a2,littleStrings);
  //
  const int *recvBuffer=a1->getConstPointer();
  DataArrayInt *descConn=DataArrayInt::New();
  descConn->alloc(tinyInfo[5],1);
  std::copy(recvBuffer,recvBuffer+tinyInfo[5],descConn->getPointer());
  DataArrayInt *descConnIndex=DataArrayInt::New();
  descConnIndex->alloc(tinyInfo[4]+1,1);
  std::copy(recvBuffer+tinyInfo[5],recvBuffer+tinyInfo[5]+tinyInfo[4]+1,descConnIndex->getPointer());
  DataArrayInt *faceConn=DataArrayInt::New();
  faceConn->alloc(tinyInfo[7],1);
  std::copy(recvBuffer+tinyInfo[5]+tinyInfo[4]+1,recvBuffer+tinyInfo[5]+tinyInfo[4]+1+tinyInfo[7],faceConn->getPointer());
  DataArrayInt *faceConnIndex=DataArrayInt::New();
  faceConnIndex->alloc(tinyInfo[6]+1,1);
  std::copy(recvBuffer+tinyInfo[5]+tinyInfo[4]+1+tinyInfo[7],
            recvBuffer+tinyInfo[5]+tinyInfo[5]+1+tinyInfo[7]+tinyInfo[6]+1,faceConnIndex->getPointer());
  setConnectivity(descConn,descConnIndex,faceConn,faceConnIndex);
  descConn->decrRef();
  descConnIndex->decrRef();
  faceConn->decrRef();
  faceConnIndex->decrRef();
  setMeshDimension(tinyInfo[2]);
}

void MEDCouplingUMeshDesc::getCellsInBoundingBox(const double *bbox, double eps, std::vector<int>& elems)
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

      for (int jface=conn_index[ielem]+1; jface<conn_index[ielem+1]; jface++)//+1 due to offset of cell type.
        {
          int iface=conn[jface];
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

void MEDCouplingUMeshDesc::getCellsInBoundingBox(const INTERP_KERNEL::DirectedBoundingBox &bbox, double eps, std::vector<int>& elems)
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

      for (int jface=conn_index[ielem]+1; jface<conn_index[ielem+1]; jface++)//+1 due to offset of cell type.
        {
          int iface=conn[jface];
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
      if (intersectsBoundingBox(bbox, elem_bb, dim, eps))
        {
          elems.push_back(ielem);
        }
    }
  delete [] elem_bb;
}

DataArrayInt *MEDCouplingUMeshDesc::mergeNodes(double precision, bool& areNodesMerged, int& newNbOfNodes)
{
  //not implemented yet.
  areNodesMerged=false;
  return 0;
}

DataArrayInt *MEDCouplingUMeshDesc::mergeNodes2(double precision, bool& areNodesMerged, int& newNbOfNodes)
{
  //not implemented yet.
  areNodesMerged=false;
  return 0;
}

void MEDCouplingUMeshDesc::tryToShareSameCoordsPermute(const MEDCouplingPointSet& other, double epsilon) throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("Not implemented yet !");
}

MEDCouplingPointSet *MEDCouplingUMeshDesc::buildPartOfMySelf(const int *start, const int *end, bool keepCoords) const
{
  //not implemented yet.
  return 0;
}

MEDCouplingPointSet *MEDCouplingUMeshDesc::buildPartOfMySelfNode(const int *start, const int *end, bool fullyIn) const
{
  //not implemented yet
  return 0;
}

MEDCouplingPointSet *MEDCouplingUMeshDesc::buildFacePartOfMySelfNode(const int *start, const int *end, bool fullyIn) const
{
  //not implemented yet
  return 0;
}

DataArrayInt *MEDCouplingUMeshDesc::simplexize(int policy) throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("MEDCouplingUMeshDesc::simplexize : Not implemented yet !");
}

void MEDCouplingUMeshDesc::findBoundaryNodes(std::vector<int>& nodes) const
{
  //not implemented yet
}

MEDCouplingPointSet *MEDCouplingUMeshDesc::buildBoundaryMesh(bool keepCoords) const
{
  //not implemented yet
  return 0;
}

MEDCouplingUMesh *MEDCouplingUMeshDesc::buildUnstructured() const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("MEDCouplingUMeshDesc::buildUnstructured : not implemented yet !");
}

void MEDCouplingUMeshDesc::renumberCells(const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("Available for UMesh desc but not implemented yet !");
}

void MEDCouplingUMeshDesc::renumberNodes(const int *newNodeNumbers, int newNbOfNodes)
{
  MEDCouplingPointSet::renumberNodes(newNodeNumbers,newNbOfNodes);
  //not implemented yet
}

MEDCouplingFieldDouble *MEDCouplingUMeshDesc::getMeasureField(bool isAbs) const
{
  //not implemented yet.
  return 0;
}

MEDCouplingFieldDouble *MEDCouplingUMeshDesc::getMeasureFieldOnNode(bool isAbs) const
{
  //not implemented yet.
  return 0;
}

MEDCouplingFieldDouble *MEDCouplingUMeshDesc::buildOrthogonalField() const
{
  if(getMeshDimension()!=2)
    throw INTERP_KERNEL::Exception("Expected a cmesh with meshDim == 2 !");
  //not implemented yet !
  return 0;
}

DataArrayInt *MEDCouplingUMeshDesc::zipCoordsTraducer()
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

void MEDCouplingUMeshDesc::checkFullyDefined() const throw(INTERP_KERNEL::Exception)
{
  if(!_desc_connec || !_desc_connec_index || !_nodal_connec_face || !_nodal_connec_face_index || !_coords)
    throw INTERP_KERNEL::Exception("full connectivity and coordinates not set in unstructured mesh.");
}

MEDCouplingMesh *MEDCouplingUMeshDesc::mergeMyselfWith(const MEDCouplingMesh *other) const
{  
  //not implemented yet.
  return 0;
}

DataArrayDouble *MEDCouplingUMeshDesc::getBarycenterAndOwner() const
{
  //not implemented yet.
  return 0;
}

int MEDCouplingUMeshDesc::getCellContainingPoint(const double *pos, double eps) const
{
  //not implemented yet.
  return -1;
}
